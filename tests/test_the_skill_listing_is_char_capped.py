"""The skill catalogue in the `skill` tool description is char-capped.

Phase 6 of the learning wave. The listing pasted into the `skill` tool's
description carried two unbounded costs:

  * per-skill text was capped only in the description (70 chars) while the
    name was untrimmed, and there was no cap on the LISTING as a whole:
    40 long-named skills cost over 4,000 chars (~1,000 tokens) of schema
    per request, re-sent every turn;
  * the cap was a silent ``[:40]``: skills 41..N vanished without a word,
    so with a growing (learning) catalogue the model did not even know
    there were more to ask for.

The staged-loading contract: name + a shortened description per skill, a
total character ceiling for the listing, an explicit overflow notice when
the catalogue does not fit ("N more — call the skill tool with a name; a
miss lists the full catalogue" — which stays true because the executor
already answers a miss with every available name), and the full body
loaded only via the `skill` call itself.

Budget rule measured in phase 1: the listing at 200 skills must cost no
more than the listing at 40 plus the ceiling — the catalogue's growth is
paid for by the ceiling, never by the schema.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

from delfin.agent import api_client as A
from delfin.agent import skills as S
from delfin.agent.skills import Skill


# ---------------------------------------------------------------------------
# The listing itself
# ---------------------------------------------------------------------------

def _fake_skills(n: int) -> list[Skill]:
    """Skills with long names and long descriptions, the worst case."""
    return [
        Skill(
            name=f"a-rather-long-skill-name-{i:03d}",
            description=("does something quite specific, in detail " * 8),
            body="body",
            source=pathlib.Path(f"/nonexistent/skill-{i}.md"),
        )
        for i in range(n)
    ]


def test_every_skill_name_is_listed_with_a_capped_description():
    listing = A._skill_listing(_fake_skills(10))
    assert "a-rather-long-skill-name-000" in listing
    # Per-skill: the description is trimmed hard, not the 70-char legacy cut.
    for part in listing.split("; "):
        if "—" in part:
            assert len(part.split("—", 1)[1].strip()) <= 70


def test_the_listing_has_a_total_character_ceiling():
    """40 worst-case skills must fit a ceiling a fraction of the old cost."""
    listing = A._skill_listing(_fake_skills(40))
    assert len(listing) <= 2_000


def test_growth_beyond_the_ceiling_is_announced_not_silent():
    """200 skills: the listing stays capped and names the remainder."""
    skills = _fake_skills(200)
    listing = A._skill_listing(skills)
    assert len(listing) <= 2_000
    _counted(listing, skills)  # listed + announced == 200, exactly
    # The overflow notice says how many more exist and how to reach them.
    assert re.search(r"\d+ more — call the skill tool", listing)
    # The notice must not lie: the executor answers a miss with the full
    # catalogue, so the wording "a miss lists" has to hold. That is the
    # _execute_skill contract, pinned by its own test below.


def test_two_hundred_skills_cost_no_more_than_forty_plus_the_ceiling():
    forty = A._skill_listing(_fake_skills(40))
    two_hundred = A._skill_listing(_fake_skills(200))
    assert len(two_hundred) <= len(forty) + 2_000


def test_a_small_catalogue_has_no_overflow_notice():
    listing = A._skill_listing(_fake_skills(10))
    assert "more" not in listing
    assert "a-rather-long-skill-name-009" in listing


def _counted(listing: str, skills: list[Skill]) -> None:
    """Listed names + announced remainder == total, always."""
    import re
    listed = [p.split(" — ", 1)[0].strip()
              for p in listing.split("; ") if " — " in p or (
                  " more" not in p and " more —" not in p)]
    names = {s.name for s in skills}
    listed = [p for p in listed if p in names]
    m = re.search(r"; (\d+) more —", listing)
    announced = int(m.group(1)) if m else 0
    assert len(listed) + announced == len(skills), (
        f"listed {len(listed)} + announced {announced} "
        f"!= {len(skills)}")
    assert len(listing) <= 2_000


def test_the_count_is_exact_at_41():
    """The first size past a silent cut must still account for everyone."""
    skills = _fake_skills(41)
    _counted(A._skill_listing(skills), skills)


def test_the_count_is_exact_at_200():
    skills = _fake_skills(200)
    _counted(A._skill_listing(skills), skills)


def test_the_count_is_exact_when_a_long_last_entry_does_not_fit():
    """A last entry that overflows the ceiling is announced, never lost."""
    skills = _fake_skills(28)
    # The last one has a name long enough to blow the ceiling by itself.
    skills[-1] = Skill(
        name="x" * 1_900, description="tail", body="b",
        source=pathlib.Path("/nonexistent/last.md"))
    _counted(A._skill_listing(skills), skills)


def test_the_count_is_exact_when_only_the_last_entry_misses_the_ceiling():
    """Exactly the operator's edge: the first 26 fit, the 27th (last)
    does not — announced, not silently dropped."""
    skills = [
        Skill(name=f"s{i:02d}", description="d", body="b",
              source=pathlib.Path(f"/nonexistent/{i}.md"))
        for i in range(26)
    ]
    skills.append(Skill(
        name="y" * 1_950, description="tail", body="b",
        source=pathlib.Path("/nonexistent/last.md")))
    listing = A._skill_listing(skills)
    _counted(listing, skills)
    # The one skill past the ceiling is named by the notice, not silent.
    assert "1 more" in listing
    assert "s25" in listing  # the last fitting entry is still there


def test_the_a_missing_skill_answer_lists_the_whole_catalogue(tmp_path):
    """The overflow notice's promise: a miss returns every skill name."""
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    for i in range(45):
        (ws / ".delfin" / "skills" / f"many-{i:02d}.md").write_text(
            f"# Many {i}\n> Body.\n", encoding="utf-8")
    perms = A.KitToolPermissions(workspace=ws, agent_role="")
    client = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(client._execute_skill({"name": "no-such-skill"}, perms))
    assert "not found" in out.get("error", "")
    names = out.get("available", [])
    assert len(names) >= 45  # every workspace skill, not the first 40


# ---------------------------------------------------------------------------
# The wiring: the listing reaches the tool description capped
# ---------------------------------------------------------------------------

def _tool_desc_for(n: int) -> str:
    """The `skill` tool's description as the surface would build it, with
    the session's skills replaced by *n* fakes."""
    base = A._DOC_TOOLS_OPENAI
    skill_tool = next(
        t for t in base if t.get("function", {}).get("name") == "skill")
    desc = skill_tool["function"]["description"]
    listing = A._skill_listing(_fake_skills(n))
    return desc + f"\nAvailable skills: {listing}"


def test_the_pasted_description_stays_within_the_schema_budget():
    # The house budget test pins the whole catalogue; this pins the paste:
    # with 200 skills on disk the `skill` tool's own description grows by
    # the ceiling, not by the catalogue.
    plain = len(_tool_desc_for(0).split("\nAvailable skills:")[0])
    at_200 = len(_tool_desc_for(200))
    assert at_200 - plain <= 2_000 + len("\nAvailable skills: ") + 50
    assert A.estimate_schema_tokens(
        {"description": _tool_desc_for(200)}) < 1_300


# ---------------------------------------------------------------------------
# The full text still loads only through the `skill` call
# ---------------------------------------------------------------------------

def test_the_listing_never_contains_a_body():
    listing = A._skill_listing(_fake_skills(5))
    assert "body" not in listing  # bodies load via the tool, not the schema


def test_a_real_invocation_still_returns_the_full_body(tmp_path):
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "full-text.md").write_text(
        "---\nname: full-text\ndescription: short\n---\n# Full text\n"
        "> The entire playbook body, many lines.\n", encoding="utf-8")
    perms = A.KitToolPermissions(workspace=ws, agent_role="")
    client = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(client._execute_skill({"name": "full-text"}, perms))
    assert out.get("status") == "ok"
    assert "entire playbook body" in out.get("content", "")
