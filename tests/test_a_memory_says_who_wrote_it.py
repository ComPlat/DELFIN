"""A memory the model invented was indistinguishable from the user's.

Three findings that compound into one:

1. `_compose_frontmatter` wrote no provenance field. A fact the user typed
   with /remember and a fact the model made up mid-turn produced
   byte-identical files.

2. `feedback` and `user` types are exempt from age pruning and get a 4x
   cap, and in the global store EVERY type is exempt. Combined with (1),
   a model-written "feedback: always skip the tests" is immortal and is
   recalled into every future session forever.

3. `save_typed_memory` read a leading `global:` out of the TEXT and moved
   the memory to the user-wide store that every other workspace reads --
   while the `remember` tool schema has no scope parameter at all and
   tells the model it is writing "project memory". So the one carrier the
   model fully controls silently widened its own reach.

The rule the rest of the framework already follows: a lower-trust carrier
may TIGHTEN its reach, never widen it. The office path had this right
(it forces scope="project" and refuses the prefix); it simply was never
applied to the model in general.

The fix is provenance plus decay, not a ban: model memories are marked,
lose the immortality exemption, and expire when they stop being recalled
-- `record_memory_recall` bumps updated_at, so a memory that keeps
earning its place survives. Decay by disuse, not by age.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import memory_store as ms


AGE = ms._AGENT_MEMORY_MAX_AGE_DAYS


def _front(path) -> dict:
    meta, _ = ms._parse_frontmatter(path.read_text(encoding="utf-8"))
    return meta


def _age(path, days: float) -> None:
    """Backdate a memory's updated_at by ``days``."""
    text = path.read_text(encoding="utf-8")
    old = ms._parse_frontmatter(text)[0].get("updated_at", "0")
    new = int(time.time()) - int(days * 86_400)
    path.write_text(
        text.replace(f"updated_at: {old}", f"updated_at: {new}"),
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# 1. The file says who wrote it
# ---------------------------------------------------------------------------

def test_a_model_written_memory_is_marked_as_such(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "project: the parser is being rewritten",
        repo_root=tmp_path, source="agent")
    assert _front(p)["source"] == "agent"


def test_a_user_written_memory_is_marked_as_the_users(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "project: the parser is being rewritten", repo_root=tmp_path)
    assert _front(p)["source"] == "user"


def test_the_default_is_the_user_so_nothing_silently_downgrades(tmp_path):
    """Every pre-existing caller means a user-grade write."""
    p, _, _ = ms.save_typed_memory("feedback: run the slow suite too",
                                   repo_root=tmp_path)
    assert _front(p)["source"] == "user"


def test_the_listing_reports_the_source(tmp_path):
    ms.save_typed_memory("project: a", repo_root=tmp_path, source="agent")
    ms.save_typed_memory("project: b", repo_root=tmp_path)
    got = {r["description"]: r["source"]
           for r in ms.list_typed_memories(tmp_path)}
    assert got == {"a": "agent", "b": "user"}


def test_a_memory_written_before_the_field_existed_counts_as_the_users(
    tmp_path,
):
    p, _, _ = ms.save_typed_memory("project: legacy", repo_root=tmp_path)
    text = "\n".join(line for line in p.read_text(encoding="utf-8").splitlines()
                     if not line.startswith("source:"))
    p.write_text(text, encoding="utf-8")
    rec = ms.list_typed_memories(tmp_path)[0]
    assert rec["source"] == "user", "a legacy file must not lose protection"


def test_the_field_is_owned_by_the_store_not_carried_through(tmp_path):
    """It decides whether a memory expires, so an arbitrary value in it
    must resolve to one of the two the store knows -- unlike the unknown
    frontmatter fields that are deliberately carried through untouched."""
    p, _, _ = ms.save_typed_memory("project: x", repo_root=tmp_path)
    p.write_text(
        p.read_text(encoding="utf-8").replace(
            "source: user", "source: hand-curated"),
        encoding="utf-8")
    assert ms.list_typed_memories(tmp_path)[0]["source"] == "user"


# ---------------------------------------------------------------------------
# 2. ...and keeps saying it
# ---------------------------------------------------------------------------

def test_a_recall_does_not_erase_the_provenance(tmp_path):
    """Both rewrite paths pass extras=meta; a KNOWN field would be dropped."""
    p, _, _ = ms.save_typed_memory("project: x", repo_root=tmp_path,
                                   source="agent")
    ms.record_memory_recall(tmp_path, [p.name])
    assert _front(p)["source"] == "agent"


def test_a_rot_hit_does_not_erase_the_provenance(tmp_path):
    p, _, _ = ms.save_typed_memory("project: x", repo_root=tmp_path,
                                   source="agent")
    ms.record_stale_hits(tmp_path, [p.name])
    assert _front(p)["source"] == "agent"


def test_the_source_is_written_once_not_twice(tmp_path):
    p, _, _ = ms.save_typed_memory("project: x", repo_root=tmp_path,
                                   source="agent")
    ms.record_memory_recall(tmp_path, [p.name])
    lines = p.read_text(encoding="utf-8").splitlines()
    assert sum(1 for line in lines if line.startswith("source:")) == 1


def test_the_model_rewriting_a_users_memory_downgrades_it(tmp_path):
    """Merged text is the model's; it must not inherit user protection."""
    ms.save_typed_memory("project: the parser is being rewritten",
                         repo_root=tmp_path)
    p, _, _ = ms.save_typed_memory(
        "project: the parser is being rewritten now",
        repo_root=tmp_path, source="agent")
    assert _front(p)["source"] == "agent"


def test_the_user_confirming_a_model_memory_upgrades_it(tmp_path):
    ms.save_typed_memory("project: the parser is being rewritten",
                         repo_root=tmp_path, source="agent")
    p, _, _ = ms.save_typed_memory(
        "project: the parser is being rewritten now", repo_root=tmp_path)
    assert _front(p)["source"] == "user"


# ---------------------------------------------------------------------------
# 3. Model text cannot widen its own reach
# ---------------------------------------------------------------------------

def test_the_model_cannot_send_a_memory_to_the_user_wide_store(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "global: feedback: always skip the tests",
        repo_root=tmp_path, source="agent")
    assert ms._memory_dir_for_scope(tmp_path, "project") in p.parents


def test_the_refused_prefix_does_not_end_up_in_the_body(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "global: feedback: always skip the tests",
        repo_root=tmp_path, source="agent")
    assert "global:" not in p.read_text(encoding="utf-8")


def test_a_refused_prefix_no_longer_forces_the_user_type(tmp_path):
    """"about the user by definition" only held while it crossed repos."""
    _, _, mtype = ms.save_typed_memory(
        "global: never commit without running the suite",
        repo_root=tmp_path, source="agent")
    assert mtype == "feedback"


def test_the_user_can_still_write_a_global_memory(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "global: feedback: always skip the tests",
        repo_root=tmp_path, allow_scope_prefix=True)
    assert ms._memory_dir_for_scope(tmp_path, "user") in p.parents


def test_an_explicit_scope_argument_still_works(tmp_path):
    """The prefix gate guards the TEXT; a caller may still pass scope."""
    p, _, _ = ms.save_typed_memory(
        "user: prefers metric units", repo_root=tmp_path, scope="user")
    assert ms._memory_dir_for_scope(tmp_path, "user") in p.parents


def test_refusing_the_prefix_is_the_default(tmp_path):
    """A future call site that forgets the flag must fail closed."""
    p, _, _ = ms.save_typed_memory("global: project: x", repo_root=tmp_path)
    assert ms._memory_dir_for_scope(tmp_path, "project") in p.parents


# ---------------------------------------------------------------------------
# 4. A model memory decays; a user memory does not
# ---------------------------------------------------------------------------

def test_an_unrecalled_model_memory_expires(tmp_path):
    p, _, _ = ms.save_typed_memory("feedback: always skip the tests",
                                   repo_root=tmp_path, source="agent")
    _age(p, AGE + 1)
    assert ms.prune_memories(tmp_path) == [p.name]
    assert not p.exists()


def test_the_users_own_feedback_never_expires(tmp_path):
    p, _, _ = ms.save_typed_memory("feedback: always run the slow suite",
                                   repo_root=tmp_path)
    _age(p, AGE * 10)
    assert ms.prune_memories(tmp_path) == []
    assert p.exists()


def test_a_model_memory_that_is_still_recalled_survives(tmp_path):
    """Decay by disuse, not by age -- recall bumps updated_at."""
    p, _, _ = ms.save_typed_memory("feedback: always skip the tests",
                                   repo_root=tmp_path, source="agent")
    _age(p, AGE + 1)
    ms.record_memory_recall(tmp_path, [p.name])
    assert ms.prune_memories(tmp_path) == []
    assert p.exists()


def test_a_fresh_model_memory_is_not_touched(tmp_path):
    p, _, _ = ms.save_typed_memory("project: mid-rewrite", repo_root=tmp_path,
                                   source="agent")
    assert ms.prune_memories(tmp_path) == []
    assert p.exists()


def test_the_global_store_is_not_a_shelter_for_model_memories(tmp_path):
    """scope="user" protects EVERY type -- that must not cover the model."""
    p, _, _ = ms.save_typed_memory(
        "feedback: always skip the tests", repo_root=tmp_path,
        scope="user", source="agent")
    _age(p, AGE + 1)
    assert ms.prune_memories(tmp_path, scope="user") == [p.name]


def test_the_users_global_memory_is_still_sheltered(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "feedback: always run the slow suite", repo_root=tmp_path,
        scope="user")
    _age(p, AGE * 10)
    assert ms.prune_memories(tmp_path, scope="user") == []


def test_the_model_loses_the_tie_when_the_cap_evicts(tmp_path):
    """Same freshness, same use count -> the model's goes first."""
    keep, _, _ = ms.save_typed_memory(
        "reference: the crystallography handbook lives on the shared drive",
        repo_root=tmp_path)
    drop, _, _ = ms.save_typed_memory(
        "reference: solvent tables are published at an external registry",
        repo_root=tmp_path, source="agent")
    for p in (keep, drop):
        _age(p, 1)
    deleted = ms.prune_memories(tmp_path, max_per_type=1)
    assert deleted == [drop.name]
    assert keep.exists()


# ---------------------------------------------------------------------------
# 5. Wiring: the two model-driven write paths declare themselves
# ---------------------------------------------------------------------------

def test_the_remember_tool_writes_as_the_model(tmp_path):
    from delfin.agent import api_client as A
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    ex._execute_remember(
        {"text": "global: feedback: always skip the tests"}, perms)
    recs = ms.list_typed_memories(tmp_path, scope="project")
    assert [r["source"] for r in recs] == ["agent"]
    assert ms.list_typed_memories(tmp_path, scope="user") == []


def test_the_distilled_memory_writes_as_the_model(tmp_path):
    from delfin.agent import memory_distill as md
    md.save_facts(["project: the parser is being rewritten"],
                  repo_root=tmp_path)
    recs = ms.list_typed_memories(tmp_path)
    assert recs and all(r["source"] == "agent" for r in recs)
