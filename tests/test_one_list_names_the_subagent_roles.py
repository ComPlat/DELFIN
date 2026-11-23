"""The roles the model may call and the roles it is told about are one list.

Input: the preset registry, and the prompt as the loader BUILDS it (not
the file on disk -- a rule is a rule when it is in the built prompt).

Output measured 2026-11-22, before this change: the ``subagent_type``
enum offered 13 specialist roles and the hand-maintained table in
``pack/agents/solo_agent.md`` named 7 of them. ``irc-checker``,
``reactivity-department``, ``metal-smiles-builder``, ``swarm-lead``,
``ts-finder`` and ``xyz-editor`` were valid enum values the dispatch
accepts, with no "when to use" cue anywhere in the prompt. In 79,056
tool calls the subagent tool was dispatched 50 times and 44 of those
were ``explore`` -- the first row of the table, with the broadest cue.
A role the model is never told the purpose of is not an offer.

Three consumers now read one dict: the tool-schema enum
(``subagent_type_names``), the dispatch (``SUBAGENT_PRESETS``
membership) and the decision table (``decision_table``). This module
pins that they cannot drift apart again.
"""
from __future__ import annotations

import re

import pytest

from delfin.agent import subagents
from delfin.agent.prompt_loader import PromptLoader

GENERATED = re.compile(
    r"^<!--\s*generated:([a-z_]+)\s*-->.*?^<!--\s*/generated\s*-->",
    re.MULTILINE | re.DOTALL)
TABLE_HEADER = "| `subagent_type` | For |"


def _built() -> str:
    return PromptLoader().load_role_prompt("solo_agent")


# --- the drift that was measured ---------------------------------------

def test_every_offered_role_is_named_in_the_built_prompt():
    """Every name the enum offers appears in the prompt the model reads.

    This is the test that was red before the table was generated: six
    specialists were callable and unnamed.
    """
    built = _built()
    missing = [n for n in subagents.subagent_type_names()
               if f"`{n}`" not in built]
    assert not missing, (
        "offered by the tool schema and not named in the built prompt: "
        f"{sorted(missing)}")


def test_every_offered_role_carries_a_cue_not_just_a_name():
    """A bare name is not an offer: each row says what the role is for."""
    built = _built()
    for name in subagents.subagent_type_names():
        row = next((ln for ln in built.splitlines()
                    if ln.startswith(f"| `{name}` |")), None)
        assert row is not None, f"{name} has no table row"
        cue = row.split("|")[2].strip()
        assert len(cue) >= 12, f"{name} has no usable cue: {row!r}"


def test_the_table_and_the_enum_are_the_same_set():
    """The two lists are one list, in both directions."""
    table = {m.group(1) for m in re.finditer(
        r"^\| `([a-z0-9-]+)` \|", subagents.decision_table(), re.MULTILINE)}
    assert table == set(subagents.subagent_type_names())


def test_the_dispatch_accepts_every_name_the_table_names():
    """Membership of the dict the dispatch validates against is what the
    table is rendered from, so this holds by construction -- and the
    test says so, because that is the property that broke."""
    for name in re.findall(r"^\| `([a-z0-9-]+)` \|",
                           subagents.decision_table(), re.MULTILINE):
        assert name in subagents.SUBAGENT_PRESETS


# --- the drift cannot come back ----------------------------------------

def test_a_preset_added_at_runtime_reaches_the_built_prompt(monkeypatch):
    """A preset registered after import is offered AND presented. The
    enum already followed the registry; the table now does too."""
    preset = subagents.SubagentPreset(
        name="late-arrival",
        description="Read-only probe added during the session. Second "
                    "sentence that must not reach the table.",
        system_prompt="probe")
    monkeypatch.setitem(subagents.SUBAGENT_PRESETS, "late-arrival", preset)
    PromptLoader._cached_read.cache_clear() if hasattr(
        PromptLoader._cached_read, "cache_clear") else None
    built = PromptLoader().load_role_prompt("solo_agent")
    assert "| `late-arrival` |" in built
    assert "Read-only probe added during the session." in built
    assert "must not reach the table" not in built, (
        "the cue is the first clause; the rest belongs to /subagents")


def test_no_pack_file_keeps_a_hand_written_copy_of_the_table():
    """A second copy is how this broke. Any table header inside the pack
    must sit inside a generated block."""
    pack = PromptLoader().agent_dir / "agents"
    for path in sorted(pack.glob("*.md")):
        text = path.read_text(errors="replace")
        if TABLE_HEADER not in text:
            continue
        outside = GENERATED.sub("", text)
        assert TABLE_HEADER not in outside, (
            f"{path.name} holds a subagent table outside a generated "
            "block: that is a hand-maintained copy and it will drift")


def test_the_dashboard_role_deliberately_has_no_table():
    """The dashboard operator explains, researches documents and runs UI
    actions; delegation is not one of its jobs, so it gets no table.

    Pinned as a decision rather than left as an accident: if the
    dashboard is ever meant to delegate, it gets the generated block,
    never a second hand-written copy.
    """
    built = PromptLoader().load_role_prompt("dashboard_agent")
    assert TABLE_HEADER not in built


# --- the prompt budget --------------------------------------------------

def test_a_cue_stays_within_the_budget_it_is_paid_from():
    """The table sits in the system prompt, so every row is paid for by
    every turn. A long description is cut at its first clause."""
    for name in subagents.subagent_type_names():
        cue = subagents._cue(
            getattr(subagents.SUBAGENT_PRESETS[name], "description", ""))
        assert len(cue) <= subagents._CUE_LIMIT + 2, (name, len(cue), cue)


@pytest.mark.parametrize("description,expected", [
    ("Read-only research agent for locating code.",
     "Read-only research agent for locating code."),
    ("First clause. Second sentence.", "First clause."),
    ("Routes a question to the right expert — the ts-finder for a TS.",
     "Routes a question to the right expert."),
    ("One job; another job.", "One job."),
    ("", ""),
])
def test_the_cue_is_the_first_clause_of_the_description(description, expected):
    """One source for the cue: the description the preset already has.
    No second frontmatter field, because two fields for one purpose is
    how the table and the enum drifted apart."""
    assert subagents._cue(description) == expected


# --- the mechanism itself ----------------------------------------------

def test_a_generated_name_nobody_owns_leaves_the_block_alone():
    """A prompt that keeps a stale section is better than a prompt with a
    hole where a section was."""
    loader = PromptLoader()
    text = ("before\n<!-- generated:not_a_real_block -->\nkept\n"
            "<!-- /generated -->\nafter\n")
    assert loader._fill_generated_blocks(text) == text


def test_a_file_without_markers_is_returned_unchanged():
    loader = PromptLoader()
    text = "# A prompt\n\nNo markers here.\n"
    assert loader._fill_generated_blocks(text) is text
