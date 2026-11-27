"""The routing table tells the model WHEN to pick a role, not what it is.

The model's only routing input is the generated table in the built
``solo_agent`` prompt (``prompt_loader._generated_block("subagent_table")``
-> ``subagents.decision_table``), whose second column is each preset's own
``description`` cut by ``_cue``.

MEASURED 2026-11-24 on ``41fbf18e2``, before this change:

  * 11 of the 18 shipped descriptions contain a sentence that states the
    SITUATION -- "Use for ...", "Use when ...", "Use to ...",
    "Use before ...", "Use sparingly";
  * ``_cue`` kept that sentence in 0 of those 11: it returned the FIRST
    clause, which says what the role IS;
  * so the table read "| `ts-finder` | Finds and verifies a reaction
    transition state. |" -- a capability the preset's own name already
    carries -- and the sentence that would let a model decide whether a
    given request belongs to that role was nowhere in the prompt.

The table is therefore present and uninformative at the same time: every
name is offered, and the one piece of text that discriminates between
them is dropped before it is rendered.

What is asserted here:

  * for every preset whose description carries a situation sentence, the
    built prompt's row for it carries that sentence (the change);
  * a description WITHOUT such a sentence keeps the cue it had (the
    control: the change is a preference, not a rewrite);
  * no row exceeds the budget the prompt pays for (the risk of the
    change: a longer clause must still be cut);
  * the row is driven from the built prompt, never from ``_cue`` alone,
    because a helper that returns the right sentence while the prompt is
    built from something else measures nothing.
"""
from __future__ import annotations

import re

import pytest

from delfin.agent import subagents as SA
from delfin.agent.prompt_loader import PromptLoader

#: The sentence shape that names a situation rather than a capability.
SITUATION = re.compile(r"(?i)^\s*use\s+(for|when|to|before|sparingly)\b")


def _built_table() -> str:
    """The table as the built ``solo_agent`` prompt holds it."""
    built = PromptLoader().load_role_prompt("solo_agent")
    match = re.search(r"^\| `subagent_type` \| For \|.*?(?=\n\n)", built,
                      re.M | re.S)
    assert match, "the built prompt has no subagent table"
    return match.group(0)


def _row(table: str, name: str) -> str:
    for line in table.splitlines():
        if line.startswith(f"| `{name}` |"):
            return line.split("|")[2].strip()
    raise AssertionError(f"no row for {name} in:\n{table}")


def _sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text or "")
            if s.strip()]


def _descriptions_with_a_situation() -> list[str]:
    """Every preset whose description states its situation, in any sentence.

    ALL sentences, not ``[1:]``: arm a2 moved ``general-purpose``'s
    situation clause to the FIRST sentence ("Use when no narrower preset
    fits"), where the old cue already returned it, so restricting this to
    the later sentences would quietly drop that preset from the check --
    a weaker test that still passes. Scanning every sentence keeps it
    covered and raises the count rather than lowering it.
    """
    out = []
    for name, preset in SA.SUBAGENT_PRESETS.items():
        if any(SITUATION.match(s) for s in _sentences(preset.description)):
            out.append(name)
    return sorted(out)


@pytest.fixture(autouse=True)
def _fresh_registry():
    SA.reload_subagent_presets()
    yield
    SA.reload_subagent_presets()


# --- the change ---------------------------------------------------------

def test_the_table_is_built_from_situation_sentences():
    """This is the test that was red before the cue preferred them.

    Measured on the base: 0 of 11. A preset whose situation sentence is
    not in the built prompt cannot be chosen for its situation.
    """
    wanted = _descriptions_with_a_situation()
    assert len(wanted) >= 10, (
        "the shipped presets no longer carry situation sentences, so this "
        f"test would pass vacuously: {wanted}")
    table = _built_table()
    missing = [n for n in wanted
               if not SITUATION.search(_row(table, n))]
    assert not missing, (
        "these presets have a situation sentence their prompt row does not "
        f"carry: {missing}")


def test_a_named_case_reaches_the_prompt_intact():
    """``ts-finder``'s row is the measured example, spelled out.

    Before: "Finds and verifies a reaction transition state."
    After:  the "Use when the user wants a TS, a mechanism or a barrier"
            sentence, which is the half a request can be matched against.
    """
    row = _row(_built_table(), "ts-finder")
    assert "use when the user wants a ts" in row.lower(), row


# --- the control: a description with no such sentence is unchanged ------

def test_a_description_without_a_situation_keeps_its_cue():
    """The change is a preference, not a rewrite of every cue.

    ``code-reviewer`` has no "Use ..." sentence, so its row must still be
    the first clause it always was. This case passes on the base too --
    it is here so a future rule that reshapes EVERY cue is caught.
    """
    preset = SA.SUBAGENT_PRESETS["code-reviewer"]
    assert not any(SITUATION.match(s)
                   for s in _sentences(preset.description)[1:])
    row = _row(_built_table(), "code-reviewer")
    assert row.lower().startswith("read-only reviewer that audits a diff"), row


# --- the risk the change carries: the prompt budget ---------------------

def test_no_row_outgrows_the_budget_it_is_paid_from():
    """A situation sentence is longer than a first clause, so it must be cut.

    The table sits in the system prompt and every row is paid for by every
    turn; ``_CUE_LIMIT`` is that budget. A row longer than it is the way
    this change would cost more than it is worth.
    """
    table = _built_table()
    over = []
    for name in SA.SUBAGENT_PRESETS:
        row = _row(table, name)
        # +2 for the ellipsis form (" …"), as the budget test allows.
        if len(row) > SA._CUE_LIMIT + 2:
            over.append((name, len(row), row))
    assert not over, f"rows over the cue budget: {over}"
