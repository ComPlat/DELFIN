"""The plan-before-act check could not tell "before" from "after".

The task is named for the order: say what you are about to do, then do
it. Its expectation was a line starting with ``1.`` / ``step 1`` /
``schritt 1`` / ``plan`` / ``todo``, matched anywhere in the answer.

Two things follow, and both were happening in the same five samples.

It missed the plan the agent actually writes. The prompt asks for "a
brief 1-line plan, then ACTIONs in order, then /done", the agent obeys
with "Ich führe alle 5 Schritte aus:", and the pattern scores it wrong
because the line does not begin with a digit. The agent was following
its instruction and being marked down for the wording.

And it accepted a summary written AFTERWARDS. The one sample that passed
did so on "1. ✅ functional auf BP86 gesetzt" -- a numbered list printed
after the actions had already run. That is the opposite of the behaviour
the task exists to measure, and it was the only thing keeping the score
above zero.

So the check rewarded reporting and punished planning, in a task called
plan-before-act. The order is part of the pattern: something has to be
said before the first ACTION line, and nothing said afterwards can
substitute for it.

2026-08-13: the expectation was REMOVED from workflow_plan_before_act,
deliberately — that task's prompt is a numbered list the user wrote, and
an agent that restates it before executing is paying tokens to repeat the
question. See the note in tasks.yaml.

The pattern is kept here, and so are its tests. It is correct, it was
expensive to get right, and plan-before-act still applies where the user
did NOT enumerate the steps. When that task gets written, the pattern is
proven and ready instead of re-derived — and these tests are what proves
it. They no longer read tasks.yaml: a test that follows a file it does
not own turns a deliberate product decision into a red suite, which is
what happened when the expectation came out.
"""

from __future__ import annotations

import re

import pytest


#: Something of substance on a line of its own before the first ACTION.
#: ``\A`` anchors it to the start of the answer, so a summary printed
#: after the actions cannot satisfy it.
PLAN_BEFORE_ACT_PATTERN = r'(?sm)\A\s*(?!ACTION:)\S[^\n]{4,}\n[\s\S]*?^ACTION:'


def _matches(answer: str) -> bool:
    return bool(re.search(PLAN_BEFORE_ACT_PATTERN, answer))


_ACTIONS = (
    "ACTION: /control key functional BP86\n"
    "ACTION: /control key main_basisset def2-TZVP\n"
    "ACTION: /control key disp_corr D3BJ\n"
    "ACTION: /tab submit\n"
    "ACTION: /done\n"
)


# ---------------------------------------------------------------------------
# A plan said before acting counts, whatever its wording
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("opening", [
    "Ich führe alle 5 Schritte aus:",
    "Ich führe alle 5 Schritte jetzt korrekt aus:",
    "I'll set the functional, the basis and the dispersion, then submit.",
    "1. functional, 2. basis, 3. dispersion, 4. submit",
    "Plan: functional, basis, dispersion, submit.",
])
def test_an_announcement_before_the_actions_counts(opening):
    assert _matches(f"{opening}\n\n{_ACTIONS}"), opening


def test_the_german_line_the_agent_actually_writes_counts():
    """Observed in 4 of 5 samples, and scored wrong in every one."""
    assert _matches("Ich führe alle 5 Schritte aus:\n\n" + _ACTIONS)


# ---------------------------------------------------------------------------
# ...and acting without a word does not
# ---------------------------------------------------------------------------

def test_actions_with_nothing_before_them_fail():
    assert not _matches(_ACTIONS)


def test_a_summary_written_afterwards_does_not_rescue_it():
    """The one passing sample passed on exactly this, which is the
    behaviour the task exists to catch."""
    after = (_ACTIONS + "\nAlle 5 Schritte wurden ausgeführt:\n"
             "1. ✅ **functional** auf `BP86` gesetzt\n"
             "2. ✅ **basis** auf `def2-TZVP` gesetzt\n")
    assert not _matches(after)


def test_a_numbered_summary_after_a_real_plan_is_still_fine():
    """Reporting afterwards is good practice; it just cannot stand in for
    the plan."""
    both = ("Ich führe alle 5 Schritte aus:\n\n" + _ACTIONS
            + "\n1. ✅ functional gesetzt\n")
    assert _matches(both)


def test_an_empty_answer_fails():
    assert not _matches("")


def test_a_blank_line_is_not_a_plan():
    assert not _matches("\n   \n" + _ACTIONS)


def test_a_single_word_is_not_a_plan():
    """Guard against the pattern degenerating into "any character wins"."""
    assert not _matches("Ok\n" + _ACTIONS)
