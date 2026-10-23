"""Before removing what looks wrong, read why it is there.

The other half of the control rule. That one says: prove the defect
before you fix it. This one says: find out whether it IS a defect.

Both halves were used within an hour of each other on 2026-09-19, and
they pointed opposite ways.

  * 1071 directories named ``delfin-gfnff-topo-*`` in /tmp. Nothing said
    why. Driving the code that made them showed the folder outliving
    every session, and the leak was real.

  * 2152 files named ``kit_bg_*`` in /tmp, the same shape of evidence.
    The module's own docstring said they are kept ON PURPOSE -- the
    output of a killed job is how one finds out why it was killed -- and
    named the single cleanup point, a seven-day prune. Measured: nothing
    older than seven days. Removing them would have taken from somebody
    the one record that explains a job they were investigating.

The difference was not in the number, which was similar. It was that one
of them had its reason written down and the other did not. An agent that
tidies what it has not understood is the most expensive kind of helpful.

So the anchor says it, and this file pins the budget it cost -- because a
ratchet is only a ratchet while somebody checks it.
"""

from __future__ import annotations

import pytest

from delfin.agent import prompt_loader as PL


WRITE_RULES = PL._CRITICAL_RULES["write"]


def test_the_anchor_says_to_read_the_reason_first():
    joined = " ".join(WRITE_RULES).lower()
    assert "why it is there" in joined or "read why" in joined, (
        "nothing tells the agent to find out whether a thing is a defect "
        "before removing it:\n  " + "\n  ".join(WRITE_RULES))


def test_the_control_rule_is_still_there():
    """The two halves are one discipline; neither may push out the other."""
    joined = " ".join(WRITE_RULES).lower()
    assert "unfixed" in joined
    assert "covers the edit" in joined


def test_the_budget_was_paid_and_the_reason_recorded():
    """A ratchet is raised by the MEASURED cost, with the reason — and
    then it holds at the new number.

    Before this rule: 5 rules, 517 characters. After it: 6 and 642. The
    rule is one sentence, the anchor is where recency puts it in front of
    the model, and there was no weaker line left to fold it into — so the
    ceiling moves by what it actually costs. 642, not 640 and not 700:
    the first draft of this test guessed 640 and was two characters
    short, which is the whole argument for measuring.
    """
    assert len(WRITE_RULES) <= 6, f"{len(WRITE_RULES)} rules"
    size = sum(len(r) for r in WRITE_RULES)
    assert size <= 642, f"{size} characters in the write anchor"
    assert size > 517, (
        "if the anchor did not grow, this ceiling was raised for nothing "
        "and should come back down")


@pytest.mark.parametrize("role", ["write", "review", "plan"])
def test_every_role_still_has_its_anchor(role):
    rules = PL._CRITICAL_RULES[role]
    assert rules and all(isinstance(r, str) and r.strip() for r in rules)
