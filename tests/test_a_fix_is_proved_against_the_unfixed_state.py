"""The agent is told to prove a fix, not just to run tests after it.

The critical anchor already said "after a code edit, run the tests that
cover it". That catches a fix which BREAKS something. It says nothing
about the other half, which is the one that matters more: a test written
alongside a fix and never run against the unfixed code proves nothing at
all. It passes because the fix is there, and it would have passed without
the defect ever having existed.

Every real finding in this codebase came from the other order. The undo
journal leaking another run's paths, the shell variable that walked out of
the workspace, the quoted directory that was granted and unreachable —
each was a claim until the same check was run against the unchanged code
and came back red. Two of them turned out to be the instrument rather
than the product, and only the control said so.

So it belongs where the model actually reads it: the critical anchor at
the END of the prompt, which is repeated for recency. And the budget is a
ratchet, so this pays from its own text — the line it replaces said less
in more words.
"""

from __future__ import annotations

import pytest

from delfin.agent import prompt_loader as PL


WRITE_RULES = PL._CRITICAL_RULES["write"]


def test_the_write_roles_are_told_to_run_the_control_first():
    joined = " ".join(WRITE_RULES).lower()
    assert "unfixed" in joined or "before the fix" in joined or (
        "fails first" in joined or "red first" in joined), (
        "a test never run against the unfixed code proves nothing, and "
        "nothing in the anchor says so:\n  " + "\n  ".join(WRITE_RULES))


def test_it_still_says_to_run_the_tests_afterwards():
    """The new half must not push out the old one: a control proves the
    defect was real, the run afterwards proves nothing else broke."""
    joined = " ".join(WRITE_RULES).lower()
    assert "covers the edit" in joined, (
        "the control proves the defect was real; running what covers the "
        "edit proves nothing else broke. Both halves or neither:\n  "
        + "\n  ".join(WRITE_RULES))
    assert "unrelated failures" in joined, (
        "and the agent must still not wander off fixing what it did not "
        "break")


@pytest.mark.parametrize("role", ["write", "review", "plan"])
def test_every_role_still_has_its_anchor(role):
    rules = PL._CRITICAL_RULES[role]
    assert rules and all(isinstance(r, str) and r.strip() for r in rules)


def test_the_anchor_did_not_grow_without_paying():
    """A ratchet: a new capability pays from its own text. The write
    anchor may not get longer than it was for adding this."""
    # The control rule itself paid: it replaced the weaker half of the
    # test rule rather than being added beside it, leaving 5 rules and
    # 517 characters. The ceiling was raised once afterwards, to 6 and
    # 642, by "read why it is there before removing it" — see
    # test_read_why_before_you_change_it.py, which records what that one
    # cost and why no line was left to fold it into.
    assert len(WRITE_RULES) <= 6, (
        f"{len(WRITE_RULES)} rules — the anchor grew instead of paying")
    assert sum(len(r) for r in WRITE_RULES) <= 642, (
        f"{sum(len(r) for r in WRITE_RULES)} characters in the write anchor")


def test_the_rule_reaches_the_built_prompt():
    """A rule is a rule when it is IN the prompt. shared/ is not
    automatically shared, and this has been found the hard way before."""
    import inspect
    src = inspect.getsource(PL)
    assert "_CRITICAL_RULES" in src
    # The anchor is selected by role category; solo and builder are write.
    assert PL._ROLE_TO_RULE_CATEGORY.get("solo_agent") == "write"
    assert PL._ROLE_TO_RULE_CATEGORY.get("builder_agent") == "write"
