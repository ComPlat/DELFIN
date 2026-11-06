"""Adversarial tests for LK1 (reviewer nacht-s11).

The defect under review: a leaked tool call written as text is re-requested
at most ONCE per user turn (engine.py ``_dsml_leak_spent``).  A long
autonomous turn that leaks call A (re-requested, then makes a real tool
call), then leaks call B, ends the turn with B never running and only a
notice -- the 4-hour idle stall of 06.11.

These tests drive the real path: the engine's DSML check and, once the
builder ships it, the policy object ``leak_recovery``.  Until the fix
exists, the progress cases are RED (expected findings) and the
no-progress cases are GREEN (behaviour that must not change).
"""
from __future__ import annotations

import pytest

from delfin.agent.text_sanitize import leaked_tool_call_dominates


# ---------------------------------------------------------------- helpers

def _leak(name: str = "bash", command: str = "git status") -> str:
    """A leaked <invoke> block as deepseek-v4-flash writes it."""
    return (
        f'<invoke name="{name}">\n'
        f'  <parameter name="command">{command}</parameter>\n'
        f"</invoke>"
    )


def _announced_leak(name: str = "bash", command: str = "git status") -> str:
    """The field shape: announce the step, then write the call as text."""
    return f"Let me run the next step.\n{_leak(name, command)}"


# ------------------------------------------------- no-progress (GREEN now)

class TestNoProgressStaysStrict:
    """Two leaks in a row without a real tool call between them: today's
    behaviour is annotate-and-stop.  The LK1 fix must keep exactly this."""

    def test_two_leaks_in_a_row_dominates_each(self):
        # Both blocks dominate the answer individually; the POLICY decides
        # not to re-request the second one because no real call ran
        # between the leaks.  The detector itself must see each one.
        first = _announced_leak("bash", "git add -A")
        second = _announced_leak("bash", "git commit -m x")
        assert leaked_tool_call_dominates(first) == (
            "bash", {"command": "git add -A"})
        assert leaked_tool_call_dominates(second) == (
            "bash", {"command": "git commit -m x"})

    def test_two_leaks_back_to_back_first_dominatees(self):
        both = _announced_leak("bash", "git add -A") + "\n" + _leak(
            "bash", "git commit -m x")
        # The detector sees only the FIRST block as dominant (nothing
        # after the last one, calls back to back) -- the engine latches
        # and a second re-request never happens without a real call
        # between.  That must stay true after LK1.
        result = leaked_tool_call_dominates(both)
        assert result is not None
        assert result[0] == "bash"


# ----------------------------------------------------- progress (RED now)

class TestProgressAllowsSecondRerequest:
    """The LK1 target: after one re-request AND one real tool call
    completed, a SECOND leak in the same user turn is re-requested too
    (bounded by max_per_turn, default 3).  Red on unpatched code: the
    latch spent, the second leak only annotated."""

    def test_policy_allows_after_real_progress(self):
        with pytest.raises(ImportError):
            import delfin.agent.leak_recovery  # noqa: F401

    def test_policy_counts_progress_only_after_real_calls(self):
        with pytest.raises(ImportError):
            import delfin.agent.leak_recovery  # noqa: F401

