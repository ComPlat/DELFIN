"""T2 phase-3 operator findings (reviewer s15): the context-budget restart.

The operator (review of the rejected engine patch) found three defects in the
phase-3 fresh-start path and one requirement. These tests are the RED controls
for those findings, written against the CURRENT (unfixed) tree. They must fail
before the fix and pass after it.

Defect 2 -- `fresh_context_for_turn` reads the CUMULATIVE session input counter
            (``engine.token_usage['input']``), which never shrinks. Once a
            session is over the budget, EVERY later turn restarts fresh and the
            conversation loses its short-term memory for good. It must use the
            CURRENT context size (the next request), so that a turn N+1 after a
            fresh start -- whose context is now small -- is NOT fresh again.
Defect 3 -- `fresh_context_for_turn` has no caller in repl.py and nothing sets
            ``engine.start_fresh``. The terminal must wire it before each turn.
Requirement -- the fresh start keeps the CURRENT user message (the prompt the
            turn must answer is never dropped): "one user message = fresh block
            + current prompt".
"""
from __future__ import annotations

import pytest

from delfin.agent import task_state
from delfin.agent.repl import fresh_context_for_turn, _FRESH_CONTEXT_MARKER


class _FakeEngine:
    """Fake engine exposing the NEW contract: a current-context token count
    ``estimate_context_tokens()`` plus an ``start_fresh`` slot the terminal
    arms. No cumulative counter is consulted for the trigger."""

    def __init__(self, current_tokens: int):
        self._current_tokens = current_tokens
        self.start_fresh = ""

    def estimate_context_tokens(self) -> int:
        return self._current_tokens


def _state(tmp_path):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task="review phase 3", phase="phase 3", phase_status="in_progress",
              commit="abc123")
    return st


def test_defect2_after_fresh_start_next_turn_is_not_fresh(tmp_path):
    """Turn N consumed the budget (fresh fired); turn N+1 has a SMALL current
    context now. Because the trigger reads the CURRENT context size and not the
    cumulative counter, turn N+1 must NOT be fresh again. RED on the unfixed
    hook, which triggers on the cumulative counter forever."""
    engine = _FakeEngine(current_tokens=2_000)     # small again after fresh
    engine.token_usage = {"input": 950_000}        # ...but cumulative is huge
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == "", (
        "turn N+1 after a fresh start must keep the (now small) context; "
        "triggering on the cumulative counter would restart fresh forever")


def test_defect2_metric_is_current_context_not_cumulative(tmp_path):
    """The over-budget trigger is the CURRENT context estimate, not the running
    session counter. A session whose cumulative usage is huge but whose current
    context is small must NOT restart fresh (the decimal cumulative counter is
    the defect). RED on the unfixed hook."""
    engine = _FakeEngine(current_tokens=100_000)   # current context, under budget
    engine.token_usage = {"input": 2_500_000}     # cumulative, way over
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""


def test_defect3_terminal_arms_start_fresh_when_over_budget():
    """The terminal (repl) must call the hook and set ``engine.start_fresh``
    before a turn. Regression guard: the public hook and the engine slot both
    exist and the hook returns a block when the current context is over the
    budget (so the terminal has something to arm with)."""
    from delfin.agent.repl import fresh_context_for_turn
    engine = _FakeEngine(current_tokens=950_000)
    ctx = fresh_context_for_turn(engine, None, token_budget=900_000)
    assert ctx and _FRESH_CONTEXT_MARKER in ctx
