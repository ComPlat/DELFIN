"""T2 phase-3 operator findings (reviewer s15): the context-budget restart.

The operator (review of the rejected engine patch) found three defects in the
phase-3 fresh-start path and one requirement. These tests are the RED/GREEN
controls for those findings.

Defect 2 -- `fresh_context_for_turn` must judge the budget on the CURRENT
             context size (``engine._estimate_context_tokens`` -- the messages
             about to be sent), not the session's cumulative input-token
             counter, which never shrinks. Once a session passes the budget
             the cumulative counter stays over it, so EVERY later turn restarts
             fresh and the conversation loses its short-term memory for good.
             A turn N+1 after a fresh start has a small CURRENT context and
             must NOT be fresh again.
Defect 3 -- `fresh_context_for_turn` must return a block over budget so the
             terminal can arm ``engine.start_fresh`` (it returned "" before
             the fix because it only read the cumulative token_usage).
Requirement -- the fresh start keeps the CURRENT user message (the prompt the
             turn must answer is never dropped): "one user message = fresh
             block + current prompt".

The fake engine exposes the real engine surface that the fixed hook reads
(``_estimate_context_tokens``) AND a settable cumulative ``token_usage``
counter that the defective hook read, so each test is red against the
unfixed hook and green against the fixed one.
"""
from __future__ import annotations

from delfin.agent import task_state
from delfin.agent.repl import fresh_context_for_turn, _FRESH_CONTEXT_MARKER


class _FakeEngine:
    """Duck-typed engine for the hook: a CURRENT-context estimate (the real
    method ``_estimate_context_tokens``) plus a settable cumulative
    ``token_usage`` counter that the defective pre-fix hook consulted."""

    def __init__(self, current_tokens: int, cumulative_input: int = 0):
        self._current_tokens = current_tokens
        # token_usage['input'] is the CUMULATIVE session counter. The fixed
        # hook never reads it; the unfixed one did -- kept so the two can be
        # told apart.
        self.token_usage = {"input": cumulative_input, "output": 0}

    def _estimate_context_tokens(self) -> int:
        return int(self._current_tokens)


def _state(tmp_path):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task="review phase 3", phase="phase 3", phase_status="in_progress",
              commit="abc123")
    return st


def test_defect2_after_fresh_start_next_turn_is_not_fresh(tmp_path):
    """Turn N consumed the budget (fresh fired); turn N+1 has a SMALL current
    context now. Because the trigger reads the CURRENT context size and not the
    cumulative counter, turn N+1 must NOT be fresh again."""
    # Small current context, yet a huge cumulative counter that the defective
    # hook kept reading forever.
    engine = _FakeEngine(current_tokens=2_000, cumulative_input=950_000)
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == "", (
        "turn N+1 after a fresh start must keep the (now small) context; "
        "triggering on the cumulative counter would restart fresh forever")


def test_defect2_metric_is_current_context_not_cumulative(tmp_path):
    """The over-budget trigger is the CURRENT context estimate, not the running
    session counter. A session whose cumulative usage is huge but whose current
    context is small must NOT restart fresh."""
    engine = _FakeEngine(current_tokens=100_000, cumulative_input=2_500_000)
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""


def test_defect3_hook_returns_a_block_over_budget():
    """The hook must return a bounded fresh block once the CURRENT context is
    over the budget, so the terminal has something to arm with."""
    # Current context over the budget -> a non-empty block with the marker.
    engine = _FakeEngine(current_tokens=950_000, cumulative_input=950_000)
    ctx = fresh_context_for_turn(engine, None, token_budget=900_000)
    assert ctx and _FRESH_CONTEXT_MARKER in ctx
