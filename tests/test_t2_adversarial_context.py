"""T2 phase-3 adversarial review (reviewer s15): fresh_context_for_turn edges.

The builder's 9 control tests cover under/over budget, naming the open
phase, bounded output, working-state carry, secret scrub, determinism, and
broken-engine/state degradation. These adversarial cases probe two edges the
control suite does not hit:

1. the EXACT boundary ``current context == token_budget``, which the control
   tests skip (they use strictly-under 100_000 vs strictly-over 900_001) —
   the hook's contract is "at or under the budget -> keep full history";
2. an OVER-budget turn with ``task_state=None`` and no working-state block —
   the control suite only reduces to the marker via a working_state block,
   so the marker-only path (no task_state, no ws) is untested; it must still
   degrade to a usable block, not raise.

Both are expected GREEN on the builder's hook (review evidence, not findings).

The fake engine exposes the CURRENT context estimate via
``_estimate_context_tokens()``, the method ``_context_tokens`` (repl.py:419)
reads; its ``token_usage`` cumulative counter must not drive the trigger.
"""
from __future__ import annotations

import pytest

from delfin.agent import task_state
from delfin.agent.repl import fresh_context_for_turn, _FRESH_CONTEXT_MARKER


class _FakeEngine:
    def __init__(self, current_context: int):
        self._current_context = current_context
        # The cumulative session counter, kept for observability: the fixed
        # hook must NOT read it (a large cumulative counter with a small
        # current context must not fire a fresh start).
        self.token_usage = {"input": current_context, "output": 0}

    def _estimate_context_tokens(self) -> int:
        return self._current_context


def _state(tmp_path):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task="review phase 3", phase="phase 3", phase_status="in_progress",
              commit="abc123")
    return st


def test_boundary_exactly_at_budget_keeps_full_history(tmp_path):
    """current context == token_budget is AT the budget -> keep history ("")."""
    engine = _FakeEngine(current_context=900_000)
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""


def test_boundary_one_over_budget_fires_fresh_context(tmp_path):
    """One token past the budget -> the fresh context fires."""
    engine = _FakeEngine(current_context=900_001)
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert _FRESH_CONTEXT_MARKER in ctx


def test_over_budget_with_no_state_and_no_working_state_degrades_usable(tmp_path):
    """Over budget with task_state=None and no working_state: the hook must
    still hand back a usable marker-only block, never raise and never ""
    ("" would silently keep the huge history the caller asked to cut)."""
    engine = _FakeEngine(current_context=950_000)
    ctx = fresh_context_for_turn(engine, None, token_budget=900_000)
    assert ctx and _FRESH_CONTEXT_MARKER in ctx


def test_over_budget_blank_working_state_block_is_dropped(tmp_path):
    """A whitespace-only working_state_block must not emit a blank section."""
    engine = _FakeEngine(current_context=950_000)
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000,
                                 working_state_block="   \n\t  ")
    assert "\n\n\n" not in ctx  # no double-blank from the trimmed empty block
    assert _FRESH_CONTEXT_MARKER in ctx
