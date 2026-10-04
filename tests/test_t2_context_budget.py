"""T2 phase-3 tests: the repl context-budget restart hook.

A long session resumes with >900k input tokens and re-sends the whole
history every turn. Phase 3 gives the terminal a restart hook in ``repl.py``
—— ``fresh_context_for_turn`` —— that, once input-token usage is above a
token budget, hands back ONE small fresh context (``task_state.render()``,
which names the open phase, plus an optional working-state block) instead of
the full history. Below the budget it returns ``""`` and the caller keeps the
history. These tests drive the hook through the public path with a fake
engine and a real ``TaskState``.
"""
from __future__ import annotations

import pytest

from delfin.agent import task_state
from delfin.agent.repl import fresh_context_for_turn, _FRESH_CONTEXT_MARKER

_SECRET = "sk-proj-" + "phase-three-never-leaks-this"


class _FakeEngine:
    """Duck-typed engine: only ``.token_usage`` is read by the hook."""

    def __init__(self, input_tokens: int):
        self.token_usage = {"input": input_tokens, "output": 0}


def _state(tmp_path, *, task="build task_state.py", phase="phase 3"):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task=task, phase=phase, phase_status="in_progress",
              commit="abc123")
    return st


def test_under_budget_keeps_full_history(tmp_path):
    engine = _FakeEngine(input_tokens=100_000)  # a long but under-budget turn
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""  # caller keeps the full history


def test_at_budget_still_keeps_full_history(tmp_path):
    engine = _FakeEngine(input_tokens=900_000)  # exactly the budget
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""


def test_over_budget_returns_fresh_context(tmp_path):
    engine = _FakeEngine(input_tokens=900_001)  # past the budget
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx
    assert _FRESH_CONTEXT_MARKER in ctx


def test_fresh_context_names_task_and_open_phase(tmp_path):
    engine = _FakeEngine(input_tokens=950_000)
    ctx = fresh_context_for_turn(engine, _state(tmp_path, task="fix the gap",
                                                phase="phase 3"),
                                 token_budget=900_000)
    assert "fix the gap" in ctx       # the task survives
    assert "phase 3" in ctx           # the OPEN phase is named, not the whole list
    assert "phase 2" not in ctx       # no stale phase bleeds in


def test_fresh_context_is_bounded_not_the_full_history(tmp_path):
    """A long-session state must yield a small block, not a replay of
    every commit (the very failure this package exists to avoid)."""
    engine = _FakeEngine(input_tokens=950_000)
    st = _state(tmp_path)
    for i in range(200):
        st.commit(task="bulk", phase="phase 3", phase_status="in_progress",
                  commit=f"bulk{i:06x}")
    ctx = fresh_context_for_turn(engine, st, token_budget=900_000)
    assert len(ctx) <= 2200, "the fresh context must stay within the budget"
    assert "bulk000000" not in ctx  # oldest dropped by the render cap


def test_working_state_block_carries_into_fresh_context(tmp_path):
    engine = _FakeEngine(input_tokens=950_000)
    ws = "- Last test outcomes:\n  tests/test_t2_context_budget.py 6 passed"
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000,
                                 working_state_block=ws)
    assert "6 passed" in ctx  # the caller's working-state block is kept


def test_fresh_context_never_leaks_secret(tmp_path):
    engine = _FakeEngine(input_tokens=950_000)
    st = _state(tmp_path)
    st.add_finding(finding="render() leaks " + _SECRET)
    ctx = fresh_context_for_turn(engine, st, token_budget=900_000)
    assert _SECRET not in ctx


def test_fresh_context_is_deterministic(tmp_path):
    engine = _FakeEngine(input_tokens=950_000)
    st = _state(tmp_path)
    assert fresh_context_for_turn(engine, st, token_budget=900_000) == \
        fresh_context_for_turn(engine, st, token_budget=900_000)


def test_broken_engine_or_missing_state_degrades_to_empty(tmp_path):
    class _NoUsage:
        """A fake engine with no token_usage at all — the safe fallback."""
        pass

    # No usage info: treat as under budget and keep the full history.
    assert fresh_context_for_turn(_NoUsage(), None, token_budget=900_000) == ""
    # A state whose render() throws must not break the caller either.
    class _Boom:
        def render(self):
            raise RuntimeError("boom")
    engine = _FakeEngine(input_tokens=950_000)
    ctx = fresh_context_for_turn(engine, _Boom(), token_budget=900_000)
    assert _FRESH_CONTEXT_MARKER in ctx  # still a usable fresh context
