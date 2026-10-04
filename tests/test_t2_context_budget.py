"""T2 phase-3 tests: the repl context-budget restart hook.

A long session resumes with >900k input tokens and re-sends the whole
history every turn. Phase 3 gives the terminal a restart hook in ``repl.py``
—— ``fresh_context_for_turn`` —— that, once the CURRENT context is above a
token budget, hands back ONE small fresh context (``task_state.render()``,
which names the open phase, plus an optional working-state block) instead of
the full history. Below the budget it returns ``""`` and the caller keeps the
history. These tests drive the hook through the public path with a fake
engine and a real ``TaskState``.

The budget is judged on the CURRENT context (``_estimate_context_tokens``),
not the session's cumulative input-token usage: a cumulative number never
shrinks, so once it passed the budget every following turn would restart
fresh and the session would lose its memory for good. The current-context
estimate drops back under the budget the moment the fresh start happens, so
turn N+1 after a fresh start must NOT be fresh again.
"""
from __future__ import annotations

import pytest

from delfin.agent import task_state
from delfin.agent.repl import fresh_context_for_turn, _FRESH_CONTEXT_MARKER

_SECRET = "sk-proj-" + "phase-three-never-leaks-this"


class _FakeEngine:
    """Duck-typed engine exposing the CURRENT-context estimate the hook reads.

    ``messages`` mirror the engine transcript; ``_estimate_context_tokens``
    tallies them (plus an optional fixed ``system_prompt_chars``), so the
    estimate shrinks when the transcript is cut, exactly like the real
    engine's ``_estimate_context_tokens``. An explicit ``context_tokens``
    override lets a test set the current-context size directly without
    building a huge transcript.
    """

    def __init__(self, messages=None, system_prompt_chars: int = 0,
                 context_tokens=None):
        self.messages = list(messages) if messages else []
        self._system_prompt_chars = int(system_prompt_chars)
        self._context_override = context_tokens

    def _estimate_context_tokens(self) -> int:
        if self._context_override is not None:
            return max(1, int(self._context_override))
        total = 0
        for m in self.messages:
            if not isinstance(m, dict):
                continue
            total += len(m.get("content", "") or "")
        return max(1, total // 4 + self._system_prompt_chars // 4)

    @classmethod
    def long_session(cls):
        """A long session: current context reads above the 900k budget."""
        return cls(context_tokens=950_000)


def _state(tmp_path, *, task="build task_state.py", phase="phase 3"):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task=task, phase=phase, phase_status="in_progress",
              commit="abc123")
    return st


def test_under_budget_keeps_full_history(tmp_path):
    engine = _FakeEngine(context_tokens=100_000)  # long but under budget
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""  # caller keeps the full history


def test_at_budget_still_keeps_full_history(tmp_path):
    engine = _FakeEngine(context_tokens=900_000)  # exactly the budget
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx == ""


def test_over_budget_returns_fresh_context(tmp_path):
    engine = _FakeEngine.long_session()  # current context past the budget
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000)
    assert ctx
    assert _FRESH_CONTEXT_MARKER in ctx


def test_fresh_context_names_task_and_open_phase(tmp_path):
    engine = _FakeEngine.long_session()
    ctx = fresh_context_for_turn(engine, _state(tmp_path, task="fix the gap",
                                                 phase="phase 3"),
                                  token_budget=900_000)
    assert "fix the gap" in ctx       # the task survives
    assert "phase 3" in ctx           # the OPEN phase is named, not the whole list
    assert "phase 2" not in ctx       # no stale phase bleeds in


def test_fresh_context_is_bounded_not_the_full_history(tmp_path):
    """A long-session state must yield a small block, not a replay of
    every commit (the very failure this package exists to avoid)."""
    engine = _FakeEngine.long_session()
    st = _state(tmp_path)
    for i in range(200):
        st.commit(task="bulk", phase="phase 3", phase_status="in_progress",
                  commit=f"bulk{i:06x}")
    ctx = fresh_context_for_turn(engine, st, token_budget=900_000)
    assert len(ctx) <= 2200, "the fresh context must stay within the budget"
    assert "bulk000000" not in ctx  # oldest dropped by the render cap


def test_working_state_block_carries_into_fresh_context(tmp_path):
    engine = _FakeEngine.long_session()
    ws = "- Last test outcomes:\n  tests/test_t2_context_budget.py 6 passed"
    ctx = fresh_context_for_turn(engine, _state(tmp_path), token_budget=900_000,
                                 working_state_block=ws)
    assert "6 passed" in ctx  # the caller's working-state block is kept


def test_fresh_context_never_leaks_secret(tmp_path):
    engine = _FakeEngine.long_session()
    st = _state(tmp_path)
    st.add_finding(finding="render() leaks " + _SECRET)
    ctx = fresh_context_for_turn(engine, st, token_budget=900_000)
    assert _SECRET not in ctx


def test_fresh_context_is_deterministic(tmp_path):
    engine = _FakeEngine.long_session()
    st = _state(tmp_path)
    assert fresh_context_for_turn(engine, st, token_budget=900_000) == \
        fresh_context_for_turn(engine, st, token_budget=900_000)


def test_after_fresh_start_next_turn_is_not_fresh_again(tmp_path):
    """Defect #2: the budget must gate the CURRENT context, so once a fresh
    start shrinks it below the budget, turn N+1 resumes a normal history
    instead of restarting fresh forever."""
    st = _state(tmp_path)
    # N: long session — over budget, so a fresh start is signalled.
    engine = _FakeEngine.long_session()
    assert fresh_context_for_turn(engine, st, token_budget=900_000)
    # The fresh start cuts the transcript; the current context collapses.
    engine.messages = [{"role": "user", "content": "what next?"}]
    engine._context_override = None
    # N+1: now tiny, so the caller keeps the full history.
    assert fresh_context_for_turn(engine, st, token_budget=900_000) == ""


def test_broken_engine_or_missing_state_degrades_to_empty(tmp_path):
    class _NoMessages:
        """A fake engine with no messages / estimate at all — safe fallback."""
        def _estimate_context_tokens(self):
            raise AttributeError("not wired in a fake")

    # An engine that cannot size its context degrades to under-budget and
    # the caller keeps the full history.
    assert fresh_context_for_turn(_NoMessages(), None, token_budget=900_000) == ""
    # A state whose render() throws must not break the caller either.
    class _Boom:
        def render(self):
            raise RuntimeError("boom")
    engine = _FakeEngine.long_session()
    ctx = fresh_context_for_turn(engine, _Boom(), token_budget=900_000)
    assert _FRESH_CONTEXT_MARKER in ctx  # still a usable fresh context


# ---------------------------------------------------------------------------
# End-to-end through TerminalAgent.turn() with a fake engine (operator defect
# #3: the hook is wired into the turn, sets start_fresh, and fuses the block
# with the prompt the model must answer).
# ---------------------------------------------------------------------------

class _TurnEngine:
    """Enough of AgentEngine for a TerminalAgent.turn() run, reproducing the
    engine patch's ``start_fresh`` surgery so the whole route is covered: the
    terminal fuses the fresh block + prompt and arms ``start_fresh``; the
    engine keeps that one message, drops the rest, clears the flag and drops
    the token floor (so turn N+1 is judged on its shrunken context)."""

    token_usage = {"input": 0, "output": 0}
    session_id = "test-sid"

    def __init__(self, context_override=None, messages=None):
        self.messages = list(messages) if messages else []
        self._context_override = context_override
        self.start_fresh = False
        self.last_user_message = ""
        self.client = type("C", (), {"model": "m"})()
        self.fresh_contexts = 0

    def _estimate_context_tokens(self):
        if self._context_override is not None:
            return self._context_override
        total = sum(len(m.get("content", "") or "") for m in self.messages)
        return max(1, total // 4)

    def _refusal_entries(self):
        return []

    def stream_response(self, user_message="", **kw):
        self.last_user_message = user_message
        if self.start_fresh:
            self.fresh_contexts += 1
            # Mimics the engine patch: keep only the fused user message,
            # drop the accumulated history, clear the one-shot flag and the
            # token floor so the next estimate is small again.
            self.messages = [{"role": "user", "content": user_message}]
            self.start_fresh = False
            self._context_override = None
        else:
            self.messages.append({"role": "user", "content": user_message})
            self.messages.append({"role": "assistant", "content": "ok"})
        return "ok"

    def get_status(self):
        return {}

    def request_stop(self):
        pass

    def clear_stop(self):
        pass


def _turn_agent(engine, *, budget=0):
    import io
    from delfin.agent import repl
    out, err = io.StringIO(), io.StringIO()
    # color="never" makes Transcript pick a plain (no ANSI) theme by itself
    # via repl_render.theme_for, so nothing more needs forcing.
    agent = repl.TerminalAgent(
        engine, out=out, err=err,
        opts=repl.ReplOptions(cwd=__import__("pathlib").Path("."),
                              context_budget=budget, color="never"))
    return agent


def test_turn_over_budget_arms_fresh_restart_and_fuses_prompt(tmp_path):
    """Defect #3: the wired hook sets ``start_fresh`` and delivers the fresh
    block fused with the prompt the model must answer — as ONE message."""
    engine = _TurnEngine(context_override=950_000)  # over the 900k budget
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    result = agent.turn("continue the build")
    assert result.text == "ok"
    assert _FRESH_CONTEXT_MARKER in engine.last_user_message   # fresh block present
    assert "continue the build" in engine.last_user_message    # current prompt kept (defect #1)
    assert engine.last_user_message.endswith("continue the build")
    assert engine.start_fresh is False                         # one-shot consumed
    assert engine.fresh_contexts == 1


def test_turn_under_budget_keeps_history_and_does_not_arm(tmp_path):
    engine = _TurnEngine(context_override=50_000)  # well under the budget
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    result = agent.turn("is_it_fresh")
    assert result.text == "ok"
    assert engine.last_user_message == "is_it_fresh"   # no fresh block
    assert _FRESH_CONTEXT_MARKER not in engine.last_user_message
    assert getattr(engine, "start_fresh", False) is False


def test_disabled_budget_never_fresh(tmp_path):
    engine = _TurnEngine(context_override=950_000)  # would be over, but hook off
    agent = _turn_agent(engine, budget=0)           # disabled (the default)
    agent._task_state = _state(tmp_path)
    result = agent.turn("plain")
    assert result.text == "ok"
    assert engine.last_user_message == "plain"
    assert _FRESH_CONTEXT_MARKER not in engine.last_user_message


def test_turn_after_a_fresh_start_is_not_fresh_again(tmp_path):
    """Defect #2 end-to-end: the fresh start shrinks the CURRENT context, so
    the very next turn is judged on its small history and is NOT fresh."""
    engine = _TurnEngine(context_override=950_000)
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    first = agent.turn("big question")          # N: over budget → fresh
    assert _FRESH_CONTEXT_MARKER in engine.last_user_message
    assert engine.fresh_contexts == 1
    second = agent.turn("what next?")           # N+1: context now tiny
    assert second.text == "ok"
    assert engine.last_user_message == "what next?"
    assert _FRESH_CONTEXT_MARKER not in engine.last_user_message
    assert engine.fresh_contexts == 1           # no second fresh restart
