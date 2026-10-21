"""T2 phase-3 operator review (reviewer s15): the real caller + archive.

The operator's review of the rejected phase-3 engine patch set four
requirements for the corrected fix, three of which live in THIS package's
owned file (repl.py):

  - "one user message = fresh block + current prompt" -- the prompt the
    current turn must answer is NEVER dropped;
  - the trigger is the CURRENT context size, not the cumulative session
    input counter (a turn N+1 after a fresh start is NOT fresh again);
  - a REAL caller in repl.py with an end-to-end test -- TerminalAgent.turn()
    wires ``fresh_context_for_turn`` and arms ``engine.start_fresh``;
  - the cut history is ARCHIVED (retrievable), not silently discarded.

These tests drive the REAL public turn path TerminalAgent.turn() (repl.py:1007)
with a fake engine that reproduces the engine patch's ``start_fresh`` surgery
so every part of the route is exercised. They are RED against the pre-fix
turn() (which had no caller: start_fresh was never armed and nothing was
archived) and GREEN against the fixed one.
"""
from __future__ import annotations

import io
from pathlib import Path

from delfin.agent import repl, task_state
from delfin.agent.repl import _FRESH_CONTEXT_MARKER


class _TurnEngine:
    """Enough of AgentEngine for a TerminalAgent.turn() run, exposing the
    real surface the fixed repl hook reads (``_estimate_context_tokens``,
    ``start_fresh``, ``_refusal_entries``) and reproducing the engine patch's
    ``start_fresh`` surgery: when armed, the accumulated history is ARCHIVED
    (``cut_history()``) and only the fused fresh+prompt message is kept."""

    token_usage = {"input": 0, "output": 0}
    session_id = "test-sid"

    def __init__(self, current_context: int | None = None, *, budget_over: bool = False):
        self.messages: list[dict] = []
        self._context_override = current_context
        self.start_fresh = False
        self.last_user_message = ""
        self._archived: list[dict] = []
        self.fresh_contexts = 0
        self.client = type("C", (), {"model": "m"})()

    def _estimate_context_tokens(self) -> int:
        if self._context_override is not None:
            # Simulate the CURRENT-context estimate the fixed hook reads.
            if self._context_override == "SHRUNK_AFTER_FRESH":
                # After a fresh start the messages are tiny again.
                total = sum(len(m.get("content", "") or "") for m in self.messages)
                return max(1, total // 4)
            return self._context_override
        total = sum(len(m.get("content", "") or "") for m in self.messages)
        return max(1, total // 4)

    def _refusal_entries(self):
        return []

    def cut_history(self) -> list:
        """The history an over-budget turn had to drop, archived."""
        return list(self._archived)

    def stream_response(self, user_message="", **kw):
        self.last_user_message = user_message
        if self.start_fresh:
            self.fresh_contexts += 1
            # Engine-patch surgery: the accumulated history is the session
            # record -- ARCHIVE it before dropping, and keep only the fused
            # fresh-block + current-prompt message.
            self._archived = self._archived + list(self.messages)
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
    out, err = io.StringIO(), io.StringIO()
    return repl.TerminalAgent(
        engine, out=out, err=err,
        opts=repl.ReplOptions(cwd=Path("."), context_budget=budget, color="never"))


def _state(tmp_path):
    st = task_state.open(tmp_path / "task_state.json")
    st.begin()
    st.commit(task="review phase 3", phase="phase 3", phase_status="in_progress",
              commit="abc123")
    return st


def test_end_to_end_over_budget_turn_arms_fresh_and_fuses_prompt(tmp_path):
    """A REAL caller: an over-budget CURRENT context makes turn() arm
    engine.start_fresh and deliver the fresh block fused with the prompt as
    ONE message — the prompt is never dropped. pre-fix turn() had no caller,
    so start_fresh was never armed (RED)."""
    engine = _TurnEngine(current_context=1_200_000)  # over the 900k budget
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    result = agent.turn("continue the build")
    assert result.text == "ok"
    assert engine.last_user_message.endswith("continue the build")  # prompt kept
    assert _FRESH_CONTEXT_MARKER in engine.last_user_message        # block present
    assert engine.start_fresh is False                              # one-shot consumed
    assert engine.fresh_contexts == 1


def test_end_to_end_fresh_turn_archives_the_cut_history(tmp_path):
    """The history a fresh start drops is ARCHIVED (retrievable via
    engine.cut_history()), not silently lost. RED on the pre-fix tree where
    nothing armed start_fresh and nothing archived."""
    engine = _TurnEngine(current_context=1_200_000)
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    # A session with real accumulated history before this turn.
    engine.messages = [{"role": "user", "content": "old q"},
                       {"role": "assistant", "content": "old a"}]
    agent.turn("current prompt")
    assert engine.cut_history(), (
        "the dropped history must be archived, not silently lost")


def test_end_to_end_under_budget_cumulative_huge_not_fresh(tmp_path):
    """The trigger is the CURRENT context size, not the cumulative session
    counter: a session with a huge cumulative input count but a small CURRENT
    context must NOT arm a fresh start — the current prompt goes through
    unchanged."""
    engine = _TurnEngine(current_context=100_000)   # small CURRENT context
    engine.token_usage = {"input": 2_500_000, "output": 0}  # cumulative huge
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    agent.turn("plain prompt")
    assert engine.start_fresh is False
    assert engine.last_user_message == "plain prompt"     # full history kept
    assert _FRESH_CONTEXT_MARKER not in engine.last_user_message


def test_end_to_end_turn_after_fresh_start_is_not_fresh_again(tmp_path):
    """Defect #3: after a fresh start the CURRENT context shrinks, so the very
    next turn is judged on its small history and is NOT fresh again."""
    engine = _TurnEngine(current_context=1_200_000)
    agent = _turn_agent(engine, budget=900_000)
    agent._task_state = _state(tmp_path)
    first = agent.turn("big question")   # N: over budget -> fresh
    assert _FRESH_CONTEXT_MARKER in engine.last_user_message
    assert engine.fresh_contexts == 1
    # N+1: the fresh start shrank the context; the estimate is now tiny.
    engine._context_override = "SHRUNK_AFTER_FRESH"
    second = agent.turn("what next?")
    assert second.text == "ok"
    assert engine.last_user_message == "what next?"   # not fresh again
    assert _FRESH_CONTEXT_MARKER not in engine.last_user_message
    assert engine.fresh_contexts == 1
