"""T2 phase-3 operator review (reviewer s15): the real caller + archive.

The operator reviewed the rejected phase-3 engine patch and set the bar for
the corrected fix coming out of repl.py (the restart hook lives in this
package's owned file):

  - "one user message = fresh block + current prompt" -- the prompt the
    current turn must answer is NEVER dropped;
  - the trigger is the CURRENT context size, not the cumulative session
    input counter (a turn N+1 after a fresh start is NOT fresh again);
  - a REAL caller in repl.py with an end-to-end test -- the hook is wired
    into the turn path, not left dead;
  - the cut history is ARCHIVED (retrievable), not silently discarded.

These tests pin that bar end-to-end through ``run_turn`` -- the public turn
path that ``TerminalAgent.turn()`` delegates to (repl.py:927). They are RED
on the current tree, where ``run_turn`` sends the full history and never
consults ``fresh_context_for_turn`` nor arms a fresh-start flag.
"""
from __future__ import annotations

from delfin.agent.repl import run_turn


class _RunEngine:
    """Duck-typed engine driving ``run_turn``: records what the turn path
    actually hands to ``stream_response`` and exposes the NEW contract the
    corrected hook must arm (``start_fresh``, ``cut_history``)."""

    def __init__(self, current_tokens: int, history: list | None = None):
        self.token_usage = {"input": current_tokens, "output": 0}
        self.messages = list(history) if history is not None else []
        self.start_fresh = ""          # armed by the corrected repl hook
        self._archived: list = []
        self.sent_user_message: str | None = None
        self.sent_history_len: int | None = None

    def estimate_context_tokens(self) -> int:
        """CURRENT context estimate (the next request), not cumulative."""
        return int(self.token_usage["input"])

    def cut_history(self) -> list:
        """The history an over-budget turn had to drop, archived."""
        return list(self._archived)

    def stream_response(self, **kwargs) -> str:
        self.sent_user_message = kwargs.get("user_message", "")
        self.sent_history_len = len(kwargs.get("history", []))
        return "answer"


def _sink(_item):
    pass


def test_end_to_end_over_budget_turn_arms_a_fresh_start():
    """An over-budget CURRENT context must make the turn path arm a fresh
    start (set ``engine.start_fresh`` to a non-empty block) so the engine can
    ship the small fresh context instead of the full history. RED on the
    current tree: run_turn sends the whole history and never arms it."""
    engine = _RunEngine(current_tokens=1_200_000)  # current context over budget
    run_turn(engine, "the-prompt", sink=_sink)
    assert engine.start_fresh, (
        "an over-budget turn must arm engine.start_fresh (real caller); "
        "it is never set on the current tree")


def test_end_to_end_fresh_turn_archives_the_cut_history():
    """When a fresh start drops the old history, that history is ARCHIVED
    (retrievable via ``engine.cut_history()``), not silently discarded.
    RED on the current tree: nothing archives anything."""
    engine = _RunEngine(
        current_tokens=1_200_000,
        history=[{"role": "user", "content": "old q"},
                 {"role": "assistant", "content": "old a"}])
    run_turn(engine, "current-prompt", sink=_sink)
    assert engine.cut_history(), (
        "the dropped history must be archived, not silently lost; "
        "cut_history() is empty on the current tree")


def test_end_to_end_fresh_turn_keeps_the_current_user_message():
    """After a fresh start the CURRENT user message (the prompt this turn
    answers) still reaches the model. The corrected patch must keep "one user
    message = fresh block + current prompt". This is a green guard on the
    currently-shipping path (which never drops the prompt) AND must stay true
    once the fresh-start swap lands."""
    engine = _RunEngine(current_tokens=1_200_000)
    run_turn(engine, "this-exact-prompt", sink=_sink)
    assert "this-exact-prompt" in (engine.sent_user_message or "")


def test_end_to_end_budget_is_current_context_not_cumulative():
    """The end-to-end trigger is the CURRENT context size, not the cumulative
    session counter: a session whose current context is small after a fresh
    restart must NOT arm another fresh start (else every later turn restarts
    and short-term memory is lost)."""
    engine = _RunEngine(current_tokens=100_000)   # small current context
    engine.token_usage = {"input": 2_500_000, "output": 0}  # cumulative huge
    run_turn(engine, "prompt", sink=_sink)
    assert engine.start_fresh == "", (
        "a small CURRENT context must not arm another fresh start even when "
        "the cumulative session counter is huge")
