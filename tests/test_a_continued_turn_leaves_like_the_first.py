"""A self-continued turn leaves the way the first one does.

``TerminalAgent.run`` ends a turn that stopped at the token ceiling
without a tool call by re-sending "continue", at most twice
(``_continue_after_length``). The FIRST turn sits inside a
``try/except KeyboardInterrupt`` that returns 130; the continuation
turns, added later, were called bare inside ``while continuation:``.

Both interrupt roads end in the same place, so the asymmetry is real:

* **Esc** arrives as a keystroke: cbreak keeps ISIG (repl_keys.py's
  stated design), Esc maps to ``rk.INTERRUPT``, and ``_on_key`` calls
  ``_stop_engine`` -- a request, never an exception.
* **Ctrl+C** arrives as SIGINT: ``_on_sigint`` counts. The first press
  requests the stop, the second warns, and the THIRD raises
  ``KeyboardInterrupt`` out of ``_pump`` -- out of ``self.turn`` -- on
  the main thread, regardless of which turn (first or continued) is
  running.

So a user pressing Ctrl+C three times during a CONTINUED turn crashed
out of ``run()`` entirely -- through the finally, past every
``except KeyboardInterrupt`` -- instead of leaving the way the same
three presses during the first turn leave: status 130, the session
saved for a resume.

Judged here:

  the escape hatch    a ``KeyboardInterrupt`` out of the continued
                      turn makes ``run()`` return 130 and say
                      "(interrupt)" on screen, exactly as the first
                      turn does -- never raise out of ``run()``

  why the first turn's handler is not enough
                      the handler covers ``self.turn(pending)`` only;
                      the continuation loop calls ``self.turn`` outside
                      it. The fix moves the hatch to cover both, not a
                      second ladder.
"""

from __future__ import annotations

from pathlib import Path

import pytest


class _LengthThenInterruptEngine:
    """Turn 1 ends at ``length`` (no tool call); turn 2 is the
    continuation. The interrupt is delivered by the agent subclass
    below, on the main thread, exactly where the third Ctrl+C raises."""

    session_id = "length-interrupt-0001"
    token_usage = {"input": 0, "output": 0}
    last_turn_stop_reason = ""

    def __init__(self):
        self.messages = []
        self.turns = 0

    def stream_response(self, user_message="", max_tokens=0, **kw):
        self.turns += 1
        self.messages.append({"role": "user", "content": user_message})
        self.messages.append({"role": "assistant", "content": "half"})
        self.last_turn_stop_reason = "length" if self.turns == 1 \
            else "end_turn"
        return "half"

    def get_status(self):
        return {}

    def export_state(self):
        return {"engine_messages": list(self.messages),
                "token_usage": dict(self.token_usage)}


def _agent_class_that_interrupts_the_continuation(repl):
    """A third Ctrl+C during the CONTINUED turn: ``turn()`` raises
    ``KeyboardInterrupt`` on the main thread, which is the one place
    the ladder's last step can surface (the worker thread swallows
    BaseException, so the raise has to come from the pump's side)."""

    class _Agent(repl.TerminalAgent):
        interrupt_on_turn = 2          # 1-based: the continuation

        def turn(self, prompt):
            if getattr(self, "_turn_call", 0) + 1 >= self.interrupt_on_turn:
                raise KeyboardInterrupt
            self._turn_call = getattr(self, "_turn_call", 0) + 1
            return super().turn(prompt)

    return _Agent


def _run(tmp_path, interrupt_on_turn):
    import io
    from delfin.agent import repl as R

    engine = _LengthThenInterruptEngine()
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    lines = ["do the big edit"]

    def _read_line(_prompt=""):
        if lines:
            return lines.pop(0)
        raise KeyboardInterrupt

    agent_cls = _agent_class_that_interrupts_the_continuation(R)
    agent = agent_cls(
        engine, opts=R.ReplOptions(cwd=tmp_path, max_tokens=0),
        out=out, err=err, read_line=_read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    agent._idle_interrupts = 1           # one more interrupt leaves (130)
    agent.interrupt_on_turn = interrupt_on_turn
    code = agent.run()
    return code, engine, out


def test_an_interrupt_in_the_continued_turn_leaves_with_130(tmp_path):
    """The continued turn is a turn like the first: the third Ctrl+C
    during it returns 130 from run(), never raises out of it."""
    code, engine, _out = _run(tmp_path, interrupt_on_turn=2)
    assert code == 130, (
        "a third Ctrl+C during a self-continued turn must leave the "
        "session the way the same presses during the first turn do "
        "(return 130); raising out of run() skips the prompt-side "
        "handler and every checkpoint after it")


def test_the_interrupt_lands_in_the_continuation_not_the_first_turn(
        tmp_path):
    """Guard for the test above: the raise must arrive in the
    CONTINUATION -- after one length-ended turn completed -- or the
    first test would pass for the old, already-covered reason."""
    code, engine, _out = _run(tmp_path, interrupt_on_turn=2)
    # One turn completed fully (its stream_response ran); the raise was
    # in the second turn() call, which _turn_call counted as done.
    assert engine.turns == 1, (
        "the fixture's interrupt fires in the second turn() call; if "
        "the first turn never completed, the test above is testing the "
        "first-turn handler, which was never broken")
    assert engine.last_turn_stop_reason == "length", (
        "the first turn must end at the ceiling -- that is what makes "
        "the second call a continuation")
