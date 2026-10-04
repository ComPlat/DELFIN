"""A line typed during a turn goes before any automatic follow-up turn.

A supervisor (or the user) typed /exit while a turn ran; the turn ended with
"Let me ..." and the announced-work follow-up ran first, so the session kept
working past its own /exit. Typed input -- above all /exit -- always wins
over every automatic continuation.
"""
from __future__ import annotations

import io

from delfin.agent import repl, repl_render as rr


class _Perms:
    def __init__(self):
        self.mode = "default"


class _Engine:
    def __init__(self):
        self.kit_permissions = _Perms()
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.pending_turn_continuation = ""

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode

    def request_stop(self):
        self.stopped = True

    def clear_turn_continuation(self):
        self.pending_turn_continuation = ""


class _Tty(io.StringIO):
    def isatty(self):
        return True


def _agent(engine):
    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=_Tty())
    agent.transcript.theme = rr.Theme(enabled=False)
    return agent


def test_a_queued_exit_ends_the_session_before_a_follow_up_turn():
    engine = _Engine()
    agent = _agent(engine)
    reads = iter(["do the task"])
    agent._read_line = lambda _p: next(reads)
    turns: list[str] = []

    def turn(prompt):
        turns.append(prompt)
        # The turn announces more work, and /exit is typed meanwhile.
        engine.pending_turn_continuation = "continue the announced step"
        agent.queued.append("/exit")
        return repl.TurnResult()

    agent.turn = turn
    assert agent.run() == 0
    assert turns == ["do the task"], turns


def test_without_typed_input_the_follow_up_still_runs():
    engine = _Engine()
    agent = _agent(engine)
    reads = iter(["do the task", "/exit"])
    agent._read_line = lambda _p: next(reads)
    turns: list[str] = []

    def turn(prompt):
        turns.append(prompt)
        if len(turns) == 1:
            engine.pending_turn_continuation = "continue the announced step"
        return repl.TurnResult()

    agent.turn = turn
    agent.run()
    assert len(turns) == 2 and turns[1] == "continue the announced step", turns
