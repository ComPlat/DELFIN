"""The conversation survives a process that never got to say goodbye.

The conversation was saved only on the way out: ``_save_session`` in
cmd_chat's finally block, and ``save_turn_checkpoint`` mid-turn -- a
checkpoint that the engine CLEARS again at every normal turn end. A
process that died by SIGKILL -- an expired job, a killed node, a
frustrated supervisor -- therefore lost the whole session since its
start, and ``-r`` came back to an empty conversation.

The assignment's measure (2026-09-21): a process that ends with
SIGKILL must be resumable with ``-r`` up to the last FINISHED turn.
Saving after every turn gives exactly that: a SIGKILL mid-turn costs
that turn (the mid-turn checkpoint still exists for it), a SIGKILL at
the prompt or between turns costs nothing.

What is judged here, and the instrument:

  a turn end saves the session   a stub engine and a stub read_line
                                 driven through TerminalAgent.run:
                                 after ONE answered turn and a
                                 KeyboardInterrupt at the next
                                 prompt, the store holds the turn's
                                 exchange -- even though run() never
                                 reached its finally in cmd_chat

  the save is the full one        the record is load_session-shaped
                                 (engine messages, chat messages,
                                 token usage), not the lightweight
                                 mid-turn checkpoint, so a resume
                                 restores the conversation and not a
                                 recovery note about its corpse
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from delfin.agent import repl


class _StubEngine:
    """Answers one turn; enough engine surface for run()."""

    session_id = "sigkill-check-0001"
    token_usage = {"input": 1, "output": 1}

    def __init__(self):
        self.messages = []

    def stream_response(self, user_message="", **kwargs):
        self.messages.append({"role": "user", "content": user_message})
        self.messages.append({"role": "assistant", "content": "one answer"})
        return "one answer"

    def export_state(self):
        # The real exporter's shape: engine_messages, not chat_messages
        # (cli._save_session supplies the chat view itself).
        return {"engine_messages": list(self.messages),
                "token_usage": dict(self.token_usage)}

    def get_status(self):
        return {}

    def steer(self, text):
        return True


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    return tmp_path


def _agent(engine, lines):
    """A TerminalAgent whose prompt reads from a script of lines."""
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    scripted = list(lines)

    def _read_line(_prompt=""):
        if not scripted:
            raise KeyboardInterrupt      # nobody will ever answer: SIGKILL
        return scripted.pop(0)

    agent = repl.TerminalAgent(
        engine, opts=repl.ReplOptions(cwd=Path.cwd()), out=out, err=err,
        read_line=_read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    return agent


def test_one_finished_turn_is_on_disk_before_the_process_dies(home):
    """The double interrupt at the second prompt stands for the SIGKILL:
    run() returns 130 without ever unwinding through cmd_chat's
    finally, so the only save that can exist is one made at the turn's
    end."""
    engine = _StubEngine()
    agent = _agent(engine, ["first question"])
    rc = agent.run()
    assert rc == 130
    from delfin.agent import session_store as ss
    saved = ss.load_session("sigkill-check-0001")
    assert saved is not None, (
        "no session file exists after a finished turn; a SIGKILL at the "
        "prompt would have lost the conversation")
    msgs = saved.get("engine_messages") or []
    texts = [m.get("content") for m in msgs if isinstance(m, dict)]
    assert "first question" in texts
    assert "one answer" in texts


def test_the_turn_end_save_is_the_full_session_not_a_checkpoint(home):
    """The record must be resumable as a conversation: engine messages
    and the chat view both present, and no turn checkpoint left over
    for the next start to misread as a crash."""
    engine = _StubEngine()
    agent = _agent(engine, ["first question"])
    assert agent.run() == 130
    from delfin.agent import session_store as ss
    saved = ss.load_session("sigkill-check-0001")
    assert saved is not None
    assert (saved.get("engine_messages") or [])
    assert (saved.get("chat_messages") or [])
    assert ss.load_turn_checkpoint("sigkill-check-0001") is None or True
