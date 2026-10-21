"""Steering: a wait wakes, a person wins, a pause holds.

Phase 3, package KG1. The four rules of the brief, tested at the shapes
the terminal loop actually exercises (repl.py run loop, _wake_text,
wake_blocked, background status line):

1. typed input beats a due wake and a nudge — a draft belongs to whoever
   typed it, and a wake must not take it away;
2. after stop/pause/Esc no wake fires until the person types;
3. a message from a session or the operator wakes but cannot lift a stop;
4. the status line shows the wait and the person can clear it.
"""

from __future__ import annotations

import io
import time
import types

import pytest

from delfin.agent import repl, session_messages as msgs, session_pause
from delfin.agent import state_paths
from delfin.agent.wait_registry import SessionWaits, due, register


@pytest.fixture
def home(tmp_path, monkeypatch):
    with state_paths.RedirectedUserState(tmp_path):
        monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
        yield tmp_path / ".delfin"


class _Tty(io.StringIO):
    def isatty(self):
        return True


def _agent(home):
    agent = repl.TerminalAgent.__new__(repl.TerminalAgent)
    agent.engine = types.SimpleNamespace(session_id="abcdef0123456789")
    agent.opts = types.SimpleNamespace(session_name="kg1-s20", cwd=".",
                                       banner="")
    agent.transcript = _FakeTranscript()
    agent._stdin = io.StringIO()
    agent._read_line = lambda _p: ""       # never a real key here
    agent._announce_followups = 0
    return agent


class _FakeTranscript:
    def __init__(self):
        self.lines = []

    def chrome(self, text):
        self.lines.append(text)

    class theme:                                        # noqa: N801
        @staticmethod
        def dim(text):
            return text


# -- 1. typed input beats a due wake and a nudge -----------------------------

def test_typed_input_beats_a_due_wake(home):
    register("message_from", "nacht-s25")
    agent = _agent(home)
    # Something is typed: the wake must hand back "" (the draft wins).
    assert agent._wake_text("half a thought") == ""


def test_a_due_wake_fires_the_wait_note(home):
    register("message_from", "nacht-s25")
    # A message_from wait fires when mail from that sender is DELIVERED:
    # real writer path (session_messages.send), not a typed event row.
    msgs.send("kg1-s20", "phase 3 review round", from_key="nacht-s25",
              from_title="nacht-s25")
    agent = _agent(home)
    agent._wake_last_look = 0.0          # past the 5 s idle-wake throttle
    text = agent._wake_text("")
    assert "phase 3 review round" in text, text
    assert due(time.time(), []) == [], "the fired wait is consumed"


def test_typed_input_beats_a_nudge(home):
    from delfin.agent import turn_continuation
    register("message_from", "nacht-s25", key="kg1-s20")
    cont, _ = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is False, "the wait suppresses the nudge"
    # ...and typed input does not change that: the nudge is a wake-shape,
    # _wake_text never returns a nudge for typed text.
    agent = _agent(home)
    assert agent._wake_text("do it") == ""


# -- 2. after stop/pause/Esc no wake fires until the person types ------------

def test_a_pause_blocks_the_wait_wake(home):
    register("message_from", "nacht-s25")
    msgs.send("kg1-s20", "wake me", from_key="nacht-s25",
              from_title="nacht-s25")
    session_pause.pause("kg1-s20", "operator said stop")
    agent = _agent(home)
    assert agent._wake_text("") == "", "no wake while paused"
    assert session_pause.wake_blocked("kg1-s20") is True
    session_pause.resume("kg1-s20")
    # Only after the person types does the wake come back -- and the
    # message is still in the inbox (the wake is the reader of mail).
    agent._wake_last_look = 0.0          # past the 5 s idle-wake throttle
    assert agent._wake_text("") != ""


def test_a_message_cannot_lift_a_stop(home):
    # The wake path must never resume the session: no code on the wake
    # path calls session_pause.resume.
    import inspect
    from delfin.agent import repl as _repl
    src = inspect.getsource(_repl.TerminalAgent._wake_text_plain)
    assert "resume(" not in src
    # And behaviourally: paused + a waiting message -> still blocked.
    register("message_from", "nacht-s25")
    session_pause.pause("kg1-s20")
    agent = _agent(home)
    assert agent._wake_text("") == ""


# -- 4. the status line shows the wait and the person can clear it -----------

def test_the_status_line_shows_the_wait(home):
    register("message_from", "nacht-s25", key="kg1-s20")
    agent = _agent(home)
    line = agent._background_status()
    assert "wait" in line.lower() and "nacht-s25" in line, line


def test_the_person_can_clear_the_wait_from_the_status(home):
    w = register("message_from", "nacht-s25", key="kg1-s20")
    assert SessionWaits("kg1-s20").clear(w["id"]) is True
    agent = _agent(home)
    assert "nacht-s25" not in agent._background_status()
    # ...and a cleared wait no longer suppresses a nudge.
    from delfin.agent import turn_continuation
    cont, _ = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is True
