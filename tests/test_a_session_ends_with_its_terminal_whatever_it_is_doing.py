"""A session ends with the terminal it runs in, whatever it is doing.

Seen while supervising, night of 2026-09-20: five sessions ran in tmux
windows, each under `script` for a recording. Overnight the tmux server
went away. In the morning all five agents were still running, eleven
hours on, and so were the five `script` processes above them; they ended
only when each was sent SIGTERM by hand.

Measured the next morning, without DELFIN and without a model: `script`
blocks SIGHUP, so it outlives its terminal, and all it passes down is a
single ctrl+d byte. At the idle prompt that ends a session cleanly. All
five sat in an approval dialog, which reads one key, ignores anything it
does not offer -- ctrl+d included -- and waits on. The turn loop dropped
ctrl+d the same way. And the reader could not tell a terminal that had
hung up from one where nobody was typing: both came back as "".

The rule: whatever the session is doing -- at the prompt, in a turn, in a
dialog -- a terminal that has gone ends it. There are three ways to learn
that it went: ctrl+d, a read that finds the line hung up, and SIGHUP or
SIGTERM. All three leave through the same door: a question in flight is
refused, the turn is stopped and not waited out, the loop returns through
its cleanup (so the session is saved and can be resumed), and a leave
that itself hangs is cut short.
"""

from __future__ import annotations

import io
import os
import signal
import threading
import time

import pytest

from delfin.agent import repl, repl_keys as rk, repl_render as rr
from delfin.agent import terminal_confirm as tc


PLAIN = rr.Theme(enabled=False)


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Engine:
    def __init__(self):
        self.kit_permissions = type("P", (), {"mode": "plan"})()
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.stopped = False

    def get_status(self):
        return {"input_tokens": 0, "output_tokens": 0, "cost_usd": 0.0}

    def request_stop(self):
        self.stopped = True

    def clear_stop(self):
        pass


class _Keys:
    """Scripted keystrokes; running out means the reader waited on."""

    active = True

    def __init__(self, keys):
        self._keys = list(keys)

    def read_ready(self, timeout):
        if not self._keys:
            raise AssertionError("still waiting for a key after ctrl+d")
        return self._keys.pop(0)


def _agent(engine=None, broker=None):
    engine = engine if engine is not None else _Engine()
    out, err = io.StringIO(), _Tty()
    agent = repl.TerminalAgent(engine, out=out, err=err, broker=broker)
    agent.transcript.theme = PLAIN
    return agent, engine, err


@pytest.fixture
def line():
    """A real pty whose far end the test can hang up."""
    pytest.importorskip("termios")
    pty = pytest.importorskip("pty")
    master, slave = pty.openpty()
    stream = os.fdopen(slave, "rb", buffering=0)
    held = {"master": master}
    try:
        yield stream, held
    finally:
        stream.close()
        if held["master"] is not None:
            os.close(held["master"])


# ---------------------------------------------------------------------------
# The reader tells a hangup from quiet
# ---------------------------------------------------------------------------

class TestTheReaderTellsAHangupFromQuiet:
    def test_nobody_typing_is_not_a_hangup(self, line):
        stream, _held = line
        with rk.RawMode(stream) as raw:
            assert raw.active
            assert raw.read_ready(0.05) == ""

    def test_a_line_that_hung_up_has_left(self, line):
        stream, held = line
        with rk.RawMode(stream) as raw:
            os.close(held["master"])
            held["master"] = None
            with pytest.raises(rk.TerminalLeft):
                raw.read_ready(0.5)


# ---------------------------------------------------------------------------
# ctrl+d is leaving, wherever it arrives
# ---------------------------------------------------------------------------

class TestCtrlDIsLeaving:
    def test_in_an_approval_dialog(self):
        agent, _engine, _err = _agent()
        with pytest.raises(rk.TerminalLeft):
            agent._read_key(_Keys(["\x04"]), {"y", "n", "\x1b"})

    def test_in_a_numbered_question(self):
        agent, _engine, _err = _agent()
        with pytest.raises(rk.TerminalLeft):
            agent._read_key(_Keys(["\x04"]), {"1", "2", "\x1b"})

    def test_during_a_turn(self):
        agent, _engine, _err = _agent()
        with pytest.raises(rk.TerminalLeft):
            agent._on_key(rk.KeyEvent(rk.EOF), rk.KeyDecoder())

    def test_a_key_the_dialog_offers_still_answers_it(self):
        agent, _engine, _err = _agent()
        assert agent._read_key(_Keys(["y"]), {"y", "n", "\x1b"}) == "y"


# ---------------------------------------------------------------------------
# Leaving unwinds everything
# ---------------------------------------------------------------------------

class TestLeavingUnwinds:
    def test_the_question_in_flight_is_refused(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
        broker = tc.TerminalConfirmBroker(session_key="t")
        agent, _engine, _err = _agent(broker=broker)
        req = broker._enqueue(tc.ConfirmRequest(
            kind=tc.CONFIRM, tool="bash", args={"command": "ls"},
            preview="x"))
        taken = broker.take()
        with pytest.raises(rk.TerminalLeft):
            agent._answer(taken, _Keys(["\x04"]))
        assert req.resolved is True
        assert req.decision is False

    def test_the_turn_is_stopped_and_not_waited_out(self, monkeypatch):
        release = threading.Event()

        def _a_tool_call_that_does_not_end(engine, prompt, *, sink,
                                           max_tokens=0):
            release.wait(30)
            sink(repl.RenderItem("done"))
            return repl.TurnResult()

        monkeypatch.setattr(repl, "run_turn", _a_tool_call_that_does_not_end)
        agent, engine, _err = _agent()

        def _hangup(_worker):
            raise rk.TerminalLeft("hangup", 129)

        agent._pump = _hangup
        t0 = time.monotonic()
        try:
            with pytest.raises(rk.TerminalLeft):
                agent.turn("go")
            assert time.monotonic() - t0 < 5
            assert engine.stopped
        finally:
            release.set()

    def test_nothing_asked_after_leaving_waits(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
        broker = tc.TerminalConfirmBroker(session_key="t")
        agent, _engine, _err = _agent(broker=broker)

        def _gone(_prompt):
            raise rk.TerminalLeft("hangup", 129)

        agent._read_line = _gone
        monkeypatch.setattr(agent, "_arm_leave_deadline", lambda code: None)
        agent.run()
        assert broker.aborted is True

    def test_the_loop_returns_through_its_cleanup(self, monkeypatch):
        agent, _engine, _err = _agent()

        def _gone(_prompt):
            raise rk.TerminalLeft("hangup", 129)

        agent._read_line = _gone
        withdrawn: list[int] = []
        armed: list[int] = []
        monkeypatch.setattr(agent, "_withdraw_presence",
                            lambda: withdrawn.append(1))
        monkeypatch.setattr(agent, "_arm_leave_deadline", armed.append)
        assert agent.run() == 129
        assert withdrawn == [1]
        assert armed == [129]


# ---------------------------------------------------------------------------
# SIGHUP and SIGTERM leave through the same door
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,code", [("SIGHUP", 129), ("SIGTERM", 143)])
def test_the_signal_leaves_through_the_same_door(name, code):
    sig = getattr(signal, name, None)
    if sig is None:
        pytest.skip(f"{name} does not exist here")
    agent, _engine, _err = _agent()
    before = signal.getsignal(sig)
    agent._install_sigint()
    try:
        handler = signal.getsignal(sig)
        assert callable(handler) and handler != before
        with pytest.raises(rk.TerminalLeft) as left:
            handler(sig, None)
        assert left.value.code == code
    finally:
        agent._restore_sigint()
        agent._restore_sigwinch()
    assert signal.getsignal(sig) == before


# ---------------------------------------------------------------------------
# A leave that hangs is cut short
# ---------------------------------------------------------------------------

def test_a_leave_that_hangs_is_cut_short(monkeypatch):
    ended = threading.Event()
    codes: list[int] = []

    def _exit(code):
        codes.append(code)
        ended.set()

    monkeypatch.setattr(repl, "_LEAVE_DEADLINE_S", 0.05)
    monkeypatch.setattr(os, "_exit", _exit)
    agent, _engine, _err = _agent()
    agent._arm_leave_deadline(129)
    assert ended.wait(5)
    assert codes == [129]
