"""What the terminal layer promises when it takes something over.

Every item here is a case where the screen or the tty was left in a state
the code did not admit to:

THE TTY IS GIVEN BACK PER PROMPT. An approval raised outside a turn needs
cbreak for one keystroke. Entering raw mode without ever leaving it puts
the next idle ``input()`` in cbreak with echo off, and stacks one atexit
hook per approval for the life of the process.

A JOIN THAT CANNOT END IS NOT AN ESCAPE. The last step of the interrupt
ladder unwinds into the settle, so an unbounded join there waits on the
very tool call the user is escaping.

ONE ROW STAYS ONE ROW. The bottom row is erased with a single ``\\r\\x1b[K``,
which reaches exactly one screen line — so anything written there has to
fit in one, typed text included.

A CHUNK BOUNDARY IS NOT A KEYPRESS. A read returns at most 1024 bytes, so
a paste can be cut immediately after an ESC byte; reading that as Escape
ends the turn on nothing but a buffer size.

THE SCREEN STATES EFFECTS THAT HAPPENED. An answer to a request that
already expired is discarded by the broker, so the line printed after it
has to follow what the broker returned.

A REFUSAL HAS TO BE A KEY THE CALLER OFFERED. With no terminal to read
from, the fallback keystroke goes back into the caller's own option chain.
"""

from __future__ import annotations

import io
import os
import threading
import time

import pytest

from delfin.agent import repl, repl_keys as rk, repl_render as rr
from delfin.agent import terminal_confirm as tc


PLAIN = rr.Theme(enabled=False)


class _Tty(io.StringIO):
    """A stream that claims to be a terminal, so cursor control is allowed."""

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
    """Stands in for the raw-mode reader: hands over scripted keystrokes."""

    active = True

    def __init__(self, keys):
        self._keys = list(keys)

    def read_ready(self, timeout):
        return self._keys.pop(0) if self._keys else "n"


def _agent(engine=None, broker=None):
    engine = engine if engine is not None else _Engine()
    out, err = io.StringIO(), _Tty()
    agent = repl.TerminalAgent(engine, out=out, err=err, broker=broker)
    agent.transcript.theme = PLAIN
    return agent, engine, err


@pytest.fixture
def tty_stdin():
    """A real pty, because raw mode is the thing under test here."""
    pytest.importorskip("termios")
    pty = pytest.importorskip("pty")
    master, slave = pty.openpty()
    stream = os.fdopen(slave, "rb", buffering=0)
    try:
        yield stream
    finally:
        stream.close()
        os.close(master)


# ---------------------------------------------------------------------------
# The tty is given back per prompt
# ---------------------------------------------------------------------------

def test_an_approval_outside_a_turn_gives_the_terminal_back(
        tty_stdin, monkeypatch):
    """`!cmd` can raise a prompt, and the reader for it has to be released.

    Two failures if it is not: the next idle `input()` runs in cbreak with
    echo off, and RawMode.__enter__ registers an atexit hook that only
    restore() takes back — so one hook accumulates per approval.
    """
    created: list = []

    class _Recorded(rk.RawMode):
        def __init__(self, stream=None):
            super().__init__(stream)
            created.append(self)

    class _Hooks:
        """Counts LIVE hooks. atexit._ncallbacks() counts emptied slots too."""

        def __init__(self):
            self.live: list = []

        def register(self, func):
            self.live.append(func)
            return func

        def unregister(self, func):
            self.live[:] = [f for f in self.live if f != func]

    hooks = _Hooks()
    monkeypatch.setattr(rk, "RawMode", _Recorded)
    monkeypatch.setattr(rk, "atexit", hooks)

    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker=broker)
    agent._stdin = tty_stdin

    entered: list[bool] = []

    def _answer(req, raw):
        entered.append(bool(getattr(raw, "active", False)))
        broker.resolve(req, True)

    agent._answer = _answer

    def _work():
        return [broker.callback("bash", {"command": "ls"}, "$ ls"),
                broker.callback("bash", {"command": "id"}, "$ id")]

    assert agent._off_thread(_work) == [True, True]

    assert entered == [True, True], (
        "a keystroke prompt needs the terminal out of canonical mode")
    assert created, "the prompt reader is built through repl_keys.RawMode"
    assert not any(raw.active for raw in created), (
        "a reader left entered leaves the next idle input() without echo")
    assert hooks.live == [], (
        "every __enter__ registers a hook that only restore() takes back, so "
        "an unreleased reader leaves one behind per approval")
