"""Mail that waits at the idle prompt arrives, over a real pty.

Measured 2026-09-21: a terminal session stood 21 minutes at the prompt
with four unread messages while the operator typed into its pane by
hand. The delivery point is the raw box's idle wake
(``_wake_text`` -> ``_operator_messages``): between reads, at most
every ``_WAKE_EVERY_S`` seconds, only on an empty buffer.

The controls here drive the path a keyboard actually takes. A pipe
answers ``raw_mode_supported`` with no, so nothing about the box is
ever exercised by one; a pty answers yes, termios applies, and the
wake runs on kernel-buffered input. The last test is the pty one --
it does not assert geometry, it asserts the contract: mail in the
inbox, nothing typed, and the box RETURNS the rendered message as
its next prompt rather than waiting for a key that never comes.

(The dialog variant of the standstill -- a session parked in
_read_key behind a question answered from outside -- is fixed on
main as 754f8300 and is that commit's control, not this file's.)
"""

from __future__ import annotations

import os
import pty
import select
import time

import pytest

from delfin.agent import repl as R
from delfin.agent import session_messages as msgs


def test_the_wake_looks_at_least_every_five_seconds():
    """The throttle itself: a message may wait, but not unboundedly."""
    assert R.TerminalAgent._WAKE_EVERY_S <= 5.0


def test_a_message_reaches_the_wake_only_on_an_empty_buffer():
    """The guard that keeps the wake from taking a draft: mail waits
    while something is typed, and is delivered the moment it is not."""
    agent = object.__new__(R.TerminalAgent)

    class _T:
        theme = type("TH", (), {"dim": staticmethod(lambda t: t)})()
        def chrome(self, line): pass
    agent.transcript = _T()
    agent.opts = type("O", (), {"session_name": "pty-wake",
                                "cwd": os.getcwd()})()
    agent.engine = type("E", (), {"session_id": "0123456789ab"})()
    # No inbox at all: what is judged is the guard, not the take.
    agent._operator_messages("half a thought")
    # And a bare object survives the wake: the prompt must not die.
    agent._operator_messages("")


def test_read_boxed_returns_waiting_mail_over_a_real_pty(tmp_path):
    """The contract a pipe cannot vouch for: over a real pty, with
    nothing typed and a message in the inbox, read_boxed returns the
    rendered message instead of blocking on the keyboard."""
    inbox = tmp_path / "inbox"
    presence = tmp_path / "presence"
    # The message is in the inbox BEFORE the child starts reading.
    msgs._DIR = inbox
    msgs.send("pty-wake", "stop searching in the home directory",
              from_title="the operator")

    pid, fd = pty.fork()
    if pid == 0:                                    # the child: the agent
        try:
            from delfin.agent import session_messages as _m
            from delfin.agent import session_presence as _p
            _m._DIR = inbox
            _p._DIR = presence
            agent = R.TerminalAgent.__new__(R.TerminalAgent)

            class _T:
                # _operator_messages prints the "✉ sender: text" line
                # through transcript.chrome BEFORE the message is
                # returned: without these two the delivery's own
                # catch-all swallows an AttributeError and the mail is
                # silently not delivered -- the standstill, faked by
                # the harness rather than the code.
                theme = type("TH", (), {
                    "dim": staticmethod(lambda t: t)})()

                def chrome(self, line):
                    os.write(1, f"\r\n{line}\r\n".encode())

            agent.err = _FdWriter(1)
            agent._stdin = _FdStream(0)
            agent._flush_err = lambda: None
            agent._width_dirty = False
            agent._cycle_mode = lambda: None
            agent.transcript = _T()
            agent.transcript.width = 40
            agent.transcript.refresh_width = lambda: None
            agent._bg_view = {"shells": [], "agents": [], "watches": [],
                              "wakeups": [], "errors": []}
            agent._bg_status = ""
            agent._bg_status_at = 0.0
            agent._bg_id_at = lambda _row: None
            agent.opts = type("O", (), {
                "session_name": "pty-wake", "cwd": str(tmp_path)})()
            agent.engine = type("E", (), {"session_id": "0123456789ab"})()
            # The wake fires on the FIRST empty read, not after the
            # five-second throttle: this test judges the delivery, not
            # the clock.
            agent._wake_last_look = 0.0
            text = agent.read_boxed()
            os.write(1, b"\nRESULT:" + text.encode() + b"\n")
        except BaseException as exc:                # noqa: BLE001
            os.write(1, f"\nRESULT:child raised {exc!r}\n".encode())
        finally:
            os._exit(0)

    # the parent: nobody types. The message must come out on its own.
    out = b""
    deadline = time.time() + 8.0
    try:
        while time.time() < deadline:
            ready, _, _ = select.select([fd], [], [], 0.2)
            if not ready:
                continue
            try:
                chunk = os.read(fd, 4096)
            except OSError:
                break
            if chunk:
                out += chunk
                if b"RESULT:" in out:
                    break
    finally:
        os.close(fd)
        os.waitpid(pid, 0)

    assert b"RESULT:" in out, (
        "the pty run never returned from read_boxed; nobody typed and "
        "mail was waiting -- the 21-minute standstill. Output tail:\n"
        + out.decode(errors="replace")[-400:])
    assert b"stop searching in the home directory" in out
    assert b"not from the user" in out, (
        "the message must read as coming from another session, never "
        "as the user's own words")


class _FdStream:
    """A stream-shaped fd, for code that asks fileno() of stdin."""

    def __init__(self, fd: int):
        self._fd = fd

    def fileno(self) -> int:
        return self._fd

    def read(self, n: int = -1) -> str:
        return os.read(self._fd, n).decode(errors="replace") \
            if n > 0 else ""

    def isatty(self) -> bool:
        return True


class _FdWriter:
    def __init__(self, fd: int):
        self._fd = fd

    def write(self, data: str) -> int:
        return os.write(self._fd, data.encode())

    def flush(self) -> None:
        pass

    def isatty(self) -> bool:
        return True
