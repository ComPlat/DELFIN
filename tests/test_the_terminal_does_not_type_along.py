"""While a turn runs, one thing draws what is being typed. Not two.

The key loop reads keystrokes itself and paints them on the bottom row,
which is only correct while the terminal is NOT also echoing them. That
property is currently a side effect of ``tty.setcbreak``, and the module
docstring states it as a fact — so it is worth holding to, because
nothing else in the suite would notice it going away. cbreak's echo
behaviour has already been argued over once upstream (3.12 rewrote the
function around ``cfmakecbreak``), and DELFIN runs on whichever
interpreter a cluster happens to offer, not only on the suite's 3.11.

Asserted on a real pty, since this is a question about a terminal: the
flag itself, before, during and after — not the call that implies it.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import repl_keys as rk

termios = pytest.importorskip("termios")
pty = pytest.importorskip("pty")


@pytest.fixture()
def terminal():
    master, slave = pty.openpty()
    stream = os.fdopen(slave, "r", closefd=False)
    try:
        yield stream, slave
    finally:
        stream.close()
        for fd in (slave, master):
            try:
                os.close(fd)
            except OSError:
                pass


def _echo_on(fd: int) -> bool:
    return bool(termios.tcgetattr(fd)[3] & termios.ECHO)


def test_the_terminal_stops_echoing_for_the_turn(terminal):
    stream, fd = terminal
    assert _echo_on(fd), "the pty starts out echoing, as terminals do"

    with rk.RawMode(stream) as raw:
        assert raw.active, "the mode was entered at all"
        assert not _echo_on(fd), (
            "the terminal would echo keystrokes while the key loop draws "
            "them too -- every character twice, the second copy wherever "
            "the cursor happens to stand")


def test_it_is_given_back_afterwards(terminal):
    stream, fd = terminal
    with rk.RawMode(stream):
        pass
    assert _echo_on(fd), "a terminal left without echo ruins the session"


def test_the_prompt_reader_gives_it_back_too(terminal):
    """The approval prompt enters the same mode by hand and releases it
    through ``restore``; both paths matter, and the pairing is what keeps
    the next idle prompt visible."""
    stream, fd = terminal
    raw = rk.RawMode(stream)
    raw.__enter__()
    assert not _echo_on(fd)
    raw.restore()
    assert _echo_on(fd)


def test_the_interrupt_still_arrives_as_a_signal(terminal):
    """cbreak rather than raw, so Ctrl+C stays a signal and the interrupt
    ladder in repl.py keeps working. ISIG is the flag that says so."""
    stream, fd = terminal
    with rk.RawMode(stream):
        lflag = termios.tcgetattr(fd)[3]
        assert not (lflag & termios.ICANON), "keys arrive without Enter"
        assert lflag & termios.ISIG, "Ctrl+C must not become a plain byte"
