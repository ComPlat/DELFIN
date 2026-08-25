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
import signal
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


# ---------------------------------------------------------------------------
# The last step of the interrupt ladder
# ---------------------------------------------------------------------------

class _StuckEngine(_Engine):
    """A turn that does not come back until the test lets it."""

    def __init__(self, release: threading.Event):
        super().__init__()
        self._release = release
        self.entered = threading.Event()

    def stream_response(self, **kwargs):
        self.entered.set()
        self._release.wait(timeout=30)
        return "late answer"


def test_the_third_interrupt_does_not_wait_on_the_call_it_escapes():
    """It raises KeyboardInterrupt, which unwinds through turn()'s finally.

    That finally settles the worker, and an unbounded join there waits for
    the tool call the interrupt exists to walk away from.
    """
    release = threading.Event()
    engine = _StuckEngine(release)
    agent, _engine, _err = _agent(engine)

    def _interrupted_pump(worker):
        engine.entered.wait(timeout=5)
        raise KeyboardInterrupt

    agent._pump = _interrupted_pump
    raised: list = []

    def _drive():
        try:
            agent.turn("the long one")
        except KeyboardInterrupt:
            raised.append("interrupt")
        except BaseException as exc:                # noqa: BLE001
            raised.append(exc)

    driver = threading.Thread(target=_drive, daemon=True)
    began = time.monotonic()
    driver.start()
    driver.join(timeout=5)
    took = time.monotonic() - began
    stuck = driver.is_alive()
    release.set()
    driver.join(timeout=10)

    assert not stuck, (
        "the settle joined with no timeout, so leaving waited on exactly the "
        "call being left")
    assert took < 4, f"the escape took {took:.1f}s"
    assert raised == ["interrupt"]


def test_the_ladder_promises_only_what_the_third_interrupt_does():
    """The step before it describes the step after it, so it has to be true."""
    agent, _engine, _err = _agent()
    agent._turn_active.set()
    agent._on_sigint(signal.SIGINT, None)
    agent._on_sigint(signal.SIGINT, None)

    said = []
    while not agent._q.empty():
        said.append(agent._q.get_nowait().text)
    text = " ".join(said).lower()

    assert "leaves this session" in text, (
        "what it does is return 130 from the loop, not kill the process")
    assert "abandons the process" not in text


# ---------------------------------------------------------------------------
# One row stays one row
# ---------------------------------------------------------------------------

def test_the_typed_row_is_cut_to_one_screen_line():
    """Past COLUMNS the terminal wraps and the row becomes two.

    `_clear_bottom` erases with one `\\r\\x1b[K`, which reaches the last of
    them — so the first half stays in the transcript, in the middle of
    whatever is printed next. The status row is already cut for this
    reason; the typed row was not.
    """
    agent, _engine, err = _agent()
    agent.transcript.width = 40
    agent._draw_input_line("x" * 200 + "the end")

    row = err.getvalue().rsplit("\x1b[K", 1)[-1]
    assert len(row) <= 39, (
        f"the typed row is {len(row)} columns wide and wraps at 40")
    assert row.endswith("the end"), (
        "the END has to survive the cut: that is where the cursor is")


def test_a_typed_row_that_fits_is_left_alone():
    agent, _engine, err = _agent()
    agent.transcript.width = 40
    agent._draw_input_line("stop and check the tests")
    assert err.getvalue().rsplit("\x1b[K", 1)[-1] == "» stop and check the tests"


# ---------------------------------------------------------------------------
# A chunk boundary is not a keypress
# ---------------------------------------------------------------------------

def _interrupts(events) -> list:
    return [e for e in events if e.kind == rk.INTERRUPT]


def test_a_paste_cut_after_an_escape_does_not_end_the_turn():
    """read_ready takes at most 1024 bytes, and a paste is longer than that.

    Cut exactly after an ESC byte, the sequence's own head becomes the last
    byte of the chunk — and reading that as Escape ends the turn because of
    a buffer size and nothing else.
    """
    decoder = rk.KeyDecoder()
    head = "paste " + "x" * 1017 + "\x1b"
    assert len(head) == 1024, "the read size is what makes the cut land here"

    assert _interrupts(decoder.feed(head)) == [], (
        "an ESC at a read boundary may still be the head of a sequence")
    assert _interrupts(decoder.feed("[Bmore text")) == []
    assert decoder.buffer.endswith("more text"), (
        "the rest of the paste still has to reach the line")


def test_a_posture_key_split_across_two_reads_is_still_one_key():
    decoder = rk.KeyDecoder()
    assert _interrupts(decoder.feed("typed\x1b")) == []
    assert decoder.feed("[Z") == [rk.KeyEvent(rk.CYCLE_MODE)]


def test_an_escape_with_nothing_behind_it_still_ends_the_turn():
    """The behaviour the hold must not cost: Esc alone interrupts."""
    assert rk.KeyDecoder().feed("\x1b") == [rk.KeyEvent(rk.INTERRUPT)]

    decoder = rk.KeyDecoder()
    assert _interrupts(decoder.feed("stop\x1b")) == []
    assert decoder.feed("") == [rk.KeyEvent(rk.INTERRUPT)], (
        "the pump feeds every read, empty ones included — that is the "
        "deadline that says nothing followed the ESC")


# ---------------------------------------------------------------------------
# The screen states effects that happened
# ---------------------------------------------------------------------------

def _expired(req):
    """A request in the state a timeout leaves it in."""
    req.expired = True
    req.resolved = True
    req.decision = tc.TerminalConfirmBroker._refusal_for(req)
    return req


def test_a_keystroke_that_lands_after_the_request_expired_says_so():
    """`resolve` discards a late answer and returns False saying it did.

    Printing "refused" anyway credits this keystroke with an effect the
    broker did not apply: the refusal was the timeout's, and the model has
    already been told.
    """
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, err = _agent(broker=broker)

    req = _expired(tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="write_file", args={"path": "x.py"},
        preview="--- a/x.py\n+++ b/x.py\n+x = 1\n"))
    agent._answer(req, _Keys(["n"]))

    text = err.getvalue()
    assert "too late" in text
    assert "  refused" not in text, (
        "that refusal belonged to the timeout, not to this keystroke")


def test_a_late_plan_approval_does_not_announce_a_posture():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, err = _agent(broker=broker)

    req = _expired(tc.ConfirmRequest(kind=tc.PLAN, preview="1. read\n2. edit"))
    agent._answer(req, _Keys(["d"]))

    assert req.decision == {"approved": False, "new_mode": "plan"}, (
        "the expiry already answered it; the key must not overwrite that")
    text = err.getvalue()
    assert "approval → default" not in text
    assert "too late" in text
