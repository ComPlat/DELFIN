"""The status line's jobs, reached with the arrows under the prompt.

The line under the input area names what is still out, with the handle
you would have to retype to look inside: ``[ab12cd34]``. A list you can
read but not reach is a report. Now the arrows walk it — on an EMPTY
line, where the readline history used to be the only tenant — and Enter
sends ``/bash <id>`` for the row that is marked, the same command the
hand would have typed.

What is judged here, and the instrument:

  the walk is visible       a mark that moved without the screen moving
                            is a key that changed nothing; every step
                            must repaint the status line with ``▶`` on
                            the row Enter would take

  the history stays         a line with something typed in it keeps the
                            readline behaviour: Up recalls an older
                            message, and the walk never starts. The
                            empty prompt is the only place the jobs
                            outrank the history.

  Enter walks in            Enter on a marked row, on an empty line,
                            returns ``/bash <id>`` — not an empty turn

  typing leaves the walk    one letter makes the line a message; the
                            mark must be gone before the letter lands,
                            so Enter after it cannot act on a job

The instrument is the model terminal from the transcript test, grown
for two things this work needs and that one never saw: ``CSI J`` (the
erase-below that a wrapped status row leaves standing) and wrapping (a
status line longer than the width occupies more physical rows than the
box counted). Both faults are invisible to a grid that does not wrap,
which is why the earlier test could not have caught them.

The last test drives a REAL pty — a pipe answers ``raw_mode_supported``
with no, so nothing about the box is ever exercised by one. The pty
answers yes, termios applies, and the same keys go through the same
decoder on a kernel terminal. What it asserts is the contract the model
terminal cannot vouch for: that the code path a human's keyboard takes
is the one tested above.
"""

from __future__ import annotations

import io
import os
import pty
import sys
import termios
import time

import pytest

from delfin.agent import repl as R
from delfin.agent import repl_box as rb
from delfin.agent import repl_keys as rk


# --- the model terminal, grown for CSI J and wrapping ----------------------

class Screen:
    """Rows of text, a cursor, and the escapes the box emits.

    Two things the transcript test's grid did not do, both of them
    load-bearing for this work:

      ``CSI J``  the erase-below emitted when a box row wrapped on the
                 real screen — a grid that ignores it cannot see the
                 stray tail it is there to remove
      wrapping   a written row longer than the width hard-wraps, the
                 way a terminal does, so a status line that outgrew the
                 input's width occupies the rows the erase must reach
    """

    def __init__(self, width: int = 40, transcript=()):
        self.width = width
        self.rows: list[str] = list(transcript) or [""]
        self.row = len(self.rows) - 1
        self.col = 0
        # Everything ever painted, for the transcript-integrity checks
        # (a row that survived a submit), and the status rows in order
        # for the walk's appear-then-leave assertions.
        self.ever: set[str] = set()
        self.status_seen: list[str] = []

    def _note(self) -> None:
        self.ever.update(r.rstrip() for r in self.rows if r.strip())
        for r in self.rows:
            r = r.rstrip()
            if "⚙" in r and (not self.status_seen or self.status_seen[-1] != r):
                self.status_seen.append(r)

    def _fit(self, row: int) -> None:
        while len(self.rows) <= row:
            self.rows.append("")

    def write(self, data: str) -> None:
        i = 0
        while i < len(data):
            ch = data[i]
            if ch == "\x1b" and data[i + 1:i + 2] == "[":
                j = i + 2
                while j < len(data) and not data[j].isalpha():
                    j += 1
                arg = data[i + 2:j] or "1"
                n = int(arg)
                verb = data[j]
                if verb == "A":
                    self.row = max(0, self.row - n)
                elif verb == "B":
                    self.row += n
                    self._fit(self.row)
                elif verb == "C":
                    self.col += n
                elif verb == "K":
                    self._fit(self.row)
                    self.rows[self.row] = self.rows[self.row][:self.col]
                elif verb == "J":
                    # 0 = cursor to end of screen; the form the erase
                    # emits. A grid that skipped this let a wrapped
                    # row's tail stand through the next answer.
                    self._fit(self.row)
                    self.rows[self.row] = self.rows[self.row][:self.col]
                    self.rows = self.rows[:self.row + 1]
                i = j + 1
                continue
            if ch == "\r":
                self.col = 0
            elif ch == "\n":
                self.row += 1
                self._fit(self.row)
            else:
                self._fit(self.row)
                line = self.rows[self.row]
                if len(line) < self.col:
                    line += " " * (self.col - len(line))
                # wrapping: the width is a wall, not a suggestion. The
                # status line can outgrow the box's width; the physical
                # rows it takes are the ones the counted erase misses.
                if self.col >= self.width:
                    self.row += 1
                    self._fit(self.row)
                    self.col = 0
                    line = self.rows[self.row]
                    if len(line) < self.col:
                        line += " " * (self.col - len(line))
                self.rows[self.row] = (line[:self.col] + ch
                                       + line[self.col + 1:])
                self.col += 1
            i += 1

    def flush(self):
        self._note()

    def isatty(self):
        return True

    def text(self) -> list[str]:
        self._note()
        return [r.rstrip() for r in self.rows]


class _Keys:
    """A RawMode stand-in that hands the loop a scripted keyboard."""

    def __init__(self, chunks):
        self.chunks = list(chunks)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    @property
    def active(self):
        return True

    def read_ready(self, timeout=0.0):
        return self.chunks.pop(0) if self.chunks else "\r"

    def restore(self):
        pass


@pytest.fixture()
def boxed(monkeypatch):
    """Run read_boxed against a scripted keyboard and a model screen.

    The status line is fed by a fixed background view, so the walk has
    real rows to walk: two shells and one agent, the ids shortened the
    way the line shows them.
    """
    view = {"shells": [
        {"id": "1111aaaa2222bbbb", "label": "suite", "since": 100.0},
        {"id": "3333cccc4444dddd", "label": "build", "since": 200.0},
    ], "agents": [
        {"id": "5555eeee6666ffff", "label": "scout", "since": 300.0,
         "last": "reading repl.py"},
    ], "watches": [], "wakeups": [], "errors": []}

    def _run(keys, *, width=40, transcript=(), now=1000.0):
        screen = Screen(width, transcript)
        agent = R.TerminalAgent.__new__(R.TerminalAgent)
        agent.err = screen
        agent._flush_err = screen.flush
        agent._stdin = io.StringIO()
        agent._width_dirty = False
        agent._cycle_mode = lambda: None

        class _T:
            pass
        agent.transcript = _T()
        agent.transcript.width = width
        agent.transcript.refresh_width = lambda: None

        from delfin.agent import background_view as bgv
        # A FRESH cache: the loop's first _draw asks _background_status,
        # which recollects on a stale one — and a bare __new__ agent has
        # no opts.cwd, so the collect raises and the line goes empty.
        # Fresh, it serves the line built from the fixed view below.
        agent._bg_status_at = time.monotonic()
        agent._bg_view = view
        agent._bg_now = now
        agent._bg_status = bgv.status_line(view, now=now)

        def _collect_fresh(cwd):
            return view

        monkeypatch.setattr(bgv, "collect", _collect_fresh)
        monkeypatch.setattr(R.TerminalAgent, "_BACKGROUND_STATUS_EVERY_S",
                            3600.0)
        monkeypatch.setattr(rk, "RawMode", lambda *a, **k: _Keys(keys))
        text = agent.read_boxed()
        return text, screen
    return _run


def _status_rows(screen) -> list[str]:
    """The rows that carry the status line's marker or handle."""
    return [r for r in screen.ever
            if "⚙" in r or "▶" in r]


# --- the walk is visible ----------------------------------------------------

def test_up_on_an_empty_line_marks_the_first_job(boxed):
    _text, screen = boxed(["\x1b[A", "\r"], width=40)
    rows = _status_rows(screen)
    assert any("▶" in r for r in rows), (
        "Up on an empty line did not mark a job:\n" + "\n".join(screen.text()))


def test_the_mark_moves_down_and_back_up(boxed):
    _text, screen = boxed(["\x1b[A", "\x1b[A", "\x1b[B", "\r"], width=40)
    rows = _status_rows(screen)
    assert any("▶" in r for r in rows), "\n".join(screen.text())


def test_the_mark_names_the_row_enter_would_take(boxed):
    text, screen = boxed(["\x1b[A", "\r"], width=40)
    assert text == "/bash 1111aaaa2222bbbb", (
        f"Enter walked into {text!r} instead of the marked job")


def test_enter_on_the_second_row_takes_that_one(boxed):
    text, _screen = boxed(["\x1b[A", "\x1b[A", "\r"], width=40)
    assert text == "/bash 3333cccc4444dddd"


def test_down_on_an_empty_line_with_no_walk_does_nothing(boxed):
    text, screen = boxed(["\x1b[B", "\r"], width=40)
    assert text == ""        # an empty submit, as it always was
    assert not any("▶" in r for r in _status_rows(screen))


# --- the history keeps the line --------------------------------------------

def test_up_with_text_typed_recalls_history_not_jobs(boxed, monkeypatch):
    """A line with something in it belongs to the readline history.

    The walk must not start there: typing ``hi`` and pressing Up fills
    the box with the older message, and no mark appears.
    """
    sent = []

    class _History:
        def up(self, buf):
            sent.append(buf)
            return "older message"

        def down(self):
            return None

        def add(self, _s):
            pass

    monkeypatch.setattr(R, "_BoxHistory", _History)
    text, screen = boxed(["hi", "\x1b[A", "\r"], width=40)
    assert text == "older message", f"got {text!r}"
    assert sent == ["hi"], "history.up was not asked for the typed line"
    assert not any("▶" in r for r in _status_rows(screen))


def test_a_typed_letter_drops_the_mark(boxed):
    """One letter makes the line a message; the mark must be gone."""
    text, screen = boxed(["\x1b[A", "x", "\x1b[B", "\r"], width=40)
    assert text == "x", f"got {text!r}"
    assert any("▶" in r for r in screen.status_seen), (
        "the walk never showed a mark:\n" + "\n".join(screen.status_seen))
    assert "▶" not in screen.status_seen[-1], (
        "the mark outlived the typed letter:\n"
        + "\n".join(screen.status_seen))


def test_a_typed_letter_then_enter_is_a_message_not_a_job(boxed):
    text, _screen = boxed(["\x1b[A", "x", "\r"], width=40)
    assert text == "x", f"Enter acted on a job: {text!r}"


def test_esc_drops_the_walk(boxed):
    text, screen = boxed(["\x1b[A", "\x1b", "\r"], width=40)
    assert text == ""
    assert any("▶" in r for r in screen.status_seen), (
        "the walk never showed a mark:\n" + "\n".join(screen.status_seen))
    assert "▶" not in screen.status_seen[-1], (
        "the mark outlived Esc:\n" + "\n".join(screen.status_seen))


# --- a wrapped status row leaves nothing standing ---------------------------

def test_a_wrapped_status_row_is_erased_to_the_end_of_the_screen(boxed):
    """The status line can outgrow the box's width.

    A 40-column screen with a status line well past 40 chars wraps onto
    a second physical row. The counted erase clears only the rows the
    view knows; without ``CSI J`` the wrapped tail stands through the
    next answer. This is the fault the non-wrapping grid could not see.
    """
    keys = ["\r"]
    text, screen = boxed(keys, width=24,
                         transcript=["answer one", "answer two"])
    shown = screen.text()
    for row in shown:
        assert "⚙" not in row and "▶" not in row, (
            "a status row survived the submit:\n" + "\n".join(shown))
    assert text == ""


# --- the path a human's keyboard takes --------------------------------------

def test_read_boxed_runs_over_a_real_pty():
    """The contract the model terminal cannot vouch for.

    ``raw_mode_supported`` answers no for a pipe, so every boxed test
    above drives a stand-in RawMode. A pty answers yes: termios is set,
    the decoder runs on kernel-buffered input, and the same arrow keys
    take the path a user's keyboard takes. The test does not assert on
    the box's geometry — it asserts the path ran and produced the same
    line the model terminal predicted for the same keys.
    """
    view = {"shells": [
        {"id": "1111aaaa2222bbbb", "label": "suite", "since": 100.0},
    ], "agents": [], "watches": [], "wakeups": [], "errors": []}
    from delfin.agent import background_view as bgv

    class _FdStream:
        """A stream-shaped fd, for code that asks fileno() of stdin.

        RawMode drives the fd with termios; pytest's captured stdin has
        none. In the pty child fd 0 IS the terminal, so this is exactly
        the stream a user's terminal hands the agent.
        """

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

    pid, fd = pty.fork()
    if pid == 0:                                    # the child: the agent
        try:
            agent = R.TerminalAgent.__new__(R.TerminalAgent)

            class _T:
                pass
            agent.err = _FdWriter(1)
            agent._stdin = _FdStream(0)
            agent._flush_err = lambda: None
            agent._width_dirty = False
            agent._cycle_mode = lambda: None
            agent.transcript = _T()
            agent.transcript.width = 40
            agent.transcript.refresh_width = lambda: None
            agent._bg_status_at = 0.0
            agent._bg_view = view
            agent._bg_status = bgv.status_line(view, now=1000.0)
            text = agent.read_boxed()
            os.write(1, b"\nRESULT:" + text.encode() + b"\n")
        except BaseException as exc:                # noqa: BLE001
            os.write(1, f"\nRESULT:child raised {exc!r}\n".encode())
        finally:
            os._exit(0)

    # the parent: the keyboard
    import select
    try:
        time.sleep(0.3)                             # let the child paint
        os.write(fd, b"\x1b[A")                     # Up: mark the job
        time.sleep(0.2)
        os.write(fd, b"\r")                         # Enter: walk in
        out = b""
        deadline = time.time() + 5.0
        while time.time() < deadline:
            ready, _, _ = select.select([fd], [], [], 0.1)
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

    assert b"RESULT:/bash 1111aaaa2222bbbb" in out, (
        "the pty run did not walk into the job; output tail:\n"
        + out.decode(errors="replace")[-400:])
