"""The framed prompt is drawn by escapes. So it is judged on a screen.

The box was built and merged with no test that ever ran ``_draw``: the
pure renderer had eleven cases, the wiring had eleven more, and between
them nothing applied a single escape sequence to anything. Three faults
lived in that gap, and every one of them destroys the user's transcript:

  the first paint climbed   ``\\x1b[{rows-1}A`` ran unconditionally, so
                            the very first box was painted over the last
                            rows above it — the banner at start-up, the
                            end of the answer after a turn

  the redraw mixed two      it walked DOWN by the old box's height and
  geometries                UP by the new one's, so a line that wrapped
                            ate one more transcript row and a line that
                            unwrapped left a stray border

  the erase was sized       ``_clear_box`` built its view from the
  from the wrong box        decoder, which had already dropped its
                            buffer on submit — so it erased an empty
                            box's worth while a taller one was on screen

What follows is a model terminal: it applies exactly the escapes the
drawing emits (CSI A/B/C/K, CR, LF) to a grid of rows, and the tests
read the grid. That is the instrument the work was missing — the claims
above are its output, not an argument about the code.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import repl as R
from delfin.agent import repl_keys as rk


# --- a terminal, as far as these escapes are concerned ---------------------

class Screen:
    """Rows of text, a cursor, and the six operations the box uses."""

    def __init__(self, width: int = 40, transcript=()):
        self.width = width
        self.rows: list[str] = list(transcript)
        self.row = len(self.rows) - 1 if self.rows else 0
        self.col = 0
        if not self.rows:
            self.rows = [""]

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
                n = int(data[i + 2:j] or 1)
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
                self.rows[self.row] = line[:self.col] + ch + line[self.col + 1:]
                self.col += 1
            i += 1

    def flush(self):
        pass

    def isatty(self):
        return True

    def text(self) -> list[str]:
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
    """Run read_boxed against a scripted keyboard and a model screen."""
    def _run(keys, *, width=40, transcript=()):
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

        monkeypatch.setattr(rk, "RawMode", lambda *a, **k: _Keys(keys))
        text = agent.read_boxed()
        return text, screen
    return _run


TRANSCRIPT = ["answer line one", "answer line two", "answer line three",
              "answer line four", ""]


# --- what the missing instrument shows -------------------------------------

def test_the_first_box_does_not_paint_over_the_transcript(boxed):
    text, screen = boxed(["hi", "\r"], transcript=list(TRANSCRIPT))
    assert text == "hi"
    shown = screen.text()
    for line in TRANSCRIPT[:4]:
        assert line in shown, (
            f"the box painted over {line!r}:\n" + "\n".join(shown))


def test_a_line_that_wraps_does_not_eat_another_row(boxed):
    long = "x" * 80                       # two content rows at width 40
    text, screen = boxed([long, "\r"], width=40, transcript=list(TRANSCRIPT))
    assert text == long
    shown = screen.text()
    for line in TRANSCRIPT[:4]:
        assert line in shown, (
            f"wrapping ate {line!r}:\n" + "\n".join(shown))


def test_a_line_that_unwraps_leaves_no_stray_border(boxed):
    # grow to two content rows, then delete back to one
    keys = ["y" * 80] + ["\x7f"] * 60 + ["\r"]
    _text, screen = boxed(keys, width=40, transcript=list(TRANSCRIPT))
    shown = [r for r in screen.text() if r.strip()]
    tops = [r for r in shown if r.startswith("╭")]
    assert len(tops) <= 1, "a border of the old box was left:\n" + "\n".join(shown)


def test_a_submitted_wrapped_line_leaves_nothing_behind(boxed):
    long = "z" * 80
    _text, screen = boxed([long, "\r"], width=40, transcript=list(TRANSCRIPT))
    shown = [r for r in screen.text() if r.strip()]
    assert not [r for r in shown if r.startswith(("╭", "╰", "│"))], (
        "the box outlived its submit:\n" + "\n".join(shown))


def test_the_transcript_is_untouched_after_a_submit(boxed):
    _text, screen = boxed(["done", "\r"], transcript=list(TRANSCRIPT))
    shown = screen.text()
    assert shown[:4] == TRANSCRIPT[:4], "\n".join(shown)
