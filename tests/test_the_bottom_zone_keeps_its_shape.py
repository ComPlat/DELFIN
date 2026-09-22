"""While a turn runs, the place you type keeps its shape.

The idle prompt draws two rules with the text between them. During a
turn that whole area used to vanish: one bare row carried the spinner or
what was being typed, and the input area was simply not there for the
length of the turn — which is the half of a session a user most wants to
know they can still reach.

The bottom zone is two rows now, a rule and the line under it, in the
same vocabulary as the prompt above. What that costs is care with the
erase: two rows to clear, and the cursor has to come back to where the
rule began, or a transcript line lands under a rule nobody removed and
the box starts eating the conversation. That fault is not hypothetical —
it happened to the framed prompt in its first week, three ways.

So it is judged on a screen: a model terminal applying the escapes the
painter emits (CSI A/K, CR, LF) to a grid of rows.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl as R


class Screen:
    """Rows of text, a cursor, and the operations the painter uses."""

    def __init__(self, width: int = 40, transcript=()):
        self.width = width
        self.rows: list[str] = list(transcript) or [""]
        self.row = len(self.rows) - 1
        self.col = 0

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

    def text(self):
        return [r.rstrip() for r in self.rows]


TRANSCRIPT = ["answer one", "answer two", "answer three", ""]


@pytest.fixture()
def agent():
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    screen = Screen(40, list(TRANSCRIPT))
    a.err = screen
    a._flush_err = screen.flush
    a._bottom = ""
    a._input_line = ""

    class _T:
        width = 40
    a.transcript = _T()
    a._can_redraw = lambda: True
    a.screen = screen
    return a


def test_the_zone_is_a_rule_and_a_line(agent):
    agent._set_bottom("⠴ 12s  esc to interrupt")
    shown = [r for r in agent.screen.text() if r.strip()]
    assert shown[-1].startswith("⠴"), shown
    assert set(shown[-2]) == {"─"}, "a rule sits above the line"


def test_it_does_not_paint_over_the_transcript(agent):
    agent._set_bottom("⠴ 12s")
    shown = agent.screen.text()
    for line in TRANSCRIPT[:3]:
        assert line in shown, f"{line!r} was painted over:\n" + "\n".join(shown)


def test_clearing_takes_both_rows(agent):
    agent._set_bottom("⠴ 12s")
    agent._clear_bottom()
    left = [r for r in agent.screen.text() if r.strip()]
    assert not [r for r in left if set(r) == {"─"}], (
        "the rule outlived the line it framed:\n" + "\n".join(left))
    assert "⠴ 12s" not in "\n".join(left)


def test_the_transcript_is_whole_after_a_repaint_cycle(agent):
    """Set, clear, set again — the shape a streaming turn draws many
    times a second. Nothing above may move."""
    for _ in range(5):
        agent._set_bottom("⠴ working")
        agent._clear_bottom()
    agent._set_bottom("⠴ working")
    shown = agent.screen.text()
    assert shown[:3] == TRANSCRIPT[:3], "\n".join(shown)


def test_clearing_nothing_writes_nothing(agent):
    before = list(agent.screen.text())
    agent._clear_bottom()
    assert agent.screen.text() == before
