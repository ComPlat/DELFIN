"""The zone under the prompt must erase exactly what it drew.

Two faults reported from a real terminal on 2026-09-19, both from the
two-row bottom zone shipped that morning, and both invisible to the model
terminal that was written with it — because that model had no line
wrapping.

FIRST, during a long first turn the rules piled up:

    ──────────────────────────────────────────
    ! ⏳ First turn on kit.glm-5.3: the endpoint is building its cache …
    ──────────────────────────────────────────
    ──────────────────────────────────────────
    ──────────────────────────────────────────      (× 60)

``_set_bottom`` writes a rule, a newline and the text; ``_clear_bottom``
erases the current row, moves up one and erases that. Exactly two rows.
But a text wider than the terminal occupies two rows or three, so the
erase came up short and every repaint left one more rule standing.

SECOND, on the way out the prompt was left on screen and the shell wrote
into it:

    > (.venv) [user@host]$ ^C
    (.venv) [user@host]$ ^C────────────────────────────
    (.venv) [user@host]$  mode · /help

Nothing tore the zone down before the loop returned.

So the model terminal here WRAPS, which is the whole point of it, and the
zone counts the rows it actually put on the screen.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl as R


class WrappingScreen:
    """Rows, a cursor, and — unlike the last one — a right edge.

    A write that reaches the last column moves to the next row, the way a
    terminal in its default (autowrap) mode does.
    """

    def __init__(self, width: int = 40, rows=()):
        self.width = width
        self.rows: list[str] = list(rows) or [""]
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
                    self.col = min(self.width - 1, self.col + n)
                elif verb == "K":
                    self._fit(self.row)
                    self.rows[self.row] = self.rows[self.row][:self.col]
                elif verb == "J":
                    # Erase from the cursor to the end of the screen.
                    self._fit(self.row)
                    self.rows[self.row] = self.rows[self.row][:self.col]
                    del self.rows[self.row + 1:]
                i = j + 1
                continue
            if ch == "\r":
                self.col = 0
            elif ch == "\n":
                self.row += 1
                self.col = 0
                self._fit(self.row)
            else:
                self._fit(self.row)
                line = self.rows[self.row]
                if len(line) < self.col:
                    line += " " * (self.col - len(line))
                self.rows[self.row] = line[:self.col] + ch + line[self.col + 1:]
                self.col += 1
                if self.col >= self.width:       # the right edge
                    self.col = 0
                    self.row += 1
                    self._fit(self.row)
            i += 1

    def flush(self):
        pass

    def isatty(self):
        return True

    def text(self):
        return [r.rstrip() for r in self.rows]

    def rules(self):
        return [r for r in self.text() if r and set(r) == {"─"}]


TRANSCRIPT = ["answer one", "answer two", ""]

LONG = ("! First turn on this model: the endpoint is building its cache "
        "for this prompt — usually about 200 s here, and on a busy day "
        "its queue adds more.")


@pytest.fixture()
def agent():
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    screen = WrappingScreen(40, list(TRANSCRIPT))
    a.err = screen
    a._flush_err = screen.flush
    a._bottom = ""
    a._bottom_rows = 0
    a._input_line = ""

    class _T:
        width = 40
    a.transcript = _T()
    a._can_redraw = lambda: True
    a.screen = screen
    return a


def test_a_text_wider_than_the_screen_leaves_one_rule(agent):
    """The reported fault, in its smallest form."""
    for _ in range(6):
        agent._set_bottom(LONG)
        agent._clear_bottom()
    agent._set_bottom("⠴ 12s")
    assert len(agent.screen.rules()) == 1, (
        f"{len(agent.screen.rules())} rules left standing:\n"
        + "\n".join(agent.screen.text()))


def test_a_wrapping_text_does_not_eat_the_transcript(agent):
    agent._set_bottom(LONG)
    agent._clear_bottom()
    shown = agent.screen.text()
    for line in TRANSCRIPT[:2]:
        assert line in shown, "\n".join(shown)


def test_the_zone_is_one_rule_and_one_line_however_long_the_text(agent):
    """A bottom line that wraps is wrong on its own: the place you type
    must keep its shape. Long text is cut, not folded."""
    agent._set_bottom(LONG)
    body = [r for r in agent.screen.text() if r.strip()]
    assert set(body[-2]) == {"─"}, body
    assert len(body[-1]) <= agent.screen.width, (
        f"the line wrapped: {body[-1]!r}")


def test_short_text_still_draws_two_rows(agent):
    agent._set_bottom("⠴ 12s")
    body = [r for r in agent.screen.text() if r.strip()]
    assert body[-1].startswith("⠴")
    assert set(body[-2]) == {"─"}


def test_clearing_leaves_no_rule_and_no_line(agent):
    agent._set_bottom(LONG)
    agent._clear_bottom()
    left = [r for r in agent.screen.text() if r.strip()]
    assert not [r for r in left if set(r) == {"─"}], "\n".join(left)
    assert not any("First turn" in r for r in left), "\n".join(left)


def test_a_repaint_cycle_does_not_grow_the_screen(agent):
    before = len(agent.screen.rows)
    for _ in range(10):
        agent._set_bottom(LONG)
        agent._clear_bottom()
    assert len(agent.screen.rows) <= before + 2, (
        f"the screen grew by {len(agent.screen.rows) - before} rows")


# -- and the way out --------------------------------------------------------

def test_leaving_tears_the_zone_down(agent):
    """The shell's prompt lands on a clean row, not inside a rule."""
    agent._set_bottom("⠴ working")
    agent._teardown_screen()
    left = [r for r in agent.screen.text() if r.strip()]
    assert not [r for r in left if set(r) == {"─"}], "\n".join(left)


def test_tearing_down_twice_is_harmless(agent):
    agent._set_bottom("⠴ working")
    agent._teardown_screen()
    agent._teardown_screen()


def test_every_exit_from_the_loop_goes_through_the_teardown():
    """The faults were not in the drawing but in the leaving: three
    returns and an exception path, and none of them erased anything."""
    import inspect
    src = inspect.getsource(R.TerminalAgent.run)
    finally_block = src.split("finally:")[-1]
    assert "_teardown_screen" in finally_block, (
        "the teardown must sit in the finally, or the path that skips it "
        "is the one a user meets")


# -- the live composer -----------------------------------------------------

def test_input_and_status_stay_below_a_streaming_answer():
    """The real working-prompt shape, driven on one shared terminal.

    A state-only assertion can say the composer exists while its escape
    sequence has actually erased the answer.  This screen applies the
    movements and proves all three things coexist, then proves the next
    streamed delta resumes at the answer rather than inside the prompt.
    """
    import time

    class Engine:
        client = type("Client", (), {"model": "kit.glm-5.3"})()
        kit_permissions = type("Permissions", (), {"mode": "default"})()
        token_usage = {}

        @staticmethod
        def get_status():
            return {}

    screen = WrappingScreen(40, ["earlier output", ""])
    agent = R.TerminalAgent(
        Engine(), out=screen, err=screen, opts=R.ReplOptions(color="never"))
    agent.transcript.width = 40
    agent._turn_t0 = time.monotonic()
    agent._turn_active.set()

    agent.transcript.answer("Antwort")
    agent._repaint_bottom(force=True)
    shown = screen.text()
    answer_at = shown.index("Antwort")
    rules_at = [i for i, row in enumerate(shown) if row and set(row) == {"─"}]
    input_at = next(i for i, row in enumerate(shown) if row.startswith(">"))
    status_at = next(i for i, row in enumerate(shown) if "esc to interrupt" in row)
    hint_at = next(i for i, row in enumerate(shown) if row.startswith("  default"))
    assert (status_at, rules_at, input_at, hint_at) == (
        answer_at + 1, [answer_at + 2, answer_at + 4],
        answer_at + 3, answer_at + 5)

    agent._render_around_bottom(R.RenderItem("text", text=" bleibt stehen"))
    agent._draw_input_line("naechste Nachricht")
    shown = screen.text()
    assert "Antwort bleibt stehen" in shown
    assert any(row == "> naechste Nachricht" for row in shown)
    input_at = next(i for i, row in enumerate(shown)
                    if row == "> naechste Nachricht")
    assert "esc to interrupt" in shown[input_at - 2], (
        "live progress belongs above the stable writing box")
    assert shown[input_at + 2].startswith("  default"), (
        "the same key hint as the idle prompt belongs below the box")
    assert len(screen.rules()) == 2

    agent._clear_bottom()
    assert "Antwort bleibt stehen" in screen.text()
    assert not screen.rules()
