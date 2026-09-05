"""Editing while the agent works, and a terminal that changed size.

Three gaps, all of the same shape: something the idle prompt gets from
readline for free, and the in-turn reader never had.

The decoder had a buffer and no cursor, so every edit happened at the
end — a typo four words back cost four words. The width was read once at
startup, so a resized window truncated to the size it used to be for the
rest of the session. And a line typed during a turn never entered the
history, which is most of what gets typed in a busy session.
"""

from __future__ import annotations

import io
import signal

from delfin.agent import repl_keys as rk
from delfin.agent.repl import ReplOptions, TerminalAgent, Transcript


def _typed(text: str) -> rk.KeyDecoder:
    d = rk.KeyDecoder()
    d.feed(text)
    return d


LEFT, RIGHT = "\x1b[D", "\x1b[C"
HOME, END = "\x1b[H", "\x1b[F"
ALT_B, ALT_F = "\x1bb", "\x1bf"
CTRL_A, CTRL_E, CTRL_W, CTRL_K, CTRL_U = "\x01", "\x05", "\x17", "\x0b", "\x15"


# ---------------------------------------------------------------------------
# The cursor
# ---------------------------------------------------------------------------

def test_a_fresh_decoder_has_its_cursor_at_the_start():
    assert rk.KeyDecoder().cursor == 0


def test_typing_moves_the_cursor_with_it():
    d = _typed("hello")
    assert (d.buffer, d.cursor) == ("hello", 5)


def test_a_character_typed_after_moving_back_lands_there():
    """The whole point. Appending regardless of the cursor is what the
    buffer-without-a-cursor did."""
    d = _typed("helo")
    d.feed(LEFT)
    d.feed("l")
    assert d.buffer == "hello"


def test_backspace_deletes_at_the_cursor_and_not_at_the_end():
    d = _typed("helllo")
    d.feed(LEFT + LEFT)          # between the third l and the o
    d.feed("\x7f")
    assert d.buffer == "hello"


def test_the_cursor_cannot_leave_the_line():
    d = _typed("ab")
    for _ in range(5):
        d.feed(LEFT)
    assert d.cursor == 0
    d.feed("\x7f")               # nothing to delete before the start
    assert d.buffer == "ab"
    for _ in range(5):
        d.feed(RIGHT)
    assert d.cursor == 2


def test_home_and_end_and_their_control_spellings():
    for key in (HOME, CTRL_A):
        d = _typed("some text")
        d.feed(key)
        assert d.cursor == 0, key
    for key in (END, CTRL_E):
        d = _typed("some text")
        d.feed(HOME)
        d.feed(key)
        assert d.cursor == len("some text"), key


def test_both_spellings_of_an_arrow_are_understood():
    """Applications-cursor mode sends O where normal mode sends [, and
    which mode a terminal is in is not this layer's decision."""
    for key in ("\x1b[D", "\x1bOD"):
        d = _typed("abc")
        d.feed(key)
        assert d.cursor == 2, key


# ---------------------------------------------------------------------------
# Words
# ---------------------------------------------------------------------------

def test_ctrl_w_deletes_the_word_behind_the_cursor():
    d = _typed("read the wrong file")
    d.feed(CTRL_W)
    assert d.buffer == "read the wrong "


def test_ctrl_w_after_a_space_still_deletes_a_word():
    """Stopping at the first space would spend a keystroke deleting
    nothing, and after a space is where the cursor usually is."""
    d = _typed("one two ")
    d.feed(CTRL_W)
    assert d.buffer == "one "


def test_ctrl_w_in_the_middle_keeps_the_tail():
    d = _typed("alpha beta gamma")
    d.feed(HOME)
    for _ in range(len("alpha beta")):
        d.feed(RIGHT)
    d.feed(CTRL_W)
    assert d.buffer == "alpha  gamma"


def test_ctrl_k_cuts_to_the_end():
    d = _typed("keep this drop that")
    d.feed(HOME)
    for _ in range(len("keep this ")):
        d.feed(RIGHT)
    d.feed(CTRL_K)
    assert d.buffer == "keep this "


def test_alt_b_and_alt_f_move_by_word():
    d = _typed("alpha beta gamma")
    d.feed(ALT_B)
    assert d.buffer[d.cursor:] == "gamma"
    d.feed(ALT_B)
    assert d.buffer[d.cursor:] == "beta gamma"
    d.feed(ALT_F)
    assert d.buffer[d.cursor:] == "gamma"


def test_alt_f_at_the_end_stays_at_the_end():
    d = _typed("one two")
    d.feed(ALT_F)
    assert d.cursor == len("one two")


def test_moving_and_deleting_agree_about_where_a_word_starts():
    """One rule, two keys. If they disagreed, Alt+B then Ctrl+W would
    delete a different word than the one just stepped over."""
    for text in ("alpha beta", "alpha beta ", "a  b", "single"):
        moved = _typed(text)
        moved.feed(ALT_B)
        deleted = _typed(text)
        deleted.feed(CTRL_W)
        # Where the move landed IS where the delete cut, and the text
        # that survives is everything up to it.
        assert moved.cursor == deleted.cursor, text
        assert deleted.buffer == text[:moved.cursor], text


# ---------------------------------------------------------------------------
# The cursor survives what resets the line
# ---------------------------------------------------------------------------

def test_enter_resets_the_cursor_with_the_buffer():
    d = _typed("send this")
    d.feed(HOME)
    d.feed("\r")
    assert (d.buffer, d.cursor) == ("", 0)


def test_ctrl_g_resets_the_cursor_with_the_buffer():
    d = _typed("steer this")
    d.feed(HOME)
    d.feed("\x07")
    assert (d.buffer, d.cursor) == ("", 0)


def test_ctrl_u_resets_the_cursor_with_the_buffer():
    d = _typed("scrap this")
    d.feed(CTRL_U)
    assert (d.buffer, d.cursor) == ("", 0)


def test_a_paste_lands_where_the_cursor_is():
    """Appending regardless would put pasted text where the user is not
    looking, which is the same defect as typing did."""
    d = _typed("before  after")
    d.feed(HOME)
    for _ in range(len("before ")):
        d.feed(RIGHT)
    d.feed("\x1b[200~PASTED\x1b[201~")
    assert d.buffer == "before PASTED after"


# ---------------------------------------------------------------------------
# The terminal changed size
# ---------------------------------------------------------------------------

class _Tty(io.StringIO):
    def isatty(self):
        return True


def test_an_explicit_width_is_never_re_asked(monkeypatch):
    """The tests need a stable width; the session needs a live one."""
    t = Transcript(_Tty(), _Tty(), width=72)
    monkeypatch.setattr("delfin.agent.repl_render.terminal_width",
                        lambda default=100: 200)
    assert t.width == 72
    assert t.refresh_width() == 72


def test_a_live_width_follows_the_terminal(monkeypatch):
    sizes = iter([80, 120])
    monkeypatch.setattr("delfin.agent.repl_render.terminal_width",
                        lambda default=100: next(sizes))
    t = Transcript(_Tty(), _Tty())
    assert t.width == 80
    assert t.refresh_width() == 120
    assert t.width == 120, "and it sticks until the next resize"


def test_the_handler_records_and_does_not_paint(monkeypatch):
    """A signal handler runs between bytecodes of whatever the main
    thread was doing — including halfway through a write to the same
    terminal. Painting there interleaves with the write it interrupted."""
    err = _Tty()
    agent = TerminalAgent(_engine(), out=_Tty(), err=err,
                          opts=ReplOptions(color="never"))
    before = err.getvalue()
    agent._on_sigwinch(getattr(signal, "SIGWINCH", 28), None)
    assert agent._width_dirty is True
    assert err.getvalue() == before, "the handler wrote to the terminal"


def test_the_loop_is_what_repaints():
    import inspect
    from delfin.agent import repl
    src = inspect.getsource(repl.TerminalAgent._pump)
    assert "_width_dirty" in src


def test_a_platform_without_sigwinch_still_starts(monkeypatch):
    monkeypatch.delattr(signal, "SIGWINCH", raising=False)
    agent = TerminalAgent(_engine(), out=_Tty(), err=_Tty(),
                          opts=ReplOptions(color="never"))
    agent._install_sigwinch()        # must not raise
    agent._restore_sigwinch()


# ---------------------------------------------------------------------------
# History
# ---------------------------------------------------------------------------

def _engine():
    return type("E", (), {
        "session_id": "sid", "messages": [], "kit_permissions": None,
        "mode": "solo", "provider": "kit",
        "client": type("C", (), {"model": "m"})(),
        "steer": staticmethod(lambda t: True),
        "get_status": staticmethod(lambda: {}),
    })()


def test_a_line_typed_during_a_turn_is_recallable(monkeypatch):
    """Only the idle prompt goes through readline, so everything typed
    while the agent worked was missing from the up-arrow afterwards."""
    added: list[str] = []
    import readline
    monkeypatch.setattr(readline, "add_history", added.append)

    agent = TerminalAgent(_engine(), out=_Tty(), err=_Tty(),
                          opts=ReplOptions(color="never"))
    agent._recall_later("run the failing test")
    agent._steer("and read the log")
    assert added == ["run the failing test", "and read the log"]


def test_recalling_is_not_remembering():
    """Two things a user might call remembering.

    `_remember` writes a durable MEMORY for `#note`. Conflating the names
    would have turned every queued line into a stored memory.
    """
    import inspect
    from delfin.agent import repl
    doc = inspect.getdoc(repl.TerminalAgent._remember) or ""
    assert "memory" in doc.lower()
    assert repl.TerminalAgent._remember is not repl.TerminalAgent._recall_later


def test_an_empty_line_is_not_recalled(monkeypatch):
    added: list[str] = []
    import readline
    monkeypatch.setattr(readline, "add_history", added.append)
    agent = TerminalAgent(_engine(), out=_Tty(), err=_Tty(),
                          opts=ReplOptions(color="never"))
    agent._recall_later("")
    assert added == []


def test_only_this_sessions_lines_are_appended():
    """Writing the whole buffer means two terminals each load the file at
    start and write their own copy at exit, so whichever quits last
    erases the other's session."""
    import inspect
    from delfin.agent import repl
    src = inspect.getsource(repl.TerminalAgent._save_history)
    assert "append_history_file" in src
    assert "_hist_base" in src
    load = inspect.getsource(repl.TerminalAgent._load_history)
    assert "_hist_base" in load, "the baseline has to be taken after loading"
