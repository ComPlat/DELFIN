"""A pasted block is one message, and the answer is styled on stdout only.

Two things that were named as missing and are now built.

**Paste.** During a turn every newline is a submit, because that is what a
newline means when someone is typing. Pasted text contains newlines that
mean nothing of the kind — a stack trace pasted into a running turn queued
one message per line. The only way a terminal program can tell the two
apart is bracketed paste: the terminal is asked to wrap pasted text in
markers, and everything between them is content.

**Markdown on stdout.** The chrome theme is a decision about stderr. The
answer goes to stdout, and the two are redirected independently — so the
answer's colour has to be decided about its own stream, or
``delfin-agent -p '…' > answer.txt`` run from a terminal would write
escape codes into the file.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import repl_keys as rk
from delfin.agent.repl import Transcript


PASTE_START = "\x1b[200~"
PASTE_END = "\x1b[201~"


def _kinds(events):
    return [e.kind for e in events]


# ---------------------------------------------------------------------------
# The whole point
# ---------------------------------------------------------------------------

def test_a_pasted_block_is_one_buffer_not_three_submits():
    d = rk.KeyDecoder()
    pasted = "line one\nline two\nline three"
    events = d.feed(PASTE_START + pasted + PASTE_END)

    assert rk.SUBMIT not in _kinds(events), (
        "a newline inside a paste is content, not a key")
    assert d.buffer == pasted
    assert "\n" in d.buffer


def test_the_same_text_typed_really_does_submit():
    """The control. Without it this could pass by disabling newlines."""
    d = rk.KeyDecoder()
    events = d.feed("line one\nline two\n")
    assert _kinds(events).count(rk.SUBMIT) == 2


def test_enter_after_a_paste_submits_the_whole_thing():
    d = rk.KeyDecoder()
    d.feed(PASTE_START + "a\nb" + PASTE_END)
    events = d.feed("\r")
    submits = [e for e in events if e.kind == rk.SUBMIT]
    assert len(submits) == 1
    assert submits[0].text == "a\nb"


# ---------------------------------------------------------------------------
# Read boundaries, which is where this actually gets decided
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cut", list(range(2, 40)))
def test_a_split_anywhere_yields_the_same_buffer(cut):
    src = PASTE_START + "alpha\nbeta\ngamma" + PASTE_END
    if cut >= len(src):
        pytest.skip("beyond the input")
    whole = rk.KeyDecoder()
    whole.feed(src)

    split = rk.KeyDecoder()
    split.feed(src[:cut])
    split.feed(src[cut:])
    assert split.buffer == whole.buffer, f"cut at {cut}"


def test_a_read_that_returns_only_the_markers_escape_reads_as_an_interrupt():
    """The one split this does not survive, named rather than hidden.

    A chunk consisting of nothing but ESC is the signature of the user
    pressing Esc: the read returned one byte because the terminal had one
    byte. Claiming it for the paste marker instead would mean a key that
    ends a running turn stops working in order to survive a read boundary
    landing at a one-byte offset of a burst the terminal writes at once.

    Written down as the accepted cost, and as the reason the prefix test
    starts at two characters.
    """
    d = rk.KeyDecoder()
    first = d.feed(PASTE_START[:1])
    assert _kinds(first) == [rk.INTERRUPT]


def test_a_start_marker_cut_in_half_is_still_a_start_marker():
    d = rk.KeyDecoder()
    d.feed(PASTE_START[:3])
    d.feed(PASTE_START[3:] + "x\ny" + PASTE_END)
    assert d.buffer == "x\ny"


def test_a_half_written_end_marker_is_not_pasted_as_text():
    d = rk.KeyDecoder()
    d.feed(PASTE_START + "done" + PASTE_END[:4])
    assert d.buffer == "done", d.buffer
    d.feed(PASTE_END[4:])
    assert d.buffer == "done"


def test_typing_before_and_after_a_paste_in_one_read():
    d = rk.KeyDecoder()
    d.feed("see: " + PASTE_START + "X\nY" + PASTE_END + "!")
    assert d.buffer == "see: X\nY!"


# ---------------------------------------------------------------------------
# Keys are not keys inside a paste
# ---------------------------------------------------------------------------

def test_an_escape_inside_a_paste_does_not_end_the_turn():
    """Pasting a log full of escape codes must not interrupt the agent."""
    d = rk.KeyDecoder()
    events = d.feed(PASTE_START + "before\x1bafter" + PASTE_END)
    assert rk.INTERRUPT not in _kinds(events)
    assert "\x1b" in d.buffer


def test_a_lone_escape_outside_a_paste_still_ends_the_turn():
    """The control for the one above."""
    d = rk.KeyDecoder()
    assert _kinds(d.feed("\x1b")) == [rk.INTERRUPT]


def test_ctrl_o_inside_a_paste_is_a_character():
    d = rk.KeyDecoder()
    events = d.feed(PASTE_START + "a\x0fb" + PASTE_END)
    assert rk.EXPAND not in _kinds(events)


# ---------------------------------------------------------------------------
# Asking the terminal for it, and stopping
# ---------------------------------------------------------------------------

def test_the_terminal_is_asked_and_then_told_to_stop(monkeypatch):
    """Left enabled, a later shell in the same terminal receives the
    markers as literal text — so the disable is as load-bearing as the
    enable, and belongs on the restore path that already runs come what
    may."""
    written: list[str] = []

    class _Sink(io.StringIO):
        def write(self, s):
            written.append(s)
            return len(s)

    def _fake_open(path, mode="r", *a, **kw):
        assert path == "/dev/tty"
        sink = _Sink()
        sink.close = lambda: None       # type: ignore[method-assign]
        return sink

    monkeypatch.setattr("builtins.open", _fake_open)
    raw = rk.RawMode.__new__(rk.RawMode)
    raw._set_paste_mode(True)
    raw._set_paste_mode(False)
    assert written == ["\x1b[?2004h", "\x1b[?2004l"]


def test_no_controlling_terminal_is_not_an_error(monkeypatch):
    def _boom(*a, **kw):
        raise OSError("no tty")
    monkeypatch.setattr("builtins.open", _boom)
    raw = rk.RawMode.__new__(rk.RawMode)
    raw._set_paste_mode(True)           # must not raise


def test_the_idle_prompt_pins_it_too():
    """Whether readline brackets a paste is a library-version accident.

    Recent GNU readline defaults it on; libedit does not. Inheriting that
    means the same paste behaves differently on two machines for a reason
    nobody can see.
    """
    import inspect
    from delfin.agent import repl
    src = inspect.getsource(repl.TerminalAgent._install_completer)
    assert "enable-bracketed-paste" in src


# ---------------------------------------------------------------------------
# Markdown reaches the answer stream, and only when that stream wants it
# ---------------------------------------------------------------------------

class _Tty(io.StringIO):
    def isatty(self):
        return True


def test_a_redirected_stdout_gets_the_model_bytes_exactly():
    out, err = io.StringIO(), _Tty()
    t = Transcript(out, err, color="auto")
    for delta in ("# Head", "ing\n", "- one\n", "**b**\n"):
        t.answer(delta)
    t.finish()
    assert out.getvalue() == "# Heading\n- one\n**b**\n"


def test_a_terminal_stdout_gets_the_styling():
    out, err = _Tty(), _Tty()
    t = Transcript(out, err, color="auto")
    t.answer("# Heading\n")
    t.finish()
    assert "\x1b[1m" in out.getvalue()
    assert "#" not in out.getvalue()


def test_colour_never_overrides_the_stream_it_is_written_to():
    """`--color never` is a decision about output, not about one stream."""
    out, err = _Tty(), _Tty()
    t = Transcript(out, err, color="never")
    t.answer("**bold**\n")
    t.finish()
    assert out.getvalue() == "**bold**\n"


def test_a_terminal_stderr_does_not_style_a_piped_stdout():
    """The exact shape of the defect this split exists to prevent."""
    out, err = io.StringIO(), _Tty()
    t = Transcript(out, err, color="auto")
    t.answer("`code`\n")
    t.finish()
    assert "\x1b[" not in out.getvalue()


def test_the_line_open_flag_follows_the_model_not_the_styling():
    """An SGR reset is not a newline.

    If the flag were taken from the styled text, a delta ending in a
    reset would read as a closed line and the status row would repaint
    over the sentence still being written.
    """
    out, err = _Tty(), _Tty()
    t = Transcript(out, err, color="auto")
    t.answer("**bold**")
    assert out.getvalue().endswith("\x1b[0m")
    assert t.answer_open is True


def test_finish_releases_what_the_renderer_held_back():
    out, err = _Tty(), _Tty()
    t = Transcript(out, err, color="auto")
    t.answer("value *")
    t.finish()
    assert out.getvalue().count("*") == 1


def test_finish_is_idempotent():
    out, err = _Tty(), _Tty()
    t = Transcript(out, err, color="auto")
    t.answer("x")
    t.finish()
    before = out.getvalue()
    t.finish()
    assert out.getvalue() == before


# ---------------------------------------------------------------------------
# And the message the model receives still has its lines
# ---------------------------------------------------------------------------

def test_a_pasted_block_keeps_its_lines_all_the_way_to_the_model():
    """The half of bracketed paste that the decoder does not decide.

    The paste arrives as ONE message — the decoder's job, done. But the
    line then goes through `expand_at_references`, which split the whole
    text on whitespace and re-joined it on single spaces: a pasted stack
    trace reached the model as a paragraph, and pasted code reached it
    without the indentation that decides what it means.
    """
    from pathlib import Path

    from delfin.agent import repl_commands as rc

    src = ('Traceback (most recent call last):\n'
           '  File "calc.py", line 3, in add\n'
           '    return a + b\n'
           'TypeError: unsupported operand')
    out = rc.expand_at_references(src, Path("/tmp"))
    assert out == src, "every line, and every indent inside a line"


def test_the_fence_keeps_its_lines_too():
    """`read_block` joins a `\"\"\"` block with newlines, and the very next
    hop discarded them — a multi-line form that could not send one."""
    from pathlib import Path

    from delfin.agent import repl_commands as rc

    assert rc.expand_at_references("one\ntwo\nthree", Path("/tmp")) == \
        "one\ntwo\nthree"


def test_an_at_reference_is_still_annotated_inside_a_block(tmp_path):
    (tmp_path / "calc.py").write_text("x = 1\n")
    from delfin.agent import repl_commands as rc

    out = rc.expand_at_references("look at\n  @calc.py\nplease", tmp_path)
    assert "\n  calc.py\n" in out, "the @ goes, the layout stays"
    assert "Files referenced: calc.py" in out


def test_a_missing_reference_is_still_flagged(tmp_path):
    from delfin.agent import repl_commands as rc

    out = rc.expand_at_references("check\n@nope.py", tmp_path)
    assert "not found" in out


def test_trailing_blank_lines_go_but_interior_ones_stay():
    from pathlib import Path

    from delfin.agent import repl_commands as rc

    out = rc.expand_at_references("first\n\nsecond\n\n\n", Path("/tmp"))
    assert out == "first\n\nsecond"


def test_the_editor_draft_and_a_typed_line_are_treated_alike():
    """One entry point behaving differently from another is a surprise
    with no upside; the editor route bypassed the expansion only because
    the expansion used to destroy what the editor exists to produce."""
    import inspect

    from delfin.agent import repl_commands as rc

    src = inspect.getsource(rc._editor_round_trip)
    assert "expand_at_references" in src
