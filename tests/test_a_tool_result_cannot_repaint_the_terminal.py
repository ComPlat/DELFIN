"""Tool output is arbitrary file content, and a terminal renders what it gets.

The dashboard got this for free: everything went through html.escape, so a
file full of escape sequences was just text. A terminal has no such layer.
`read_file` on a log with ANSI in it, `bash` running something that forces
colour, a repository that ships a file with a title-setting OSC string —
each of those can move the cursor, clear the screen, or overwrite the line
that was just written, which is how output hides itself.

Also in here: the formatting rules that keep a transcript readable when
the model streams a sentence and then calls a tool in the middle of it.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import repl, repl_render as rr


PLAIN = rr.Theme(enabled=False)


# ---------------------------------------------------------------------------
# Control characters
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("payload, expected", [
    ("ok\x1b[2Jdone", "okdone"),                     # clear screen
    ("ok\x1b[1;1Hdone", "okdone"),                   # cursor home
    ("ok\x1b]0;pwned\x07done", "okdone"),            # OSC title, BEL
    ("ok\x1b]0;pwned\x1b\\done", "okdone"),          # OSC title, ST
    ("ok\rdone", "okdone"),                          # overwrite the line
    ("ok\x08done", "okdone"),                        # backspace
    ("ok\x00done", "okdone"),                        # NUL
    ("ok\x1b[31mred\x1b[0m", "okred"),               # colour
])
def test_an_escape_sequence_never_survives(payload, expected):
    assert rr.strip_control(payload) == expected


def test_the_characters_a_transcript_needs_are_kept():
    assert rr.strip_control("a\tb\nc") == "a\tb\nc"


def test_a_tool_result_line_carries_no_escape():
    line = rr.tool_result_line(
        "read_file", "\x1b[31mred\x1b[0m\n\x1b]0;title\x07", theme=PLAIN)
    assert "\x1b" not in line


def test_a_tool_headline_carries_no_escape():
    line = rr.tool_headline(
        "bash", {"command": "echo \x1b[2J"}, theme=PLAIN)
    assert "\x1b" not in line


# ---------------------------------------------------------------------------
# Naming a call by the argument that matters
# ---------------------------------------------------------------------------

def test_the_headline_names_the_interesting_argument():
    assert rr.tool_headline(
        "bash", '{"command": "pytest -x tests/"}', theme=PLAIN
    ).endswith("pytest -x tests/")
    line = rr.tool_headline(
        "read_file", {"file_path": "/a/b.py", "offset": 0}, theme=PLAIN)
    assert "/a/b.py" in line and "offset" not in line


def test_the_engines_json_string_is_accepted_as_given():
    """stream_response hands over tool_input as a JSON STRING, not a dict."""
    assert "ls -la" in rr.tool_headline("bash", '{"command": "ls -la"}',
                                        theme=PLAIN)


def test_a_heredoc_becomes_one_line():
    line = rr.tool_headline("bash", {"command": "a\nb\nc"}, theme=PLAIN)
    assert "\n" not in line
    assert "a b c" in line


def test_a_long_argument_cannot_grow_the_line():
    line = rr.tool_headline("bash", {"command": "x" * 500}, width=80,
                            theme=PLAIN)
    assert len(line) <= 80


def test_both_ends_of_a_long_command_survive():
    """A shell command means one thing at the front and another at the back.

    Cutting the tail hides which file is about to be written, which is the
    single most important word in the line.
    """
    cmd = "rm -rf " + "a" * 200 + "/the_important_target"
    line = rr.tool_headline("bash", {"command": cmd}, width=60, theme=PLAIN)
    assert "rm -rf" in line
    assert "target" in line


def test_unparseable_input_does_not_raise():
    line = rr.tool_headline("weird_tool", "not json at all", theme=PLAIN)
    assert "weird_tool" in line


def test_an_unknown_tool_with_args_does_not_dump_them():
    line = rr.tool_headline("mystery", {"alpha": "x" * 400, "beta": 2},
                            theme=PLAIN)
    assert "mystery" in line
    assert "xxxx" not in line


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def test_a_blocked_result_says_blocked_and_why():
    line = rr.tool_result_line(
        "write_file", "", meta={"ok": False, "error": "path is on the deny-list"},
        theme=PLAIN)
    assert "blocked" in line
    assert "deny-list" in line


def test_a_truncated_result_says_it_was_truncated():
    line = rr.tool_result_line(
        "read_file", "head", meta={"ok": True, "truncated": True,
                                   "chars": 100_000}, theme=PLAIN)
    assert "truncated" in line
    assert "kB" in line or "MB" in line


def test_a_count_from_a_head_slice_is_reported_as_a_floor():
    """Without meta the numbers come from a 2000-char slice.

    Stating them as exact would be an estimate wearing the clothes of a
    measurement, which is the shape of several defects already fixed here.
    """
    line = rr.tool_result_line("read_file", "a\nb\nc", theme=PLAIN)
    assert "+" in line


# ---------------------------------------------------------------------------
# Colour
# ---------------------------------------------------------------------------

def test_colour_is_off_when_the_output_is_not_a_terminal():
    assert rr.theme_for(io.StringIO()).enabled is False


def test_no_color_is_honoured(monkeypatch):
    class _Tty(io.StringIO):
        def isatty(self):
            return True

    monkeypatch.setenv("NO_COLOR", "1")
    assert rr.theme_for(_Tty()).enabled is False
    monkeypatch.delenv("NO_COLOR")
    monkeypatch.setenv("TERM", "dumb")
    assert rr.theme_for(_Tty()).enabled is False


def test_colour_can_be_forced_either_way():
    assert rr.theme_for(io.StringIO(), "always").enabled is True

    class _Tty(io.StringIO):
        def isatty(self):
            return True

    assert rr.theme_for(_Tty(), "never").enabled is False


# ---------------------------------------------------------------------------
# The transcript
# ---------------------------------------------------------------------------

def test_a_tool_line_never_lands_inside_a_sentence():
    out, err = io.StringIO(), io.StringIO()
    t = repl.Transcript(out, err, theme=PLAIN)
    t.answer("Let me check the tests")          # no trailing newline
    t.render(repl.RenderItem("tool_use", name="bash",
                             text='{"command": "pytest"}'))
    t.finish()

    assert out.getvalue() == "Let me check the tests\n"
    assert err.getvalue().startswith("\n"), (
        "the break that closes the answer line belongs on stderr")
    assert "pytest" in err.getvalue()


def test_the_answer_stream_carries_only_the_answer():
    out, err = io.StringIO(), io.StringIO()
    t = repl.Transcript(out, err, theme=PLAIN)
    t.answer("hello ")
    for item in (
        repl.RenderItem("notice", text="retrying"),
        repl.RenderItem("tool_use", name="bash", text='{"command": "ls"}'),
        repl.RenderItem("tool_result", name="bash", text="a\nb"),
        repl.RenderItem("denied", name="write_file"),
        repl.RenderItem("error", text="boom"),
    ):
        t.render(item)
    t.answer("world")
    t.finish()

    assert out.getvalue() == "hello world\n"
    for expected in ("retrying", "ls", "write_file", "boom"):
        assert expected in err.getvalue()


def test_thinking_is_off_unless_it_was_asked_for():
    out, err = io.StringIO(), io.StringIO()
    t = repl.Transcript(out, err, theme=PLAIN, show_thinking=False)
    t.render(repl.RenderItem("thinking", text="hmm"))
    assert err.getvalue() == ""

    t2 = repl.Transcript(out, err, theme=PLAIN, show_thinking=True)
    t2.render(repl.RenderItem("thinking", text="hmm"))
    assert "hmm" in err.getvalue()
