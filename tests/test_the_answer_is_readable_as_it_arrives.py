"""Markdown in a stream that cannot be taken back.

The model writes markdown. On a terminal it arrived as literal asterisks,
backticks and hash marks, so a heading looked like punctuation and a code
block looked like prose with fences around it.

The hard part is not the syntax, it is the medium: deltas arrive in
whatever chunks the provider sends, and by the time a closing marker is
seen the opening text is already on screen. Every test here is written
against that — the same input is fed in one piece, in single characters,
and split at the worst place, and all three must agree.

The other rule the tests enforce is the stream contract: with colour off
the bytes must come out exactly as the model produced them, because
`delfin-agent -p '…' > answer.txt` is a deliverable, not a display.
"""

from __future__ import annotations

import pytest

from delfin.agent.repl_render import MarkdownStream, Theme


ON = Theme(enabled=True)
OFF = Theme(enabled=False)


def _whole(text: str, theme: Theme = ON) -> str:
    md = MarkdownStream(theme)
    return md.feed(text) + md.flush()


def _by_char(text: str, theme: Theme = ON) -> str:
    md = MarkdownStream(theme)
    return "".join(md.feed(c) for c in text) + md.flush()


def _split_at(text: str, i: int, theme: Theme = ON) -> str:
    md = MarkdownStream(theme)
    return md.feed(text[:i]) + md.feed(text[i:]) + md.flush()


def _plain(styled: str) -> str:
    """The visible characters, with every SGR sequence removed."""
    import re
    return re.sub(r"\x1b\[[0-9;]*m", "", styled)


# ---------------------------------------------------------------------------
# The contract that outranks everything else
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "# Heading\n\n- one\n- two\n\n```python\nx = 1\n```\n**bold** and `code`\n",
    "no markdown at all",
    "",
    "*",
    "```",
])
def test_with_colour_off_nothing_is_touched(text):
    assert _whole(text, OFF) == text
    assert _by_char(text, OFF) == text


def test_the_visible_characters_are_never_invented(  # noqa: D103
):
    src = "# Title\n- one\n**b** `c`\n```\nraw\n```\n"
    seen = _plain(_whole(src))
    # Markers are consumed, never replaced by content of our own — except
    # the bullet glyph, which is a deliberate substitution.
    for fragment in ("Title", "one", "b", "c", "raw"):
        assert fragment in seen
    assert "**" not in seen
    assert "```" not in seen
    assert "#" not in seen


# ---------------------------------------------------------------------------
# Chunk boundaries — the reason this class exists at all
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("src", [
    "**bold** tail\n",
    "# Heading\n",
    "- item\n",
    "```py\ncode\n```\ndone\n",
    "`inline` after\n",
    "a * b * c\n",
    "text with ** unmatched\nnext line\n",
])
def test_the_split_does_not_change_the_result(src):
    reference = _whole(src)
    assert _by_char(src) == reference, "character by character"
    for i in range(1, len(src)):
        assert _split_at(src, i) == reference, f"split at {i}"


def test_a_fence_arriving_one_backtick_at_a_time_is_still_a_fence():
    md = MarkdownStream(ON)
    out = md.feed("`") + md.feed("`") + md.feed("`") + md.feed("\nx\n```\n")
    out += md.flush()
    assert "\x1b[2m" in out, "the block is dimmed"
    assert "`" not in _plain(out)


def test_a_single_backtick_at_line_start_is_still_inline_code():
    out = _whole("`x` y\n")
    assert _plain(out) == "x y\n"
    assert "\x1b[36m" in out


# ---------------------------------------------------------------------------
# What is styled
# ---------------------------------------------------------------------------

def test_a_heading_is_bold():
    out = _whole("## Results\n")
    assert out.startswith("\x1b[1m")
    assert _plain(out) == "Results\n"


def test_a_hash_that_is_not_a_heading_is_left_alone():
    assert _plain(_whole("#nothashtag\n")) == "#nothashtag\n"
    assert _plain(_whole("####### seven\n")) == "####### seven\n"


def test_a_bullet_becomes_a_glyph():
    assert _plain(_whole("- one\n")) == "• one\n"
    assert _plain(_whole("* two\n")) == "• two\n"


def test_a_minus_sign_is_not_a_bullet():
    assert _plain(_whole("-5 degrees\n")) == "-5 degrees\n"


def test_bold_opens_and_closes():
    out = _whole("say **this** now\n")
    assert _plain(out) == "say this now\n"
    assert out.count("\x1b[1m") == 1


def test_inline_code_is_coloured():
    out = _whole("run `pytest -q` first\n")
    assert _plain(out) == "run pytest -q first\n"
    assert out.count("\x1b[36m") == 1


def test_markers_inside_a_fence_are_literal():
    """A code block is code. Styling its contents would be a lie about it."""
    out = _whole("```\n**not bold** and `not code`\n```\n")
    assert "**not bold**" in _plain(out)
    assert "`not code`" in _plain(out)


def test_a_heading_inside_a_fence_stays_a_hash():
    out = _whole("```sh\n# a shell comment\n```\n")
    assert "# a shell comment" in _plain(out)


# ---------------------------------------------------------------------------
# The documented failure, bounded
# ---------------------------------------------------------------------------

def test_an_unmatched_marker_costs_one_line_and_no_more():
    """Stated in the class docstring, and enforced here.

    A stream cannot look ahead, so an opening `**` with no closer styles
    what follows. The newline is where that stops.
    """
    out = _whole("**never closed\nplain again\n")
    body = out.split("\n")
    assert "\x1b[0m" in body[0], "the style is dropped at the newline"
    assert "\x1b[1m" not in body[1]


def test_flush_closes_an_open_style():
    md = MarkdownStream(ON)
    first = md.feed("**still open")
    assert "\x1b[1m" in first
    assert md.flush().endswith("\x1b[0m")


def test_flush_emits_a_held_back_tail():
    """A two-character tail that looked like a marker must still print."""
    md = MarkdownStream(ON)
    md.feed("value ")
    held = md.feed("*")
    assert _plain(held) == "", "held, because it may be the start of **"
    assert _plain(md.flush()) == "*", "and released when nothing follows"


def test_an_unclosed_fence_does_not_leave_the_terminal_dim():
    md = MarkdownStream(ON)
    md.feed("```python\nx = 1\n")
    assert md.flush().endswith("\x1b[0m")


# ---------------------------------------------------------------------------
# Volume
# ---------------------------------------------------------------------------

def test_a_long_answer_round_trips(monkeypatch):
    src = "".join(
        f"## Section {i}\n- point one\n- point two\n\n```py\nf({i})\n```\n\n"
        f"Plain text with **emphasis** and `code`.\n\n"
        for i in range(50)
    )
    assert _plain(_whole(src)) == _plain(_by_char(src))
    assert _whole(src, OFF) == src
