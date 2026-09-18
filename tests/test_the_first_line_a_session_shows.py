"""The first line of a session is the one impression nobody re-reads.

Reported from a live session on kit.glm-5.3, the line under the very
first "Hallo" read:

    ! > ⏳ First turn on kit.glm-5.3: the endpoint is building its cache
      for this prompt — usually about 200 s here ...

Two markers, one meaning. The notice is composed once and read on two
surfaces: the dashboard renders it as markdown, where a leading "> "
sets the line apart, and the terminal has "! " for exactly that. Neither
is wrong on its own; printed together they read like a shell prompt got
loose in the output.

The caveats share that text, so this is fixed where the terminal writes
its own marker, not in the sentence.
"""

from __future__ import annotations

from delfin.agent import repl_render as R
from delfin.agent import verify_guard


def _plain(text: str) -> str:
    return R.strip_ansi(text) if hasattr(R, "strip_ansi") else text


def test_the_cold_start_notice_shows_one_marker():
    said = verify_guard.cold_start_notice("kit.glm-5.3", 200.0)
    assert said, "the notice is produced at all"
    line = _plain(R.notice_line(said)).splitlines()[0]
    assert line.startswith("! "), line
    assert "! >" not in line, (
        "a markdown blockquote marker reached the terminal: " + line)
    assert "First turn" in line


def test_a_quoted_caveat_keeps_its_words():
    """Only the marker goes; the sentence is untouched."""
    out = _plain(R.notice_line("> ⚠️ This answer states 42 results"))
    assert out == "! ⚠️ This answer states 42 results"


def test_every_line_of_a_block_is_unquoted():
    out = _plain(R.notice_line("> first\n> second"))
    assert out.splitlines() == ["! first", "  second"]


def test_a_greater_than_inside_a_line_is_not_a_quote():
    """Harness speech talks about shell redirection and comparisons."""
    out = _plain(R.notice_line("delfin-agent > answer.txt collects it"))
    assert "> answer.txt" in out


def test_only_the_first_marker_is_taken():
    """A line that really begins with two of them keeps the second: the
    fix removes formatting, it does not rewrite the sentence."""
    out = _plain(R.notice_line(">> deeper"))
    assert out == "! > deeper"
