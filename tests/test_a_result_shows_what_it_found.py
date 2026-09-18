"""A count is not a result.

A finished tool call read

    ⏺ bash  python -m pytest -q
      ⎿ 34 lines, 2.1 kB

which says a command ran and nothing about what it found. To learn
whether the suite passed, the reader had to ask the agent what its own
tool had just told it.

The last few lines are where the answer usually is: a pytest tally, the
bottom of a traceback, the path a file was written to. Three of them,
because a thousand-line run must not push the conversation off the
screen — the same discipline the headline already follows, where a
400-line heredoc collapses to one line.
"""

from __future__ import annotations

import pytest

from delfin.agent.repl_render import tool_result_line


SUITE = ("collected 900 items\n\n"
         "tests/test_a.py ......  [ 40%]\n"
         "tests/test_b.py ..F..   [ 80%]\n"
         "2 failed, 17209 passed in 1908.18s")


def _plain(text: str) -> list[str]:
    import re
    return [re.sub(r"\x1b\[[0-9;]*m", "", ln) for ln in text.splitlines()]


def test_the_last_lines_are_shown(the_meta=None):
    lines = _plain(tool_result_line("bash", SUITE,
                                    meta={"chars": len(SUITE), "ok": True},
                                    width=70))
    assert "2 failed, 17209 passed" in lines[-1], lines
    assert len(lines) == 4, "one summary and three lines of it"


def test_the_summary_is_still_there():
    lines = _plain(tool_result_line("bash", SUITE,
                                    meta={"chars": len(SUITE)}, width=70))
    assert "lines" in lines[0] and "B" in lines[0]


def test_a_long_run_does_not_take_the_screen():
    out = "\n".join(f"line {i}" for i in range(1000))
    lines = _plain(tool_result_line("bash", out,
                                    meta={"chars": len(out)}, width=70))
    assert len(lines) == 4, f"{len(lines)} rows for a thousand-line run"
    assert "line 999" in lines[-1]


def test_blank_lines_are_not_spent():
    out = "result\n\n\n\n"
    lines = _plain(tool_result_line("bash", out, meta={"chars": len(out)},
                                    width=70))
    assert len(lines) == 2 and "result" in lines[1]


def test_a_long_line_is_cut_to_the_width():
    out = "x" * 500
    lines = _plain(tool_result_line("bash", out, meta={"chars": len(out)},
                                    width=70))
    assert all(len(ln) <= 70 for ln in lines), [len(ln) for ln in lines]


def test_no_output_says_so_and_shows_nothing():
    assert _plain(tool_result_line("bash", "", meta={})) == ["  ⎿ (no output)"]


def test_a_blocked_call_is_unchanged():
    """A refusal is its own line and gains no excerpt — there is nothing
    to excerpt, and the reason is the result."""
    lines = _plain(tool_result_line(
        "bash", "", meta={"ok": False, "error": "outside the workspace"},
        width=70))
    assert len(lines) == 1
    assert "blocked" in lines[0] and "outside the workspace" in lines[0]
