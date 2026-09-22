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

import json

import pytest

from delfin.agent.repl_render import tool_result_line


SUITE = ("collected 900 items\n\n"
         "tests/test_a.py ......  [ 40%]\n"
         "tests/test_b.py ..F..   [ 80%]\n"
         "2 failed, 17209 passed in 1908.18s")

EDIT = """Edited delfin/agent/repl.py (1 replacement(s)):

--- a/delfin/agent/repl.py
+++ b/delfin/agent/repl.py
@@ -1336,5 +1336,5 @@
         below = rows - 1 - cursor_row
-        if below:
-            out.append("move up")
+        if below or cursor_col != natural_col:
+            out.append("move when needed")
         self.err.write("".join(out))
"""

PATCH = """--- a/delfin/agent/repl.py
+++ b/delfin/agent/repl.py
@@ -10,3 +10,3 @@
 keep
-old
+new
 keep too
"""


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


def test_an_edit_shows_the_file_counts_and_relevant_diff_lines():
    lines = _plain(tool_result_line(
        "edit_file", EDIT, meta={"ok": True, "chars": len(EDIT)},
        width=78))

    assert lines[0] == (
        "  ⎿ Edited delfin/agent/repl.py (+2 -2)"), lines
    assert any("1337 -" in line and "if below" in line for line in lines)
    assert any("1337 +" in line and "natural_col" in line for line in lines)
    assert not any("lines," in line for line in lines), lines
    assert all(len(line) <= 78 for line in lines)


def test_apply_patch_uses_the_matching_calls_diff():
    output = json.dumps({
        "status": "ok", "files_touched": ["delfin/agent/repl.py"],
    })
    lines = _plain(tool_result_line(
        "mcp__kit-coding__apply_patch", output,
        meta={"ok": True, "chars": len(output)},
        tool_input={"diff": PATCH}, width=78))

    assert lines[0] == (
        "  ⎿ Edited delfin/agent/repl.py (+1 -1)"), lines
    assert any("11 -old" in line for line in lines)
    assert any("11 +new" in line for line in lines)


def test_a_failed_patch_never_claims_that_it_edited_the_file():
    output = json.dumps({"status": "check_failed", "error": "does not apply"})
    lines = _plain(tool_result_line(
        "apply_patch", output, meta={"ok": True, "chars": len(output)},
        tool_input={"diff": PATCH}, width=78))
    assert not any("Edited" in line for line in lines)


def test_a_large_edit_preview_is_bounded_but_says_what_was_left_out():
    old = "\n".join(f"-old {i}" for i in range(80))
    new = "\n".join(f"+new {i}" for i in range(80))
    diff = ("Edited many.py:\n\n--- a/many.py\n+++ b/many.py\n"
            "@@ -1,80 +1,80 @@\n" + old + "\n" + new)
    lines = _plain(tool_result_line("edit_file", diff, width=70))

    assert len(lines) == 26, "header + 24 diff rows + omission notice"
    assert "more diff lines" in lines[-1]


def test_an_edit_diff_cannot_smuggle_terminal_control_sequences():
    dirty = PATCH.replace("+new", "+new\x1b[2J\x1b]0;pwned\x07")
    rendered = tool_result_line(
        "apply_patch", '{"status": "ok", "files_touched": ["x.py"]}',
        meta={"ok": True}, tool_input={"diff": dirty})
    assert "\x1b" not in rendered
