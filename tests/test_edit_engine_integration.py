"""What api_client's edit tools must do once they route through
edit_engine.

These are RED against today's api_client (they xfail strictly), because
today:

* a whitespace-drifted old_string is applied by the silent fuzzy
  fallback (api_client.py:15146) instead of being diagnosed,
* a no-match error names no line number or near-miss text (s7),
* an edit that leaves Python unparseable is applied without a word
  (s6).

After the integration (routing _execute_edit_file / _execute_multi_edit
through edit_engine.apply_edit / apply_multi_edit, keeping the
permission/baseline/write apparatus unchanged) each test must turn
green. strict=True: if one starts passing without the integration, or
keeps failing after it, that is a finding, not noise.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
)


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _read_file(perms, path):
    return _doc_executor.execute("read_file", {"path": str(path)}, perms)


def _edit_file(perms, path, **args):
    args["path"] = str(path)
    return _doc_executor.execute("edit_file", args, perms)


@pytest.mark.xfail(
    reason="api_client still applies drifted matches via the silent "
           "fuzzy fallback; edit_engine diagnoses instead",
    strict=True,
)
def test_a_drifted_match_is_diagnosed_not_applied(workspace):
    """The s7/silent-fuzzy replacement: the engine's near-miss report
    (line number, actual text) reaches the tool result, and the file is
    untouched."""
    t = workspace / "y.py"
    t.write_text("def m():\n    if cond:\n        return 1\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t,
                     old_string="if cond:\n    return 1",
                     new_string="if cond:\n    return 2")
    assert "error" in out
    assert "line 2" in out            # the near miss, with its line
    assert t.read_text() == "def m():\n    if cond:\n        return 1\n"


@pytest.mark.xfail(
    reason="api_client's no-match error carries no line information",
    strict=True,
)
def test_a_no_match_error_names_the_near_miss_line(workspace):
    t = workspace / "big.py"
    t.write_text("def a():\n    return 1\n\n\ndef b():\n    return 2\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # The block stands at a different indent — a near miss.
    out = _edit_file(perms, t,
                     old_string="def b():\nreturn 2",
                     new_string="def b():\nreturn 3")
    assert "error" in out
    assert "line 5" in out
    assert "def b():" in out          # the actual text there, verbatim


@pytest.mark.xfail(
    reason="api_client reports the match count but not the lines",
    strict=True,
)
def test_an_ambiguous_match_names_the_lines(workspace):
    t = workspace / "a.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="a", new_string="A")
    assert "error" in out
    assert "lines 1, 3" in out


@pytest.mark.xfail(
    reason="api_client applies syntax-destroying edits without a word",
    strict=True,
)
def test_a_syntax_regression_is_reported_in_the_result(workspace):
    """The s6 case: the edit lands (that stays the caller's policy),
    but the result says the file no longer parses, with the line."""
    t = workspace / "z.py"
    t.write_text("def a():\n    x = 1\n    return x\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t,
                     old_string="    x = 1\n    return x",
                     new_string="    x = 1\n  return x")
    assert "Edited" in out            # still applied
    assert "syntax" in out.lower()    # but flagged
    assert "line 3" in out
