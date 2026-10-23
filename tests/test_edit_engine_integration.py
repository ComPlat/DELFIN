"""What api_client's edit tools do since they route through edit_engine.

The whitespace-tolerant fallback (editblock.fuzzy_replace) still runs
first and applies a UNIQUE drifted match. What edit_engine adds is what
happens when there is no unique match: every near miss with its line
and verbatim text, the lines of every exact match, and a warning when
an applied edit leaves a Python file unparseable.
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


def test_two_drifted_near_misses_are_named_not_guessed(workspace):
    """The fallback gives up when the drift fits two places; the error
    names both, with line and text, and the file is untouched."""
    t = workspace / "y.py"
    body = ("def m():\n    if cond:\n        return 1\n"
            "def n():\n    if cond:\n        return 1\n")
    t.write_text(body)
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t,
                     old_string="if cond:\n    return 1",
                     new_string="if cond:\n    return 2")
    assert "error" in out
    assert "line 2" in out and "line 5" in out
    assert t.read_text() == body


def test_a_multi_edit_ambiguity_names_the_lines(workspace):
    t = workspace / "big.py"
    body = "def a():\n    return 1\n\n\ndef b():\n    return 1\n"
    t.write_text(body)
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _doc_executor.execute("multi_edit", {"path": str(t), "edits": [
        {"old_string": "def a():", "new_string": "def a2():"},
        {"old_string": "return 1\n", "new_string": "return 9\n",
         "replace_all": False},
    ]}, perms)
    assert "edit #2" in out
    assert "lines 2, 6" in out
    assert t.read_text() == body


def test_an_ambiguous_match_names_the_lines(workspace):
    t = workspace / "a.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="a", new_string="A")
    assert "error" in out
    assert "lines 1, 3" in out


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
