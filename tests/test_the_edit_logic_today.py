"""How edit_file / multi_edit behave TODAY, pinned by test.

Characterisation of ``_DocToolExecutor.execute("edit_file"|"multi_edit",
...)`` in delfin/agent/api_client.py, exercised the same way
tests/test_editblock_and_bundle.py does it: a real executor, a real
KitToolPermissions in "default" mode, a file under tmp_path.

What is pinned here is deliberate input for the edit_engine work:

* the special cases the engine must keep answering identically
  (empty old_string, old == new, replace_all, ambiguity, the
  read-before-edit baseline),
* one behavior the engine intentionally REPLACES: the whitespace-tolerant
  fallback that silently applies a drifted match (api_client.py:15146).
  The engine must diagnose a near miss instead of applying it,
* one gap the engine closes: a no-match error that names neither a line
  number nor the text that is nearly there (the s7 report),
* one regression the engine adds: a Python file that parsed before the
  edit and does not after is applied today without a word (the s6
  report),
* and one stale belief it corrects: multi_edit IS atomic on today's
  code (the in-memory loop at api_client.py:15268 writes only after all
  edits validated). The s9 report described a state that no longer
  exists; this file is the evidence.
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


def _multi_edit(perms, path, edits):
    return _doc_executor.execute("multi_edit",
                                 {"path": str(path), "edits": edits}, perms)


PY_BLOCK = "def m():\n    if cond:\n        print('hi')\n        return 1\n"


# ---------------------------------------------------------------------------
# Special cases the engine must keep answering identically
# ---------------------------------------------------------------------------

def test_empty_old_string_is_rejected(workspace):
    t = workspace / "a.py"
    t.write_text("x = 1\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="", new_string="y")
    assert "error" in out
    assert "old_string is required" in out


def test_old_equals_new_is_rejected(workspace):
    t = workspace / "a.py"
    t.write_text("x = 1\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="x = 1", new_string="x = 1")
    assert "error" in out
    assert "must differ" in out


def test_exact_unique_match_is_applied_once(workspace):
    t = workspace / "a.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # "a\nb" occurs once (the second "a" has no "b" after it).
    out = _edit_file(perms, t, old_string="a\nb", new_string="A\nB")
    assert "Edited" in out
    assert t.read_text() == "A\nB\na\n"  # only the first occurrence


def test_replace_all_replaces_every_occurrence(workspace):
    t = workspace / "a.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="a", new_string="A",
                     replace_all=True)
    assert "Edited" in out
    assert t.read_text() == "A\nb\nA\n"


def test_ambiguous_match_without_replace_all_is_rejected(workspace):
    t = workspace / "a.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="a", new_string="A")
    assert "error" in out
    assert "matches 2 times" in out
    # Today's message names the count but no line numbers.
    assert t.read_text() == "a\nb\na\n"


def test_edit_requires_a_read_baseline(workspace):
    t = workspace / "a.py"
    t.write_text("x = 1\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = _edit_file(perms, t, old_string="x = 1", new_string="x = 2")
    assert "error" in out
    assert "read_file" in out


# ---------------------------------------------------------------------------
# What the engine replaces / closes. Each of these documents TODAY'S gap
# with the observed behavior; edit_engine must do better (see
# test_edit_engine_regressions.py).
# ---------------------------------------------------------------------------

def test_today_a_drifted_match_is_applied_silently(workspace):
    """s7's cousin: the whitespace-tolerant fallback (api_client.py:15146)
    APPLIES a near miss without asking. The engine must diagnose, not
    apply."""
    t = workspace / "y.py"
    t.write_text(PY_BLOCK)
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # old_string at column-0, file has it indented — a drift, not a match.
    out = _edit_file(perms, t,
                     old_string="if cond:\n    print('hi')\n    return 1",
                     new_string="if cond:\n    print('bye')\n    return 7")
    assert "Edited" in out
    assert "fuzzy" in out.lower()  # it went ahead and applied it


def test_today_a_no_match_error_names_no_line(workspace):
    """s7: 'old_string not found' on a big file, with no hint WHERE the
    text nearly stands. The engine must report line number + actual text."""
    t = workspace / "big.py"
    t.write_text("def a():\n    return 1\n\n\ndef b():\n    return 2\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # Nothing matches, not even the whitespace-tolerant fallback (the
    # text simply is not there): today's error names no line number.
    out = _edit_file(perms, t, old_string="def c():\n    return 3",
                     new_string="def c():\n    return 4")
    assert "error" in out
    msg = out
    assert "not found" in msg
    # No line number anywhere in today's message:
    assert "line" not in msg.lower()


def test_today_a_syntax_destroying_edit_is_applied_without_a_word(workspace):
    """s6: the edit matches exactly but destroys the following block's
    indentation. Today it is written without any syntax check."""
    t = workspace / "z.py"
    t.write_text(
        "def a():\n"
        "    x = 1\n"
        "    return x\n"
        "\n"
        "def b():\n"
        "    y = 2\n"
        "    return y\n"
    )
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # Exact unique match; new_string uses 2-space indent where the file
    # uses 4 — the result mixes indents inside one block and no longer
    # parses. Today it is written without any syntax check.
    out = _edit_file(perms, t,
                     old_string="    return x",
                     new_string="  return x")
    assert "Edited" in out  # applied...
    import ast
    with pytest.raises(SyntaxError):
        ast.parse(t.read_text())  # ...and left the file unparseable


def test_today_multi_edit_is_atomic(workspace):
    """s9 said a failed batch left earlier edits applied. On today's code
    (in-memory loop, single write at the end) that is NOT the case —
    this test is the evidence, and the engine must keep the behavior."""
    t = workspace / "m.py"
    t.write_text("a\nb\nc\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _multi_edit(perms, t, [
        {"old_string": "a", "new_string": "A"},          # would apply
        {"old_string": "zzz", "new_string": "Z"},        # never matches
        {"old_string": "c", "new_string": "C"},          # would apply
    ])
    assert "error" in out
    assert "edit #2" in out
    assert t.read_text() == "a\nb\nc\n"  # nothing applied


def test_today_multi_edit_names_the_failing_index(workspace):
    t = workspace / "m.py"
    t.write_text("a\nb\na\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _multi_edit(perms, t, [
        {"old_string": "a", "new_string": "A"},          # ambiguous here
        {"old_string": "b", "new_string": "B"},
    ])
    assert "error" in out
    assert "edit #1" in out
    assert "matches 2 times" in out
    assert t.read_text() == "a\nb\na\n"


def test_today_line_endings_survive_an_edit(workspace):
    """text_files.write_text_file(like=shape) already preserves CRLF; the
    engine must not regress this (edit point 5)."""
    t = workspace / "w.csv"
    t.write_bytes(b"alpha\r\nbeta\r\ngamma\r\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    _edit_file(perms, t, old_string="beta", new_string="BETA")
    assert t.read_bytes() == b"alpha\r\nBETA\r\ngamma\r\n"
