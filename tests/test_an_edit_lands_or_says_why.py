"""How edit_file / multi_edit behave, pinned by test.

Written as a characterisation of the pre-edit_engine code (2026-09-25,
runde3-s4) and kept as the contract after edit_engine was wired in:

* the special cases answered identically before and after (empty
  old_string, old == new, replace_all, ambiguity, the read-before-edit
  baseline, CRLF survival),
* the whitespace-tolerant fallback STAYS: a unique drifted match is
  applied and the result says "fuzzy match" (it was never silent -- the
  first draft of edit_engine wanted to drop it; kept because weaker
  models lean on it and it reports itself),
* a Python file that parsed before an edit and does not after is still
  written, but the result now carries a WARNING with the line (s6),
* multi_edit is atomic (in-memory loop, one write at the end); the s9
  report described an older state.
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
# Where edit_engine changed or kept the behavior.
# ---------------------------------------------------------------------------

def test_a_unique_drifted_match_is_applied_and_says_fuzzy(workspace):
    """The whitespace-tolerant fallback applies a UNIQUE near miss and
    says so in the result. Deliberate, kept after edit_engine."""
    t = workspace / "y.py"
    t.write_text(PY_BLOCK)
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    # old_string at column-0, file has it indented — a drift, not a match.
    out = _edit_file(perms, t,
                     old_string="if cond:\n    print('hi')\n    return 1",
                     new_string="if cond:\n    print('bye')\n    return 7")
    assert "Edited" in out
    assert "fuzzy" in out.lower()  # applied, and it says how


def test_a_plain_miss_says_nothing_was_replaced(workspace):
    """Nothing near the old_string: no near miss to name, and the error
    says the file is untouched."""
    t = workspace / "big.py"
    t.write_text("def a():\n    return 1\n\n\ndef b():\n    return 2\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    out = _edit_file(perms, t, old_string="def c():\n    return 3",
                     new_string="def c():\n    return 4")
    assert "error" in out
    assert "not found" in out
    assert "Nothing was replaced" in out
    assert t.read_text() == "def a():\n    return 1\n\n\ndef b():\n    return 2\n"


def test_a_syntax_destroying_edit_is_applied_with_a_warning(workspace):
    """s6: the edit matches exactly but destroys the block's indentation.
    It is still written (a rewrite may pass through a broken state), and
    the result says the file no longer parses, with the line."""
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
    # parses.
    out = _edit_file(perms, t,
                     old_string="    return x",
                     new_string="  return x")
    assert "Edited" in out  # applied...
    assert "no longer parses" in out and "line 3" in out  # ...and said so
    import ast
    with pytest.raises(SyntaxError):
        ast.parse(t.read_text())


def test_multi_edit_is_atomic(workspace):
    """s9 said a failed batch left earlier edits applied. It does not
    (in-memory loop, single write at the end) -- this is the evidence."""
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


def test_multi_edit_names_the_failing_index(workspace):
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


def test_line_endings_survive_an_edit(workspace):
    """text_files.write_text_file(like=shape) already preserves CRLF; the
    engine must not regress this (edit point 5)."""
    t = workspace / "w.csv"
    t.write_bytes(b"alpha\r\nbeta\r\ngamma\r\n")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _read_file(perms, t)
    _edit_file(perms, t, old_string="beta", new_string="BETA")
    assert t.read_bytes() == b"alpha\r\nBETA\r\ngamma\r\n"
