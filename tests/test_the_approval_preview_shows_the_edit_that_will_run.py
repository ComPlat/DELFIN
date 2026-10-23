"""The approval preview shows the edit that will actually be applied.

When old_string is not found verbatim, edit_file falls back to a
whitespace-tolerant match and applies it. The preview built for the
dialog replaced exactly and found nothing: on 2026-09-26 the
self-modification guard asked a person to approve "(no changes)" for an
edit that then rewrote engine.py. The preview now takes the same path.
"""

from __future__ import annotations

from delfin.agent.api_client import _DocToolExecutor

SOURCE = "def f():\n    if x:\n        return 1\n    return 2\n"


def _preview(tmp_path, old, new):
    path = tmp_path / "m.py"
    path.write_text(SOURCE)
    return _DocToolExecutor()._build_change_preview(
        "edit_file", {"path": str(path), "old_string": old,
                      "new_string": new}, path)


def test_a_whitespace_tolerant_edit_is_previewed(tmp_path):
    # Indentation differs from the file: the executor's fuzzy fallback
    # applies this; the dialog must show it.
    out = _preview(tmp_path, "if x:\n    return 1", "if x:\n    return 3")
    assert "(no changes)" not in out
    assert "+" in out and "return 3" in out


def test_an_exact_edit_is_previewed_as_before(tmp_path):
    out = _preview(tmp_path, "return 2", "return 4")
    assert "-    return 2" in out and "+    return 4" in out


def test_an_edit_that_matches_nothing_still_says_so(tmp_path):
    assert _preview(tmp_path, "no such text", "x") == "(no changes)"
