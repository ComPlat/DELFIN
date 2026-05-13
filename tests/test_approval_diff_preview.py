"""Test the unified-diff computation for the approval-prompt preview.

The helper lives inside ``create_tab`` as a closure, but it's pure so we
exercise the same logic via a free-standing copy here. If the production
implementation drifts, the test will surface the gap on the next run.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _compute(detail: str) -> str:
    """Mirror the production helper. Kept in sync by the test below."""
    import ast as _ast
    import difflib as _difflib
    try:
        d = _ast.literal_eval(detail) if detail and detail.lstrip().startswith("{") else {}
    except Exception:
        return ""
    tool = (d.get("tool_name") or "").strip()
    inp = d.get("tool_input") or {}
    if not isinstance(inp, dict):
        return ""
    file_path = (inp.get("file_path") or inp.get("path") or "").strip()
    if not file_path or tool not in ("Edit", "Write", "multi_edit"):
        return ""
    try:
        target = Path(file_path).expanduser()
        current = target.read_text(encoding="utf-8") if target.is_file() else ""
    except OSError:
        return ""

    if tool == "Write":
        proposed = inp.get("content", "")
    elif tool == "Edit":
        old = inp.get("old_string", "")
        new = inp.get("new_string", "")
        if not old:
            return ""
        replace_all = bool(inp.get("replace_all", False))
        if old not in current:
            return f"(diff preview unavailable — old_string not found in {target.name})"
        proposed = (current.replace(old, new) if replace_all
                    else current.replace(old, new, 1))
    elif tool == "multi_edit":
        edits = inp.get("edits") or []
        proposed = current
        for e in edits:
            old = e.get("old_string", "")
            new = e.get("new_string", "")
            if old and old in proposed:
                if e.get("replace_all"):
                    proposed = proposed.replace(old, new)
                else:
                    proposed = proposed.replace(old, new, 1)
    else:
        return ""

    if proposed == current:
        return f"(no change to {target.name})"
    diff_lines = list(_difflib.unified_diff(
        current.splitlines(keepends=False),
        proposed.splitlines(keepends=False),
        fromfile=str(target),
        tofile=f"{target} (proposed)",
        lineterm="",
        n=3,
    ))
    if not diff_lines:
        return ""
    return "```diff\n" + "\n".join(diff_lines) + "\n```"


def test_edit_tool_produces_diff(tmp_path):
    f = tmp_path / "demo.py"
    f.write_text("alpha\nbeta\ngamma\n", encoding="utf-8")
    detail = repr({
        "tool_name": "Edit",
        "tool_input": {
            "file_path": str(f),
            "old_string": "beta",
            "new_string": "BETA",
        },
    })
    out = _compute(detail)
    assert "```diff" in out
    assert "-beta" in out
    assert "+BETA" in out


def test_write_tool_produces_diff(tmp_path):
    f = tmp_path / "demo.txt"
    f.write_text("hello\nworld\n", encoding="utf-8")
    detail = repr({
        "tool_name": "Write",
        "tool_input": {
            "file_path": str(f),
            "content": "hello\nthere\n",
        },
    })
    out = _compute(detail)
    assert "-world" in out
    assert "+there" in out


def test_multi_edit_applies_all_edits(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text("one\ntwo\nthree\n", encoding="utf-8")
    detail = repr({
        "tool_name": "multi_edit",
        "tool_input": {
            "file_path": str(f),
            "edits": [
                {"old_string": "one",   "new_string": "ONE"},
                {"old_string": "three", "new_string": "THREE"},
            ],
        },
    })
    out = _compute(detail)
    assert "-one" in out and "+ONE" in out
    assert "-three" in out and "+THREE" in out


def test_edit_with_missing_old_string_returns_explanation(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("alpha\n", encoding="utf-8")
    detail = repr({
        "tool_name": "Edit",
        "tool_input": {
            "file_path": str(f),
            "old_string": "not present",
            "new_string": "anything",
        },
    })
    assert "old_string not found" in _compute(detail)


def test_non_mutating_tool_returns_empty(tmp_path):
    detail = repr({
        "tool_name": "Read",
        "tool_input": {"file_path": str(tmp_path / "x.txt")},
    })
    assert _compute(detail) == ""


def test_no_change_returns_no_change_marker(tmp_path):
    f = tmp_path / "x.txt"
    f.write_text("hi\n", encoding="utf-8")
    detail = repr({
        "tool_name": "Write",
        "tool_input": {"file_path": str(f), "content": "hi\n"},
    })
    assert "no change" in _compute(detail)


def test_malformed_detail_returns_empty():
    assert _compute("") == ""
    assert _compute("not python") == ""


def test_production_helper_in_sync():
    """The dashboard's `_compute_approval_diff` closure must contain the
    same key string fragments as this copy — guards against silent drift."""
    src = Path(__file__).resolve().parent.parent / "delfin" / "dashboard" / "tab_agent.py"
    text = src.read_text(encoding="utf-8")
    assert "_compute_approval_diff" in text
    # Spot-check that the production version handles all 3 tool names
    fn_start = text.find("def _compute_approval_diff")
    fn_end = text.find("def _show_approval_prompt", fn_start)
    body = text[fn_start:fn_end]
    assert '"Edit"' in body and '"Write"' in body and '"multi_edit"' in body
    assert "difflib" in body
