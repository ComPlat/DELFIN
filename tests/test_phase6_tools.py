"""Tests for Phase 6 tools: run_tests, apply_patch, find_definition/references."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.agent import code_nav as CN
from delfin.agent import patch_apply as PA
from delfin.agent import test_runner as TR
from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


# ---- run_tests ------------------------------------------------------------


def _make_minimal_repo() -> Path:
    """Create a tmp dir with a passing pytest test."""
    d = Path(tempfile.mkdtemp(prefix="run_tests_"))
    (d / "test_demo.py").write_text(
        "def test_addition():\n"
        "    assert 1 + 1 == 2\n"
        "def test_failing():\n"
        "    assert 'a' == 'b'\n",
        encoding="utf-8",
    )
    return d


def test_run_tests_reports_passes_and_fails():
    repo = _make_minimal_repo()
    try:
        result = TR.run_tests(repo, target="test_demo.py")
        assert result["status"] in ("failed", "ok")  # we expect failed
        assert result["framework"] == "pytest"
        s = result["summary"]
        assert s["passed"] == 1
        assert s["failed"] == 1
        # The failing test must show up in the failures list
        ids = [f["node_id"] for f in result["failures"]]
        assert any("test_failing" in i for i in ids)
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_run_tests_target_filter():
    repo = _make_minimal_repo()
    try:
        result = TR.run_tests(repo, target="test_demo.py::test_addition")
        s = result["summary"]
        assert s["passed"] == 1
        assert s["failed"] == 0
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_run_tests_missing_workspace():
    result = TR.run_tests("/no/such/path")
    assert result["status"] == "error"


def test_run_tests_tool_dispatch():
    repo = _make_minimal_repo()
    try:
        perms = KitToolPermissions(workspace=repo, mode="default")
        out = _DocToolExecutor().execute(
            "run_tests", {"target": "test_demo.py::test_addition"},
            permissions=perms,
        )
        payload = json.loads(out)
        assert payload["summary"]["passed"] == 1
    finally:
        shutil.rmtree(repo, ignore_errors=True)


# ---- apply_patch ----------------------------------------------------------


_SIMPLE_DIFF = """\
--- a/foo.txt
+++ b/foo.txt
@@ -1,3 +1,3 @@
 line1
-line2
+LINE2
 line3
"""


def test_patch_check_only_does_not_mutate():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "foo.txt").write_text("line1\nline2\nline3\n")
        result = PA.apply_patch(ws, _SIMPLE_DIFF, check_only=True, prefer="py")
        assert result["status"] == "ok"
        assert result["check_only"] is True
        # Content unchanged
        assert (ws / "foo.txt").read_text() == "line1\nline2\nline3\n"


def test_patch_applies_simple_change():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "foo.txt").write_text("line1\nline2\nline3\n")
        result = PA.apply_patch(ws, _SIMPLE_DIFF, prefer="py")
        assert result["status"] == "ok"
        assert (ws / "foo.txt").read_text() == "line1\nLINE2\nline3\n"
        assert "foo.txt" in result["files_touched"]


def test_patch_rejects_invalid_diff():
    with tempfile.TemporaryDirectory() as d:
        result = PA.apply_patch(d, "not a diff at all", prefer="py")
        assert result["status"] == "check_failed"


def test_patch_rejects_context_mismatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        # Wrong content — patch expects 'line2' but file has 'XXXXX'
        (ws / "foo.txt").write_text("line1\nXXXXX\nline3\n")
        result = PA.apply_patch(ws, _SIMPLE_DIFF, prefer="py")
        assert result["status"] == "check_failed"
        # Original unchanged
        assert "XXXXX" in (ws / "foo.txt").read_text()


def test_patch_atomic_failure_does_not_partially_write():
    """Two-file diff where the second hunk fails — first file untouched."""
    diff = """\
--- a/a.txt
+++ b/a.txt
@@ -1,1 +1,1 @@
-old_a
+new_a
--- a/b.txt
+++ b/b.txt
@@ -1,1 +1,1 @@
-MISSING
+new_b
"""
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "a.txt").write_text("old_a\n")
        (ws / "b.txt").write_text("real_b\n")
        result = PA.apply_patch(ws, diff, prefer="py")
        assert result["status"] == "check_failed"
        # Both files unchanged because we pre-checked all hunks
        assert (ws / "a.txt").read_text() == "old_a\n"
        assert (ws / "b.txt").read_text() == "real_b\n"


def test_patch_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "foo.txt").write_text("line1\nline2\nline3\n")
        perms = KitToolPermissions(workspace=ws, mode="default")
        out = _DocToolExecutor().execute(
            "apply_patch",
            {"diff": _SIMPLE_DIFF},
            permissions=perms,
        )
        payload = json.loads(out)
        assert payload["status"] == "ok"


def test_patch_blocked_in_plan_mode_unless_check_only():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "foo.txt").write_text("line1\nline2\nline3\n")
        perms = KitToolPermissions(workspace=ws, mode="plan")
        out = _DocToolExecutor().execute(
            "apply_patch", {"diff": _SIMPLE_DIFF}, permissions=perms,
        )
        payload = json.loads(out)
        assert "plan mode" in payload["error"]
        # check_only allowed
        out2 = _DocToolExecutor().execute(
            "apply_patch",
            {"diff": _SIMPLE_DIFF, "check_only": True},
            permissions=perms,
        )
        payload2 = json.loads(out2)
        assert payload2["status"] == "ok"


# ---- code_nav ------------------------------------------------------------


def test_find_definition_python():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "lib.py").write_text(
            "def compute_mean(values):\n"
            "    return sum(values) / len(values)\n"
            "\n"
            "def use_it():\n"
            "    return compute_mean([1, 2, 3])\n",
            encoding="utf-8",
        )
        result = CN.find_definition(ws, "compute_mean")
        assert result["matches"]
        first = result["matches"][0]
        assert first["path"] == "lib.py"
        assert first["line"] == 1
        assert "def compute_mean" in first["preview"]


def test_find_references_python():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "lib.py").write_text(
            "def compute_mean(values):\n"
            "    return sum(values) / len(values)\n"
            "\n"
            "def use_it():\n"
            "    return compute_mean([1, 2, 3])\n"
            "\n"
            "x = compute_mean([4, 5])\n",
            encoding="utf-8",
        )
        result = CN.find_references(ws, "compute_mean")
        # At least 3 hits: definition + 2 callsites
        assert len(result["matches"]) >= 3


def test_find_definition_invalid_symbol():
    with tempfile.TemporaryDirectory() as d:
        result = CN.find_definition(d, "not a valid identifier!")
        assert "error" in result


def test_find_definition_grep_fallback_for_js():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "app.js").write_text(
            "function helper() { return 1; }\n"
            "const x = helper() + helper();\n",
            encoding="utf-8",
        )
        result = CN.find_definition(ws, "helper", language="any")
        assert result["backend"] == "grep"
        # At least one match (the definition line)
        assert any("function helper" in m["preview"]
                   for m in result["matches"])


def test_find_definition_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "x.py").write_text("def foo(): return 1\n")
        perms = KitToolPermissions(workspace=ws, mode="default")
        out = _DocToolExecutor().execute(
            "find_definition", {"symbol": "foo"}, permissions=perms,
        )
        payload = json.loads(out)
        assert payload["matches"]


def test_find_references_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / "x.py").write_text(
            "def foo(): return 1\nbar = foo() + foo()\n"
        )
        perms = KitToolPermissions(workspace=ws, mode="default")
        out = _DocToolExecutor().execute(
            "find_references", {"symbol": "foo"}, permissions=perms,
        )
        payload = json.loads(out)
        assert len(payload["matches"]) >= 1
