"""Tests for delfin.agent.project_memory (CLAUDE.md auto-loader)."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent.project_memory import (
    discover_project_memory_files,
    load_project_memory,
)


def test_discover_finds_claude_md_in_cwd():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "CLAUDE.md").write_text("hello", encoding="utf-8")
        files = discover_project_memory_files(root)
        assert any(p.name == "CLAUDE.md" and p.parent == root for p in files)


def test_discover_walks_up_to_parent():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        sub = root / "a" / "b"
        sub.mkdir(parents=True)
        (root / "AGENTS.md").write_text("root rules", encoding="utf-8")
        (sub / "CLAUDE.md").write_text("sub rules", encoding="utf-8")
        files = discover_project_memory_files(sub)
        names_paths = [(p.parent, p.name) for p in files]
        assert (sub, "CLAUDE.md") in names_paths
        assert (root, "AGENTS.md") in names_paths
        # closest first
        assert names_paths.index((sub, "CLAUDE.md")) < names_paths.index((root, "AGENTS.md"))


def test_load_concatenates_with_headers():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        sub = root / "x"
        sub.mkdir()
        (root / "CLAUDE.md").write_text("ROOTRULE", encoding="utf-8")
        (sub / "DELFIN.md").write_text("SUBRULE", encoding="utf-8")
        text = load_project_memory(sub)
        assert "ROOTRULE" in text
        assert "SUBRULE" in text
        assert "Project memory:" in text


def test_load_empty_when_no_files():
    with tempfile.TemporaryDirectory() as d:
        text = load_project_memory(d)
        assert text == ""


def test_load_truncates_long_files():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        big = "x" * 50_000
        (root / "CLAUDE.md").write_text(big, encoding="utf-8")
        text = load_project_memory(root, max_chars=2000, per_file_cap=1500)
        assert "[... truncated ...]" in text
        assert len(text) <= 3000


def test_load_with_extra_roots():
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        r1 = Path(d1).resolve()
        r2 = Path(d2).resolve()
        (r1 / "CLAUDE.md").write_text("primary", encoding="utf-8")
        (r2 / "AGENTS.md").write_text("extra", encoding="utf-8")
        text = load_project_memory(r1, extra_roots=[r2])
        assert "primary" in text
        assert "extra" in text


def test_recognises_all_three_filenames():
    for name in ("CLAUDE.md", "AGENTS.md", "DELFIN.md"):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            (root / name).write_text(f"FROM-{name}", encoding="utf-8")
            text = load_project_memory(root)
            assert f"FROM-{name}" in text


def test_skips_empty_files():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        (root / "CLAUDE.md").write_text("", encoding="utf-8")
        (root / "AGENTS.md").write_text("real", encoding="utf-8")
        text = load_project_memory(root)
        assert "real" in text
