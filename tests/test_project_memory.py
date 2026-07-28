"""Tests for delfin.agent.project_memory (DELFIN.MD auto-loader)."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent.project_memory import (
    discover_project_memory_files,
    load_project_memory,
)


def test_discover_finds_delfin_md_in_cwd():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "DELFIN.MD").write_text("hello", encoding="utf-8")
        files = discover_project_memory_files(root)
        assert any(p.name == "DELFIN.MD" and p.parent == root for p in files)


def test_discover_walks_up_to_parent():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        sub = root / "a" / "b"
        sub.mkdir(parents=True)
        (root / "AGENTS.md").write_text("root rules", encoding="utf-8")
        (sub / "DELFIN.MD").write_text("sub rules", encoding="utf-8")
        files = discover_project_memory_files(sub)
        names_paths = [(p.parent, p.name) for p in files]
        assert (sub, "DELFIN.MD") in names_paths
        assert (root, "AGENTS.md") in names_paths
        # closest first
        assert names_paths.index((sub, "DELFIN.MD")) < names_paths.index((root, "AGENTS.md"))


def test_load_concatenates_with_headers():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        sub = root / "x"
        sub.mkdir()
        (root / "DELFIN.MD").write_text("ROOTRULE", encoding="utf-8")
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
        (root / "DELFIN.MD").write_text(big, encoding="utf-8")
        text = load_project_memory(root, max_chars=2000, per_file_cap=1500)
        assert "[... truncated ...]" in text
        assert len(text) <= 3000


def test_load_with_extra_roots():
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        r1 = Path(d1).resolve()
        r2 = Path(d2).resolve()
        (r1 / "DELFIN.MD").write_text("primary", encoding="utf-8")
        (r2 / "AGENTS.md").write_text("extra", encoding="utf-8")
        text = load_project_memory(r1, extra_roots=[r2])
        assert "primary" in text
        assert "extra" in text


def test_recognises_all_delfin_filenames():
    for name in ("DELFIN.MD", "DELFIN.md", "AGENTS.md"):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            (root / name).write_text(f"FROM-{name}", encoding="utf-8")
            text = load_project_memory(root)
            assert f"FROM-{name}" in text


def test_ignores_other_frameworks_instruction_files():
    """DELFIN is standalone: instruction files of other agent tools must
    NOT leak into its prompt."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        (root / "CLAUDE.md").write_text("FOREIGN-RULES", encoding="utf-8")
        (root / "GEMINI.md").write_text("FOREIGN-RULES", encoding="utf-8")
        assert load_project_memory(root) == ""
        (root / "DELFIN.MD").write_text("OWN-RULES", encoding="utf-8")
        text = load_project_memory(root)
        assert "OWN-RULES" in text
        assert "FOREIGN-RULES" not in text


def test_no_duplicate_when_both_case_variants_resolve_to_same_file():
    """On a case-insensitive filesystem DELFIN.MD and DELFIN.md are one
    file — the loader must not inject it twice (inode dedupe)."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        (root / "DELFIN.MD").write_text("ONCE", encoding="utf-8")
        files = discover_project_memory_files(root)
        idents = {p.stat().st_ino for p in files}
        assert len(files) == len(idents)


def test_skips_empty_files():
    with tempfile.TemporaryDirectory() as d:
        root = Path(d).resolve()
        (root / "DELFIN.MD").write_text("", encoding="utf-8")
        (root / "AGENTS.md").write_text("real", encoding="utf-8")
        text = load_project_memory(root)
        assert "real" in text
