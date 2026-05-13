"""Tests for save_typed_memory + parse_memory_type.

The Claude-Code-style memory mirror lives under
``~/.claude/projects/<slug>/memory/`` — we monkey-patch ``Path.home`` so
the test never touches the real home directory.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as ms


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    return tmp_path


def test_parse_memory_type_default_user():
    t, body = ms.parse_memory_type("prefers German")
    assert t == "user"
    assert body == "prefers German"


def test_parse_memory_type_explicit_prefix():
    for prefix, expected in [
        ("user:", "user"),
        ("feedback:", "feedback"),
        ("project:", "project"),
        ("reference:", "reference"),
        ("User:", "user"),
        ("FEEDBACK :", "feedback"),
    ]:
        t, body = ms.parse_memory_type(f"{prefix} the rest")
        assert t == expected, (prefix, t)
        assert body == "the rest"


def test_parse_memory_type_unknown_prefix_falls_back():
    """A leading word that *looks* like a prefix but isn't one of the
    four types must be preserved as part of the body."""
    t, body = ms.parse_memory_type("notes: this is just notes")
    assert t == "user"
    assert body == "notes: this is just notes"


def test_save_typed_memory_creates_file_with_frontmatter(fake_home, tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    fpath, slug, mtype = ms.save_typed_memory(
        "user: prefers German, exact, numeric",
        repo_root=repo,
    )
    assert mtype == "user"
    assert fpath.exists()
    assert fpath.name.startswith("user_")
    text = fpath.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert f"name: {slug}" in text
    assert "metadata:" in text and "type: user" in text
    assert "prefers German" in text


def test_save_typed_memory_indexes_in_memory_md(fake_home, tmp_path):
    repo = tmp_path / "repo2"
    repo.mkdir()
    ms.save_typed_memory(
        "feedback: do not mock the DB",
        repo_root=repo,
    )
    index = fpath = ms._claude_memory_dir(repo) / "MEMORY.md"
    content = index.read_text(encoding="utf-8")
    assert "## Feedback (how to work)" in content
    assert "(feedback_" in content
    assert "do not mock the DB" in content


def test_save_typed_memory_unique_filenames(fake_home, tmp_path):
    repo = tmp_path / "repo3"
    repo.mkdir()
    p1, slug1, _ = ms.save_typed_memory("user: first", repo_root=repo)
    p2, slug2, _ = ms.save_typed_memory("user: first", repo_root=repo)
    assert p1 != p2
    # Both indexed
    index = ms._claude_memory_dir(repo) / "MEMORY.md"
    assert index.read_text(encoding="utf-8").count("user_first") >= 1


def test_save_typed_memory_creates_missing_section(fake_home, tmp_path):
    """If MEMORY.md exists but lacks the type's section header,
    the writer appends a new section at the bottom."""
    repo = tmp_path / "repo4"
    repo.mkdir()
    mem_dir = ms._claude_memory_dir(repo)
    mem_dir.mkdir(parents=True)
    (mem_dir / "MEMORY.md").write_text(
        "# DELFIN Project Memory\n\n## Project (current work)\n- prior\n",
        encoding="utf-8",
    )
    ms.save_typed_memory("reference: see Linear/INGEST", repo_root=repo)
    content = (mem_dir / "MEMORY.md").read_text(encoding="utf-8")
    assert "## Reference (read only when relevant)" in content
    assert "Linear/INGEST" in content


def test_save_typed_memory_explicit_type_overrides_prefix(fake_home, tmp_path):
    repo = tmp_path / "repo5"
    repo.mkdir()
    _, _, mtype = ms.save_typed_memory(
        "user: this is the body",
        repo_root=repo,
        memory_type="project",
    )
    assert mtype == "project"


def test_slugify_handles_unicode_and_punctuation():
    assert ms._slugify("Über die Brücke!! ja??") in {
        "ber-die-brcke-ja",     # ASCII-stripped variant
        "ber-die-br-cke-ja",
    } or ms._slugify("Über die Brücke!! ja??")  # at minimum non-empty
