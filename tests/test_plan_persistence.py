"""Tests for plan-mode plan-file persistence via memory_store.save_plan."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as ms


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    return tmp_path


def test_save_plan_writes_file_with_frontmatter(fake_home, tmp_path):
    repo = tmp_path / "myrepo"
    repo.mkdir()
    body = (
        "# Plan: rewire flux capacitor\n\n"
        "## Context\nThe capacitor needs 1.21 GW.\n\n"
        "## Implementation\n1. Wire X\n2. Wire Y\n"
    )
    fpath = ms.save_plan(body, repo_root=repo)
    assert fpath.exists()
    assert fpath.suffix == ".md"
    text = fpath.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert "name:" in text
    assert "description:" in text
    assert "created_at:" in text
    assert "rewire flux capacitor" in text


def test_save_plan_lands_under_per_project_plans_dir(fake_home, tmp_path):
    repo = tmp_path / "repo2"
    repo.mkdir()
    fpath = ms.save_plan("# Tiny plan\n\nDo X.", repo_root=repo)
    plans_dir = ms._claude_plans_dir(repo)
    assert fpath.parent == plans_dir
    assert plans_dir.exists()


def test_save_plan_does_not_clobber_existing(fake_home, tmp_path):
    repo = tmp_path / "repo3"
    repo.mkdir()
    p1 = ms.save_plan("# Plan A\nA body", repo_root=repo)
    p2 = ms.save_plan("# Plan A\nA body again", repo_root=repo)
    assert p1 != p2
    assert p1.exists() and p2.exists()


def test_save_plan_uses_explicit_title_when_given(fake_home, tmp_path):
    repo = tmp_path / "repo4"
    repo.mkdir()
    fpath = ms.save_plan(
        "## Random heading\n\nthe rest",
        repo_root=repo,
        title="explicit title for slug",
    )
    assert "explicit-title-for-slug" in fpath.name


def test_save_plan_empty_body_raises(fake_home, tmp_path):
    repo = tmp_path / "repo5"
    repo.mkdir()
    with pytest.raises(ValueError):
        ms.save_plan("   \n   ", repo_root=repo)


def test_claude_plans_dir_distinct_from_memory_dir(fake_home, tmp_path):
    repo = tmp_path / "repo6"
    repo.mkdir()
    mem_dir = ms._claude_memory_dir(repo)
    plans_dir = ms._claude_plans_dir(repo)
    assert mem_dir != plans_dir
    assert mem_dir.parent == plans_dir.parent
    assert mem_dir.name == "memory"
    assert plans_dir.name == "plans"
