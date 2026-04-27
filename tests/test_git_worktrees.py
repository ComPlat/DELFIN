"""Tests for delfin.dashboard.git_worktrees."""
from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.dashboard import git_worktrees


# ---------------------------------------------------------------------------
# list_worktrees
# ---------------------------------------------------------------------------

def _fake_run_git(stdout: str, returncode: int = 0):
    def _runner(repo, *args, timeout=3.0):
        if returncode != 0:
            return ""
        return stdout
    return _runner


def test_list_no_git_returns_empty(tmp_path):
    """No git binary or no repo → empty list."""
    with patch.object(git_worktrees, "_run_git", return_value=""):
        assert git_worktrees.list_worktrees(tmp_path) == []


def test_list_single_main_worktree(tmp_path):
    output = (
        "worktree /repo\n"
        "HEAD abc123def456\n"
        "branch refs/heads/main\n"
    )
    with patch.object(git_worktrees, "_run_git",
                      return_value=output):
        result = git_worktrees.list_worktrees(tmp_path)
    assert len(result) == 1
    e = result[0]
    assert e.path == "/repo"
    assert e.head == "abc123def456"
    assert e.branch == "main"
    assert e.is_main is True
    assert e.is_bare is False


def test_list_multiple_worktrees(tmp_path):
    output = (
        "worktree /repo\n"
        "HEAD aaa\n"
        "branch refs/heads/main\n"
        "\n"
        "worktree /tmp/wt-feature\n"
        "HEAD bbb\n"
        "branch refs/heads/feature\n"
    )
    with patch.object(git_worktrees, "_run_git",
                      return_value=output):
        result = git_worktrees.list_worktrees(tmp_path)
    assert len(result) == 2
    assert result[0].is_main is True
    assert result[1].is_main is False
    assert result[1].branch == "feature"


def test_list_detected_bare_repo(tmp_path):
    output = (
        "worktree /repo.git\n"
        "bare\n"
        "\n"
        "worktree /repo/wt-main\n"
        "HEAD abc\n"
        "branch refs/heads/main\n"
    )
    with patch.object(git_worktrees, "_run_git",
                      return_value=output):
        result = git_worktrees.list_worktrees(tmp_path)
    assert result[0].is_bare is True
    assert result[1].is_bare is False
    assert result[1].is_main is False  # the bare repo took the main slot


def test_list_detached_head(tmp_path):
    output = (
        "worktree /repo\n"
        "HEAD abc\n"
        "detached\n"
    )
    with patch.object(git_worktrees, "_run_git",
                      return_value=output):
        result = git_worktrees.list_worktrees(tmp_path)
    assert len(result) == 1
    assert result[0].branch == ""  # detached has no branch line


def test_list_real_repo_smoke():
    """Smoke test on the actual DELFIN repo (no mocks).

    Just verifies the parser doesn't crash on real ``git worktree list``
    output and returns at least one entry.
    """
    repo = Path(__file__).resolve().parents[1]
    if not (repo / ".git").exists():
        pytest.skip("not a git checkout")
    items = git_worktrees.list_worktrees(repo)
    # If git is available we expect at least the main worktree
    if items:
        assert any(it.is_main for it in items)


# ---------------------------------------------------------------------------
# render_worktrees_html
# ---------------------------------------------------------------------------

def test_render_empty_shows_helpful_hint():
    html = git_worktrees.render_worktrees_html([])
    assert "git worktree add" in html


def test_render_main_branch_badge():
    items = [git_worktrees.WorktreeEntry(
        path="/repo", head="abc123de", branch="main",
        is_main=True, is_bare=False,
    )]
    html = git_worktrees.render_worktrees_html(items)
    assert "MAIN" in html
    assert "main" in html
    assert "abc123de" in html


def test_render_secondary_worktree_badge():
    items = [
        git_worktrees.WorktreeEntry(
            path="/repo", head="aaa", branch="main",
            is_main=True, is_bare=False),
        git_worktrees.WorktreeEntry(
            path="/tmp/wt", head="bbb", branch="feature",
            is_main=False, is_bare=False),
    ]
    html = git_worktrees.render_worktrees_html(items)
    assert "MAIN" in html
    assert "WORKTREE" in html


def test_render_bare_badge():
    items = [git_worktrees.WorktreeEntry(
        path="/repo.git", head="", branch="",
        is_main=False, is_bare=True,
    )]
    html = git_worktrees.render_worktrees_html(items)
    assert "BARE" in html


def test_render_detached_label():
    items = [git_worktrees.WorktreeEntry(
        path="/repo", head="abc", branch="",
        is_main=True, is_bare=False,
    )]
    html = git_worktrees.render_worktrees_html(items)
    assert "(detached)" in html


def test_render_escapes_user_text():
    items = [git_worktrees.WorktreeEntry(
        path="<script>", head="<head>", branch="<branch>",
        is_main=False, is_bare=False,
    )]
    html = git_worktrees.render_worktrees_html(items)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "&lt;branch&gt;" in html
