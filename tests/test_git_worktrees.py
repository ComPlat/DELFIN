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


# ---------------------------------------------------------------------------
# A4 — grouping by path prefix
# ---------------------------------------------------------------------------

def _make(path: str, *, is_main: bool = False, is_bare: bool = False):
    return git_worktrees.WorktreeEntry(
        path=path, head="abc1234", branch="GUPPY",
        is_main=is_main, is_bare=is_bare,
    )


def test_classify_main_worktree():
    e = _make("/home/u/repo", is_main=True)
    assert git_worktrees._classify_worktree(e) == "main"


def test_classify_bare_repo_is_main_bucket():
    e = _make("/home/u/repo.git", is_bare=True)
    assert git_worktrees._classify_worktree(e) == "main"


def test_classify_commit_sweep():
    e = _make("/home/u/agent_workspace/commit_sweep/worktrees/abc123")
    assert git_worktrees._classify_worktree(e) == "commit_sweep"


def test_classify_bisect_variants():
    bisect_paths = [
        "/tmp/bvs_abc", "/tmp/vf_abc", "/tmp/vsf_abc", "/tmp/vs51_abc",
        "/tmp/sweep_wt_abc", "/tmp/sd_abc", "/tmp/sw_abc",
        "/tmp/r20_abc", "/tmp/r20m_abc", "/tmp/feschiff_abc",
        "/tmp/main_wt", "/tmp/wt_pre_failfast", "/tmp/delfin_bisect_wt",
        "/tmp/delfin-main-cherry",
    ]
    for p in bisect_paths:
        assert git_worktrees._classify_worktree(_make(p)) == "bisect", p


def test_classify_unknown_path_is_other():
    e = _make("/home/u/some/random/spot")
    assert git_worktrees._classify_worktree(e) == "other"


def test_group_worktrees_partitions_correctly():
    items = [
        _make("/home/u/repo", is_main=True),
        _make("/home/u/agent_workspace/commit_sweep/worktrees/a"),
        _make("/home/u/agent_workspace/commit_sweep/worktrees/b"),
        _make("/tmp/bvs_xyz"),
        _make("/tmp/vf_abc"),
        _make("/home/u/random"),
    ]
    g = git_worktrees._group_worktrees(items)
    assert len(g["main"]) == 1
    assert len(g["commit_sweep"]) == 2
    assert len(g["bisect"]) == 2
    assert len(g["other"]) == 1


def test_render_uses_grouped_layout():
    items = [
        _make("/home/u/repo", is_main=True),
        _make("/home/u/agent_workspace/commit_sweep/worktrees/a"),
        _make("/tmp/bvs_xyz"),
    ]
    html = git_worktrees.render_worktrees_html(items)
    assert "1 main" in html
    assert "1 commit_sweep" in html
    assert "1 bisect" in html
    assert "<details" in html
    assert "Main worktrees" in html
    assert "External: commit_sweep" in html
    assert "External: bisect" in html


def test_render_main_group_open_external_groups_collapsed():
    items = [
        _make("/home/u/repo", is_main=True),
        _make("/tmp/bvs_xyz"),
    ]
    html = git_worktrees.render_worktrees_html(items)
    # Walk back from each section title to its <details ...> open tag.
    main_idx = html.find("Main worktrees")
    bisect_idx = html.find("External: bisect")
    main_open_idx = html.rfind("<details", 0, main_idx)
    bisect_open_idx = html.rfind("<details", 0, bisect_idx)
    main_tag = html[main_open_idx:html.find(">", main_open_idx) + 1]
    bisect_tag = html[bisect_open_idx:html.find(">", bisect_open_idx) + 1]
    assert " open" in main_tag, main_tag
    assert " open" not in bisect_tag, bisect_tag


def test_render_skips_empty_groups():
    items = [_make("/home/u/repo", is_main=True)]
    html = git_worktrees.render_worktrees_html(items)
    assert "Main worktrees" in html
    assert "commit_sweep" not in html
    assert "External: bisect" not in html
