"""Tests for delfin.agent.worktree (git worktree isolation)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.agent import worktree as W
from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _git_init(path: Path) -> None:
    """Init a fresh git repo with one commit."""
    subprocess.run(["git", "init", "--quiet", "--initial-branch=main"],
                   cwd=str(path), check=True)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=str(path), check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=str(path), check=True)
    subprocess.run(["git", "config", "commit.gpgsign", "false"],
                   cwd=str(path), check=True)
    (path / "README.md").write_text("hi\n")
    subprocess.run(["git", "add", "."], cwd=str(path), check=True)
    subprocess.run(["git", "commit", "--quiet", "-m", "init"],
                   cwd=str(path), check=True)


def _git_available() -> bool:
    return shutil.which("git") is not None


pytestmark = pytest.mark.skipif(
    not _git_available(), reason="git not available",
)


def test_enter_creates_worktree():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        info = W.enter_worktree(repo)
        try:
            assert info.path.is_dir()
            assert info.branch.startswith("agent/")
            assert (info.path / "README.md").is_file()
        finally:
            W.exit_worktree(info, keep_if_changed=False)


def test_exit_cleans_when_no_changes():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        info = W.enter_worktree(repo)
        wt_path = info.path
        W.exit_worktree(info, keep_if_changed=False)
        assert not wt_path.is_dir()
        assert info.cleaned_up
        assert not info.had_changes
        # branch should be gone too
        out = subprocess.check_output(
            ["git", "branch", "--list"], cwd=str(repo), text=True,
        )
        assert info.branch not in out


def test_exit_keeps_when_changes_detected():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        info = W.enter_worktree(repo)
        # introduce a change
        (info.path / "new.txt").write_text("changed")
        W.exit_worktree(info, keep_if_changed=True)
        assert info.had_changes
        assert info.final_path is not None


def test_session_context_manager():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        with W.worktree_session(repo, keep_if_changed=False) as info:
            assert info.path.is_dir()
            wt = info.path
        assert not wt.is_dir()


def test_enter_rejects_non_repo():
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(W.WorktreeError):
            W.enter_worktree(d)


def test_enter_worktree_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        perms = KitToolPermissions(workspace=repo, mode="default")
        out = _DocToolExecutor().execute(
            "enter_worktree", {"repo_dir": str(repo)}, permissions=perms,
        )
        payload = json.loads(out)
        try:
            assert payload["status"] == "ok"
            assert payload["branch"].startswith("agent/")
            wt_path = Path(payload["path"])
            assert wt_path.is_dir()
            # The worktree must be auto-registered as an extra dir
            assert wt_path in perms.extra_workspace_dirs
        finally:
            _DocToolExecutor().execute(
                "exit_worktree",
                {"path": payload["path"], "keep_if_changed": False},
                permissions=perms,
            )


def test_exit_worktree_tool_dispatch():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d) / "repo"; repo.mkdir()
        _git_init(repo)
        perms = KitToolPermissions(workspace=repo, mode="default")
        out = _DocToolExecutor().execute(
            "enter_worktree", {"repo_dir": str(repo)}, permissions=perms,
        )
        wt_path = json.loads(out)["path"]
        out2 = _DocToolExecutor().execute(
            "exit_worktree",
            {"path": wt_path, "keep_if_changed": False},
            permissions=perms,
        )
        payload = json.loads(out2)
        assert payload["status"] == "ok"
        assert not payload["had_changes"]
        assert not payload["kept"]


def test_exit_worktree_missing_path():
    perms = KitToolPermissions(
        workspace=Path("/tmp"), mode="default",
    )
    out = _DocToolExecutor().execute(
        "exit_worktree", {"path": "/no/such/dir"}, permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload
