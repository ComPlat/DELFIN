"""worktree.diff_summary — review what an isolated subagent changed."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import worktree as wt


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=str(repo), check=True,
                   capture_output=True)


def _repo(tmp_path):
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    _git(r, "config", "user.email", "t@t")
    _git(r, "config", "user.name", "t")
    (r / "a.txt").write_text("hello\n")
    _git(r, "add", "-A")
    _git(r, "commit", "-qm", "init")
    return r


def test_diff_summary_shows_tracked_and_untracked(tmp_path):
    r = _repo(tmp_path)
    info = wt.enter_worktree(r)
    (Path(info.path) / "a.txt").write_text("hello\nworld\n")     # tracked edit
    (Path(info.path) / "new.txt").write_text("brand new\n")      # untracked
    info = wt.exit_worktree(info, keep_if_changed=True)
    assert info.had_changes
    ds = wt.diff_summary(info)
    assert "a.txt" in ds
    assert "new.txt" in ds and "Untracked" in ds


def test_diff_summary_empty_when_unchanged(tmp_path):
    r = _repo(tmp_path)
    info = wt.enter_worktree(r)
    info = wt.exit_worktree(info, keep_if_changed=True)
    assert not info.had_changes
    assert wt.diff_summary(info) == ""


def test_subagent_surfaces_diff_summary_wired():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "subagents.py").read_text(encoding="utf-8")
    assert "diff_summary" in src and 'worktree_summary["diff_summary"]' in src
