"""worktree.merge_worktree — safely bring a worktree's changes back.

Clean-apply-or-nothing: a clean patch lands in the target working tree; a
conflicting one leaves the target untouched.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

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


def test_merge_clean_applies_to_target(tmp_path):
    r = _repo(tmp_path)
    info = wt.enter_worktree(r)
    (Path(info.path) / "a.txt").write_text("hello\nworld\n")   # tracked edit
    (Path(info.path) / "new.txt").write_text("brand new\n")    # untracked add

    res = wt.merge_worktree(info)                              # target == source
    assert res.ok and res.applied
    assert "a.txt" in res.files and "new.txt" in res.files
    # Changes landed in the source working tree, uncommitted.
    assert (r / "a.txt").read_text() == "hello\nworld\n"
    assert (r / "new.txt").read_text() == "brand new\n"
    status = subprocess.run(["git", "status", "--porcelain"], cwd=str(r),
                            text=True, capture_output=True).stdout
    assert "new.txt" in status and "a.txt" in status


def test_merge_no_changes(tmp_path):
    r = _repo(tmp_path)
    info = wt.enter_worktree(r)
    res = wt.merge_worktree(info)
    assert res.ok and not res.applied
    assert "no changes" in res.message.lower()


def test_merge_conflict_leaves_target_untouched(tmp_path):
    r = _repo(tmp_path)
    info = wt.enter_worktree(r)
    (Path(info.path) / "a.txt").write_text("hello\nWORKTREE\n")  # worktree edit
    # Diverge the target working tree so the worktree patch cannot apply.
    (r / "a.txt").write_text("MAIN-ONLY\n")

    res = wt.merge_worktree(info)
    assert not res.ok and not res.applied
    assert "conflict" in res.message.lower() or "untouched" in res.message.lower()
    # Target file is preserved exactly — never half-merged.
    assert (r / "a.txt").read_text() == "MAIN-ONLY\n"


def test_merge_tool_registered():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert '"name": "worktree_merge"' in src
    assert 'name == "worktree_merge"' in src
    assert "_execute_worktree_merge" in src
