"""A worktree merge writes only what a write_file would be allowed to write.

Found running wave 14 (2026-10-09): worktree_merge applied a worktree's
diff into the target with git, file by file unchecked. A file the session
may not write -- outside its write scope, or a protected one behind the
Self-Modification Guard -- could be edited in a worktree (which carries no
write scope of its own) and merged back. Every path the merge would add,
change, delete or rename now goes through the target's write gate first; one
refusal and nothing is applied.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _git(cwd, *args):
    subprocess.run(["git", *args], cwd=str(cwd), check=True,
                   capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    (r / "allowed").mkdir(parents=True)
    (r / "other").mkdir()
    (r / "delfin" / "agent").mkdir(parents=True)
    (r / "allowed" / "a.txt").write_text("a\n")
    (r / "other" / "b.txt").write_text("b\n")
    (r / "delfin" / "agent" / "api_client.py").write_text("x = 1\n")
    _git(r, "init", "-q")
    _git(r, "-c", "user.email=t@t", "-c", "user.name=t", "add", "-A")
    _git(r, "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "base")
    wt = tmp_path / "wt"
    _git(r, "worktree", "add", "-q", "-b", "side", str(wt))
    return r, wt


def _merge(repo, wt, **perm_kw):
    perms = KitToolPermissions(workspace=str(repo), mode="bypassPermissions",
                               extra_workspace_dirs=(str(wt),), **perm_kw)
    out = _DocToolExecutor()._execute_worktree_merge({"path": str(wt)}, perms)
    return json.loads(out)


def test_a_file_outside_the_write_scope_is_not_merged(repo):
    r, wt = repo
    (wt / "other" / "b.txt").write_text("changed outside the scope\n")
    out = _merge(r, wt, write_allow_globs=("allowed/*",))
    assert out.get("applied") is not True, out
    assert (r / "other" / "b.txt").read_text() == "b\n"


def test_a_protected_file_is_not_merged_without_an_approval(repo):
    r, wt = repo
    (wt / "delfin" / "agent" / "api_client.py").write_text("x = 2  # sneaked in\n")
    out = _merge(r, wt)
    assert out.get("applied") is not True, out
    assert (r / "delfin" / "agent" / "api_client.py").read_text() == "x = 1\n"


def test_a_deletion_outside_the_scope_is_not_merged(repo):
    r, wt = repo
    (wt / "other" / "b.txt").unlink()
    out = _merge(r, wt, write_allow_globs=("allowed/*",))
    assert out.get("applied") is not True, out
    assert (r / "other" / "b.txt").exists()


def test_an_allowed_change_still_merges(repo):
    r, wt = repo
    (wt / "allowed" / "a.txt").write_text("a2\n")
    out = _merge(r, wt, write_allow_globs=("allowed/*",))
    assert out.get("applied") is True, out
    assert (r / "allowed" / "a.txt").read_text() == "a2\n"
