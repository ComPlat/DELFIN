"""`enter_worktree` turned a refusal into a writable copy.

The tool takes a `repo_dir` from the model, creates a git worktree of it,
and registers that worktree as a writable root so later edits and shell
calls land without a separate grant. Every part of that is right for the
repository the user is working in. `repo_dir` was not checked against
anything.

So a path the read gate refuses — outside the workspace, no confirm
callback, "add the directory via remember_permission" — became reachable
in one call: the worktree holds the same files, the tool grants write
access to it, and `worktree_merge` puts changes back into the original
working tree. Measured 2026-09-08 with two scratch repositories: reading
`elsewhere/crown_jewels.py` was denied before the call and allowed after
it.

The other two verbs have the same shape. `worktree_merge` takes a
`target_dir` it writes INTO, and `exit_worktree` takes a path it removes.
Both were unchecked.

Containment, not confirmation: a worktree of a repository the user
already granted is exactly what the tool is for, and the grant that
follows is the point of it. What it may not do is choose a different
repository.
"""

from __future__ import annotations

import json
import subprocess

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _repo(path):
    path.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    (path / "crown_jewels.py").write_text("secret business logic\n")
    subprocess.run(["git", "-C", str(path), "add", "."], check=True)
    subprocess.run(["git", "-C", str(path), "-c", "user.email=t@t",
                    "-c", "user.name=t", "commit", "-qm", "init"], check=True)
    return path


@pytest.fixture
def two_repos(tmp_path):
    ws = _repo(tmp_path / "workspace")
    other = _repo(tmp_path / "elsewhere")
    perms = KitToolPermissions(workspace=ws)
    return _DocToolExecutor(), perms, ws, other


def test_a_repository_outside_the_workspace_is_refused(two_repos):
    ex, perms, ws, other = two_repos
    out = json.loads(ex._execute_enter_worktree({"repo_dir": str(other)}, perms))
    assert "error" in out, out
    assert str(other) in out["error"]
    assert perms.extra_workspace_dirs == (), (
        "a writable root was granted for a repository that was refused")


def test_the_refusal_does_not_leave_a_worktree_behind(two_repos):
    ex, perms, ws, other = two_repos
    ex._execute_enter_worktree({"repo_dir": str(other)}, perms)
    listing = subprocess.check_output(
        ["git", "-C", str(other), "worktree", "list"], text=True)
    assert listing.strip().count("\n") == 0, listing


def test_the_workspace_itself_still_works(two_repos):
    ex, perms, ws, other = two_repos
    out = json.loads(ex._execute_enter_worktree({}, perms))
    assert out.get("status") == "ok", out
    assert out["path"] in [str(p) for p in perms.extra_workspace_dirs]


def test_naming_the_workspace_explicitly_still_works(two_repos):
    ex, perms, ws, other = two_repos
    out = json.loads(ex._execute_enter_worktree({"repo_dir": str(ws)}, perms))
    assert out.get("status") == "ok", out


def test_a_granted_directory_is_a_repository_it_may_use(two_repos):
    """The check is containment, not identity: once the user has granted
    a directory, a worktree of it is ordinary work."""
    ex, perms, ws, other = two_repos
    perms.add_extra_dir(other)
    out = json.loads(ex._execute_enter_worktree({"repo_dir": str(other)}, perms))
    assert out.get("status") == "ok", out


def test_a_merge_cannot_choose_where_the_changes_land(two_repos):
    ex, perms, ws, other = two_repos
    made = json.loads(ex._execute_enter_worktree({}, perms))
    assert made.get("status") == "ok", made
    out = json.loads(ex._execute_worktree_merge(
        {"path": made["path"], "target_dir": str(other)}, perms))
    assert "error" in out, out
    assert str(other) in out["error"]


def test_a_merge_cannot_take_a_worktree_from_anywhere(two_repos):
    ex, perms, ws, other = two_repos
    out = json.loads(ex._execute_worktree_merge({"path": str(other)}, perms))
    assert "error" in out, out


def test_leaving_cannot_remove_a_directory_it_was_never_given(two_repos):
    ex, perms, ws, other = two_repos
    out = json.loads(ex._execute_exit_worktree({"path": str(other)}, perms))
    assert "error" in out, out
    assert (other / "crown_jewels.py").is_file()


def test_without_a_permissions_object_nothing_is_derived_from_a_workspace(
        two_repos):
    """No permissions means no workspace to be contained by. The tool
    already required an explicit repo_dir there; it must not become the
    way round the check."""
    ex, _perms, ws, other = two_repos
    out = json.loads(ex._execute_enter_worktree({"repo_dir": str(other)}, None))
    assert "error" in out, out
