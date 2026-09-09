"""Three tools that reported something other than what happened.

Found by driving tools no model has ever called — 33 of the 72 advertised
have zero calls across 2509 recorded benchmark runs, so their tests are
the only thing that has ever exercised them, and a test written beside a
tool shares the blind spots of whoever wrote it.

**A grant outlived the thing it was granted for.** enter_worktree adds
the worktree to ``extra_workspace_dirs``, which every write gate treats
as inside the sandbox. Both ways out remove the directory — exit_worktree
tears it down, worktree_merge consumes it — and neither gave the root
back, because ``add_extra_dir`` had no inverse at all: the writable
roots of a session could only ever grow. So the ordinary path left a
session writing-enabled on a path that no longer existed.

**A dry run was recorded as a write.** apply_patch with check_only is
``git apply --check``: it names the files the diff WOULD touch and writes
none of them. The audit recorded one write per named file, so
list_changes_made — the one tool that answers "what did you do" — listed
a file that was never modified, and marked it "NOT undoable (no
pre-image)", which invites worry about an unrecoverable change that does
not exist.

**A watch was accepted on a job that cannot report.** watch_job took a
bash job id the registry does not know and answered ``status: watching``,
with a note telling the model to end its turn and wait for a notification
that can never arrive. bash_kill refuses the same id outright.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def git_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    (repo / "f.txt").write_text("x\n")
    subprocess.run(["git", "-C", str(repo), "add", "f.txt"], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.email=a@b",
                    "-c", "user.name=t", "commit", "-qm", "init"], check=True)
    return repo


def _perms(root, session="wt"):
    perms = KitToolPermissions(workspace=str(root))
    perms.mode = "acceptEdits"
    perms.task_session_id = session
    return perms


def _commit_in(worktree):
    (Path(worktree) / "f.txt").write_text("x\nchanged\n")
    subprocess.run(["git", "-C", worktree, "add", "f.txt"], check=True)
    subprocess.run(["git", "-C", worktree, "-c", "user.email=a@b",
                    "-c", "user.name=t", "commit", "-qm", "w"], check=True)


# ---------------------------------------------------------------------------
# The writable root
# ---------------------------------------------------------------------------

def test_exit_gives_the_writable_root_back(git_repo):
    perms = _perms(git_repo)
    ex = _DocToolExecutor()
    wt = json.loads(ex.execute("enter_worktree",
                               {"repo_dir": str(git_repo)}, perms))["path"]
    assert list(perms.extra_workspace_dirs), (
        "premise gone: entering no longer grants a root")
    ex.execute("exit_worktree", {"path": wt}, perms)
    assert not list(perms.extra_workspace_dirs), (
        "the session can still write to a worktree that is gone")


def test_merge_gives_the_writable_root_back_too(git_repo):
    """The likelier path: work, merge, carry on. It removes the worktree
    just as exit does."""
    perms = _perms(git_repo, "merge")
    ex = _DocToolExecutor()
    wt = json.loads(ex.execute("enter_worktree",
                               {"repo_dir": str(git_repo)}, perms))["path"]
    _commit_in(wt)
    ex.execute("worktree_merge",
               {"path": wt, "target_dir": str(git_repo)}, perms)
    assert not Path(wt).exists(), "premise gone: merge no longer consumes it"
    assert not list(perms.extra_workspace_dirs)


def test_a_root_still_in_use_is_never_taken_away(tmp_path):
    """The rule is one sentence: a grant for a directory that is gone is
    dead already. A directory that still exists is somebody's."""
    keep = tmp_path / "project"
    keep.mkdir()
    perms = _perms(tmp_path)
    perms.add_extra_dir(keep)
    assert keep.resolve() in perms.extra_workspace_dirs
    assert perms.drop_extra_dir(keep) is False
    assert keep.resolve() in perms.extra_workspace_dirs


def test_the_workspace_itself_is_never_dropped(tmp_path):
    perms = _perms(tmp_path)
    assert perms.drop_extra_dir(tmp_path) is False
    assert perms.workspace == Path(tmp_path).resolve()


def test_dropping_a_root_that_was_never_granted_is_a_no_op(tmp_path):
    perms = _perms(tmp_path)
    assert perms.drop_extra_dir(tmp_path / "never" / "granted") is False


# ---------------------------------------------------------------------------
# The dry run
# ---------------------------------------------------------------------------

def _own_lines(report: str, workspace) -> list[str]:
    """Only this test's paths — the audit log outlives one process."""
    return [ln.strip() for ln in report.splitlines() if str(workspace) in ln]


def test_a_check_only_patch_is_not_reported_as_a_change(tmp_path):
    target = tmp_path / "m.py"
    target.write_text("a = 1\n")
    perms = _perms(tmp_path, "check")
    ex = _DocToolExecutor()
    diff = "--- a/m.py\n+++ b/m.py\n@@ -1 +1 @@\n-a = 1\n+a = 2\n"

    out = json.loads(ex.execute(
        "apply_patch", {"diff": diff, "check_only": True}, perms))
    assert out["status"] == "ok" and out["files_touched"] == ["m.py"]
    assert target.read_text() == "a = 1\n", "a check wrote to the file"
    assert not _own_lines(ex.execute("list_changes_made", {}, perms), tmp_path), (
        "the changes report lists a file the check never modified")


def test_a_real_patch_is_still_reported(tmp_path):
    """The other half — a fix that silenced the audit entirely would be
    worse than the phantom write."""
    target = tmp_path / "m.py"
    target.write_text("a = 1\n")
    perms = _perms(tmp_path, "real")
    ex = _DocToolExecutor()
    ex.execute("apply_patch", {
        "diff": "--- a/m.py\n+++ b/m.py\n@@ -1 +1 @@\n-a = 1\n+a = 2\n"}, perms)
    assert target.read_text() == "a = 2\n"
    lines = _own_lines(ex.execute("list_changes_made", {}, perms), tmp_path)
    assert any("m.py" in ln and "apply_patch" in ln for ln in lines), lines


# ---------------------------------------------------------------------------
# The watch
# ---------------------------------------------------------------------------

def test_watching_an_unknown_bash_job_is_refused(tmp_path):
    perms = _perms(tmp_path, "watch")
    out = json.loads(_DocToolExecutor().execute(
        "watch_job", {"job_id": "does-not-exist"}, perms))
    assert "error" in out, out
    assert "unknown job_id" in out["error"]


def test_bash_kill_and_watch_job_agree_about_an_unknown_id(tmp_path):
    """They disagreed: one refused, the other answered 'watching'."""
    perms = _perms(tmp_path, "watch")
    ex = _DocToolExecutor()
    killed = json.loads(ex.execute("bash_kill", {"job_id": "nope"}, perms))
    watched = json.loads(ex.execute("watch_job", {"job_id": "nope"}, perms))
    assert killed.get("status") == "error"
    assert "error" in watched


def test_a_scheduler_id_is_still_accepted(tmp_path):
    """A SLURM id cannot be checked without asking the scheduler, so the
    refusal is deliberately limited to bash jobs this process owns."""
    perms = _perms(tmp_path, "watch")
    out = json.loads(_DocToolExecutor().execute(
        "watch_job", {"job_id": "123456"}, perms))
    assert out.get("status") == "watching", out
    assert out.get("kind") == "slurm"
