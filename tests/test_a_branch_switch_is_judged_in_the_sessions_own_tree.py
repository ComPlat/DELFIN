"""A branch switch is routine only when it is one in the session's own tree.

Input: routine_actions._existing_branch_routine (GT1), reached through
action_classifier.screen. It looked for the branch and for a file of the
same name in the PROCESS's cwd, while the command runs in the session's
workspace. In the dashboard those differ: with a file ``main`` (edited,
uncommitted) in the workspace and a branch ``main`` in the cwd's repository,
``git checkout main`` was waved through, and in the workspace git reads it
as a path checkout and throws the edits away.

Now the classifier context carries the workspace; the branch and the file
are looked up there; without a workspace the shape is not routine.
"""
from __future__ import annotations

import subprocess

from delfin.agent import action_classifier as ac


def _repo(path, *, branch=None, file=None):
    path.mkdir()
    run = lambda *a: subprocess.run(["git", *a], cwd=path, check=True,
                                    capture_output=True)
    run("init", "-q", "-b", "main")
    run("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q",
        "--allow-empty", "-m", "x")
    if branch:
        run("branch", branch)
    if file:
        (path / file).write_text("edited\n")
    return path


def _verdict(cmd, workspace):
    ctx = {"task": "switch back"}
    if workspace is not None:
        ctx["workspace"] = str(workspace)
    return ac.screen({"tool": "bash", "command": cmd}, ctx)


def test_a_file_in_the_workspace_wins_over_a_branch_in_the_cwd(tmp_path,
                                                               monkeypatch):
    cwd_repo = _repo(tmp_path / "cwd", branch="notes")
    ws = _repo(tmp_path / "ws", file="notes")
    monkeypatch.chdir(cwd_repo)
    assert _verdict("git checkout notes", ws) is None


def test_a_branch_of_the_workspace_is_routine_from_anywhere(tmp_path,
                                                            monkeypatch):
    ws = _repo(tmp_path / "ws", branch="feature")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    v = _verdict("git switch feature", ws)
    assert v is not None and v.decision == "allow"


def test_without_a_workspace_a_switch_is_judged(tmp_path, monkeypatch):
    monkeypatch.chdir(_repo(tmp_path / "cwd", branch="feature"))
    assert _verdict("git switch feature", None) is None
