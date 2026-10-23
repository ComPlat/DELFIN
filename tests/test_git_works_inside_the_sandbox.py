"""A session in its own worktree can still use git.

The isolation confines a command to its workspace. A git worktree checks
out under ``<repo>/.delfin/worktrees/<name>`` but keeps its metadata and
object store under the main repo's ``.git`` -- outside the workspace -- so
every git command failed with "not a git repository" (driven 2026-09-17).
The sandbox now grants git's own directories, so commit and branch
work, while a hook, the repository config and info/ stay read-only (bubblewrap and
Seatbelt) or caught by the command gate.
"""
import json
import os
import subprocess
import sys

import pytest

from delfin.agent import api_client as A
from delfin.agent import socket_guard as SG


@pytest.fixture
def worktree(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    env = {**os.environ, "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null"}
    run = lambda *a: subprocess.run(["git", *a], cwd=repo, env=env, check=True,
                                    capture_output=True)
    run("init", "-q")
    run("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "--allow-empty", "-m", "init")
    run("worktree", "add", "-q", ".delfin/worktrees/s1", "-b", "session-s1")
    return repo, (repo / ".delfin" / "worktrees" / "s1")


def test_the_git_roots_are_the_common_dir_and_the_worktree_gitdir(worktree):
    repo, ws = worktree
    roots = A._git_metadata_roots(str(ws))
    assert str((repo / ".git").resolve()) in roots
    assert str((repo / ".git" / "worktrees" / "s1").resolve()) in roots


def test_a_plain_directory_has_no_git_roots(tmp_path):
    (tmp_path / "plain").mkdir()
    assert A._git_metadata_roots(str(tmp_path / "plain")) == []


def _commit(ws, perms):
    cmd = ("printf hi > f.txt && git add f.txt && "
           "git -c user.email=t@t -c user.name=t commit -q -m work && echo COMMIT-OK; "
           "git checkout -q -b feature && echo BRANCH-OK")
    argv = A._bash_isolation_argv(cmd, str(ws), perms)
    return subprocess.run(argv, cwd=str(ws), capture_output=True, text=True, timeout=90)


@pytest.mark.skipif(not A._bwrap_functional(), reason="needs bubblewrap")
def test_git_commits_under_bubblewrap(worktree, monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    _repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    out = _commit(ws, perms)
    assert "COMMIT-OK" in out.stdout and "BRANCH-OK" in out.stdout, (out.stdout, out.stderr)


@pytest.mark.skipif(not (A._landlock_functional() and SG.available()),
                    reason="needs Landlock and the socket guard")
def test_git_commits_under_landlock(worktree, monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    _repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    out = _commit(ws, perms)
    assert "COMMIT-OK" in out.stdout and "BRANCH-OK" in out.stdout, (out.stdout, out.stderr)


def test_bubblewrap_binds_the_hooks_dir_read_only(worktree, monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    argv = A._bash_isolation_argv("git status", str(ws), perms)
    hooks = str((repo / ".git" / "hooks").resolve())
    assert "--ro-bind-try" in argv
    i = argv.index(hooks)
    assert argv[i - 1] == "--ro-bind-try"
    # bound writable first, re-bound read-only after
    assert argv.index("--bind", 0, i) < i


def test_the_command_gate_refuses_a_literal_hook_write(worktree, monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    _repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), mode="bypassPermissions")
    out = json.loads(A._doc_executor._execute_bash(
        {"command": "echo x > .git/hooks/pre-commit"}, perms))
    assert "error" in out and ".git/hooks" in out["error"]


def test_git_config_is_not_hidden_from_the_sandbox():
    from delfin.agent.sandbox import _HOME_SECRET_DIRS
    assert ".config/git" not in _HOME_SECRET_DIRS
    assert ".config/gh" in _HOME_SECRET_DIRS      # the token store stays hidden


@pytest.mark.parametrize("name", ["config", "config.worktree", "info"])
def test_bubblewrap_binds_the_git_config_read_only(worktree, monkeypatch, name):
    """The repository's config, the per-worktree config and info/ (where
    info/attributes lives) are what git reads to decide which programs it
    runs. A command must not be able to change them for the git commands
    that run outside the sandbox later -- DELFIN's own and the user's."""
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    argv = A._bash_isolation_argv("git status", str(ws), perms)
    for gitdir in (repo / ".git", repo / ".git" / "worktrees" / "s1"):
        target = str((gitdir / name).resolve())
        assert target in argv, (name, gitdir)
        i = argv.index(target)
        assert argv[i - 1] == "--ro-bind-try"
        assert argv.index("--bind", 0, i) < i


@pytest.mark.skipif(not A._bwrap_functional(), reason="needs bubblewrap")
def test_a_sandboxed_command_leaves_the_git_config_as_it_was(worktree, monkeypatch):
    """Driven, not read from the argv: a git command and a script that both
    try to change the repository config leave it byte for byte as it was,
    and a commit in the same command still lands."""
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    repo, ws = worktree
    config = repo / ".git" / "config"
    before = config.read_bytes()
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    cmd = ("git config example.key value; "
           "python3 -c \"import sys; open(sys.argv[1], 'a').write('[example]\\\\n')\" "
           f"{config}; "
           "printf hi > g.txt && git add g.txt && "
           "git -c user.email=t@t -c user.name=t commit -q -m more && echo COMMIT-OK")
    argv = A._bash_isolation_argv(cmd, str(ws), perms)
    out = subprocess.run(argv, cwd=str(ws), capture_output=True, text=True, timeout=90)
    assert config.read_bytes() == before, out.stderr
    assert "COMMIT-OK" in out.stdout, (out.stdout, out.stderr)


def test_seatbelt_denies_writing_the_git_config(worktree, monkeypatch):
    """macOS: the profile denies writes below hooks/, config,
    config.worktree and info/ of every git dir -- read here from the
    profile text, since Seatbelt does not run on this host."""
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    from delfin.agent import seatbelt as SB
    seen = {}
    monkeypatch.setattr(SB, "argv", lambda plain, prof: seen.setdefault("p", prof) and plain)
    repo, ws = worktree
    perms = A.KitToolPermissions(workspace=str(ws), lock_workspace=True)
    A._seatbelt_argv(["true"], perms)
    deny = [ln for ln in seen["p"].splitlines() if ln.startswith("(deny file-write* (subpath")]
    assert deny, seen["p"]
    for gitdir in (repo / ".git", repo / ".git" / "worktrees" / "s1"):
        for name in ("hooks", "config", "config.worktree", "info"):
            assert str((gitdir / name).resolve()) in deny[-1], (gitdir, name)
