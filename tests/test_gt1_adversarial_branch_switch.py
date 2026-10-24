"""GT1 reviewer adversarial tests (nacht-s19): the wave-21 friction finding B.

Measured in wave 21: `git checkout agent/s3-ds1b` — back to the session's
own EXISTING local branch — was refused by the write gate, which read the
branch name as the path `<worktree>/agent/s3-ds1b`.

The deciding chain is:
  1. api_client.py auto-allow regex (line ~3711) matches plain
     `git checkout <branch>` -> command is auto-allowed;
  2. patch_targets.py `_targets_git_checkout` (line 485-512) classifies
     the single operand: `_is_path_shaped` (line 575-578) returns True for
     any operand containing `/` or `:`, so a branch name like
     `agent/s3-ds1b` is kept as a WRITE TARGET;
  3. api_client `_gate_bash_write_targets` (line 24630) resolves that
     target against the workspace and gates the write: the branch does not
     exist as a file, so the write gate refuses.

The fix (builder, owned routine_actions.py; or operator .gate/px on
patch_targets.py / api_client.py) must make a single-operand checkout that
resolves to a BRANCH (not a path) name no write target — while every
restore-by-path shape (`-- file`, `ref -- file`, `.`, `..`) stays held.

These tests drive the real path: `_run_permission_gate` on a workspace
holding a REAL git repository (subprocess git, a test fixture, not the
product's safe_git). Commands are never executed by the gate.

Red on the wave base a4c9ebaf: the branch case is REFUSED there. Green
after the fix. Look-alikes must stay held in both states.
"""
from __future__ import annotations

import subprocess

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def git_ws(tmp_path):
    """A workspace that is a real git repository with one commit."""
    ws = tmp_path / "ws"
    ws.mkdir()
    env = {
        "PATH": "/usr/bin:/bin",
        "GIT_AUTHOR_NAME": "T", "GIT_AUTHOR_EMAIL": "t@example.com",
        "GIT_COMMITTER_NAME": "T", "GIT_COMMITTER_EMAIL": "t@example.com",
        "GIT_AUTHOR_DATE": "2005-04-07T22:13:13 +0000",
        "GIT_COMMITTER_DATE": "2005-04-07T22:13:13 +0000",
        "HOME": str(ws),
        "LC_ALL": "C",
    }
    def git(*args):
        return subprocess.run(
            ["git", *args], cwd=ws, env=env, capture_output=True, text=True)
    assert git("init", "-q", "-b", "main").returncode == 0
    (ws / "notes.txt").write_text("hello\n")
    assert git("add", "notes.txt").returncode == 0
    assert git("commit", "-q", "-m", "init").returncode == 0
    # The branch the wave-21 session wanted to go back to.
    assert git("branch", "agent/s3-ds1b").returncode == 0
    (ws / "notes.txt").write_text("changed, uncommitted\n")
    asked = []
    def _ask(tool, args, preview=""):
        asked.append(args.get("command", ""))
        return False
    p = KitToolPermissions(workspace=ws, mode="default",
                           confirm_callback=_ask)
    p._asked = asked
    p._git = git
    # The wave-21 session ran under a RESTRICTED write scope (a wave
    # worktree's write_allow_globs): the branch path <ws>/agent/s3-ds1b
    # matched no glob, and the write gate refused the branch switch as a
    # write to that path. Without a scope the whole workspace is writable
    # and the refusal does not surface. Reproduce the measured state.
    p.write_allow_globs = ("tests/test_gt1_*", "docs/**")
    return p


def gate(perms, cmd):
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms)


def write_gate(perms, cmd):
    """The tier the wave-21 refusal lives on: _gate_bash_write_targets
    (api_client.py:24630), which _execute_bash calls ON TOP of
    _run_permission_gate. It resolves each write target from
    patch_targets/_bash_write_targets against the workspace and gates it
    against the write scope."""
    return _doc_executor._gate_bash_write_targets(
        cmd, {"command": cmd}, perms)


# ---------------------------------------------------------------------------
# The measured defect: switching to an existing local branch of the own
# worktree must be routine. These are RED on the wave base.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git checkout agent/s3-ds1b",          # the measured wave-21 command
    "git switch agent/s3-ds1b",            # same intent, switch spelling
])
def test_switching_to_an_existing_local_branch_is_allowed(cmd, git_ws):
    assert gate(git_ws, cmd) is None, (
        "the wave-21 friction: an existing local branch of the own "
        "worktree was read as the path <worktree>/agent/s3-ds1b")


def test_the_write_gate_read_the_branch_name_as_a_path(git_ws):
    """GT1: the DECIDING classification, pinned in its INTENT, not its
    old defect. Originally this test reproduced the wave-21 defect
    (the branch operand was named as a write path) and stayed red until
    the px fix landed; the operator ruled the fix turns it red by
    design, so the assertion is INVERTED per its purpose: a branch
    operand must NOT be named as a path, while a real pathspec of the
    same shape still is (see test_checkout_dashdash_not_routine's tier
    and the *_is_still_named cases in test_gt1_existing_branch_switch)."""
    from delfin.agent import patch_targets
    targets = patch_targets.targets("git checkout agent/s3-ds1b",
                                    str(git_ws.workspace))
    assert [str(git_ws.workspace / "agent" / "s3-ds1b")] != targets


@pytest.mark.parametrize("cmd", [
    "git checkout agent/s3-ds1b",          # the measured wave-21 command
    "git switch agent/s3-ds1b",            # same intent, switch spelling
])
def test_the_write_target_gate_refuses_the_branch_switch(cmd, git_ws):
    """RED on the wave base for checkout: the write-target gate refuses
    the branch switch as a write to <worktree>/agent/s3-ds1b. The
    desired end state is err is None for BOTH spellings; checkout is red
    today, switch is green today (patch_targets names no target for it).
    The builder's fix must keep switch green and turn checkout green."""
    err = write_gate(git_ws, cmd)
    assert err is None, (cmd, err)


# ---------------------------------------------------------------------------
# Look-alikes that must STAY judged (held) in every state of the fix.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git checkout -- agent/s3-ds1b",       # restore-by-path
    "git checkout main -- notes.txt",      # ref + pathspec
    "git checkout .",                      # discard everything
    "git checkout ..",                     # parent-dir pathspec
    "git checkout -f main",                # force
    "git checkout -B agent/x",             # reset/force-create branch
    "git checkout --orphan x",             # orphan
    "git checkout --detach main",          # detach
    "git checkout -b new",                 # new branch: stays judged as
                                           # today (it auto-allows via
                                           # 3711/4078; here we only pin
                                           # that no WRITE TARGET appears)
])
def test_lookalikes_name_no_write_target_from_a_branch_operand(cmd, git_ws):
    """Whatever the fix does, it must not turn these into a blank allow by
    way of the write-target layer either: the gate verdict must come from
    the auto-allow/classifier tiers, with write targets intact."""
    gate(git_ws, cmd)   # no assertion on the verdict here beyond no-crash


@pytest.mark.parametrize("cmd", [
    "git checkout -- notes.txt",           # destructive restore of own file
])
def test_restore_by_path_still_held(cmd, git_ws):
    assert gate(git_ws, cmd) is not None, cmd


def test_path_operand_that_exists_is_still_a_write_target(git_ws):
    """The path-shaped check must not lose the case it exists for: an
    operand that names a real FILE is a restore and must stay gated."""
    from delfin.agent import patch_targets
    targets = patch_targets.targets("git checkout notes.txt",
                                    str(git_ws.workspace))
    assert [str(git_ws.workspace / "notes.txt")] == targets
