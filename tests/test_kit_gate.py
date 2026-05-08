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


def test_bash_no_match_default_fails_clearly(workspace):
    """In default mode, an unmatched command must fail with a clear hint —
    no confirm dialog popped, no silent allow."""
    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        bash_auto_allow_patterns=(r"^\s*pytest\b",),
        confirm_callback=lambda *_: True,  # would have approved if asked
    )
    err = _gate(perms, "bash",
                {"command": "rsync -a foo bar", "description": "sync"})
    assert err is not None
    assert "auto-allow" in err.lower()
    assert "remember_permission" in err.lower() or "bypass" in err.lower()


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
    """Reading absolute paths is allowed (the agent should be able to scout
    other dirs), but secret-style files must be refused."""
    fake_ssh = tmp_path / ".ssh"
    fake_ssh.mkdir()
    (fake_ssh / "id_rsa").write_text("FAKE-KEY")

    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = _doc_executor._execute_read_file(
        {"path": str(fake_ssh / "id_rsa")}, perms=perms,
    )
    assert "deny-glob" in out or "secret" in out.lower() or ".ssh" in out


def test_read_file_outside_workspace_allowed(workspace, tmp_path):
    """Plain file outside workspace + not on deny list → readable."""
    other = tmp_path / "other.txt"
    other.write_text("hello world\n")

    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = _doc_executor._execute_read_file(
        {"path": str(other)}, perms=perms,
    )
    assert "hello world" in out
