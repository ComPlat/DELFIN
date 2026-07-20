"""apply_patch is gated with the SAME per-path write policy as write_file.

Historically apply_patch only ran the self-mod guard, so a unified diff could
create ``.git/hooks/*`` (code execution on the next git op), ``.env``, or
overwrite a stored calc / read-only archive file — everything ``write_file``
hard-denies. The gate now runs deny-list, read-only-archive, self-mod-guard and
calc-confirm for every file the diff touches, natively AND over MCP.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / ".git" / "hooks").mkdir(parents=True)
    return ws


def _diff_create(path: str, content: str = "x") -> str:
    return f"--- /dev/null\n+++ b/{path}\n@@ -0,0 +1 @@\n+{content}\n"


def _apply_gate(perms, diff, **extra):
    return _doc_executor._run_permission_gate(
        "apply_patch", {"diff": diff, **extra}, perms)


def _write_gate(perms, path):
    return _doc_executor._run_permission_gate(
        "write_file", {"path": path, "content": "x"}, perms)


# ---------------------------------------------------------------------------
# Native gate: deny-listed / read-only / normal paths behave like write_file.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", [
    ".git/hooks/pre-commit",   # arbitrary code exec on next git op
    ".env",                     # secrets
])
def test_apply_patch_blocks_deny_listed_paths(workspace, path):
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    err = _apply_gate(perms, _diff_create(path))
    assert err is not None, f"apply_patch to {path} was not blocked"
    assert "deny-list" in err.lower() or "refus" in err.lower()


def test_apply_patch_allows_normal_workspace_path(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    assert _apply_gate(perms, _diff_create("src/module.py")) is None


@pytest.mark.parametrize("path", [
    ".git/config", ".env", "src/module.py", "a/b/c.txt",
])
def test_apply_patch_gate_matches_write_file_gate(workspace, path):
    """For every path, apply_patch blocks iff write_file blocks — the tiers
    are now identical (equivalence is the whole point of the fix)."""
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    apply_blocked = _apply_gate(perms, _diff_create(path)) is not None
    write_blocked = _write_gate(perms, path) is not None
    assert apply_blocked == write_blocked, (
        f"{path}: apply={apply_blocked} write={write_blocked}")


def test_apply_patch_multi_file_diff_blocks_on_any_bad_file(workspace):
    """A diff that touches one safe and one deny-listed file is refused."""
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    diff = _diff_create("ok.py") + _diff_create(".git/hooks/pre-push")
    assert _apply_gate(perms, diff) is not None


def test_apply_patch_self_mod_guard_headless_refused(workspace):
    """A diff editing the agent's own safety layer is refused head-less."""
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    diff = _diff_create("delfin/agent/engine.py")
    err = _apply_gate(perms, diff)
    assert err is not None and ("safety layer" in err or "self" in err.lower())


def test_apply_patch_check_only_not_gated(workspace):
    """check_only (git apply --check) is read-only → the write gate is not
    applied (so it stays usable, incl. in plan mode via the executor)."""
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    # Even a deny-listed target is not write-gated when only checking.
    assert _apply_gate(perms, _diff_create(".env"), check_only=True) is None


def test_apply_patch_blocked_in_plan_mode(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="plan")
    err = _apply_gate(perms, _diff_create("src/x.py"))
    assert err is not None and "plan mode" in err.lower()


# ---------------------------------------------------------------------------
# End-to-end through the executor: the file is never written.
# ---------------------------------------------------------------------------

def test_execute_apply_patch_does_not_write_denied_file(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    out = _doc_executor.execute(
        "apply_patch", {"diff": _diff_create(".env", "SECRET=1")}, perms)
    payload = json.loads(out)
    assert "error" in payload
    assert not (workspace / ".env").exists()


# ---------------------------------------------------------------------------
# MCP path: mcp__*__apply_patch is gated identically (was fully ungated).
# ---------------------------------------------------------------------------

def test_mcp_apply_patch_gated_like_native(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    diff = _diff_create(".git/hooks/pre-commit")
    mcp = _doc_executor._gate_mcp_tool(
        "mcp__kit-coding__apply_patch", {"diff": diff}, perms)
    native = _apply_gate(perms, diff)
    assert mcp == native
    assert mcp is not None  # blocked


def test_mcp_apply_patch_allows_normal_path(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               confirm_callback=None)
    assert _doc_executor._gate_mcp_tool(
        "mcp__kit-coding__apply_patch",
        {"diff": _diff_create("src/x.py")}, perms) is None


def test_mcp_apply_patch_refused_without_permissions():
    err = _doc_executor._gate_mcp_tool(
        "mcp__kit-coding__apply_patch",
        {"diff": _diff_create(".env")}, None)
    assert err is not None and "permission" in err.lower()
