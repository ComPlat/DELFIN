"""Tests for exit_plan_mode tool dispatch + permissions integration."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _make_perms(mode: str = "plan", **kwargs) -> KitToolPermissions:
    tmp = tempfile.mkdtemp(prefix="planmode_")
    return KitToolPermissions(workspace=Path(tmp), mode=mode, **kwargs)


def _exec_plan(plan: str, perms: KitToolPermissions) -> dict:
    out = _DocToolExecutor().execute(
        "exit_plan_mode", {"plan": plan}, permissions=perms,
    )
    return json.loads(out)


def test_exit_plan_no_callback_does_not_self_approve():
    """No approval channel → plan mode is a HARD human gate: the agent must
    NOT self-approve. Mode stays 'plan' so edits/writes remain blocked; the
    call returns 'awaiting_approval' and the agent should stop and wait.
    (Regression guard for the self-exit bug fixed 2026-07-14.)"""
    perms = _make_perms(mode="plan")
    result = _exec_plan("1. step one\n2. step two", perms)
    assert result["status"] == "awaiting_approval"
    assert perms.mode == "plan"                    # NOT flipped to 'default'
    assert not perms.last_approved_plan            # nothing was approved


def test_exit_plan_blocked_outside_plan_mode():
    perms = _make_perms(mode="default")
    result = _exec_plan("anything", perms)
    assert "error" in result
    assert "plan" in result["error"].lower()


def test_exit_plan_empty_plan_rejected():
    perms = _make_perms(mode="plan")
    result = _exec_plan("   ", perms)
    assert "error" in result
    assert perms.mode == "plan"


def test_exit_plan_with_callback_approve():
    seen = []

    def cb(plan):
        seen.append(plan)
        return {"approved": True, "new_mode": "acceptEdits"}

    perms = _make_perms(mode="plan", plan_approval_callback=cb)
    result = _exec_plan("a plan", perms)
    assert result["status"] == "approved"
    assert result["new_mode"] == "acceptEdits"
    assert perms.mode == "acceptEdits"
    assert seen == ["a plan"]


def test_approval_result_anchors_on_the_plan():
    """On approval the result re-states the plan as the authoritative, most
    recent instruction — so execution anchors on it, not on stale context.
    Hardening for bug 2026-06-25: a session whose context still held a prior
    Tetris build presented a correct spreadsheet plan, got it approved, then
    built Tetris from the leftover context."""
    plan = "# Plan: Mini-Spreadsheet\n1. grid.py\n2. formula.py\n3. tests"
    result = _exec_plan(plan, _make_perms(
        mode="plan", plan_approval_callback=lambda _p: {"approved": True}))
    assert result["status"] == "approved"
    instr = result.get("instruction", "")
    assert "grid.py" in instr and "formula.py" in instr   # full plan echoed
    assert "EXACTLY this approved plan" in instr
    assert "IGNORE it" in instr                            # ignore earlier task


def test_exit_plan_with_callback_reject():
    perms = _make_perms(
        mode="plan",
        plan_approval_callback=lambda _p: {"approved": False},
    )
    result = _exec_plan("rejected plan", perms)
    assert result["status"] == "rejected"
    assert perms.mode == "plan"


def test_exit_plan_callback_invalid_mode():
    perms = _make_perms(
        mode="plan",
        plan_approval_callback=lambda _p: {
            "approved": True, "new_mode": "garbage",
        },
    )
    result = _exec_plan("plan", perms)
    assert "error" in result
    assert perms.mode == "plan"


def test_exit_plan_callback_exception():
    def boom(_):
        raise RuntimeError("nope")
    perms = _make_perms(mode="plan", plan_approval_callback=boom)
    result = _exec_plan("p", perms)
    assert "error" in result
    assert "nope" in result["error"]


def test_plan_mode_blocks_writes_via_gate():
    """Sanity check that plan mode still blocks edit/write."""
    perms = _make_perms(mode="plan")
    out = _DocToolExecutor().execute(
        "write_file", {"path": "x.txt", "content": "hi"},
        permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload
    assert "plan mode" in payload["error"]


def test_plan_mode_lifecycle_unblocks_writes():
    """After exit_plan_mode → default, write_file should be possible
    (subject to confirm_callback / sandbox)."""
    confirms: list = []

    def confirm(name, args, preview):
        confirms.append(name)
        return True

    perms = _make_perms(mode="plan", confirm_callback=confirm,
                        plan_approval_callback=lambda _p: {"approved": True})
    _exec_plan("step 1", perms)
    assert perms.mode == "default"
    out = _DocToolExecutor().execute(
        "write_file",
        {"path": "x.txt", "content": "hello"},
        permissions=perms,
    )
    # Write succeeded (returns plain text on success, JSON error on
    # failure). The key invariant: plan-mode block is gone.
    assert "plan mode" not in out


def test_exit_plan_timeout_is_awaiting_not_rejected():
    """A plan-approval TIMEOUT must report 'awaiting_approval', not 'rejected',
    so the agent stops instead of resubmitting the plan (which re-blocks for the
    whole window). Regression for the 2026-06-25 21-min double-hang."""
    perms = _make_perms(
        mode="plan",
        plan_approval_callback=lambda _p: {"approved": False, "timed_out": True},
    )
    result = _exec_plan("a plan", perms)
    assert result["status"] == "awaiting_approval"
    assert "resubmit" in result["message"].lower()   # tells the agent to stop
    assert perms.mode == "plan"                        # still read-only (no flip)
