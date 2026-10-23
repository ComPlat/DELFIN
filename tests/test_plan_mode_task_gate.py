"""Plan-mode must not let the agent self-start a task list (bug 20260708-092217).

A dashboard agent running with perms=plan (READ-ONLY) created a task and moved
it to in_progress, and the per-turn open-tasks reminder then auto-drove it into
"↻ auto-continue → next task" with no user acceptance. Root cause: task tools
bypassed the plan-mode gate, and the reminder fired even in plan mode. These
tests pin the fix:

  * task_create / task_update(->in_progress|completed) are refused in plan mode,
  * metadata-only task_update and read-only task_list stay allowed,
  * the per-turn open-tasks reminder is suppressed in plan mode,
  * the dashboard agent's source-inspection tools are hard-denied.
"""

from __future__ import annotations

import json
import tempfile
import types
from pathlib import Path

from delfin.agent.api_client import _doc_executor, KitToolPermissions


def _perms(mode: str, ws: Path) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode=mode)
    p.task_session_id = "sess"
    return p


# --- executor gate ----------------------------------------------------------

def test_task_create_blocked_in_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        perms = _perms("plan", Path(d))
        out = json.loads(_doc_executor._execute_task_create(
            {"subject": "do the thing"}, perms))
        assert "error" in out
        assert "plan mode" in out["error"].lower()
        # Nothing was actually created.
        listed = json.loads(_doc_executor._execute_task_list({}, perms))
        assert listed["count"] == 0


def test_task_update_to_in_progress_blocked_in_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        perms = _perms("default", ws)
        created = json.loads(_doc_executor._execute_task_create(
            {"subject": "seeded"}, perms))
        tid = created["task"]["id"]
        # Now flip to plan mode and try to START the task.
        perms.mode = "plan"
        for status in ("in_progress", "completed"):
            out = json.loads(_doc_executor._execute_task_update(
                {"task_id": tid, "status": status}, perms))
            assert "error" in out, status
            assert "plan mode" in out["error"].lower(), status


def test_task_update_metadata_only_allowed_in_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        perms = _perms("default", ws)
        created = json.loads(_doc_executor._execute_task_create(
            {"subject": "seeded"}, perms))
        tid = created["task"]["id"]
        perms.mode = "plan"
        # Editing the subject (no status escalation) is harmless → allowed.
        out = json.loads(_doc_executor._execute_task_update(
            {"task_id": tid, "subject": "renamed"}, perms))
        assert out.get("status") == "updated"
        assert out["task"]["subject"] == "renamed"


def test_task_create_allowed_outside_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        perms = _perms("acceptEdits", Path(d))
        out = json.loads(_doc_executor._execute_task_create(
            {"subject": "real work"}, perms))
        assert out.get("status") == "created"


# --- per-turn reminder ------------------------------------------------------

def _open_tasks_block(perms) -> str:
    from delfin.agent.engine import AgentEngine
    stub = types.SimpleNamespace(kit_permissions=perms)
    return AgentEngine._build_open_tasks_block(stub)


def test_open_tasks_reminder_suppressed_in_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        seed = _perms("default", ws)
        _doc_executor._execute_task_create({"subject": "one"}, seed)
        # Same workspace/session, but now in plan mode → no reminder pressure.
        plan = _perms("plan", ws)
        assert _open_tasks_block(plan) == ""


def test_open_tasks_reminder_present_outside_plan_mode():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        perms = _perms("default", ws)
        _doc_executor._execute_task_create({"subject": "one"}, perms)
        block = _open_tasks_block(perms)
        assert "Open tasks" in block
        assert "one" in block


# --- dashboard source-inspection scope --------------------------------------

def test_dashboard_source_deny_covers_source_tools_not_docsearch():
    """Asked of the layer that refuses, not of a list that describes it.

    This used to read a second copy of the ban in engine.py, checked
    against the tool_use event — which is a report of a call already
    made. The question that matters is whether the executor refuses it,
    and it takes the namespaced name because that is how it arrives.
    """
    from delfin.agent.api_client import _tool_denied_for_role

    for tool in ("read_file", "grep_file", "list_files", "glob_files"):
        assert _tool_denied_for_role("dashboard_agent",
                                     f"mcp__delfin-docs__{tool}"), tool
        assert _tool_denied_for_role("dashboard_agent", tool), tool
    # … while the dashboard's legitimate research tools go through.
    for tool in ("search_docs", "read_section"):
        assert not _tool_denied_for_role("dashboard_agent",
                                         f"mcp__delfin-docs__{tool}"), tool
