"""Agent safety CONTRACT — behaviour-level invariants enforced END-TO-END.

Unlike the per-fix unit tests, this pins the load-bearing agent-safety
properties THROUGH the real code paths (the executor's central authorization
gate; the prompt composition). A regression — a removed gate call, a prompt
edit that re-adds coding scaffolding to a guide role, a new tool that skips
scope — fails CI here. One canonical file so the safety surface is visible at
a glance. (Bugs 20260708-092217 & follow-ups.)
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
from delfin.agent.prompt_loader import PromptLoader


# ---------------------------------------------------------------------------
# 1. Dashboard guide cannot EXECUTE source/shell/edit tools (through execute())
# ---------------------------------------------------------------------------

def _dashboard_perms(tmp_path):
    return KitToolPermissions(
        workspace=tmp_path, mode="plan", agent_role="dashboard_agent")


def _is_scope_denied(result: str) -> bool:
    try:
        err = str(json.loads(result).get("error", ""))
    except Exception:
        return False
    return "not available to the 'dashboard_agent'" in err


@pytest.mark.parametrize("tool,args", [
    ("read_file", {"path": "delfin/agent/engine.py"}),
    ("grep_file", {"pattern": "token", "path": "delfin"}),
    ("list_files", {"path": "delfin"}),
    ("bash", {"command": "cat delfin/agent/engine.py"}),
    ("write_file", {"path": "x.py", "content": "x"}),
    ("edit_file", {"path": "x.py", "old_string": "a", "new_string": "b"}),
    ("apply_patch", {"patch": "..."}),
])
def test_dashboard_guide_denied_end_to_end(tmp_path, tool, args):
    ex = _DocToolExecutor()
    res = ex.execute(tool, args, _dashboard_perms(tmp_path))
    assert _is_scope_denied(res), f"{tool} not denied through execute(): {res[:200]}"


def test_dashboard_guide_research_tool_not_scope_denied(tmp_path):
    # search_docs is on the guide's allow-list — it must NOT hit the scope gate
    # (it may return an empty/no-index message, which is fine).
    ex = _DocToolExecutor()
    res = ex.execute("search_docs", {"query": "PBE0"}, _dashboard_perms(tmp_path))
    assert not _is_scope_denied(res)


def test_coding_role_not_scope_denied_end_to_end(tmp_path):
    # A coding role has no execution allow-list at this layer — the central gate
    # must NOT deny it (its plan-mode / per-tool gates apply elsewhere).
    ex = _DocToolExecutor()
    perms = KitToolPermissions(
        workspace=tmp_path, mode="default", agent_role="solo_agent")
    res = ex.execute("grep_file", {"pattern": "x", "path": str(tmp_path)}, perms)
    assert not _is_scope_denied(res)


# ---------------------------------------------------------------------------
# 1b. The SAME scope deny reaches MCP-dispatched tools. MCP tools skip
# execute() (they route straight to the registry), so _gate_mcp_tool must
# re-apply the role allow-list — else the guide reaches source/shell/edit
# through an MCP backend (bug 20260708-092217, MCP dispatch path).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,args", [
    ("mcp__kit-coding__bash", {"command": "cat delfin/agent/engine.py"}),
    ("mcp__delfin-docs__read_file", {"path": "delfin/agent/engine.py"}),
    ("mcp__delfin-docs__grep_file", {"pattern": "token"}),
    ("mcp__delfin-docs__list_files", {"pattern": "**/*.py"}),
    ("mcp__kit-coding__write_file", {"path": "x.py", "content": "x"}),
    ("mcp__kit-coding__edit_file",
     {"path": "x.py", "old_string": "a", "new_string": "b"}),
])
def test_dashboard_guide_denied_for_mcp_tools(tmp_path, name, args):
    ex = _DocToolExecutor()
    err = ex._gate_mcp_tool(name, args, _dashboard_perms(tmp_path))
    assert err is not None and "dashboard_agent" in err, \
        f"{name} not role-denied via MCP dispatch: {err!r}"


def test_dashboard_guide_mcp_research_tool_not_scope_denied(tmp_path):
    # search_docs is on the guide's allow-list → not role-denied even via MCP.
    ex = _DocToolExecutor()
    err = ex._gate_mcp_tool(
        "mcp__delfin-docs__search_docs", {"query": "PBE0"},
        _dashboard_perms(tmp_path))
    assert err is None


# ---------------------------------------------------------------------------
# 2. Prompt hygiene: guide role vs coding role
# ---------------------------------------------------------------------------

def _prompt(role, route, mode):
    return PromptLoader().build_system_prompt(
        role_id=role, mode_id=mode, route=route, task_text="open the submit tab")


_CODING_SCAFFOLD = (
    "Collaboration protocol",
    "automated multi-agent pipeline",
    "Efficiency rules",
    "read DELFIN source code, calculation data, archives (anywhere)",
    "Self-reflection",
)


def test_guide_prompt_free_of_coding_scaffolding():
    p = _prompt("dashboard_agent", ["dashboard_agent"], "dashboard")
    for frag in _CODING_SCAFFOLD:
        assert frag not in p, f"guide prompt regressed — carries {frag!r}"


def test_coding_prompt_keeps_scaffolding():
    p = _prompt("builder_agent", ["session_manager", "builder_agent"], "solo")
    assert "Collaboration protocol" in p
    assert "Efficiency rules" in p
