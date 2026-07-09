"""Central deny-by-default execution authorization (bug 20260708-092217).

The dashboard guide reached read_file / grep_file / list_files / bash /
task_create through the KIT backend because the engine's per-role tool
whitelist exempts every ``mcp__…``-namespaced tool. The fix re-checks a
per-role execution allow-list at the ONE choke point every tool passes
through (``_DocToolExecutor.execute`` → ``_tool_denied_for_role``), so the
leak is impossible regardless of the advertising / namespacing path.
"""

from __future__ import annotations

from delfin.agent.api_client import (
    _DASHBOARD_AGENT_ALLOWED_TOOLS,
    _tool_denied_for_role,
)


# --- dashboard_agent: source-reading / shell / edit tools are denied --------

def test_dashboard_role_denies_source_and_shell_tools():
    for tool in ("read_file", "grep_file", "list_files", "glob_files",
                 "bash", "bash_background", "write_file", "edit_file",
                 "multi_edit", "apply_patch", "run_tests"):
        assert _tool_denied_for_role("dashboard_agent", tool), tool


def test_dashboard_role_allows_its_research_and_ui_tools():
    for tool in _DASHBOARD_AGENT_ALLOWED_TOOLS:
        assert not _tool_denied_for_role("dashboard_agent", tool), tool


def test_dashboard_role_denies_mcp_namespaced_leak():
    # The exact vector from the bug: namespaced coding/doc tools.
    assert _tool_denied_for_role("dashboard_agent", "mcp__delfin-docs__read_file")
    assert _tool_denied_for_role("dashboard_agent", "mcp__delfin-docs__grep_file")
    assert _tool_denied_for_role("dashboard_agent", "mcp__kit-coding__bash")
    assert _tool_denied_for_role("dashboard_agent", "mcp__kit-coding__task_create") is False
    # …but a namespaced allowed tool still passes.
    assert not _tool_denied_for_role("dashboard_agent", "mcp__delfin-docs__search_docs")


def test_always_allowed_meta_tools_pass_for_dashboard():
    for tool in ("exit_plan_mode", "ask_user_question", "subagent_result"):
        assert not _tool_denied_for_role("dashboard_agent", tool)


# --- coding roles are unrestricted at THIS layer ----------------------------

def test_coding_roles_are_not_restricted_here():
    for role in ("solo_agent", "builder_agent", "test_agent", ""):
        for tool in ("read_file", "bash", "write_file", "edit_file"):
            assert not _tool_denied_for_role(role, tool), (role, tool)


def test_unknown_role_is_not_restricted():
    assert not _tool_denied_for_role("some_future_role", "bash")
