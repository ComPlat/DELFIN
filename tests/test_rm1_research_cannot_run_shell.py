# RM1 finding (s21): the research role must deny shell execution.
"""The dashboard promises research has "no shell commands" (tab_agent.py
description, README.md:523). The enforcement is the per-role denylist in
api_client.py; this control pins it at the real gate site,
_tool_denied_for_role — the function the executor calls before a tool
runs (api_client.py:7370). Red before the patch: the denylist for
research_agent is _WRITE_TOOL_NAMES, which does not contain bash.
"""

from delfin.agent.api_client import _tool_denied_for_role


def test_research_role_denies_bash():
    assert _tool_denied_for_role("research_agent", "bash")


def test_research_role_denies_bash_background():
    assert _tool_denied_for_role("research_agent", "bash_background")


def test_research_role_denies_namespaced_bash():
    # The gate strips the MCP namespace (api_client.py:7378), so a
    # namespaced call must be denied by the same entry.
    assert _tool_denied_for_role("research_agent", "mcp__kit-coding__bash")


def test_research_role_keeps_ordinary_reads():
    # Narrow the surface, do not disable the role.
    assert not _tool_denied_for_role("research_agent", "read_file")
    assert not _tool_denied_for_role("research_agent", "search_docs")


def test_write_tools_stay_denied_for_research():
    assert _tool_denied_for_role("research_agent", "write_file")
    assert _tool_denied_for_role("research_agent", "edit_file")


def test_other_roles_keep_their_shell():
    # The change must touch only the research role.
    assert not _tool_denied_for_role("solo_agent", "bash")
    assert not _tool_denied_for_role("builder_agent", "bash")
