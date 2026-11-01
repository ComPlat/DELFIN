"""The result_checks tool stays out of a non-chemistry prompt.

The tool rides on the delfin-ops MCP server (``mcp__delfin-ops__*``),
which the workspace toolsets map to the ``chemistry`` group.  A plain
(non-DELFIN) workspace never offers that group, so ``drop_hidden``
removes the tool from the advertised schema — the built prompt for a
non-chemistry task is identical to what a tree without this tool
offers.  These tests pin that, plus the positive case for a DELFIN
workspace.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.toolsets import (
    CHEMISTRY,
    WorkspaceToolsets,
    drop_hidden,
    tool_group,
)


def _plain_repo(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "package.json").write_text(
        '{"name": "app", "scripts": {"test": "vitest"}}')
    (base / "src").mkdir()
    (base / "src" / "index.ts").write_text("export const x = 1;\n")
    return base


def _delfin_workspace(base: Path) -> Path:
    (base / "delfin" / "agent").mkdir(parents=True)
    (base / "delfin" / "__init__.py").write_text("")
    (base / "delfin" / "agent" / "__init__.py").write_text("")
    return base


def _tool_entry(name: str) -> dict:
    return {"type": "function",
            "function": {"name": name, "description": "d" * 3000}}


def test_result_checks_is_a_chemistry_group_tool():
    assert tool_group("mcp__delfin-ops__result_checks") == CHEMISTRY


def test_a_plain_workspace_hides_the_new_tool(tmp_path):
    ts = WorkspaceToolsets(_plain_repo(tmp_path / "app"))
    assert ts.offered(CHEMISTRY) is False
    hidden = ts.hidden_groups()
    assert CHEMISTRY in hidden
    kept = drop_hidden(
        [_tool_entry("mcp__delfin-ops__result_checks"),
         _tool_entry("mcp__delfin-ops__parse_orca_output"),
         _tool_entry("core_read_file")],
        hidden)
    names = [t["function"]["name"] for t in kept]
    assert "mcp__delfin-ops__result_checks" not in names
    assert "mcp__delfin-ops__parse_orca_output" not in names
    assert "core_read_file" in names


def test_a_delfin_workspace_keeps_the_new_tool(tmp_path):
    ts = WorkspaceToolsets(_delfin_workspace(tmp_path / "d"),
                           is_delfin_workspace=True)
    assert ts.offered(CHEMISTRY) is True
    kept = drop_hidden([_tool_entry("mcp__delfin-ops__result_checks")],
                       ts.hidden_groups())
    names = [t["function"]["name"] for t in kept]
    assert "mcp__delfin-ops__result_checks" in names
