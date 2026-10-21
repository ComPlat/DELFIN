# Probe (NOT committed): dashboard <-> research mode switch.
#
# Research runs its engine in ctx.agent_dir (the fixed agent scratch
# folder, the same folder Dashboard uses) while its mode name is
# "research". _mode_workspace_differs groups by workspace: office /
# dashboard / launch. With research grouped under "launch", a switch
# Dashboard -> Research (or back) does NOT drop the engine, so the next
# send reuses or rebuilds an engine bound to the WRONG workspace and the
# previous mode's role.
#
# Red on the current branch: research groups as "launch", so the two
# switch directions report False and the tests fail.

from __future__ import annotations

import ast
import inspect
from pathlib import Path

from delfin.dashboard import tab_agent as T

_SRC = Path(inspect.getfile(T)).read_text(encoding="utf-8")


def test_switch_dashboard_research_drops_the_engine():
    assert T._mode_workspace_differs("dashboard", "research")
    assert T._mode_workspace_differs("research", "dashboard")


def test_switch_research_code_drops_the_engine():
    # research works in the fixed scratch folder, Code follows the launch
    # directory: different folders, the engine must be rebuilt.
    assert T._mode_workspace_differs("research", "solo")


def test_switch_within_launch_group_keeps_the_engine():
    # A switch that changes nothing must not discard a live session: the
    # launch group keeps its existing pairs.
    assert not T._mode_workspace_differs("solo", "pipeline")
    assert not T._mode_workspace_differs("pipeline", "solo")


def test_research_group_is_named_in_the_workspace_differs_docstring():
    # The grouping is the contract: the docstring lists which modes sit in
    # which group, so a new mode cannot silently fall into "launch" again.
    tree = ast.parse(_SRC)
    found = False
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "_mode_workspace_differs":
            doc = ast.get_docstring(node) or ""
            assert "research" in doc, (
                "_mode_workspace_differs docstring does not name research; "
                "the grouping is undocumented and will drift again"
            )
            found = True
    assert found, "_mode_workspace_differs not found"
