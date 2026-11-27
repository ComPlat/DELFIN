"""The tool-reach table names routes the product really has.

``delfin/tool_reach.py`` says, per CONTROL module, which registered building
blocks and which agent tools reach it — the count the measurement package
(wave 53) needs.  A count is only worth something if the routes in it resolve,
so this test checks both halves against the code:

* every ``steps`` name must be a building block the registry holds;
* every ``agent_tools`` name must be in the agent's tool catalogue.

A name that stopped existing fails here instead of inflating a reachability
figure. The module that is reachable *only* by a CONTROL run is asserted too —
as unreachable through the tool surface — because "every module is reachable"
and "we counted three different things under one word" are not the same
statement.

No ORCA, no xTB, no network: the registry and the catalogue are Python data.
"""

from __future__ import annotations

import pytest

from delfin import tool_reach
from delfin.tools._registry import list_steps


def _catalogue_names() -> set[str]:
    import delfin.api as api

    return {entry["name"] for entry in api._TOOL_CATALOG}


def test_every_building_block_the_table_names_exists():
    present = set(list_steps())
    missing = {
        route.module: [s for s in route.steps if s not in present]
        for route in tool_reach.TOOL_ROUTES
        if any(s not in present for s in route.steps)
    }
    assert not missing, (
        "the tool-reach table names building blocks the registry does not hold: "
        f"{missing}"
    )


def test_every_agent_tool_the_table_names_is_in_the_catalogue():
    catalogue = _catalogue_names()
    missing = {
        route.module: [t for t in route.agent_tools if t not in catalogue]
        for route in tool_reach.TOOL_ROUTES
        if any(t not in catalogue for t in route.agent_tools)
    }
    assert not missing, (
        "the tool-reach table names agent tools the catalogue does not offer: "
        f"{missing}"
    )


def test_the_report_resolves_every_route_it_claims():
    rows = tool_reach.module_tool_reach()
    assert rows, "the report is empty"
    broken = [
        r for r in rows
        if r["steps_missing"] or (r["agent_tools"] and not r["reachable"])
    ]
    assert not broken, f"routes that do not resolve: {broken}"


def test_thermodynamics_is_reported_as_reached_by_a_control_run_only():
    """The honest gap: no building block, no agent tool — a CONTROL run only.

    Asserted as *unreachable through the tool surface* so that a later change
    which adds a tool shows up as a deliberate edit of this test, not as a
    silently improved number.
    """
    rows = {r["module"]: r for r in tool_reach.module_tool_reach()}
    thdy = next(r for name, r in rows.items() if name.startswith("thermodynamics"))
    assert thdy["steps"] == []
    assert thdy["agent_tools"] == []
    assert thdy["reachable"] is False, (
        "thermodynamics became reachable through a tool; update the measurement "
        "and this assertion together"
    )


@pytest.mark.parametrize("module_prefix", ["ESD", "IMAG", "OCCUPIER", "CO2"])
def test_the_modules_the_package_names_have_a_route(module_prefix):
    """ESD, IMAG, OCCUPIER and CO2 each have at least one route that resolves."""
    rows = {r["module"]: r for r in tool_reach.module_tool_reach()}
    row = next((r for name, r in rows.items() if name.startswith(module_prefix)), None)
    assert row is not None, f"no tool-reach entry for {module_prefix}"
    assert row["reachable"], f"{module_prefix} has no reachable route: {row}"
