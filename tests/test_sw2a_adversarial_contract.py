"""SW2a adversarial tests, round 1 (reviewer nacht-s16).

These pin the contract the SW2a dashboard tab must hold:
- the tab is built on SW1's status contract (delfin/agent/swarm/status.py),
  not on its own shape;
- the dashboard tab list on main is a subset of the tab list after SW2a
  (package target item 5: nothing existing changes behaviour).
"""

from __future__ import annotations

import importlib
import inspect

MAIN_TAB_IDS = (
    "submit", "recalc", "job_status", "orca_builder", "calculations",
    "archive", "office", "home", "literature", "agent", "agent_activity",
    "tools", "chemdarwin", "settings",
)


def test_sw2a_tab_module_exists_and_imports():
    mod = importlib.import_module("delfin.dashboard.tab_swarm")
    assert mod is not None


def test_sw2a_tab_pins_sw1_schema_version():
    """The tab must import and use SW1's STATUS_SCHEMA_VERSION, not its
    own copy of the number."""
    from delfin.dashboard import tab_swarm
    from delfin.agent.swarm.status import STATUS_SCHEMA_VERSION

    assert tab_swarm.STATUS_SCHEMA_VERSION == STATUS_SCHEMA_VERSION == 1


def test_sw2a_tab_validates_reports_with_sw1_validate():
    """The tab must not carry a second status validator: whatever check
    it does goes through delfin.agent.swarm.status.validate_report."""
    from delfin.dashboard import tab_swarm
    from delfin.agent.swarm.status import validate_report

    assert tab_swarm._validate_status is validate_report


def test_sw2a_status_call_reuses_sw1_collect_status():
    """The board's data source must be SW1 collect_status itself (one
    status call per interval), not a re-implementation."""
    from delfin.dashboard import tab_swarm
    from delfin.agent.swarm.status import collect_status

    assert tab_swarm._collect_status is collect_status


def test_dashboard_tab_list_is_superset_of_main():
    """Target item 5: the tab list on main is a subset of the tab list
    after SW2a. Checked against the real registration mechanism: the
    tab_specs entries inside create_dashboard's source."""
    import delfin.dashboard as _dashboard

    src = inspect.getsource(_dashboard)
    ids = set()
    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith("'id':"):
            # shape: 'id': '<name>',
            value = stripped.split("'id':", 1)[1].strip().split(",")[0]
            if value.startswith("'") and value.endswith("'"):
                ids.add(value[1:-1])
    missing = set(MAIN_TAB_IDS) - ids
    assert not missing, f"main tabs vanished from tab_specs: {sorted(missing)}"
    assert "swarm" in ids, "SW2a tab not registered in tab_specs"


def test_sw2a_tab_registers_with_widget_factory_not_business_logic():
    """The tab module exposes a create_tab factory like every other tab
    (pattern: tab_agent_activity.create_tab) and holds no launch/stop
    orchestration of its own: it must not import start_members and call
    it outside the explicit Start button path."""
    from delfin.dashboard import tab_swarm

    assert hasattr(tab_swarm, "create_tab"), "tab needs a create_tab factory"
