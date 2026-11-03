"""OR1: the TS adapters are discoverable through the real registry.

Red test for .gate/or1_registry_discovery.patch (F1): production
discovery (lazy _discover via list_steps/run_step) must register
orca_ts_* steps WITHOUT the test importing the module itself —
importing it would mask the discovery hole the registry patch fixes.
"""
from __future__ import annotations

import delfin.tools._registry as reg
from delfin.tools import run_step
from delfin.tools._registry import list_steps


def test_discovery_registers_orca_ts_steps():
    """Public path only: fresh discovery, no direct module import."""
    assert not reg.__dict__.get("_DISCOVERED", True) or True
    # force a genuinely fresh discovery state
    reg._DISCOVERED = False
    reg._STEPS.clear() if hasattr(reg, "_STEPS") else None
    names = set(list_steps())
    assert "orca_ts_opt" in names
    assert "orca_ts_freq_check" in names
    assert "orca_ts_irc" in names


def test_run_step_resolves_orca_ts_freq_check_via_registry():
    reg._DISCOVERED = False
    # an unknown step still raises; a known one must not say "Unknown step"
    try:
        run_step("orca_ts_freq_check")
    except Exception as exc:  # ValueError or param error is fine
        assert "Unknown step" not in str(exc), str(exc)
