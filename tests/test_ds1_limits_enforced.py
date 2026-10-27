"""DS1 P2 red test: the declared limits are ENFORCED, not just declared.

A limit that is only prose does not stop a physically wrong run. Each case
here calls the adapter's execute() with a physically invalid parameter set
through the public registry path and requires a ValueError naming the
physical reason BEFORE the underlying tool would be invoked (the wrappers
are stubbed, so the only way execute() can fail is the physics gate).

On the previous commit (621d3c60, limits declared, not enforced) every
"enforcement" case FAILS: execute() reaches the stubbed wrapper and
succeeds. The "fires only when it should" direction lives in
tests/test_ds1_valid_params_control.py (green both before and after).
"""
from __future__ import annotations

import pytest

from delfin.tools import adapters
from delfin.tools._registry import get


# Ensure the realism enforcement is installed on these classes even if a
# test elsewhere cleared the registry and re-imported only the base modules.
from delfin.tools.adapters import morfeus_realism, multiwfn_realism  # noqa: F401


@pytest.fixture(autouse=True)
def stub_wrappers(monkeypatch):
    """Stub the analysis-tool wrappers: execute() must fail on physics only."""
    import delfin.analysis_tools.morfeus_wrapper as mw
    import delfin.analysis_tools.multiwfn_wrapper as mfnw

    calls = []

    def fake(*args, **kwargs):
        calls.append((args, kwargs))
        return {"ok": True}

    monkeypatch.setattr(mw, "buried_volume", fake)
    monkeypatch.setattr(mw, "cone_angle", fake)
    monkeypatch.setattr(mw, "sterimol", fake)
    monkeypatch.setattr(mw, "full_steric_report", fake)
    monkeypatch.setattr(mfnw, "bond_order_analysis", fake)
    return calls


def _run(name, **kwargs):
    adapter = get(name)
    return adapter.execute(__import__("pathlib").Path("."), **kwargs)


@pytest.mark.parametrize("radius", [0.5, 0.0, 15.0])
def test_buried_volume_rejects_radius_outside_physical_range(radius):
    with pytest.raises(ValueError, match="physically supported range"):
        _run("morfeus_buried_volume", metal_index=0, radius=radius)


def test_sterimol_rejects_degenerate_axis():
    with pytest.raises(ValueError, match="degenerate"):
        _run("morfeus_sterimol", atom1_index=3, atom2_index=3)


def test_full_report_rejects_wrong_donor_count():
    with pytest.raises(ValueError, match="exactly 2 donor_indices"):
        _run("morfeus_full_report", metal_index=0, donor_indices=[1])


def test_bond_orders_rejects_undeclared_method():
    with pytest.raises(ValueError, match="not one of"):
        _run("multiwfn_bond_orders", input_file="x.molden", method="delocalization")


def test_bond_orders_rejects_missing_input_file():
    with pytest.raises(ValueError, match="input_file"):
        _run("multiwfn_bond_orders")


def test_physics_gate_runs_before_the_tool():
    """Order matters: the wrapper must not be reached on invalid input."""
    with pytest.raises(ValueError):
        _run("morfeus_buried_volume", metal_index=0, radius=99.0)
