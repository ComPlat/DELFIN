"""DS1 P1 red test: descriptor adapters must declare physical-validity limits.

Wave-21 rule (physical realism is correctness): a tool step must declare what
it is valid for, not merely describe it in prose. The descriptor adapters
(morfeus steric, Multiwfn bond orders, cclib parsing) declared parameters only
— no limits anywhere. The limits are attached by the realism modules
(delfin/tools/adapters/{morfeus,multiwfn,cclib}_realism.py) which the adapter
directory-walk discovery executes; this test drives the real path
(``_registry.get`` -> discovery), never a .gate copy.

On the previous commit every case here FAILS: no adapter carries
``validity_limits``.
"""
from __future__ import annotations

import pytest

from delfin.tools._registry import get, list_steps


DESCRIPTOR_STEPS = [
    "morfeus_buried_volume",
    "morfeus_cone_angle",
    "morfeus_sterimol",
    "morfeus_full_report",
    "multiwfn_bond_orders",
    "cclib_parse",
    "cclib_energies",
    "cclib_vibrations",
]


@pytest.mark.parametrize("name", DESCRIPTOR_STEPS)
def test_descriptor_step_declares_validity_limits(name):
    """Every descriptor step declares machine-readable validity limits."""
    assert name in list_steps(), f"step {name} not registered"
    limits = getattr(get(name), "validity_limits", None)
    assert limits, (
        f"{name} declares no validity_limits — physical realism must be "
        "declared, not described"
    )
    assert set(limits) >= {"applies_to", "not_valid_for"}, (
        f"{name} limits lack the required keys 'applies_to' and "
        f"'not_valid_for': {sorted(limits)}"
    )
    assert limits["applies_to"], "applies_to must not be empty"
    assert limits["not_valid_for"], "not_valid_for must not be empty"


def test_limits_come_from_the_real_discovery_path():
    """The realism modules load through the registry walk, not a side import.

    Red control for the wave-20 OR1 defect class: an install module nothing
    imports changes nothing. Discovery must be the thing that executes the
    realism modules, so a fresh get() carries the limits.
    """
    # list_steps() triggers discovery; a fresh process reaches here first.
    assert get("morfeus_buried_volume").validity_limits is not None


def test_morfeus_buried_volume_limits_are_correct():
    """Spot-check content, not just presence: %Vbur is a metal-complex
    descriptor valid for a metal centre plus its ligand shell, not for
    isolated organic molecules; the sphere radius is bounded.
    """
    bv = get("morfeus_buried_volume")
    limits = bv.validity_limits
    assert "metal" in " ".join(limits["applies_to"])
    assert any(
        "organic" in s or "no metal" in s for s in limits["not_valid_for"]
    )
    assert limits["radius_range"] == (1.0, 10.0)


def test_multiwfn_bond_order_limits_are_correct():
    """Spot-check content: bond-order analysis needs a wavefunction file from
    a converged QM run and is meaningless for a geometry alone.
    """
    bo = get("multiwfn_bond_orders")
    limits = bo.validity_limits
    assert any("wavefunction" in s for s in limits["applies_to"])
    assert limits["requires_converged"] is True


def test_cclib_parse_limits_are_correct():
    """cclib parses QM program output; it cannot parse a bare geometry."""
    parse = get("cclib_parse")
    limits = parse.validity_limits
    assert any("output" in s or "log" in s for s in limits["applies_to"])
