"""DS1 P1 opposite-direction control: valid params still validate cleanly.

Red-control discipline, both directions: the phase-2 enforcement test must
FIRE on invalid params and NOT FIRE on valid ones. This file pins the
"does not fire" direction against the unchanged original behavior — an
ordinary, physically sensible parameter set must keep passing validation.
"""
from __future__ import annotations

import pytest

from delfin.tools._registry import get


def test_buried_volume_valid_params_pass():
    bv = get("morfeus_buried_volume")
    bv.validate_params(metal_index=0, radius=3.5)


def test_full_report_valid_params_pass():
    fr = get("morfeus_full_report")
    fr.validate_params(metal_index=0, donor_indices=[1, 2], radius=3.5)


def test_sterimol_valid_params_pass():
    st = get("morfeus_sterimol")
    st.validate_params(atom1_index=1, atom2_index=4)


def test_bond_orders_valid_params_pass():
    bo = get("multiwfn_bond_orders")
    bo.validate_params(input_file="x.molden")
    bo.validate_params(input_file="x.molden", method="wiberg")


@pytest.mark.parametrize("name", [
    "cclib_parse", "cclib_energies", "cclib_vibrations",
])
def test_cclib_valid_params_pass(name):
    get(name).validate_params(filepath="job.out")
