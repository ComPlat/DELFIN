"""Tests for delfin.fffree.oxoanion_vsepr_template -- iter-32a-3 re-activation.

Validate the build-time oxoanion VSEPR template:
  - Default OFF byte-identical (no env flag set).
  - Graph-based oxoanion detector (NO3-/ClO4-/SO4^2-/PO4^3-).
  - Td template projects perchlorate / sulfate / phosphate to 109.47 deg.
  - D3h template projects nitrate to 120 deg.
  - Rigid central atom: X never moved.
  - Donor oxygen invariant: when one O is a donor, its position is
    preserved (M-O distance invariant); when 2+ O are donors the
    motif is skipped.
  - Symmetric env flags (per-fix + master).
  - Determinism: same input -> identical output bytes (no random ordering).
"""
import os

import numpy as np
import pytest


_FLAGS = (
    "DELFIN_FFFREE_OXOANION_VSEPR",
    "DELFIN_FFFREE_CONSTRUCTION_FIX_ALL",
)


def _reset_env():
    for k in _FLAGS:
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


# ---- helpers ----

def _o_o_angle(P_c, P_a, P_b):
    """O-X-O angle in degrees at centre P_c."""
    va = P_a - P_c
    vb = P_b - P_c
    va /= np.linalg.norm(va)
    vb /= np.linalg.norm(vb)
    return float(np.degrees(np.arccos(np.clip(np.dot(va, vb), -1.0, 1.0))))


def _max_dev(P_c, oxys, ideal):
    devs = []
    n = len(oxys)
    for i in range(n):
        for j in range(i + 1, n):
            devs.append(abs(_o_o_angle(P_c, oxys[i], oxys[j]) - ideal))
    return max(devs) if devs else 0.0


# ---- 1. default-OFF byte identity ----


def test_default_off_byte_identical_nitrate():
    """No env flag: nitrate left untouched, n_fixed == 0."""
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    # severely distorted nitrate: O-N-O = 90 deg pair (broken)
    syms = ["N", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],         # N
        [1.245, 0.0, 0.0],       # O1
        [0.0, 1.245, 0.0],       # O2 (90 deg from O1)
        [-1.245, 0.0, 0.0],      # O3 (180 deg from O1)
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


def test_default_off_byte_identical_perchlorate():
    """No env flag: perchlorate left untouched."""
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["Cl", "O", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.45, 0.0, 0.0],
        [0.0, 1.45, 0.0],
        [0.0, 0.0, 1.45],
        [-1.45, -1.45, -1.45],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


# ---- 2. nitrate (D3h) ----


def test_nitrate_projected_to_120_deg():
    """Distorted nitrate -> O-N-O = 120 deg."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    # Build a clearly-broken nitrate: O-N-O pair at 90 deg
    syms = ["N", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.245, 0.0, 0.0],
        [0.0, 1.245, 0.0],
        [-1.245, -1.245, 0.0],
    ], dtype=float)
    dev_before = _max_dev(P[0], P[1:], 120.0)
    assert dev_before > 5.0
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1
    dev_after = _max_dev(P_new[0], P_new[1:], 120.0)
    assert dev_after < 1.0, f"nitrate angle dev still {dev_after:.2f} deg"


def test_nitrate_central_n_never_moves():
    """Central N atom must remain at its input coordinate (rigid pivot)."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["N", "O", "O", "O"]
    P = np.array([
        [1.5, 2.5, -0.7],
        [2.745, 2.5, -0.7],
        [1.5, 3.745, -0.7],
        [0.255, 1.255, -0.7],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1
    np.testing.assert_allclose(P_new[0], P[0], atol=1e-9)


# ---- 3. perchlorate / sulfate / phosphate (Td) ----


def test_perchlorate_projected_to_109_5_deg():
    """Distorted perchlorate -> O-Cl-O = 109.47 deg."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["Cl", "O", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.45, 0.0, 0.0],
        [0.0, 1.45, 0.0],
        [-1.45, 0.0, 0.0],
        [0.0, -1.45, 0.0],
    ], dtype=float)  # planar -- 90 deg pairs everywhere
    dev_before = _max_dev(P[0], P[1:], 109.471)
    assert dev_before > 5.0
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1
    dev_after = _max_dev(P_new[0], P_new[1:], 109.471)
    assert dev_after < 5.0, f"ClO4 dev still {dev_after:.2f} deg"


def test_sulfate_projected_to_109_5_deg():
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["S", "O", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.49, 0.0, 0.0],
        [0.0, 1.49, 0.0],
        [-1.49, 0.0, 0.0],
        [0.0, -1.49, 0.0],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1
    dev_after = _max_dev(P_new[0], P_new[1:], 109.471)
    assert dev_after < 5.0


def test_phosphate_projected_to_109_5_deg():
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["P", "O", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.54, 0.0, 0.0],
        [0.0, 1.54, 0.0],
        [-1.54, 0.0, 0.0],
        [0.0, -1.54, 0.0],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1
    dev_after = _max_dev(P_new[0], P_new[1:], 109.471)
    assert dev_after < 5.0


# ---- 4. M-D invariant ----


def test_donor_oxygen_position_preserved():
    """When one oxygen is a donor, its position is preserved exactly
    (M-O distance invariant by construction).

    Geometry: Cu at origin; nitrate at N = (2.7, 0, 0); donor O1 closer
    to the metal, O2/O3 distorted with WIDE inter-O separation so the
    motif detector sees a clean N-O-only neighbourhood (no
    spurious O-O bond perception).
    """
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["Cu", "O", "N", "O", "O"]
    # Distorted nitrate: O1 in the M-N axis, O2/O3 at ~90 deg from O1,
    # placed >2 A apart from each other so they don't form a spurious
    # O-O bond.  Bond N-O distances all ~1.245 A.
    P = np.array([
        [0.0, 0.0, 0.0],         # metal
        [1.5, 0.0, 0.0],         # donor O1  (M-O1 = 1.5)
        [2.745, 0.0, 0.0],       # central N
        [2.745, 1.245, 0.0],     # O2  (90 deg from O1 -- distorted)
        [2.745, 0.0, 1.245],     # O3  (90 deg from O1 -- distorted, z-axis)
    ], dtype=float)
    metal_idx = 0
    donor_idxs = (1,)
    md_before = float(np.linalg.norm(P[1] - P[0]))
    P_new, n = enforce_oxoanion_vsepr(syms, P, metal_idx=metal_idx,
                                       donor_idxs=donor_idxs)
    assert n == 1
    md_after = float(np.linalg.norm(P_new[1] - P_new[0]))
    assert abs(md_before - md_after) < 1e-6, (
        f"M-O changed: {md_before:.6f} -> {md_after:.6f}"
    )


def test_skipped_when_two_donor_oxygens():
    """When 2+ oxygens are donors (kappa^2), the template is skipped
    (geometry already constrained by chelate)."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["Cu", "O", "N", "O", "O"]
    # same well-separated-O nitrate as the donor test
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.745, 0.0, 0.0],
        [2.745, 1.245, 0.0],
        [2.745, 0.0, 1.245],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P, metal_idx=0,
                                       donor_idxs=(1, 3))
    assert n == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


# ---- 5. master env flag ----


def test_master_flag_enables():
    """Master flag DELFIN_FFFREE_CONSTRUCTION_FIX_ALL enables the fix too."""
    os.environ["DELFIN_FFFREE_CONSTRUCTION_FIX_ALL"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import (
        enforce_oxoanion_vsepr, _flag_active,
    )
    assert _flag_active()

    syms = ["N", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.245, 0.0, 0.0],
        [0.0, 1.245, 0.0],
        [-1.245, -1.245, 0.0],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 1


# ---- 6. detector ----


def test_find_oxoanions_detects_NO3_ClO4_SO4_PO4():
    """All four oxoanions detected at correct (x_idx, [o_idx], geom)."""
    from delfin.fffree.oxoanion_vsepr_template import find_oxoanions

    # NO3
    syms_n = ["N", "O", "O", "O"]
    P_n = np.array([
        [0.0, 0.0, 0.0], [1.245, 0.0, 0.0],
        [-0.62, 1.08, 0.0], [-0.62, -1.08, 0.0],
    ], dtype=float)
    out = find_oxoanions(syms_n, P_n)
    assert len(out) == 1 and out[0][0] == 0 and out[0][2] == "D3h"

    # ClO4
    syms_c = ["Cl", "O", "O", "O", "O"]
    a = 1.45 / np.sqrt(3.0)
    P_c = np.array([
        [0.0, 0.0, 0.0],
        [a, a, a], [a, -a, -a],
        [-a, a, -a], [-a, -a, a],
    ], dtype=float)
    out = find_oxoanions(syms_c, P_c)
    assert len(out) == 1 and out[0][0] == 0 and out[0][2] == "Td"

    # SO4
    syms_s = ["S", "O", "O", "O", "O"]
    a = 1.49 / np.sqrt(3.0)
    P_s = np.array([
        [0.0, 0.0, 0.0],
        [a, a, a], [a, -a, -a],
        [-a, a, -a], [-a, -a, a],
    ], dtype=float)
    out = find_oxoanions(syms_s, P_s)
    assert len(out) == 1 and out[0][2] == "Td"

    # PO4
    syms_p = ["P", "O", "O", "O", "O"]
    a = 1.54 / np.sqrt(3.0)
    P_p = np.array([
        [0.0, 0.0, 0.0],
        [a, a, a], [a, -a, -a],
        [-a, a, -a], [-a, -a, a],
    ], dtype=float)
    out = find_oxoanions(syms_p, P_p)
    assert len(out) == 1 and out[0][2] == "Td"


def test_carbonate_not_treated_as_oxoanion():
    """CO3 (carbonate) is NOT a nitrate / phosphate -- C is not in the
    central-atom template set, so it must be skipped."""
    from delfin.fffree.oxoanion_vsepr_template import find_oxoanions

    syms = ["C", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0], [1.28, 0.0, 0.0],
        [-0.64, 1.11, 0.0], [-0.64, -1.11, 0.0],
    ], dtype=float)
    assert find_oxoanions(syms, P) == []


def test_nitro_attached_to_C_not_treated():
    """An R-NO2 nitro group (N has C+O+O neighbours) must not match the
    NO3- template (which requires ALL 3 N-neighbours to be O)."""
    from delfin.fffree.oxoanion_vsepr_template import find_oxoanions

    syms = ["C", "N", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],      # C
        [1.50, 0.0, 0.0],     # N
        [2.25, 1.08, 0.0],    # O1
        [2.25, -1.08, 0.0],   # O2
    ], dtype=float)
    assert find_oxoanions(syms, P) == []


# ---- 7. determinism ----


def test_determinism_byte_identical_two_runs():
    """Same syms/P twice -> byte-identical output (no random ordering)."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    syms = ["Cl", "O", "O", "O", "O"]
    P = np.array([
        [0.5, -0.3, 0.7],
        [1.95, -0.3, 0.7],
        [0.5, 1.15, 0.7],
        [-0.95, -0.3, 0.7],
        [0.5, -1.75, 0.7],
    ], dtype=float)
    P_a, n_a = enforce_oxoanion_vsepr(syms, P.copy())
    P_b, n_b = enforce_oxoanion_vsepr(syms, P.copy())
    assert n_a == n_b == 1
    # Byte identity of the float arrays
    assert P_a.tobytes() == P_b.tobytes()


def test_accept_if_better_gate():
    """If the input is already perfect, n_fixed == 0 (no spurious change)."""
    os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"
    from delfin.fffree.oxoanion_vsepr_template import enforce_oxoanion_vsepr

    # Perfect Td perchlorate
    syms = ["Cl", "O", "O", "O", "O"]
    a = 1.45 / np.sqrt(3.0)
    P = np.array([
        [0.0, 0.0, 0.0],
        [a, a, a], [a, -a, -a],
        [-a, a, -a], [-a, -a, a],
    ], dtype=float)
    P_new, n = enforce_oxoanion_vsepr(syms, P)
    assert n == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)
