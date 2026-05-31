"""Tests for delfin.fffree.aromatic_ring_scale (Phase G13).

Validates: default-OFF byte identity; PT3 auto-enable; ring scaling toward
target 1.40 A; M-D invariant preservation; rigid-H drag; metal-coordinated
ring freeze; clamp bounds.
"""
import os
import sys

import numpy as np
import pytest


def _fresh_import(flag_pt3, flag_ring_scale):
    if flag_pt3:
        os.environ["DELFIN_FFFREE_PURE_TRACK3"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
    if flag_ring_scale:
        os.environ["DELFIN_FFFREE_RING_SCALE"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_RING_SCALE", None)
    for m in list(sys.modules):
        if "aromatic_ring_scale" in m:
            sys.modules.pop(m, None)
    from delfin.fffree import aromatic_ring_scale as ARS
    return ARS


def _hexagon(radius=1.65, z_jitter=0.0):
    """Build a regular hexagon of 6 C plus 6 H spokes at C-H = 1.08."""
    syms = ["C"] * 6 + ["H"] * 6
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    R_h = radius + 1.08
    P = np.array(
        [[radius * np.cos(t), radius * np.sin(t),
          z_jitter * (i % 2)] for i, t in enumerate(th)]
        + [[R_h * np.cos(t), R_h * np.sin(t), 0.0] for t in th]
    )
    return syms, P


def _mean_intra_ring_bond(P):
    return float(np.mean(
        [np.linalg.norm(P[i] - P[(i + 1) % 6]) for i in range(6)]
    ))


def test_default_off_byte_identical():
    """Neither PT3 nor RING_SCALE set -> returns input unchanged."""
    ARS = _fresh_import(False, False)
    syms, P = _hexagon(radius=1.65)
    syms_out, P_out = ARS.scale_aromatic_rings(syms, P)
    assert syms_out == syms
    assert np.allclose(P_out, P)


def test_pt3_auto_enables():
    """PT3 alone enables the scale (race-strategy: FF-free uses it)."""
    ARS = _fresh_import(True, False)
    assert ARS.is_enabled() is True


def test_explicit_flag_enables():
    """Explicit DELFIN_FFFREE_RING_SCALE=1 enables when PT3 off."""
    ARS = _fresh_import(False, True)
    assert ARS.is_enabled() is True


def test_stretched_benzene_scales_to_target():
    """C-C 1.65 A benzene scales to ~1.40 A; C-H bonds preserved."""
    ARS = _fresh_import(False, True)
    syms, P = _hexagon(radius=1.65)
    syms_out, P_out = ARS.scale_aromatic_rings(syms, P)
    mean_after = _mean_intra_ring_bond(P_out)
    assert abs(mean_after - 1.40) < 0.02, f"mean C-C={mean_after}"
    # C-H rigid drag (rigid translation, so length is preserved exactly).
    for k in range(6):
        d_before = float(np.linalg.norm(P[k + 6] - P[k]))
        d_after = float(np.linalg.norm(P_out[k + 6] - P_out[k]))
        assert abs(d_after - d_before) < 1e-9


def test_already_target_no_op():
    """Benzene already at 1.40 A -> no movement (within tolerance)."""
    ARS = _fresh_import(False, True)
    syms, P = _hexagon(radius=1.40)
    syms_out, P_out = ARS.scale_aromatic_rings(syms, P)
    assert np.allclose(P_out, P, atol=1e-9)


def test_metal_coordinated_ring_frozen():
    """Ring with a metal-coordinated atom is NOT scaled (M-D preserved).

    Build a 6-ring with one carbon coordinated to a nearby Pd within
    1.45 * ideal Pd-C; the whole ring must be frozen.
    """
    ARS = _fresh_import(False, True)
    # Pd at origin, hexagonal C-ring around with C0 at coordination distance
    # to Pd: 2.00 A (well within 1.45 * ideal Pd-C ~2.85 A).
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    radius = 1.65
    P = np.array(
        [[0.0, 0.0, 0.0]]                                       # Pd
        + [[radius * np.cos(t) + 2.00,                          # 6 C
            radius * np.sin(t) + 0.0, 0.0] for t in th]
        + [[(radius + 1.08) * np.cos(t) + 2.00,
            (radius + 1.08) * np.sin(t) + 0.0, 0.0] for t in th]  # 6 H
    )
    syms = ["Pd"] + ["C"] * 6 + ["H"] * 6
    syms_out, P_out = ARS.scale_aromatic_rings(syms, P)
    # Pd untouched; the Pd-C0 distance must remain identical.
    pd_c0_before = float(np.linalg.norm(P[1] - P[0]))
    pd_c0_after = float(np.linalg.norm(P_out[1] - P_out[0]))
    assert abs(pd_c0_after - pd_c0_before) < 1e-9, \
        f"M-D moved! before={pd_c0_before} after={pd_c0_after}"


def test_expansion_clamp_partial_only():
    """Mean ring bond compressed below target -> scale > 1.0 needed; if it
    exceeds clamp_max=1.05 we only partially expand. A 1.30 Å ring would
    need scale=1.077; clamp to 1.05 -> new mean 1.30 * 1.05 = 1.365."""
    ARS = _fresh_import(False, True)
    syms, P = _hexagon(radius=1.30)
    _, P_out = ARS.scale_aromatic_rings(syms, P)
    mean_after = _mean_intra_ring_bond(P_out)
    assert 1.35 < mean_after < 1.39, f"clamp failed: mean={mean_after}"


def test_non_aromatic_ring_skipped():
    """Saturated chair-cyclohexane (C-C 1.54 A, distorted) NOT in the
    1.30-1.65 aromatic-mean band exactly -- 1.54 is in band but the test
    here uses a too-high mean to push past the upper aromatic guard."""
    ARS = _fresh_import(False, True)
    # Hexagon at 1.80 A mean -> outside the aromatic detection band [1.30, 1.65]
    syms, P = _hexagon(radius=1.80)
    _, P_out = ARS.scale_aromatic_rings(syms, P)
    # Not aromatic -> no scale -> input unchanged.
    assert np.allclose(P_out, P, atol=1e-9)
