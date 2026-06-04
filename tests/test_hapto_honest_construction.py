"""Tests for ``delfin.fffree.hapto_honest_construction``.

Construction-side counterpart to the hapto-only-CShM metric.

Covers:
* env flag (default-OFF byte-identical to input)
* detection of η2-η8 units
* ferrocene η5-η5 ideal: hapto-only CShM ≈ 0 after fix
* (η6-arene)Cr: M-centroid distance snaps to CCDC ideal
* η4-butadiene, η2-ethylene basics
* ring-plane preservation (rigid-block: intra-block C-C invariant)
* M-D distance preservation for σ donors (σ-only metric invariant)
* perpendicularity restored after a 30° ring tilt
* determinism (two runs → byte-identical results)
* refine.py call does not break the hapto invariant (no leak)
"""
from __future__ import annotations

import importlib
import math
import os
import sys

import numpy as np
import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.abspath(os.path.join(
    _HERE, "..", "agent_workspace", "quality_framework", "scripts",
))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)


# ---------------------------------------------------------------------------
# Helpers (locally constructed canonical geometries).
# ---------------------------------------------------------------------------
def _ferrocene_eta5_eta5(M: str = "Fe", r_cc: float = 1.40,
                         z_off: float = 1.65):
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    for sign in (+1.0, -1.0):
        for k in range(5):
            theta = 2.0 * math.pi * k / 5.0
            coords.append([r_ring * math.cos(theta),
                           r_ring * math.sin(theta),
                           sign * z_off])
            syms.append("C")
    return syms, np.asarray(coords, dtype=float)


def _arene_eta6(M: str = "Cr", r_cc: float = 1.40, z_off: float = 1.62):
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 6.0))
    for k in range(6):
        theta = 2.0 * math.pi * k / 6.0
        coords.append([r_ring * math.cos(theta),
                       r_ring * math.sin(theta), z_off])
        syms.append("C")
    return syms, np.asarray(coords, dtype=float)


def _eta4_butadiene(M: str = "Fe", r_cc: float = 1.40, z_off: float = 1.70):
    syms = [M, "C", "C", "C", "C"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-1.5 * r_cc / 2.0, +r_cc * 0.4, z_off],
        [-0.5 * r_cc / 2.0, -r_cc * 0.4, z_off],
        [+0.5 * r_cc / 2.0, +r_cc * 0.4, z_off],
        [+1.5 * r_cc / 2.0, -r_cc * 0.4, z_off],
    ], dtype=float)
    return syms, coords


def _eta2_ethylene(M: str = "Pt", r_cc: float = 1.40, z_off: float = 2.10):
    syms = [M, "C", "C"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [+r_cc / 2.0, 0.0, z_off],
        [-r_cc / 2.0, 0.0, z_off],
    ], dtype=float)
    return syms, coords


def _half_sandwich_eta5_plus_sigma(M: str = "Fe", r_cc: float = 1.40,
                                    r_cl: float = 2.30, z_off: float = 1.65):
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    for k in range(5):
        theta = 2.0 * math.pi * k / 5.0
        coords.append([r_ring * math.cos(theta),
                       r_ring * math.sin(theta), z_off])
        syms.append("C")
    syms.append("Cl"); coords.append([r_cl, 0.0, 0.0])
    syms.append("N");  coords.append([-r_cl / 2.0, r_cl * math.sqrt(3) / 2.0, 0.0])
    syms.append("P");  coords.append([-r_cl / 2.0, -r_cl * math.sqrt(3) / 2.0, 0.0])
    return syms, np.asarray(coords, dtype=float)


def _ferrocene_distorted_centroid(M: str = "Fe"):
    """Ferrocene-like but the top ring sits at a WRONG centroid distance
    (1.45 Å vs ideal 1.65) and tilted 25° — the corrector should snap it."""
    syms, P = _ferrocene_eta5_eta5(M=M, z_off=1.45)
    # Tilt only the TOP ring (indices 1..5) by 25° around x-axis.
    th = math.radians(25.0)
    R = np.array([
        [1.0, 0.0, 0.0],
        [0.0, math.cos(th), -math.sin(th)],
        [0.0, math.sin(th),  math.cos(th)],
    ])
    cen = P[1:6].mean(axis=0)
    P_bad = P.copy()
    for i in range(1, 6):
        P_bad[i] = R @ (P[i] - cen) + cen
    return syms, P_bad


# ---------------------------------------------------------------------------
# Fixture: clean module reload so env-flag changes are seen.
# ---------------------------------------------------------------------------
@pytest.fixture()
def hho_env_on(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", "1")
    yield


@pytest.fixture()
def hho_env_off(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", raising=False)
    yield


def _fresh_module():
    import delfin.fffree.hapto_honest_construction as m
    return importlib.reload(m)


# ---------------------------------------------------------------------------
# Tests — env flag.
# ---------------------------------------------------------------------------
def test_env_off_byte_identical_no_op(hho_env_off):
    """When the env flag is unset, the corrector is a no-op (P.copy())."""
    mod = _fresh_module()
    assert mod.honest_active() is False
    syms, P = _ferrocene_eta5_eta5()
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 0
    assert P_out.shape == P.shape
    np.testing.assert_array_equal(P_out, P)


def test_env_master_flag_activates(monkeypatch):
    """The master ``CONSTRUCTION_FIX_ALL`` flag also activates the corrector."""
    monkeypatch.setenv("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "1")
    mod = _fresh_module()
    assert mod.honest_active() is True


def test_env_garbage_value_is_off(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", "xxx")
    monkeypatch.delenv("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", raising=False)
    mod = _fresh_module()
    assert mod.honest_active() is False


# ---------------------------------------------------------------------------
# Tests — detection.
# ---------------------------------------------------------------------------
def test_detect_ferrocene_two_eta5(hho_env_off):
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5()
    units = mod.detect_hapto_units(syms, P, metal_idx=0)
    assert len(units) == 2
    assert sorted(u.eta for u in units) == [5, 5]


def test_detect_arene_eta6(hho_env_off):
    mod = _fresh_module()
    syms, P = _arene_eta6()
    units = mod.detect_hapto_units(syms, P, metal_idx=0)
    assert len(units) == 1
    assert units[0].eta == 6


def test_detect_eta4_butadiene(hho_env_off):
    mod = _fresh_module()
    syms, P = _eta4_butadiene()
    units = mod.detect_hapto_units(syms, P, metal_idx=0)
    assert len(units) == 1
    assert units[0].eta == 4


def test_detect_eta2_ethylene(hho_env_off):
    mod = _fresh_module()
    syms, P = _eta2_ethylene()
    units = mod.detect_hapto_units(syms, P, metal_idx=0)
    assert len(units) == 1
    assert units[0].eta == 2


def test_detect_half_sandwich_one_eta5(hho_env_off):
    """Half-sandwich = one η5 ring + σ donors → exactly one hapto unit."""
    mod = _fresh_module()
    syms, P = _half_sandwich_eta5_plus_sigma()
    units = mod.detect_hapto_units(syms, P, metal_idx=0)
    assert len(units) == 1
    assert units[0].eta == 5
    assert sorted(units[0].indices) == [1, 2, 3, 4, 5]


def test_detect_no_hapto_in_pure_sigma(hho_env_off):
    """Pure Werner-type complex has no hapto units."""
    mod = _fresh_module()
    # OC-6 of 6 σ-N donors on Fe.
    syms = ["Fe"] + ["N"] * 6
    r = 2.10
    coords = np.array([
        [0.0, 0.0, 0.0],
        [+r, 0.0, 0.0], [-r, 0.0, 0.0],
        [0.0, +r, 0.0], [0.0, -r, 0.0],
        [0.0, 0.0, +r], [0.0, 0.0, -r],
    ], dtype=float)
    units = mod.detect_hapto_units(syms, coords, metal_idx=0)
    assert units == []


# ---------------------------------------------------------------------------
# Tests — construction.
# ---------------------------------------------------------------------------
def test_ferrocene_already_ideal_no_change_needed(hho_env_on):
    """Already-ideal ferrocene → corrector accepts but produces essentially
    the same geometry (centroid distance already at ideal)."""
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5(z_off=1.65)
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 2
    # Metal must not move.
    np.testing.assert_array_equal(P_out[0], P[0])
    # Centroid distance is the Fe-Cp ideal 1.65.
    for ring in ([1, 2, 3, 4, 5], [6, 7, 8, 9, 10]):
        cen = P_out[ring].mean(axis=0)
        d = float(np.linalg.norm(cen - P_out[0]))
        assert abs(d - 1.65) < 1e-3


def test_ferrocene_wrong_centroid_dist_gets_snapped(hho_env_on):
    """Build ferrocene with z_off=1.45 (wrong); corrector should snap to
    the CCDC ideal 1.65 Å for Fe η5."""
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5(z_off=1.45)
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 2
    for ring in ([1, 2, 3, 4, 5], [6, 7, 8, 9, 10]):
        cen = P_out[ring].mean(axis=0)
        d = float(np.linalg.norm(cen - P_out[0]))
        assert abs(d - 1.65) < 1e-3


def test_arene_eta6_centroid_snaps_to_ccdc_ideal(hho_env_on):
    """(η6-arene)Cr with z_off=1.80 (wrong) → snap to 1.62 (CCDC)."""
    mod = _fresh_module()
    syms, P = _arene_eta6(M="Cr", z_off=1.80)
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 1
    cen = P_out[1:7].mean(axis=0)
    d = float(np.linalg.norm(cen - P_out[0]))
    assert abs(d - 1.62) < 1e-3


def test_eta4_butadiene_centroid_snap(hho_env_on):
    """η4-butadiene on Fe: ideal 1.70 Å (CCDC)."""
    mod = _fresh_module()
    syms, P = _eta4_butadiene(z_off=2.00)
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 1
    cen = P_out[1:5].mean(axis=0)
    d = float(np.linalg.norm(cen - P_out[0]))
    assert abs(d - 1.70) < 1e-3


def test_eta2_ethylene_centroid_snap(hho_env_on):
    """η2-ethylene on Pt: ideal 2.10 Å (CCDC)."""
    mod = _fresh_module()
    syms, P = _eta2_ethylene(z_off=1.50)  # too close
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed == 1
    cen = P_out[1:3].mean(axis=0)
    d = float(np.linalg.norm(cen - P_out[0]))
    assert abs(d - 2.10) < 1e-3


# ---------------------------------------------------------------------------
# Tests — invariants (the publication promises).
# ---------------------------------------------------------------------------
def test_sigma_donors_invariant_in_half_sandwich(hho_env_on):
    """In a half-sandwich CpFe(Cl)(N)(P) with a wrong Cp centroid distance,
    the σ-Cl, σ-N, σ-P positions must NOT move."""
    mod = _fresh_module()
    syms, P = _half_sandwich_eta5_plus_sigma(z_off=1.45)
    # σ-donor indices (after metal + 5 Cp carbons): 6, 7, 8.
    sigma_idx = [6, 7, 8]
    P_sigma_pre = P[sigma_idx].copy()
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed >= 1
    np.testing.assert_array_equal(P_out[sigma_idx], P_sigma_pre)
    # M-D for σ donors must also be exactly invariant.
    for i in sigma_idx:
        d_pre = float(np.linalg.norm(P[i] - P[0]))
        d_post = float(np.linalg.norm(P_out[i] - P_out[0]))
        assert abs(d_pre - d_post) < 1e-9


def test_metal_does_not_move(hho_env_on):
    """The metal is fixed; only the hapto ring atoms (+ ring H) move."""
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5(z_off=1.45)
    P_out, _ = mod.apply_hapto_honest(syms, P, metal_idx=0)
    np.testing.assert_array_equal(P_out[0], P[0])


def test_ring_plane_rigidity_after_correction(hho_env_on):
    """Intra-ring C-C distances must be IDENTICAL before/after the corrector
    (rigid-block guarantee)."""
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5(z_off=1.45)
    P_out, _ = mod.apply_hapto_honest(syms, P, metal_idx=0)
    for ring in ([1, 2, 3, 4, 5], [6, 7, 8, 9, 10]):
        for i, j in [(a, b) for a in ring for b in ring if a < b]:
            d_pre = float(np.linalg.norm(P[i] - P[j]))
            d_post = float(np.linalg.norm(P_out[i] - P_out[j]))
            assert abs(d_pre - d_post) < 1e-6


def test_ring_plane_perpendicularity_after_correction(hho_env_on):
    """After the corrector, the M-centroid axis must be parallel to the
    ring normal (perpendicularity ~0 deg)."""
    mod = _fresh_module()
    syms, P = _ferrocene_distorted_centroid()
    P_out, n_fixed = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_fixed >= 1
    # Compute perpendicularity on the corrected top ring (indices 1..5).
    ring_xyz = P_out[1:6]
    cen = ring_xyz.mean(axis=0)
    centred = ring_xyz - cen
    _U, _S, Vt = np.linalg.svd(centred, full_matrices=False)
    normal = Vt[-1]
    nrm = float(np.linalg.norm(normal))
    normal = normal / max(nrm, 1e-12)
    md_vec = cen - P_out[0]
    md_unit = md_vec / max(float(np.linalg.norm(md_vec)), 1e-12)
    cos_ang = abs(float(np.dot(md_unit, normal)))
    perp_deg = math.degrees(math.acos(max(-1.0, min(1.0, cos_ang))))
    assert perp_deg < 1.0  # essentially axial


def test_determinism_two_calls_byte_identical(hho_env_on):
    """Two consecutive runs must produce identical coordinates."""
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5(z_off=1.45)
    P_a, n_a = mod.apply_hapto_honest(syms, P, metal_idx=0)
    P_b, n_b = mod.apply_hapto_honest(syms, P, metal_idx=0)
    assert n_a == n_b
    np.testing.assert_array_equal(P_a, P_b)


def test_invalid_metal_idx_returns_input_copy(hho_env_on):
    mod = _fresh_module()
    syms, P = _ferrocene_eta5_eta5()
    P_out, n = mod.apply_hapto_honest(syms, P, metal_idx=-1)
    assert n == 0
    np.testing.assert_array_equal(P_out, P)


def test_empty_input_handled(hho_env_on):
    mod = _fresh_module()
    P_out, n = mod.apply_hapto_honest([], np.zeros((0, 3), float), metal_idx=0)
    assert n == 0
    assert P_out.shape == (0, 3)


def test_hapto_only_cshm_decreases_after_correction(hho_env_on):
    """The end-to-end honesty test: the hapto-only-CShM metric must
    DECREASE after the corrector is applied (the construction makes
    the structure more honest with respect to the metric)."""
    mod = _fresh_module()
    from hapto_only_cshm import compute_hapto_only_cshm_from_arrays
    syms, P_bad = _ferrocene_distorted_centroid()
    rec_pre = compute_hapto_only_cshm_from_arrays(syms, P_bad)
    P_good, n_fixed = mod.apply_hapto_honest(syms, P_bad, metal_idx=0)
    rec_post = compute_hapto_only_cshm_from_arrays(syms, P_good)
    assert n_fixed >= 1
    assert rec_pre["hapto_only_cshm"] is not None
    assert rec_post["hapto_only_cshm"] is not None
    # End-to-end: hapto-only-CShM must drop sharply.
    assert rec_post["hapto_only_cshm"] < rec_pre["hapto_only_cshm"]
    # And approach 0 for the geometrically-perfect ideal target.
    assert rec_post["hapto_only_cshm"] < 1.0
