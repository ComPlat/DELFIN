"""Unit tests for Welle-5j A — Cp piano-stool hapticity refinement.

Validates ``delfin.manta._cp_piano_stool``:

* synthetic ferrocene-like η⁵ ring + off-axis Fe -> Fe snaps to axis at
  ideal η⁵ distance
* benzene-like 6-ring + Cr -> NOT detected as Cp (no false-positive on
  η⁶)
* 5-ring far from metal (>3 Å) -> NOT detected (outside cutoff)
* σ-only Pt(NH₃)₄ analog -> nothing moves (no 5-ring)
* M-D invariant rollback -> shift is reverted when another donor H/
  donor would drift > 0.05 Å
* env-flag default-OFF: dispatch helper in smiles_converter is
  bit-exact when ``DELFIN_5J_A_CP_PIANO_STOOL`` is unset / 0
"""
from __future__ import annotations

import os
from typing import List, Tuple

import numpy as np
import pytest

cp = pytest.importorskip(
    "delfin.manta._cp_piano_stool",
    reason="_cp_piano_stool module not present",
)


def _ferrocene_like(metal_sym: str = "Fe",
                   metal_pos: Tuple[float, float, float] = (0.0, 0.0, 1.65),
                   ) -> Tuple[List[str], np.ndarray]:
    """Build a planar Cp ring (5 C in z=0) + 5 H + metal at given pos."""
    syms = ["C"] * 5 + ["H"] * 5 + [metal_sym]
    pts: List[np.ndarray] = []
    r_c = 1.20
    r_h = r_c + 1.08
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_c * np.cos(a), r_c * np.sin(a), 0.0]))
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_h * np.cos(a), r_h * np.sin(a), 0.0]))
    pts.append(np.array(metal_pos, dtype=float))
    return syms, np.vstack(pts)


def test_ferrocene_off_axis_metal_snaps_to_eta5_axis() -> None:
    """Fe placed off the Cp axis is snapped to (0, 0, +1.65 Å)."""
    syms, pts = _ferrocene_like(metal_pos=(0.30, 0.15, 1.55))
    pts_work = pts.copy()
    moved = cp.refine_cp_piano_stool(syms, pts_work)
    assert moved == 1
    centroid = pts_work[:5].mean(axis=0)
    d = float(np.linalg.norm(pts_work[10] - centroid))
    assert abs(d - 1.65) < 1e-6


def test_benzene_arene_not_detected_as_cp() -> None:
    """A 6-ring of C does not trigger the η⁵ Cp snap."""
    syms = ["C"] * 6 + ["H"] * 6 + ["Cr"]
    pts: List[List[float]] = []
    r_c = 1.39
    r_h = r_c + 1.08
    for k in range(6):
        a = 2 * np.pi * k / 6
        pts.append([r_c * np.cos(a), r_c * np.sin(a), 0.0])
    for k in range(6):
        a = 2 * np.pi * k / 6
        pts.append([r_h * np.cos(a), r_h * np.sin(a), 0.0])
    pts.append([0.0, 0.0, 1.61])  # Cr on η⁶ axis
    arr = np.array(pts, dtype=float)
    pre = arr.copy()
    moved = cp.refine_cp_piano_stool(syms, arr)
    assert moved == 0
    assert np.allclose(arr, pre)


def test_5ring_far_from_metal_not_detected() -> None:
    """A 5-ring at 5 Å from the metal is outside the coordination cutoff."""
    syms, pts = _ferrocene_like(metal_pos=(0.0, 0.0, 5.0))
    pts_work = pts.copy()
    moved = cp.refine_cp_piano_stool(syms, pts_work)
    assert moved == 0
    assert np.allclose(pts_work, pts)


def test_sigma_only_complex_does_nothing() -> None:
    """A square-planar Pt(NH₃)₄ analog has no 5-ring -> no shift."""
    syms = ["Pt"] + ["N"] * 4 + ["H"] * 12
    pts: List[List[float]] = [[0.0, 0.0, 0.0],
                              [2.05, 0.0, 0.0], [-2.05, 0.0, 0.0],
                              [0.0, 2.05, 0.0], [0.0, -2.05, 0.0]]
    for nx in [(2.05, 0, 0), (-2.05, 0, 0), (0, 2.05, 0), (0, -2.05, 0)]:
        for k in range(3):
            a = 2 * np.pi * k / 3
            pts.append([nx[0] + 0.9 * np.cos(a),
                        nx[1] + 0.9 * np.sin(a), 0.3])
    arr = np.array(pts, dtype=float)
    pre = arr.copy()
    moved = cp.refine_cp_piano_stool(syms, arr)
    assert moved == 0
    assert np.allclose(arr, pre)


def test_md_drift_allowed_when_bond_stays_intact() -> None:
    """The M shift naturally moves M relative to non-ring donors by
    several tenths of an Å — those drifts are allowed as long as each
    bonded donor remains within the bonded cutoff (1.30 · Σr_cov)."""
    syms = ["C"] * 5 + ["H"] * 5 + ["Fe", "Cl"]
    pts: List[np.ndarray] = []
    r_c = 1.20
    r_h = r_c + 1.08
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_c * np.cos(a), r_c * np.sin(a), 0.0]))
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_h * np.cos(a), r_h * np.sin(a), 0.0]))
    # Fe off-axis but not catastrophically: simulates UFF-distorted Cp
    # where one ring carbon has drifted ~0.3 Å from the ideal cone.
    pts.append(np.array([0.40, 0.20, 1.55]))
    # Cl ~2.30 Å from this Fe (bonded), will stay bonded after snap.
    pts.append(np.array([2.30, 0.20, 1.55]))
    arr = np.vstack(pts)
    moved = cp.refine_cp_piano_stool(syms, arr)
    assert moved == 1
    post_d = float(np.linalg.norm(arr[10] - arr[11]))
    bonded_cutoff = cp._MD_INTACT_FACTOR * (cp._COV_RADII["Fe"]
                                            + cp._COV_RADII["Cl"])
    assert post_d <= bonded_cutoff


def test_md_invariant_blocks_dissociating_donor() -> None:
    """If the snap would pull a bonded donor past the 1.30·Σr_cov bonded
    cutoff (i.e. dissociate it), the whole shift is reverted."""
    syms = ["C"] * 5 + ["H"] * 5 + ["Fe", "Cl"]
    pts: List[np.ndarray] = []
    r_c = 1.20
    r_h = r_c + 1.08
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_c * np.cos(a), r_c * np.sin(a), 0.0]))
    for k in range(5):
        a = 2 * np.pi * k / 5
        pts.append(np.array([r_h * np.cos(a), r_h * np.sin(a), 0.0]))
    # Fe near the Cp axis (small offset so the candidate IS detected);
    # Cl placed so that after the snap the Fe-Cl bond would exceed the
    # 1.30·(1.32+1.02) = 3.04 Å bonded cutoff.
    pts.append(np.array([0.30, 0.20, 1.55]))
    # Fe-Cl current = sqrt(2.7² + 0² + 0²) = 2.70 Å (bonded, snapshot
    # records it).  After snap Fe moves to (0,0,1.65), new Fe-Cl =
    # sqrt(3.0² + 0.2² + 0.10²) = 3.01 — still inside cutoff.  Push Cl
    # further so post-snap dist > 3.04.
    pts.append(np.array([3.30, 0.20, 1.55]))
    arr = np.vstack(pts)
    pre = arr.copy()
    moved = cp.refine_cp_piano_stool(syms, arr)
    # Pre-shift Fe-Cl = 3.00 Å (inside 3.04 cutoff -> snapshotted).
    # Post-shift Fe-Cl = sqrt(3.30² + 0.20² + 0.10²) = 3.31 Å (past cutoff
    # -> dissociated).  Snap must revert.
    assert moved == 0
    assert np.allclose(arr, pre)


def test_correct_xyz_round_trip_on_non_metal_string() -> None:
    """A non-metal XYZ string is returned bit-exact unchanged."""
    xyz = "3\nbenzene fragment\nC 0.0 0.0 0.0\nC 1.4 0.0 0.0\nC 0.7 1.2 0.0\n"
    out = cp.correct_xyz(xyz)
    assert out == xyz


def test_correct_results_signature_parity_with_b3_b4() -> None:
    """``correct_results(mol, [(xyz, label)])`` preserves the tuple shape
    matching ``_coord_angle_corrector`` / ``_pi_h_projector`` contract."""
    xyz = "3\nfrag\nC 0.0 0.0 0.0\nO 1.2 0.0 0.0\nH 1.9 0.0 0.0\n"
    results = [(xyz, "label_1"), (xyz, "label_2")]
    out = cp.correct_results(None, results)
    assert len(out) == 2
    assert out[0][1] == "label_1"
    assert out[1][1] == "label_2"


def test_env_flag_default_off_dispatch_is_passthrough() -> None:
    """``_apply_5j_a_cp_piano_stool_if_enabled`` returns input unchanged
    when ``DELFIN_5J_A_CP_PIANO_STOOL`` is unset (default OFF)."""
    sc = pytest.importorskip("delfin.smiles_converter")
    dispatcher = getattr(sc, "_apply_5j_a_cp_piano_stool_if_enabled", None)
    assert dispatcher is not None
    saved = os.environ.pop("DELFIN_5J_A_CP_PIANO_STOOL", None)
    try:
        xyz = "1\nx\nFe 0.0 0.0 0.0\n"
        results = [(xyz, "lbl")]
        out = dispatcher(None, results, dual_parse_done=False)
        # Bit-exact identity when env-flag is OFF.
        assert out is results
    finally:
        if saved is not None:
            os.environ["DELFIN_5J_A_CP_PIANO_STOOL"] = saved


def test_env_flag_dispatch_skips_inner_dual_parse() -> None:
    """Even with the env-flag ON, the dispatcher must skip the inner
    dual-parse call (heavy-atom signature dedup invariant)."""
    sc = pytest.importorskip("delfin.smiles_converter")
    dispatcher = sc._apply_5j_a_cp_piano_stool_if_enabled
    os.environ["DELFIN_5J_A_CP_PIANO_STOOL"] = "1"
    try:
        xyz = "1\nx\nFe 0.0 0.0 0.0\n"
        results = [(xyz, "lbl")]
        out = dispatcher(None, results, dual_parse_done=True)
        assert out is results  # identity (no copy made)
    finally:
        os.environ.pop("DELFIN_5J_A_CP_PIANO_STOOL", None)


def test_module_constants_have_expected_values() -> None:
    """Sanity: the η=5 target distances cover the common ferrocene metals
    and the fallback is in the documented 1.85 Å range."""
    for m in ("Fe", "Ru", "Co", "Ni", "Cr", "Mn", "Ti", "Zr"):
        assert m in cp._ETA5_TARGET_MC
    assert 1.50 < cp._ETA5_TARGET_MC_FALLBACK < 2.30
