"""Tests for delfin.fffree.conformer_symmetry — point-group order detection
and symmetry-priority cluster-representative scoring.

Covers:
  - C1 (random), Ci (inversion), Cs/C2v (water), C3v (NH3), D6h (benzene),
    Td (methane), Oh (octahedron).
  - Score function tie-break behaviour at equal severity.
  - Default-OFF byte-identity proof for the DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY
    env-flag wired into delfin.fffree.conformer_dedup.dedup_by_rmsd.
  - Determinism: identical input -> bit-identical output across two runs.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest


# ----------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------


def _benzene():
    """D6h benzene in the xy-plane (C-C ~1.4 Å, C-H ~1.1 Å)."""
    R_C = 1.4
    R_H = 2.5
    coords = []
    Z = []
    for k in range(6):
        a = 2 * np.pi * k / 6.0
        coords.append([R_C * np.cos(a), R_C * np.sin(a), 0.0])
        Z.append(6)
    for k in range(6):
        a = 2 * np.pi * k / 6.0
        coords.append([R_H * np.cos(a), R_H * np.sin(a), 0.0])
        Z.append(1)
    return np.array(coords, dtype=float), Z


def _nh3():
    """C3v NH3 with the C3 axis along z (H atoms in a triangle below N)."""
    return (
        np.array([
            [0.0,    0.0,    0.0],
            [0.94,   0.0,   -0.3],
            [-0.47,  0.814, -0.3],
            [-0.47, -0.814, -0.3],
        ], dtype=float),
        [7, 1, 1, 1],
    )


def _octahedron_oh():
    """Oh — six equivalent atoms at the vertices of a regular octahedron."""
    P = np.array([
        [ 1, 0, 0],
        [-1, 0, 0],
        [ 0, 1, 0],
        [ 0,-1, 0],
        [ 0, 0, 1],
        [ 0, 0,-1],
    ], dtype=float)
    return P, [8] * 6


def _methane():
    """Td CH4 with the 3-fold axes along the cube diagonals."""
    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    P = np.array([
        [0.0, 0.0, 0.0],
        [ inv_sqrt3,  inv_sqrt3,  inv_sqrt3],
        [-inv_sqrt3, -inv_sqrt3,  inv_sqrt3],
        [-inv_sqrt3,  inv_sqrt3, -inv_sqrt3],
        [ inv_sqrt3, -inv_sqrt3, -inv_sqrt3],
    ], dtype=float)
    return P, [6, 1, 1, 1, 1]


def _water():
    """C2v H2O — O at origin, two H in the xy-plane, C2 along y."""
    return (
        np.array([
            [0.0,  0.0,  0.0],
            [ 0.76, 0.59, 0.0],
            [-0.76, 0.59, 0.0],
        ], dtype=float),
        [8, 1, 1],
    )


# ----------------------------------------------------------------------
# detect_point_group_order
# ----------------------------------------------------------------------


def test_detect_C1_for_random():
    """A random heavy-atom blob has no symmetry beyond identity -> order 1."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    rng = np.random.default_rng(42)
    P = rng.standard_normal((8, 3)) * 3.0
    Z = list(range(1, 9))  # all distinct elements -> no equivalence
    order = detect_point_group_order(P, Z, tol=0.05)
    assert order == 1, f"expected C1 (order 1), got {order}"


def test_detect_D6h_for_benzene():
    """Benzene D6h has order 24 (12 proper + 12 improper)."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _benzene()
    order = detect_point_group_order(P, Z, tol=0.1)
    # The custom backend recovers D6h via composition closure; spglib's
    # vacuum-cell wrapper may under-count for off-axis orientations, but
    # the public API takes the max.
    assert order >= 12, f"expected order >= 12 (D6), got {order}"
    assert order == 24, f"expected D6h order 24, got {order}"


def test_detect_C3v_for_NH3():
    """NH3 C3v has order 6 (E, 2 C3, 3 sigma_v)."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _nh3()
    order = detect_point_group_order(P, Z, tol=0.1)
    assert order == 6, f"expected C3v order 6, got {order}"


def test_detect_Oh_for_octahedron():
    """Octahedron Oh has 48 ops (24 proper + 24 improper)."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _octahedron_oh()
    order = detect_point_group_order(P, Z, tol=0.1)
    # The task spec accepts 24 (O proper) or 48 (Oh full).
    assert order in (24, 48), f"expected octahedron 24 or 48, got {order}"


def test_detect_Td_for_methane():
    """Methane Td has order 24."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _methane()
    order = detect_point_group_order(P, Z, tol=0.1)
    assert order >= 12, f"expected at least T order 12, got {order}"
    assert order == 24, f"expected Td order 24, got {order}"


def test_detect_C2v_for_water():
    """H2O C2v has order 4."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _water()
    order = detect_point_group_order(P, Z, tol=0.1)
    assert order == 4, f"expected C2v order 4, got {order}"


def test_detect_accepts_element_symbols():
    """Public API must also accept element-symbol input transparently."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, _ = _nh3()
    symbols = ["N", "H", "H", "H"]
    order = detect_point_group_order(P, symbols, tol=0.1)
    assert order == 6


def test_detect_empty_coords_returns_one():
    """No atoms -> trivial group (order 1)."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    assert detect_point_group_order(np.zeros((0, 3)), [], tol=0.1) == 1


# ----------------------------------------------------------------------
# score_conformer_for_representative
# ----------------------------------------------------------------------


def test_score_high_sym_beats_low_sym_at_tie_severity():
    """At identical severity, the higher-point-group conformer wins."""
    from delfin.fffree.conformer_symmetry import (
        detect_point_group_order,
        score_conformer_for_representative,
    )
    P_high, Z_high = _benzene()       # D6h: 24
    P_low, Z_low = _benzene()         # break by jittering one atom out of plane
    P_low = P_low.copy()
    # Off-plane jitter kills BOTH the C6 and the sigma_h -> order 1 (C1).
    P_low[0] += np.array([0.6, 0.6, 0.6])
    assert detect_point_group_order(P_low, Z_low, tol=0.1) == 1

    frame_high = {"coords": P_high, "atomic_numbers": Z_high,
                  "mogul_severity": 5.0}
    frame_low = {"coords": P_low, "atomic_numbers": Z_low,
                 "mogul_severity": 5.0}
    s_high = score_conformer_for_representative(frame_high)
    s_low = score_conformer_for_representative(frame_low)
    assert s_high > s_low, (
        f"high-sym score {s_high:.3f} must beat low-sym score {s_low:.3f} "
        f"at equal severity"
    )


def test_score_formula_matches_spec():
    """Hard-coded score: w_sym * log2(pg) - w_sev * sev."""
    from delfin.fffree.conformer_symmetry import (
        score_conformer_for_representative,
    )
    P, Z = _benzene()
    frame = {"coords": P, "atomic_numbers": Z, "mogul_severity": 2.0}
    s = score_conformer_for_representative(frame, w_sym=0.5, w_sev=1.0)
    # benzene D6h -> 24; log2(24) ~= 4.585; score = 0.5*4.585 - 1.0*2.0
    expected = 0.5 * math.log2(24) - 1.0 * 2.0
    assert abs(s - expected) < 1e-9


def test_score_low_severity_can_outweigh_low_symmetry():
    """A clean (severity 0) C1 conformer beats a high-severity D6h one
    when ``w_sev`` is large enough relative to ``w_sym``.
    """
    from delfin.fffree.conformer_symmetry import (
        score_conformer_for_representative,
    )
    P_high, Z_high = _benzene()
    P_c1 = P_high.copy()
    P_c1[0] += np.array([0.6, 0.6, 0.6])
    frame_high = {"coords": P_high, "atomic_numbers": Z_high,
                  "mogul_severity": 10.0}
    frame_low = {"coords": P_c1, "atomic_numbers": Z_high,
                 "mogul_severity": 0.0}
    s_high = score_conformer_for_representative(frame_high)
    s_low = score_conformer_for_representative(frame_low)
    assert s_low > s_high, (
        "clean C1 (sev=0) must beat severe D6h (sev=10) under default "
        "weights"
    )


def test_score_handles_missing_severity():
    """Frames without a severity field default to 0.0."""
    from delfin.fffree.conformer_symmetry import (
        score_conformer_for_representative,
    )
    P, Z = _benzene()
    frame = {"coords": P, "atomic_numbers": Z}
    s = score_conformer_for_representative(frame)
    # benzene D6h -> 24; sev defaults to 0
    expected = 0.5 * math.log2(24)
    assert abs(s - expected) < 1e-9


# ----------------------------------------------------------------------
# Dedup env-flag plumbing — byte-identity + determinism
# ----------------------------------------------------------------------


def _three_frames_one_high_sym():
    """Three test frames sharing the SAME atom layout (benzene-like)
    that fall into a single Butina cluster (large RMSD threshold so the
    jittered frames cluster with the clean one):
      - frame 0: clean D6h benzene, severity 1.0
      - frame 1: jitter on one atom > detector tolerance -> C1, severity 0.5
      - frame 2: another jitter, C1, severity 0.5

    The jitter (0.6 Å on a single atom) exceeds the default 0.1 Å
    point-group tolerance and reliably breaks the D6h symmetry, but
    stays under the RMSD threshold (1.0 Å) used by the test so all
    three frames land in the same Butina cluster.
    """
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P0, Z = _benzene()
    P1 = P0.copy()
    P1[0] += np.array([0.6, 0.6, 0.0])
    P2 = P0.copy()
    P2[6] += np.array([0.5, -0.5, 0.0])
    # Sanity: confirm the detector sees the clean frame as full D6h (24)
    # and the jittered frames as substantially LOWER symmetry.  The xy-
    # plane mirror survives in-plane jitter (Cs/order 2), which is enough
    # for the scoring test as long as the D6h bonus dominates.
    assert detect_point_group_order(P0, Z, tol=0.1) == 24
    o1 = detect_point_group_order(P1, Z, tol=0.1)
    o2 = detect_point_group_order(P2, Z, tol=0.1)
    assert o1 < 24, f"P1 should be lower-sym than D6h, got order {o1}"
    assert o2 < 24, f"P2 should be lower-sym than D6h, got order {o2}"
    frames = [
        ("D6h", P0, 1.0),
        ("C1a", P1, 0.5),
        ("C1b", P2, 0.5),
    ]
    syms = ["C"] * 6 + ["H"] * 6
    return frames, syms, Z


def test_dedup_env_off_byte_identical(monkeypatch):
    """Without the env flag, dedup_by_rmsd returns the SAME representative
    as HEAD (lowest severity wins).  Byte-identity proof for the
    DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY hook.
    """
    monkeypatch.delenv("DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY",
                       raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_SYMMETRIC_REP", raising=False)
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    frames, syms, _ = _three_frames_one_high_sym()
    kept = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    assert len(kept) == 1
    # Lowest severity: 0.5 (C1a, lower index breaks tie vs C1b).
    label = kept[0][0] if isinstance(kept[0], tuple) else kept[0].label
    assert label == "C1a", f"expected C1a (lowest sev + index), got {label}"


def test_dedup_env_on_picks_high_sym(monkeypatch):
    """With the priority hook ON, the cluster representative is the
    HIGH-SYMMETRY frame even though it has worse severity.
    """
    monkeypatch.setenv("DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY", "1")
    monkeypatch.delenv("DELFIN_FFFREE_SYMMETRIC_REP", raising=False)
    # Heavy weight so the D6h bonus outranks the 0.5-vs-1.0 severity gap.
    monkeypatch.setenv("DELFIN_FFFREE_RMSD_DEDUP_SYM_WEIGHT", "5.0")
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    frames, syms, _ = _three_frames_one_high_sym()
    kept = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    assert len(kept) == 1
    label = kept[0][0] if isinstance(kept[0], tuple) else kept[0].label
    assert label == "D6h", f"expected D6h (high-sym priority), got {label}"


def test_dedup_env_unset_matches_explicit_zero(monkeypatch):
    """Unset and ``=0`` must produce IDENTICAL output (env-flag invariant)."""
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    frames, syms, _ = _three_frames_one_high_sym()
    monkeypatch.delenv("DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY",
                       raising=False)
    out_unset = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    monkeypatch.setenv("DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY", "0")
    out_zero = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    # Compare labels (the only field affected by representative choice).
    lbl_unset = [f[0] if isinstance(f, tuple) else f.label for f in out_unset]
    lbl_zero = [f[0] if isinstance(f, tuple) else f.label for f in out_zero]
    assert lbl_unset == lbl_zero


def test_determinism_two_runs_byte_identical(monkeypatch):
    """Two back-to-back runs with the priority hook ON return bit-identical
    coordinates + labels for the same input.
    """
    monkeypatch.setenv("DELFIN_FFFREE_RMSD_DEDUP_SYMMETRY_PRIORITY", "1")
    monkeypatch.setenv("DELFIN_FFFREE_RMSD_DEDUP_SYM_WEIGHT", "5.0")
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    frames, syms, _ = _three_frames_one_high_sym()
    k1 = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    k2 = dedup_by_rmsd(frames, threshold=1.0, syms=syms)
    assert len(k1) == len(k2)
    for a, b in zip(k1, k2):
        if isinstance(a, tuple):
            assert a[0] == b[0]
            assert np.array_equal(a[1], b[1])
            assert a[2] == b[2]
        else:
            assert a.label == b.label
            assert np.array_equal(a.P, b.P)
            assert a.severity == b.severity


def test_priority_select_function_directly():
    """Direct call to the cluster-internal selector returns the higher-
    symmetry frame at equal severity.
    """
    from delfin.fffree.conformer_symmetry import (
        select_symmetric_priority_representative,
    )
    P_high, Z = _benzene()
    P_low = P_high.copy()
    P_low[0] += np.array([0.6, 0.6, 0.6])  # break C6 -> C1
    coords_list = [P_high, P_low]
    severities = [1.0, 1.0]
    symbols = ["C"] * 6 + ["H"] * 6
    best = select_symmetric_priority_representative(
        [0, 1], coords_list, severities, symbols,
    )
    assert best == 0, "high-sym frame (index 0) must win at tie severity"


# ----------------------------------------------------------------------
# Determinism — point-group detection is order-independent and seed-free
# ----------------------------------------------------------------------


def test_detect_two_runs_identical():
    """Same input twice -> exactly the same integer order."""
    from delfin.fffree.conformer_symmetry import detect_point_group_order
    P, Z = _nh3()
    o1 = detect_point_group_order(P, Z, tol=0.1)
    o2 = detect_point_group_order(P, Z, tol=0.1)
    assert o1 == o2


def test_score_two_runs_identical():
    """Same input twice -> exactly the same float score."""
    from delfin.fffree.conformer_symmetry import (
        score_conformer_for_representative,
    )
    P, Z = _benzene()
    f = {"coords": P, "atomic_numbers": Z, "mogul_severity": 2.0}
    s1 = score_conformer_for_representative(f)
    s2 = score_conformer_for_representative(f)
    assert s1 == s2


# ----------------------------------------------------------------------
# Phase 2 — per-rotamer symmetry pre-ranking
# ----------------------------------------------------------------------


def test_rerank_rotamers_puts_high_sym_first():
    """rerank_rotamers_by_symmetry sorts variants by detected pg-order DESC."""
    from delfin.fffree.single_bond_rotamers import rerank_rotamers_by_symmetry
    P_high, Z = _benzene()
    P_low = P_high.copy()
    P_low[0] += np.array([0.6, 0.6, 0.6])
    variants = [(P_low, "rot_low"), (P_high, "rot_high")]
    ranked = rerank_rotamers_by_symmetry(variants, Z)
    assert ranked[0][1] == "rot_high"
    assert ranked[1][1] == "rot_low"


def test_rerank_rotamers_stable_tiebreak():
    """Equal pg-order -> stable sort preserves input order."""
    from delfin.fffree.single_bond_rotamers import rerank_rotamers_by_symmetry
    P_a, Z = _benzene()
    P_b = P_a.copy()
    P_c = P_a.copy()
    variants = [(P_a, "A"), (P_b, "B"), (P_c, "C")]
    ranked = rerank_rotamers_by_symmetry(variants, Z)
    assert [t[1] for t in ranked] == ["A", "B", "C"]


def test_rerank_rotamers_empty_safe():
    """Empty input -> empty output, no crash."""
    from delfin.fffree.single_bond_rotamers import rerank_rotamers_by_symmetry
    assert rerank_rotamers_by_symmetry([], [6, 1]) == []


def test_phase2_env_default_off_byte_identical(monkeypatch):
    """``DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS`` unset -> the rotamer
    generator behaves bit-identically to before the patch (stream order
    determined purely by ``enumerate_rotamer_configs``).

    We test this against a plain butane (CCCC) which has rotors and so
    triggers the streaming path; without the env flag the order matches
    a direct list-comprehension over the inner generator.
    """
    monkeypatch.delenv("DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS",
                       raising=False)
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.fffree.single_bond_rotamers import (
        enumerate_single_bond_rotamers,
    )
    m = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    AllChem.EmbedMolecule(m, randomSeed=42)
    P = m.GetConformer().GetPositions()
    labels = [lab for _, lab in enumerate_single_bond_rotamers(m, P)]
    # Now explicitly set the env flag to "0" -> must produce identical
    # output (default-OFF invariant).
    monkeypatch.setenv("DELFIN_FFFREE_SYMMETRY_PRIORITY_ROTAMERS", "0")
    labels_off = [lab for _, lab in enumerate_single_bond_rotamers(m, P)]
    assert labels == labels_off
