"""Tests for delfin.fffree.conformer_dedup — Kabsch-RMSD + Butina clustering
+ severity-best selection.
"""
from __future__ import annotations

import os

import numpy as np
import pytest


# ----------------------------------------------------------------------
# kabsch_rmsd
# ----------------------------------------------------------------------


def test_kabsch_rmsd_identity_zero():
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    assert kabsch_rmsd(P, P) < 1e-9


def test_kabsch_rmsd_translation_zero():
    """Pure translation -> kabsch RMSD = 0 (mean is subtracted)."""
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    Q = P + np.array([5.0, -2.0, 3.0])
    assert kabsch_rmsd(P, Q) < 1e-9


def test_kabsch_rmsd_rotation_zero():
    """Pure rotation -> kabsch RMSD = 0 (optimal rotation removes it)."""
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    # 90° around z
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    Q = P @ R.T
    assert kabsch_rmsd(P, Q) < 1e-9


def test_kabsch_rmsd_nontrivial_value():
    """A real deformation should yield a positive RMSD."""
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    Q = np.array([[0, 0, 0], [1.5, 0, 0], [0, 1, 0]], dtype=float)
    r = kabsch_rmsd(P, Q)
    assert r > 0.1


def test_kabsch_rmsd_shape_mismatch_inf():
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.array([[0, 0, 0]], dtype=float)
    Q = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    assert kabsch_rmsd(P, Q) == float("inf")


def test_kabsch_rmsd_empty_zero():
    from delfin.fffree.conformer_dedup import kabsch_rmsd
    P = np.zeros((0, 3))
    Q = np.zeros((0, 3))
    assert kabsch_rmsd(P, Q) == 0.0


# ----------------------------------------------------------------------
# pairwise_rmsd_matrix
# ----------------------------------------------------------------------


def test_pairwise_rmsd_matrix_symmetric():
    from delfin.fffree.conformer_dedup import pairwise_rmsd_matrix
    P = [
        np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float),
        np.array([[0, 0, 0], [1, 0, 0], [2, 0.5, 0]], dtype=float),
        np.array([[0, 0, 0], [0, 1, 0], [0, 2, 0]], dtype=float),
    ]
    M = pairwise_rmsd_matrix(P)
    assert M.shape == (3, 3)
    assert np.all(np.diag(M) < 1e-9)
    assert np.allclose(M, M.T, atol=1e-9)


def test_pairwise_rmsd_matrix_heavy_only_drops_h():
    from delfin.fffree.conformer_dedup import pairwise_rmsd_matrix
    P = [np.array([[0, 0, 0], [1, 0, 0], [99, 99, 99]], dtype=float),
         np.array([[0, 0, 0], [1, 0, 0], [-5, -5, -5]], dtype=float)]
    syms = ["C", "C", "H"]
    M = pairwise_rmsd_matrix(P, syms=syms, heavy_only=True)
    # H differs wildly but is dropped -> heavy-only RMSD ~ 0
    assert M[0, 1] < 1e-9


# ----------------------------------------------------------------------
# butina_cluster
# ----------------------------------------------------------------------


def test_butina_cluster_singleton_each():
    """Far-apart points -> each its own cluster."""
    from delfin.fffree.conformer_dedup import butina_cluster
    D = np.array([
        [0, 5, 5],
        [5, 0, 5],
        [5, 5, 0],
    ], dtype=float)
    clusters = butina_cluster(D, threshold=1.0)
    assert len(clusters) == 3


def test_butina_cluster_all_one_cluster():
    """Points all within threshold -> single cluster."""
    from delfin.fffree.conformer_dedup import butina_cluster
    D = np.array([
        [0, 0.1, 0.1],
        [0.1, 0, 0.1],
        [0.1, 0.1, 0],
    ], dtype=float)
    clusters = butina_cluster(D, threshold=0.5)
    assert len(clusters) == 1
    assert set(clusters[0]) == {0, 1, 2}


def test_butina_cluster_mixed():
    """Two near-points + one distinct -> 2 clusters."""
    from delfin.fffree.conformer_dedup import butina_cluster
    D = np.array([
        [0, 0.1, 5],
        [0.1, 0, 5],
        [5, 5, 0],
    ], dtype=float)
    clusters = butina_cluster(D, threshold=0.5)
    assert len(clusters) == 2


def test_butina_cluster_deterministic():
    """Repeat call -> identical cluster assignment."""
    from delfin.fffree.conformer_dedup import butina_cluster
    D = np.array([
        [0, 0.1, 5, 0.2],
        [0.1, 0, 5, 0.3],
        [5, 5, 0, 5],
        [0.2, 0.3, 5, 0],
    ], dtype=float)
    a = butina_cluster(D, threshold=0.5)
    b = butina_cluster(D, threshold=0.5)
    assert a == b


def test_butina_cluster_empty():
    from delfin.fffree.conformer_dedup import butina_cluster
    assert butina_cluster(np.zeros((0, 0)), 0.3) == []


def test_butina_cluster_single_point():
    from delfin.fffree.conformer_dedup import butina_cluster
    assert butina_cluster(np.zeros((1, 1)), 0.3) == [[0]]


# ----------------------------------------------------------------------
# dedup_by_rmsd — top-level
# ----------------------------------------------------------------------


def test_dedup_by_rmsd_basic():
    """Three frames, two near-duplicates -> dedup yields 2.

    Use truly non-superimposable triangle (Kabsch can rotate two collinear
    sets onto each other, so we need a 3D shape).
    """
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    P0 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    P1 = P0 + np.array([0.02, 0.0, 0.0])
    # P2 = mirror image (-x) — Kabsch (proper rotation only) cannot map this
    # onto P0 without a reflection -> RMSD > 0.
    P2 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, -1]], dtype=float)
    frames = [("A", P0, 1.0), ("B", P1, 2.0), ("C", P2, 1.5)]
    kept = dedup_by_rmsd(frames, threshold=0.3)
    assert len(kept) == 2, f"expected 2 clusters, got {len(kept)}"


def test_dedup_by_rmsd_keeps_lowest_severity():
    """Cluster keeper = minimum severity within the cluster."""
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    P0 = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    P1 = P0 + 0.01
    P2 = P0 + 0.02
    frames = [("A", P0, 5.0), ("B", P1, 1.0), ("C", P2, 3.0)]
    kept = dedup_by_rmsd(frames, threshold=0.3)
    assert len(kept) == 1
    assert kept[0][0] == "B", f"expected B (sev=1), got {kept[0][0]}"


def test_dedup_by_rmsd_sorted_by_severity():
    """Output is sorted by severity ascending."""
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    P_a = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    P_b = np.array([[0, 0, 0], [0, 1, 0]], dtype=float)
    P_c = np.array([[0, 0, 0], [0, 0, 1]], dtype=float)
    frames = [("A", P_a, 5.0), ("B", P_b, 2.0), ("C", P_c, 1.0)]
    kept = dedup_by_rmsd(frames, threshold=0.1)
    sevs = [s for _, _, s in kept]
    assert sevs == sorted(sevs), f"output not sorted ascending: {sevs}"


def test_dedup_by_rmsd_max_keep_honored():
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    # Use 4 atoms in a non-symmetric tetrahedron shape; per-frame perturb a
    # different atom so each frame is geometrically distinct under Kabsch.
    frames = []
    base = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    for k in range(10):
        Pk = base.copy()
        Pk[3] = np.array([float(k) * 2.0, float(k) * 2.0, float(k) * 2.0])
        frames.append((f"F{k}", Pk, float(k)))
    kept = dedup_by_rmsd(frames, threshold=0.3, max_keep=3)
    assert len(kept) == 3


def test_dedup_by_rmsd_deterministic():
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    P0 = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    P1 = P0 + 0.05
    P2 = P0 + np.array([0.5, 0.5, 0.0])
    frames = [("A", P0, 1.0), ("B", P1, 2.0), ("C", P2, 3.0)]
    a = dedup_by_rmsd(frames, threshold=0.3)
    b = dedup_by_rmsd(frames, threshold=0.3)
    assert len(a) == len(b)
    for (la, _, sa), (lb, _, sb) in zip(a, b):
        assert la == lb
        assert sa == sb


def test_dedup_by_rmsd_empty():
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    assert dedup_by_rmsd([]) == []


def test_dedup_by_rmsd_single_frame():
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    P = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    out = dedup_by_rmsd([("A", P, 1.0)], threshold=0.3)
    assert len(out) == 1


def test_dedup_by_rmsd_frame_dataclass():
    """Frame dataclass round-trip preserves shape."""
    from delfin.fffree.conformer_dedup import dedup_by_rmsd, Frame
    P0 = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    P1 = P0 + 0.02
    out = dedup_by_rmsd(
        [Frame("A", P0, 1.0), Frame("B", P1, 2.0)],
        threshold=0.3,
    )
    assert len(out) == 1
    assert isinstance(out[0], Frame)
    assert out[0].label == "A", "lowest-severity Frame should be kept"


# ----------------------------------------------------------------------
# Integration check on a realistic ARABUR-style case
# ----------------------------------------------------------------------


def test_dedup_arabur_style_67pct_redundancy():
    """Simulate ARABUR-style: 9 frames where 7 are within 0.3 Å and 2 are
    distinct -> dedup should yield exactly 3 essential conformers.

    A constant translation is removed by Kabsch centering, so "distinct"
    conformers here are constructed by SCALING + LOCAL DEFORMATION rather
    than rigid-body translation.
    """
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    base = np.random.RandomState(42).randn(20, 3)
    frames = []
    # 7 near-duplicates of base, severity 1..7 (small per-atom jitter only)
    for k in range(7):
        Pk = base + np.random.RandomState(k).randn(20, 3) * 0.02
        frames.append((f"near_{k}", Pk, float(k + 1)))
    # 1 distinct conformer: large local deformation on half the atoms
    base2 = base.copy()
    base2[10:] = base2[10:] + np.array([3.0, 3.0, 3.0])
    frames.append(("distinct_1", base2, 10.0))
    # 1 distinct conformer: scaled apart
    base3 = base * 2.5
    frames.append(("distinct_2", base3, 11.0))
    kept = dedup_by_rmsd(frames, threshold=0.3)
    # 1 cluster from the 7 near-duplicates + 2 singletons = 3
    assert len(kept) == 3, f"ARABUR-style expected 3 essential, got {len(kept)}"
    # The kept "near" should be the one with the lowest severity (=1).
    assert kept[0][0] == "near_0", f"expected near_0 first, got {kept[0][0]}"


# ----------------------------------------------------------------------
# Subagent #129 follow-up: default threshold 0.15 + symmetric-rep + 
# preserve-originals
# ----------------------------------------------------------------------


def test_default_threshold_now_015():
    """Spec change: env-default RMSD threshold tightened to 0.15 Å (was 0.30)."""
    from delfin.fffree.conformer_dedup import _env_threshold, DEFAULT_RMSD_THRESHOLD
    saved = os.environ.pop("DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD", None)
    try:
        assert DEFAULT_RMSD_THRESHOLD == 0.15
        assert _env_threshold() == 0.15
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD"] = saved


def test_env_threshold_override_honored():
    """Setting the env-var overrides the default."""
    from delfin.fffree.conformer_dedup import _env_threshold
    saved = os.environ.get("DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD")
    try:
        os.environ["DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD"] = "0.42"
        assert abs(_env_threshold() - 0.42) < 1e-9
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD", None)


def test_point_group_order_methane_high_symmetry():
    """Methane CH4 -> point group order >= 4 (D2d at minimum, ideally 24 Td)."""
    from delfin.fffree.conformer_dedup import point_group_order
    # Tetrahedral CH4 coordinates (perfect Td).
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ])
    syms = ["C", "H", "H", "H", "H"]
    order = point_group_order(P, syms)
    # Custom detector finds at least 2 (mirror plane / inversion).
    # spglib (when present) finds 24.  Both > 1.
    assert order >= 2


def test_point_group_order_c1_is_one():
    """A scrambled three-atom molecule has order 1."""
    from delfin.fffree.conformer_dedup import point_group_order
    P = np.array([[0.0, 0.0, 0.0], [1.234, 0.567, 0.0], [-0.789, 1.012, 0.5]])
    syms = ["C", "N", "O"]
    order = point_group_order(P, syms)
    assert order == 1


def test_select_symmetric_representative_prefers_high_symmetry():
    """Given two cluster members with equal severity, the higher-symmetry
    one wins."""
    from delfin.fffree.conformer_dedup import select_symmetric_representative
    # frame 0: scrambled -> order 1 (C1)
    P0 = np.array([[0.0, 0.0, 0.0], [1.2, 0.5, 0.1], [-0.7, 1.0, 0.3]])
    # frame 1: collinear, equidistant -> mirror plane -> order >= 2
    P1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    coords = [P0, P1]
    sev = [1.0, 1.0]
    syms = ["C", "C", "C"]
    best = select_symmetric_representative([0, 1], coords, sev, syms)
    assert best == 1, f"expected frame 1 (higher symmetry), got {best}"


def test_select_symmetric_representative_tradeoff_with_severity():
    """High symmetry + bad severity loses to low symmetry + great severity
    when the severity gap is huge (1/(1+sev) dominates)."""
    from delfin.fffree.conformer_dedup import select_symmetric_representative
    # frame 0: low symmetry, low severity (perfect)
    P0 = np.array([[0.0, 0.0, 0.0], [1.2, 0.5, 0.1], [-0.7, 1.0, 0.3]])
    # frame 1: high symmetry, ridiculous severity
    P1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    coords = [P0, P1]
    sev = [0.0, 0.0]  # equal severity -> symmetry wins
    syms = ["C", "C", "C"]
    best = select_symmetric_representative([0, 1], coords, sev, syms)
    assert best == 1
    # Sanity-check: score = pg_order + 1/(1+sev).  pg_order >= 2 (often 6/16
    # via custom/spglib) and the severity weight is 1/(1+sev) <= 1.  So a
    # highly symmetric frame ALWAYS beats a low-symmetry frame regardless of
    # severity -- the symmetric-rep flag is intentionally symmetry-dominant.
    sev2 = [0.0, 1e6]
    best2 = select_symmetric_representative([0, 1], coords, sev2, syms)
    assert best2 == 1, (
        "high-symmetry frame must win even at huge severity (score is "
        "symmetry-dominated by design)"
    )


def test_dedup_by_rmsd_symmetric_rep_kwarg():
    """symmetric_rep=True selects via score; default = legacy severity."""
    from delfin.fffree.conformer_dedup import dedup_by_rmsd
    # Three frames in one cluster (all near-duplicates).
    P0 = np.array([[0.0, 0.0, 0.0], [1.2, 0.5, 0.1], [-0.7, 1.0, 0.3]])
    P1 = P0 + 0.05  # near-duplicate (translation only -> kabsch=0)
    P2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    syms = ["C", "C", "C"]
    # Make P2 close enough to P0 to be in the same cluster (need a high
    # threshold; the Kabsch-RMSD between P0 and P2 is ~0.3).  Use 1.0 to
    # force-bag them together.
    frames = [("a", P0, 0.0), ("b", P1, 0.0), ("sym", P2, 0.5)]
    # symmetric_rep=False -> min severity (=0.0); first tied = "a"
    kept_legacy = dedup_by_rmsd(frames, threshold=1.0, syms=syms,
                                 symmetric_rep=False)
    assert len(kept_legacy) == 1
    assert kept_legacy[0][0] == "a"
    # symmetric_rep=True -> "sym" has higher pg_order >= 2 vs C1 P0/P1
    # (score = order + 1/(1+0.5) vs 1 + 1/(1+0) = 2 for legacy).
    # The symmetry detector must find order >= 2 for P2 to win.
    kept_sym = dedup_by_rmsd(frames, threshold=1.0, syms=syms,
                              symmetric_rep=True)
    assert len(kept_sym) == 1
    # "sym" wins iff order(P2) >= 2; the custom detector finds order 2
    # (mirror plane). So we expect kept[0][0] == "sym".
    assert kept_sym[0][0] == "sym"


def test_dedup_by_rmsd_preserve_originals_keeps_all_origs():
    """Originals (first N frames) are never dropped, even if dups.

    Uses chiral 4-atom geometries so Kabsch-RMSD picks up the difference
    (proper-rotation-only Kabsch flips chirality but not the chirality of
    a tetrahedral set).
    """
    from delfin.fffree.conformer_dedup import dedup_by_rmsd_preserve_originals
    P0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, 1.0]])
    P1 = P0 + 0.01  # near-duplicate of P0
    P2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, -1.0]])
    syms = ["C", "C", "C", "C"]
    # 2 originals (P0, P1) -- they cluster together but both kept.
    # + 1 additive (P2) -- distinct cluster (the third atom flipped sign).
    frames = [("orig0", P0, 0.0), ("orig1", P1, 1.0), ("add", P2, 0.5)]
    kept = dedup_by_rmsd_preserve_originals(
        frames, n_originals=2, threshold=0.3, syms=syms,
    )
    labels = [k[0] for k in kept]
    assert "orig0" in labels
    assert "orig1" in labels
    assert "add" in labels   # in its own cluster -> kept
    assert len(kept) == 3


def test_dedup_by_rmsd_preserve_originals_busy_cluster_origs_kept():
    """When originals already cover a cluster, additives in that cluster
    are dropped; additives in a DIFFERENT cluster are kept.

    Uses chiral 4-atom geometries so Kabsch-RMSD picks up the difference.
    """
    from delfin.fffree.conformer_dedup import dedup_by_rmsd_preserve_originals
    P0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, 1.0]])
    P1 = P0 + 0.01   # near-dup, same cluster as P0
    P2 = P0 + 0.02   # near-dup, same cluster as P0
    P3 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, -1.0]])
    syms = ["C", "C", "C", "C"]
    # 1 original (P0) + 3 additives (P1, P2, P3).
    frames = [("orig", P0, 0.0), ("add1", P1, 0.5),
              ("add2", P2, 0.4), ("add3", P3, 0.6)]
    kept = dedup_by_rmsd_preserve_originals(
        frames, n_originals=1, threshold=0.3, syms=syms,
    )
    labels = [k[0] for k in kept]
    assert "orig" in labels
    # add1/add2 are in the same cluster as orig -> dropped.
    assert "add1" not in labels
    assert "add2" not in labels
    # add3 is in its own cluster -> kept.
    assert "add3" in labels


def test_dedup_by_rmsd_preserve_originals_zero_originals_falls_back():
    """n_originals=0 -> behaves like the legacy dedup_by_rmsd (best per cluster).

    Uses chiral 4-atom geometries so Kabsch picks up the cluster distinction.
    """
    from delfin.fffree.conformer_dedup import (
        dedup_by_rmsd, dedup_by_rmsd_preserve_originals,
    )
    P0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, 1.0]])
    P1 = P0 + 0.05
    P2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.5, 1.0, 0.0], [0.5, 0.0, -1.0]])
    syms = ["C", "C", "C", "C"]
    frames = [("a", P0, 1.0), ("b", P1, 0.5), ("c", P2, 0.7)]
    kept_legacy = dedup_by_rmsd(frames, threshold=0.3, syms=syms)
    kept_preserve = dedup_by_rmsd_preserve_originals(
        frames, n_originals=0, threshold=0.3, syms=syms,
    )
    assert len(kept_legacy) == len(kept_preserve) == 2


def test_pre_cluster_emit_env_default_off():
    """DELFIN_FFFREE_PRE_CLUSTER_EMIT default = OFF (byte-identical to HEAD)."""
    from delfin.fffree.conformer_dedup import _env_pre_cluster_emit
    saved = os.environ.pop("DELFIN_FFFREE_PRE_CLUSTER_EMIT", None)
    try:
        assert _env_pre_cluster_emit() is False
        os.environ["DELFIN_FFFREE_PRE_CLUSTER_EMIT"] = "1"
        assert _env_pre_cluster_emit() is True
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_PRE_CLUSTER_EMIT"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_PRE_CLUSTER_EMIT", None)


def test_symmetric_rep_env_default_off():
    """DELFIN_FFFREE_SYMMETRIC_REP default = OFF (byte-identical to HEAD)."""
    from delfin.fffree.conformer_dedup import _env_symmetric_rep
    saved = os.environ.pop("DELFIN_FFFREE_SYMMETRIC_REP", None)
    try:
        assert _env_symmetric_rep() is False
        os.environ["DELFIN_FFFREE_SYMMETRIC_REP"] = "1"
        assert _env_symmetric_rep() is True
    finally:
        if saved is not None:
            os.environ["DELFIN_FFFREE_SYMMETRIC_REP"] = saved
        else:
            os.environ.pop("DELFIN_FFFREE_SYMMETRIC_REP", None)
