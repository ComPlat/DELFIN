"""Tests for delfin.fffree.mddir_alignment

Validate per-donor M-D direction alignment by axis rotation:
  - Default OFF byte-identical.
  - Aromatic donor face-on -> 90 deg rotation around M-D -> ring perpendicular.
  - sp3-anti donor (M on substituent side) -> 180 deg rotation -> M anti.
  - donor_linearized SKIPPED (n_skipped reported).
  - M-D invariant exactly preserved (axis rotation through donor).
  - No donor-element atom: no-op.
"""
import math
import os
import numpy as np
import pytest


def _reset_env():
    for k in (
        "DELFIN_FFFREE_MDDIR_ALIGN",
        "DELFIN_FFFREE_CONSTRUCTION_FIX_ALL",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


def test_default_off_byte_identical():
    """No env flag set: coordinates unchanged."""
    from delfin.fffree.mddir_alignment import align_donor_lonepairs

    # Pd + aromatic-N donor (pyridine-like)
    syms = ["Pd", "N", "C", "C", "C", "C", "C", "H", "H", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Pd
        [2.10, 0.0, 0.0],    # N (donor)
        [2.80, 1.20, 0.0],   # C neighbour 1 in plane
        [2.80, -1.20, 0.0],  # C neighbour 2 in plane
        [4.20, 1.20, 0.0],   # ring C
        [4.20, -1.20, 0.0],  # ring C
        [4.90, 0.0, 0.0],    # ring C
        [4.75, 2.10, 0.0],   # H
        [4.75, -2.10, 0.0],  # H
        [5.95, 0.0, 0.0],    # H
        [2.80, 0.0, 1.0],    # extra
    ], dtype=float)
    P_new, n_fixed, n_skip = align_donor_lonepairs(syms, P, 0, (1,))
    assert n_fixed == 0
    assert n_skip == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


def test_aromatic_face_on_rotated_when_active():
    """When M is on the ring-normal axis (face-on), rotate 90 deg around M-D.

    Pre-rotation: M->D direction is in z, ring lies in xy with normal in z =>
    face-on (theta = 0).
    Post-rotation: ring should now have its normal perpendicular to M-D
    (i.e. ring in xz plane), and the metric_md_direction face-on test should
    no longer fire.
    """
    os.environ["DELFIN_FFFREE_MDDIR_ALIGN"] = "1"
    import sys
    sys.path.insert(0, "/home/qmchem_max/agent_workspace/quality_framework/scripts")
    from delfin.fffree.mddir_alignment import align_donor_lonepairs
    from metric_md_direction import md_direction_violations  # type: ignore

    # Pd at +z, ring in xy with N at origin (donor).  Ring substituents are
    # the 2 ring carbons in xy.  M->D vector points -z; ring normal is +z =>
    # face-on (theta 0).
    syms = ["Pd", "N", "C", "C", "C", "C", "C"]
    P = np.array([
        [0.0, 0.0, 2.10],     # Pd (above ring)
        [0.0, 0.0, 0.0],      # N donor at ring vertex
        [1.20, 0.7, 0.0],     # C neighbour 1 in xy
        [-1.20, 0.7, 0.0],    # C neighbour 2 in xy
        [1.20, 2.10, 0.0],    # ring C
        [-1.20, 2.10, 0.0],   # ring C
        [0.0, 2.80, 0.0],     # ring C
    ], dtype=float)
    # confirm pre-rotation: face-on detected
    pre = md_direction_violations(syms, P)
    assert any(v["mode"] == "aromatic_face" for v in pre), (
        f"pre-rotation should flag aromatic_face, got {pre}"
    )
    P_new, n_fixed, _ = align_donor_lonepairs(syms, P, 0, (1,))
    assert n_fixed == 1
    # M-D distance EXACTLY preserved (axis rotation through donor)
    d_old = float(np.linalg.norm(P[1] - P[0]))
    d_new = float(np.linalg.norm(P_new[1] - P_new[0]))
    assert abs(d_old - d_new) < 1e-8
    # Post: face-on no longer triggers
    post = md_direction_violations(syms, P_new)
    assert not any(v["mode"] == "aromatic_face" for v in post), (
        f"post-rotation should NOT flag aromatic_face, got {post}"
    )


def test_md_invariant_preserved_aromatic():
    """Pure axis rotation preserves M-D distance exactly."""
    os.environ["DELFIN_FFFREE_MDDIR_ALIGN"] = "1"
    from delfin.fffree.mddir_alignment import align_donor_lonepairs

    syms = ["Pd", "N", "C", "C", "C", "C", "C"]
    P = np.array([
        [0.0, 0.0, 2.10],
        [0.0, 0.0, 0.0],
        [1.20, 0.7, 0.0],
        [-1.20, 0.7, 0.0],
        [1.20, 2.10, 0.0],
        [-1.20, 2.10, 0.0],
        [0.0, 2.80, 0.0],
    ], dtype=float)
    P_new, _, _ = align_donor_lonepairs(syms, P, 0, (1,))
    np.testing.assert_allclose(P_new[1], P[1], atol=1e-9)
    np.testing.assert_allclose(P_new[0], P[0], atol=1e-9)


def test_donor_linearized_skipped():
    """If max(M-D-X) > 168, the donor is sp not sp3 -> skip, count in n_skip."""
    os.environ["DELFIN_FFFREE_MDDIR_ALIGN"] = "1"
    from delfin.fffree.mddir_alignment import align_donor_lonepairs

    # M-O-C linear at ~179.9 deg
    syms = ["Cd", "O", "C", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Cd
        [2.40, 0.0, 0.0],    # O donor
        [3.60, 0.0, 0.0],    # C linear with M-O
        [4.80, 0.0, 0.0],    # next C
    ], dtype=float)
    _, n_fixed, n_skip = align_donor_lonepairs(syms, P, 0, (1,))
    assert n_fixed == 0
    assert n_skip == 1


def test_sp3_anti_180_rotation():
    """sp3 donor with M on substituent side -> 180 rotation around M-D."""
    os.environ["DELFIN_FFFREE_MDDIR_ALIGN"] = "1"
    import sys
    sys.path.insert(0, "/home/qmchem_max/agent_workspace/quality_framework/scripts")
    from delfin.fffree.mddir_alignment import align_donor_lonepairs
    from metric_md_direction import md_direction_violations  # type: ignore

    # Detector definition: v_md = D->M; v_sub = D->centroid;
    # phi = angle(v_md, v_sub).  FLAG when phi < 90 (M on substituent side).
    # Make v_md and v_sub aligned: M at +z, substituents at +z (both above D).
    # Use Co (transition metal) -> Co-N ideal ~2.20 A -> place M at +1.95 z.
    # Substituents at distance ~1.47 A from N (within geom bond detection).
    syms = ["Co", "N", "C", "C", "C"]
    P = np.array([
        [0.0, 0.0, 1.95],     # Co (above N, within M-D detection)
        [0.0, 0.0, 0.0],      # N donor
        [0.80, 0.46, 0.40],   # C (above N, same side as Co)
        [-0.80, 0.46, 0.40],  # C (above N)
        [0.0, -0.92, 0.40],   # C (above N)
    ], dtype=float)
    pre = md_direction_violations(syms, P)
    # sp3_anti should be the mode flagged (phi < 90)
    assert any(v["mode"] == "sp3_anti" for v in pre), pre
    P_new, n_fixed, _ = align_donor_lonepairs(syms, P, 0, (1,))
    assert n_fixed == 1
    # M-D invariant
    d_old = float(np.linalg.norm(P[1] - P[0]))
    d_new = float(np.linalg.norm(P_new[1] - P_new[0]))
    assert abs(d_old - d_new) < 1e-8
    # After 180 rotation around M-D = z axis through donor:
    # the substituents flip to z<0 (anti to M which is at +z)
    post = md_direction_violations(syms, P_new)
    assert not any(v["mode"] == "sp3_anti" for v in post)


def test_no_donor_element_no_op():
    """Donor list with only non-donor elements (e.g. C/H) -> no action."""
    os.environ["DELFIN_FFFREE_MDDIR_ALIGN"] = "1"
    from delfin.fffree.mddir_alignment import align_donor_lonepairs

    syms = ["Pd", "C", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 1.0, 0.0],
        [3.0, -1.0, 0.0],
    ], dtype=float)
    P_new, n_fixed, n_skip = align_donor_lonepairs(syms, P, 0, (1,))
    assert n_fixed == 0
    assert n_skip == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


def test_master_flag_enables():
    os.environ["DELFIN_FFFREE_CONSTRUCTION_FIX_ALL"] = "1"
    from delfin.fffree.mddir_alignment import _flag_active
    assert _flag_active() is True
