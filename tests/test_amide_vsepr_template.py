"""Tests for delfin.fffree.amide_vsepr_template

Validate the build-time amide planarity template:
  - Default OFF byte-identical (no env flag set).
  - Graph-based amide motif detection (C(=O)-N with 3 N substituents).
  - SVD projection of pyramidal amide-N onto neighbour plane.
  - Rigid-H drag for H attached to the moved N.
  - M-D invariant: donor / metal never moved; if N IS a donor, skipped.
  - Symmetric env flags (per-fix + master).
"""
import os
import numpy as np
import pytest


def _reset_env():
    for k in (
        "DELFIN_FFFREE_AMIDE_VSEPR",
        "DELFIN_FFFREE_CONSTRUCTION_FIX_ALL",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


def test_default_off_byte_identical():
    """No env flag set: coordinates unchanged, n_fixed == 0."""
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity

    # synthetic amide: O=C-N(C)(C) pyramidal N (0.3 A off plane)
    syms = ["O", "C", "N", "C", "C"]
    P = np.array([
        [0.0, 1.21, 0.0],    # O (C=O at 1.21 A)
        [0.0, 0.0, 0.0],     # carbonyl C
        [1.30, -0.75, 0.0],  # amide N (in plane = below carbonyl)
        [2.20, -0.10, 0.30], # alpha-C1 (above plane)
        [1.30, -2.05, 0.0],  # alpha-C2 (in plane)
    ], dtype=float)
    P_new, n_fixed = enforce_amide_planarity(syms, P)
    assert n_fixed == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


def test_amide_n_projected_to_plane_when_active():
    """With env flag on, an amide-N >0.10 A off the substituent plane is
    projected back into that plane.  Carbonyl C at (0,0,0) with =O at
    (0,1.21,0); N at (1.30,-0.75,0.40) (pyramidal); 2 N substituents in z=0.
    """
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity

    syms = ["O", "C", "N", "C", "C"]
    P = np.array([
        [0.0, 1.21, 0.0],    # O
        [0.0, 0.0, 0.0],     # carbonyl C  (N nbr 1)
        [1.30, -0.75, 0.40], # amide N (off plane by ~0.40 A)
        [2.40, -1.20, 0.0],  # alpha-C1 (N nbr 2)
        [1.30, -2.10, 0.0],  # alpha-C2 (N nbr 3)
    ], dtype=float)
    P_new, n_fixed = enforce_amide_planarity(syms, P)
    assert n_fixed == 1
    # Compute the plane of the 3 N-substituents and verify N now lies on it
    nbr_idxs = [1, 3, 4]  # carbonyl-C + 2 alpha-C
    centroid = np.mean([P_new[k] for k in nbr_idxs], axis=0)
    Q = np.array([P_new[k] - centroid for k in nbr_idxs])
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    normal = Vt[-1] / np.linalg.norm(Vt[-1])
    off = abs(float(np.dot(P_new[2] - centroid, normal)))
    assert off < 1e-6, f"amide-N still off plane: {off}"


def test_master_flag_enables():
    """Master flag DELFIN_FFFREE_CONSTRUCTION_FIX_ALL enables the fix too."""
    os.environ["DELFIN_FFFREE_CONSTRUCTION_FIX_ALL"] = "1"
    from delfin.fffree.amide_vsepr_template import (
        enforce_amide_planarity, _flag_active,
    )
    assert _flag_active() is True
    syms = ["O", "C", "N", "C", "C"]
    P = np.array([
        [0.0, 1.21, 0.0],
        [0.0, 0.0, 0.0],
        [1.30, -0.75, 0.40],
        [2.40, -1.20, 0.0],
        [1.30, -2.10, 0.0],
    ], dtype=float)
    _, n_fixed = enforce_amide_planarity(syms, P)
    assert n_fixed == 1


def test_donor_amide_n_skipped():
    """If the amide-N is in donor_idxs (rare amidate coordination), skip it.

    M-D invariant must be exactly preserved.
    """
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity

    # M (Pd at origin) + amide-N donor + carbonyl-C + 2 alpha-C + O
    syms = ["Pd", "N", "C", "O", "C", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Pd (M)
        [2.10, 0.0, 0.30],   # amide-N (DONOR; off plane by 0.30 A)
        [3.30, 0.75, 0.0],   # carbonyl C
        [3.30, 1.96, 0.0],   # O
        [3.30, -0.10, 0.0],  # alpha-C1 (N nbr extension)
        [2.10, -1.20, 0.0],  # alpha-C2
    ], dtype=float)
    P_new, n_fixed = enforce_amide_planarity(
        syms, P, metal_idx=0, donor_idxs=(1,)
    )
    # N at idx 1 should NOT move (donor)
    np.testing.assert_allclose(P_new[1], P[1], atol=1e-12)


def test_md_invariant_preserved():
    """When the amide-N is NOT a donor but the metal+donors are present,
    M-D distances must be preserved exactly.
    """
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity

    # Pd at origin + Cl donor at +x + amide ligand (O=C-N-CC) elsewhere
    syms = ["Pd", "Cl", "O", "C", "N", "C", "C"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Pd
        [2.30, 0.0, 0.0],    # Cl donor
        [-2.0, 1.21, 0.0],   # O
        [-2.0, 0.0, 0.0],    # carbonyl C
        [-0.70, -0.75, 0.40],# amide N off plane (NOT a donor here)
        [0.40, -1.20, 0.0],  # alpha-C1
        [-0.70, -2.05, 0.0], # alpha-C2
    ], dtype=float)
    d_md_old = float(np.linalg.norm(P[1] - P[0]))
    P_new, n_fixed = enforce_amide_planarity(
        syms, P, metal_idx=0, donor_idxs=(1,)
    )
    d_md_new = float(np.linalg.norm(P_new[1] - P_new[0]))
    assert abs(d_md_new - d_md_old) < 1e-8
    # N actually moved (and was below the plane)
    assert n_fixed == 1
    assert not np.allclose(P_new[4], P[4])


def test_no_amide_no_change():
    """A molecule without any amide motif must produce 0 fixes."""
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity

    # Ethylamine C-C-N-H2 — no carbonyl C neighbour to N, NOT an amide
    syms = ["C", "C", "N", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.54, 0.0, 0.0],
        [2.30, 1.20, 0.50],
        [3.30, 1.20, 0.50],
        [2.30, 1.20, 1.50],
    ], dtype=float)
    P_new, n_fixed = enforce_amide_planarity(syms, P)
    assert n_fixed == 0
    np.testing.assert_allclose(P_new, P, atol=1e-12)


def test_find_amide_nitrogens_motif():
    """Pure detection: find_amide_nitrogens returns (N_idx, C_idx) pairs."""
    from delfin.fffree.amide_vsepr_template import find_amide_nitrogens

    # 2 amides + 1 amine in one molecule
    syms = ["O", "C", "N", "C", "C", "O", "C", "N", "C", "C", "N", "C", "C", "C"]
    # Amide 1 (idx 2): O-C-N-C-C
    # Amide 2 (idx 7): O-C-N-C-C
    # Amine (idx 10): N-C-C-C (no carbonyl)
    P = np.array([
        [0.0, 1.21, 0.0],   # 0 O
        [0.0, 0.0, 0.0],    # 1 carbonyl-C
        [1.30, -0.75, 0.0], # 2 amide-N
        [2.20, -0.10, 0.0], # 3 alpha-C1
        [1.30, -2.05, 0.0], # 4 alpha-C2
        [5.0, 1.21, 0.0],   # 5 O
        [5.0, 0.0, 0.0],    # 6 carbonyl-C
        [6.30, -0.75, 0.0], # 7 amide-N
        [7.20, -0.10, 0.0], # 8 alpha-C1
        [6.30, -2.05, 0.0], # 9 alpha-C2
        [10.0, 0.0, 0.0],   # 10 amine N
        [11.30, 0.0, 0.0],  # 11
        [10.0, 1.30, 0.0],  # 12
        [10.0, -1.30, 0.0], # 13
    ], dtype=float)
    amides = find_amide_nitrogens(syms, P)
    # both N at idx 2 and idx 7 detected, amine at idx 10 NOT
    assert (2, 1) in amides
    assert (7, 6) in amides
    assert all(n_idx != 10 for n_idx, _ in amides)
