"""Tests for delfin.fffree.amide_vsepr_template — F24-plane mode (surgical fix 3).

Validate:
  - Default OFF: enforce_amide_planarity uses the SVD plane (byte-identical
    to pre-fix HEAD).
  - With DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE=1 (and AMIDE_VSEPR active):
    the projector uses ``cross(N->nbr0, N->nbr1)`` as the plane normal
    (F24-detector's plane) instead of the substituent-SVD normal.
  - For coplanar substituents the two plane definitions give the same
    result (since SVD on 3 points yields the unique 3-point plane normal).
  - Master flag enables both AMIDE_VSEPR + F24 plane.
"""
import os
import numpy as np
import pytest


def _reset_env():
    for k in (
        "DELFIN_FFFREE_AMIDE_VSEPR",
        "DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE",
        "DELFIN_FFFREE_F24_INTERLIG_FIX_ALL",
        "DELFIN_FFFREE_CONSTRUCTION_FIX_ALL",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


def _formamide_like(oop: float = 0.4):
    """Build a minimal amide-N motif.

    Atom layout (6 atoms):
      0 C_carbonyl   at (-0.5, 0, 0)         (sp2 C, partner of amide N)
      1 O_carbonyl   at (-1.2, 0.95, 0)      (C=O)
      2 R_carbonyl   at (-1.2, -0.95, 0)     (3rd C-substituent on carbonyl)
      3 N_amide      at (0.7, 0, oop)        (pyramidal: OOP from substituent plane)
      4 C_alpha1     at (1.4, 1.0, 0)        (sp3 C, sticking out of plane in
                                              real life but here flat for test)
      5 C_alpha2     at (1.4, -1.0, 0)       (sp3 C, also flat for test)
    """
    syms = ["C", "O", "C", "N", "C", "C"]
    P = np.array([
        [-0.5, 0.0, 0.0],
        [-1.2, 0.95, 0.0],
        [-1.2, -0.95, 0.0],
        [0.7, 0.0, oop],
        [1.4, 1.0, 0.0],
        [1.4, -1.0, 0.0],
    ], dtype=float)
    return syms, P


def test_f24_plane_default_off():
    """Without the F24 flag, _f24_plane_active is False."""
    from delfin.fffree.amide_vsepr_template import _f24_plane_active
    assert _f24_plane_active() is False


def test_f24_plane_env_on_per_fix():
    """Per-fix flag activates."""
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE"] = "1"
    from delfin.fffree.amide_vsepr_template import _f24_plane_active
    assert _f24_plane_active() is True


def test_f24_plane_env_on_master_f24():
    """F24-interlig master flag activates."""
    os.environ["DELFIN_FFFREE_F24_INTERLIG_FIX_ALL"] = "1"
    from delfin.fffree.amide_vsepr_template import _f24_plane_active
    assert _f24_plane_active() is True


def test_env_off_returns_pcopy_byte_identical():
    """With AMIDE_VSEPR OFF the projector returns (P_copy, 0) and the
    coords are byte-identical."""
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity
    syms, P = _formamide_like(oop=0.4)
    P_new, n = enforce_amide_planarity(syms, P, metal_idx=-1, donor_idxs=())
    assert n == 0
    np.testing.assert_array_equal(P_new, P)
    assert P_new.tobytes() == P.tobytes()


def test_f24_plane_byte_identical_to_svd_when_substituents_coplanar():
    """For 3 coplanar substituents AND N already lying in the substituent
    plane (no OOP at all), the cross-product plane normal and the SVD
    plane normal are parallel and projection is a no-op for both modes.
    """
    # N at (0.7, 0, 0) -- already coplanar with substituents on z=0.
    syms, P = _formamide_like(oop=0.0)
    # Run with SVD plane (default).
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity
    P_svd, n_svd = enforce_amide_planarity(
        syms, P, metal_idx=-1, donor_idxs=(),
    )
    # Run with F24 plane.
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE"] = "1"
    P_f24, n_f24 = enforce_amide_planarity(
        syms, P, metal_idx=-1, donor_idxs=(),
    )
    # Both find an amide-N but both are below the OOP threshold so both
    # return the input unchanged.
    assert n_svd == n_f24 == 0
    np.testing.assert_allclose(P_svd, P_f24, atol=1e-9)


def test_f24_plane_diverges_from_svd_when_N_pyramidal():
    """When N is pyramidal (OOP > 0) the two plane definitions point in
    different directions: SVD normal = substituent-plane-normal (z-axis
    here), F24 normal = cross(N->nbr0, N->nbr1) which depends on N.
    Projection results therefore differ -- this is the whole point of
    surgical fix 3 (template-vs-detector plane mismatch).
    """
    syms, P = _formamide_like(oop=0.4)
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity
    P_svd, n_svd = enforce_amide_planarity(
        syms, P, metal_idx=-1, donor_idxs=(),
    )
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE"] = "1"
    P_f24, n_f24 = enforce_amide_planarity(
        syms, P, metal_idx=-1, donor_idxs=(),
    )
    assert n_svd >= 1
    assert n_f24 >= 1
    # The two results must DIFFER -- SVD projects onto substituent plane,
    # F24 projects onto plane (substituent centroid, normal=cross(N->nbr0,
    # N->nbr1)).
    assert not np.allclose(P_svd, P_f24, atol=1e-6)


def test_f24_plane_matches_detector_definition():
    """After F24-mode projection, N lies on a plane whose normal direction
    matches the F24 detector's plane normal (computed at the PRE-projection
    N coordinates).

    The F24 detector's plane is ``plane(N, nbr0, nbr1)`` with normal
    ``cross(nbr0 - N, nbr1 - N)``.  When ``DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE
    =1`` is set the projector uses precisely that normal (rather than the
    SVD-fit of the 3 substituent points) so the projection direction
    matches the detector's plane direction.

    NOTE: an EXACT zero of the post-projection F24 metric is mathematically
    impossible because moving N changes the (N, nbr0, nbr1) plane.  The
    spec's contract is plane-DIRECTION agreement, which this test verifies.
    """
    syms, P = _formamide_like(oop=0.4)
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE"] = "1"
    from delfin.fffree.amide_vsepr_template import (
        enforce_amide_planarity, find_amide_nitrogens,
    )
    P_new, n = enforce_amide_planarity(
        syms, P, metal_idx=-1, donor_idxs=(),
    )
    assert n >= 1
    amides = find_amide_nitrogens(syms, P_new)
    assert amides
    ni, _ = amides[0]
    from delfin._bond_decollapse import _ideal_bond, _is_metal
    nbrs = []
    for j in range(len(syms)):
        if j == ni or syms[j] == "H" or _is_metal(syms[j]):
            continue
        d = float(np.linalg.norm(P_new[ni] - P_new[j]))
        if d < 1.30 * _ideal_bond(syms[ni], syms[j]):
            nbrs.append(j)
    assert len(nbrs) == 3
    # F24 detector's plane normal at PRE-projection coordinates.
    v1 = P[nbrs[0]] - P[ni]
    v2 = P[nbrs[1]] - P[ni]
    nrm = np.cross(v1, v2)
    nl = float(np.linalg.norm(nrm))
    assert nl > 1e-9
    nrm /= nl
    # The movement vector (N_new - N_old) should be PARALLEL to the F24
    # detector's plane-normal direction (the projection moves N along the
    # normal axis exactly).
    move = P_new[ni] - P[ni]
    move_n = float(np.linalg.norm(move))
    assert move_n > 1e-9
    move_unit = move / move_n
    # Parallel means |dot| ~ 1.
    cos_align = float(abs(np.dot(move_unit, nrm)))
    assert cos_align > 1.0 - 1e-9, (
        f"projection direction not parallel to F24 normal: cos={cos_align}"
    )


def test_master_flag_enables_amide_and_f24():
    """DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1 enables BOTH amide-vsepr +
    F24-plane mode."""
    os.environ["DELFIN_FFFREE_F24_INTERLIG_FIX_ALL"] = "1"
    from delfin.fffree.amide_vsepr_template import (
        _flag_active, _f24_plane_active,
    )
    assert _flag_active() is True
    assert _f24_plane_active() is True


def test_global_master_flag_enables():
    """DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1 also enables F24-plane mode."""
    os.environ["DELFIN_FFFREE_CONSTRUCTION_FIX_ALL"] = "1"
    from delfin.fffree.amide_vsepr_template import (
        _flag_active, _f24_plane_active,
    )
    assert _flag_active() is True
    assert _f24_plane_active() is True


def test_md_invariant_preserved_with_f24_plane():
    """With donors specified, M-D distances are preserved after projection."""
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR"] = "1"
    os.environ["DELFIN_FFFREE_AMIDE_VSEPR_F24_PLANE"] = "1"
    from delfin.fffree.amide_vsepr_template import enforce_amide_planarity
    syms, P = _formamide_like(oop=0.4)
    # Pretend the carbonyl C (idx 0) is a donor coordinated to a metal
    # we'll fake at index 99 (not in the array, so no effect).
    # Use a real test: M at idx 0, donor at idx 3 (the amide-N itself).
    syms2 = ["Au"] + syms
    P2 = np.vstack([np.array([[-3.0, 0.0, 0.0]]), P])
    # Amide-N is at original idx 3 -> new idx 4.  Make idx 1 (was C
    # carbonyl) the donor to the metal at idx 0.
    P_new, n = enforce_amide_planarity(
        syms2, P2, metal_idx=0, donor_idxs=[1],
    )
    md_before = float(np.linalg.norm(P2[1] - P2[0]))
    md_after = float(np.linalg.norm(P_new[1] - P_new[0]))
    assert abs(md_before - md_after) < 1e-9
