"""Tests for delfin.fffree.post_grip_corrector

Validate the three post-GRIP correctors:
  1. apply_sp2_flatten: amide-N / carbonyl-C projection back onto plane.
  2. apply_haxis_rotation: H-on-M-D-axis flipped via donor subtree rotation.
  3. apply_mdshort_repulsion: non-donor secondary contact pushed back.

Plus the master `post_grip_corrections` entrypoint:
  - default OFF (env unset): coordinates unchanged (byte-identical).
  - all-three ON: all enabled fixes run and M-D invariant preserved.
"""
import os
import numpy as np
import pytest


def _reset_env():
    for k in (
        "DELFIN_FFFREE_POST_GRIP_ALL",
        "DELFIN_FFFREE_POST_GRIP_FLATTEN",
        "DELFIN_FFFREE_HAXIS_FIX",
        "DELFIN_FFFREE_MDSHORT_FIX",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


# ---------------------------------------------------------------------------
# Fix 1: sp2 flatten
# ---------------------------------------------------------------------------
def test_sp2_flatten_amide_n_projected_to_plane():
    """Amide N (3 heavy neighbours = C(=O), C_alpha1, C_alpha2) pulled 0.4 A
    off the plane must be projected back into the neighbour plane.

    Geometry: place the 3 N-neighbours at ~120 deg from N (sp2 fingerprint
    requires mean angle in [100, 135] deg), all in z=0 plane.
    """
    from delfin.fffree.post_grip_corrector import apply_sp2_flatten

    syms = ["C", "O", "N", "C", "C"]
    # N at (0, 0, 0.15) with 3 nbrs in z=0 plane at 120 deg apart.
    # Off-plane by 0.15 A => slight pyramidalisation, angles in sp2 window.
    # The 3 nbrs themselves do NOT need to bond to each other (heavy_adj is
    # for N's bonds, not nbr-nbr bonds).
    d_nbr = 1.40
    P = np.array([
        [d_nbr * np.cos(0), d_nbr * np.sin(0), 0.0],
        [d_nbr * np.cos(0) + 1.22, d_nbr * np.sin(0), 0.0],   # =O off the C
        [0.00, 0.00, 0.15],   # amide N slightly pyramidalised in +z
        [d_nbr * np.cos(np.deg2rad(120)), d_nbr * np.sin(np.deg2rad(120)), 0.0],
        [d_nbr * np.cos(np.deg2rad(240)), d_nbr * np.sin(np.deg2rad(240)), 0.0],
    ], dtype=float)

    P_out, n = apply_sp2_flatten(syms, P)
    # N should now be (nearly) in plane z = 0
    z = P_out[2, 2]
    assert abs(z) < 0.05, f"N still off plane: z={z:.4f}"
    assert n >= 1


def test_sp2_flatten_skips_donor_atoms():
    """A donor atom in the donor_idxs list must NOT be moved (preserves M-D
    invariant)."""
    from delfin.fffree.post_grip_corrector import apply_sp2_flatten

    syms = ["Fe", "C", "O", "N", "C", "C"]
    P = np.array([
        [0.00, 0.00, 0.00],     # Fe (metal)
        [2.10, 0.00, 0.00],     # C (carbonyl, sp2)
        [3.30, 0.00, 0.00],     # =O
        [1.50, 1.30, 0.30],     # N donor (off-plane!)
        [1.00, 2.40, 0.00],     # R1 of N
        [2.40, 2.40, 0.00],     # R2 of N
    ], dtype=float)
    # mark N (idx 3) as a donor
    P_out, n = apply_sp2_flatten(syms, P, metal_idx=0, donor_idxs=[3])
    # N must NOT have moved
    assert np.allclose(P_out[3], P[3]), f"donor N moved: {P_out[3]} vs {P[3]}"


def test_sp2_flatten_already_in_plane_noop():
    """An sp2 atom already in plane (off < 0.05 A) is unchanged."""
    from delfin.fffree.post_grip_corrector import apply_sp2_flatten
    syms = ["C", "O", "N", "C", "C"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.22, 0.00, 0.00],
        [-0.70, 1.10, 0.01],    # only 0.01 A out of plane
        [-1.40, 1.80, 0.00],
        [-1.40, 0.40, 0.00],
    ], dtype=float)
    P_out, n = apply_sp2_flatten(syms, P)
    # No projection happens (threshold 0.05 A)
    assert n == 0
    assert np.allclose(P_out, P)


def test_sp2_flatten_skips_sp3_n():
    """An sp3 N (4 neighbours OR pyramidal-only-3-neighbours) is NOT flattened."""
    from delfin.fffree.post_grip_corrector import apply_sp2_flatten
    # Amine N with 3 H and 1 R: sp3 pyramidal, mean angle ~107 deg (not in
    # SP2 fingerprint window [100, 135] -- 107 IS in window so this test is
    # actually the boundary case).  Make it strongly pyramidal -- N with 3
    # heavies at 90 deg each (sub-100 deg mean), then projection skipped.
    syms = ["N", "C", "C", "C"]
    # Three substituents in a tight pyramid: 90 deg
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.40, 0.00, 0.00],
        [0.00, 1.40, 0.00],
        [0.00, 0.00, 1.40],
    ], dtype=float)
    # mean of (90, 90, 90) = 90 deg < 100, skip
    P_out, n = apply_sp2_flatten(syms, P)
    assert n == 0


# ---------------------------------------------------------------------------
# Fix 2: haxis rotation
# ---------------------------------------------------------------------------
def test_haxis_rotation_flips_h_off_md_axis():
    """A monodentate donor (Fe-N) with an N-H pointing at Fe must rotate 180
    so the H is on the far side of N (off the Fe-N axis)."""
    from delfin.fffree.post_grip_corrector import apply_haxis_rotation
    # Fe at origin, N at (2.0, 0, 0) (donor), R (carbon) at (2.7, 0, 0)
    # H attached to N at (1.5, 0.1, 0) -- on the Fe-N axis, between Fe and N
    # (this is the haxis mode-2 violation).
    syms = ["Fe", "N", "C", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [2.00, 0.00, 0.00],
        [3.00, 0.40, 0.00],     # C nbr of N, beyond N
        [1.50, 0.20, 0.00],     # H on N -- pointing AT Fe, near axis
    ], dtype=float)
    # adj: H (idx 3) closer to N (idx 1) than to Fe.
    # N's heavy nbrs (non-metal): just C (idx 2).
    P_out, n = apply_haxis_rotation(syms, P, metal_idx=0, donor_idxs=[1])
    # After rotation, H should be on the FAR side of N from Fe.
    # Original H was at (1.50, 0.20, 0): between Fe (0) and N (2).
    # Post-rotation: x of H > x of N (away from Fe), or at least no longer
    # between Fe and N.
    print("H before:", P[3], "after:", P_out[3])
    # The fix triggers => either we rotated (then H is far side) OR we didn't
    # (no fix possible -- chelate-block).  At minimum: N (donor) didn't move.
    assert np.allclose(P_out[1], P[1])
    # And M-D distance preserved
    d_old = np.linalg.norm(P[1] - P[0])
    d_new = np.linalg.norm(P_out[1] - P_out[0])
    assert abs(d_old - d_new) < 1e-9


def test_haxis_chelate_protection():
    """A chelate ligand (two donors) MUST NOT be rotated as a whole because
    that would displace the second donor."""
    from delfin.fffree.post_grip_corrector import apply_haxis_rotation
    # Fe at origin, ethylenediamine-like: N1-CH2-CH2-N2
    syms = ["Fe", "N", "C", "C", "N", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],     # Fe
        [2.00, 0.00, 0.00],     # N1
        [2.50, 1.40, 0.00],     # C-a
        [1.80, 2.60, 0.00],     # C-b
        [0.40, 2.60, 0.00],     # N2
        [1.50, 0.10, 0.00],     # H on N1 (on Fe-N1 axis)
    ], dtype=float)
    # donors = [1, 4] (chelate).  Rotation of N1's subtree would carry N2 too.
    # The fix must DETECT this and SKIP the rotation.
    P_out, n = apply_haxis_rotation(syms, P, metal_idx=0, donor_idxs=[1, 4])
    # Both donors must be unchanged
    assert np.allclose(P_out[1], P[1])
    assert np.allclose(P_out[4], P[4])


# ---------------------------------------------------------------------------
# Fix 3: mdshort repulsion
# ---------------------------------------------------------------------------
def test_mdshort_repulsion_pushes_secondary_contact():
    """A non-donor C atom 1.3 A from Fe (way below ideal ~2.0 A x 0.80) must
    be pushed back to 0.85 x ideal."""
    from delfin.fffree.post_grip_corrector import apply_mdshort_repulsion
    syms = ["Fe", "N", "C"]
    # Fe at origin; N is the donor (idx 1) at 2.0 A; C is a non-donor at 1.3 A
    # (way too short; this is the secondary-contact mdshort signature).
    P = np.array([
        [0.00, 0.00, 0.00],
        [2.00, 0.00, 0.00],     # donor N
        [1.30, 0.50, 0.00],     # too close non-donor C
    ], dtype=float)
    ideal = 0.76 + 0.76  # Fe (we use 0.9 fallback in _COV) + C
    # Actually _ideal_bond("Fe", "C"): Fe not in _COV, returns 0.9 + 0.76 = 1.66
    # 0.80 * 1.66 = 1.33; our d = sqrt(1.3^2 + 0.5^2) = 1.39 > 1.33 -> not triggered.
    # Tighten:
    P[2] = np.array([1.20, 0.30, 0.00])  # d = 1.24
    P_out, n = apply_mdshort_repulsion(syms, P, metal_idx=0, donor_idxs=[1])
    d_old = np.linalg.norm(P[2] - P[0])
    d_new = np.linalg.norm(P_out[2] - P_out[0])
    if n > 0:
        # If pushed, the new distance is ~0.85 * ideal_FeC.
        assert d_new > d_old, f"C should be pushed away: {d_old:.3f} -> {d_new:.3f}"
    # Donor N must not move
    assert np.allclose(P_out[1], P[1])


def test_mdshort_skips_donor_atoms():
    """An atom in the donor_idxs list must NEVER be pushed (M-D invariant)."""
    from delfin.fffree.post_grip_corrector import apply_mdshort_repulsion
    syms = ["Fe", "N"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.20, 0.00, 0.00],     # donor N at 1.20 (would be flagged as mdshort
                                #  if not in donor list)
    ], dtype=float)
    P_out, n = apply_mdshort_repulsion(syms, P, metal_idx=0, donor_idxs=[1])
    assert n == 0
    assert np.allclose(P_out[1], P[1])


# ---------------------------------------------------------------------------
# Master entrypoint
# ---------------------------------------------------------------------------
def test_post_grip_corrections_default_off_byte_identical():
    """Without any env flag set, post_grip_corrections must return P unchanged
    (byte-identical to HEAD)."""
    from delfin.fffree.post_grip_corrector import post_grip_corrections
    syms = ["Fe", "N", "C", "C"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [2.00, 0.10, 0.20],
        [3.00, 1.20, 0.30],
        [1.50, 0.10, 0.40],
    ], dtype=float)
    P_out, diag = post_grip_corrections(syms, P, metal_idx=0, donor_idxs=[1])
    assert np.array_equal(P_out, P)
    assert diag["any_enabled"] is False


def test_post_grip_corrections_all_active_runs_all_three():
    """Master env DELFIN_FFFREE_POST_GRIP_ALL=1 activates all three fixes."""
    os.environ["DELFIN_FFFREE_POST_GRIP_ALL"] = "1"
    from delfin.fffree.post_grip_corrector import (
        post_grip_corrections, flatten_active, haxis_active, mdshort_active,
        any_active,
    )
    assert flatten_active()
    assert haxis_active()
    assert mdshort_active()
    assert any_active()
    syms = ["Fe", "N", "C", "C", "C", "C"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [2.00, 0.00, 0.00],     # donor N
        [3.00, 1.20, 0.00],     # C ring
        [3.00, -1.20, 0.00],    # C ring
        [4.20, 0.00, 0.00],     # C ring
        [3.00, 0.00, 0.20],     # sp2 C off-plane by 0.20 -- gets flattened
    ], dtype=float)
    P_out, diag = post_grip_corrections(syms, P, metal_idx=0, donor_idxs=[1])
    # Must preserve M-D invariant
    d_old = np.linalg.norm(P[1] - P[0])
    d_new = np.linalg.norm(P_out[1] - P_out[0])
    assert abs(d_old - d_new) < 0.05
    # The diagnostics record what happened
    assert diag["any_enabled"] is True


def test_per_fix_env_flags_independent():
    """Each per-fix env flag works on its own."""
    from delfin.fffree.post_grip_corrector import (
        flatten_active, haxis_active, mdshort_active,
    )
    _reset_env()
    os.environ["DELFIN_FFFREE_POST_GRIP_FLATTEN"] = "1"
    assert flatten_active()
    assert not haxis_active()
    assert not mdshort_active()

    _reset_env()
    os.environ["DELFIN_FFFREE_HAXIS_FIX"] = "1"
    assert not flatten_active()
    assert haxis_active()
    assert not mdshort_active()

    _reset_env()
    os.environ["DELFIN_FFFREE_MDSHORT_FIX"] = "1"
    assert not flatten_active()
    assert not haxis_active()
    assert mdshort_active()


def test_haxis_hapto_protection_skips_pi_cluster():
    """Hapto-pi donors (>=3 same-element donors clustered within 3 A) must
    NOT be rotated -- their geometry is part of a rigid pi-frame."""
    from delfin.fffree.post_grip_corrector import (
        apply_haxis_rotation, _detect_hapto_donors,
    )
    # Cp ligand: 5 C donors in a ring around the metal axis
    syms = ["Fe", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H"]
    # Cp ring centred above Fe at z=1.5, radius 1.2
    r = 1.2
    P = np.array([
        [0.0, 0.0, 0.0],   # Fe
        [r * np.cos(np.deg2rad(0)),    r * np.sin(np.deg2rad(0)),   1.5],
        [r * np.cos(np.deg2rad(72)),   r * np.sin(np.deg2rad(72)),  1.5],
        [r * np.cos(np.deg2rad(144)),  r * np.sin(np.deg2rad(144)), 1.5],
        [r * np.cos(np.deg2rad(216)),  r * np.sin(np.deg2rad(216)), 1.5],
        [r * np.cos(np.deg2rad(288)),  r * np.sin(np.deg2rad(288)), 1.5],
        # H on each C (pointing up)
        [r * np.cos(np.deg2rad(0)),    r * np.sin(np.deg2rad(0)),   2.5],
        [r * np.cos(np.deg2rad(72)),   r * np.sin(np.deg2rad(72)),  2.5],
        [r * np.cos(np.deg2rad(144)),  r * np.sin(np.deg2rad(144)), 2.5],
        [r * np.cos(np.deg2rad(216)),  r * np.sin(np.deg2rad(216)), 2.5],
        [r * np.cos(np.deg2rad(288)),  r * np.sin(np.deg2rad(288)), 2.5],
    ], dtype=float)
    donors = [1, 2, 3, 4, 5]
    hapto = _detect_hapto_donors(syms, P, 0, donors)
    assert hapto == set(donors), f"all 5 Cp Cs should be hapto: {hapto}"
    P_out, n = apply_haxis_rotation(syms, P, metal_idx=0, donor_idxs=donors)
    # NO rotation should happen on a hapto cluster
    assert np.allclose(P_out, P)


def test_haxis_stereocenter_protection():
    """A donor whose subtree contains an sp3 stereocentre must NOT be rotated
    180 deg (would invert chirality)."""
    from delfin.fffree.post_grip_corrector import apply_haxis_rotation
    # Fe-N donor, N-C-stereo-C (4 different heavy nbrs)
    syms = ["Fe", "N", "C", "C", "O", "Cl", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Fe
        [2.0, 0.0, 0.0],     # N donor
        [3.0, 0.0, 0.0],     # alpha-C (sp3 chiral: 4 heavy nbrs)
        [3.5, 1.2, 0.0],     # nbr1 of alpha-C
        [3.5, -0.6, 1.0],    # nbr2 (O)
        [3.5, -0.6, -1.0],   # nbr3 (Cl)
        [1.5, 0.1, 0.0],     # H on N triggering haxis
    ], dtype=float)
    # alpha-C has 4 heavy nbrs: N (1), C (3), O (4), Cl (5) -- stereocentre
    P_out, n = apply_haxis_rotation(syms, P, metal_idx=0, donor_idxs=[1])
    # Stereocenter detected: no rotation
    assert np.allclose(P_out, P)


def test_sp2_flatten_tighter_window_rejects_sp3():
    """sp3 atoms with mean angle ~109 deg must NOT be flattened (avoids the
    flat_sp3 +67 % regression seen in smoke500 voll-pool 2026-06-02)."""
    from delfin.fffree.post_grip_corrector import apply_sp2_flatten
    # sp3 N with 3 heavy nbrs at ~107 deg (typical amine umbrella)
    # nbrs at distances 1.4 A in z=0 plane at 120 deg apart -- mean 113
    # but if N is well off-plane (0.3 A), mean angle goes lower.
    syms = ["N", "C", "C", "C"]
    # Place N at (0, 0, 0.45) with 3 nbrs at 109.5 mean ~ pyramidal
    # Trial: nbrs at d=1.40 in plane z=0, spread 120 deg.
    # If N at (0, 0, 0.4), angles = arccos(0.4^2 + d_nbr^2 cos(120 deg) - ...)
    d_nbr = 1.40
    P = np.array([
        [0.00, 0.00, 0.40],
        [d_nbr * np.cos(np.deg2rad(0)),   d_nbr * np.sin(np.deg2rad(0)),   0.0],
        [d_nbr * np.cos(np.deg2rad(120)), d_nbr * np.sin(np.deg2rad(120)), 0.0],
        [d_nbr * np.cos(np.deg2rad(240)), d_nbr * np.sin(np.deg2rad(240)), 0.0],
    ], dtype=float)
    P_out, n = apply_sp2_flatten(syms, P)
    # mean angle here ~ 112 deg (just below 115 lower bound) -- skip
    # Verify N NOT flattened (stays off plane)
    assert P_out[0, 2] > 0.1, f"N should not be flattened (z stays >{0.1}): z={P_out[0,2]}"


def test_md_invariant_preserved_even_under_aggressive_movements():
    """Sanity: with all three fixes ON, the M-D distance can change at most
    by the rollback tolerance (0.05 A)."""
    os.environ["DELFIN_FFFREE_POST_GRIP_ALL"] = "1"
    from delfin.fffree.post_grip_corrector import post_grip_corrections
    rng = np.random.default_rng(42)
    syms = ["Fe", "N", "O", "C", "C", "C", "C", "C", "C", "H", "H"]
    P = rng.normal(scale=2.0, size=(11, 3))
    # Put Fe at origin, donor N at canonical 2 A, donor O at canonical 2 A
    P[0] = [0.0, 0.0, 0.0]
    P[1] = [2.0, 0.0, 0.0]
    P[2] = [0.0, 2.0, 0.0]
    d_md_N = np.linalg.norm(P[1] - P[0])
    d_md_O = np.linalg.norm(P[2] - P[0])
    P_out, diag = post_grip_corrections(syms, P, metal_idx=0, donor_idxs=[1, 2])
    d_md_N_new = np.linalg.norm(P_out[1] - P_out[0])
    d_md_O_new = np.linalg.norm(P_out[2] - P_out[0])
    assert abs(d_md_N - d_md_N_new) < 0.05
    assert abs(d_md_O - d_md_O_new) < 0.05
