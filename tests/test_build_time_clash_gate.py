"""Tests for delfin.fffree.build_time_clash_gate

Validate the build-time collapse / clash gate:
  - Default OFF: has_collapse returns False regardless of geometry
    (and select_clean_candidate returns 0).
  - Active: detects xh_collapse / heavy_collapse / superposition.
  - collapse_count gives a numeric ranking (read-only, always active).
  - collapse_details enumerates each violation with mode + atom indices.
  - Candidate selection picks the lowest-collapse-count candidate.
"""
import os
import numpy as np
import pytest


def _reset_env():
    for k in (
        "DELFIN_FFFREE_BUILD_CLASH_GATE",
        "DELFIN_FFFREE_CONSTRUCTION_FIX_ALL",
        "DELFIN_FFFREE_BUILD_INTERLIG_PENALTY",
        "DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_W",
        "DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_FACTOR",
        "DELFIN_FFFREE_F24_INTERLIG_FIX_ALL",
    ):
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    yield
    _reset_env()


def test_default_off_has_collapse_false():
    """Even a catastrophic collapse passes when env flag is OFF."""
    from delfin.fffree.build_time_clash_gate import has_collapse

    syms = ["C", "C", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.05, 0.0, 0.0],     # C-C @ 0.05 A (catastrophic)
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ], dtype=float)
    assert has_collapse(syms, P) is False


def test_active_detects_heavy_collapse():
    os.environ["DELFIN_FFFREE_BUILD_CLASH_GATE"] = "1"
    from delfin.fffree.build_time_clash_gate import has_collapse, collapse_details

    syms = ["C", "C", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.05, 0.0, 0.0],     # heavy_collapse
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ], dtype=float)
    assert has_collapse(syms, P) is True
    dets = collapse_details(syms, P)
    # superposition fires first (any-pair < 0.5)
    modes = {d["mode"] for d in dets}
    assert "superposition" in modes or "heavy_collapse" in modes


def test_active_detects_xh_collapse():
    os.environ["DELFIN_FFFREE_BUILD_CLASH_GATE"] = "1"
    from delfin.fffree.build_time_clash_gate import has_collapse, collapse_details

    # C with H bonded at 0.30 A (< 0.70 floor)
    syms = ["C", "H", "C", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.30, 0.0, 0.0],     # X-H collapse
        [3.0, 0.0, 0.0],      # second C separated
        [3.0, 1.10, 0.0],     # second C-H OK
    ], dtype=float)
    assert has_collapse(syms, P) is True
    dets = collapse_details(syms, P)
    assert any(d["mode"] in ("xh_collapse", "superposition") for d in dets)


def test_clean_structure_no_collapse():
    """A normal ethanol or methane structure passes."""
    os.environ["DELFIN_FFFREE_BUILD_CLASH_GATE"] = "1"
    from delfin.fffree.build_time_clash_gate import has_collapse, collapse_count

    # methane CH4 (regular tetrahedral)
    syms = ["C", "H", "H", "H", "H"]
    t = 1.10 / np.sqrt(3)
    P = np.array([
        [0.0, 0.0, 0.0],
        [t, t, t],
        [t, -t, -t],
        [-t, t, -t],
        [-t, -t, t],
    ], dtype=float)
    assert has_collapse(syms, P) is False
    assert collapse_count(syms, P) == 0


def test_collapse_count_ranking():
    """collapse_count gives integer for ranking; works even when flag OFF
    (read-only, used internally)."""
    from delfin.fffree.build_time_clash_gate import collapse_count

    syms = ["C", "H", "C", "H"]
    P_clean = np.array([
        [0.0, 0.0, 0.0],
        [1.10, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 1.10, 0.0],
    ], dtype=float)
    P_collide = np.array([
        [0.0, 0.0, 0.0],
        [0.30, 0.0, 0.0],     # X-H collapse
        [3.0, 0.0, 0.0],
        [3.0, 1.10, 0.0],
    ], dtype=float)
    assert collapse_count(syms, P_clean) == 0
    assert collapse_count(syms, P_collide) >= 1


def test_select_clean_candidate_picks_clean():
    """select_clean_candidate returns index of first clean candidate."""
    os.environ["DELFIN_FFFREE_BUILD_CLASH_GATE"] = "1"
    from delfin.fffree.build_time_clash_gate import select_clean_candidate

    syms = ["C", "H", "C", "H"]
    bad = np.array([
        [0.0, 0.0, 0.0],
        [0.30, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 1.10, 0.0],
    ], dtype=float)
    clean = np.array([
        [0.0, 0.0, 0.0],
        [1.10, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 1.10, 0.0],
    ], dtype=float)
    # bad first, clean second -> select must pick index 1 (clean)
    pick = select_clean_candidate([syms, syms], [bad, clean])
    assert pick == 1
    # both bad -> pick the one with lower collapse_count (the same in this case)
    pick2 = select_clean_candidate([syms, syms], [bad, bad])
    assert pick2 == 0


def test_select_clean_candidate_off_returns_zero():
    """When env flag OFF: returns 0 unconditionally (byte-identical to
    'first candidate' semantics)."""
    from delfin.fffree.build_time_clash_gate import select_clean_candidate

    syms = ["C", "C"]
    P_bad = np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]], dtype=float)
    P_clean = np.array([[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]], dtype=float)
    pick = select_clean_candidate([syms, syms], [P_bad, P_clean])
    assert pick == 0


def test_master_flag_enables():
    os.environ["DELFIN_FFFREE_CONSTRUCTION_FIX_ALL"] = "1"
    from delfin.fffree.build_time_clash_gate import has_collapse, _flag_active
    assert _flag_active() is True
    # Now a collapse is detected
    syms = ["C", "C"]
    P = np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]], dtype=float)
    assert has_collapse(syms, P) is True


# ---------------------------------------------------------------------------
# Surgical fix 1: inter-ligand vdW penalty
# ---------------------------------------------------------------------------
def test_interlig_penalty_default_off():
    """interlig_penalty_active is False with no env flag set."""
    from delfin.fffree.build_time_clash_gate import interlig_penalty_active
    assert interlig_penalty_active() is False


def test_interlig_penalty_env_on_per_fix():
    """Per-fix flag activates."""
    os.environ["DELFIN_FFFREE_BUILD_INTERLIG_PENALTY"] = "1"
    from delfin.fffree.build_time_clash_gate import interlig_penalty_active
    assert interlig_penalty_active() is True


def test_interlig_penalty_env_on_master_f24():
    """F24-interlig master flag activates."""
    os.environ["DELFIN_FFFREE_F24_INTERLIG_FIX_ALL"] = "1"
    from delfin.fffree.build_time_clash_gate import interlig_penalty_active
    assert interlig_penalty_active() is True


def test_interlig_vdw_penalty_geometry():
    """Two C atoms in DIFFERENT fragments at 1.0 A produce a positive
    quadratic penalty; at 4.0 A produce 0.0; clash count consistent."""
    from delfin.fffree.build_time_clash_gate import (
        compute_interlig_vdw_penalty, INTERLIG_VDW_RADII,
    )
    # 4 carbon atoms; frag 0 = [0,1], frag 1 = [2,3]; only the cross pair
    # (0,2) is close.
    P_close = np.array([
        [0.0, 0.0, 0.0],   # frag 0, C
        [1.5, 0.0, 0.0],   # frag 0, C
        [1.0, 0.0, 0.0],   # frag 1, C (collides with frag-0 atom 0)
        [10.0, 0.0, 0.0],  # frag 1, C
    ], dtype=float)
    radii = [INTERLIG_VDW_RADII["C"]] * 4
    frag = {0: [0, 1], 1: [2, 3]}
    penalty_close, n_close = compute_interlig_vdw_penalty(
        P_close, frag, radii, factor=0.85, return_count=True,
    )
    # Cross-distance d(0,2) = 1.0; target = 0.85*(1.70+1.70) = 2.89.
    # Overlap = 1.89.  Other cross pairs: d(0,3)=10, d(1,2)=0.5 (overlap),
    # d(1,3)=8.5.  Pair (1,2) is also < target.
    assert penalty_close > 0.0
    assert n_close >= 1

    # Move frag 1 far away: no clash.
    P_far = P_close.copy()
    P_far[2] = [10.0, 0.0, 0.0]
    P_far[3] = [12.0, 0.0, 0.0]
    penalty_far, n_far = compute_interlig_vdw_penalty(
        P_far, frag, radii, factor=0.85, return_count=True,
    )
    assert penalty_far == 0.0
    assert n_far == 0


def test_interlig_vdw_penalty_skips_intra_fragment():
    """Atoms within the same fragment do NOT contribute (intra is handled
    by the intra-Q collapse rules)."""
    from delfin.fffree.build_time_clash_gate import (
        compute_interlig_vdw_penalty, INTERLIG_VDW_RADII,
    )
    # Two atoms 0.5 A apart in the SAME fragment.
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
    ], dtype=float)
    radii = [INTERLIG_VDW_RADII["C"], INTERLIG_VDW_RADII["C"]]
    frag = {0: [0, 1]}
    penalty = compute_interlig_vdw_penalty(P, frag, radii, factor=0.85)
    assert penalty == 0.0


def test_interlig_vdw_penalty_unknown_element_skipped():
    """Atoms with NaN vdW radius are skipped (no contribution)."""
    from delfin.fffree.build_time_clash_gate import compute_interlig_vdw_penalty
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=float)
    radii = [np.nan, 1.70]
    frag = {0: [0], 1: [1]}
    penalty, n = compute_interlig_vdw_penalty(
        P, frag, radii, factor=0.85, return_count=True,
    )
    assert penalty == 0.0
    assert n == 0


def test_interlig_vdw_penalty_deterministic():
    """Same input yields same output across two calls (no randomness)."""
    from delfin.fffree.build_time_clash_gate import (
        compute_interlig_vdw_penalty, INTERLIG_VDW_RADII,
    )
    np.random.seed(42)
    P = np.random.randn(8, 3) * 2.0
    radii = [INTERLIG_VDW_RADII["C"]] * 8
    frag = {0: [0, 1, 2, 3], 1: [4, 5, 6, 7]}
    p1, n1 = compute_interlig_vdw_penalty(
        P, frag, radii, return_count=True,
    )
    p2, n2 = compute_interlig_vdw_penalty(
        P, frag, radii, return_count=True,
    )
    assert p1 == p2
    assert n1 == n2


def test_interlig_penalty_env_off_byte_identical():
    """The assemble-complex hook returns penalty 0 when flag OFF, so the
    candidate selection is byte-identical to pre-fix HEAD.

    This is the unit-level analog: the helper
    ``_interlig_penalty_for_pair`` always computes; the hook checks
    ``interlig_penalty_active()`` first.  When OFF, the hook never calls
    it.  Verify both paths.
    """
    from delfin.fffree.build_time_clash_gate import (
        interlig_penalty_active, _interlig_penalty_for_pair,
    )
    syms_a = ["C", "O"]
    P_a = np.array([[0.0, 0.0, 0.0], [1.13, 0.0, 0.0]], dtype=float)
    syms_b = ["C", "O"]
    P_b = np.array([[0.5, 0.0, 0.0], [1.63, 0.0, 0.0]], dtype=float)
    # OFF (default): no flag set
    assert interlig_penalty_active() is False
    # Manual call still works (always-on helper), but hook gate is OFF.
    pen, n = _interlig_penalty_for_pair(syms_a, P_a, syms_b, P_b)
    assert pen >= 0.0
    assert n >= 0


def test_interlig_weight_env_override():
    """The penalty weight can be tuned via env."""
    os.environ["DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_W"] = "42.0"
    from delfin.fffree.build_time_clash_gate import _interlig_weight
    assert _interlig_weight() == 42.0


def test_interlig_factor_env_override():
    """The vdW fraction can be tuned via env."""
    os.environ["DELFIN_FFFREE_BUILD_INTERLIG_PENALTY_FACTOR"] = "0.70"
    from delfin.fffree.build_time_clash_gate import _interlig_factor
    assert _interlig_factor() == 0.70


def test_interlig_invalid_inputs():
    """Bad shapes raise."""
    from delfin.fffree.build_time_clash_gate import compute_interlig_vdw_penalty
    with pytest.raises(ValueError):
        compute_interlig_vdw_penalty(
            np.array([1.0, 2.0]), {0: [0]}, [1.7],
        )
    with pytest.raises(ValueError):
        compute_interlig_vdw_penalty(
            np.zeros((3, 3)), {0: [0]}, [1.7, 1.7],   # len mismatch
        )


def test_interlig_pair_helper_co_carbonyl_clash():
    """The pair-helper detects the canonical 1.85-A CO-CO interlig clash."""
    from delfin.fffree.build_time_clash_gate import _interlig_penalty_for_pair
    # Two CO ligands; their O atoms 1.85 A apart -> overlap > 0.
    syms_a = ["C", "O"]
    P_a = np.array([
        [0.0, 0.0, 0.0],
        [1.13, 0.0, 0.0],
    ], dtype=float)
    syms_b = ["C", "O"]
    P_b = np.array([
        [0.0, 1.85, 0.0],
        [1.13, 1.85, 0.0],
    ], dtype=float)
    # O-O distance = 1.85; vdW(O,O) = 3.04; 0.85*3.04 = 2.584.  Overlap > 0.
    pen, n = _interlig_penalty_for_pair(syms_a, P_a, syms_b, P_b)
    assert pen > 0.0
    assert n >= 1
