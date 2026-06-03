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
