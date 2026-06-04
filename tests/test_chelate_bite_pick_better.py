"""Unit tests for the CHELATE_BITE per-structure pick-better gate.

Task #55-56 (2026-06-04): the σ-tail CHELATE_BITE constrained-per-chelate
embed was VALIDATED on ABUMET (11.97 -> 0.055 CShM) but smoke MIXED
(median down, mean up, catastrophic regressions D-YOMHIR 0.6 -> 36.1).
The fix is a per-structure pick-better gate: compute CShM for BOTH
constrained and unconstrained embeds, keep the lower-CShM one.  This
eliminates the catastrophic regressions while keeping the wins.

The gate is env-double-gated:
  * ``DELFIN_FFFREE_CHELATE_BITE=1`` keeps the constrained per-chelate
    embed active (existing behaviour),
  * ``DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER=1`` additionally computes
    the unconstrained embed and picks per-structure.

Default-OFF byte-identical: when either flag is unset, the helper's
``pick_better_active()`` returns False and assemble_from_config takes
the legacy path (no alt embed computed, no swap).

These tests exercise the public surface of
``delfin.fffree.chelate_bite_pick_better``; the wiring in
``assemble_from_config`` is exercised indirectly by HEAD's existing
smoke/coverage tests when both env flags are off (which is the default
in the test environment).
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

from delfin.fffree import chelate_bite_pick_better as PB


# ---------------------------------------------------------------------------
# Determinism + env-gate
# ---------------------------------------------------------------------------
def _clear_env(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_CHELATE_BITE", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER", raising=False)


def test_pick_better_active_default_off(monkeypatch):
    """Default env -> gate OFF (byte-identical legacy path)."""
    _clear_env(monkeypatch)
    assert PB.pick_better_active() is False


def test_pick_better_active_only_bite(monkeypatch):
    """Only CHELATE_BITE set -> pick-better STILL OFF (legacy bite behaviour)."""
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FFFREE_CHELATE_BITE", "1")
    assert PB.pick_better_active() is False


def test_pick_better_active_only_pick(monkeypatch):
    """Only PICK_BETTER set -> still OFF (bite is the umbrella gate)."""
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER", "1")
    assert PB.pick_better_active() is False


def test_pick_better_active_both_set(monkeypatch):
    _clear_env(monkeypatch)
    monkeypatch.setenv("DELFIN_FFFREE_CHELATE_BITE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER", "1")
    assert PB.pick_better_active() is True


# ---------------------------------------------------------------------------
# cshm_local: degenerate input never crashes
# ---------------------------------------------------------------------------
def test_cshm_local_none_returns_inf():
    assert PB.cshm_local(None, [0, 1], "OC-6 octahedron") == float("inf")


def test_cshm_local_too_few_donors_returns_inf():
    P = np.array([[1.0, 0, 0], [0, 1, 0]])
    assert PB.cshm_local(P, [], "OC-6 octahedron") == float("inf")
    assert PB.cshm_local(P, [0], "OC-6 octahedron") == float("inf")


def test_cshm_local_perfect_octahedron_near_zero():
    """6 donors on the ideal octahedron -> CShM should be very small."""
    P = np.array([
        [2.0, 0, 0], [-2.0, 0, 0],
        [0, 2.0, 0], [0, -2.0, 0],
        [0, 0, 2.0], [0, 0, -2.0],
    ])
    s = PB.cshm_local(P, list(range(6)), "OC-6 octahedron")
    assert s < 1e-6


def test_cshm_local_distorted_octahedron_finite_and_positive():
    """Perturb one donor and CShM must be finite and > 0."""
    P = np.array([
        [2.0, 0, 0], [-2.0, 0, 0],
        [0, 2.0, 0], [0, -2.0, 0],
        [0, 0, 2.0], [0.5, 0.5, -1.8],   # perturbed
    ])
    s = PB.cshm_local(P, list(range(6)), "OC-6 octahedron")
    assert np.isfinite(s)
    assert s > 0.0


# ---------------------------------------------------------------------------
# pick_better core contract
# ---------------------------------------------------------------------------
def _two_donor_block(d1, d2):
    """Build a 2-donor block with metal at origin and 2 donor positions."""
    return np.array([d1, d2], dtype=float)


def test_pick_better_constrained_wins_when_lower_cshm():
    """ABUMET-like case: constrained = ideal octahedral cis (90°),
    unconstrained = distorted (acute) -> constrained wins."""
    md = 2.0
    # constrained donors: cis-octahedral 90° apart on the x,y axes
    Qc = _two_donor_block([md, 0, 0], [0, md, 0])
    # unconstrained: 60° apart (over-stretched ring buckled)
    Qu = _two_donor_block([md, 0, 0],
                          [md * math.cos(math.radians(60)),
                           md * math.sin(math.radians(60)), 0])
    coords, src, cc, cu = PB.pick_better(
        Qc, Qu, donor_idxs=[0, 1], geometry="OC-6 octahedron",
    )
    # The geometry only sees TWO donors but cshm vs OC-6 needs 6 -> returns 0.0
    # (polyhedra.cshm degenerate guard).  To get a meaningful comparison we
    # invoke a CN2 geometry instead.
    # → fix: use L-2 linear with 180° donors.  Re-run with proper case:
    # Constrained: 180° apart on x-axis (linear)
    Qc = _two_donor_block([md, 0, 0], [-md, 0, 0])
    # Unconstrained: 100° apart (heavily bent)
    Qu = _two_donor_block(
        [md, 0, 0],
        [md * math.cos(math.radians(100)), md * math.sin(math.radians(100)), 0],
    )
    coords, src, cc, cu = PB.pick_better(
        Qc, Qu, donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert src == "constrained"
    assert cc <= cu
    assert np.allclose(coords, Qc)


def test_pick_better_catastrophic_unconstrained_wins():
    """D-YOMHIR-class: constrained forced-90° geometry blows up CShM,
    unconstrained natural-bite is cleaner -> unconstrained wins."""
    md = 2.0
    # constrained: heavily distorted (forced-90° broke the rigid ring backbone,
    # so the donors actually landed far off the targets -> we model that by
    # giving constrained a *bad* geometry vs CN2 linear).
    Qc = _two_donor_block(
        [md, 0, 0],
        [md * math.cos(math.radians(60)), md * math.sin(math.radians(60)), 0],
    )
    # unconstrained: cleanly linear (closer to ideal L-2 = 180°)
    Qu = _two_donor_block([md, 0, 0], [-md, 0, 0])
    coords, src, cc, cu = PB.pick_better(
        Qc, Qu, donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert src == "unconstrained"
    assert cu <= cc
    assert np.allclose(coords, Qu)


def test_pick_better_tie_constrained_wins():
    """Exact tie -> deterministic rule = constrained wins."""
    md = 2.0
    Q = _two_donor_block([md, 0, 0], [-md, 0, 0])
    # Same coordinates: tie
    coords, src, cc, cu = PB.pick_better(
        Q, Q.copy(), donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert src == "constrained"
    assert cc == cu


def test_pick_better_none_unconstrained_keeps_constrained():
    """Unconstrained embed failed -> keep the constrained one."""
    md = 2.0
    Qc = _two_donor_block([md, 0, 0], [-md, 0, 0])
    coords, src, cc, cu = PB.pick_better(
        Qc, None, donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert src == "constrained"
    assert not math.isfinite(cu)
    assert np.allclose(coords, Qc)


def test_pick_better_none_constrained_keeps_unconstrained():
    """Constrained embed failed -> fall back to unconstrained."""
    md = 2.0
    Qu = _two_donor_block([md, 0, 0], [-md, 0, 0])
    coords, src, cc, cu = PB.pick_better(
        None, Qu, donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert src == "unconstrained"
    assert not math.isfinite(cc)
    assert np.allclose(coords, Qu)


def test_pick_better_both_none_returns_none():
    coords, src, cc, cu = PB.pick_better(
        None, None, donor_idxs=[0, 1], geometry="L-2 linear",
    )
    assert coords is None
    assert src == "none"
    assert not math.isfinite(cc) and not math.isfinite(cu)


# ---------------------------------------------------------------------------
# M-D invariant defence-in-depth
# ---------------------------------------------------------------------------
def test_md_invariant_preserved_true_when_donors_unchanged():
    A = np.array([[2.0, 0, 0], [-2.0, 0, 0], [1.0, 1.0, 0.5]])
    B = A.copy()
    B[2] = [3.0, 0, 0]   # non-donor moved
    assert PB.md_invariant_preserved(A, B, donor_idxs=[0, 1])


def test_md_invariant_preserved_false_when_donor_moved():
    A = np.array([[2.0, 0, 0], [-2.0, 0, 0]])
    B = A.copy()
    B[0] = [3.5, 0, 0]   # M-D distance changed 2.0 -> 3.5
    assert not PB.md_invariant_preserved(A, B, donor_idxs=[0, 1], tol=0.05)


def test_md_invariant_shape_mismatch_returns_false():
    A = np.array([[2.0, 0, 0], [-2.0, 0, 0]])
    B = np.array([[2.0, 0, 0]])
    assert not PB.md_invariant_preserved(A, B, donor_idxs=[0])


# ---------------------------------------------------------------------------
# Determinism: same input -> same pick across calls
# ---------------------------------------------------------------------------
def test_pick_better_is_deterministic():
    md = 2.0
    Qc = _two_donor_block([md, 0, 0], [-md, 0, 0])
    Qu = _two_donor_block(
        [md, 0, 0],
        [md * math.cos(math.radians(100)), md * math.sin(math.radians(100)), 0],
    )
    out1 = PB.pick_better(Qc, Qu, donor_idxs=[0, 1], geometry="L-2 linear")
    out2 = PB.pick_better(Qc, Qu, donor_idxs=[0, 1], geometry="L-2 linear")
    assert out1[1] == out2[1]
    assert out1[2] == out2[2]
    assert out1[3] == out2[3]
    assert np.array_equal(out1[0], out2[0])


# ---------------------------------------------------------------------------
# Public API smoke (pick_better signature + return shape)
# ---------------------------------------------------------------------------
def test_pick_better_return_shape():
    md = 2.0
    Qc = _two_donor_block([md, 0, 0], [-md, 0, 0])
    Qu = _two_donor_block([md, 0, 0], [-md, 0, 0])
    out = PB.pick_better(Qc, Qu, donor_idxs=[0, 1], geometry="L-2 linear")
    assert len(out) == 4
    coords, src, cc, cu = out
    assert isinstance(coords, np.ndarray)
    assert src in {"constrained", "unconstrained", "none"}
    assert isinstance(cc, float)
    assert isinstance(cu, float)
