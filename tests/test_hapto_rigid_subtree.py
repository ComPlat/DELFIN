"""Unit tests for delfin.fffree.hapto_rigid_subtree (Mission A3, Phase G13b).

Validates:
  1. Determinism (same input -> identical output across calls)
  2. OFF byte-identity (env-flag unset -> input == output exactly)
  3. ON byte-difference (env-flag set -> rigid extension applied)
  4. Rigid invariant: intra-block pairwise distances preserved exactly
  5. Substituent subtree drag: methyl-C and methyl-H on Cp ride with ring
  6. Boundary respect: metal + σ donors NEVER move
  7. Multi-hapto (ferrocene-like): two Cp rings each get their own block,
     the BFS does NOT bridge through the metal
  8. Fused ring boundary: BFS stops at ``other_ring_atoms``
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest


# Determinism
os.environ.setdefault("PYTHONHASHSEED", "0")


@pytest.fixture
def env_off(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_HAPTO_RIGID_SUBTREE", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", raising=False)


@pytest.fixture
def env_on(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_RIGID_SUBTREE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", "1")
    monkeypatch.delenv("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", raising=False)


def _build_methyl_cp_iron(z_ring: float = 1.95, r_ring: float = 1.21):
    """Construct ferrocene-like fragment: Fe + Cp (5 C + 5 H) +
    one methyl substituent on ring atom 0 (3 H + 1 C).

    Layout indices:
      0 : Fe
      1-5: ring C
      6-10: ring H
      11 : methyl C
      12-14: methyl H
    """
    symbols = (
        ["Fe"]
        + ["C"] * 5
        + ["H"] * 5
        + ["C", "H", "H", "H"]
    )
    P = np.zeros((len(symbols), 3))
    P[0] = [0.0, 0.0, 0.0]
    for i in range(5):
        th = i * 2.0 * math.pi / 5.0
        P[1 + i] = [r_ring * math.cos(th), r_ring * math.sin(th), z_ring]
        # H sits 1.08 Å outward from ring C in the ring plane
        P[6 + i] = [(r_ring + 1.08) * math.cos(th),
                    (r_ring + 1.08) * math.sin(th), z_ring]
    # Methyl C attached to ring C atom 1 (index in P = 1), in-plane offset
    P[11] = P[1] + np.array([1.50, 0.0, 0.0])
    P[12] = P[11] + np.array([1.08, 0.0, 0.0])
    P[13] = P[11] + np.array([-0.36, 1.02, 0.0])
    P[14] = P[11] + np.array([-0.36, -0.51, 0.88])
    return symbols, P


# --------------------------------------------------------------------------- #
# 1. Determinism / OFF byte-identity                                          #
# --------------------------------------------------------------------------- #
def test_module_imports(env_off):
    from delfin.fffree.hapto_rigid_subtree import (
        collect_rigid_subtree_atoms,
        apply_rigid_subtree_translation,
        rigid_subtree_active,
    )
    assert rigid_subtree_active() is False


def test_off_byte_identity_translation(env_off):
    """When the env flag is OFF, apply_rigid_subtree_translation returns
    coords unchanged."""
    from delfin.fffree.hapto_rigid_subtree import apply_rigid_subtree_translation
    rng = np.random.RandomState(7)
    P = rng.randn(15, 3)
    P0 = P.copy()
    P_out = apply_rigid_subtree_translation(P, [1, 2, 3], np.array([5.0, 0.0, 0.0]))
    assert np.array_equal(P_out, P0)


def test_off_byte_identity_honest(env_off):
    """When NEITHER honest_construction nor rigid_subtree is set,
    apply_hapto_honest returns input copy byte-identical."""
    from delfin.fffree.hapto_honest_construction import apply_hapto_honest
    symbols, P = _build_methyl_cp_iron()
    P0 = P.copy()
    P_out, n = apply_hapto_honest(symbols, P, metal_idx=0)
    assert n == 0
    assert np.array_equal(P_out, P0)


def test_determinism(env_on):
    """Same input twice → bit-identical outputs."""
    from delfin.fffree.hapto_rigid_subtree import apply_rigid_subtree_translation
    rng = np.random.RandomState(11)
    P = rng.randn(20, 3)
    delta = np.array([1.5, -0.7, 0.3])
    P_a = apply_rigid_subtree_translation(P, [2, 4, 6], delta)
    P_b = apply_rigid_subtree_translation(P, [2, 4, 6], delta)
    assert np.array_equal(P_a, P_b)


# --------------------------------------------------------------------------- #
# 2. Rigid invariant                                                          #
# --------------------------------------------------------------------------- #
def test_rigid_invariant_translation(env_on):
    """All pairwise distances in the rigid block are preserved exactly
    after translation."""
    from delfin.fffree.hapto_rigid_subtree import apply_rigid_subtree_translation
    rng = np.random.RandomState(13)
    P = rng.randn(12, 3)
    block = [2, 3, 4, 7]
    delta = np.array([3.14, -2.71, 0.5])
    P_out = apply_rigid_subtree_translation(P, block, delta)
    for a in block:
        for b in block:
            if a >= b:
                continue
            d_in = float(np.linalg.norm(P[a] - P[b]))
            d_out = float(np.linalg.norm(P_out[a] - P_out[b]))
            assert abs(d_in - d_out) < 1e-12, (a, b, d_in, d_out)


def test_rigid_invariant_rotation(env_on):
    """Pairwise distances preserved after rigid rotation about a pivot."""
    from delfin.fffree.hapto_rigid_subtree import apply_rigid_subtree_transform
    rng = np.random.RandomState(17)
    P = rng.randn(10, 3)
    block = [1, 4, 5, 8, 9]
    # 30° rotation about z
    th = math.pi / 6
    R = np.array([[math.cos(th), -math.sin(th), 0.0],
                  [math.sin(th), math.cos(th), 0.0],
                  [0.0, 0.0, 1.0]])
    pivot = P[block].mean(axis=0)
    P_out = apply_rigid_subtree_transform(P, block, pivot=pivot, rotation=R)
    for a in block:
        for b in block:
            if a >= b:
                continue
            d_in = float(np.linalg.norm(P[a] - P[b]))
            d_out = float(np.linalg.norm(P_out[a] - P_out[b]))
            assert abs(d_in - d_out) < 1e-9, (a, b, d_in, d_out)


# --------------------------------------------------------------------------- #
# 3. Substituent subtree drag (the heart of the mission)                      #
# --------------------------------------------------------------------------- #
def test_collect_includes_methyl_subtree(env_on):
    """For a methyl-Cp, the rigid set for the ring must include the
    methyl carbon AND its three H atoms."""
    from delfin.fffree.hapto_rigid_subtree import collect_rigid_subtree_atoms
    symbols, P = _build_methyl_cp_iron()
    rigid = collect_rigid_subtree_atoms(
        symbols, P,
        ring_atom_indices=[1, 2, 3, 4, 5],
        metal_idx=0,
    )
    # Ring 1-5, ring H 6-10, methyl 11, methyl H 12-14.
    expected = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}
    assert set(rigid) == expected


def test_honest_drag_methyl_with_ring(env_on):
    """When honest correction snaps the Cp centroid, the methyl carbon
    AND its H atoms MUST travel by the same delta.  The C(ring)-C(methyl)
    bond is then preserved EXACTLY."""
    from delfin.fffree.hapto_honest_construction import apply_hapto_honest
    symbols, P = _build_methyl_cp_iron()
    d_bond_in = float(np.linalg.norm(P[11] - P[1]))
    P_out, n = apply_hapto_honest(symbols, P, metal_idx=0)
    assert n == 1
    ring_delta = P_out[1] - P[1]
    methyl_delta = P_out[11] - P[11]
    # rigid drag: methyl C delta == ring C[1] delta
    assert np.allclose(ring_delta, methyl_delta, atol=1e-12), (
        ring_delta, methyl_delta
    )
    # methyl H atoms ride too
    for h in (12, 13, 14):
        h_delta = P_out[h] - P[h]
        assert np.allclose(ring_delta, h_delta, atol=1e-12), (h, ring_delta, h_delta)
    # Attachment bond preserved exactly
    d_bond_out = float(np.linalg.norm(P_out[11] - P_out[1]))
    assert abs(d_bond_in - d_bond_out) < 1e-12, (d_bond_in, d_bond_out)


def test_honest_off_does_not_drag(env_off):
    """When DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION=1 is set but the
    rigid-subtree flag is OFF, the legacy behaviour (ring + ring-H only)
    is preserved -- methyl substituent stays behind.  This proves the
    extension is gated cleanly behind the flag."""
    os.environ["DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION"] = "1"
    os.environ.pop("DELFIN_FFFREE_HAPTO_RIGID_SUBTREE", None)
    os.environ.pop("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", None)
    try:
        from delfin.fffree.hapto_honest_construction import apply_hapto_honest
        symbols, P = _build_methyl_cp_iron()
        P_out, n = apply_hapto_honest(symbols, P, metal_idx=0)
        assert n == 1
        # ring moved
        ring_delta = float(np.linalg.norm(P_out[1] - P[1]))
        assert ring_delta > 0.01
        # methyl did NOT move (legacy behaviour)
        methyl_delta = float(np.linalg.norm(P_out[11] - P[11]))
        assert methyl_delta < 1e-9, methyl_delta
    finally:
        os.environ.pop("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", None)


# --------------------------------------------------------------------------- #
# 4. Boundary respect                                                         #
# --------------------------------------------------------------------------- #
def test_metal_never_moves(env_on):
    """Metal stays at the origin no matter what."""
    from delfin.fffree.hapto_honest_construction import apply_hapto_honest
    symbols, P = _build_methyl_cp_iron()
    P_out, n = apply_hapto_honest(symbols, P, metal_idx=0)
    assert n == 1
    assert np.allclose(P_out[0], P[0], atol=1e-12)


def test_bfs_stops_at_metal(env_on):
    """The BFS does NOT cross a metal edge.  For a ferrocene-like
    [Cp-Fe-Cp], the BFS started from ring A must NOT touch ring B
    atoms (they reach only via the metal)."""
    from delfin.fffree.hapto_rigid_subtree import (
        collect_rigid_subtree_atoms, build_geometric_adjacency,
    )
    # Build ferrocene: Fe + 2 Cp rings (top + bottom) -- no methyl
    symbols = ["Fe"] + ["C"] * 5 + ["C"] * 5
    P = np.zeros((11, 3))
    P[0] = [0.0, 0.0, 0.0]
    r = 1.21
    for i in range(5):
        th = i * 2.0 * math.pi / 5.0
        P[1 + i] = [r * math.cos(th), r * math.sin(th), 1.65]
        P[6 + i] = [r * math.cos(th), r * math.sin(th), -1.65]
    # Both rings have their own atom indices; both reach the metal only.
    rigid_top = collect_rigid_subtree_atoms(
        symbols, P,
        ring_atom_indices=[1, 2, 3, 4, 5],
        metal_idx=0,
        other_ring_atoms=[6, 7, 8, 9, 10],
    )
    # Top ring's rigid set must equal {1,...,5}; the bottom ring (6-10) is
    # not reachable except through the metal (which is excluded).
    assert set(rigid_top) == {1, 2, 3, 4, 5}


def test_frozen_extra_boundary(env_on):
    """An atom passed in ``frozen_extra`` is NOT pulled into the rigid
    set even when it is graph-bonded to a ring atom (e.g. a σ donor on
    a chelating phosphine arm)."""
    from delfin.fffree.hapto_rigid_subtree import collect_rigid_subtree_atoms
    symbols, P = _build_methyl_cp_iron()
    rigid = collect_rigid_subtree_atoms(
        symbols, P,
        ring_atom_indices=[1, 2, 3, 4, 5],
        metal_idx=0,
        frozen_extra={11},  # mark methyl carbon as frozen σ donor
    )
    # Methyl C is blocked; its H atoms can't be reached without going
    # through the methyl C (which is in frozen) so they stay out too.
    assert 11 not in rigid
    assert 12 not in rigid and 13 not in rigid and 14 not in rigid
    # Ring + ring H still in
    assert {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}.issubset(set(rigid))


# --------------------------------------------------------------------------- #
# 5. Sanity: detect-collect-translate round-trip                              #
# --------------------------------------------------------------------------- #
def test_translation_then_rotation_invariant(env_on):
    """A full translate-then-rotate of a rigid block preserves all
    intra-block pairwise distances exactly."""
    from delfin.fffree.hapto_rigid_subtree import (
        apply_rigid_subtree_translation,
        apply_rigid_subtree_transform,
    )
    rng = np.random.RandomState(23)
    P = rng.randn(15, 3)
    block = [3, 4, 5, 6, 7, 8]
    delta = np.array([1.0, 2.0, -0.5])
    P1 = apply_rigid_subtree_translation(P, block, delta)
    R = np.array([
        [math.cos(0.4), -math.sin(0.4), 0.0],
        [math.sin(0.4),  math.cos(0.4), 0.0],
        [0.0, 0.0, 1.0],
    ])
    pivot = P1[block].mean(axis=0)
    P2 = apply_rigid_subtree_transform(P1, block, pivot=pivot, rotation=R)
    for a in block:
        for b in block:
            if a >= b:
                continue
            d_in = float(np.linalg.norm(P[a] - P[b]))
            d_out = float(np.linalg.norm(P2[a] - P2[b]))
            assert abs(d_in - d_out) < 1e-9, (a, b, d_in, d_out)
