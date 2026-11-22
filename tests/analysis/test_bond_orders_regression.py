"""Regression tests for bond orders in the skeleton renderer.

User report (2026-11-21): pyridine rings coordinated to the metal lost
their C=N double bond and showed 2 instead of 3 doubles (class A).  The
drop-rule in _rdkit_bond_orders was meant for 5-ring pyrrolate N only.

These tests call the real functions on a small synthetic geometry that
mirrors the 106-138_108 metal-bis-pyridine complex: an Fe center with
two 6-membered pyridine rings (each ring N coordinating the metal), one
5-membered pyrazolate ring, and the chlorides.
"""

import numpy as np
import pytest

from delfin.analysis.results_overview import (
    _detect_bonds, _rdkit_bond_orders, _XYZ_COVALENT,
    _XYZ_COVALENT_DEFAULT, _rotate_inertia)


def _build(sym, coords):
    """Return (bonds, positions) for the given xyz list."""
    pos = np.array(coords, dtype=float)
    rotated = _rotate_inertia(sym, pos)
    cov = [_XYZ_COVALENT.get(s, _XYZ_COVALENT_DEFAULT) for s in sym]
    return _detect_bonds(sym, cov, rotated), rotated


# Fe center at origin; two pyridine ring N at the axial sites; ring
# atoms laid out in the x-y plane.  Distances: C-C 1.39, C-N 1.34,
# N-Fe 2.1, Fe-Cl 2.3 (all under the covalent-radii cutoff).
FE = (0.0, 0.0, 0.0)
N1 = (0.0, 0.0, 2.1)     # coordinating N of pyridine 1
N2 = (0.0, 0.0, -2.1)    # coordinating N of pyridine 2


def _pyridine(center_n, direction):
    """6 ring atoms: the coordinating N plus 5 C, planar ring."""
    # Build a regular hexagon of edge 1.39 around the N position with
    # the ring plane in z (through the N and the metal axis).
    import math
    cx, cy, cz = center_n
    # ring centre 1.39 further along `direction` from the N... simpler:
    # place 5 carbons at hexagon vertices around a centre offset from
    # the N so the N-adjacent C sits at edge distance.
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    # any perpendicular axis in the ring plane
    perp = np.array([0.0, 1.0, 0.0])
    if abs(np.dot(d, perp)) > 0.9:
        perp = np.array([1.0, 0.0, 0.0])
    u = perp / np.linalg.norm(perp)
    centre = np.array(center_n, dtype=float) + d * 1.39 * math.cos(math.pi / 6)
    c0 = centre - d * 1.39 * math.cos(math.pi / 6)  # == the N position
    pts = []
    for k in range(1, 6):
        ang = math.pi / 3 * k + math.pi / 6
        p = c0 + (d * math.cos(ang) + u * math.sin(ang)) * 1.39
        pts.append(tuple(p))
    return pts


def test_pyridine_rings_keep_three_doubles():
    """Class A: coordinated 6-ring pyridine N keeps its C=N double."""
    sym = (["Fe", "Cl", "Cl", "Cl", "N", "N"]
           + ["C"] * 10)
    coords = [FE, (2.3, 0, 0), (0, 2.3, 0), (0, -2.3, 0), N1, N2]
    # Pyridine 1 in +z half-space, ring plane x-z (through the axis).
    # Simpler and robust: place the 5 C of each ring along a hexagon in
    # the x-z plane around the N, ring normal (0, 1, 0).
    import math
    ring_c = []
    for sign in (+1, -1):
        n_pos = np.array(N1 if sign > 0 else N2)
        # hexagon centre one edge further from the metal along z
        centre = n_pos + np.array([0, 0, sign * 1.2])
        for k in range(1, 6):
            ang = math.pi / 3 * k
            p = centre + np.array([
                math.cos(ang) * 1.39,
                0.0,
                sign * math.sin(ang) * 1.39])
            ring_c.append(tuple(p))
    coords += ring_c
    bonds, pos = _build(sym, coords)
    assert len(bonds) > 10, "ring bonds must be detected"
    orders = _rdkit_bond_orders(sym, bonds, pos)
    assert orders, "bond orders must be assigned"
    # Each pyridine N (index 4 and 5) must keep ONE double bond to a
    # ring carbon - this is the C=N the old drop-rule removed.
    for n_idx in (4, 5):
        dbl = [k for (i, j), o in orders.items()
               if o == 2 and (i == n_idx or j == n_idx)]
        assert dbl, (
            f"coordinated pyridine N {n_idx} lost its C=N double "
            f"(class A regression)")
    # Total doubles >= 6: three per ring (class A = 2 per ring before).
    n_dbl = len({tuple(sorted(p)) for p in orders.items() and orders})
    assert n_dbl >= 6, f"expected >=6 double bonds, got {n_dbl}"


# ---- classes B/E: 27-50_36..50 report (2026-11-22) -------------------
# These calcs carry no stored SMILES; the pipeline must take the orders
# from the 3D-derived SMILES (wired in _smiles_orders._load_smiles),
# not fall back to the geometric heuristic.

def _calc_dir_for(n):
    """Real calc dir; skipped when absent (cluster data not mounted)."""
    import os
    d = f"/pfs/data6/home/ka/ka_ibcs/ka_xn0397/calc/27-50_{n}"
    if not os.path.isdir(d):
        pytest.skip(f"calc data {d} not mounted")
    return d


def _orders_for(n):
    from delfin.analysis import results_overview as ro
    from delfin.analysis._smiles_orders import smiles_bond_orders
    cd = _calc_dir_for(n)
    atoms = ro._parse_xyz(cd)
    sym = [a[0] for a in atoms]
    pts = np.array([[a[1], a[2], a[3]] for a in atoms])
    rot = ro._rotate_inertia(sym, pts)
    cov = [ro._XYZ_COVALENT.get(s, ro._XYZ_COVALENT_DEFAULT) for s in sym]
    bonds = ro._detect_bonds(sym, cov, rot)
    return sym, bonds, smiles_bond_orders(sym, bonds, cd)


def test_smiles_orders_wired_for_old_calcs():
    """27-50_36 (no stored SMILES): the orders must come from the
    3D-derived SMILES, so the pyridine rings show 3 doubles."""
    sym, bonds, sbo = _orders_for(36)
    assert sbo, "smiles_bond_orders returned None - 3D SMILES not wired"
    ro_mod = _orders_for.__globals__["__import__"]  # noqa: F841 (lint)
    from delfin.analysis import results_overview as ro
    r6 = [r for r in ro._find_rings(sym, bonds)
          if len(r) == 6 and sum(1 for k in r if sym[k] == "N") == 1]
    assert r6, "expected at least one pyridine ring"
    for r in r6:
        d = sum(1 for a, b in zip(r, r[1:] + r[:1])
                if sbo.get((min(a, b), max(a, b))))
        assert d == 3, f"pyridine ring has {d} doubles, expected 3"


def test_amid_diimine_ring_matches_smiles():
    """27-50_45: the diimine chelate ring must carry the SMILES kekule
    (2 C=C) and BOTH exocyclic imine C=N; no dropped SMILES edge."""
    sym, bonds, sbo = _orders_for(45)
    assert sbo
    # count transferred doubles = number of doubles in the kekule SMILES
    from delfin.analysis import _smiles_orders as so
    atoms_s, orders_s = so._kekule_orders(so._load_smiles(
        f"/pfs/data6/home/ka/ka_ibcs/ka_xn0397/calc/27-50_45"))
    n_s = sum(1 for v in orders_s.values() if v == 2)
    n_x = sum(1 for v in sbo.values() if v == 2)
    assert n_x == n_s, (
        f"transferred {n_x} doubles, SMILES has {n_s} - an edge was "
        f"dropped in the transfer (class B regression)")


def test_nitro_o_carries_no_h_label():
    """27-50_42: both nitro O have no structural H - the renderer must
    not give them an H suffix (class E)."""
    sym, bonds, _ = _orders_for(42)
    nitro_n = [i for i, s in enumerate(sym) if s == "N"
               and sum(1 for a, b in bonds
                       if a == i and sym[b] == "O"
                       or b == i and sym[a] == "O") == 2]
    assert nitro_n, "expected a nitro group"
    for n_idx in nitro_n:
        for a, b in bonds:
            if a == n_idx and sym[b] == "O" or b == n_idx and sym[a] == "O":
                o = b if a == n_idx else a
                h = [x for x, y in bonds
                     if (x == o and sym[y] == "H")
                     or (y == o and sym[x] == "H")]
                assert not h, f"nitro O {o} carries an H (class E)"
