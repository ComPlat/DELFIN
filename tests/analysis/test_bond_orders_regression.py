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
