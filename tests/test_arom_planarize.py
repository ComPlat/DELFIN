"""Unit tests for delfin.manta._arom_planarize (Iter-33 universal aromatic
ring-SYSTEM planarisation).

Geometry-only (no SMILES build), so these run fast and deterministically.
We construct a puckered benzene + a coordinated heteroaromatic donor ring
and assert: pendant rings flatten, coordinated donors are anchor-preserved
(M-D distance unchanged), output is finite + deterministic, and an already-flat
/ non-aromatic input is returned bit-identically.
"""
from __future__ import annotations

import math

import numpy as np

from delfin.manta._arom_planarize import correct_xyz, correct_results, _ring_oop
from delfin.manta._pi_h_projector import _parse_xyz, _build_geometric_adjacency


def _benzene(z_pucker: float = 0.0, r: float = 1.39) -> str:
    """Flat or puckered benzene (alternating +/- z) as a 6-atom XYZ block."""
    lines = ["6", "benzene"]
    for k in range(6):
        ang = 2 * math.pi * k / 6
        x = r * math.cos(ang)
        y = r * math.sin(ang)
        z = z_pucker * (1 if k % 2 == 0 else -1)
        lines.append(f"C {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


def test_pendant_benzene_flattens():
    xyz = _benzene(z_pucker=0.12)  # ~0.12 Å pucker, in actionable band
    syms, pts, _ = _parse_xyz(xyz)
    ring = tuple(range(6))
    oop_before = _ring_oop(pts, ring)
    assert oop_before > 0.10

    out = correct_xyz(xyz)
    syms2, pts2, _ = _parse_xyz(out)
    oop_after = _ring_oop(pts2, ring)
    assert oop_after < 0.02, (oop_before, oop_after)
    assert np.all(np.isfinite(pts2))


def test_already_flat_is_bit_identical():
    xyz = _benzene(z_pucker=0.0)
    assert correct_xyz(xyz) == xyz


def test_non_aromatic_is_bit_identical():
    # A short methane-ish blob: no 5/6 aromatic ring → untouched.
    xyz = (
        "5\nmethane\n"
        "C 0.000000 0.000000 0.000000\n"
        "H 0.629000 0.629000 0.629000\n"
        "H -0.629000 -0.629000 0.629000\n"
        "H -0.629000 0.629000 -0.629000\n"
        "H 0.629000 -0.629000 -0.629000\n"
    )
    assert correct_xyz(xyz) == xyz


def test_saturated_ring_untouched():
    # Cyclohexane-like (mean bond ~1.54) puckered ring must NOT be flattened.
    xyz = _benzene(z_pucker=0.20, r=1.54 / (2 * math.sin(math.pi / 6)))
    # bond length here ~1.54 → above aromatic gate → no flatten
    assert correct_xyz(xyz) == xyz


def test_md_distance_preserved_for_coordinated_donor():
    """A pyridine-like C5N ring coordinated to a metal: the N donor must not
    move (M-N distance preserved), even while the ring is flattened."""
    # Build a slightly puckered 6-ring (5 C + 1 N at index 5) and put a metal
    # 2.10 Å from the N along -y.
    r = 1.39
    atoms = []
    coords = []
    for k in range(6):
        ang = 2 * math.pi * k / 6
        x = r * math.cos(ang)
        y = r * math.sin(ang)
        z = 0.13 * (1 if k % 2 == 0 else -1)
        atoms.append("N" if k == 5 else "C")
        coords.append([x, y, z])
    # place metal below the N (atom 5)
    nx, ny, nz = coords[5]
    metal = [nx, ny - 2.10, nz]
    atoms.append("Fe")
    coords.append(metal)
    lines = [str(len(atoms)), "pyridine-Fe"]
    for s, (x, y, z) in zip(atoms, coords):
        lines.append(f"{s} {x:12.6f} {y:12.6f} {z:12.6f}")
    xyz = "\n".join(lines) + "\n"

    syms, pts, _ = _parse_xyz(xyz)
    m_idx, n_idx = 6, 5
    d_before = float(np.linalg.norm(pts[m_idx] - pts[n_idx]))

    out = correct_xyz(xyz)
    syms2, pts2, _ = _parse_xyz(out)
    d_after = float(np.linalg.norm(pts2[m_idx] - pts2[n_idx]))
    # M-N distance preserved within the M-D invariant tolerance
    assert abs(d_after - d_before) < 0.05, (d_before, d_after)
    assert np.all(np.isfinite(pts2))


def test_deterministic():
    xyz = _benzene(z_pucker=0.15)
    a = correct_xyz(xyz)
    b = correct_xyz(xyz)
    assert a == b


def test_correct_results_passthrough_on_empty():
    assert correct_results(None, []) == []
    res = [(_benzene(0.0), "flat")]
    assert correct_results(None, res) == res


def _puckered_toluene():
    """Puckered benzene whose ring C0 carries a methyl lifted 0.30 A out of the ring plane, with its
    three hydrogens at 1.09 A (tetrahedral)."""
    lines, r = [], 1.39
    for k in range(6):
        ang = 2 * math.pi * k / 6
        z = 0.12 * (1 if k % 2 == 0 else -1)
        lines.append(("C", np.array([r * math.cos(ang), r * math.sin(ang), z])))
    c0 = lines[0][1]
    cm = c0 + np.array([1.47, 0.0, 0.30])
    u = (cm - c0) / np.linalg.norm(cm - c0)
    a = np.cross(u, [0.0, 0.0, 1.0]); a /= np.linalg.norm(a)
    b = np.cross(u, a)
    for k in range(3):
        phi = 2 * math.pi * k / 3
        d = math.cos(math.radians(70.5)) * u + math.sin(math.radians(70.5)) * (math.cos(phi) * a + math.sin(phi) * b)
        lines.append(("H", cm + 1.09 * d))
    lines.insert(6, ("C", cm))
    body = "\n".join(f"{s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for s, p in lines)
    return f"{len(lines)}\ntoluene\n{body}\n"


def _ch(pts):
    return [float(np.linalg.norm(pts[h] - pts[6])) for h in (7, 8, 9)]


def test_dragged_methyl_leaves_its_hydrogens_without_the_switch(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_AROM_PLANARIZE_CARRY_H", raising=False)
    xyz = _puckered_toluene()
    _s, p0, _ = _parse_xyz(xyz)
    _s, p1, _ = _parse_xyz(correct_xyz(xyz))
    assert abs(p1[6][2] - p0[6][2]) > 0.05                  # the methyl C was dragged
    assert np.abs(p1[7:10] - p0[7:10]).max() == 0.0         # its hydrogens were not
    assert max(abs(d - 1.09) for d in _ch(p1)) > 0.05       # -> stretched / squeezed C-H


def test_dragged_methyl_carries_its_hydrogens_with_the_switch(monkeypatch):
    xyz = _puckered_toluene()
    monkeypatch.delenv("DELFIN_FFFREE_AROM_PLANARIZE_CARRY_H", raising=False)
    _s, poff, _ = _parse_xyz(correct_xyz(xyz))
    monkeypatch.setenv("DELFIN_FFFREE_AROM_PLANARIZE_CARRY_H", "1")
    _s, pon, _ = _parse_xyz(correct_xyz(xyz))
    assert np.abs(pon[:7] - poff[:7]).max() < 1e-9          # ring + methyl C exactly as without the switch
    assert max(abs(d - 1.09) for d in _ch(pon)) < 1e-6      # every C-H kept
