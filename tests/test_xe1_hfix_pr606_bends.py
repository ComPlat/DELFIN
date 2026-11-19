"""fix_hydrogens, operator real-run round 7 on #606.

Measured on the PR head e7bc84eec: the H on a two-valent O or S with one heavy
neighbour was placed on the bond axis (methanol, phenol C-O-H 180.0; CH3SH
C-S-H 180.0), and ethanol's CH2 got its two H in the C-C-O plane, one of them
0.42 A from the O. Real values: C-O-H 108.5, C-S-H 96.5; an sp3 CH2 has one H
on each side of the bond plane (H-C-H 109.5).
"""
from __future__ import annotations

import math

import pytest

from delfin.structure_edit import fix_hydrogens


def _xyz(atoms):
    return f"{len(atoms)}\n\n" + "\n".join(f"{s} {x} {y} {z}" for s, x, y, z in atoms) + "\n"


def _atoms(text):
    return [(ln.split()[0], tuple(map(float, ln.split()[1:4])))
            for ln in text.strip().splitlines()[2:]]


def _angle(a, b, c):
    v1 = [a[k] - b[k] for k in range(3)]
    v2 = [c[k] - b[k] for k in range(3)]
    cos = sum(v1[k] * v2[k] for k in range(3)) / (math.dist(a, b) * math.dist(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def _hydrogens_of(atoms, i, cutoff=1.45):
    return [j for j, (s, p) in enumerate(atoms) if s == "H" and math.dist(p, atoms[i][1]) < cutoff]


PHENOL = [("C", 1.39 * math.cos(math.radians(a)), 1.39 * math.sin(math.radians(a)), 0.0)
          for a in range(0, 360, 60)] + [("O", 2.75, 0.0, 0.0)]


@pytest.mark.parametrize("name, skeleton, y, x, angle", [
    ("methanol", [("C", 0, 0, 0), ("O", 1.43, 0, 0)], 1, 0, 108.5),
    ("CH3SH", [("C", 0, 0, 0), ("S", 1.82, 0, 0)], 1, 0, 96.5),
    ("phenol", PHENOL, 6, 0, 108.5),
])
def test_one_h_on_a_two_valent_atom_is_bent(name, skeleton, y, x, angle):
    atoms = _atoms(fix_hydrogens(_xyz(skeleton))["xyz"])
    hs = _hydrogens_of(atoms, y)
    assert len(hs) == 1, name
    assert _angle(atoms[x][1], atoms[y][1], atoms[hs[0]][1]) == pytest.approx(angle, abs=0.5)


def test_ethanol_ch2_is_tetrahedral_and_clear_of_the_oxygen():
    skeleton = [("C", 0, 0, 0), ("C", 1.52, 0, 0), ("O", 2.05, 1.35, 0)]
    out = fix_hydrogens(_xyz(skeleton))
    atoms = _atoms(out["xyz"])
    ch2 = [j for j in _hydrogens_of(atoms, 1, cutoff=1.2)]
    assert len(ch2) == 2
    a, b = (atoms[j][1] for j in ch2)
    assert _angle(a, atoms[1][1], b) == pytest.approx(109.47, abs=1.0)
    # one on each side of the C-C-O plane (z = 0)
    assert a[2] * b[2] < 0
    # no H closer than 1.8 A to a heavy atom it is not bonded to
    for j, (s, p) in enumerate(atoms):
        if s != "H":
            continue
        for k, (t, q) in enumerate(atoms):
            if t != "H" and math.dist(p, q) > 1.2:
                assert math.dist(p, q) > 1.8
    assert not any("WARNING" in c for c in out["checks"])
