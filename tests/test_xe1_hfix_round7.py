"""XE1 round-7: aromatic rings (benzene skeleton), pyridine N, phenol
C-O, amide N planarity. Red controls for the operator's round-6 rerun
findings; every earlier real-run case is regression-pinned."""
from __future__ import annotations

import math

from delfin.structure_edit import fix_hydrogens


def _parse(out_xyz):
    atoms = []
    for ln in out_xyz.splitlines():
        p = ln.split()
        if len(p) == 4:
            try:
                atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
            except ValueError:
                pass
    return atoms


def _angle(a, b, c):
    u = [b[k] - a[k] for k in range(3)]
    v = [c[k] - a[k] for k in range(3)]
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    return math.degrees(math.acos(
        max(-1.0, min(1.0, sum(u[k] * v[k] for k in range(3)) / (nu * nv)))
    ))


def _h_of(atoms, vertex):
    return [i for i, (s, p) in enumerate(atoms)
            if s == "H" and math.dist(p, atoms[vertex][1]) < 1.3]


BENZENE = (
    "6\n\nC  1.39  0.0  0.0\nC  0.695  1.204  0.0\nC  -0.695  1.204  0.0\n"
    "C  -1.39  0.0  0.0\nC  -0.695  -1.204  0.0\nC  0.695  -1.204  0.0\n"
)
# Pyridine: N replaces the bottom-left ring C (C-N ~1.34 aromatic).
PYRIDINE = (
    "6\n\nC  1.39  0.0  0.0\nC  0.695  1.204  0.0\nC  -0.695  1.204  0.0\n"
    "C  -1.39  0.0  0.0\nN  -0.695  -1.204  0.0\nC  0.695  -1.204  0.0\n"
)
# Phenol-like: ring C1 carries O at 1.36 A pointing outward from the ring
# (real phenol topology — the C bonded to O keeps NO H, missing = 4 - 1.5
# - 1.5 - 1 = 0).
PHENOL = (
    "7\n\nC  1.39  0.0  0.0\nC  0.695  1.204  0.0\nC  -0.695  1.204  0.0\n"
    "C  -1.39  0.0  0.0\nC  -0.695  -1.204  0.0\nC  0.695  -1.204  0.0\n"
    "O  1.375  2.382  0.0\n"
)
# Formamide skeleton: C with O (1.23) and N (1.34, next to the C=O).
AMIDE = (
    "3\n\nC  0.0  0.0  0.0\nO  1.23  0.0  0.0\nN  -0.68  1.17  0.0\n"
)


def test_benzene_gets_exactly_six_h_in_plane():
    out = fix_hydrogens(BENZENE)
    atoms = _parse(out["xyz"])
    n_c = sum(1 for s, _ in atoms if s == "C")
    n_h = sum(1 for s, _ in atoms if s == "H")
    assert n_c == 6 and n_h == 6, f"got {n_c} C, {n_h} H"
    for ci in range(n_c):
        hs = _h_of(atoms, ci)
        assert len(hs) == 1, f"C{ci} carries {len(hs)} H"
        # In-plane: H z matches its C's z.
        assert abs(atoms[hs[0]][1][2] - atoms[ci][1][2]) < 0.2
    # No H-H below 1.4 A anywhere.
    hh = [math.dist(atoms[a][1], atoms[b][1])
          for a in range(len(atoms)) for b in range(a + 1, len(atoms))
          if atoms[a][0] == "H" and atoms[b][0] == "H"]
    assert min(hh) >= 1.4, f"min H-H {min(hh):.2f}"


def test_pyridine_n_carries_no_h():
    out = fix_hydrogens(PYRIDINE)
    atoms = _parse(out["xyz"])
    ni = [i for i, (s, _) in enumerate(atoms) if s == "N"][0]
    assert _h_of(atoms, ni) == [], "pyridine N must carry no H"
    assert sum(1 for s, _ in atoms if s == "H") == 5


def test_phenol_ring_c_bonded_to_o_gets_no_h():
    out = fix_hydrogens(PHENOL)
    atoms = _parse(out["xyz"])
    oi = [i for i, (s, _) in enumerate(atoms) if s == "O"][0]
    # The C bonded to O carries NO H (real phenol: valence saturated by
    # the two aromatic bonds + the C-O single bond). Its ring neighbours
    # keep exactly one H each.
    ci = min(range(len(atoms)),
             key=lambda i: math.dist(atoms[i][1], atoms[oi][1])
             if atoms[i][0] == "C" else 99.0)
    assert len(_h_of(atoms, ci)) == 0, "C bonded to O must carry no H"
    oH = sum(1 for s, _ in atoms if s == "H")
    # 6 ring C: 5 with one H each (C-O carbon none) + the O gets 1 H.
    assert oH == 6, f"expected 6 H (5 ring + 1 OH), got {oH}"


def test_amide_n_is_planar():
    out = fix_hydrogens(AMIDE)
    atoms = _parse(out["xyz"])
    ni = [i for i, (s, _) in enumerate(atoms) if s == "N"][0]
    hs = _h_of(atoms, ni)
    assert len(hs) == 2
    a, b = hs
    ang = _angle(atoms[ni][1], atoms[a][1], atoms[b][1])
    assert abs(ang - 120.0) < 6.0, f"H-N-H {ang:.1f} (planar amide)"


def test_regression_all_previous_cases():
    # Every earlier real-run case, pinned (operator round-5 demand).
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nC  1.20  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    for c in range(2):
        hs = _h_of(atoms, c)
        assert len(hs) == 1
        assert abs(_angle(atoms[c][1], atoms[1 - c][1],
                          atoms[hs[0]][1]) - 180.0) < 3.0
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nC  1.33  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    assert sum(1 for s, _ in atoms if s == "H") == 4
    out = fix_hydrogens("1\n\nO  0.0  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    hs = _h_of(atoms, 0)
    assert abs(_angle(atoms[0][1], atoms[hs[0]][1],
                      atoms[hs[1]][1]) - 104.5) < 3.0
