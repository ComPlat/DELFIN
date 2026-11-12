"""XE1 round-5: lone-atom X-H lengths and amine N hybridisation.
Red controls for the operator's two round-4 findings, measuring at the
vertex atom."""
from __future__ import annotations

import math

from delfin.structure_edit import fix_hydrogens


def _parse(out_xyz):
    lines = out_xyz.splitlines()
    atoms = []
    for ln in lines:
        p = ln.split()
        if len(p) == 4:
            try:
                atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
            except ValueError:
                pass
    return atoms


def _angle(a, b, c):
    # Both legs FROM the vertex b (the earlier u = b - a mixed an
    # incoming with an outgoing leg and distorted the angle readouts;
    # operator round-4 finding applies here too).
    u = [a[k] - b[k] for k in range(3)]
    v = [c[k] - b[k] for k in range(3)]
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    return math.degrees(math.acos(
        max(-1.0, min(1.0, sum(u[k] * v[k] for k in range(3)) / (nu * nv)))
    ))


def test_lone_c_gets_real_ch_bond_length():
    out = fix_hydrogens("1\n\nC  0.0  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    hs = [p for s, p in atoms if s == "H"]
    assert len(hs) == 4
    for p in hs:
        d = math.dist(p, atoms[0][1])
        assert abs(d - 1.09) < 0.02, f"C-H {d:.3f}"  # was 1.888


def test_lone_n_and_o_lengths():
    out = fix_hydrogens("1\n\nN  0.0  0.0  0.0\n")
    for p in (a for a in _parse(out["xyz"]) if a[0] == "H"):
        assert abs(math.dist(p[1], (0.0, 0.0, 0.0)) - 1.01) < 0.02
    out = fix_hydrogens("1\n\nO  0.0  0.0  0.0\n")
    for p in (a for a in _parse(out["xyz"]) if a[0] == "H"):
        assert abs(math.dist(p[1], (0.0, 0.0, 0.0)) - 0.96) < 0.02


def test_methylamine_amine_n_is_pyramidal_not_planar():
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nN  1.47  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    n = next(i for i, (s, _) in enumerate(atoms) if s == "N")
    c = next(i for i, (s, _) in enumerate(atoms) if s == "C")
    hn = [i for i, (s, p) in enumerate(atoms)
          if s == "H" and math.dist(p, atoms[n][1]) < 1.3]
    assert len(hn) == 2, f"N carries {len(hn)} H"
    cnh = [_angle(atoms[c][1], atoms[n][1], atoms[h][1]) for h in hn]
    hnh = _angle(atoms[hn[0]][1], atoms[n][1], atoms[hn[1]][1])
    # Pyramidal amine: C-N-H ~110, H-N-H ~106 (NOT planar 120/120).
    for a in cnh:
        assert abs(a - 110.0) < 4.0, f"C-N-H {a:.1f}"
    assert abs(hnh - 106.7) < 4.0, f"H-N-H {hnh:.1f}"


def test_imine_n_double_bond_gets_one_h():
    # C=N double bond: N valence 3, order 2 -> 1 free valence -> 1 N-H
    # (real imine R-N=CH-R' carries the N-H); NOT 0.
    out = fix_hydrogens("2\n\nN  0.0  0.0  0.0\nC  1.29  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    n = next(i for i, (s, _) in enumerate(atoms) if s == "N")
    hn = [i for i, (s, p) in enumerate(atoms)
          if s == "H" and math.dist(p, atoms[n][1]) < 1.3]
    assert len(hn) == 1, "C=N imine N carries exactly 1 H"
