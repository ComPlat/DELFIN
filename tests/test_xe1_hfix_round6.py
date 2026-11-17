"""XE1 round-6: sp carbon (acetylene) and the H direction at a carbon
with several heavy neighbours (formamide). Red controls for the
operator's round-5 rerun findings; every earlier real-run case is
regression-pinned so one fix cannot break another again."""
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


def _h_of(atoms, vertex):
    return [i for i, (s, p) in enumerate(atoms)
            if s == "H" and math.dist(p, atoms[vertex][1]) < 1.3]


def test_acetylene_sp_h_is_linear_not_trans_bent():
    # REGRESSION (round 4 had it linear; round 5 broke it): C#C 1.20 A
    # is a triple bond -> sp -> X-C-H 180, NOT the 120-deg imine ring.
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nC  1.20  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    c0, c1 = 0, 1
    for c, other in ((c0, c1), (c1, c0)):
        hs = _h_of(atoms, c)
        assert len(hs) == 1, f"C{c} carries {len(hs)} H"
        ang = _angle(atoms[other][1], atoms[c][1], atoms[hs[0]][1])
        assert abs(ang - 180.0) < 3.0, f"C-C-H {ang:.1f}"


def test_formamide_c_h_points_away_from_the_heavy_neighbours():
    # C with 2+ heavy neighbours: H along -(sum of unit vectors to the
    # neighbours), so it points AWAY from them (operator round 5: the H
    # sat 0.265 A from N).
    xyz = ("3\n\nC  0.0  0.0  0.0\nO  1.23  0.0  0.0\n"
           "N  -0.68  1.17  0.0\n")
    out = fix_hydrogens(xyz)
    atoms = _parse(out["xyz"])
    c = 0
    hs = _h_of(atoms, c)
    assert len(hs) == 1, f"C carries {len(hs)} H"
    # No clash with the neighbours: distance H..O and H..N well above
    # the covalent-contact range.
    for other in (1, 2):
        d = math.dist(atoms[hs[0]][1], atoms[other][1])
        assert d > 2.0, f"H..{atoms[other][0]} {d:.3f} A (clash)"
    # C-H still a proper bond length.
    assert abs(math.dist(atoms[hs[0]][1], atoms[c][1]) - 1.09) < 0.03


def test_regression_all_previous_realrun_cases():
    # Methanol: 1 H on O, 3 on C, no coincident H.
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nO  1.42  0.0  0.0\n")
    a = _parse(out["xyz"])
    assert len(_h_of(a, 1)) == 1 and len(_h_of(a, 0)) == 3
    # HCN: 1 H on C, 0 on N.
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nN  1.16  0.0  0.0\n")
    a = _parse(out["xyz"])
    assert len(_h_of(a, 0)) == 1 and len(_h_of(a, 1)) == 0
    # Ethylene: 2 H per C, planar, ~120 deg.
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nC  1.33  0.0  0.0\n")
    a = _parse(out["xyz"])
    h0 = _h_of(a, 0)
    assert len(h0) == 2
    for h in h0:
        assert abs(_angle(a[1][1], a[0][1], a[h][1]) - 120.0) < 3.0
    # Methylamine: pyramidal NH2.
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nN  1.47  0.0  0.0\n")
    a = _parse(out["xyz"])
    n = next(i for i, (s, _) in enumerate(a) if s == "N")
    hn = _h_of(a, n)
    assert len(hn) == 2
    assert abs(_angle(a[hn[0]][1], a[n][1], a[hn[1]][1]) - 106.7) < 4.0
    # Lone CH4: C-H 1.09, H-C-H 109.47.
    out = fix_hydrogens("1\n\nC  0.0  0.0  0.0\n")
    a = _parse(out["xyz"])
    hs = [i for i, (s, _) in enumerate(a) if s == "H"]
    assert len(hs) == 4
    assert abs(_angle(a[hs[0]][1], a[0][1], a[hs[1]][1]) - 109.47) < 2.0
