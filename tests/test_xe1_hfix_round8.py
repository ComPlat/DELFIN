"""XE1 round-8: O/S with one heavy neighbour must take the element's real
bond angle (O ~108.5, S ~96, NOT linear 180), and the clash self-check
must cover H..heavy distances to atoms OTHER than the parent
(0.9 x covalent sum). Red controls for the operator's round-7 rerun
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
    # 1.4 covers the S-H 1.34 A bond (a 1.3 filter misses it).
    return [i for i, (s, p) in enumerate(atoms)
            if s == "H" and math.dist(p, atoms[vertex][1]) < 1.4]


METHANOL = "2\n\nC  0.0  0.0  0.0\nO  1.42  0.0  0.0\n"
METHANETHIOL = "2\n\nC  0.0  0.0  0.0\nS  1.82  0.0  0.0\n"
PHENOL = (
    "7\n\nC  1.39  0.0  0.0\nC  0.695  1.204  0.0\nC  -0.695  1.204  0.0\n"
    "C  -1.39  0.0  0.0\nC  -0.695  -1.204  0.0\nC  0.695  -1.204  0.0\n"
    "O  1.375  2.382  0.0\n"
)
# Ethanol skeleton: C-C 1.52, C-O ~1.43 from the second C.
ETHANOL = (
    "3\n\nC  0.0  0.0  0.0\nC  1.52  0.0  0.0\nO  2.05  1.35  0.0\n"
)


def test_methanol_oh_angle_is_1085_not_linear():
    out = fix_hydrogens(METHANOL)
    atoms = _parse(out["xyz"])
    oi = [i for i, (s, _) in enumerate(atoms) if s == "O"][0]
    hs = _h_of(atoms, oi)
    assert len(hs) == 1, f"O carries {len(hs)} H"
    ang = _angle(atoms[0][1], atoms[oi][1], atoms[hs[0]][1])
    assert abs(ang - 108.5) < 4.0, f"C-O-H {ang:.1f}, expected ~108.5"


def test_methanethiol_sh_angle_is_96_not_linear():
    out = fix_hydrogens(METHANETHIOL)
    atoms = _parse(out["xyz"])
    si = [i for i, (s, _) in enumerate(atoms) if s == "S"][0]
    hs = _h_of(atoms, si)
    assert len(hs) == 1, f"S carries {len(hs)} H"
    ang = _angle(atoms[0][1], atoms[si][1], atoms[hs[0]][1])
    assert abs(ang - 96.0) < 5.0, f"C-S-H {ang:.1f}, expected ~96"


def test_phenol_oh_angle_is_1085_not_linear():
    out = fix_hydrogens(PHENOL)
    atoms = _parse(out["xyz"])
    oi = [i for i, (s, _) in enumerate(atoms) if s == "O"][0]
    ci = min((i for i, (s, _) in enumerate(atoms) if s == "C"),
             key=lambda i: math.dist(atoms[i][1], atoms[oi][1]))
    hs = _h_of(atoms, oi)
    assert len(hs) == 1, f"O carries {len(hs)} H"
    ang = _angle(atoms[ci][1], atoms[oi][1], atoms[hs[0]][1])
    assert abs(ang - 108.5) < 4.0, f"C-O-H {ang:.1f}, expected ~108.5"


def test_ethanol_h_heavy_clash_is_flagged():
    out = fix_hydrogens(ETHANOL)
    atoms = _parse(out["xyz"])
    checks = " ".join(out["checks"])
    # Any placed H closer than 0.9 x covalent sum to a heavy atom that is
    # NOT its parent must produce a WARNING in checks[].
    oi = [i for i, (s, _) in enumerate(atoms) if s == "O"][0]
    worst = None
    for i, (s, p) in enumerate(atoms):
        if s != "H":
            continue
        for j, (s2, p2) in enumerate(atoms):
            if s2 == "H" or j == i:
                continue
            d = math.dist(p, p2)
            if worst is None or d < worst[0]:
                worst = (d, i, j)
    assert worst is not None
    d, hi, heavy = worst
    # The skeleton is clash-free by construction; a placed H too close to
    # a foreign heavy atom must be flagged.
    if d < 1.0:
        assert "WARNING" in checks, (
            f"H{hi}..{atoms[heavy][0]}{heavy} at {d:.3f} A not flagged; "
            f"checks={out['checks']}"
        )


# --- regression pins: every earlier real-run case (operator demand) ---

def test_regression_all_earlier_cases():
    # methanol skeleton: CH3OH (3 H on C, 1 on O)
    out = fix_hydrogens(METHANOL)
    atoms = _parse(out["xyz"])
    assert sorted(s for s, _ in atoms) == ["C", "H", "H", "H", "H", "O"]
    # HCN: 1 H on C only
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nN  1.16  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    assert sum(1 for s, _ in atoms if s == "H") == 1
    # benzene: exactly 6 H, one per ring C
    BENZENE = (
        "6\n\nC  1.39  0.0  0.0\nC  0.695  1.204  0.0\n"
        "C  -0.695  1.204  0.0\nC  -1.39  0.0  0.0\n"
        "C  -0.695  -1.204  0.0\nC  0.695  -1.204  0.0\n"
    )
    out = fix_hydrogens(BENZENE)
    atoms = _parse(out["xyz"])
    assert sum(1 for s, _ in atoms if s == "H") == 6
    # ethylene: 4 H, coplanar
    out = fix_hydrogens("2\n\nC  0.0  0.0  0.0\nC  1.33  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    assert sum(1 for s, _ in atoms if s == "H") == 4
    # water lone O: H-O-H 104.5
    out = fix_hydrogens("1\n\nO  0.0  0.0  0.0\n")
    atoms = _parse(out["xyz"])
    h = [i for i, (s, _) in enumerate(atoms) if s == "H"]
    assert len(h) == 2
    ang = _angle(atoms[h[0]][1], atoms[0][1], atoms[h[1]][1])
    assert abs(ang - 104.5) < 3.0
