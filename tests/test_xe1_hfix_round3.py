"""XE1 real-run round-3 controls for fix_hydrogens (operator FAIL report,
head 3ffde73ae): ethylene sp2 hydrogens must point away from the C=C bond
at ~121.5 deg (not 70.5 deg inward), a lone O must give a water-like
H-O-H angle (~104.5 deg, not 141 deg), min H-H distance >= 1.4 A, no
placed X-H may point INTO a heavy-heavy bond, and a blank input comment
line must not make the output comment repeat the first atom line.
"""
from __future__ import annotations

import math

from delfin.structure_edit import fix_hydrogens


def _parse(xyz: str):
    lines = xyz.strip().splitlines()
    n = int(lines[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append(tuple(float(p) for p in parts[1:4]))
    return syms, coords


ETHYLENE = """2

C  0.000  0.000  0.000
C  1.330  0.000  0.000
"""


def _angle(p1, p2, p3):
    """Angle at p2 in degrees."""
    import math as m

    def sub(a, b):
        return [a[i] - b[i] for i in range(3)]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def norm(a):
        return m.sqrt(dot(a, a)) or 1e-12

    u, v = sub(p1, p2), sub(p3, p2)
    return m.degrees(
        m.acos(max(-1.0, min(1.0, dot(u, v) / norm(u) / norm(v))))
    )


def test_ethylene_sp2_h_away_from_partner():
    syms, coords = _parse(fix_hydrogens(ETHYLENE)["xyz"])
    assert sorted(syms) == ["C", "C", "H", "H", "H", "H"]
    # each C gets 2 H; C-C-H angle ~121.5 deg (sp2), never < 100 deg
    for ci, cj in ((0, 1), (1, 0)):
        hs = [k for k, s in enumerate(syms) if s == "H" and
              math.dist(coords[k], coords[ci]) < 1.2]
        assert len(hs) == 2, f"C{ci} carries {len(hs)} H"
        for h in hs:
            ang = _angle(coords[cj], coords[ci], coords[h])
            assert ang > 100.0, f"C-C-H angle {ang:.1f} deg (want ~121.5)"


def test_lone_o_water_angle():
    O = "1\n\nO  0.0  0.0  0.0\n"
    out = fix_hydrogens(O)
    syms, coords = _parse(out["xyz"])
    assert sorted(syms) == ["H", "H", "O"]
    hs = [coords[i] for i, s in enumerate(syms) if s == "H"]
    ang = _angle(hs[0], coords[0], hs[1])
    assert 95.0 <= ang <= 115.0, f"H-O-H {ang:.1f} deg (want ~104.5)"


def test_min_hh_distance():
    for xyz in (ETHYLENE, "1\n\nO  0.0  0.0  0.0\n"):
        syms, coords = _parse(fix_hydrogens(xyz)["xyz"])
        hs = [c for c, s in zip(coords, syms) if s == "H"]
        for a in range(len(hs)):
            for b in range(a + 1, len(hs)):
                d = math.dist(hs[a], hs[b])
                assert d >= 1.4, f"H-H {d:.3f} A < 1.4"


def test_no_h_pointing_into_a_heavy_bond():
    syms, coords = _parse(fix_hydrogens(ETHYLENE)["xyz"])
    for k, s in enumerate(syms):
        if s != "H":
            continue
        # bonded heavy atom
        ci = min(
            (i for i, ss in enumerate(syms) if ss == "C"),
            key=lambda i: math.dist(coords[k], coords[i]),
        )
        d_h_c = math.dist(coords[k], coords[ci])
        cj = [i for i, ss in enumerate(syms) if ss == "C" and i != ci][0]
        # angle H-C(=partner) must be obtuse-ish: H must not sit between
        # the two carbons (that would be < 90 deg with a short H..Cj)
        d_h_cj = math.dist(coords[k], coords[cj])
        ang = _angle(coords[ci], coords[k], coords[cj])
        assert not (ang < 90.0 and d_h_cj < 1.8), (
            f"H points into the C=C bond: angle {ang:.1f}, d(H..Cj) {d_h_cj:.2f}"
        )


def test_blank_comment_not_repeated():
    xyz = "1\n\nO  0.0  0.0  0.0\n"
    out = fix_hydrogens(xyz)
    lines = out["xyz"].strip().splitlines()
    # comment line (index 1) must not be an atom line
    parts = lines[1].split()
    assert not (parts and parts[0] in ("O", "H", "C", "N")), (
        f"comment line repeats an atom line: {lines[1]!r}"
    )
