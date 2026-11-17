"""XE1 round-4 red controls: fix_hydrogens angle correctness measured AT
THE VERTEX atom (H-X-H between placed H, X-X-H between a heavy neighbour,
the vertex and an H), as the operator's real-run round 4 demands.

Skeletons: ethylene C2H4 (sp2: C-C-H ~121.5, H-C-H ~117, coplanar),
lone N -> NH3 (H-N-H ~106.7), lone C -> CH4 (H-C-H 109.47, H-H ~1.78),
and the clash self-check must flag any H-H < 1.4 A.
"""
from __future__ import annotations

import math

import pytest

from delfin.structure_edit import _clash_report, fix_hydrogens

ETHYLENE = "2\n\nC  0.000  0.000  0.000\nC  1.330  0.000  0.000\n"
LONE_N = "1\n\nN  0.000  0.000  0.000\n"
LONE_C = "1\n\nC  0.000  0.000  0.000\n"


def _parse(out_xyz):
    lines = [ln for ln in out_xyz.splitlines() if ln.strip()]
    start = 2 if lines[0].strip().isdigit() else 0
    atoms = []
    for ln in lines[start:]:
        p = ln.split()
        atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
    return atoms


def _angle(p, v1, v2):
    a = [v1[k] - p[k] for k in range(3)]
    b = [v2[k] - p[k] for k in range(3)]
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    return math.degrees(math.acos(
        max(-1.0, min(1.0, sum(a[k] * b[k] for k in range(3)) / (na * nb)))
    ))


def _vec(atoms, i, j):
    return atoms[j][1]


def test_ethylene_sp2_angles_at_the_vertex():
    atoms = _parse(fix_hydrogens(ETHYLENE)["xyz"])
    c0, c1 = 0, 1
    # Attribute H by bond distance (an H belongs to the vertex it sits
    # ~1.09 A from; the <150 deg angular filter alone double-counts the
    # far carbon's H, which subtend 120 deg too).
    h0 = [i for i, (s, p) in enumerate(atoms)
          if s == "H" and math.dist(p, atoms[c0][1]) < 1.3]
    h1 = [i for i, (s, p) in enumerate(atoms)
          if s == "H" and math.dist(p, atoms[c1][1]) < 1.3]
    assert len(h0) == 2 and len(h1) == 2
    # C-C-H at each vertex ~120 (planar sp2 ideal; the 121.5 real value
    # carries small pyramidalisation this guess does not model).
    for c, other, hs in ((c0, c1, h0), (c1, c0, h1)):
        for h in hs:
            ang = _angle(atoms[c][1], atoms[other][1], atoms[h][1])
            assert abs(ang - 120.0) < 3.0, f"C-C-H {ang:.1f}"
    # H-C-H ~120 (not 141.1).
    a, b = h0
    ang = _angle(atoms[c0][1], atoms[a][1], atoms[b][1])
    assert abs(ang - 120.0) < 3.0, f"H-C-H {ang:.1f}"
    # Coplanar: all four H and both C share z ~ 0.
    for i in (c0, c1, *h0, *h1):
        assert abs(atoms[i][1][2]) < 1e-6, f"atom {i} off-plane"


def test_lone_n_gives_nh3_angle_106_7():
    atoms = _parse(fix_hydrogens(LONE_N)["xyz"])
    hs = [i for i, (s, _) in enumerate(atoms) if s == "H"]
    assert len(hs) == 3
    for a in range(3):
        for b in range(a + 1, 3):
            ang = _angle(atoms[0][1], atoms[hs[a]][1], atoms[hs[b]][1])
            assert abs(ang - 106.7) < 2.5, f"H-N-H {ang:.1f}"


def test_lone_c_gives_methane_tetrahedral():
    atoms = _parse(fix_hydrogens(LONE_C)["xyz"])
    hs = [i for i, (s, _) in enumerate(atoms) if s == "H"]
    assert len(hs) == 4
    for a in range(4):
        for b in range(a + 1, 4):
            ang = _angle(atoms[0][1], atoms[hs[a]][1], atoms[hs[b]][1])
            assert abs(ang - 109.47) < 2.5, f"H-C-H {ang:.1f}"
    hh = [math.dist(atoms[hs[a]][1], atoms[hs[b]][1])
          for a in range(4) for b in range(a + 1, 4)]
    assert min(hh) >= 1.4, f"min H-H {min(hh):.2f}"


def test_clash_report_flags_hh_below_1_4():
    flagged = _clash_report(
        ["H", "H"], [(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)]
    )
    assert flagged is not None and "1.20" in flagged


def test_clash_report_silent_for_normal_hh():
    assert _clash_report(
        ["H", "H"], [(0.0, 0.0, 0.0), (1.8, 0.0, 0.0)]
    ) is None
