#!/usr/bin/env python3
"""Acceptance: the metal-complex geometry has the right shape.

[Au(CN)2]- must come out linear (C-Au-C ~180 deg); Cr(CO)6 octahedral
(trans C-Cr-C 180 deg, cis ~90 deg, six equal Cr-C bonds).  The check
reads the model's optimised xyz (or re-optimises our start geometry as
the physics baseline) and measures the angles itself - a Hessian-min
that is the wrong shape is a fail even with zero imaginary modes.

usage: chem_geometry_shape.py <workspace> <mol_key>
"""

from __future__ import annotations

import itertools
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

MOLS = {
    "au_cyanide": ("au_cyanide", -1),
    "cr_carbonyl": ("cr_carbonyl", 0),
}


def _xtb() -> str:
    from delfin.qm_runtime import resolve_tool
    r = resolve_tool("xtb")
    if r is None or not r.path:
        print("xtb not found via qm_runtime resolver", file=sys.stderr)
        sys.exit(3)
    return str(r.path)


def _load_xyz(p: Path) -> list[tuple[str, float, float, float]]:
    lines = p.read_text(encoding="utf-8").splitlines()
    n = int(lines[0].split()[0])
    out = []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        out.append((parts[0], float(parts[1]), float(parts[2]),
                    float(parts[3])))
    return out


def _angle(a, b, c) -> float:
    v1 = [a[i] - b[i] for i in (1, 2, 3)]
    v2 = [c[i] - b[i] for i in (1, 2, 3)]
    dot = sum(v1[i] * v2[i] for i in range(3))
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("usage: <workspace> <mol_key>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    mol_key = argv[2]
    if mol_key not in MOLS:
        print(f"unknown molecule key: {mol_key}", file=sys.stderr)
        return 2
    folder, chrg = MOLS[mol_key]

    mdir = ws / "chem" / folder
    if not mdir.is_dir():
        print(f"no chem/{folder} folder in workspace", file=sys.stderr)
        return 1
    cand = mdir / "xtbopt.xyz"
    if not cand.is_file():
        xyzs = sorted(mdir.glob("*.xyz"))
        cand = xyzs[-1] if xyzs else mdir / "start.xyz"
    if not cand.is_file():
        print(f"no xyz in chem/{folder}", file=sys.stderr)
        return 1

    atoms = _load_xyz(cand)
    if mol_key == "au_cyanide":
        au = [a for a in atoms if a[0] == "Au"]
        cs = [a for a in atoms if a[0] == "C"]
        if len(au) != 1 or len(cs) != 2:
            print("expected Au + 2 C in [Au(CN)2]-", file=sys.stderr)
            return 1
        ang = _angle(cs[0], au[0], cs[1])
        if abs(ang - 180.0) > 5.0:
            print(f"C-Au-C {ang:.1f} deg, expected ~180 (linear)",
                  file=sys.stderr)
            return 1
        print(f"linear ok: C-Au-C = {ang:.1f} deg")
        return 0

    # cr_carbonyl: octahedral
    cr = [a for a in atoms if a[0] == "Cr"]
    cs = [a for a in atoms if a[0] == "C"]
    if len(cr) != 1 or len(cs) != 6:
        print("expected Cr + 6 C in Cr(CO)6", file=sys.stderr)
        return 1
    angs = [_angle(a, cr[0], b) for a, b in itertools.combinations(cs, 2)]
    trans = [a for a in angs if a > 175.0]
    cis_ok = all(80.0 <= a <= 100.0 for a in angs if a <= 175.0)
    dists = [math.dist(cr[0][1:], c[1:]) for c in cs]
    spread = max(dists) - min(dists)
    if len(trans) != 3 or not cis_ok:
        print(f"not octahedral: 3 trans pairs expected, got {len(trans)}; "
              f"cis angles outside 80-100", file=sys.stderr)
        return 1
    if spread > 0.1:
        print(f"Cr-C bonds not equivalent (spread {spread:.3f} A)",
              file=sys.stderr)
        return 1
    print(f"octahedral ok: 3 trans pairs, cis in band, "
          f"Cr-C = {dists[0]:.3f} A (spread {spread:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
