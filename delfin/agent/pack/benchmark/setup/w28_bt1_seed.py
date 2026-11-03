#!/usr/bin/env python3
"""Seed the butane rotation task (wave 28, bt1).

Seeds <workspace>/chem/butane/start.xyz: hand-built anti-periplanar
(staggered) n-butane from z-matrix arithmetic in Angstrom
(C-C 1.526, C-H 1.093, CCC 112.5 deg, HCC ~109.5 deg), phi(C1-C2-C3-C4)
= 180.  No CSD-derived coordinates.  Refuses to overwrite an existing
chem/butane tree, like w20_or1_seed.py does.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path


def _anti_butane() -> str:
    c_cc, c_ch = 1.526, 1.093
    ccc = math.radians(112.5)
    # Chain along +z: C1-C2-C3-C4 with phi(C1-C2-C3-C4)=180 (anti).
    r23 = c_cc * math.sin(ccc)
    dz = c_cc * math.cos(ccc)
    # C1 at origin, C2 along z, C3 in the x-z plane.
    c1 = (0.0, 0.0, 0.0)
    c2 = (0.0, 0.0, c_cc)
    c3 = (r23, 0.0, c_cc + dz)
    # C4: torsion 180 deg about the C2-C3 axis from C1 -> anti position.
    hcc = math.radians(109.5)
    r_ch_xz = c_ch * math.sin(hcc)
    # two methyls on C1 and C3 plus two mid CH2 hydrogens; approximate
    # staggered arrangement, refined by the pre-optimisation the task asks.
    hydrogens = []
    # methyl on C1 (staggered w.r.t. C2-C3 bond direction)
    for phi in (90.0, 210.0, 330.0):
        p = math.radians(phi)
        hydrogens.append((c1[0] + r_ch_xz * math.cos(p),
                          c1[1] + r_ch_xz * math.sin(p),
                          c1[2] - c_ch * math.cos(hcc)))
    # methyl on C4 (anti: continues the zigzag in -x direction)
    c4 = (c3[0] - r23 * math.cos(-ccc), 0.0,
          c3[2] + c_cc * math.cos(ccc))
    # recompute C4 exactly as anti continuation: rotate C1 frame by 180
    c4 = (2 * c3[0] - c2[0], 0.0, 2 * c3[2] - c2[2])
    for phi in (270.0, 30.0, 150.0):
        p = math.radians(phi)
        hydrogens.append((c4[0] + r_ch_xz * math.cos(p),
                          c4[1] + r_ch_xz * math.sin(p),
                          c4[2] + c_ch * math.cos(hcc)))
    # two CH2 hydrogens on C2 and C3 (anti pair per carbon)
    for center, phi_pair, zsign in ((c2, (0.0, 180.0), -1),
                                    (c3, (0.0, 180.0), +1)):
        for phi in phi_pair:
            p = math.radians(phi)
            hydrogens.append((center[0] + r_ch_xz * math.cos(p),
                              center[1] + r_ch_xz * math.sin(p) * zsign,
                              center[2]))
    lines = ["10", "n-butane anti-periplanar start (hand-built)"]
    for el, xyz in ([("C", c1), ("C", c2), ("C", c3), ("C", c4)]
                    + [("H", h) for h in hydrogens]):
        lines.append(f"{el} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f}")
    return "\n".join(lines) + "\n"


README = """\
Butane rotation benchmark start geometry
=======================================
Molecule: n-butane, C4H10, neutral, singlet (charge 0, mult 1).
start.xyz: anti-periplanar (staggered) conformation, phi(C1-C2-C3-C4) = 180.
Atoms 1-4 are the carbons in chain order; the scan dihedral is D 0 1 2 3.
Hand-built from z-matrix arithmetic (C-C 1.526 A, C-H 1.093 A);
pre-optimise (xtb or MLP) before the DFT scan.
"""


def main(argv) -> int:
    if len(argv) < 2:
        print("usage: w28_bt1_seed.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    folder = ws / "chem" / "butane"
    if folder.exists():
        print(f"ERROR: {folder} already exists; refusing to overwrite",
              file=sys.stderr)
        return 1
    folder.mkdir(parents=True)
    (folder / "start.xyz").write_text(_anti_butane(), encoding="utf-8")
    (folder / "README.md").write_text(README, encoding="utf-8")
    print(f"seeded {folder}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
