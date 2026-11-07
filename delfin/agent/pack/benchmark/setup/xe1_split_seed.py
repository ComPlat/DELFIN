#!/usr/bin/env python3
"""Setup for w28_xe1_split_a_dimer.

Places in ``struc/`` a ``dimer.xyz`` of two IDENTICAL CO molecules placed
far apart (C-O bond ~1.13 A, molecules separated by ~7 A, so there is no
covalent bond between them). The correct action is to split the dimer into
its two molecules and write each to a NEW file ``fragment_a.xyz`` and
``fragment_b.xyz`` (each a single CO, 2 atoms). The verify script reads the
two written files and asserts each is one CO of the right geometry.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from pathlib import Path

DIMER_XYZ = """4
dimer: two CO molecules, far apart
C    -3.500    0.000    0.000
O    -4.630    0.000    0.000
C     3.500    0.000    0.000
O     4.630    0.000    0.000
"""


def build(workspace: str) -> Path:
    root = Path(workspace)
    struc = root / "struc"
    struc.mkdir(parents=True, exist_ok=True)
    (struc / "dimer.xyz").write_text(DIMER_XYZ, encoding="utf-8")
    from delfin.tools.geometry_edit import parse_xyz

    g = parse_xyz(DIMER_XYZ)
    if len(g) != 4:  # C O C O
        raise AssertionError("split fixture does not have 4 atoms")
    return root


def main(build_fn=build) -> int:
    if len(sys.argv) >= 2:
        build_fn(sys.argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
