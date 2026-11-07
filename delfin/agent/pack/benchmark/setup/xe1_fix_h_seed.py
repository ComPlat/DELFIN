#!/usr/bin/env python3
"""Setup for w28_xe1_fix_missing_h_on_metal_complex.

Places in ``struc/`` a hand-built zinc amide ``complex.xyz`` in which the
amido nitrogen N1 carries only ONE hydrogen (an amido N must carry two,
N-H ~1.01 A, Zn-N ~1.95 A). The correct action is to add the missing H
(geometry_fix_hydrogens) and write the corrected structure to a NEW file
``complex_fixed.xyz``. The verify script asserts the new file has exactly
one more H than the input, the same heavy atoms, and that N1 is now
satisfied (two H).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Zn -- N1(amido) with only ONE H on N1 (missing H5), plus a spectator
# NH2 group N2 with both H. Hand-built with Zn-N 1.95, N-H 1.01, H-N-H ~109.
COMPLEX_XYZ = """6
zinc amide, N1 missing one H
Zn    0.000    0.000    0.000
N    -1.950    0.000    0.000
H    -2.430    0.700    0.000
N     1.950    0.000    0.000
H     2.430    0.700    0.000
H     2.430   -0.700    0.000
"""


def build(workspace: str) -> Path:
    root = Path(workspace)
    struc = root / "struc"
    struc.mkdir(parents=True, exist_ok=True)
    (struc / "complex.xyz").write_text(COMPLEX_XYZ, encoding="utf-8")
    from delfin.tools.geometry_edit import parse_xyz

    g = parse_xyz(COMPLEX_XYZ)
    if len(g) != 6:  # Zn, N, H, N, H, H
        raise AssertionError("fix-H fixture does not have 6 atoms")
    return root


def main(build_fn=build) -> int:
    if len(sys.argv) >= 2:
        build_fn(sys.argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
