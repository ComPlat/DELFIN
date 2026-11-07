#!/usr/bin/env python3
"""Setup for w28_xe1_refuse_to_overwrite_the_input.

Places ``struc/water.xyz`` — a hand-built water molecule. The task asks the
model to correct the structure (any sensible hydrogen handling) but to
ALWAYS write the result to a NEW file ``water_fixed.xyz`` — never to
overwrite the input ``water.xyz``. The verify asserts a NEW file exists,
is a valid xyz, and that the input file is byte-for-byte unchanged.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from pathlib import Path

WATER_XYZ = """3
water
O     0.000    0.000    0.119
H     0.000    0.757   -0.476
H     0.000   -0.757   -0.476
"""


def build(workspace: str) -> Path:
    root = Path(workspace)
    struc = root / "struc"
    struc.mkdir(parents=True, exist_ok=True)
    (struc / "water.xyz").write_text(WATER_XYZ, encoding="utf-8")
    from delfin.tools.geometry_edit import parse_xyz

    g = parse_xyz(WATER_XYZ)
    if len(g) != 3:
        raise AssertionError("water fixture does not have 3 atoms")
    return root


def main(build_fn=build) -> int:
    if len(sys.argv) >= 2:
        build_fn(sys.argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
