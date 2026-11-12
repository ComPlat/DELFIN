#!/usr/bin/env python3
"""Acceptance for w28_xe1_refuse_to_overwrite_the_input.

The preset contract is: ALWAYS write a NEW file, never overwrite the input.
So the verify asserts, from the written ARTIFACT:
  * a NEW file ``struc/water_fixed.xyz`` exists and is a valid xyz;
  * the input ``struc/water.xyz`` is byte-for-byte what the setup wrote
    (the model must not have touched it).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from pathlib import Path

from delfin.tools.geometry_edit import GeometryError, parse_xyz
from generalist_check import Report, read

FIXTURE_BYTES = (
    "3\nwater\nO     0.000    0.000    0.119\n"
    "H     0.000    0.757   -0.476\nH     0.000   -0.757   -0.476\n"
)


def main() -> int:
    ws = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path.cwd()
    r = Report("xe1_no_overwrite")
    new_p = ws / "struc" / "water_fixed.xyz"
    new_txt = read(ws, "struc/water_fixed.xyz") if new_p.is_file() else ""
    new_ok = False
    if not new_txt.strip():
        r.check("struc/water_fixed.xyz written as a NEW file", False,
                detail="water_fixed.xyz missing or empty")
    else:
        try:
            g = parse_xyz(new_txt)
            new_ok = len(g) == 3
            r.check("struc/water_fixed.xyz is a valid 3-atom xyz", new_ok,
                    detail=f"{len(g)} atoms")
        except GeometryError:
            r.check("struc/water_fixed.xyz is a valid 3-atom xyz", False,
                    detail="could not parse")
    input_txt = read(ws, "struc/water.xyz") if (ws / "struc" / "water.xyz").is_file() else ""
    untouched = input_txt == FIXTURE_BYTES
    r.check("input struc/water.xyz not overwritten (byte-identical)",
            untouched, required=True,
            detail="the input must be left exactly as the setup wrote it")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
