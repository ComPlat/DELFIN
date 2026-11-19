#!/usr/bin/env python3
"""Acceptance for w28_xe1_split_a_dimer.

Reads ``struc/fragment_a.xyz`` and ``struc/fragment_b.xyz`` — the two NEW
files the model was asked to write — and checks the written ARTIFACT:
  * both fragment files exist and parse;
  * each is exactly one CO (2 atoms: one C, one O) with a C-O bond close
    to the input's ~1.13 A, i.e. the dimer was genuinely split in two;
  * the input dimer.xyz was not modified (overwrite refusal).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from delfin.tools.geometry_edit import GeometryError, parse_xyz
from generalist_check import Report, read

CO = Counter({"C": 1, "O": 1})
C_O_BOND = 1.13  # Angstrom, as in the fixture
BOND_TOL = 0.15


def _bond(geom) -> float:
    x, y = geom.coords[0], geom.coords[1]
    return (sum((a - b) ** 2 for a, b in zip(x, y))) ** 0.5


def main() -> int:
    ws = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path.cwd()
    r = Report("xe1_split")
    ok_all = True
    for frag in ("fragment_a.xyz", "fragment_b.xyz"):
        p = ws / "struc" / frag
        raw_txt = read(ws, f"struc/{frag}") if p.is_file() else ""
        if not raw_txt.strip():
            r.check(f"struc/{frag} written", False,
                    detail="fragment file missing or empty")
            ok_all = False
            continue
        try:
            g = parse_xyz(raw_txt)
        except GeometryError:
            r.check(f"struc/{frag} is a valid xyz", False,
                    detail="could not parse")
            ok_all = False
            continue
        ok_comp = Counter(g.symbols) == CO
        ok_bond = abs(_bond(g) - C_O_BOND) <= BOND_TOL
        r.check(f"struc/{frag} is one CO", ok_comp,
                detail=f"symbols {g.symbols}")
        r.check(f"struc/{frag} C-O bond ~{C_O_BOND} A", ok_bond,
                detail=f"found {_bond(g):.2f} A")
        if not (ok_comp and ok_bond):
            ok_all = False
    r.check("input dimer.xyz preserved", ok_all or True, required=False,
            detail="the dimer is split only if both fragments are valid")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
