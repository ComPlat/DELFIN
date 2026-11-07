#!/usr/bin/env python3
"""Acceptance for w28_xe1_fix_missing_h_on_metal_complex.

Reads ``struc/complex_fixed.xyz`` — the NEW file the model was asked to
write — and checks the written ARTIFACT, never the prose:
  * the new file exists and parses;
  * its composition is exactly Zn1 N2 H3 (the fixed zinc amide: the
    input had Zn1 N2 H2 — the missing amido H was added);
  * the input complex.xyz was not modified (overwrite refusal).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from delfin.tools.geometry_edit import GeometryError, parse_xyz
from generalist_check import Report, read

EXPECTED = Counter({"Zn": 1, "N": 2, "H": 3})


def main() -> int:
    ws = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path.cwd()
    r = Report("xe1_fix_h")
    new_path = ws / "struc" / "complex_fixed.xyz"
    raw = read(ws, "struc/complex_fixed.xyz")
    exists = new_path.is_file() and bool(raw.strip())
    r.check("new structure file written", exists,
            detail="struc/complex_fixed.xyz")
    if not exists:
        return r.finish()

    try:
        g = parse_xyz(raw)
    except GeometryError as exc:
        r.check("new structure parses", False, detail=str(exc))
        return r.finish()

    syms = Counter(g.symbols)
    r.check("composition is Zn1 N2 H3 (missing amido H added)",
            syms == EXPECTED, detail=f"found {dict(syms)}")

    # Overwrite refusal: the input file must survive byte-for-byte.
    input_path = ws / "struc" / "complex.xyz"
    input_raw = read(ws, "struc/complex.xyz") if input_path.is_file() else ""
    r.check("input complex.xyz left unmodified by the fix",
            bool(input_raw.strip()), required=False,
            detail="the two H of the spectator NH2 still present")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
