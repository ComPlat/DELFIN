#!/usr/bin/env python3
"""Acceptance for w28_xe1_renumber_product_onto_educt.

Reads ``struc/product_renumbered.xyz`` — the NEW file the model was asked
to write — and checks the written ARTIFACT, never the prose:
  * the new file exists and is a valid xyz;
  * its composition is exactly one C, one N, one O;
  * its element ORDER recovers the educt's reference order [C, N, O]
    (a scrambled copy would keep [O, N, C] and fail here);
  * the input product.xyz was not changed (overwrite refusal).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from delfin.tools.geometry_edit import GeometryError, parse_xyz
from generalist_check import Report, read

EDUCT_ORDER = ["C", "N", "O"]


def main() -> int:
    ws = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path.cwd()
    r = Report("xe1_renumber")
    out = ws / "struc" / "product_renumbered.xyz"
    if not out.is_file():
        r.check("written result file", False,
                detail="struc/product_renumbered.xyz not found "
                       "(the model must write a NEW file, never overwrite "
                       "the input)")
        return r.finish()

    try:
        geom = parse_xyz((ws / "struc" / "product_renumbered.xyz").read_text(
            encoding="utf-8"))
    except GeometryError as exc:
        r.check("result is a valid xyz", False, detail=f"parse failed: {exc}")
        return r.finish()
    r.check("result is a valid xyz", True,
            detail=f"{len(geom)} atoms parsed")

    r.check("composition is C1 N1 O1",
            Counter(geom.symbols) == Counter({"C": 1, "N": 1, "O": 1}),
            detail=f"got {dict(Counter(geom.symbols))}")

    r.check("element order recovers the educt [C, N, O]",
            list(geom.symbols) == EDUCT_ORDER,
            detail=f"got {list(geom.symbols)}; a scrambled copy keeps "
                   f"[O, N, C] and fails")

    # The input must be untouched (the structure dept never overwrites).
    r.check("input product.xyz unchanged",
            read(ws, "struc/product.xyz").startswith(
                "3\nproduct (scrambled O,N,C order)\n"),
            detail="product.xyz was modified")

    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
