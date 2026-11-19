#!/usr/bin/env python3
"""Setup for w28_xe1_renumber_product_onto_educt.

Places in ``struc/`` a renumberable pair of hand-built xyz structures:
``educt.xyz`` (reference atom order) and ``product.xyz`` (the same atoms
in a scrambled order). Both are the same three carbon atoms of a bent C-C-C
skeleton; only the line ordering differs, so the correct answer (reorder
the product onto the educt mapping) is a permutation that the verify script
can assert against the educt's element order plus positions.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

# educt: bent C-N-O triatomic, atom order C(0), N(1), O(2).
EDUCT_XYZ = """3
educt C-N-O reference order
C     0.000   0.000   0.000
N     1.200   0.000   0.000
O     2.000   0.400   0.000
"""

# product: the SAME three atoms, listed in a scrambled order O, N, C
# (same positions, only the element/atom ordering changed).
PRODUCT_XYZ = """3
product (scrambled O,N,C order)
O     2.000   0.400   0.000
N     1.200   0.000   0.000
C     0.000   0.000   0.000
"""


def build(workspace: str) -> Path:
    root = Path(workspace)
    struc = root / "struc"
    if struc.exists():
        print(f"refusing to overwrite {struc}", file=sys.stderr)
        sys.exit(1)
    struc.mkdir(parents=True, exist_ok=True)
    (struc / "educt.xyz").write_text(EDUCT_XYZ, encoding="utf-8")
    (struc / "product.xyz").write_text(PRODUCT_XYZ, encoding="utf-8")
    # Precondition: both fixtures are parseable as xyz and share a formula.
    from delfin.tools.geometry_edit import parse_xyz

    a = parse_xyz(EDUCT_XYZ)
    b = parse_xyz(PRODUCT_XYZ)
    if Counter(a.symbols) != Counter(b.symbols):
        raise AssertionError("renumber fixtures do not share a formula")
    return root


if __name__ == "__main__":
    from generalist_repo import main

    main(build)
