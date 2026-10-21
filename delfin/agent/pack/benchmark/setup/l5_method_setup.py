#!/usr/bin/env python3
"""Fixture for l5_name_method_with_basis_pairing.

Input: the workspace root (argv[1]). Output: l5_method_notes/ with

* brief.txt -- asks for ONE recommended DFT-functional/basis-set combination
  for a closed-shell H,C,N,O molecule, no solvent, no dispersion, in ORCA.

The answer is the file l5_method_notes/choice.txt, judged by
accept/l5_method_verify.py: one functional, one basis set, no dispersion
correction, and an ORCA keyword line naming both.

Exit 1 without writing anything when l5_method_notes/ already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_BRIEF = ("Recommend ONE density-functional/basis-set combination for a "
          "closed-shell molecule with only H, C, N and O. No solvent, no "
          "dispersion correction. Name the functional, the basis set, and "
          "the keyword line you would use in ORCA, and write them to "
          "l5_method_notes/choice.txt.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_method_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "l5_method_notes"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    if "functional" not in _BRIEF:
        print("brief lost the functional", file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "brief.txt").write_text(_BRIEF, encoding="utf-8")
    print(f"prepared {d} -- one recommended functional/basis pair")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
