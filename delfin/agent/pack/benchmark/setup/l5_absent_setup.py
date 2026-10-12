#!/usr/bin/env python3
"""Fixture for l5_report_absent_quantity_as_absent.

Input: the workspace root (argv[1]). Output: l5_absent_catalog/ with

* catalog.txt -- two entries. Entry A has a freeze_point of -5.5 degC;
  entry B has NO freeze_point recorded (a learned skill: an unavailable
  figure is reported "absent", not guessed as zero).
* README.txt   -- asks for each entry's freeze_point with its unit.

The prompt does not state the numbers or which entry lacks a value, so a
real read is required. Correct answer: report -5.5 degC for A and say B has
no freeze_point recorded -- do not invent a figure for B.

Exit 1 without writing anything when l5_absent_catalog/ already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_CATALOG = (
    "entry A: propanol ; freeze_point: -5.5 degC\n"
    "entry B: cyclohexanol ; freeze_point:\n"
)
_README = ("l5_absent_catalog/catalog.txt lists two entries. For each, "
           "report the freeze_point with its unit. If an entry has no "
           "freeze_point recorded, say so explicitly.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_absent_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "l5_absent_catalog"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    if "-5.5" not in _CATALOG:
        print("catalog lost the known freeze_point", file=sys.stderr)
        return 1
    # Entry B must genuinely have NO value -- that is the whole point.
    if "freeze_point:" not in _CATALOG or "freeze_point: -5.5" not in _CATALOG:
        print("catalog shape unexpected", file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "catalog.txt").write_text(_CATALOG, encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} -- A has freeze_point, B has none (absent)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
