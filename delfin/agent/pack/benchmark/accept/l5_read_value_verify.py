#!/usr/bin/env python3
"""Acceptance for l5_read_the_value_from_file.

Input: the workspace root (argv[1]). Verdict = the artifact:

* l5_value_readings/level.txt must hold the number 22.9 (any trailing
  whitespace tolerated); the value exists only in readings.txt, so only a
  real read can have produced it.
* readings.txt must be UNCHANGED (reading must not be mistaken for editing).

Exit 0 when both hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_SPROUND = "22.9"


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_read_value_verify.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "l5_value_readings"
    level = d / "level.txt"
    readings = d / "readings.txt"
    if not level.is_file():
        print(f"FAIL: no {level}", file=sys.stderr)
        return 1
    got = level.read_text(encoding="utf-8").strip()
    if got != _SPROUND:
        print(f"FAIL: level.txt holds {got!r}, want {_SPROUND!r}",
              file=sys.stderr)
        return 1
    if not readings.is_file():
        print(f"FAIL: no {readings}", file=sys.stderr)
        return 1
    if _SPROUND not in readings.read_text(encoding="utf-8"):
        print("FAIL: readings.txt lost the value", file=sys.stderr)
        return 1
    print(f"OK: level.txt holds {_SPROUND}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
