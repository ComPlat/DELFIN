#!/usr/bin/env python3
"""Acceptance for l5_read_the_value_from_file.

Input: the workspace root (argv[1]). Verdict = the artifact:

* l5_value_readings/level.txt holds exactly one number, optionally
  followed by the unit ``m``, within 0.01 of 22.9 (``22.9``, ``22.90`` and
  ``22.9 m`` all pass); the value exists only in readings.txt.
* readings.txt is byte-for-byte what the setup wrote (reading is not
  editing).

Exit 0 when both hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_EXPECTED = 22.9
_TOLERANCE = 0.01
# Must equal l5_read_value_setup._READINGS.
_READINGS = "level: 22.9 m\nrecorded at 2026-10-12T09:00Z on probe L5-A\n"
_NUMBER_WITH_UNIT = re.compile(r"\s*([-+]?\d+(?:\.\d+)?)\s*(?:m)?\s*")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_read_value_verify.py <workspace>", file=sys.stderr)
        return 2
    d = Path(argv[1]).resolve() / "l5_value_readings"
    level = d / "level.txt"
    readings = d / "readings.txt"
    if not level.is_file():
        print(f"FAIL: no {level}", file=sys.stderr)
        return 1
    got = level.read_text(encoding="utf-8")
    m = _NUMBER_WITH_UNIT.fullmatch(got)
    if not m:
        print(f"FAIL: level.txt holds {got.strip()!r}, not one number "
              "(optionally with the unit m)", file=sys.stderr)
        return 1
    value = float(m.group(1))
    if abs(value - _EXPECTED) > _TOLERANCE:
        print(f"FAIL: level.txt holds {value}, want {_EXPECTED} "
              f"+/- {_TOLERANCE}", file=sys.stderr)
        return 1
    if (not readings.is_file()
            or readings.read_text(encoding="utf-8") != _READINGS):
        print("FAIL: readings.txt was changed or removed", file=sys.stderr)
        return 1
    print(f"OK: level.txt holds {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
