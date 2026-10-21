#!/usr/bin/env python3
"""Fixture for l5_read_the_value_from_file.

Input: the workspace root (argv[1]). Output: l5_value_readings/ with

* readings.txt -- the value 22.9 m; the prompt does not state it, so only a
  real read can learn it.
* README.txt   -- names readings.txt and asks for the number in level.txt.

A learning comparison (l5) runs this SAME task twice; a skill/memory that
recalls the exact path + unit lets the second run skip the read-and-guess
calls, which is what learning_verdict measures.

Exit 1 without writing anything when l5_value_readings/ already exists or a
precondition fails: readings.txt holds the value.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_VALUE = "22.9"
_UNIT = "m"
_READINGS = f"level: {_VALUE} {_UNIT}\nrecorded at 2026-10-12T09:00Z on probe L5-A\n"
_README = ("l5_value_readings/readings.txt holds the latest level. Read it "
           "with your reading tools, report the value (with its unit), and "
           "write just the number to l5_value_readings/level.txt.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_read_value_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "l5_value_readings"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    if _VALUE not in _READINGS:
        print("readings lost the value", file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "readings.txt").write_text(_READINGS, encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} -- readings.txt holds {_VALUE} {_UNIT}; the prompt "
          "does not state it, so only a real read learns it")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
