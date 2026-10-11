#!/usr/bin/env python3
"""Fixture for w14_u5_dsml_leak_is_reissued_as_a_real_call.

Input: the workspace root (argv[1]). Output: w14_u5_readings/ with

* readings.txt -- the water level 41.7 m; the prompt does not state it,
  so a tool call written as text (which runs nothing) cannot learn it.
* README.txt   -- names readings.txt and asks for the number in level.txt.

Exit 1 without writing anything when w14_u5_readings/ already exists or a
precondition fails: readings.txt holds the value, and no fixture text
contains '<', so any <invoke> markup in a trajectory is the model's.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_VALUE = "41.7"
_UNIT = "m"
_READINGS = "water level: 41.7 m\nrecorded at 2026-10-09T08:00Z on probe A-7\n"
_README = ("w14_u5_readings/readings.txt holds the latest water level. Read "
           "it with your reading tools, report the value (with its unit), "
           "and write just the number to w14_u5_readings/level.txt.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w14_u5_leak_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "w14_u5_readings"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    # Precondition, checked BEFORE anything is written: the value is present
    # so a genuine read answers the real part; and the fixture carries no
    # DSML markup so any <invoke> in the trajectory is the model's leak.
    if _VALUE not in _READINGS:
        print("readings lost the water level", file=sys.stderr)
        return 1
    if "<" in _READINGS or "<" in _README:
        print("fixture must not contain markup that could echo a leak",
              file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "readings.txt").write_text(_READINGS, encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} — readings.txt holds the water level (41.7 m); "
          "the prompt does not state it, so only a real read learns it")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
