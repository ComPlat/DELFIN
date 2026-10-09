#!/usr/bin/env python3
"""Fixture for the Wave-14 U5 benchmark task: a tool call the model writes
as text (a leaked DSML <invoke> block) must be re-issued as a real call.

w14_u5_readings/ holds two files after setup:

* readings.txt — the only place the real value (water level 41.7 m) lives.
  The task PROMPT does not state the value, so a model that leaks a
  <invoke>...cat readings.txt...</invoke> block as text (the wave-13 DSML
  leak) never actually reads it and therefore never learns 41.7: it either
  stands still or fabricates. Only a REAL read_file/grep_file/bash call gets
  the value.
* README.txt — a short line naming the file to read and what to write to
  level.txt (folder-handling instructions only; no value).

After the turn the accept script (accept/w14_u5_leak_verify.py) checks the
artifact: level.txt must carry the EXACT real value 41.7 (proving a real
read+write happened — text never runs) and readings.txt must be byte
identical (the source was not overwritten to smuggle the value in).

Precondition (checked BEFORE writing anything): the notes below still state
the value so the real part is answerable via a genuine read; the fixture
must NOT contain DSML markup, so a forbidden-signal hit is always the
model's leak, never a fixture echo.
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
