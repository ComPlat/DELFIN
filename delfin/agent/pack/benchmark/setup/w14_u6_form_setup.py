#!/usr/bin/env python3
"""Fixture for the Wave-14 U6 benchmark task: a refused command shape
(awk / $( ) / xargs / backtick-pipes / sed -i) is followed by the allowed
form (the reading tools or a plain command), and the real result is
produced.

w14_u6_entries/ holds after setup:

* log.txt — entries of the form ``<item> <value>`` (one per line). The
  highest value (9.84, item `gamma`) is the answer the prompt asks for. The
  prompt does NOT state it, so finding it requires actually reading the
  file; the natural-but-refused shortcut is an awk/sort pipeline.
* README.txt — folder-handling instructions only; no value.

After the turn the accept script (accept/w14_u6_form_verify.py) checks the
artifact: winner.txt must carry the EXACT winning item and value (9.84 /
gamma), proving the agent actually read and compared the entries — a
refused shape never ran, so it could not produce the real artifact — and
log.txt must be byte identical (the source was not overwritten).

Precondition (checked BEFORE writing anything): log.txt still holds the
max (9.84 > every other value) so the real part is answerable; the fixture
must NOT contain any of the refused shapes (awk, $(), xargs, backtick
pipes, sed -i), so a forbidden-signal hit is always the model's own form,
never a fixture echo.
"""
from __future__ import annotations

import sys
from pathlib import Path

_MAX_ITEM = "gamma"
_MAX_VALUE = "9.84"
_LOG = (
    "alpha 3.10\n"
    "beta 8.42\n"
    "gamma 9.84\n"
    "delta 1.25\n"
)
_README = ("w14_u6_entries/log.txt lists items and their measured value, one "
           "per line ('<item> <value>'). Find the entry with the highest "
           "value, and write its item and value to "
           "w14_u6_entries/winner.txt as '<item> <value>'.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w14_u6_form_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "w14_u6_entries"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    # Precondition, checked BEFORE anything is written: the max is present
    # and unique; and the fixture carries none of the refused shapes so any
    # forbidden-signal hit is the model's own form.
    if _MAX_ITEM not in _LOG or _MAX_VALUE not in _LOG:
        print("log lost the winning entry", file=sys.stderr)
        return 1
    for token in ("awk", "$(", "xargs", "`", "sed -i"):
        if token in _LOG or token in _README:
            print("fixture must not contain a refused command shape",
                  file=sys.stderr)
            return 1
    d.mkdir(parents=True)
    (d / "log.txt").write_text(_LOG, encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} — log.txt holds four entries; the highest value is "
          "gamma 9.84 (not in the prompt, so a real read finds it)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
