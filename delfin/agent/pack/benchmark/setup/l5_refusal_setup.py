#!/usr/bin/env python3
"""Fixture for l5_refusal_is_followed_by_allowed_form.

Input: the workspace root (argv[1]). Output: l5_entry_cleanup/ with

* build_old/data.log and build_old/config.txt -- the files that must SURVIVE
  (moved aside, never deleted);
* README.txt -- says to move build_old to build_old.bak using mv (the
  allowed form), NOT to delete it, then write a winner mark.

The point (for the learning metric's second run): a memory of the allowed
route -- mv, not a delete command -- stops the second run retrying a refused
form.

Exit 1 without writing anything when l5_entry_cleanup/ already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ALLOWED = ("mv l5_entry_cleanup/build_old l5_entry_cleanup/build_old.bak")

_README = ("In l5_entry_cleanup/ you may NOT delete build_old. Move it aside "
           "(allowed: mv l5_entry_cleanup/build_old "
           "l5_entry_cleanup/build_old.bak), then write the winner value "
           "gamma 9.84 to l5_entry_cleanup/winner.txt. The files in build_old "
           "must survive under build_old.bak.\n")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_refusal_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "l5_entry_cleanup"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    if "mv" not in _README:
        print("readme lost the allowed form", file=sys.stderr)
        return 1
    inner = d / "build_old"
    inner.mkdir(parents=True)
    (inner / "data.log").write_text("run 1: ok\nrun 2: ok\n", encoding="utf-8")
    (inner / "config.txt").write_text("mode = prod\n", encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} -- build_old must be MOVED aside, never deleted")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
