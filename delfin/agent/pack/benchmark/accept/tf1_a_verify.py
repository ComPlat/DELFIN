#!/usr/bin/env python3
"""Acceptance: w28_tf1_a_pick_correct_ts (ts-finder task A).

The seeded candidate (tf1/a/ts.xyz + tf1/a/freq.hess, from setup/tf1_a_seed.py)
is a genuine transition state: exactly one imaginary mode and that mode
stretches the intended H-transfer pair (0,1)=H-C. The ts-finder expert must
VERIFY it (verdict.json = "found and verified") and record the reactive pairs
[[0, 1]]. The grader re-derives the verdict from the seeded files with the
real delfin.ts_finder checks and requires verdict.json to agree.

usage: tf1_a_verify.py <workspace>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from delfin.ts_finder import (
    check_imaginary_mode_connects_bonds,
    CheckResult,
)


def main(argv) -> int:
    if len(argv) != 2:
        print("usage: tf1_a_verify.py <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1])
    folder = root / "tf1" / "a"
    hess = folder / "freq.hess"
    xyz = folder / "ts.xyz"
    verdict_path = folder / "verdict.json"

    preconditions = []

    # 1. The expert produced a verdict.json with the honest label.
    if verdict_path.is_file():
        try:
            verdict = json.loads(verdict_path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            verdict = None
        label = verdict.get("verdict") if isinstance(verdict, dict) else None
        preconditions.append(
            (label == "found and verified",
             "verdict.json says 'found and verified'"))
        pairs = list(verdict.get("pairs", [])) if isinstance(verdict, dict) else []
        preconditions.append(
            (pairs == [[0, 1]],
             f"verdict.json reactive pairs record H-C as [[0, 1]] (got {pairs!r})"))
    else:
        preconditions.append((False, "verdict.json exists"))

    # 2. Re-derive from the files with the real check.
    if hess.is_file() and xyz.is_file():
        r = check_imaginary_mode_connects_bonds(str(hess), xyz.read_text(
            encoding="utf-8"), [(0, 1)])
        preconditions.append(
            (isinstance(r, CheckResult) and r.ok,
             "the real mode-check VERIFIES the candidate (ok=True)"))
    else:
        preconditions.append((False, "freq.hess and ts.xyz exist"))

    ok_all = True
    met = 0
    for ok, label in preconditions:
        print(("PASS" if ok else "FAIL") + ": " + label)
        ok_all = ok_all and ok
        met += 1 if ok else 0
    print(f"QUALITY: {met / len(preconditions):.2f}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
