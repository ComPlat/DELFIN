#!/usr/bin/env python3
"""Acceptance: w28_tf1_c_reject_mixed_levels (ts-finder task C).

The seeded record (tf1/c/levels.txt from setup/tf1_c_seed.py) assigns the TS a
basis (def2-TZVP) different from the other stages (def2-SVP), so the barrier is
NOT at one level of theory and must be rejected. The grader parses levels.txt
into the per-stage dict delfin.ts_finder.check_single_level_of_theory expects
and requires ok=False, and requires verdict.json to say "mixed levels of
theory".

usage: tf1_c_verify.py <workspace>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from delfin.ts_finder import (
    check_single_level_of_theory,
    CheckResult,
)


def _levels_dict(text: str) -> dict:
    levels = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        key, _, value = line.partition(" ")
        if key:
            levels[key] = value
    return levels


def main(argv) -> int:
    if len(argv) != 2:
        print("usage: tf1_c_verify.py <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1])
    folder = root / "tf1" / "c"
    levels_path = folder / "levels.txt"
    verdict_path = folder / "verdict.json"

    preconditions = []

    # 1. The expert wrote an honest "mixed levels" verdict.
    if verdict_path.is_file():
        try:
            verdict = json.loads(verdict_path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            verdict = None
        label = verdict.get("verdict") if isinstance(verdict, dict) else None
        preconditions.append(
            (label == "mixed levels of theory",
             "verdict.json says 'mixed levels of theory'"))
    else:
        preconditions.append((False, "verdict.json exists"))

    # 2. Re-derive from the record: one-level check must reject.
    if levels_path.is_file():
        r = check_single_level_of_theory(
            _levels_dict(levels_path.read_text(encoding="utf-8")))
        preconditions.append(
            (isinstance(r, CheckResult) and not r.ok,
             "the real one-level check REJECTS the mixed record (ok=False)"))
    else:
        preconditions.append((False, "levels.txt exists"))

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
