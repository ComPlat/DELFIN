#!/usr/bin/env python3
"""Acceptance: w28_tf1_b_reject_wrong_mode (ts-finder task B).

The seeded candidate (tf1/b/ts.xyz + tf1/b/freq.hess, from setup/tf1_b_seed.py)
has a single imaginary mode that stretches the SPECTATOR C-N bond only (N is
off the H-C axis, so C moving along y does not change the intended H-C length).
The intended reactive pair is H-C (0,1); the mode does NOT act on it, so the
real delfin.ts_finder check returns ok=False and the honest verdict is
"found but wrong mode". The grader re-derives that from the files and requires
verdict.json to agree.

usage: tf1_b_verify.py <workspace>
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
        print("usage: tf1_b_verify.py <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1])
    folder = root / "tf1" / "b"
    hess = folder / "freq.hess"
    xyz = folder / "ts.xyz"
    verdict_path = folder / "verdict.json"

    preconditions = []

    # 1. The expert wrote an honest "wrong mode" verdict.
    if verdict_path.is_file():
        try:
            verdict = json.loads(verdict_path.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            verdict = None
        label = verdict.get("verdict") if isinstance(verdict, dict) else None
        # A wrong imaginary mode must never be called "found and verified".
        preconditions.append(
            (label == "found but wrong mode",
             "verdict.json says 'found but wrong mode'"))
        preconditions.append(
            (label != "found and verified",
             "verdict.json does NOT claim 'found and verified'"))
    else:
        preconditions.append((False, "verdict.json exists"))

    # 2. Re-derive from the files: the mode-check must NOT verify H-C.
    if hess.is_file() and xyz.is_file():
        r = check_imaginary_mode_connects_bonds(str(hess), xyz.read_text(
            encoding="utf-8"), [(0, 1)])
        preconditions.append(
            (isinstance(r, CheckResult) and not r.ok,
             "the real mode-check REJECTS the candidate for H-C (ok=False)"))
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
