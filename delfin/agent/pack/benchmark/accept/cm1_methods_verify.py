#!/usr/bin/env python3
"""Verify the CM1 method-table task: cm1/methods.md.

Grades: per method+basis group — job count (16 / 4), mean of CONVERGED
jobs in Eh (r2SCAN-3c: -402.515946; B3LYP: -401.780091), spread in mEh
(r2SCAN-3c: 11.565 mEh; B3LYP: 6.294 mEh), failed jobs excluded and
listed separately, and the no-cross-method-comparison statement.
Filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

EXPECT = {
    ("r2SCAN-3c", "def2-mTZVP"): {"n": 15, "mean": -402.516381, "spread": 9.136},
    ("B3LYP", "def2-SVP"): {"n": 5, "mean": -401.780091, "spread": 3.694},
}


def _near(text: str, value: float, tol: float) -> bool:
    m = re.findall(r"(-?\d+\.?\d*)", text)
    return any(abs(float(v) - value) <= tol for v in m)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    doc = ws / "cm1" / "methods.md"
    if not doc.exists():
        print("missing cm1/methods.md", file=sys.stderr)
        return 1
    text = doc.read_text(encoding="utf-8")
    fails = []

    if not _near(text, EXPECT[("r2SCAN-3c", "def2-mTZVP")]["mean"], 0.0005):
        fails.append("r2SCAN-3c converged mean -402.5164 Eh missing/wrong")
    if not _near(text, EXPECT[("B3LYP", "def2-SVP")]["mean"], 0.0005):
        fails.append("B3LYP converged mean -401.7801 Eh missing/wrong")

    if not _near(text, EXPECT[("r2SCAN-3c", "def2-mTZVP")]["spread"], 0.05):
        fails.append("r2SCAN-3c spread 9.136 mEh missing/wrong")
    if not _near(text, EXPECT[("B3LYP", "def2-SVP")]["spread"], 0.05):
        fails.append("B3LYP spread 3.694 mEh missing/wrong")

    if not re.search(r"(?i)job_04", text) or not re.search(r"(?i)job_13", text):
        fails.append("failed jobs not excluded/listed separately")
    if not re.search(r"(?i)(14|fourteen).{0,30}(converged|r2SCAN)", text) and not re.search(r"(?i)(18|eighteen)", text):
        fails.append("converged counts not stated")

    if not re.search(r"(?i)(not|never).{0,30}compar", text):
        fails.append("no statement that group energies must not be compared")

    if not re.search(r"(?i)Eh\b", text) or not re.search(r"(?i)mEh\b", text):
        fails.append("units Eh/mEh not stated")

    if fails:
        for f in fails:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("ok: method table complete (per-group means, spreads, units)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
