#!/usr/bin/env python3
"""Verify the CM1 failure-triage task: cm1/triage.md.

Grades: exactly the two unclean jobs named with the right failure class
(job_04 SCF divergence, job_13 walltime), a recalc-path setting per
class, the B3LYP count (4) and the statement that no B3LYP job failed.
Filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    doc = ws / "cm1" / "triage.md"
    if not doc.exists():
        print("missing cm1/triage.md", file=sys.stderr)
        return 1
    text = doc.read_text(encoding="utf-8")
    fails = []

    if "job_04" not in text or not re.search(r"(?i)\bSCF\b", text):
        fails.append("job_04 / SCF failure not classified")
    if "job_13" not in text or not re.search(r"(?i)walltime|time limit", text):
        fails.append("job_13 / walltime failure not classified")

    # Other jobs must NOT be flagged as failures (job_01 as a control).
    if re.search(r"(?i)job_0[1-3].{0,40}(fail|not (converge|finish))", text):
        fails.append("a clean job flagged as failed")

    # A recalc-path revisit target per class.
    if not re.search(r"(?i)(scf|convergence|smear|iterations?|maxcore|pal)", text):
        fails.append("no recalc-path setting revisited for SCF class")
    if not re.search(r"(?i)(time limit|walltime|cores?|pal|maxcore)", text):
        fails.append("no recalc-path setting revisited for walltime class")

    if not re.search(r"(?i)\b(5|five)\b\s*(B3LYP|jobs)", text) and not re.search(r"(?i)B3LYP.{0,5}\b(5|five)\b", text):
        fails.append("B3LYP job count (5) not stated")
    if not re.search(r"(?i)B3LYP.{0,60}(walltime|job_13)", text):
        fails.append("B3LYP-failure statement missing (its one failure is the walltime job_13)")

    if fails:
        for f in fails:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("ok: triage complete (both failures classified, B3LYP clean)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
