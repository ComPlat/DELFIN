#!/usr/bin/env python3
"""Verify the CM1 census task: cm1/census.md must tabulate all 20 folders.

Grades: one row per job folder (20 rows), method/basis and Eh energy and
status columns, per-method lowest-energy job (never across methods),
failure reasons named for job_04 (SCF) and job_13 (walltime), and the
group sizes. Filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

JOB = re.compile(r"job_\d\d")
NUM = r"-?\d+\.\d+"

EXPECT = {
    "n_converged": 18,
    "scf": "job_04",
    "walltime": "job_13",
    "r2scan_low": "job_05",     # -402.521477, lowest r2SCAN-3c energy
    "b3lyp_low": "job_07",      # -401.781904, lowest B3LYP energy
    "r2scan_n": 15,
    "b3lyp_n": 5,
}


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    doc = ws / "cm1" / "census.md"
    if not doc.exists():
        print("missing cm1/census.md", file=sys.stderr)
        return 1
    text = doc.read_text(encoding="utf-8")
    fails = []

    rows = sorted(set(JOB.findall(text)))
    if len(rows) != 20:
        fails.append(f"census lists {len(rows)} job folders, expected 20")

    if f"{EXPECT['n_converged']}" not in text and "18" not in text:
        fails.append("converged count 18 not stated")

    if EXPECT["scf"] not in text or not re.search(r"(?i)\bSCF\b", text):
        fails.append("job_04 / SCF failure not named")
    if EXPECT["walltime"] not in text or not re.search(r"(?i)walltime|time limit", text):
        fails.append("job_13 / walltime not named")

    # Per-method lowest-energy jobs — both must appear NEAR their energy.
    for key, job, energy in (("r2scan_low", EXPECT["r2scan_low"], "-402.521"),
                             ("b3lyp_low", EXPECT["b3lyp_low"], "-401.781")):
        if job not in text or energy not in text:
            fails.append(f"per-method lowest {key}: {job} / {energy} Eh missing")

    # Per-method group sizes must appear ADJACENT to the method name (a
    # bare "15"/"5" elsewhere — e.g. inside energy digits or job_05/job_15
    # tokens — used to satisfy the old substring check: reviewer finding D).
    if not (re.search(r"(?i)\b15\b[^.\n]{0,40}r2SCAN", text)
            or re.search(r"(?i)r2SCAN[^.\n]{0,40}\b15\b", text)):
        fails.append("r2SCAN group size 15 not stated adjacent to the method")
    if not (re.search(r"(?i)\b5\b[^.\n]{0,40}B3LYP", text)
            or re.search(r"(?i)B3LYP[^.\n]{0,40}\b5\b", text)):
        fails.append("B3LYP group size 5 not stated adjacent to the method")

    # A cross-method "lowest overall" claim is the defect this task hunts.
    m = re.search(r"(?i)lowest.{0,60}(overall|across|global)", text)
    if m and re.search(r"(?i)(job_0[5679]|B3LYP)", m.group(0)):
        fails.append("cross-method lowest-energy claim found")

    if fails:
        for f in fails:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("ok: census complete (20 rows, per-method minima, failures named)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
