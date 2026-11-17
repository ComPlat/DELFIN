#!/usr/bin/env python3
"""Fixture for w25-sl1-03 (failure classification): the recorded
failure text, host-neutral.

Writes <workspace>/sl1/failure/ with the CONTROL.txt resource snippet
and the job-log excerpts the prompt quotes. Refuses to overwrite an
existing sl1/ tree. The numbers are the recorded wave-21/22 values
(TIMEOUT at the 48 h limit; the sacct TotalCPU artefact on walltime
kills is operator-documented). No scheduler is contacted.
"""
from __future__ import annotations

import sys
from pathlib import Path

CONTROL_SNIPPET = """\
# CONTROL.txt resource lines of the failed run (excerpt)
PAL = 8
maxcore = 4000
time_limit = 24:00:00
"""

ERR_TAIL = """\
slurmstepd: error: *** JOB 7924102 CANCELLED AT 2026-11-01T03:12 \
DUE TO TIME LIMIT ***
"""

OUT_TAIL = """\
...  (optimisation cycle 214, energy converged at step 211)
TOTAL ENERGY SUMMARY — cycle 214 continues, no error termination
"""


def build(ws: Path) -> None:
    root = ws / "sl1" / "failure"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        raise SystemExit(1)
    root.mkdir(parents=True)
    (root / "CONTROL_snippet.txt").write_text(CONTROL_SNIPPET,
                                              encoding="utf-8")
    (root / "job.err").write_text(ERR_TAIL, encoding="utf-8")
    (root / "job.out").write_text(OUT_TAIL, encoding="utf-8")
    (root / "README.txt").write_text(
        "Recorded failure fixture, host-neutral: walltime kill at the "
        "48:00:00 limit (time_limit was 24:00:00 — the backend doubles "
        "nothing; the 48 h is the operator's account default). sacct: "
        "State TIMEOUT, Elapsed 48:00:00, Timelimit 48:00:00, TotalCPU "
        "00:00:00 (artefact on walltime kills — the job computed).\n",
        encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    build(Path(argv[1]).resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
