#!/usr/bin/env python3
"""Setup for the four wave-25 MN1 fixtures (package MN1, nacht-s9).

Builds four sibling fixture folders under <workspace>/mn1_w25/:

* wait_line/    -- run_slow_job.sh: sleeps a FIXED 15 s, appends to
  job.log ("PHASE ONE ...", after 7 s "PHASE TWO REACHED <reason>"),
  so the agent must start it in the background and wait on a LINE.
* appear_file/  -- make_report.py: sleeps a FIXED 10 s, then writes
  report.json with fixed numbers (n=5, mean_kj_per_mol=41.73).
* fail_job/     -- flaky_job.sh: sleeps 5 s, writes progress lines,
  then fails with a fixed reason line ("FATAL: disk quota exceeded
  while writing checkpoint.bin") and exit 3. The agent must notice,
  NOT retry, and report the reason.
* two_jobs/     -- fast_job.sh (2 s -> fast.log) and slow_job.sh
  (12 s -> slow.log); the agent must start both and act on whichever
  finishes first.

Every sleep is a fixed constant written into the fixture: no machine
speed is probed and no test depends on wall-clock timing beyond the
task's max_duration_s. Exit non-zero on any failure: a task whose
precondition was not built is reported as unmeasured, never as a model
failure. With no argument it prints usage and exits non-zero (the
rule the benchmark harness tests setups by).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT_NAME = "mn1_w25"

WAIT_LINE_SH = """#!/bin/sh
# Fake long job: two phases, the second one is the event to wait for.
sleep 7
echo "PHASE ONE started at cycle 1" >> job.log
sleep 8
echo "PHASE TWO REACHED: the optimizer left the first basin" >> job.log
"""

APPEAR_PY = '''"""Fake report job: sleeps 10 s, then writes report.json."""
import json
import time

time.sleep(10)
with open("report.json", "w", encoding="utf-8") as fh:
    json.dump({"n": 5, "mean_kj_per_mol": 41.73}, fh)
print("report written")
'''

FLAKY_SH = """#!/bin/sh
# Fails half-way with a fixed, findable reason.
echo "progress 1/3" >> job.log
sleep 5
echo "progress 2/3" >> job.log
echo "FATAL: disk quota exceeded while writing checkpoint.bin" >> job.log
echo "job aborted, exit 3" >> job.log
exit 3
"""

FAST_SH = """#!/bin/sh
# Finishes first (2 s).
sleep 2
echo "fast job done: 3 iterations, residual 1.2e-6" >> fast.log
"""

SLOW_SH = """#!/bin/sh
# Finishes last (12 s).
sleep 12
echo "slow job done: full sweep" >> slow.log
"""


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def build(workspace: Path) -> None:
    root = workspace / ROOT_NAME
    wait = root / "wait_line"
    wait.mkdir(parents=True, exist_ok=True)
    _write(wait / "run_slow_job.sh", WAIT_LINE_SH)
    (wait / "run_slow_job.sh").chmod(0o755)

    appear = root / "appear_file"
    appear.mkdir(parents=True, exist_ok=True)
    _write(appear / "make_report.py", APPEAR_PY)

    fail = root / "fail_job"
    fail.mkdir(parents=True, exist_ok=True)
    _write(fail / "flaky_job.sh", FLAKY_SH)
    (fail / "flaky_job.sh").chmod(0o755)

    two = root / "two_jobs"
    two.mkdir(parents=True, exist_ok=True)
    _write(two / "fast_job.sh", FAST_SH)
    _write(two / "slow_job.sh", SLOW_SH)
    (two / "fast_job.sh").chmod(0o755)
    (two / "slow_job.sh").chmod(0o755)

    for folder, need in (
            (wait, ["run_slow_job.sh"]),
            (appear, ["make_report.py"]),
            (fail, ["flaky_job.sh"]),
            (two, ["fast_job.sh", "slow_job.sh"])):
        missing = [n for n in need if not (folder / n).is_file()]
        if missing:
            raise SystemExit(f"setup incomplete: {folder} misses {missing}")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: mn1_w25_setup.py <workspace>")
        return 2
    build(Path(sys.argv[1]))
    print("mn1_w25 fixtures built")
    return 0


if __name__ == "__main__":
    sys.exit(main())
