#!/usr/bin/env python3
"""A2 of the capability census, through ``bash``.

Question: ``b2/`` holds three instrument runs as CSV. What is the mean
cycle time of each run, and which run is fastest?

Why the answer differs when the capability works: the means are not in
any file. They must be computed. The three runs have 5, 5 and 3 rows, so
the pooled mean (116.00 s) equals NONE of the three per-run means
(120.20, 90.80, 151.00) -- an answer that averages everything gives a
number that is not an answer to any of the questions asked. ``run_c``
also carries a ``warmup_seconds`` column, so summing every column named
"seconds"-ish inflates it. A reasoning agent that cannot run code guesses
or reports the raw list; one with a working shell states the real means.

Exit 0 = pass.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

MEANS = {"a": 120.20, "b": 90.80, "c": 151.00}
POOLED = 116.00
TOL = 0.02


def _numbers(text: str) -> list[float]:
    return [float(m.group()) for m in
            re.finditer(r"\d+(?:\.\d+)?", text)]


def _states(nums: list[float], want: float) -> bool:
    return any(abs(n - want) <= TOL for n in nums)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    answer = ws / "b2" / "answer.md"
    if not answer.is_file():
        print("missing b2/answer.md", file=sys.stderr)
        return 1
    text = answer.read_text(encoding="utf-8", errors="replace")
    nums = _numbers(text)
    fails: list[str] = []

    for run, mean in MEANS.items():
        if not _states(nums, mean):
            fails.append(f"the mean cycle time of run_{run} ({mean:.2f} s) "
                         "is not stated with the right value")

    # The pooled mean is the signature of not separating the runs. It is
    # only a defect when it is offered AS a run's mean, so the check is
    # that all three per-run means are present (above) and that 116.00 is
    # not passed off as one of them.
    if _states(nums, POOLED) and not all(_states(nums, m) for m in MEANS.values()):
        fails.append(
            f"the answer states {POOLED:.2f} s, the pooled mean over all 13 "
            "rows; none of the three per-run means is that number")

    # "b" is fastest. Naming it is required; naming a or c is wrong.
    fast = re.search(r"(?i)\brun[_\s-]?b\b", text)
    if not fast:
        fails.append("the fastest run (run_b) is not named")
    for wrong in ("a", "c"):
        if re.search(rf"(?i)run[_\s-]?{wrong}\b\s*(is|was)?\s*(the\s*)?"
                     rf"(fastest|quickest)", text):
            fails.append(f"run_{wrong} is called the fastest; run_b is")
    if re.search(r"(?i)(fastest|quickest)[^.\n]{0,40}\brun[_\s-]?([ac])\b", text):
        fails.append("a run other than run_b is presented as the fastest")

    if fails:
        print("FAIL: " + "; ".join(fails), file=sys.stderr)
        return 1
    print("OK: all three per-run means stated and run_b identified as "
          "fastest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
