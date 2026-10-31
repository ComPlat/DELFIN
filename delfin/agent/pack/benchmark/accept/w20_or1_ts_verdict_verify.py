#!/usr/bin/env python3
"""Acceptance: or1_ts_candidate_verdict.

Hard precondition: a frequency output under chem/ts_candidate was
produced by running the verification (any mode count is acceptable —
the task is the VERDICT, and the prose is graded by the runner).

usage: w20_or1_ts_verdict_verify.py <workspace>
"""
from __future__ import annotations

import sys

import w20_or1_common as C


def main(argv) -> int:
    ws = C.workspace(argv)
    folder = ws / "chem" / "ts_candidate"
    if not folder.is_dir():
        print(f"FAIL: no task folder {folder}")
        print("QUALITY: 0.00")
        return 1
    outs = C.find_freq_outputs(folder)
    return C.finish([
        (bool(outs), "a frequency output exists (the verification ran)"),
    ])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
