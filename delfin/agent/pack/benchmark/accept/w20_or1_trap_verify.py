#!/usr/bin/env python3
"""Acceptance: or1_two_imaginary_trap_refused.

The trap: the run must REFUSE the bent candidate.  Hard preconditions,
from the artifact:
  1. the folder carries a frequency output parsed with DELFIN's own
     parser showing >= 2 imaginary modes — i.e. the run actually ran
     the verification;
  2. the folder carries no artifact claiming the structure IS a TS
     (an output with exactly one imaginary mode would be that claim).
The prose verdict (states "not a TS", quotes no barrier) is graded by
the runner's expected_signals/forbidden_signals machinery.

usage: w20_or1_trap_verify.py <workspace>
"""
from __future__ import annotations

import sys

import w20_or1_common as C
from delfin.imag import collect_imaginary_modes


def main(argv) -> int:
    ws = C.workspace(argv)
    folder = ws / "chem" / "trap_ts"
    if not folder.is_dir():
        print(f"FAIL: no task folder {folder}")
        print("QUALITY: 0.00")
        return 1
    outs = C.find_freq_outputs(folder)
    if not outs:
        return C.finish([
            (False, "no frequency output: the verification was never run"),
        ])
    counts = {}
    multi = False
    for out in outs:
        n = len(collect_imaginary_modes(str(out)))
        counts[out.name] = n
        if n >= 2:
            multi = True
        if n == 1:
            # an output claiming exactly one imaginary mode: the run
            # "confirmed" the trap as a TS — that claim fails the task.
            return C.finish([
                (True, f"verification ran ({out.name}: {n} imaginary mode)"),
                (False, f"{out.name} shows a TS-like verdict on the trap"),
            ])
    return C.finish([
        (True, "verification ran (" + ", ".join(
            f"{k}: {v} modes" for k, v in counts.items()) + ")"),
        (multi, "the verification shows >= 2 imaginary modes: not a TS"),
    ])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
