#!/usr/bin/env python3
"""Acceptance: or1_hcn_hnc_ts_and_irc.

Hard preconditions, from the artifact the run left behind:
  1. a frequency output in chem/hcn_hnc shows EXACTLY one imaginary
     mode (the TS verification);
  2. the prose states a barrier in kcal/mol inside the cited band
     45-53 (HCN<->HNC isomerisation, literature band, grading only).
Exit 0 = both hold; stdout is the evidence.

usage: w20_or1_ts_hcn_verify.py <workspace>
"""
from __future__ import annotations

import sys
from pathlib import Path

import w20_or1_common as C


def main(argv) -> int:
    ws = C.workspace(argv)
    folder = ws / "chem" / "hcn_hnc"
    if not folder.is_dir():
        print(f"FAIL: no task folder {folder}")
        print("QUALITY: 0.00")
        return 1
    ok_mode, ev = C.one_imaginary_mode(folder)
    prose = ""
    for md in folder.rglob("*.md"):
        prose += md.read_text(encoding="utf-8", errors="ignore")
    # F-QS-1: the band is graded HERE on the artifact (a computed
    # barrier statement over the freq output), not left to the runner's
    # permissive expected_values judge, which a cited band also passes.
    ok_band = C.barrier_with_artifact(folder, prose, 45, 53)
    return C.finish([
        (ok_mode, f"TS carries exactly one imaginary mode ({ev})"),
        (ok_band, "computed barrier stated in 45-53 kcal/mol over the artifact"),
    ])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
