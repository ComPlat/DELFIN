#!/usr/bin/env python3
"""Acceptance: or1_sn2_identity_barrier.

Hard precondition, from the artifact: a frequency output under
chem/sn2 with EXACTLY one imaginary mode (the symmetric TS); the
barrier band (15 +/- 4 kcal, literature 12-20) is graded on the prose
by the runner's expected_values machinery.

usage: w20_or1_ts_sn2_verify.py <workspace>
"""
from __future__ import annotations

import sys

import w20_or1_common as C


def main(argv) -> int:
    ws = C.workspace(argv)
    folder = ws / "chem" / "sn2"
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
    ok_band = C.barrier_with_artifact(folder, prose, 12, 20)
    return C.finish([
        (ok_mode, f"SN2 TS carries exactly one imaginary mode ({ev})"),
        (ok_band, "computed barrier stated in 12-20 kcal/mol over the artifact"),
    ])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
