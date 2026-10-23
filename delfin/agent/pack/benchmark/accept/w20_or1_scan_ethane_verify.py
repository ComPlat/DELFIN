#!/usr/bin/env python3
"""Acceptance: or1_ethane_methyl_rotation_scan.

Hard precondition, from the artifact: a scan output exists under
chem/ethane whose input carried a %geom constraint/scan block, and the
prose graded by the runner pins the barrier band (2.9 +/- 0.6 kcal).

usage: w20_or1_scan_ethane_verify.py <workspace>
"""
from __future__ import annotations

import sys

import w20_or1_common as C


def main(argv) -> int:
    ws = C.workspace(argv)
    folder = ws / "chem" / "ethane"
    if not folder.is_dir():
        print(f"FAIL: no task folder {folder}")
        print("QUALITY: 0.00")
        return 1
    # The scan artifact: an ORCA input with a scan constraint, and an
    # output beside it.
    inputs = list(folder.rglob("*.inp"))
    scan_inputs = []
    for inp in inputs:
        try:
            if "%geom" in inp.read_text(encoding="utf-8", errors="ignore"):
                scan_inputs.append(inp)
        except OSError:
            continue
    scan_ran = any(
        (inp.with_suffix(".out")).is_file() for inp in scan_inputs)
    return C.finish([
        (bool(scan_inputs), "an input with a %geom scan/constraint block exists"),
        (scan_ran, "an output sits beside the scan input"),
    ])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
