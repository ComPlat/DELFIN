#!/usr/bin/env python3
"""Fixture for the w20 classic opt+freq workflow tasks.

Input: the workspace root (argv[1]). Output: w20_<family>_<k>/ with

* start.xyz  -- the practice molecule's start geometry (SMILES-derived;
   the runner's chemistry workspace carries the same molecule so the
   geometry here matches the instance the holdout drew)
* brief.txt  -- the task prompt, verbatim from the holdout sampler.

Exit 1 without writing anything when the folder already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Public-seed draw for w20_classic_opt_freq_1 (see bench_holdout.sample).
_MOLECULE = "nitric_oxide"
_FORMULA = "NO"
_CHARGE = 0
_MULTIPLICITY = 2


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w20_classic_opt_freq.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "w20_classic_opt_freq_1"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "start.xyz").write_text(
        f"2\n{_FORMULA} charge={_CHARGE} mult={_MULTIPLICITY}\n"
        "N 0.0 0.0 0.0\nO 0.0 0.0 1.15\n", encoding="utf-8")
    (d / "brief.txt").write_text(
        "Optimise nitric oxide in the gas phase and confirm the "
        "resulting geometry is a true minimum. Report the lowest-energy "
        "structure's thermochemical correction at 298 K and the count "
        "of imaginary frequencies.\n", encoding="utf-8")
    print(f"prepared {d} -- classic opt+freq fixture")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
