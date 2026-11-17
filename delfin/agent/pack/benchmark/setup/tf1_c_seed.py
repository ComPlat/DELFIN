#!/usr/bin/env python3
"""Seed task C for the TF1 benchmark: reject a barrier built from mixed levels
of theory.

Writes <workspace>/tf1/c/levels.txt -- a per-stage level record for a reaction
barrier whose reactant, TS and product are NOT all at one level of theory
(presented in the same key convention delfin.ts_finder.check_single_level_of_
theory reads). Here TS uses def2-TZVP while reactant/product/ZPE use def2-SVP,
so the honest verdict is "reject": a barrier assembled from single points at
different levels is not a barrier (the ts-finder preset requires ONE level of
theory for reactant, TS, product and ZPE).

The ts-finder expert must reject this and write <workspace>/tf1/c/verdict.json
with the honest label {"verdict": "mixed levels of theory"} (the library check
check_single_level_of_theory returns ok=False for this input). The level names
are DELFIN CONTROL vocabulary (functional/basis), not measurements.

usage: tf1_c_seed.py <workspace>
"""
from __future__ import annotations

import sys
from pathlib import Path


def _levels() -> str:
    return "\n".join([
        "reactant PBE0/def2-SVP",
        "ts PBE0/def2-TZVP",
        "product PBE0/def2-SVP",
        "zpe PBE0/def2-SVP",
    ]) + "\n"


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv
    if len(argv) != 2:
        print("usage: tf1_c_seed.py <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1])
    folder = root / "tf1" / "c"
    if folder.exists():
        print(f"refusing to overwrite existing {folder}", file=sys.stderr)
        return 1
    folder.mkdir(parents=True)
    (folder / "levels.txt").write_text(_levels(), encoding="utf-8")
    print(f"seeded {folder}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
