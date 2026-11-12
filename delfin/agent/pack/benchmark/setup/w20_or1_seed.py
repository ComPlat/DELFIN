#!/usr/bin/env python3
"""Seed start geometries for the OR1 TS tasks (wave 20, nacht-s14).

Each task gets its own folder under <workspace>/chem with a start
geometry (hand-built, tiny systems only) and a README stating charge
and multiplicity.  Refuses to overwrite an existing chem/ tree, like
chem_start_geometries.py does.

Hand-built linear geometries (reasonable bond lengths in Angstrom):
  HCN: H-C 1.064, C-N 1.160 (near-equilibrium linear)
  HNC: N-H 0.996, C-N 1.171 (near-equilibrium linear)
  CH3F: C-F 1.383, C-H 1.095 (tetrahedral)
  F-: single atom
  ethane: staggered, C-C 1.535, C-H 1.095
  ts_candidate: a linear HCN-H-N arrangement near the HCN/HNC TS
  trap_ts: a deliberately BENT "TS" (two soft modes -> not a TS)
"""
from __future__ import annotations

import sys
from pathlib import Path

HCN = """3
linear HCN start
H 0.0 0.0 -1.064
C 0.0 0.0 0.0
N 0.0 0.0 1.160
"""
HNC = """3
linear HNC guess
N 0.0 0.0 -1.171
C 0.0 0.0 0.0
H 0.0 0.0 0.996
"""
TS_CANDIDATE = """3
HCN/HNC transfer TS candidate (linear)
H 0.0 0.0 0.05
C 0.0 0.0 1.02
N 0.0 0.0 2.19
"""
TRAP_TS = """4
bent four-atom candidate (NOT a TS: two soft modes)
H 0.0 0.0 0.60
C 0.0 0.0 -0.50
N 0.0 1.05 -1.05
O 0.0 -1.05 -1.05
"""
CH3F = """5
fluoromethane start
C 0.0 0.0 0.0
F 0.0 0.0 1.383
H 0.941 0.0 -0.365
H -0.470 0.815 -0.365
H -0.470 -0.815 -0.365
"""
F_ANION = """1
fluoride anion
F 0.0 0.0 0.0
"""
ETHANE = """8
ethane staggered start
C 0.0 0.0 0.0
C 0.0 0.0 1.535
H 0.0 1.024 -0.362
H 0.887 -0.512 -0.362
H -0.887 -0.512 -0.362
H 0.0 -1.024 1.897
H -0.887 0.512 1.897
H 0.887 0.512 1.897
"""

FILES = {
    "hcn_hnc": [("start.xyz", HCN), ("hnc_guess.xyz", HNC)],
    "ethane": [("start.xyz", ETHANE)],
    "sn2": [("start.xyz", CH3F), ("f_anion.xyz", F_ANION)],
    "ts_candidate": [("ts_candidate.xyz", TS_CANDIDATE)],
}


def build(ws: Path) -> None:
    root = ws / "chem"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        raise SystemExit(1)
    for folder, files in FILES.items():
        d = root / folder
        d.mkdir(parents=True)
        for name, text in files:
            (d / name).write_text(text, encoding="utf-8")
        (d / "README.txt").write_text(
            f"{folder}: hand-built start geometry, charge 0 unless the "
            "task says otherwise; run everything through DELFIN's steps\n",
            encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    try:
        build(ws)
    except SystemExit as exc:
        return int(exc.code or 1)
    print(f"seeded OR1 TS start geometries under {ws / 'chem'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
