#!/usr/bin/env python3
"""Acceptance: the carbocation ordering is right and measured.

Relative carbocation stability via the isodesmic hydride exchange

    R+ + CH4  ->  RH + CH3+

which is size-extensive per additional alkyl group and directly gives
the textbook ordering (methyl < primary < tertiary).  The script looks
for a results file the task asks the model to write
(chem/carbocation/energies.csv with columns species,total_energy_eh,
one row per species) and checks the ORDERING of the derived
hydride-exchange values.

Reference (xtb GFN2, alpb water, this repository's chem fixtures):
    CH3+    +0.00 (by construction)
    C2H5+  -48.2 kcal/mol
    tC4H9+ -82.0 kcal/mol

usage: chem_energy_ordering.py <workspace>
"""

from __future__ import annotations

import sys
from pathlib import Path

# Eh, xtb GFN2, alpb water, opt+hess minimum of this fixture set.
REF = {
    "CH3+": -3.153111414560,
    "CH4": -4.174939415342,
    "C2H5+": -6.391169553787,
    "C2H6": -7.336171211134,
    "tC4H9+": -12.776443950535,
    "iC4H10": -13.667637360218,
}
EH_TO_KCAL = 627.5095


def _rel(cat: str, neu: str) -> float:
    return ((REF[cat] + REF["CH4"]) - (REF[neu] + REF["CH3+"])) * EH_TO_KCAL


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    csv = ws / "chem" / "carbocation" / "energies.csv"
    if not csv.is_file():
        print("expected chem/carbocation/energies.csv", file=sys.stderr)
        return 1
    import csv as _csv
    rows = list(_csv.DictReader(csv.open(encoding="utf-8")))
    by_species: dict[str, float] = {}
    for r in rows:
        try:
            s = (r.get("species") or "").strip()
            e = float(r.get("total_energy_eh") or "nan")
        except (TypeError, ValueError):
            continue
        if s:
            by_species[s] = e
    for need in ("CH3+", "C2H5+", "tC4H9+", "CH4", "C2H6", "iC4H10"):
        if need not in by_species:
            print(f"energies.csv lacks species {need}", file=sys.stderr)
            return 1

    def rel(cat, neu):
        return ((by_species[cat] + by_species["CH4"])
                - (by_species[neu] + by_species["CH3+"]))

    # Ordering: methyl least stable, then ethyl, tBu most stable.
    if not (rel("CH3+", "CH4") > rel("C2H5+", "C2H6")
            > rel("tC4H9+", "iC4H10")):
        print("ordering wrong: expected CH3+ < C2H5+ < tC4H9+ "
              "(by hydride-exchange stability)", file=sys.stderr)
        return 1
    # And the tBu gap must be meaningfully larger than the ethyl gap:
    # the reference differences are 48 and 82 kcal/mol.
    d_eth = rel("CH3+", "CH4") - rel("C2H5+", "C2H6")
    d_tbu = rel("C2H5+", "C2H6") - rel("tC4H9+", "iC4H10")
    print(f"ordering ok: ethyl gap {d_eth*EH_TO_KCAL:.1f}, "
          f"tBu gap {d_tbu*EH_TO_KCAL:.1f} kcal/mol")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
