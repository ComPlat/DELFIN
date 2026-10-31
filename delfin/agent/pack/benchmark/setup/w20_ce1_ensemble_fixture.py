#!/usr/bin/env python3
"""Fixture for the four w20_ce1_* tasks.

Input: the workspace root (argv[1]). Output: w20_ce1_ensemble/ with

* ensemble_default.xyz     -- 4 frames, energies gfn2-tagged in the
                              comment lines; the lowest two are 0.005760
                              Eh apart (= 3.61455... kcal/mol).
* ensemble_degenerate.xyz  -- 3 frames, frames 2 and 3 degenerate within
                              0.001 kcal/mol (identical energy).
* ensemble_mixed.xyz       -- 2 frames, one gfn2- and one r2scan-3c-
                              tagged: the two-levels trap.
* ensemble_untagged.xyz    -- 4 frames, plain CREST comment numbers.
* README.txt               -- the assignment text (threshold 4.0 kcal/mol).

Exit 1 without writing anything when w20_ce1_ensemble/ already exists or
a precondition fails. Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_FRAME = """\
3
{comment}
C 0.0 0.0 0.0
H 1.09 0.0 0.0
H -0.55 {y} 0.0
"""

_ENERGY_DEFAULT = [-12.345670, -12.339910, -12.331000, -12.325000]
_ENERGY_DEGENERATE = [-12.345670, -12.340000, -12.340000]
_ENERGY_UNTAGGED = [-12.345670, -12.339910, -12.331000, -12.320000]


def _ensemble(energies, tags) -> str:
    ys = [0.95, -0.95, 0.30, -0.30]
    blocks = []
    for i, e in enumerate(energies):
        comment = f"energy: {e:.6f} Eh"
        if tags and tags[i]:
            comment += f" {tags[i]}"
        blocks.append(_FRAME.format(comment=comment, y=ys[i % len(ys)]))
    return "".join(blocks)


_README = (
    "w20_ce1_ensemble/ holds CREST-style conformer ensembles (multi-frame\n"
    "XYZ, energy in each comment line, energies in Eh).\n\n"
    "For the threshold task (ensemble_untagged.xyz): apply the CENSO-style\n"
    "prescreening threshold of 4.0 kcal/mol on the relative energies at\n"
    "the ensemble's single level of theory; every frame whose relative\n"
    "energy is <= 4.0 kcal/mol survives into the DFT workflow.\n"
)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w20_ce1_ensemble_fixture.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "w20_ce1_ensemble"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    # Precondition: the two lowest gfn2 frames are 0.005760 Eh apart, which
    # is 3.6145... kcal/mol -- the gap task 1 expects (tolerance 0.05).
    gap_eh = _ENERGY_DEFAULT[1] - _ENERGY_DEFAULT[0]
    if abs(gap_eh * 627.5094740631 - 3.6145) > 0.05:
        print("fixture gap drifted from the expected 3.6145 kcal/mol",
              file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "ensemble_default.xyz").write_text(
        _ensemble(_ENERGY_DEFAULT, ["gfn2", "gfn2", "gfn2", "gfn2"]),
        encoding="utf-8")
    (d / "ensemble_degenerate.xyz").write_text(
        _ensemble(_ENERGY_DEGENERATE, ["gfn2", "gfn2", "gfn2"]),
        encoding="utf-8")
    (d / "ensemble_mixed.xyz").write_text(
        _ensemble([-12.345670, -12.300000], ["gfn2", "r2scan-3c"]),
        encoding="utf-8")
    (d / "ensemble_untagged.xyz").write_text(
        _ensemble(_ENERGY_UNTAGGED, None), encoding="utf-8")
    (d / "README.txt").write_text(_README, encoding="utf-8")
    print(f"prepared {d} -- 4 ensembles + README; no binary called, the "
          "energies are fixture numbers")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
