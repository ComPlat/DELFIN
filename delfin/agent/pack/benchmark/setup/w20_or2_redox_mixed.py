"""Setup for the OR2 redox trap task (w20-or2-05).

Writes two fake ORCA-style result records of ONE redox couple computed
at DIFFERENT levels of theory. The trap: an agent that subtracts them
reports a meaningless potential. No reference electrode is stated in
the task prompt — the correct answer refuses, names each level, and
demands a common level plus an explicit electrode (E_ref).

MAX_ATOMS documents the small-systems limit this package is graded on.
"""

from __future__ import annotations

import pathlib
import sys

MAX_ATOMS = 16


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(f"usage: {pathlib.Path(argv[0]).name} <workspace>",
              file=sys.stderr)
        return 2
    run(pathlib.Path(argv[1]).resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))


# Deliberately mixed: oxidised form at B3LYP/def2-SVP, reduced form at
# PBE0/def2-TZVP. The energies are plausible, not interchangeable.
OX_RECORD = (
    "oxidised.xyz_level.txt:\n"
    "charge: +1  multiplicity: 2\n"
    "method: B3LYP  basis: def2-SVP  solvent: CPCM(acetonitrile)\n"
    "final Gibbs energy (Eh): -411.982314\n"
)
RED_RECORD = (
    "reduced.xyz_level.txt:\n"
    "charge: 0  multiplicity: 1\n"
    "method: PBE0  basis: def2-TZVP  solvent: CPCM(acetonitrile)\n"
    "final Gibbs energy (Eh): -411.998902\n"
)


def run(workspace: pathlib.Path) -> dict:
    """Write the mixed-level records under <workspace>/or2/redox_mixed/."""
    d = workspace / "or2" / "redox_mixed"
    d.mkdir(parents=True, exist_ok=True)
    (d / "ox_result.txt").write_text(OX_RECORD, encoding="utf-8")
    (d / "red_result.txt").write_text(RED_RECORD, encoding="utf-8")
    return {"folder": str(d), "files": ["ox_result.txt", "red_result.txt"]}
