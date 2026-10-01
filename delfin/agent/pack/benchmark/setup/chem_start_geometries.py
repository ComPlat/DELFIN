#!/usr/bin/env python3
"""Seed start geometries for the chemistry tasks (xtb level).

Each molecule the chemistry suite asks about is written as a start.xyz
inside its own subfolder, plus a README stating charge and what the task
wants.  The geometries are RDKit embeddings (seeded, so reproducible)
with two hand-built cases RDKit valence cannot express: Cr(CO)6
(octahedral, literature Cr-C 1.918 A / C-O 1.140 A) and the tert-butyl
cation.  Optimised REFERENCE values live in the tasks file as tolerance
bands; this script only provides the starting point the model sees.
"""

from __future__ import annotations

import sys
from pathlib import Path

MOLS = {
    "acetaminophen": ("CC(=O)Nc1ccc(O)cc1", 0),
    "chloronitrobenzene": ("[Cl]c1ccccc1[N+](=O)[O-]", 0),
    "au_cyanide": ("[Au-](C#N)(C#N)", -1),
    "cat_methyl": ("[CH3+]", 1),
    "cat_ethyl": ("C[CH2+]", 1),
    "cat_tbutyl": ("C[C+](C)(C)", 1),
    "methane": ("C", 0),
    "ethane": ("CC", 0),
    "isobutane": ("CC(C)C", 0),
}

CR_CO6 = """14
Cr(CO)6 hand-built octahedral start geometry
Cr    0.0000   0.0000   0.0000
C     1.9180   0.0000   0.0000
O     3.0580   0.0000   0.0000
C    -1.9180   0.0000   0.0000
O    -3.0580   0.0000   0.0000
C     0.0000   1.9180   0.0000
O     0.0000   3.0580   0.0000
C     0.0000  -1.9180   0.0000
O     0.0000  -3.0580   0.0000
C     0.0000   0.0000   1.9180
O     0.0000   0.0000   3.0580
C     0.0000   0.0000  -1.9180
O     0.0000   0.0000  -3.0580
"""


def build(ws: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for name, (smi, q) in MOLS.items():
        d = ws / "chem" / name
        d.mkdir(parents=True, exist_ok=True)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise RuntimeError(f"SMILES for {name} failed to parse: {smi}")
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            raise RuntimeError(f"embedding failed for {name}")
        (d / "start.xyz").write_text(Chem.MolToXYZBlock(m), encoding="utf-8")
        (d / "README.txt").write_text(
            f"{name}: SMILES {smi}, charge {q}\n", encoding="utf-8")
    d = ws / "chem" / "cr_carbonyl"
    d.mkdir(parents=True, exist_ok=True)
    (d / "start.xyz").write_text(CR_CO6, encoding="utf-8")
    (d / "README.txt").write_text(
        "cr_carbonyl: Cr(CO)6, charge 0 (hand-built octahedral start)\n",
        encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "chem"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    try:
        build(ws)
    except Exception as exc:
        print(f"fixture build failed: {exc}", file=sys.stderr)
        return 1
    print(f"prepared {root} with start geometries")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
