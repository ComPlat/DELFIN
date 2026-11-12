#!/usr/bin/env python3
"""Seed start geometries for the XT1 xtb-method practice tasks.

Six systems, each in its own subfolder under <workspace>/xtb/ with a
start.xyz and a README stating SMILES, charge and multiplicity:

  water_cation      H2O+   (3 atoms,  doublet)          -- cation, open shell
  methylamine_radical CH3NH. (6 atoms, doublet, hand-built) -- neutral radical
  oh_anion          OH-    (2 atoms,  anion, singlet)   -- anion
  butanone          C4H8O  (13 atoms, neutral)          -- in-solvent opt target
  ammonium          NH4+   (5 atoms,  closed shell cation)
  methoxy_radical   CH3O.  (5 atoms,  doublet, hand-built) -- the wrong-uhf
                    trap task grades against a triplet run of this

Geometries are RDKit embeddings (seeded, reproducible); the two species
RDKit valence cannot express are written as hand-built xyz. Reference
values live as tolerance bands in tasks_w20_xt1.yaml; this script only
provides the start the model sees.
"""

from __future__ import annotations

import sys
from pathlib import Path

# key -> (SMILES, charge, multiplicity)
MOLS = {
    "water_cation": ("[OH2+]", 1, 2),
    "oh_anion": ("[OH-]", -1, 1),
    "butanone": ("CCCC=O", 0, 1),
    "ammonium": ("[NH4+]", 1, 1),
}

# RDKit cannot express the monovalent radical N center of CH3NH. as a
# normal molecule; hand-built start geometries instead.
METHYLAMINE_RADICAL = """6
methylamine radical, doublet start (hand-built)
N    0.0000   0.0000   0.0000
C    1.4700   0.0000   0.0000
H    1.0000   0.9400   0.8100
H    1.8000  -0.9400   0.0000
H    2.1000   0.9200   0.0000
H    1.8000   0.5000  -0.8600
"""

METHOXY_RADICAL = """5
methoxy radical, doublet start (hand-built)
O    0.0000   0.0000   0.0000
C    1.4300   0.0000   0.0000
H    1.8500   1.0200   0.0000
H    1.8500  -0.5100   0.8900
H    1.8500  -0.5100  -0.8900
"""


def build(ws: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    for name, (smi, q, mult) in MOLS.items():
        d = ws / "xtb" / name
        d.mkdir(parents=True, exist_ok=True)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise RuntimeError(f"SMILES for {name} failed to parse: {smi}")
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            raise RuntimeError(f"embedding failed for {name}")
        (d / "start.xyz").write_text(Chem.MolToXYZBlock(m), encoding="utf-8")
        (d / "README.txt").write_text(
            f"{name}: SMILES {smi}, charge {q}, multiplicity {mult}\n",
            encoding="utf-8")
    for name, block, q, mult in (
        ("methylamine_radical", METHYLAMINE_RADICAL, 0, 2),
        ("methoxy_radical", METHOXY_RADICAL, 0, 2),
    ):
        d = ws / "xtb" / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "start.xyz").write_text(block, encoding="utf-8")
        (d / "README.txt").write_text(
            f"{name}: hand-built start, charge {q}, multiplicity {mult}\n",
            encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "xtb"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    try:
        build(ws)
    except Exception as exc:
        print(f"fixture build failed: {exc}", file=sys.stderr)
        return 1
    print(f"prepared {root} with XT1 start geometries")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
