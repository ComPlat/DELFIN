"""delfin.fffree.tautomers — Universal tautomer enumeration for TMC ligands.

Many ligand functional groups exist as tautomers (keto/enol, amide/imidic acid,
imine/enamine, lactam/lactim). For complete enumeration, each tautomer can
present different donor atoms to the metal:

  keto/enol: C=O ↔ C-OH (carbonyl vs hydroxyl as donor)
  amide/imidic: C(=O)NHR ↔ C(OH)=NR (O vs N donor change)
  imine/enamine: C=N-H ↔ C-NH-H (N vs C donor)
  imidazole tautomers: Nδ-H ↔ Nε-H (different N donor)
  lactam/lactim: cyclic amide tautomers
  azole tautomers: triazole/tetrazole N-positions
  thione/thiol: C=S ↔ C-SH

For complete TMC isomer enumeration: each tautomer × Pólya orbit = distinct
DFT-startable isomer with potentially different bonding.

Universal SMARTS pattern detection.
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_TAUTOMERS=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_TAUTOMERS = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_TAUTOMERS", "0") == "1"


# Tautomer patterns: {SMARTS_form_A: [(SMARTS_form_B, description)]}
TAUTOMER_PAIRS = {
    # Keto / enol
    "C(=O)[CX4H]":     [("C(O)=[CX3H]", "enol")],
    # Amide / imidic acid
    "[CX3](=O)[NH1]":  [("[CX3]([OH1])=[NX2]", "imidic acid")],
    # Imine / enamine
    "[CX3]=[NH1]":     [("[CX4]-[NH2]", "enamine")],
    # Imidazole tautomers (1H vs 3H)
    "[nH]1cncc1":      [("n1c[nH]cc1", "3H tautomer")],
    "[nH]1ccnc1":      [("n1cc[nH]c1", "3H tautomer")],
    # Lactam / lactim
    "O=C1[NH1][CX4]CC1":  [("OC1=[NX2][CX4]CC1", "lactim")],
    # Triazole tautomers
    "[nH]1ncnc1":      [("n1[nH]cnc1", "2H triazole"),
                        ("n1nc[nH]c1", "4H triazole")],
    "[nH]1nncc1":      [("n1[nH]ncc1", "2H triazole")],
    # Tetrazole tautomers
    "[nH]1nnnc1":      [("n1[nH]nnc1", "2H tetrazole")],
    # Pyrazole / 4H-pyrazole
    "[nH]1nccc1":      [("n1[nH]ccc1", "2H pyrazole")],
    # Thione / thiol
    "[CX3](=S)[NH1]":  [("[CX3]([SH1])=[NX2]", "thiol-imine")],
}


def detect_tautomer_sites(mol) -> List[Dict]:
    """Detect tautomerizable functional groups in an RDKit Mol.

    Returns list of {pattern, alternatives, atom_indices, n_tautomers}.
    """
    if not _TAUTOMERS:
        return []
    try:
        from rdkit import Chem
    except ImportError:
        return []
    out = []
    for patt_smiles, alternatives in TAUTOMER_PAIRS.items():
        try:
            patt = Chem.MolFromSmarts(patt_smiles)
            if patt is None:
                continue
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                out.append({
                    "pattern": patt_smiles,
                    "atom_indices": list(match),
                    "alternatives": alternatives,
                    "n_tautomers": 1 + len(alternatives),
                })
        except Exception:
            continue
    return out


def enumerate_tautomer_count(mol) -> int:
    """Return number of distinct tautomer combinations for a molecule."""
    if not _TAUTOMERS:
        return 1
    sites = detect_tautomer_sites(mol)
    if not sites:
        return 1
    total = 1
    for s in sites:
        total *= s["n_tautomers"]
    return total


if __name__ == "__main__":
    print(f"=== Tautomer patterns: {len(TAUTOMER_PAIRS)} ===")
    for patt, alts in TAUTOMER_PAIRS.items():
        n = len(alts)
        labels = ", ".join(lbl for _, lbl in alts)
        print(f"  {patt:<30} → {n} ({labels})")

    os.environ["DELFIN_FFFREE_TAUTOMERS"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.tautomers", None)
    from delfin.fffree.tautomers import detect_tautomer_sites, enumerate_tautomer_count
    try:
        from rdkit import Chem
        cases = [
            ("CC(=O)C", "acetone (keto-enol)"),
            ("CC(=O)N", "acetamide (amide-imidic)"),
            ("c1[nH]ncn1", "1H-imidazole tautomers"),
        ]
        for smi, label in cases:
            m = Chem.MolFromSmiles(smi)
            if m is None: continue
            sites = detect_tautomer_sites(m)
            n = enumerate_tautomer_count(m)
            print(f"  {label}: {len(sites)} sites, {n} tautomer states")
    except Exception as e:
        print(f"RDKit test skipped: {e}")
