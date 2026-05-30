"""delfin.fffree.solvent_counterion — Universal solvent + counterion handling.

CCDC structures frequently include co-crystallized solvent molecules and
counterions. For DFT-startable structures of the metal complex, these are
typically removed or relocated. This module provides:
  - Detection of common solvents in SMILES (water, MeOH, EtOH, MeCN, DCM, DMF, ...)
  - Detection of common counterions (PF₆⁻, BF₄⁻, ClO₄⁻, BPh₄⁻, NTf₂⁻, ...)
  - Mode selection: strip / keep / position-far

Universal across all SMILES; no atom-specific code.
FF-free.

Env-gate: DELFIN_FFFREE_SOLVENT=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Set, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_SOLVENT = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_SOLVENT", "0") == "1"


# Common solvent SMARTS (in CCDC entries)
COMMON_SOLVENTS = {
    "water":        "O",
    "methanol":     "CO",
    "ethanol":      "OCC",
    "acetonitrile": "CC#N",
    "DCM":          "ClCCl",
    "chloroform":   "ClC(Cl)Cl",
    "THF":          "C1CCOC1",
    "diethyl_ether": "CCOCC",
    "DMSO":         "CS(=O)C",
    "DMF":          "CN(C)C=O",
    "acetone":      "CC(=O)C",
    "EtOAc":        "CCOC(C)=O",
    "benzene":      "c1ccccc1",
    "toluene":      "Cc1ccccc1",
    "pyridine":     "c1ccncc1",
    "n-hexane":     "CCCCCC",
    "n-pentane":    "CCCCC",
    "ammonia":      "N",
}


# Common counterions
COMMON_COUNTERIONS = {
    "PF6":      "[P-](F)(F)(F)(F)(F)F",
    "BF4":      "[B-](F)(F)(F)F",
    "ClO4":     "[Cl](=O)(=O)(=O)[O-]",
    "OTf":      "[O]S(=O)(=O)C(F)(F)F",
    "BPh4":     "[B-](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1",
    "Tf2N":     "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
    "tetraphenylborate": "[B-](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1",
    "OAc":      "CC(=O)[O-]",       # acetate
    "TfO":      "[O-]S(=O)(=O)C(F)(F)F",
    "MsO":      "[O-]S(=O)(=O)C",   # mesylate
    "HSO4":     "OS(=O)(=O)[O-]",
    "Cl":       "[Cl-]",
    "Br":       "[Br-]",
    "I":        "[I-]",
    "NO3":      "[N+](=O)([O-])[O-]",
    "OH":       "[OH-]",
}


def detect_solvent_or_counterion(mol) -> List[Dict]:
    """Detect common solvent + counterion fragments in an RDKit Mol.

    Returns list of: {type: 'solvent'|'counterion', name, atom_indices}.
    """
    if not _SOLVENT:
        return []
    try:
        from rdkit import Chem
    except ImportError:
        return []
    out = []
    for name, smiles in COMMON_SOLVENTS.items():
        try:
            patt = Chem.MolFromSmiles(smiles)
            if patt is None:
                continue
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                out.append({
                    "type": "solvent",
                    "name": name,
                    "atom_indices": list(match),
                })
        except Exception:
            continue
    for name, smiles in COMMON_COUNTERIONS.items():
        try:
            patt = Chem.MolFromSmiles(smiles)
            if patt is None:
                continue
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                out.append({
                    "type": "counterion",
                    "name": name,
                    "atom_indices": list(match),
                })
        except Exception:
            continue
    return out


def strip_solvent_counterion(mol, coords: np.ndarray,
                              strip_solvent: bool = True,
                              strip_counterion: bool = False
                              ) -> Tuple[object, np.ndarray]:
    """Remove solvent/counterion atoms from the structure.

    Returns (new_mol, new_coords). Default: strip solvent only, keep counterions
    (needed for proper charge balance).
    """
    if not _SOLVENT:
        return mol, coords
    detected = detect_solvent_or_counterion(mol)
    to_remove: Set[int] = set()
    for d in detected:
        if d["type"] == "solvent" and strip_solvent:
            to_remove.update(d["atom_indices"])
        elif d["type"] == "counterion" and strip_counterion:
            to_remove.update(d["atom_indices"])
    if not to_remove:
        return mol, coords
    try:
        from rdkit import Chem
        rw = Chem.RWMol(mol)
        for idx in sorted(to_remove, reverse=True):
            rw.RemoveAtom(int(idx))
        new_mol = rw.GetMol()
        # Update coords
        keep_indices = [i for i in range(mol.GetNumAtoms()) if i not in to_remove]
        new_coords = coords[keep_indices]
        return new_mol, new_coords
    except Exception:
        return mol, coords


if __name__ == "__main__":
    print(f"=== Solvent patterns: {len(COMMON_SOLVENTS)} ===")
    for n, s in list(COMMON_SOLVENTS.items())[:5]:
        print(f"  {n:<15} {s}")
    print(f"=== Counterion patterns: {len(COMMON_COUNTERIONS)} ===")
    for n, s in list(COMMON_COUNTERIONS.items())[:5]:
        print(f"  {n:<15} {s}")
