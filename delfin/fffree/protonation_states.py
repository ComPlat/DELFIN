"""delfin.fffree.protonation_states — Universal protonation-state enumeration.

For ligands with pH-dependent protonation (carboxylic acids, amines, phenols,
imidazoles, water/hydroxide/oxo), the protonation state significantly
affects coordination chemistry. This module enumerates relevant states.

  Carboxylic acid (COOH) : neutral COOH vs anionic COO⁻
  Amine (NHₓ)            : NHₓ vs deprotonated NHₓ₋₁⁻
  Phenol (Ar-OH)         : Ar-OH vs Ar-O⁻
  Imidazole              : neutral vs anionic (Nδ vs Nε protonation)
  Water/hydroxide/oxo    : H2O ↔ OH⁻ ↔ O²⁻ (pH-dependent)
  Phosphate (HxPO4)      : H3PO4 → H2PO4⁻ → HPO4²⁻ → PO4³⁻

For complete TMC enumeration, multiple protonation states = multiple distinct
isomers contributing to the Pólya orbit count.

Universal SMARTS-pattern detection.
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_PROTON=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_PROTON = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_PROTON", "0") == "1"


# Protonation-state SMARTS patterns + alternatives
# Format: pattern_protonated → list of possible deprotonated states
PROTONATION_STATES = {
    # Carboxylic acid
    "[CX3](=O)[OH]":  [("[CX3](=O)[O-]", "COO-")],
    "[CX3](=O)[OD1]": [("[CX3](=O)[O-]", "COO-")],
    # Primary amine
    "[NH2]":          [("[NH-]", "NH-")],
    "[NH3+]":         [("[NH2]", "NH2 neutral")],
    # Secondary amine
    "[NH1][CX4]":     [("[N-][CX4]", "deprotonated")],
    # Phenol
    "c[OH]":          [("c[O-]", "phenolate")],
    # Imidazole
    "[nH]1ccnc1":     [("[n-]1ccnc1", "imidazolate")],
    # Water/hydroxide/oxo
    "[OH2]":          [("[OH-]", "hydroxide"), ("[O--]", "oxo")],
    "[OH-]":          [("[O--]", "oxo")],
    # Phosphate ladder
    "OP(=O)(O)O":     [("OP(=O)(O)[O-]", "H2PO4-"),
                       ("OP(=O)([O-])[O-]", "HPO4 2-"),
                       ("[O-]P(=O)([O-])[O-]", "PO4 3-")],
    # Sulfonic acid
    "[SX4](=O)(=O)[OH]": [("[SX4](=O)(=O)[O-]", "sulfonate")],
}


def detect_protonation_sites(mol) -> List[Dict]:
    """Detect protonation sites in an RDKit Mol.

    Returns list of {pattern, atom_indices, alternatives}.
    """
    if not _PROTON:
        return []
    try:
        from rdkit import Chem
    except ImportError:
        return []
    out = []
    for patt_smiles, alternatives in PROTONATION_STATES.items():
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
                    "n_states": 1 + len(alternatives),
                })
        except Exception:
            continue
    return out


def enumerate_protonation_states(mol) -> int:
    """Return number of distinct protonation-state combinations for a Mol."""
    if not _PROTON:
        return 1
    sites = detect_protonation_sites(mol)
    if not sites:
        return 1
    total = 1
    for site in sites:
        total *= site["n_states"]
    return total


if __name__ == "__main__":
    print(f"=== Protonation patterns: {len(PROTONATION_STATES)} ===")
    for patt, alts in PROTONATION_STATES.items():
        n_alts = len(alts)
        print(f"  {patt:<30} → {n_alts} alternatives")
    print()
    os.environ["DELFIN_FFFREE_PROTON"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.protonation_states", None)
    from delfin.fffree.protonation_states import (
        detect_protonation_sites, enumerate_protonation_states
    )
    try:
        from rdkit import Chem
        # Acetic acid: 1 site, 2 states
        m = Chem.MolFromSmiles("CC(=O)O")
        sites = detect_protonation_sites(m)
        n = enumerate_protonation_states(m)
        print(f"Acetic acid: {len(sites)} sites detected, {n} total states")
        # Phosphoric acid
        m = Chem.MolFromSmiles("OP(=O)(O)O")
        sites = detect_protonation_sites(m)
        n = enumerate_protonation_states(m)
        print(f"Phosphoric acid: {len(sites)} sites, {n} total states")
        # His (imidazole)
        m = Chem.MolFromSmiles("c1ncnc1")
        sites = detect_protonation_sites(m)
        n = enumerate_protonation_states(m)
        print(f"Imidazole: {len(sites)} sites, {n} total states")
    except Exception as e:
        print(f"RDKit test skipped: {e}")
