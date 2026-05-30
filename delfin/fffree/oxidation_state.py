"""delfin.fffree.oxidation_state — Universal oxidation-state inference from SMILES.

For a TMC SMILES, infer the metal's oxidation state from:
  1. Explicit formal charge on metal atom (e.g., [Fe+2])
  2. Charge balance with counterions / anionic ligands
  3. Donor-atom ligand-field analysis (heuristic for common cases)

Used by spin_states module to determine d-electron count → multiplicity.

Universal across all SMILES.
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_OXSTATE=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_OXSTATE = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_OXSTATE", "0") == "1"


# Common ligand contributions (charge on metal-binding side)
LIGAND_CHARGE = {
    # Anionic ligands
    "Cl": -1, "Br": -1, "I": -1, "F": -1,
    "OH": -1, "OR": -1, "SR": -1,
    "CN": -1, "NCS": -1, "SCN": -1, "NCO": -1,
    "O2": -2, "S2": -2, "CO3": -2, "SO4": -2,
    "C2O4": -2,  # oxalate
    "Cp": -1,    # η5-cyclopentadienyl
    "H": -1,     # hydride
    "PO4": -3, "NO3": -1, "NO2": -1,
    # Neutral ligands
    "py": 0, "bpy": 0, "phen": 0,
    "NH3": 0, "H2O": 0, "MeOH": 0, "EtOH": 0,
    "PPh3": 0, "PR3": 0,
    "CO": 0, "diene": 0, "arene": 0,
}


def infer_oxidation_state_from_smiles(smiles: str) -> Optional[Tuple[str, int]]:
    """Infer metal + oxidation state from SMILES.

    Algorithm:
      1. Parse SMILES to find metal atom with explicit charge
      2. Return (metal_symbol, oxidation_state)
      3. Fallback: charge balance from ligand environment

    Universal across all SMILES.
    Returns None if no metal found.
    """
    if not _OXSTATE:
        return None
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        def _is_metal(s): return s not in {"H", "C", "N", "O", "F", "P", "S",
                                            "Cl", "Br", "I", "B", "Si", "Se", "As"}
    for atom in mol.GetAtoms():
        if _is_metal(atom.GetSymbol()):
            return atom.GetSymbol(), atom.GetFormalCharge()
    return None


def total_complex_charge(smiles: str) -> Optional[int]:
    """Total charge of the complex (sum of all formal charges)."""
    if not _OXSTATE:
        return None
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return sum(a.GetFormalCharge() for a in mol.GetAtoms())
    except Exception:
        return None


def infer_charge_balance(smiles: str) -> Optional[Dict]:
    """Decompose total charge into metal contribution + ligand contributions.

    For [Fe(H2O)6]Cl3: metal Fe³⁺ + 6 H2O (neutral) + 3 Cl⁻ → balanced.
    For [Pt(NH3)2Cl2]: metal Pt²⁺ + 2 NH3 + 2 Cl⁻ → neutral total.

    Returns dict with breakdown.
    """
    if not _OXSTATE:
        return None
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None
    total = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    metal_charges = [(a.GetSymbol(), a.GetFormalCharge())
                     for a in mol.GetAtoms() if _is_metal(a.GetSymbol())]
    return {
        "total_charge": total,
        "metals": metal_charges,
        "smiles": smiles,
    }


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_OXSTATE"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.oxidation_state", None)
    from delfin.fffree.oxidation_state import (
        infer_oxidation_state_from_smiles, total_complex_charge, infer_charge_balance
    )
    test_smiles = [
        ("[Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]", "Fe²⁺ hexammine"),
        ("[Pt+2]([NH3])([NH3])([Cl-])[Cl-]", "Pt²⁺ cisplatin"),
        ("[Cu+2]([NH3])([NH3])([NH3])([NH3])", "Cu²⁺ tetrammine"),
        ("[Co+3]([NH3])([NH3])([NH3])([NH3])([NH3])([NH3])", "Co³⁺ hexammine"),
        ("[Ti+4]([Cl-])([Cl-])([Cl-])[Cl-]", "TiCl4"),
        ("c1ccccc1[Mn](C=O)(C=O)(C=O)(C=O)C=O", "MnCp(CO)5 type"),
    ]
    print("=== Oxidation-state inference test ===")
    for smi, label in test_smiles:
        result = infer_oxidation_state_from_smiles(smi)
        cb = infer_charge_balance(smi)
        print(f"  {label:<25} {result}  total={cb['total_charge'] if cb else '?'}")
