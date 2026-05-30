"""delfin.fffree.linkage_isomers — Universal ambidentate / linkage-isomer enumeration.

For ambidentate ligands (multiple potential donor atoms), enumerate all possible
M-binding configurations. Universal across all known ambidentate ligands.

Common ambidentate ligands:
  SCN⁻  : κ-S (thiocyanato-S)  vs κ-N (isothiocyanato-N)
  SeCN⁻ : κ-Se vs κ-N
  CN⁻   : κ-C  vs κ-N (cyanido)
  NO₂⁻  : κ-N (nitro) vs κ-O (nitrito)
  ONO   : same as NO₂ (linkage isomer)
  N₃⁻   : κ-Nα (terminal) vs κ-Nγ
  NCO⁻  : κ-N vs κ-O (isocyanate)
  DMSO  : κ-O vs κ-S (O- or S-bound)
  SO₃²⁻ : κ-S vs κ-O (sulfito vs sulfonato)

For each ambidentate ligand, generate 2 (or more) M-binding configurations as
distinct isomers in the Pólya enumeration.

Universal: any SMILES containing an ambidentate functional group is detected
via SMARTS patterns + treated as multi-donor in the isomer count.

Env-gate: DELFIN_FFFREE_LINKAGE=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for complete linkage-isomer coverage.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_LINKAGE = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_LINKAGE", "0") == "1"


# Ambidentate ligand fragments + their possible binding atoms.
# Format: { SMARTS_pattern: [(donor_idx_a, label_a), (donor_idx_b, label_b)] }
AMBIDENTATE_LIGANDS = {
    # SCN⁻ / NCS⁻
    "[S]=C=[N-]":  [(0, "κ-S thiocyanato"), (2, "κ-N isothiocyanato")],
    "[S]C#N":      [(0, "κ-S thiocyanato"), (2, "κ-N isothiocyanato")],
    "[N]=C=[S-]":  [(0, "κ-N isothiocyanato"), (2, "κ-S thiocyanato")],
    # SeCN⁻
    "[Se]C#N":     [(0, "κ-Se selenocyanato"), (2, "κ-N isoselenocyanato")],
    # CN⁻
    "[C-]#N":      [(0, "κ-C cyanido"), (1, "κ-N isocyanido")],
    "[N-]=[C]":    [(1, "κ-C cyanido"), (0, "κ-N isocyanido")],
    # NO₂⁻
    "[N+](=O)[O-]": [(0, "κ-N nitro"), (2, "κ-O nitrito")],
    # NCO⁻
    "[N]=C=[O-]":   [(0, "κ-N isocyanato"), (2, "κ-O cyanato")],
    "[O]C#N":       [(0, "κ-O cyanato"),    (2, "κ-N isocyanato")],
    # DMSO (dimethyl sulfoxide)
    "[S](=O)([CH3])[CH3]": [(0, "κ-S dmso"), (1, "κ-O dmso")],
    # Sulfite SO₃²⁻
    "[S](=O)(=O)[O-]": [(0, "κ-S sulfito"), (3, "κ-O sulfonato")],
    # Pyridine-N-oxide (typically O-bound, but possible N)
    "c1ccncc1=O":  [(4, "κ-O"), (3, "κ-N")],
    # Azide N₃⁻
    "[N-]=[N+]=[N-]": [(0, "κ-Nα"), (2, "κ-Nγ")],
}


def detect_ambidentate_groups(mol) -> List[Dict]:
    """Detect ambidentate functional groups in an RDKit Mol.

    Returns list of dicts:
      {pattern, atom_indices, possible_donors, current_donor}

    Universal across all ligand types.
    """
    if not _LINKAGE:
        return []
    try:
        from rdkit import Chem
    except ImportError:
        return []
    out = []
    for smarts, options in AMBIDENTATE_LIGANDS.items():
        try:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                out.append({
                    "pattern": smarts,
                    "atom_indices": list(match),
                    "donor_options": [
                        (int(match[idx]), label) for idx, label in options
                    ],
                    "n_linkage_isomers": len(options),
                })
        except Exception:
            continue
    return out


def enumerate_linkage_isomer_donor_choices(mol) -> List[List[int]]:
    """Enumerate all combinations of κ-binding choices for ambidentate ligands.

    For each ambidentate group, choose one of the possible donors. Cartesian
    product across all ambidentate groups gives all linkage isomers.

    Returns: list of "donor index per ambidentate group" choice tuples.
    """
    if not _LINKAGE:
        return [[]]
    groups = detect_ambidentate_groups(mol)
    if not groups:
        return [[]]
    from itertools import product as _prod
    choice_lists = [list(range(g["n_linkage_isomers"])) for g in groups]
    return list(_prod(*choice_lists))


def get_donor_atom_for_choice(group_info: Dict, choice_idx: int) -> Tuple[int, str]:
    """Return (atom_idx, label) for a chosen ambidentate donor option."""
    return group_info["donor_options"][choice_idx]


if __name__ == "__main__":
    print("=== Ambidentate ligand patterns ===")
    for smarts, opts in AMBIDENTATE_LIGANDS.items():
        labels = " | ".join(lbl for _, lbl in opts)
        print(f"  {smarts:<25} → {labels}")
    print(f"\nTotal ambidentate patterns: {len(AMBIDENTATE_LIGANDS)}")

    # Test detection on a simple SCN-containing complex
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles("[Cl][Fe]([S]C#N)([S]C#N)([S]C#N)")
        if m is not None:
            os.environ["DELFIN_FFFREE_LINKAGE"] = "1"
            import importlib
            import sys
            sys.modules.pop("delfin.fffree.linkage_isomers", None)
            from delfin.fffree.linkage_isomers import detect_ambidentate_groups
            groups = detect_ambidentate_groups(m)
            print(f"\nDetection test on Fe(SCN)₃Cl: {len(groups)} ambidentate groups found")
            for g in groups:
                print(f"  {g['pattern']:<20} atoms={g['atom_indices']}, options={g['donor_options']}")
    except Exception as e:
        print(f"Detection test skipped: {e}")
