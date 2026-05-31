"""delfin.fffree.electron_counting — Universal 18-electron rule electron counting.

For organometallic complexes, the 18-electron rule (and 16-electron for d8 sqp)
is a fundamental valence-electron bookkeeping that determines geometry stability.

Donor electron counts (covalent / neutral-ligand-method):
  σ-donor (L-type): 2 electrons per donor (e.g. PR3, NH3, CO)
  X-type donor: 1 electron (e.g. Cl, R, OH)
  η²-alkene: 2 electrons
  η³-allyl (X-type): 3 electrons
  η⁵-Cp (X-type): 5 electrons
  η⁶-arene: 6 electrons
  η⁷-cycloheptatrienyl: 7 electrons
  η⁸-COT: 8 electrons
  μ-bridging halide: 3 electrons per metal

For each TMC, compute total valence electrons:
  TVE = d_count(metal) + Σ(ligand_contributions) - charge
  If TVE ≈ 18 (or 16 for d8 sqp): valid 18-electron complex

Universal across all organometallic complexes.
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_ELECTRON_COUNT=1 (default OFF, byte-identical when unset).
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_E_COUNT = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_ELECTRON_COUNT", "0") == "1"


# Donor electron contributions (covalent / neutral-ligand-method)
DONOR_ELECTRONS = {
    "CO": 2, "PR3": 2, "NR3": 2, "py": 2, "NH3": 2, "H2O": 2,
    "bpy": 4, "phen": 4,                # bidentate L-type
    "Cl": 1, "Br": 1, "I": 1, "F": 1,   # X-type halides
    "H": 1,                              # hydride X-type
    "Me": 1, "Et": 1, "R": 1,           # alkyl X-type
    "OH": 1, "OR": 1,                   # alkoxide
    "NCS": 1, "SCN": 1,                 # ambidentate as X
    "Cp": 5,                            # η5-cyclopentadienyl X-type
    "Cp*": 5,                           # pentamethylcyclopentadienyl
    "indenyl": 5,                       # η5-indenyl
    "arene": 6,                         # η6
    "diene": 4,                         # η4
    "alkene": 2,                        # η2
    "allyl": 3,                         # η3 X-type
    "COT": 8,                           # η8
}


def electron_count_for_complex(
    metal: str,
    oxidation_state: int,
    ligands: List[str],
    extra_charge: int = 0,
) -> Dict:
    """Compute total valence electrons for an organometallic complex.

    Parameters
    ----------
    metal : metal symbol (e.g. 'Fe', 'Pt')
    oxidation_state : metal oxidation state (used for d-count)
    ligands : list of ligand labels matching DONOR_ELECTRONS keys
    extra_charge : net charge of complex (e.g. -1 for anionic)

    Returns dict {tve, d_count, ligand_contributions, satisfies_18e}.

    Universal across all organometallic complexes.
    """
    if not _E_COUNT:
        return {"electron_count": "disabled"}
    try:
        from delfin.fffree.spin_states import d_electron_count
        d_n = d_electron_count(metal, oxidation_state)
    except Exception:
        d_n = max(0, 10 - oxidation_state)
    if d_n < 0:
        d_n = 0
    contributions = []
    total = d_n
    for lig in ligands:
        n = DONOR_ELECTRONS.get(lig, 2)  # default L-type
        contributions.append((lig, n))
        total += n
    total -= extra_charge  # extra electrons if anionic, fewer if cationic
    return {
        "metal": metal,
        "oxidation_state": oxidation_state,
        "d_count": d_n,
        "ligand_contributions": contributions,
        "tve": total,
        "satisfies_18e": abs(total - 18) <= 0,
        "satisfies_16e": abs(total - 16) <= 0,
    }


def is_organometallic_stable(tve: int, geometry: str = "oct") -> bool:
    """Check if TVE is consistent with stability rules.

    18e rule: oct, tbp, td, all common organometallic.
    16e rule: d8 sqp.
    14e rule: some 5d Pt(II)/Au(I) low-coord.
    """
    if not _E_COUNT:
        return True
    if geometry == "sqp":
        return tve == 16
    return tve == 18


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_ELECTRON_COUNT"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.electron_counting", None)
    from delfin.fffree.electron_counting import (
        electron_count_for_complex, DONOR_ELECTRONS, is_organometallic_stable
    )

    print(f"=== Donor electron table: {len(DONOR_ELECTRONS)} ===")
    for k, v in DONOR_ELECTRONS.items():
        print(f"  {k}: {v} e-")

    print("\n=== Test complexes ===")
    cases = [
        ("Fe", 0, ["CO"] * 5, 0, "Fe(CO)5"),
        ("Cr", 0, ["arene", "arene"], 0, "Cr(η6-C6H6)2"),
        ("Ni", 0, ["CO"] * 4, 0, "Ni(CO)4"),
        ("Mn", 1, ["CO"] * 5, 0, "Mn(CO)5+"),
        ("Co", 1, ["CO"] * 4, 0, "[Co(CO)4]+"),
        ("Fe", 2, ["Cp", "Cp"], 0, "ferrocene"),
        ("Mo", 0, ["CO"] * 6, 0, "Mo(CO)6"),
        ("Pt", 2, ["Cl", "Cl", "NH3", "NH3"], 0, "cisplatin"),
    ]
    for m, ox, lig, charge, label in cases:
        r = electron_count_for_complex(m, ox, lig, charge)
        print(f"  {label:<18}: d{r['d_count']}, TVE={r['tve']}, "
              f"18e={r['satisfies_18e']}, 16e={r['satisfies_16e']}")
