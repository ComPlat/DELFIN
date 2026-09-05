"""_oxidation_state.py — Oxidation state and d-electron count from the LIGAND BALANCE.

WHY (16.08.2026, user).  Three measured open work items hang on the same missing number:

  1. `smiles_converter._PREFERRED_CN4_GEOMETRY` lists `'Cu': 'TH'` -- tetrahedral.  That holds
     for Cu(I) d10.  **Cu(II) is d9 and square-planar up to JT-elongated octahedral**,
     and Cu(II) is by far the more common case.  The same row for both states.  Likewise
     `'Au': 'SQ'` -- correct for Au(III) d8, wrong for Au(I) d10, which is LINEAR.
     The builder already thinks in d-counts (the comments there say "d8 metals",
     "d6 low-spin", "d0-d5, d7, d10") -- but as a hard-wired ELEMENT LIST.
  2. JAHN-TELLER is strong exactly for unequal e_g occupation: **d4 HS, d7 LS, d9**.  The
     e_g orbitals point at the ligands, the distortion lies at 0.2-0.4 A -- measurable.
     (The t2g cases d1/d2/d4 LS/d5 LS/d6 HS/d7 HS lie at ~0.05 A and vanish in the
     noise.)  Measured on 16.08. on `_ml_prior.json`: **Cu-N skew 3.7x at n=12339,
     Cu-O 3.3x at n=7987** -- the long upper tail of the "two long" bonds.  Mn and
     Cr do NOT show it (0.9x / 1.2x), because in the CSD they occur predominantly as Mn(II) d5
     and Cr(III) d3 -- **not JT-active**.  This cannot be decided at the ELEMENT,
     only at the STATE.
  3. The M-D length key: HS and LS differ by 0.1-0.2 A (Co(II) HS ~2.1,
     Co(III) LS ~1.9).  The builder table carries ONE value per (metal, donor).

⚠ WHY NOT FROM THE FORMAL CHARGE.  DELFIN writes formal charges that the SMILES writer
sets for VALENCE SATISFACTION: `[Cu-2]`, `[Re-2]`, `[Ti-2]`, `[O-][Re]`.  These are not
physical oxidation states -- the metal carries a NEGATIVE formal charge there, while
in reality it is positive.  The SUM of all formal charges, by contrast, is correct: it is the
total charge of the complex.

THE CALCULATION is the classical ionic counting:

    OS(metal)  =  total charge(complex)  -  sum(ligand charges)

    Example [Cu(acac)2], neutral:   0 - (-1 -1) = +2  ->  Cu(II), d9   ✓
    Example [Re(=O)Cl3(PPh3)2]:     0 - (-2 -3 +0) = +5 -> Re(V), d2   ✓

    d-count = group number - OS

⚠⚠ THE MOST IMPORTANT RULE: WHEN IN DOUBT, NOTHING.  If even ONE donor cannot be classified,
this module returns `None` -- no estimated state.  A WRONG oxidation state is
worse than none, because it actively spoils the geometry choice: it would treat a Cu(I) d10 as
d9 and force a JT elongation on it that it does not have.  On `None` the caller falls
back to today's behaviour, and that is byte-identical by construction.

Likewise `None` for SEVERAL metals: the total charge can then not be distributed
unambiguously (a mixed-valence Fe(II)/Fe(III) would otherwise be Fe(2.5) twice).

Geometry-free, deterministic, without crystal, without an element table with radii.
"""
from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

_LOG = logging.getLogger(__name__)

# Group number of the transition metals -- d-count = group - OS.
_OX_GROUP: Dict[str, int] = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}

# JT-STRONG: unequal e_g occupation.  (d-count, low_spin) -> True.
# d4 HS (t2g3 eg1) | d7 LS (t2g6 eg1) | d9 (t2g6 eg3, one hole).
# d9 is spin-independent -- there is only one occupation.
_OX_JT_STRONG = {(4, False), (7, True), (9, True), (9, False)}


def _ox_donor_charge(atom, metal_idx: int) -> Optional[int]:
    """Ionic charge of ONE donor atom.  None = not classifiable -> abort.

    Classification is by the element, the heavy non-metal neighbours and the
    H count -- not by the formal charge, which is exactly what is unusable here.
    """
    sym = atom.GetSymbol()
    nh = atom.GetTotalNumHs()
    heavy = [n for n in atom.GetNeighbors()
             if n.GetSymbol() != "H" and n.GetIdx() != metal_idx]
    nheavy = len(heavy)

    if sym in ("F", "Cl", "Br", "I"):
        return -1 if (nheavy == 0 and nh == 0) else None      # halide; otherwise organically bound

    if sym == "O":
        if nheavy == 0 and nh == 0:
            return -2                                          # terminal oxo
        if nheavy == 0 and nh == 1:
            return -1                                          # hydroxide
        if nheavy == 0 and nh == 2:
            return 0                                           # aqua
        if nheavy == 1 and nh == 0:
            return -1                                          # alkoxide / carboxylate O / phenolate
        if nheavy == 2 and nh == 0:
            return 0                                           # ether
        return None

    if sym == "S":
        if nheavy == 0 and nh == 0:
            return -2                                          # sulfido
        if nheavy == 1 and nh == 0:
            return -1                                          # thiolate
        if nheavy == 2 and nh == 0:
            return 0                                           # thioether
        return None

    if sym == "N":
        if nheavy == 0 and nh == 0:
            return -3                                          # nitrido
        if atom.GetIsAromatic():
            return 0                                           # pyridine / imidazole -- L-type
        # ⚠ AMIDE VERSUS AMMINE -- the self-test exposed it (cisplatin read as
        # Pt(IV) d6 instead of Pt(II) d8).  A metal-bound NH3 and an NR2 amide BOTH have
        # "two substituents next to the metal" in the graph; the bond to the metal consumes
        # one valence.  The distinction is made at the HYDROGEN: if the N carries an H, in
        # these SMILES it is practically always a neutral amine/ammine -- a deprotonated amide
        # would be written with an explicit charge.  Only an N with TWO heavy substituents and
        # NO H is an amide.
        if nh > 0:
            return 0                                           # ammine / amine / imine -- L-type
        if nheavy == 1:
            return -2                                          # imido
        if nheavy == 2:
            return -1                                          # amide
        if nheavy == 3:
            return 0                                           # tertiary amine -- L-type
        return None

    if sym == "P":
        if nheavy + nh == 3:
            return 0                                           # phosphane -- L-type
        if nheavy + nh == 2:
            return -1                                          # phosphide
        return None

    if sym == "C":
        # Carbonyl (C with terminal O) and isonitrile are L-type; alkyl/aryl are X-type.
        for n in heavy:
            if n.GetSymbol() == "O" and len([q for q in n.GetNeighbors()
                                             if q.GetSymbol() != "H"]) == 1:
                return 0                                       # CO
        if atom.GetIsAromatic() or nheavy + nh >= 1:
            return -1                                          # alkyl / aryl / cyclopentadienyl C
        return None

    return None                                                # unknown donor -> abort


def oxidation_state(mol, metal_idx: int) -> Optional[Tuple[int, int]]:
    """(oxidation state, d-electron count) or None.

    None means EXPLICITLY "not determinable", not "zero" -- the caller must then fall
    back to its current behaviour.  Reasons for None: several metals, unknown
    donor, unknown metal, unphysical result.
    """
    if mol is None:
        return None
    try:
        m = mol.GetAtomWithIdx(int(metal_idx))
        msym = m.GetSymbol()
        if msym not in _OX_GROUP:
            return None
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP]
        if len(metals) != 1:
            return None                       # distribution over several centres is not unambiguous
        total = int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
        lig = 0
        for n in m.GetNeighbors():
            q = _ox_donor_charge(n, int(metal_idx))
            if q is None:
                return None                   # ONE unknown donor is enough -- guess nothing
            lig += q
        os_ = total - lig
        if os_ < 0 or os_ > 8:
            return None                       # unphysical -> the balance does not add up
        d = _OX_GROUP[msym] - os_
        if d < 0 or d > 10:
            return None
        return os_, d
    except Exception as exc:
        _LOG.debug("oxidation-state: nicht bestimmbar (%s)", exc)
        return None




# ⚠ A JT QUERY IS DELIBERATELY NOT HERE YET.  It would have no caller today -- the
# JT enumeration (three elongation axes per octahedron) is not built.  The reach guard
# reported exactly that, and it is right: unreachable code cannot change anything, and writing
# it anyway is the same mistake this project found three times today.
# `_OX_JT_STRONG` above carries the knowledge as DATA; the query comes with its consumer.


if __name__ == "__main__":  # pragma: no cover -- self-test on known cases
    # Every case is recomputed by hand; the expected state is fixed by the chemistry,
    # not by an earlier run of this module.  A self-test that only confirms its own
    # output from yesterday tests nothing.
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    CASES = [
        ("CC(=O)C=C(C)[O][Cu][O]C(C)=CC(C)=O", 2, 9, "Cu(acac)2 -- d9, JT-aktiv"),
        ("Cl[Pt](Cl)(N)N", 2, 8, "cis-Platin -- d8"),
        ("[O][Re](Cl)(Cl)Cl", 5, 2, "Re(V) oxo-trichlorid"),
        ("[Ni](C#[O])(C#[O])(C#[O])C#[O]", 0, 10, "Ni(CO)4 -- d10, alle L-Typ"),
        ("Cl[Fe](Cl)Cl", 3, 5, "FeCl3 -- d5"),
    ]
    ok = bad = 0
    for smi, exp_os, exp_d, note in CASES:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                except Exception:
                    pass
        if mol is None:
            print(f"  ?? {smi}  -- SMILES nicht lesbar ({note})"); bad += 1; continue
        mi = next((a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP), None)
        r = oxidation_state(mol, mi) if mi is not None else None
        if r is None:
            print(f"  -- {smi}\n     None (nicht bestimmbar) -- erwartet OS {exp_os}, d{exp_d}  [{note}]")
            bad += 1
            continue
        got_os, got_d = r
        mark = "OK" if (got_os, got_d) == (exp_os, exp_d) else "!!"
        ok += (mark == "OK"); bad += (mark != "OK")
        print(f"  {mark} {smi}\n     OS {got_os} d{got_d}   erwartet OS {exp_os} d{exp_d}   [{note}]")
    print(f"\n  {ok} richtig, {bad} nicht.  "
          f"'None' ist KEIN Fehler des Verfahrens -- es ist die Weigerung zu raten.")
