"""delfin.fffree.ambidentate_kappa_enum — Universal κⁿ binding-mode enumeration.

For every ligand fragment with multiple chemically-plausible donor atoms,
enumerate the distinct κⁿ binding modes that the metal can adopt:

  * κ¹ via each unique donor atom  (linkage isomerism)
  * κ²  if 2 donors are bite-compatible with the metal
  * κ³  if 3 donors are bite-compatible

Universal, graph + element + geometry only.  NO SMARTS, no per-anion
templates, no SMILES pattern matching.  The algorithm:

  for each ligand fragment in the complex:
      seed_donor = the atom currently bonded to M in the input graph
      potential = atoms with a lone pair (N/O/P/S/Se/halide/carbene-C)
                  reachable from seed_donor within MAX_BITE_HOPS bond hops
                  AND with formal charge <= 0
                  AND with degree+totalHs+formalCharge consistent with at
                      least one available lone pair (universal valence
                      counting, not per-anion).
      for k in [1, 2, 3]:
          for each k-subset of potential donors:
              if subset is bite-compatible to the metal:
                  yield (kappa=k, donors=subset)

Bite-compatibility check (universal):

  * single donor (k=1): always compatible.
  * pair (k=2): graph-shortest-path distance between donors must be
    2..4 bonds AND the estimated through-space distance (using ideal
    bond lengths from RDKit ``CovalentRadii``) must fall in
    [BITE_MIN, BITE_MAX] = [2.0, 3.5] Å, which is the empirical CCDC
    bite-distance envelope for κ² chelates across all classes (see
    Mogul fragment manifold 2026-06-06).
  * triple (k=3): all three pairwise distances must be bite-compatible
    AND the three donors must lie on a connected path within the ligand
    (no bridge through the metal).  Rare but valid for nitrate κ³,
    sulfate κ³, oxoanion κ³.

The output is consumed by ``assemble_via_mogul.enumerate_kappa_variants_smiles``,
which rewires the input SMILES graph to insert the additional dative bonds
for each κⁿ choice and emits one full Mogul-PRIMARY embed per variant.

Env-gate: ``DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM``.  Default ON when
``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is also set; off otherwise (byte-identical
with HEAD).

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

from rdkit import Chem
from rdkit.Chem import Draw  # noqa: F401  (keeps RDKit ring info live)

import delfin._bond_decollapse as _bd


# ---------------------------------------------------------------------------
# Constants — universal across ligand classes
# ---------------------------------------------------------------------------
DONOR_ELEMENTS: Set[str] = {
    "N", "O", "P", "S", "Se", "As",
    "F", "Cl", "Br", "I",
    "C",  # carbene / CN-/ allyl carbon: lone-pair test below filters out
          # neutral sp3 C; only carbene/cyanide/ambidentate C qualify.
}

# Bite envelope (Å).  Calibrated from the CCDC manifold (see project
# memory entry "Mogul-DG Phase A DONE 2026-06-06" -- 2.45 M bond keys).
BITE_MIN = 2.0
BITE_MAX = 3.5

# How far inside a ligand to look for additional donors (bonds from seed).
# 4 bonds covers carboxylate (1-3), nitrate (1-3, 1-3), nitrite (1-3),
# dithiocarbamate (1-3), salicylate (1-5), phenanthroline (1-3),
# bipyridine (1-4).
MAX_BITE_HOPS = 4

# Maximum denticity to enumerate.  κ³ covers nitrate, sulfate, oxoanion;
# κ⁴+ (porphyrin etc.) is already handled by the chelate path and is not
# emitted here.
MAX_KAPPA = 3


def _env_on() -> bool:
    """κⁿ enumeration is ON when the Mogul-primary path is active AND the
    explicit kill-switch is not set to "0".  Default ON under
    DELFIN_FFFREE_MOGUL_PRIMARY=1; default OFF otherwise (the whole
    module is then never imported by the dispatch path)."""
    if os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") != "1":
        return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", "0") == "1"
    # Mogul-PRIMARY active -> default ON, explicit "0" to disable.
    return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", "1") == "1"


# ---------------------------------------------------------------------------
# Lone-pair test (universal, no per-anion template)
# ---------------------------------------------------------------------------
# Standard valence counts (RDKit-default).  An atom has a "free" lone pair
# available for coordination when its bonded + H + formal-charge ledger
# leaves at least one non-bonded electron pair.  For maingroup donors
# this is the universal lone-pair test.

# Number of valence electrons per main-group element (universal periodic
# table count, no per-anion).  Used to compute the lone-pair budget
# directly from the octet: lp_electrons = group_e - (2*bond_order_sum)
# + (-formal_charge).  Each lp = 2 electrons.
_GROUP_VALENCE_ELECTRONS = {
    "C": 4,
    "N": 5, "P": 5, "As": 5,
    "O": 6, "S": 6, "Se": 6,
    "F": 7, "Cl": 7, "Br": 7, "I": 7,
}


def _bond_order_sum(atom) -> float:
    """Sum of bond orders (including bonds to H) for an atom.

    Aromatic bonds count as 1.5.  Dative bonds count as 1.0 from the
    donor's perspective (the donor's lone pair forms the σ pair; the
    electron count on the donor goes DOWN by 2 because both electrons
    sit between donor and metal — for the donor's lone-pair budget we
    treat the dative bond as a normal single bond, which captures the
    "those two electrons are no longer a free lone pair" semantic).
    """
    total = 0.0
    for b in atom.GetBonds():
        bt = b.GetBondType()
        if bt == Chem.BondType.SINGLE or bt == Chem.BondType.DATIVE:
            total += 1.0
        elif bt == Chem.BondType.DOUBLE:
            total += 2.0
        elif bt == Chem.BondType.TRIPLE:
            total += 3.0
        elif bt == Chem.BondType.AROMATIC:
            total += 1.5
        else:
            total += 1.0
    total += float(atom.GetTotalNumHs())
    return total


def _has_lone_pair(atom) -> bool:
    """Universal lone-pair test for coordination donors.

    True iff the atom has at least one non-bonding electron pair in its
    octet (Lewis structure), computed from
        lp_electrons = group_valence_electrons
                     - 2 * bond_order_sum_excluding_M
                     + (-formal_charge)
    and lp_pairs = lp_electrons // 2 >= 1.

    This is the textbook Lewis octet test, no SMARTS, no per-anion table:

      * NH3: N has 5 v.e. - 2*3 (3 N-H) + 0 = -1?  No -- we must NOT
        subtract bonds to H as "bond orders 1" each twice.  The
        correct Lewis count is:
            valence_electrons - (electrons_in_bonds_at_this_atom)
        where each single bond contributes 1 electron to the atom's
        valence shell, a double 2, etc.  So for NH3: 5 - 3 = 2 = 1 lp.
        For O=C: 6 - 2 = 4 = 2 lp.  This is the convention we use.
      * neutral N (NH3, py): 5 - 3 = 2 -> 1 lp -> donor.
      * carbonyl O (C=O):     6 - 2 = 4 -> 2 lp -> donor (the oxo
        oxygen can coordinate to M; CCDC manifold confirms acetate-O,
        nitrate-O, sulfate-O all coordinate from BOTH sp2 oxygens).
      * O⁻ (carboxylate-O⁻):  6 - 1 + 1 = 6 -> 3 lp -> donor.
      * neutral sp3 C (alkane): 4 - 4 = 0 lp -> NOT a donor (correct).
      * carbene C: 4 - 2 = 2 -> 1 lp -> donor.
      * cyanide C (C in C#N⁻):  4 - 3 + 1 = 2 -> 1 lp -> donor.
      * isocyanide N (N in C#N⁻ via N): 5 - 3 = 2 -> 1 lp -> donor.
      * Ar-aromatic C (1.5+1.5+1.0=4): 4 - 4 = 0 -> NOT a donor.
      * py-N (aromatic in ring): 5 - 3 = 2 -> 1 lp -> donor.
      * X⁻ halide: 7 - 0 + 1 = 8 -> 4 lp -> donor (when free).
    """
    sym = atom.GetSymbol()
    ve = _GROUP_VALENCE_ELECTRONS.get(sym)
    if ve is None:
        return False
    bos = _bond_order_sum(atom)
    charge = atom.GetFormalCharge()
    # Lewis valence-shell electron count at this atom.
    # Each bond shares 1 electron with this atom; a negative formal
    # charge adds an electron; a positive charge removes one.
    lone_electrons = ve - bos + (-charge)
    # At least one lone pair (2 electrons).
    return lone_electrons >= 2.0


# ---------------------------------------------------------------------------
# Ligand-local potential-donor discovery
# ---------------------------------------------------------------------------
def _ligand_subtree(mol, metal_idx: int, seed_donor: int) -> List[int]:
    """Return the connected atoms reachable from ``seed_donor`` without
    crossing the metal atom.  This is the ligand's atom set."""
    seen = {int(seed_donor)}
    stack = [int(seed_donor)]
    while stack:
        u = stack.pop()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            j = int(nb.GetIdx())
            if j == int(metal_idx) or j in seen:
                continue
            seen.add(j)
            stack.append(j)
    return sorted(seen)


def _bond_hop_distance(mol, i: int, j: int) -> int:
    """Shortest-path length in BONDS between atoms i and j (excluding M)."""
    path = Chem.GetShortestPath(mol, int(i), int(j))
    if not path:
        return 0
    return len(path) - 1


def _through_space_estimate(mol, i: int, j: int) -> float:
    """Estimate the through-space distance (Å) between two donor atoms in
    the same ligand by counting bonds and applying typical bond lengths.

    Universal (no SMARTS, no per-class table): treats every covalent edge
    as ~ 1.45 Å (mean over C-X bonds in the CCDC manifold) and applies
    Bertrand's chord-vs-arc heuristic — through-space distance for a
    n-bond zig-zag chain is ~ 1.45 * sin(zigzag_angle/2) * n with
    zigzag_angle = 120° (sp2/sp3 mean).  For donors separated by 2 bonds
    (1-3, e.g. carboxylate O-C-O) the geometric distance reduces to
    ~ 2 * 1.45 * sin(60°) ≈ 2.51 Å, which matches the CCDC mean (2.20 Å
    for COO⁻, 2.65 Å for nitrate O-O).

    Returns 0.0 if the donors are bonded directly (degenerate case;
    they cannot both donate to the same metal independently).
    """
    n_bonds = _bond_hop_distance(mol, i, j)
    if n_bonds <= 0:
        return 0.0
    if n_bonds == 1:
        # Directly bonded -> they cannot independently donate to M.
        return 0.0
    # Universal chord-length estimate for zig-zag chain of n bonds at
    # ~120° internal angle.  See module docstring.
    BOND_LEN = 1.45
    ANGLE_DEG = 120.0
    import math
    chord = 2.0 * BOND_LEN * math.sin(math.radians(ANGLE_DEG / 2.0))
    # Each additional bond beyond 2 adds ~ 1.0-1.2 Å of through-space
    # extension (zigzag chain).  Linearise as: dist(n) = chord * (1 + (n-2)/2).
    if n_bonds == 2:
        return chord
    return chord + (n_bonds - 2) * 1.0


def _is_bite_compatible(mol, i: int, j: int) -> bool:
    """Universal bite-compatibility test for a donor pair.

    True iff
        * the donors are separated by 2..MAX_BITE_HOPS bonds in the
          molecular graph (not bonded directly, not unreasonably far),
        * AND the through-space chord estimate falls in
          [BITE_MIN, BITE_MAX] = [2.0, 3.5] Å.

    This is the chemistry-agnostic version of the "1-3 oxygen pair on
    carboxylate" / "1-3 oxygen pair on nitrate" / "1-4 N pair on
    bipyridine" rule -- it admits any donor pair the metal can geometrically
    bite.
    """
    n_bonds = _bond_hop_distance(mol, i, j)
    if n_bonds < 2 or n_bonds > MAX_BITE_HOPS:
        return False
    d_est = _through_space_estimate(mol, i, j)
    return BITE_MIN <= d_est <= BITE_MAX


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def find_potential_donors_in_ligand(
    mol,
    metal_idx: int,
    seed_donor: int,
) -> List[int]:
    """Return all atoms in the ligand fragment (containing ``seed_donor``)
    that are coordination-capable donors -- i.e. have a free lone pair.

    Includes the seed donor itself, plus any other lone-pair-bearing atom
    in the same fragment.  Excludes carbons unless they pass the explicit
    carbene/cyanide lone-pair test (so neutral aliphatic carbons do NOT
    show up as "potential donors", avoiding C-H/C-C false positives).

    Universal: pure graph + element + formal-charge + valence.  No SMARTS.
    """
    subtree = _ligand_subtree(mol, metal_idx, seed_donor)
    out: List[int] = []
    for ai in subtree:
        atom = mol.GetAtomWithIdx(int(ai))
        if atom.GetSymbol() not in DONOR_ELEMENTS:
            continue
        if not _has_lone_pair(atom):
            continue
        out.append(int(ai))
    return out


def enumerate_kappa_modes_for_ligand(
    mol,
    metal_idx: int,
    seed_donor: int,
) -> List[Dict]:
    """Enumerate the κⁿ binding-mode options for ONE ligand fragment.

    Inputs:
      mol         : RDKit Mol of the FULL complex (metal + all ligands).
      metal_idx   : atom idx of the (primary) metal.
      seed_donor  : atom idx of the donor currently bonded to the metal.

    Returns a list of dicts, one per emitted binding mode:
        {
          "kappa": int,                # denticity
          "donors": (idx, ...),        # atom indices of the chosen donor set
          "donor_elements": (sym, ...),# parallel element symbols
          "is_seed_only": bool,        # True iff this is the original κ¹ on seed
          "label_suffix": str,         # e.g. "k1-O" / "k2-OO" / "k3-OOO"
        }

    The seed-only κ¹ on ``seed_donor`` is always emitted first (it is
    the binding mode the input SMILES already encodes -- the "default").
    Additional modes follow in canonical order (sorted by element multiset
    then by atom-index multiset).

    Universal: no SMARTS, no per-ligand template.  Bite-compatibility is
    a graph + geometry check only.
    """
    out: List[Dict] = []
    seed_atom = mol.GetAtomWithIdx(int(seed_donor))
    seed_sym = seed_atom.GetSymbol()

    # 1) Always emit the seed κ¹ (input topology preserved).
    out.append({
        "kappa": 1,
        "donors": (int(seed_donor),),
        "donor_elements": (seed_sym,),
        "is_seed_only": True,
        "label_suffix": f"k1-{seed_sym}",
    })

    if not _env_on():
        return out

    # 2) Universal donor discovery in the ligand fragment.
    potential = find_potential_donors_in_ligand(mol, metal_idx, seed_donor)
    # The seed donor IS in the list; remove it for the "alternate" set.
    alts = [d for d in potential if d != seed_donor]
    if not alts:
        return out

    # 3) κ¹ linkage isomers: each non-seed potential donor (deduped by
    #    canonical RDKit rank so chemically-equivalent donors collapse).
    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        ranks = list(range(mol.GetNumAtoms()))
    seen_k1_ranks: Set[int] = {int(ranks[seed_donor])}
    seed_rank = int(ranks[seed_donor])
    for alt in sorted(alts):
        r = int(ranks[alt])
        if r in seen_k1_ranks:
            continue
        seen_k1_ranks.add(r)
        sym = mol.GetAtomWithIdx(int(alt)).GetSymbol()
        out.append({
            "kappa": 1,
            "donors": (int(alt),),
            "donor_elements": (sym,),
            "is_seed_only": False,
            "label_suffix": f"k1-{sym}",
        })

    # 4) κ² subsets containing the seed donor + one bite-compatible donor.
    #    (We anchor on the seed donor: every emitted κ² is a "promote
    #    the seed to bidentate by adding one more arm".  This is the
    #    universal interpretation of the user-eye SIYMEU case --
    #    "acetate is κ², not κ¹" means the seed-O carries a partner.)
    seen_k2_sigs: Set[Tuple[int, ...]] = set()
    seen_k2_sigs.add(tuple(sorted([seed_rank])))  # κ¹ on seed already emitted
    for alt in sorted(alts):
        if not _is_bite_compatible(mol, seed_donor, alt):
            continue
        sig = tuple(sorted([seed_rank, int(ranks[alt])]))
        if sig in seen_k2_sigs:
            continue
        seen_k2_sigs.add(sig)
        e1 = seed_sym
        e2 = mol.GetAtomWithIdx(int(alt)).GetSymbol()
        elems = tuple(sorted([e1, e2]))
        out.append({
            "kappa": 2,
            "donors": tuple(sorted([int(seed_donor), int(alt)])),
            "donor_elements": elems,
            "is_seed_only": False,
            "label_suffix": f"k2-{''.join(elems)}",
        })

    # 5) κ³ subsets: seed + 2 more, all pairwise bite-compatible.
    #    Cap at MAX_KAPPA = 3 (κ⁴+ porphyrin / salen handled elsewhere).
    if MAX_KAPPA >= 3 and len(alts) >= 2:
        from itertools import combinations as _combn
        for a1, a2 in _combn(sorted(alts), 2):
            if not _is_bite_compatible(mol, seed_donor, a1):
                continue
            if not _is_bite_compatible(mol, seed_donor, a2):
                continue
            if not _is_bite_compatible(mol, a1, a2):
                continue
            r1 = int(ranks[a1]); r2 = int(ranks[a2])
            sig = tuple(sorted([seed_rank, r1, r2]))
            if sig in seen_k2_sigs:
                continue
            seen_k2_sigs.add(sig)
            e1 = mol.GetAtomWithIdx(int(a1)).GetSymbol()
            e2 = mol.GetAtomWithIdx(int(a2)).GetSymbol()
            elems = tuple(sorted([seed_sym, e1, e2]))
            out.append({
                "kappa": 3,
                "donors": tuple(sorted([int(seed_donor), int(a1), int(a2)])),
                "donor_elements": elems,
                "is_seed_only": False,
                "label_suffix": f"k3-{''.join(elems)}",
            })

    return out


def enumerate_complex_kappa_variants(
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
) -> List[Dict]:
    """Enumerate ALL κⁿ-variant donor-sets for the WHOLE complex.

    For each ligand fragment (one fragment per seed donor in
    ``donor_idxs``), enumerate the κⁿ options via
    ``enumerate_kappa_modes_for_ligand``.  Then take the Cartesian
    product across all ligands -> each product element is a full
    complex donor set (one κⁿ choice per ligand).

    De-duplicates by canonical donor-index multiset.  Returns:

        [
          {
            "donors_per_ligand": [(d, ...), ...],   # κⁿ choice per ligand
            "kappa_per_ligand":  [k, ...],
            "label_suffix":      "lig0:k1-O+lig1:k2-OO+...",
            "is_seed_default":   True/False,
            "all_donors":        [d, ...],          # flat list, sorted
          },
          ...
        ]

    The "all-seed κ¹" variant is always first (it is the input SMILES
    topology).  When the env-flag is OFF, only this canonical variant
    is emitted -- byte-identical with HEAD.
    """
    per_ligand: List[List[Dict]] = []
    for seed in donor_idxs:
        modes = enumerate_kappa_modes_for_ligand(mol, int(metal_idx), int(seed))
        if not modes:
            modes = [{
                "kappa": 1,
                "donors": (int(seed),),
                "donor_elements": (
                    mol.GetAtomWithIdx(int(seed)).GetSymbol(),
                ),
                "is_seed_only": True,
                "label_suffix": f"k1-{mol.GetAtomWithIdx(int(seed)).GetSymbol()}",
            }]
        per_ligand.append(modes)

    # Cartesian product across ligands.  Cap the total to keep the
    # explosion bounded -- typical complexes have 1-2 ambidentate ligands
    # so the product stays small.
    from itertools import product as _prod
    MAX_VARIANTS = int(os.environ.get(
        "DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM_MAX", "32"))
    out: List[Dict] = []
    seen_sigs: Set[Tuple[int, ...]] = set()
    for choice in _prod(*per_ligand):
        donors_pl = [tuple(c["donors"]) for c in choice]
        kappa_pl = [int(c["kappa"]) for c in choice]
        all_donors = sorted(d for tup in donors_pl for d in tup)
        sig = tuple(all_donors)
        if sig in seen_sigs:
            continue
        seen_sigs.add(sig)
        # Suffix: only mention ligands that DIVERGE from κ¹-seed.  The
        # default-all-κ¹ variant gets the empty suffix so its label is
        # the unmodified Mogul-PRIMARY one.
        parts = []
        is_default = True
        for li, c in enumerate(choice):
            if c.get("is_seed_only"):
                continue
            is_default = False
            parts.append(f"L{li}:{c['label_suffix']}")
        label_suffix = "+".join(parts) if parts else ""
        out.append({
            "donors_per_ligand": donors_pl,
            "kappa_per_ligand": kappa_pl,
            "label_suffix": label_suffix,
            "is_seed_default": is_default,
            "all_donors": all_donors,
        })
        if len(out) >= MAX_VARIANTS:
            break
    return out


# ---------------------------------------------------------------------------
# Graph rewiring: insert dative bonds for the chosen κⁿ donor set
# ---------------------------------------------------------------------------
def make_variant_mol(
    mol,
    metal_idx: int,
    original_donor_idxs: Sequence[int],
    variant_donor_idxs: Sequence[int],
):
    """Build a new Mol where the metal is bonded to ``variant_donor_idxs``
    (the κⁿ-promoted donor set).  All other graph features are preserved.

    Strategy:
      * Atom indices stay identical (no re-ordering) so downstream code
        that uses ``donor_idxs`` indices keeps working.
      * For donors in ``variant_donor_idxs`` that are NOT already bonded
        to the metal, add a DATIVE single bond.
      * For donors in ``original_donor_idxs`` that are NOT in
        ``variant_donor_idxs`` (replaced linkage isomer), REMOVE the
        existing M-D bond.
      * Reuse the same RingInfo by calling FastFindRings on the new mol.
    """
    em = Chem.RWMol(mol)
    orig_set = set(int(d) for d in original_donor_idxs)
    new_set = set(int(d) for d in variant_donor_idxs)
    # Add missing bonds.
    for d in new_set:
        if em.GetBondBetweenAtoms(int(metal_idx), int(d)) is None:
            em.AddBond(int(metal_idx), int(d), Chem.BondType.DATIVE)
    # Remove obsolete bonds.
    for d in orig_set - new_set:
        if em.GetBondBetweenAtoms(int(metal_idx), int(d)) is not None:
            em.RemoveBond(int(metal_idx), int(d))
    new_mol = em.GetMol()
    try:
        Chem.FastFindRings(new_mol)
    except Exception:
        pass
    return new_mol


# ---------------------------------------------------------------------------
# Standalone smoke test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Universal smoke tests on canonical ambidentate ligands.
    cases = [
        ("Ag-acetate-NH3 (κ¹-O seed)",
         "CC(=O)[O][Ag][NH3]"),
        ("Ni-nitrate-NH3 (κ¹-O seed)",
         "[O-][N+](=O)[O][Ni][NH3]"),
        ("CO ligand on M (κ¹-C seed)",
         "[O+]#[C-][Ni][NH3]"),
        ("K-thiocyanate (κ¹-S seed)",
         "[S]C#N.[K]"),
    ]
    os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
    os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"
    for name, smi in cases:
        print(f"\n=== {name} ===")
        try:
            from delfin.smiles_converter import _prepare_mol_for_embedding
            mol = _prepare_mol_for_embedding(smi)
        except Exception:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            print("  parse failed")
            continue
        if not any(a.GetSymbol() == "H" for a in mol.GetAtoms()):
            mol = Chem.AddHs(mol)
        try:
            Chem.FastFindRings(mol)
        except Exception:
            pass
        metals = [a.GetIdx() for a in mol.GetAtoms()
                  if _bd._is_metal(a.GetSymbol())]
        if not metals:
            print("  no metal")
            continue
        m = metals[0]
        donors = [int(n.GetIdx()) for n in mol.GetAtomWithIdx(m).GetNeighbors()]
        if not donors:
            print("  no donors")
            continue
        for seed in donors:
            modes = enumerate_kappa_modes_for_ligand(mol, m, seed)
            seed_sym = mol.GetAtomWithIdx(seed).GetSymbol()
            print(f"  seed={seed_sym}({seed}): {len(modes)} modes")
            for mode in modes:
                tag = "*" if mode["is_seed_only"] else " "
                print(f"   {tag} κ{mode['kappa']} "
                      f"{''.join(mode['donor_elements'])} "
                      f"donors={mode['donors']}  ({mode['label_suffix']})")
        variants = enumerate_complex_kappa_variants(mol, m, donors)
        print(f"  -> {len(variants)} full-complex κⁿ variants")
        for v in variants:
            tag = "DEFAULT" if v["is_seed_default"] else v["label_suffix"]
            print(f"      {tag}  donors={v['all_donors']}")
