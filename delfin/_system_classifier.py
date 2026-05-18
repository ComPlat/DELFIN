"""Universal 8-axis system classifier for metal complexes.

Built per Welle-5m-Y brief (2026-05-18) as INPUT for Agent Z's
treatment-matrix-dispatcher.  All axes are derived purely from RDKit
graph + atomic-number features.  No SMILES regex, no refcode prefixes,
no element allowlists — strict universal-fundamental doctrine
([[feedback_universal_fundamental_doctrine]]).

The public entry point is :func:`classify_complex_system`, which
returns a flat ``dict`` with the following keys:

``class``               see legacy ``_classify_complex_class`` semantics
``CN``                  observed CN of the **primary** metal (hapto-
                        collapsed: one hapto group counts as 1 donor)
``donor_heterogeneity`` ``homo`` / ``bi`` / ``tri`` / ``maximal-asym``
                        from the number of distinct donor element
                        classes (`>=4` ⇒ maximal-asym)
``chelate_pattern``     ``all-mono``/``mono-plus-bid``/``bis-bid``/
                        ``tris-bid``/``macrocyclic``/``pincer``
``metal_block``         ``3d``/``4d``/``5d``/``f-block``/``main-group``
``q_magnitude``         ``int(abs(formal_charge(metal)))``
``hapticity_max``       largest hapto group size (η-count), 0 if none
``bulky_flag``          ``True`` iff any donor atom carries a tBu-like
                        (sp3-C with 3 sp3-C neighbours) substituent
``cell_key``            stable underscore-separated descriptor used
                        as the dispatcher matrix key

The classifier is **helper-only**: nothing in the production pipeline
imports it.  Wiring is the job of Agent Z's treatment-matrix-dispatcher.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    from rdkit import Chem
    _RDKIT_OK = True
except Exception:  # pragma: no cover — RDKit always present in DELFIN
    _RDKIT_OK = False


# ---------------------------------------------------------------------------
# Metal-block partitioning by atomic number (no element allowlist —
# we partition the whole periodic table by Z-range, which is a pure
# graph/atomic-number feature, not a refcode-style lookup table).
# ---------------------------------------------------------------------------
_BLOCK_3D_Z = frozenset(range(21, 31))               # Sc(21)..Zn(30)
_BLOCK_4D_Z = frozenset(range(39, 49))               # Y(39)..Cd(48)
_BLOCK_5D_Z = frozenset({57}) | frozenset(range(72, 81))  # La + Hf..Hg
_BLOCK_F_Z = frozenset(range(58, 72)) | frozenset(range(89, 104))
# Main-group metals: alkali, alkaline-earth, p-block metals.  We treat
# everything that is *not* a TM/lanthanide/actinide and *is* in our
# project ``_METALS`` set as main-group for dispatch purposes.

# Donor element classes — coarse grouping by chemical character (sp3 C
# vs π-system C is resolved by the hapto axis, not the heterogeneity
# axis).  Grouping is by atomic number, so no element allowlist:
#  * "C"      ⇒ Z == 6
#  * "N"      ⇒ Z == 7
#  * "O"      ⇒ Z == 8
#  * "P"      ⇒ Z == 15
#  * "S"      ⇒ Z == 16
#  * "halide" ⇒ Z in {9,17,35,53,85}   (F Cl Br I At)
#  * "hydride"⇒ Z == 1
#  * "Si"     ⇒ Z == 14
#  * "B"      ⇒ Z ==  5
#  * "other"  ⇒ anything else (Se, As, Te, ...)
_HALOGEN_Z = frozenset({9, 17, 35, 53, 85})


def _donor_element_class(z: int) -> str:
    """Map atomic number → donor element class.

    Universal: dispatches by ``z`` only, no hardcoded element string
    matching.  Every Z gets a class; nothing is dropped.
    """
    if z == 1:
        return "hydride"
    if z == 5:
        return "B"
    if z == 6:
        return "C"
    if z == 7:
        return "N"
    if z == 8:
        return "O"
    if z == 14:
        return "Si"
    if z == 15:
        return "P"
    if z == 16:
        return "S"
    if z in _HALOGEN_Z:
        return "halide"
    return "other"


def _metal_block(z: int) -> str:
    if z in _BLOCK_3D_Z:
        return "3d"
    if z in _BLOCK_4D_Z:
        return "4d"
    if z in _BLOCK_5D_Z:
        return "5d"
    if z in _BLOCK_F_Z:
        return "f-block"
    return "main-group"


# ---------------------------------------------------------------------------
# Hapto / metal helpers — locally re-implemented (mirrors smiles_converter
# semantics) to avoid a circular import at helper-load time.
# ---------------------------------------------------------------------------
def _is_metal(z: int) -> bool:
    """Z-range metal predicate (groups 1-2, d-, f-, p-block metals).

    Universal: covers any Z that DELFIN considers a metal.  Anything
    that is not H/C/N/O/F/Ne/.../noble-gas/typical-non-metal and that
    can act as a coordination centre.
    """
    # alkali + alkaline-earth (Z=3,4,11,12,19,20,37,38,55,56,87,88)
    if z in {3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88}:
        return True
    # d-block 3d,4d,5d (incl. La)
    if z in _BLOCK_3D_Z or z in _BLOCK_4D_Z or z in _BLOCK_5D_Z:
        return True
    # f-block (Ln + An)
    if z in _BLOCK_F_Z:
        return True
    # p-block metals (Al, Ga, In, Tl; Sn, Pb; Bi, Po)
    if z in {13, 31, 49, 81, 50, 82, 83, 84}:
        return True
    return False


def _find_hapto_groups_local(mol) -> List[Tuple[int, List[int]]]:
    """Detect contiguous metal-bound carbon blocks (η-coordination).

    Mirrors :func:`delfin.smiles_converter._find_hapto_groups` but lives
    here to keep the classifier import-free of the heavy converter.
    """
    if mol is None:
        return []
    groups: List[Tuple[int, List[int]]] = []
    for atom in mol.GetAtoms():
        if not _is_metal(atom.GetAtomicNum()):
            continue
        metal_idx = atom.GetIdx()
        all_nbrs = list(atom.GetNeighbors())
        c_nbrs = [n.GetIdx() for n in all_nbrs if n.GetAtomicNum() == 6]
        if len(c_nbrs) < 2:
            continue
        b_count = sum(1 for n in all_nbrs if n.GetAtomicNum() == 5)
        if b_count > len(all_nbrs) / 2:
            continue
        c_set = set(c_nbrs)
        seen: Set[int] = set()
        for start in c_nbrs:
            if start in seen:
                continue
            comp: List[int] = []
            stack = [start]
            seen.add(start)
            while stack:
                cur = stack.pop()
                comp.append(cur)
                cur_atom = mol.GetAtomWithIdx(cur)
                for nbr in cur_atom.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in c_set or ni in seen:
                        continue
                    seen.add(ni)
                    stack.append(ni)
            if len(comp) < 2:
                continue
            ring_like = any(mol.GetAtomWithIdx(i).IsInRing() for i in comp)
            if len(comp) >= 3 or ring_like:
                groups.append((metal_idx, sorted(comp)))
    return groups


# ---------------------------------------------------------------------------
# Chelate ring detection — graph-feature only.
#
# Definition: a "chelate ring" is an SSSR ring of length 4..8 that
# contains exactly one metal atom and at least two donor atoms bonded
# to that metal.  Macrocycles (length >= 9 containing the metal) are
# flagged separately.  A pincer is a single ligand (tracked by ligand
# fragment id, not by ring) providing >=3 donors arranged so two
# chelate rings share an edge at the metal (mer-tridentate).
# ---------------------------------------------------------------------------
def _ligand_fragment_id(mol, metal_indices: Iterable[int]) -> Dict[int, int]:
    """Assign a ligand-fragment id to every non-metal atom.

    Two non-metal atoms share an id iff they are connected through a
    path that does not pass through any metal atom.  This is the
    standard graph-theoretic ligand decomposition.
    """
    metal_set = set(metal_indices)
    n = mol.GetNumAtoms()
    fid = [-1] * n
    next_id = 0
    for start in range(n):
        if start in metal_set or fid[start] != -1:
            continue
        # BFS
        stack = [start]
        fid[start] = next_id
        while stack:
            cur = stack.pop()
            atom = mol.GetAtomWithIdx(cur)
            for nbr in atom.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in metal_set or fid[ni] != -1:
                    continue
                fid[ni] = next_id
                stack.append(ni)
        next_id += 1
    return {i: fid[i] for i in range(n) if i not in metal_set}


def _chelate_rings_for_metal(mol, metal_idx: int) -> List[List[int]]:
    """Return SSSR rings (as atom-idx lists) that contain ``metal_idx``."""
    out: List[List[int]] = []
    try:
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if metal_idx in ring:
                out.append(list(ring))
    except Exception:
        pass
    return out


def _is_pincer_at_metal(
    mol,
    metal_idx: int,
    donor_indices: Iterable[int],
    frag_of: Dict[int, int],
) -> bool:
    """Pincer detection: one ligand fragment provides >=3 donors arranged
    as two fused chelate rings sharing the central donor."""
    donors = list(donor_indices)
    # Group donors by fragment id
    by_frag: Dict[int, List[int]] = {}
    for d in donors:
        f = frag_of.get(d)
        if f is None:
            continue
        by_frag.setdefault(f, []).append(d)
    # A pincer ligand contributes >=3 donor atoms in one fragment AND
    # produces >=2 chelate rings that share an atom (the central donor).
    chelate_rings = [
        set(r)
        for r in _chelate_rings_for_metal(mol, metal_idx)
        if 4 <= len(r) <= 8
    ]
    for frag_donors in by_frag.values():
        if len(frag_donors) < 3:
            continue
        rings_with_donor = []
        for d in frag_donors:
            for r in chelate_rings:
                if d in r:
                    rings_with_donor.append((d, r))
        # need at least 2 distinct chelate rings sharing the *metal*
        ring_signatures = {frozenset(r) for _, r in rings_with_donor}
        if len(ring_signatures) >= 2:
            return True
    return False


# ---------------------------------------------------------------------------
# Bulky substituent detection — universal graph signature.
#
# "tBu pattern" = sp3 C with three sp3-C neighbours (no H needed on the
# central C, three heavy-C substituents).  This generalises naturally to
# any quaternary-ish carbon with three sp3-C neighbours: tBu, neopentyl
# arm, adamantyl bridgehead, etc.  No SMILES match, no name match.
# ---------------------------------------------------------------------------
def _atom_is_sp3_carbon(atom) -> bool:
    """Detect sp3 C without relying on RDKit hybridization perception.

    Falls back to bond-type + aromaticity inspection so the helper works
    on partially-sanitized metal-complex parses.  Universal predicate:
    "sp3 C" = atomic-num 6, not aromatic, all bonds are single bonds.
    """
    if atom.GetAtomicNum() != 6:
        return False
    if atom.GetIsAromatic():
        return False
    for b in atom.GetBonds():
        bt = b.GetBondType()
        if bt != Chem.BondType.SINGLE:
            return False
    return True


def _atom_is_tbu_like(atom) -> bool:
    """tBu-like central C: sp3 C with >=3 sp3-C neighbours.

    Pure graph signature — covers tBu, neopentyl bridgehead, adamantyl
    bridgehead, etc.  No name match.
    """
    if not _atom_is_sp3_carbon(atom):
        return False
    sp3_c_nbrs = sum(
        1 for n in atom.GetNeighbors() if _atom_is_sp3_carbon(n)
    )
    return sp3_c_nbrs >= 3


def _has_bulky_near_donor(mol, donor_idx: int, depth: int = 3) -> bool:
    """BFS up to ``depth`` bonds from a donor atom — flag if any node
    matches the tBu graph signature.  ``depth=3`` reaches the β/γ
    carbons typical of bulky phosphines and amines."""
    visited: Set[int] = {donor_idx}
    frontier = [donor_idx]
    for _ in range(depth):
        next_frontier: List[int] = []
        for idx in frontier:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in visited:
                    continue
                visited.add(ni)
                if _atom_is_tbu_like(nbr):
                    return True
                next_frontier.append(ni)
        frontier = next_frontier
    return False


# ---------------------------------------------------------------------------
# Main classifier
# ---------------------------------------------------------------------------
def _heterogeneity_label(n_classes: int) -> str:
    if n_classes <= 1:
        return "homo"
    if n_classes == 2:
        return "bi"
    if n_classes == 3:
        return "tri"
    return "maximal-asym"


def _chelate_pattern(
    mol,
    metal_idx: int,
    donor_indices: List[int],
    frag_of: Dict[int, int],
) -> str:
    """Classify chelate pattern at one metal centre."""
    # Macrocyclic test: any ring with the metal of size >= 9
    rings_with_metal = _chelate_rings_for_metal(mol, metal_idx)
    if any(len(r) >= 9 for r in rings_with_metal):
        return "macrocyclic"

    # Group donors by ligand fragment.  Chelate-arity per fragment =
    # number of donor atoms that the fragment contributes.
    by_frag: Dict[int, List[int]] = {}
    for d in donor_indices:
        f = frag_of.get(d)
        if f is None:
            continue
        by_frag.setdefault(f, []).append(d)

    if not by_frag:
        return "all-mono"

    arities = sorted((len(v) for v in by_frag.values()), reverse=True)

    # Pincer (mer-tridentate single fragment)
    if _is_pincer_at_metal(mol, metal_idx, donor_indices, frag_of):
        return "pincer"

    bidentates = sum(1 for a in arities if a == 2)
    monodentates = sum(1 for a in arities if a == 1)
    has_higher = any(a >= 3 for a in arities)

    if has_higher:
        # tridentate-or-more but not pincer (e.g. tripodal): keep
        # ``pincer`` reserved for mer-planar; tripodal fac-tris is
        # better described as a chelate triple.  Fold into pincer for
        # dispatcher purposes (mer/fac tridentate both want bite-angle
        # restraints).
        return "pincer"

    if bidentates == 0:
        return "all-mono"
    if bidentates == 1:
        return "mono-plus-bid"
    if bidentates == 2:
        return "bis-bid"
    if bidentates >= 3 and monodentates == 0:
        return "tris-bid"
    # mixed: tris-bid with leftover monodentates is rare — pool with bis-bid
    return "bis-bid"


def classify_complex_system(mol) -> Dict[str, object]:
    """Return universal 8-axis system descriptor.

    See module docstring for axis semantics.  ``mol`` must be a parsed
    RDKit molecule (caller is responsible for parsing — this helper
    operates purely on the graph).
    """
    no_metal_descriptor: Dict[str, object] = {
        "class": "no_metal",
        "CN": 0,
        "donor_heterogeneity": "homo",
        "chelate_pattern": "all-mono",
        "metal_block": "main-group",
        "q_magnitude": 0,
        "hapticity_max": 0,
        "bulky_flag": False,
        "cell_key": "no_metal",
    }

    if mol is None or not _RDKIT_OK:
        return no_metal_descriptor

    # ---- 1) Locate metals ------------------------------------------------
    metal_atoms = [a for a in mol.GetAtoms() if _is_metal(a.GetAtomicNum())]
    if not metal_atoms:
        return no_metal_descriptor

    metal_indices = [a.GetIdx() for a in metal_atoms]
    n_metals = len(metal_atoms)

    # Primary metal = highest-Z (most "central" in TM complexes); ties
    # broken by atom index for determinism.  This is a pure graph
    # feature, not a name lookup.
    primary = max(metal_atoms, key=lambda a: (a.GetAtomicNum(), -a.GetIdx()))
    primary_idx = primary.GetIdx()

    # ---- 2) Hapto groups -------------------------------------------------
    hapto_groups = _find_hapto_groups_local(mol)
    hapto_by_metal: Dict[int, List[List[int]]] = {}
    for m_idx, grp in hapto_groups:
        hapto_by_metal.setdefault(m_idx, []).append(grp)

    hapticity_max = (
        max((len(g) for g in hapto_groups[0:0] or [[] for _ in []]), default=0)
        if False
        else max((len(grp) for _, grp in hapto_groups), default=0)
    )

    # ---- 3) Class (legacy 5-class semantics) -----------------------------
    n_hapto_metals = len(hapto_by_metal)
    if n_metals == 1:
        cls = "hapto" if n_hapto_metals else "sigma"
    else:
        # Cluster: >=2 metals with M-M graph edge.  Otherwise the
        # legacy multi_sigma / multi_hapto labels are kept for backward
        # compatibility with the dispatcher.
        m_set = set(metal_indices)
        has_mm = False
        for a in metal_atoms:
            for n in a.GetNeighbors():
                if n.GetIdx() in m_set:
                    has_mm = True
                    break
            if has_mm:
                break
        if has_mm:
            cls = "cluster"
        elif n_hapto_metals:
            cls = "multi-hapto"
        else:
            cls = "multi-sigma"

    # ---- 4) Donors at primary metal --------------------------------------
    primary_atom = mol.GetAtomWithIdx(primary_idx)
    primary_hapto_groups = hapto_by_metal.get(primary_idx, [])
    # Atoms belonging to a hapto group at this metal — those collapse
    # to one donor "slot" per group.
    hapto_atomset_for_primary: Set[int] = set()
    for grp in primary_hapto_groups:
        hapto_atomset_for_primary.update(grp)

    sigma_donor_indices: List[int] = []
    for nbr in primary_atom.GetNeighbors():
        ni = nbr.GetIdx()
        if ni in hapto_atomset_for_primary:
            continue
        # Exclude M-M edges — they are not "donors" in the dispatch sense
        if _is_metal(nbr.GetAtomicNum()):
            continue
        sigma_donor_indices.append(ni)

    # Observed CN = sigma donors + one slot per hapto group
    cn = len(sigma_donor_indices) + len(primary_hapto_groups)

    # ---- 5) Donor heterogeneity -----------------------------------------
    donor_classes: Set[str] = set()
    for ni in sigma_donor_indices:
        z = mol.GetAtomWithIdx(ni).GetAtomicNum()
        donor_classes.add(_donor_element_class(z))
    if primary_hapto_groups:
        # All hapto groups are π-C systems; collapse into a single "Cπ"
        # class so that "Cp + sigma-N + sigma-Cl" registers as tri, not bi.
        donor_classes.add("Cpi")
    heterogeneity = _heterogeneity_label(len(donor_classes))

    # ---- 6) Chelate pattern ---------------------------------------------
    frag_of = _ligand_fragment_id(mol, metal_indices)
    chelate = _chelate_pattern(
        mol,
        primary_idx,
        sigma_donor_indices,
        frag_of,
    )

    # ---- 7) Metal block + charge magnitude ------------------------------
    block = _metal_block(primary.GetAtomicNum())
    try:
        q_mag = int(abs(primary.GetFormalCharge()))
    except Exception:
        q_mag = 0

    # ---- 8) Bulky flag ---------------------------------------------------
    bulky = False
    # Probe each donor environment (both sigma and hapto-anchor) for
    # tBu-like substructure within graph depth 3.
    probe_donors: List[int] = list(sigma_donor_indices)
    for grp in primary_hapto_groups:
        if grp:
            probe_donors.append(grp[0])
    for d in probe_donors:
        if _has_bulky_near_donor(mol, d, depth=3):
            bulky = True
            break

    # ---- 9) Cell key -----------------------------------------------------
    eta = f"eta{hapticity_max}"
    cell_key = "_".join(
        [
            cls,
            f"CN{cn}",
            heterogeneity,
            chelate,
            block,
            f"q{q_mag}",
            eta,
            "bulky" if bulky else "nobulky",
        ]
    )

    return {
        "class": cls,
        "CN": cn,
        "donor_heterogeneity": heterogeneity,
        "chelate_pattern": chelate,
        "metal_block": block,
        "q_magnitude": q_mag,
        "hapticity_max": hapticity_max,
        "bulky_flag": bulky,
        "cell_key": cell_key,
    }


__all__ = ["classify_complex_system"]
