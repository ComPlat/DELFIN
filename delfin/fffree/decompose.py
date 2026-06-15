"""delfin.fffree.decompose — SMILES → coordination decomposition for the
metal-FF-free backend.

Scope (falls back to None otherwise): mononuclear complex with explicit
metal-donor bonds in the graph, CN 4/5/6.  Returns the metal, CN, default
geometry, and per-ligand fragments (mol + local donor indices + donor elements +
denticity) — monodentate and chelating (bi-/tridentate) ligands are both
detected; ``has_chelate`` flags whether any ligand is polydentate.  CN outside
4-6, >1 metal, >tridentate ligands, or unparseable input -> return None.
"""
from __future__ import annotations
import os
from typing import Optional, Dict, List
from rdkit import Chem
import delfin._bond_decollapse as bd

# d8 square-planar-preferring metals (else CN4 -> tetrahedral)
_D8 = {"Pt", "Pd", "Ni", "Au", "Rh", "Ir"}


def _default_geometry(metal: str, cn: int) -> Optional[str]:
    if cn == 2:
        # iter-32f (DELFIN_FFFREE_CN_EXTEND): linear two-coordinate — the canonical
        # d10 geometry (Cu(I)/Ag(I)/Au(I)/Hg(II)).  No metal-specific branch: CN2 is
        # essentially always linear.
        return "L-2 linear"
    if cn == 6:
        return "OC-6 octahedron"
    if cn == 5:
        return "TBP-5 trigonal bipyramid"
    if cn == 4:
        return "SP-4 square planar" if metal in _D8 else "T-4 tetrahedron"
    if cn == 3:
        # iter-32c (User 2026-05-28 ADUMOD: Pd CN3 was built as Td/linear).
        # d⁸ Pt/Pd/Ni/Au/Rh/Ir CN3 prefer T-shape (square-planar with one vacancy);
        # all other metals get trigonal-planar SP-3.  Both isomers can be additively
        # enumerated via DELFIN_FFFREE_DUAL_CN3=1 (mirror of dual-CN4 / dual-CN6 wiring).
        return "T-3 T-shape" if metal in _D8 else "SP-3 trigonal planar"
    if cn == 7:
        return "PB-7 pentagonal bipyramid"
    if cn == 8:
        return "SQAP-8 square antiprism"
    if cn == 9:
        return "TTP-9 tricapped trigonal prism"
    return None


def _rigid_hapto_enabled() -> bool:
    """Rigid-η construction on the FF-free path (default OFF -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_RIGID_HAPTO", "0") == "1"


def _eta_groups(mol, m: int) -> List[List[int]]:
    """Detect η (hapto) groups among the metal's carbon neighbours: maximal sets
    of ≥2 metal-bound carbons that are mutually contiguous through C-C bonds
    (a Cp / arene / diene / allyl π-face).  Graph-only, deterministic; returns a
    list of sorted carbon-index lists.  A lone metal-bound carbon (σ M-C, e.g. a
    carbonyl C or a methyl) is NOT an η group and is left out.

    Mirrors smiles_converter._find_hapto_groups but scoped to ONE metal and used
    only by the rigid-hapto FF-free path (env-gated)."""
    matom = mol.GetAtomWithIdx(m)
    c_nbrs = [n.GetIdx() for n in matom.GetNeighbors() if n.GetAtomicNum() == 6]
    cset = set(c_nbrs)
    if len(cset) < 2:
        return []
    seen: set = set()
    groups: List[List[int]] = []
    for start in c_nbrs:
        if start in seen:
            continue
        comp: List[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nb.GetIdx()
                if ni in cset and ni not in seen:
                    seen.add(ni)
                    stack.append(ni)
        if len(comp) >= 2:                            # contiguous π-face = η group
            groups.append(sorted(comp))
    return groups


def _decompose_hapto(smiles: str, mol, m: int, matom) -> Optional[Dict]:
    """Hapto-aware decomposition (DELFIN_FFFREE_RIGID_HAPTO=1 only).

    Treats each contiguous metal-bound π-face (Cp / arene / diene / allyl) as ONE
    coordination site occupying ONE polyhedron vertex (the ring centroid), so the
    effective CN is (#σ-donors) + (#η-groups) rather than the inflated η-carbon
    count.  Returns a decompose dict with per-ligand 'eta' metadata that
    assemble_complex builds as a rigid ring on the centroid vertex, or None to
    fall back to the legacy hapto path.  Universal, graph-only, deterministic."""
    egroups = _eta_groups(mol, m)
    if not egroups:
        return None                                   # no η-face -> not our case
    eta_atoms = set()
    for g in egroups:
        eta_atoms.update(g)
    nbr_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    sigma_idx = [d for d in nbr_idx if d not in eta_atoms]
    cn = len(sigma_idx) + len(egroups)                # effective coordination number
    _allowed = {4, 5, 6}
    if os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    if cn not in _allowed:
        return None
    metal = matom.GetSymbol()
    geometry = _default_geometry(metal, cn)
    if geometry is None:
        return None
    donor_elem = {d: mol.GetAtomWithIdx(d).GetSymbol() for d in sigma_idx}
    # Cleave: break every M-σ-donor bond AND every M-(η-carbon) bond, then split
    # into sanitized fragment mols (indices stay stable inside each fragment).
    em = Chem.RWMol(mol)
    for d in nbr_idx:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    mapping: List = []
    try:
        frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                 fragsMolAtomMapping=mapping)
    except Exception:
        return None
    sigma_set = set(sigma_idx)
    # map global η-atom -> its group id (for per-fragment grouping)
    eta_gid = {}
    for gid, g in enumerate(egroups):
        for a in g:
            eta_gid[a] = gid
    ligands: List[Dict] = []
    n_sites = 0
    for fmol, orig_idxs in zip(frags, mapping):
        orig = list(orig_idxs)
        if m in orig:
            continue                                  # the metal's own fragment
        f_sigma = [o for o in orig if o in sigma_set]
        f_eta_gids = sorted({eta_gid[o] for o in orig if o in eta_gid})
        if not f_sigma and not f_eta_gids:
            return None                               # bridging / spectator -> legacy
        # A fragment may carry both σ-donors and η-faces (rare); but each must map
        # cleanly onto its own vertex.  Keep v1 simple+safe: a fragment is EITHER a
        # set of σ-donors (handled like the Werner path) OR exactly ONE η-face.
        if f_eta_gids:
            if len(f_eta_gids) != 1 or f_sigma:
                return None                           # mixed/multi η in one frag -> legacy
            gid = f_eta_gids[0]
            g = egroups[gid]
            local_eta = [orig.index(o) for o in g]
            ligands.append({
                "mol": fmol,
                "is_eta": True,
                "eta_local_idxs": local_eta,          # ring carbons (local indices)
                "eta_n": len(g),                      # η hapticity (5=Cp, 6=arene, ...)
                "denticity": 1,                       # occupies ONE vertex (centroid)
                "donor_elem": "C",
            })
            n_sites += 1
        else:
            if len(f_sigma) > 3:
                return None                           # >tridentate (rare) -> legacy
            local_donors = [orig.index(o) for o in f_sigma]
            ligands.append({
                "mol": fmol,
                "is_eta": False,
                "donor_local_idx": local_donors[0],
                "donor_local_idxs": local_donors,
                "denticity": len(f_sigma),
                "donor_elem": donor_elem[f_sigma[0]],
                "donor_elems": [donor_elem[o] for o in f_sigma],
            })
            n_sites += len(f_sigma)
    if n_sites != cn:
        return None
    # Ligand-complexity gate (mirror of the Werner path), but EXEMPT η-faces: an
    # η-ring is placed rigidly as one unit, so its heavy-atom count is irrelevant.
    MAX_HEAVY_PER_DONOR = 8
    for lg in ligands:
        if lg.get("is_eta"):
            continue
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        if nheavy / max(lg["denticity"], 1) > MAX_HEAVY_PER_DONOR:
            return None
    has_chelate = any((not lg.get("is_eta")) and lg["denticity"] >= 2
                      for lg in ligands)
    has_eta = any(lg.get("is_eta") for lg in ligands)
    return {"metal": metal, "cn": cn, "geometry": geometry,
            "has_chelate": has_chelate, "has_eta": has_eta,
            "donor_elems": [lg["donor_elem"] for lg in ligands],
            "ligands": ligands}


def decompose(smiles: str) -> Optional[Dict]:
    # Reuse the converter's full organometallic mol-preparation (stk / dative-bond
    # conversion / charge+H perception) so cleaved ligands have correct chemistry
    # (NH3 stays NH3, Cl⁻ stays Cl⁻ — not HCl).  Lazy import avoids the import
    # cycle (smiles_converter -> converter_backend -> decompose).
    mol = None
    try:
        from delfin.smiles_converter import _prepare_mol_for_embedding
        mol = _prepare_mol_for_embedding(smiles)
    except Exception:
        mol = None
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    metals = [a.GetIdx() for a in mol.GetAtoms() if bd._is_metal(a.GetSymbol())]
    if len(metals) != 1:
        return None                                   # mononuclear only (v1)
    m = metals[0]
    matom = mol.GetAtomWithIdx(m)
    # Rigid-hapto path (env-gated, default OFF): if the metal carries a contiguous
    # π-face (Cp / arene / diene / allyl), collapse each face to ONE coordination
    # site so the effective CN passes the 4-6 gate and reaches the FF-free build.
    # Byte-identical when the flag is off (the branch is never entered).
    if _rigid_hapto_enabled():
        d_hap = _decompose_hapto(smiles, mol, m, matom)
        if d_hap is not None:
            return d_hap
        # no η-face (or unhandled hapto topology) -> fall through to the Werner path,
        # which is byte-identical to the flag-off behaviour for non-hapto inputs.
    donor_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    cn = len(donor_idx)
    # High-CN (7-9) support is env-gated: default coverage stays CN 4-6 (byte-
    # identical), CN 7-9 polyhedra (PB-7/SQAP-8/TTP-9) only when DELFIN_FFFREE_HIGHCN=1.
    # Low-CN (3) support is env-gated: DELFIN_FFFREE_CN3=1 enables SP-3/T-3 (iter-32c).
    # CN-coverage extension (iter-32f): DELFIN_FFFREE_CN_EXTEND=1 accepts the most
    # impactful out-of-range buckets — CN2 (linear, the big one: d10 Cu/Ag/Au/Hg) plus
    # CN7 (PB-7) and CN8 (SQAP-8) — routing them to the polyhedra that already exist.
    # Default OFF -> byte-identical (the branch is never entered).
    _allowed = set()
    _allowed.update({4, 5, 6})
    if os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    if os.environ.get("DELFIN_FFFREE_CN_EXTEND", "0") == "1":
        _allowed.update({2, 7, 8})
    if cn not in _allowed:
        return None
    metal = matom.GetSymbol()
    geometry = _default_geometry(metal, cn)
    if geometry is None:
        return None

    donor_elem = {d: mol.GetAtomWithIdx(d).GetSymbol() for d in donor_idx}
    # break metal-donor bonds (keep atom indices stable), then split into
    # sanitized fragment mols (sanitize fixes valences -> N regains its H etc.).
    em = Chem.RWMol(mol)
    for d in donor_idx:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    mapping: List = []
    try:
        frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                 fragsMolAtomMapping=mapping)
    except Exception:
        return None
    donor_set = set(donor_idx)
    ligands: List[Dict] = []
    n_chelate_bonds = 0
    for fmol, orig_idxs in zip(frags, mapping):
        orig = list(orig_idxs)
        if m in orig:
            continue                                  # the metal's own fragment
        fdonors = [o for o in orig if o in donor_set]
        if len(fdonors) == 0:
            return None                               # bridging / spectator -> legacy
        if len(fdonors) > 3:
            return None                               # >tridentate (rare) -> legacy
        local_donors = [orig.index(o) for o in fdonors]
        ligands.append({
            "mol": fmol,
            "donor_local_idx": local_donors[0],       # primary (back-compat)
            "donor_local_idxs": local_donors,         # all donors (chelate)
            "denticity": len(fdonors),
            "donor_elem": donor_elem[fdonors[0]],
            "donor_elems": [donor_elem[o] for o in fdonors],
        })
        n_chelate_bonds += len(fdonors)
    if n_chelate_bonds != cn:
        return None
    has_chelate = any(lg["denticity"] >= 2 for lg in ligands)
    # Ligand-complexity gate, measured PER DONOR ARM (heavy atoms / denticity):
    # small/simple ligands place cleanly; very large conjugated ligands need
    # conformer work the rigid placement doesn't do, so bail to None for those.
    # Per-arm (not per-ligand) so a polydentate chelate (e.g. citrate, 13 heavy
    # over 6 donors ~ 2/arm) is not unfairly rejected for its total size.
    MAX_HEAVY_PER_DONOR = 8
    for lg in ligands:
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        if nheavy / max(lg["denticity"], 1) > MAX_HEAVY_PER_DONOR:
            return None
    return {"metal": metal, "cn": cn, "geometry": geometry,
            "has_chelate": has_chelate,
            "donor_elems": [lg["donor_elem"] for lg in ligands],
            "ligands": ligands}


if __name__ == "__main__":
    for label, smi in [("cisplatin", "N[Pt](N)(Cl)Cl"),
                       ("hexammineCo", "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]"),
                       ("[CoCl3(NH3)3]", "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]")]:
        d = decompose(smi)
        if d:
            print(f"{label:<16} metal={d['metal']} CN={d['cn']} geom={d['geometry']} "
                  f"donors={d['donor_elems']}")
        else:
            print(f"{label:<16} -> None (legacy fallback)")
