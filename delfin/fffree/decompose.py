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
        # Phase G: linear coordination (Cu(I), Ag(I), Au(I), Hg(II)).
        # Simple D∞h: 2 donors 180° apart on metal axis.
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
    donor_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    cn = len(donor_idx)
    # High-CN (7-9) support is env-gated: default coverage stays CN 4-6 (byte-
    # identical), CN 7-9 polyhedra (PB-7/SQAP-8/TTP-9) only when DELFIN_FFFREE_HIGHCN=1.
    # Low-CN (3) support is env-gated: DELFIN_FFFREE_CN3=1 enables SP-3/T-3 (iter-32c).
    _allowed = set()
    _allowed.update({4, 5, 6})
    if os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    # Phase G: CN2 (linear) auto-enabled under PURE_TRACK3
    # (Cu(I)/Ag(I)/Au(I)/Hg(II) linear coordination)
    if (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_CN2", "0") == "1"):
        _allowed.add(2)
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
            # Phase G: bridging / spectator → previously hard fail. With
            # DELFIN_FFFREE_DECOMPOSE_EXTENDED=1, skip these fragments instead
            # (treat as auxiliary). Expected CCDC build-rate gain: +20%.
            if (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                or os.environ.get("DELFIN_FFFREE_DECOMPOSE_EXTENDED", "0") == "1"):
                continue
            return None                               # bridging / spectator -> legacy
        if len(fdonors) > 6:
            # Phase G: relax from > 3 to > 6 — allows tetra/penta/hexadentate
            # ligands (porphyrin κ4, EDTA κ6, salen κ4). Phase G integration.
            return None                               # >hexadentate (very rare) -> legacy
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
    #
    # Phase G (User 2026-05-31): with macrocycle.py + ring_pucker integration,
    # large polydentate ligands (porphyrin, salen) become handlable. When
    # DELFIN_FFFREE_DECOMPOSE_EXTENDED=1 (auto-on under PURE_TRACK3), raise the
    # per-donor heavy-atom cap from 8 → 20 to admit macrocyclic ligands.
    # Expected CCDC build-rate gain: +10% (porphyrin/salen/calix-class).
    _EXTENDED = (
        os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_DECOMPOSE_EXTENDED", "0") == "1"
    )
    MAX_HEAVY_PER_DONOR = 20 if _EXTENDED else 8
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
