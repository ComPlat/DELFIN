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
from typing import Optional, Dict, List
from rdkit import Chem
import delfin._bond_decollapse as bd

# d8 square-planar-preferring metals (else CN4 -> tetrahedral)
_D8 = {"Pt", "Pd", "Ni", "Au", "Rh", "Ir"}


def _default_geometry(metal: str, cn: int) -> Optional[str]:
    if cn == 6:
        return "OC-6 octahedron"
    if cn == 5:
        return "TBP-5 trigonal bipyramid"
    if cn == 4:
        return "SP-4 square planar" if metal in _D8 else "T-4 tetrahedron"
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
    if cn not in (4, 5, 6):
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
