"""delfin.fffree.decompose — SMILES → coordination decomposition for the
metal-FF-free backend.

v1 scope (safe, falls back to legacy otherwise): mononuclear complex with explicit
metal-donor bonds in the graph, all-monodentate donors, CN 4/5/6.  Returns the
metal, per-donor ligand fragments (mol + local donor index + element), CN and a
default geometry.  Chelates (ligand with >1 metal bond), dative/dot-separated
representations, hapto, and CN outside 4-6 -> return None (legacy handles them).
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
    for fmol, orig_idxs in zip(frags, mapping):
        orig = list(orig_idxs)
        if m in orig:
            continue                                  # the metal's own fragment
        fdonors = [o for o in orig if o in donor_set]
        if len(fdonors) != 1:
            return None                               # chelate / bridging -> legacy
        local_donor = orig.index(fdonors[0])
        ligands.append({"mol": fmol, "donor_local_idx": local_donor,
                        "donor_elem": donor_elem[fdonors[0]]})
    if len(ligands) != cn:
        return None
    return {"metal": metal, "cn": cn, "geometry": geometry,
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
