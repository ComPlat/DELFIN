"""Diagnostic: measure GRIP library hit rate for SIYMEU/BERTEB/ALAHEB.

Sweeps bonds/angles/impropers via the SAME fragment-detection logic that
grip_polish uses, and reports per-case hit ratios for a given library.

Universal — no per-class branching.
"""
from __future__ import annotations
import os, sys, json, argparse
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("DELFIN_FFFREE_MOGUL_PRIMARY", "1")

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

SMILES_SIYMEU = (
    "CC(=O)[O][Ag-][C+]1N(CC2=CC=C(C)C=C2)C(C2=CC=C(C(C)C)C=C2)"
    "=C(C2=CC=C(C(C)C)C=C2)N1CC1=CC=C(C)C=C1"
)
SMILES_BERTEB = "[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12"
SMILES_ALAHEB = (
    "CC(C)(C)[P+]1(C(C)(C)C)C=C2C=CC=C3C[P+](C(C)(C)C)(C(C)(C)C)"
    "[Ir-3]1([C]#[O+])[N]23"
)

CASES = [
    ("SIYMEU", SMILES_SIYMEU),
    ("BERTEB", SMILES_BERTEB),
    ("ALAHEB", SMILES_ALAHEB),
]


def _hyb_str(atom):
    """Match delfin.fffree.grip_fragment_detect._hyb_str."""
    h = atom.GetHybridization()
    from rdkit.Chem.rdchem import HybridizationType as H
    if h == H.SP:
        return "sp"
    if h == H.SP2:
        return "sp2"
    if h == H.SP3:
        return "sp3"
    if h == H.SP3D:
        return "sp3d"
    if h == H.SP3D2:
        return "sp3d2"
    return "*"


def _atom_is_aromatic(atom):
    try:
        return bool(atom.GetIsAromatic())
    except Exception:
        return False


def _ring_size_min(mol, idx):
    try:
        ri = mol.GetRingInfo()
        sizes = [len(r) for r in ri.AtomRings() if idx in r]
        return min(sizes) if sizes else -1
    except Exception:
        return -1


def _bond_ring_min(mol, a, b):
    try:
        ri = mol.GetRingInfo()
        sizes = []
        for r in ri.AtomRings():
            if a in r and b in r:
                sizes.append(len(r))
        return min(sizes) if sizes else -1
    except Exception:
        return -1


METALS = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu","Ac","Th","Pa","U","Np","Pu",
})


def measure_one(smiles, lib_path, tm_fallback=False):
    """Return dict with hit/total per category for ONE smiles + library."""
    if tm_fallback:
        os.environ["DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK", None)
    os.environ["DELFIN_GRIP_LIB_PATH"] = str(lib_path)

    # Late import so env is respected
    from delfin.fffree.grip_mogul_lookup import GripLibrary
    lib = GripLibrary(str(lib_path))

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    metal_idx = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALS:
            metal_idx = atom.GetIdx()
            break

    # ---- Bonds ----
    bond_q = 0
    bond_hit = 0
    bond_hit_ligand_internal = 0
    bond_q_ligand_internal = 0
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)
        z1, z2 = atom_a.GetSymbol(), atom_b.GetSymbol()
        if z1 == "H" and z2 == "H":
            continue
        is_md = (z1 in METALS) or (z2 in METALS)
        hyb1 = _hyb_str(atom_a)
        hyb2 = _hyb_str(atom_b)
        rmin = _bond_ring_min(mol, a, b)
        arom = _atom_is_aromatic(atom_a) and _atom_is_aromatic(atom_b)
        hit = lib.lookup_bond(z1, hyb1, z2, hyb2, rmin, arom, min_n=5)
        bond_q += 1
        if hit is not None:
            bond_hit += 1
        if not is_md:
            bond_q_ligand_internal += 1
            if hit is not None:
                bond_hit_ligand_internal += 1

    # ---- Angles ----
    angle_q = 0
    angle_hit = 0
    angle_q_li = 0
    angle_hit_li = 0
    for atom in mol.GetAtoms():
        nbrs = sorted(int(n.GetIdx()) for n in atom.GetNeighbors())
        if len(nbrs) < 2:
            continue
        bidx = atom.GetIdx()
        zb = atom.GetSymbol()
        hyb_b = _hyb_str(atom)
        rmin_b = _ring_size_min(mol, bidx)
        arom_b = _atom_is_aromatic(atom)
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                a_idx, c_idx = nbrs[i], nbrs[j]
                atom_a = mol.GetAtomWithIdx(a_idx)
                atom_c = mol.GetAtomWithIdx(c_idx)
                z1, z3 = atom_a.GetSymbol(), atom_c.GetSymbol()
                if z1 == "H" and z3 == "H" and zb == "H":
                    continue
                hyb1 = _hyb_str(atom_a)
                hyb3 = _hyb_str(atom_c)
                hit = lib.lookup_angle(z1, zb, hyb_b, z3, rmin_b, arom_b,
                                       hyb1=hyb1, hyb3=hyb3, min_n=5)
                angle_q += 1
                if hit is not None:
                    angle_hit += 1
                involves_metal = (z1 in METALS) or (zb in METALS) or (z3 in METALS)
                if not involves_metal:
                    angle_q_li += 1
                    if hit is not None:
                        angle_hit_li += 1

    # ---- Impropers ----
    imp_q = 0
    imp_hit = 0
    for atom in mol.GetAtoms():
        c_idx = atom.GetIdx()
        zc = atom.GetSymbol()
        if zc == "H" or zc in METALS:
            continue
        hyb_c = _hyb_str(atom)
        if hyb_c not in ("sp2", "sp3"):
            continue
        nbrs = sorted(int(n.GetIdx()) for n in atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        nbrs_data = []
        for ni in nbrs:
            a = mol.GetAtomWithIdx(ni)
            nbrs_data.append((a.GetSymbol(), _hyb_str(a), int(ni)))
        nbrs_data.sort(key=lambda t: (t[0], t[1]))
        zs_sorted = [t[0] for t in nbrs_data]
        hybs_sorted = [t[1] for t in nbrs_data]
        rmin = _ring_size_min(mol, c_idx)
        hit = lib.lookup_improper(zc, hyb_c, zs_sorted, ring_size_min=rmin,
                                   neighbor_hybs_sorted=hybs_sorted, min_n=5)
        imp_q += 1
        if hit is not None:
            imp_hit += 1

    return {
        "bond_hit": bond_hit, "bond_q": bond_q,
        "bond_hit_LI": bond_hit_ligand_internal, "bond_q_LI": bond_q_ligand_internal,
        "angle_hit": angle_hit, "angle_q": angle_q,
        "angle_hit_LI": angle_hit_li, "angle_q_LI": angle_q_li,
        "improper_hit": imp_hit, "improper_q": imp_q,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", required=True)
    parser.add_argument("--tm-fallback", action="store_true")
    args = parser.parse_args()

    print(f"Library: {args.lib}")
    print(f"TM-fallback env: {args.tm_fallback}")
    print()
    for name, smi in CASES:
        r = measure_one(smi, args.lib, args.tm_fallback)
        if r is None:
            print(f"{name}: FAILED to parse smiles")
            continue
        def _pct(num, den):
            return f"{100.0*num/max(den,1):.1f}%"
        print(f"--- {name} ---")
        print(f"  bond (ALL):  {r['bond_hit']:>4}/{r['bond_q']:<4} = {_pct(r['bond_hit'], r['bond_q']):>6}")
        print(f"  bond (LI):   {r['bond_hit_LI']:>4}/{r['bond_q_LI']:<4} = {_pct(r['bond_hit_LI'], r['bond_q_LI']):>6}  [non-metal-incident]")
        print(f"  angle (ALL): {r['angle_hit']:>4}/{r['angle_q']:<4} = {_pct(r['angle_hit'], r['angle_q']):>6}")
        print(f"  angle (LI):  {r['angle_hit_LI']:>4}/{r['angle_q_LI']:<4} = {_pct(r['angle_hit_LI'], r['angle_q_LI']):>6}")
        print(f"  improper:    {r['improper_hit']:>4}/{r['improper_q']:<4} = {_pct(r['improper_hit'], r['improper_q']):>6}")


if __name__ == "__main__":
    main()
