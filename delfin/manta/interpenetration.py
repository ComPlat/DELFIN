"""delfin.manta.interpenetration -- contact census of a placed legacy topology frame (Ursache 5, 27.10.2026).

A legacy topology frame is assembled from rigid ligand fragments (template conformer, Procrustes onto the
polyhedron vertices, spin about the donor axis) and only THEN relaxed by UFF.  When two placed fragments
interpenetrate, a local force field cannot pull them apart: it compresses bonds instead.  This module
measures, on the SMILES graph of ``mol`` (atom order of the frame), how deep the placed fragments overlap:

* ``inter``: the closest pair of atoms in DIFFERENT non-metal fragments, as fraction of the Bondi vdW sum,
* ``intra``: the closest pair inside one fragment at topological distance >= 3, same fraction,
* ``bond``:  the shortest non-metal SMILES bond, as fraction of its bond-order-aware covalent reference.

Pure geometry on public radii (Bondi 1964; RDKit covalent radii).  No CSD/CCDC data.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Bondi (1964) van der Waals radii, A; Se/Te/As/B/Si/Ge/Sb after Mantina et al. 2009.
BONDI_VDW: Dict[str, float] = {
    "H": 1.20, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Si": 2.10, "P": 1.80, "S": 1.80,
    "Cl": 1.75, "As": 1.85, "Se": 1.90, "Br": 1.85, "Te": 2.06, "I": 1.98, "Ge": 2.11, "Sb": 2.06,
}
_VDW_DEFAULT = 1.80


def u5_trace_enabled() -> bool:
    """DELFIN_U5_TRACE=1 (default off, output only)."""
    return os.environ.get("DELFIN_U5_TRACE", "0") == "1"


def _u5_fragments(mol, metal_set) -> Tuple[List[int], Dict[int, int]]:
    """Non-metal atoms and their connected-component label (smallest index of the component)."""
    nm = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() not in metal_set]
    lab: Dict[int, int] = {}
    for s in nm:
        if s in lab:
            continue
        stack = [s]
        while stack:
            x = stack.pop()
            if x in lab:
                continue
            lab[x] = s
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                if nb.GetSymbol() not in metal_set and nb.GetIdx() not in lab:
                    stack.append(nb.GetIdx())
    return nm, lab


def _u5_topo_within(mol, metal_set, idx: Sequence[int], maxd: int = 2) -> set:
    """Pairs (i<j) of non-metal atoms at topological distance <= maxd through the non-metal graph."""
    near = set()
    for s in idx:
        front = {s}
        seen = {s}
        for _ in range(maxd):
            nxt = set()
            for x in front:
                for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                    j = nb.GetIdx()
                    if nb.GetSymbol() in metal_set or j in seen:
                        continue
                    seen.add(j)
                    nxt.add(j)
            front = nxt
        for j in seen:
            if j != s:
                near.add((min(s, j), max(s, j)))
    return near


def contact_census(mol, coords, metal_set, donor_indices: Optional[Sequence[int]] = None) -> dict:
    """Closest inter-fragment / intra-fragment (topological distance >= 3) contact and the shortest bond."""
    from rdkit import Chem
    X = np.asarray(coords, dtype=float)
    nm, lab = _u5_fragments(mol, metal_set)
    el = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
    donors = set(donor_indices or [])
    dent: Dict[int, int] = {}
    for i in nm:
        dent[lab[i]] = dent.get(lab[i], 0) + (1 if i in donors else 0)
    ids = np.array(nm, dtype=int)
    P = X[ids]
    R = np.array([BONDI_VDW.get(el[i], _VDW_DEFAULT) for i in ids])
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    Q = D / (R[:, None] + R[None, :])
    F = np.array([lab[i] for i in ids])
    same = F[:, None] == F[None, :]
    iu = np.triu_indices(len(ids), k=1)
    out = {"inter": (9.0, ""), "intra": (9.0, ""), "bond": (9.0, "")}
    qi = np.where(same[iu], 9.0, Q[iu])
    if qi.size and float(qi.min()) < 9.0:
        k = int(np.argmin(qi))
        a, b = int(ids[iu[0][k]]), int(ids[iu[1][k]])
        out["inter"] = (float(qi[k]), "%s%d-%s%d %.2f dent %d/%d" % (
            el[a], a, el[b], b, float(D[iu[0][k], iu[1][k]]), dent[lab[a]], dent[lab[b]]))
    near = _u5_topo_within(mol, metal_set, nm, 2)
    best, btxt = 9.0, ""
    for n in range(len(iu[0])):
        r, c = iu[0][n], iu[1][n]
        if not same[r, c]:
            continue
        a, b = int(ids[r]), int(ids[c])
        if (min(a, b), max(a, b)) in near:
            continue
        q = float(Q[r, c])
        if q < best:
            best, btxt = q, "%s%d-%s%d %.2f" % (el[a], a, el[b], b, float(D[r, c]))
    out["intra"] = (best, btxt)
    pt = Chem.GetPeriodicTable()
    sc = {Chem.BondType.SINGLE: 1.0, Chem.BondType.DOUBLE: 0.87, Chem.BondType.TRIPLE: 0.78,
          Chem.BondType.AROMATIC: 0.91}
    for bd in mol.GetBonds():
        a, b = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
        if el[a] in metal_set or el[b] in metal_set:
            continue
        ref = (pt.GetRcovalent(mol.GetAtomWithIdx(a).GetAtomicNum())
               + pt.GetRcovalent(mol.GetAtomWithIdx(b).GetAtomicNum())) * sc.get(bd.GetBondType(), 1.0)
        d = float(np.linalg.norm(X[a] - X[b]))
        if d / ref < out["bond"][0]:
            out["bond"] = (d / ref, "%s%d-%s%d %.2f" % (el[a], a, el[b], b, d))
    return out


def u5_trace(tag: str, xyz, mol, metal_set, donor_indices=None, extra: str = "") -> None:
    """One [U5] line per frame and stage (DELFIN_U5_TRACE=1, output only)."""
    if not u5_trace_enabled():
        return
    try:
        if hasattr(xyz, "GetAtomPosition"):
            P = [tuple(xyz.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        elif isinstance(xyz, str):
            P = []
            for ln in xyz.strip().splitlines():
                p = ln.split()
                if len(p) >= 4:
                    P.append((float(p[1]), float(p[2]), float(p[3])))
        else:
            P = [tuple(r) for r in xyz]
        if len(P) != mol.GetNumAtoms():
            os.write(2, ("[U5] %s natoms-mismatch %s\n" % (tag, extra)).encode())
            return
        c = contact_census(mol, P, metal_set, donor_indices)
        os.write(2, ("[U5] %s inter %.2f (%s) intra %.2f (%s) bond %.2f (%s) %s\n" % (
            tag, c["inter"][0], c["inter"][1], c["intra"][0], c["intra"][1], c["bond"][0], c["bond"][1],
            extra)).encode())
    except Exception as exc:  # a trace line must never kill a build
        os.write(2, ("[U5] %s trace-failed %s\n" % (tag, str(exc)[:80])).encode())
