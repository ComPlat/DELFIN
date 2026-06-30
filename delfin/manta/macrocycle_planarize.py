"""delfin.manta.macrocycle_planarize — flatten CONJUGATED metallo-macrocycles
(porphyrin / phthalocyanine / corrole) to a plane.

The analytical/embed paths planarise INDIVIDUAL aromatic rings but not the fused
MACROCYCLE as a whole, so a porphyrin comes out saddled/ruffled (measured ~0.5 A
RMS from the mean plane; a planar porphyrin is <0.1, and plain UFF/"quick" gets
~0.36 -> MANTA must not be worse).  This pass detects the conjugated macrocyclic
core (the fused aromatic ring system carrying the metal's donor atoms) and
projects it onto its best-fit plane, then restores the core bond lengths in-plane
-- yielding a flat, bond-intact macrocycle while leaving substituents (meso-aryls)
and axial ligands untouched.

GOODHART-SAFE GATES (so it NEVER flattens a genuinely 3-D macrocycle -- crown
ether / cyclam / sepulchrate, which are correctly puckered):
  * the core must be AROMATIC/conjugated (from the RDKit mol), AND
  * the core must be ALREADY mostly-planar in the build (RMS < max_rms, default
    0.8 A) -- a would-be-planar macrocycle that merely ruffled, not a 3-D one.
Both gates must pass.  Geometry-only otherwise; deterministic; never raises.

Env-gated ``DELFIN_FFFREE_MACROCYCLE_PLANARIZE`` (default ``"0"`` -> identity).
"""
from __future__ import annotations

import os
from typing import List, Optional, Sequence

import numpy as np

_COV = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "B": 0.84, "H": 0.31}
_MAX_RMS = 0.8          # only flatten cores already this planar (would-be-planar)
_MIN_CORE = 12          # a macrocycle (porphyrin core ~23); avoids small chelates


def _is_metal(sym: str, mol, idx) -> bool:
    try:
        from delfin.manta import _bond_decollapse as _bd
        return _bd._is_metal(sym)
    except Exception:
        return mol.GetAtomWithIdx(idx).GetAtomicNum() > 20


def detect_core(mol):
    """Return (core_atom_idxs, is_aromatic) for the fused conjugated ring system
    carrying a metal's donor atoms, or (None, False) if none qualifies."""
    try:
        from rdkit import Chem
        Chem.FastFindRings(mol)
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        metals = [i for i, s in enumerate(syms)
                  if mol.GetAtomWithIdx(i).GetAtomicNum() > 20 and _is_metal(s, mol, i)]
        if not metals:
            return None, False
        donor = set()
        for m in metals:
            for nb in mol.GetAtomWithIdx(m).GetNeighbors():
                donor.add(nb.GetIdx())
        rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
        core = set()
        for r in rings:
            if r & donor:
                core |= r
        for _ in range(5):                       # fuse rings sharing >=2 atoms
            for r in rings:
                if len(r & core) >= 2:
                    core |= r
        core = sorted(a for a in core if syms[a] != "H" and a not in metals)
        if len(core) < _MIN_CORE:
            return None, False
        # aromatic/conjugated gate: most core atoms aromatic OR in a >=1.5 bond
        arom = 0
        for a in core:
            at = mol.GetAtomWithIdx(a)
            if at.GetIsAromatic() or any(b.GetBondTypeAsDouble() >= 1.5 or b.GetIsAromatic()
                                          for b in at.GetBonds()):
                arom += 1
        return core, (arom >= 0.8 * len(core))
    except Exception:
        return None, False


def _core_bonds(mol, core):
    cs = set(core)
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in cs and b.GetEndAtomIdx() in cs]


def _plan_rms(P, core):
    C = P[core] - P[core].mean(0)
    _, _, V = np.linalg.svd(C)
    return float(np.sqrt(((C @ V[2]) ** 2).mean())), V[2], P[core].mean(0)


def planarize(syms: Sequence[str], coords, mol, max_rms: float = _MAX_RMS):
    """Flatten the conjugated metallo-macrocycle core onto its plane (bond-intact).
    Returns flattened coords; identity on any anomaly / gate-fail."""
    try:
        P = np.array(coords, dtype=float)
        n = len(syms)
        if P.shape != (n, 3) or not np.all(np.isfinite(P)):
            return P
        core, is_arom = detect_core(mol)
        if core is None or not is_arom:
            return P
        rms0, nrm, cen = _plan_rms(P, core)
        if rms0 > max_rms or rms0 < 0.05:
            return P                              # not a ruffled-planar macrocycle
        metals = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() > 20
                  and _is_metal(syms[i], mol, i)]
        flat = set(core) | {m for m in metals
                            if any(np.linalg.norm(P[m] - P[c]) < 2.6 for c in core)}
        cb = _core_bonds(mol, core)
        Pf = P.copy()
        for a in flat:                            # project to plane
            Pf[a] = Pf[a] - float(np.dot(Pf[a] - cen, nrm)) * nrm
        for _ in range(150):                      # in-plane bond-restore
            F = np.zeros_like(Pf)
            for i, j in cb:
                d = Pf[j] - Pf[i]; L = float(np.linalg.norm(d))
                if L < 1e-6:
                    continue
                tgt = _COV.get(syms[i], 0.8) + _COV.get(syms[j], 0.8)
                f = 0.3 * (L - tgt) * (d / L); F[i] += f; F[j] -= f
            for a in flat:
                Pf[a] = Pf[a] + F[a]
                Pf[a] = Pf[a] - float(np.dot(Pf[a] - cen, nrm)) * nrm
        rms1, _, _ = _plan_rms(Pf, core)
        if rms1 > rms0 + 1e-6 or not np.all(np.isfinite(Pf)):
            return P                              # never-worse on planarity
        return Pf
    except Exception:
        return np.array(coords, dtype=float)
