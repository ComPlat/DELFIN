#!/usr/bin/env python3
"""conformer_enum.py — deterministic combinatorial conformer
enumeration (Confab-style, NOT stochastic global search).  Replaces SEARCH with
ENUMERATION: rotatable bonds -> discrete torsion-state product -> clash-prune ->
LOCAL MMFF relax -> RMSD/energy dedup -> keep within the relevant energy window.

This is the Layer-3 generator; completeness here is "all relevant conformers
within the energy window, up to torsion resolution" (validated by convergence:
finer grid stops adding distinct minima).  Organic ligand part only (MMFF valid);
the metal sphere is built separately.

Deterministic: ETKDG initial embed with fixed seed + numThreads=1; fixed torsion
grid; MMFF; canonical-RMSD dedup with deterministic ordering.
"""
from __future__ import annotations
import itertools
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

import os
SEED = 42
TORSIONS = (60.0, 180.0, 300.0)   # gauche+/anti/gauche- ; resolution = 120 deg
RMSD_DEDUP = 0.5                  # A; conformers closer than this are the same
ENERGY_WINDOW = 5.0               # kcal/mol above the global min = "relevant"
MAX_ROT = 7                       # cap combinatorial blow-up (3^7 = 2187)

# Track 3 pure FF-free: replace MMFF with structural-quality score.
# DELFIN_FFFREE_PURE_TRACK3=1 → use defect-count score instead of MMFF energy.
# Universal, deterministic, no force field.
_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"


def _rotatable_bonds(mol) -> List[Tuple[int, int]]:
    """Non-ring single bonds whose BOTH atoms have >=2 heavy neighbours (excludes
    terminal methyl/halide spins = no distinct conformer)."""
    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        h1 = sum(1 for n in a1.GetNeighbors() if n.GetAtomicNum() > 1)
        h2 = sum(1 for n in a2.GetNeighbors() if n.GetAtomicNum() > 1)
        if h1 >= 2 and h2 >= 2:
            out.append((a1.GetIdx(), a2.GetIdx()))
    return out


def _dihedral_ref(mol, b1: int, b2: int):
    a1 = next((n.GetIdx() for n in mol.GetAtomWithIdx(b1).GetNeighbors()
               if n.GetIdx() != b2 and n.GetAtomicNum() > 1), None)
    a3 = next((n.GetIdx() for n in mol.GetAtomWithIdx(b2).GetNeighbors()
               if n.GetIdx() != b1 and n.GetAtomicNum() > 1), None)
    return a1, a3


def _mmff_energy(mol) -> float:
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    return ff.CalcEnergy()


# Pure Track 3 FF-free: structural-quality score replaces MMFF energy
# Universal, no template, classical mechanics principles only.
_VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10, "Se": 1.90,
        "As": 1.85, "Te": 2.06}
_COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
        "Cl": 1.02, "Br": 1.20, "I": 1.39, "B": 0.84, "Si": 1.11, "Se": 1.20,
        "As": 1.19}
_CLASH_FACTOR = 0.70


def _structural_quality_score(mol) -> float:
    """Pure FF-free defect-count score (lower = better conformer).

    Counts: (a) bonded pairs with length deviation^2 from ideal covalent sum,
    (b) non-bonded close contacts < 0.7 * (vdW_i + vdW_j). Universal, deterministic,
    no force field. Used when DELFIN_FFFREE_PURE_TRACK3=1.
    """
    import numpy as _np
    n = mol.GetNumAtoms()
    conf = mol.GetConformer()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n)]
    coords = _np.array([[conf.GetAtomPosition(i).x,
                          conf.GetAtomPosition(i).y,
                          conf.GetAtomPosition(i).z] for i in range(n)])
    bonded = set()
    score = 0.0
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonded.add((min(i, j), max(i, j)))
        ideal = _COV.get(syms[i], 1.0) + _COV.get(syms[j], 1.0)
        d = float(_np.linalg.norm(coords[i] - coords[j]))
        score += (d - ideal) ** 2
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            d = float(_np.linalg.norm(coords[i] - coords[j]))
            clash_threshold = _CLASH_FACTOR * (_VDW.get(syms[i], 1.7) + _VDW.get(syms[j], 1.7))
            if d < clash_threshold:
                score += (clash_threshold - d) ** 2 * 10.0  # clash penalty
    return score


def _conformer_score(mol) -> float:
    """Dispatch based on DELFIN_FFFREE_PURE_TRACK3 env-flag."""
    if _PURE_TRACK3:
        return _structural_quality_score(mol)
    return _mmff_energy(mol)


def _optimize_or_skip(mol):
    """MMFF-optimize (Track 1/2) or skip (Track 3 pure FF-free)."""
    if _PURE_TRACK3:
        return  # rely on ETKDG starting coords; no FF post-process
    AllChem.MMFFOptimizeMolecule(mol)


def enumerate_conformers(smiles: str):
    """Return (n_combos, conformers) where conformers = list of (energy, mol) within
    the energy window, deduped by RMSD.  Deterministic."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, []
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
        return 0, []
    _optimize_or_skip(mol)
    rot = _rotatable_bonds(mol)
    refs = [( _dihedral_ref(mol, b1, b2), b1, b2) for b1, b2 in rot]
    refs = [(r, b1, b2) for (r, b1, b2) in refs if r[0] is not None and r[1] is not None]
    if len(refs) > MAX_ROT:
        refs = refs[:MAX_ROT]   # cap; note: high-flexibility -> global-opt fallback
    grids = [TORSIONS] * len(refs)
    combos = list(itertools.product(*grids)) if refs else [()]

    kept: List[Tuple[float, Chem.Mol]] = []
    for combo in combos:
        m = Chem.Mol(mol)
        conf = m.GetConformer()
        for ((a1, a3), b1, b2), ang in zip(refs, combo):
            rdMolTransforms.SetDihedralDeg(conf, a1, b1, b2, a3, ang)
        if not _PURE_TRACK3:
            if AllChem.MMFFOptimizeMolecule(m) not in (0, 1):
                continue
        e = _conformer_score(m)
        # dedup by RMSD against kept
        dup = False
        for _, mk in kept:
            try:
                if AllChem.GetBestRMS(Chem.RemoveHs(m), Chem.RemoveHs(mk)) < RMSD_DEDUP:
                    dup = True; break
            except Exception:
                pass
        if not dup:
            kept.append((e, m))
    if not kept:
        return len(combos), []
    emin = min(e for e, _ in kept)
    relevant = sorted([(e, m) for e, m in kept if e - emin <= ENERGY_WINDOW],
                      key=lambda x: x[0])
    return len(combos), relevant


if __name__ == "__main__":
    tests = {
        "ethane (no rot)": "CC",
        "butane (1 rot: anti+gauche)": "CCCC",
        "n-pentane": "CCCCC",
        "1,2-dichloroethane": "ClCCCl",
        "ethylene glycol": "OCCO",
    }
    for name, smi in tests.items():
        nc, confs = enumerate_conformers(smi)
        if confs:
            es = [round(e, 2) for e, _ in confs]
            span = round(es[-1] - es[0], 2)
            print(f"{name:<32} combos={nc:>4}  distinct-in-window={len(confs):>2}  "
                  f"ΔE={span} kcal/mol  energies={es[:6]}")
        else:
            print(f"{name:<32} combos={nc:>4}  (embed/opt failed)")
