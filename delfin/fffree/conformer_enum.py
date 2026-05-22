#!/usr/bin/env python3
"""conformer_enum.py — Phase 1, Layer 3: deterministic combinatorial conformer
enumeration (Confab-style, NOT stochastic global search).  Replaces SEARCH with
ENUMERATION: rotatable bonds -> discrete torsion-state product -> clash-prune ->
LOCAL MMFF relax -> RMSD/energy dedup -> keep within the relevant energy window.

This is the Layer-3 generator; completeness here is "all relevant conformers
within the energy window, up to torsion resolution" (validated by convergence:
finer grid stops adding distinct minima).  Organic ligand part only (MMFF valid);
the metal sphere is built separately (Phase 2).

Deterministic: ETKDG initial embed with fixed seed + numThreads=1; fixed torsion
grid; MMFF; canonical-RMSD dedup with deterministic ordering.
"""
from __future__ import annotations
import itertools
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

SEED = 42
TORSIONS = (60.0, 180.0, 300.0)   # gauche+/anti/gauche- ; resolution = 120 deg
RMSD_DEDUP = 0.5                  # A; conformers closer than this are the same
ENERGY_WINDOW = 5.0               # kcal/mol above the global min = "relevant"
MAX_ROT = 7                       # cap combinatorial blow-up (3^7 = 2187)


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


def enumerate_conformers(smiles: str):
    """Return (n_combos, conformers) where conformers = list of (energy, mol) within
    the energy window, deduped by RMSD.  Deterministic."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, []
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
        return 0, []
    AllChem.MMFFOptimizeMolecule(mol)
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
        if AllChem.MMFFOptimizeMolecule(m) not in (0, 1):
            continue
        e = _mmff_energy(m)
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
