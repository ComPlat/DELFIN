#!/usr/bin/env python3
"""coordination_mode_enum.py — enumerate the valid COORDINATION
MODES of a ligand (the hard open core).  For a ligand with multiple potential
coordination sites, enumerate: donor-selection + denticity (κ1/κ2/κ3…) + linkage
isomerism (ambidentate N-vs-O, S-vs-N, carboxylate κ1-vs-κ2).

Universal, graph-only (no SMILES-specific rules): donors identified by element +
lone-pair availability; chelate feasibility by metallacycle ring size from the
molecular graph; symmetric donors merged by RDKit canonical rank so equivalent
modes are not double-counted.  Deterministic (canonical SMILES atom order).

VALIDITY RULES (the calibratable core; refine against COD/CCDC):
  - donor eligible: N/O/P/S/halide with formal charge <= 0 and degree < 4 (lone pair)
  - chelate ring size 4..7 feasible (4 = carboxylate-type strained-but-real, 5-6
    ideal, 7 ok); 3 = not a chelate (side-on), >7 = unusual -> excluded
  - max denticity capped at MAX_KAPPA
  - (TODO refine: exclude linear-backbone chelates e.g. SCN κ2; charge/CN balance)

This is the Layer-1 enumerator; the count it returns is the Layer-1 coverage
denominator (modes the generator must build).  Tested on en/acetate/glycinate/
SCN/pyridine/terpy.
"""
from __future__ import annotations
import itertools
from typing import Dict, List, Tuple
from rdkit import Chem

DONOR_ELEMENTS = {"N", "O", "P", "S", "F", "Cl", "Br", "I", "Se", "As"}
CHELATE_MIN, CHELATE_MAX = 4, 7      # feasible metallacycle ring size (atoms incl. metal)
MAX_KAPPA = 6


def _donor_atoms(mol) -> List[int]:
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in DONOR_ELEMENTS and a.GetFormalCharge() <= 0 and a.GetDegree() < 4:
            out.append(a.GetIdx())
    return out


def _ring_size(mol, i: int, j: int) -> int:
    """Metallacycle ring size if donors i,j chelate = (atoms on shortest path) + 1
    (the metal).  Returns 0 if no path."""
    path = Chem.GetShortestPath(mol, i, j)
    return len(path) + 1 if path else 0


def _feasible_pair(mol, i: int, j: int) -> bool:
    return CHELATE_MIN <= _ring_size(mol, i, j) <= CHELATE_MAX


def _canon_ranks(mol) -> List[int]:
    return list(Chem.CanonicalRankAtoms(mol, breakTies=False))


def enumerate_modes(smiles: str) -> List[Dict]:
    """Return the distinct coordination modes of the ligand.  Each mode:
    {donors: tuple[int], kappa: int, ring_sizes: tuple[int], donor_elems: tuple[str]}.
    Symmetric modes merged via canonical ranks."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    donors = _donor_atoms(mol)
    if not donors:
        return []
    rank = _canon_ranks(mol)
    sym = {a.GetIdx(): a.GetSymbol() for a in mol.GetAtoms()}
    seen = set()
    modes: List[Dict] = []

    def add(subset: Tuple[int, ...], ring_sizes: Tuple[int, ...]):
        # canonical signature: multiset of (rank) for symmetry-merge
        sig = tuple(sorted(rank[i] for i in subset))
        if sig in seen:
            return
        seen.add(sig)
        modes.append({
            "donors": subset,
            "kappa": len(subset),
            "ring_sizes": ring_sizes,
            "donor_elems": tuple(sym[i] for i in subset),
        })

    # κ1: each donor (symmetry-merged) — captures linkage isomers (distinct ranks)
    for d in donors:
        add((d,), ())

    # κ2..MAX_KAPPA: subsets where the chelate graph is connected via feasible pairs
    for k in range(2, min(len(donors), MAX_KAPPA) + 1):
        for subset in itertools.combinations(donors, k):
            # feasibility: the subset must form a connected chelate using only
            # feasible donor-donor pairs (spanning connectivity)
            feas_edges = [(a, b) for a, b in itertools.combinations(subset, 2)
                          if _feasible_pair(mol, a, b)]
            if not feas_edges:
                continue
            # connected? union-find over subset
            parent = {d: d for d in subset}
            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]; x = parent[x]
                return x
            for a, b in feas_edges:
                parent[find(a)] = find(b)
            if len({find(d) for d in subset}) != 1:
                continue
            ring_sizes = tuple(sorted(_ring_size(mol, a, b) for a, b in feas_edges))
            add(tuple(sorted(subset)), ring_sizes)
    return modes


if __name__ == "__main__":
    tests = {
        "ethylenediamine (NCCN)": "NCCN",
        "acetate (CC(=O)[O-])": "CC(=O)[O-]",
        "glycinate": "[NH2]CC(=O)[O-]",
        "thiocyanate (ambidentate)": "[S-]C#N",
        "pyridine (monodentate)": "c1ccncc1",
        "2,2'-bipyridine": "c1ccc(-c2ccccn2)nc1",
        "terpyridine": "c1ccc(-c2cccc(-c3ccccn3)n2)nc1",
        "water": "O",
    }
    for name, smi in tests.items():
        ms = enumerate_modes(smi)
        print(f"\n{name}: {len(ms)} modes")
        for m in ms:
            print(f"   κ{m['kappa']} {''.join(m['donor_elems'])}  rings={m['ring_sizes']}")
