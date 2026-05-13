"""Symmetry detection algorithms for Baustein 6.

Implements the four detection tiers:

  * Tier A - Hungarian assignment of donors to ideal polyhedron slots
  * Tier B - Chemical equivalence via Morgan canonical ranking
             (atoms, bond pairs, angle triples)
  * Tier C - Fragment / SMARTS matching (delegated to RDKit at call site;
             this module only ships the equivalence-class helpers Tier C
             needs).
  * Tier D - Global molecular point-group detection from 3D coordinates

All RDKit / scipy / sister-module imports are *lazy* inside the public
functions so that this module can be imported in environments where the
heavy dependencies are not available (callers see an empty fallback
result instead of an ImportError).

Tolerances default to 0.3 Å which is chemically reasonable for slightly
distorted but recognisable symmetric structures.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Tier B - chemical equivalence classes
# ---------------------------------------------------------------------------

def _canonical_ranks(mol) -> Optional[List[int]]:
    """Return Morgan canonical ranks (breakTies=False) or None on failure."""
    try:
        from rdkit.Chem import CanonicalRankAtoms
    except ImportError:
        return None
    try:
        return list(CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        return None


def find_equivalent_atoms(mol) -> List[Set[int]]:
    """Tier B: chemically equivalent atom groups via Morgan canonical ranking.

    Returns
    -------
    list of sets
        Each set contains atom indices that share a canonical rank.
        Singletons (rank seen exactly once) are *excluded* - only true
        multi-atom equivalence classes are returned.
    """
    ranks = _canonical_ranks(mol)
    if ranks is None:
        return []

    groups: Dict[int, Set[int]] = defaultdict(set)
    for atom_idx, rank in enumerate(ranks):
        groups[rank].add(atom_idx)
    return [g for g in groups.values() if len(g) > 1]


def find_equivalent_bond_pairs(mol) -> List[List[Tuple[int, int]]]:
    """Tier B: bonds grouped by chemical equivalence.

    Two bonds (i, j) and (k, l) belong to the same class iff the
    *sorted* pair of canonical ranks of their endpoints matches.

    Returns
    -------
    list of lists of (i, j) tuples
        Only multi-member classes (>= 2 bonds) are returned.
    """
    ranks = _canonical_ranks(mol)
    if ranks is None:
        return []

    bond_groups: Dict[Tuple[int, int], List[Tuple[int, int]]] = defaultdict(list)
    try:
        bonds = list(mol.GetBonds())
    except Exception:
        return []

    for bond in bonds:
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        key = tuple(sorted([ranks[i], ranks[j]]))
        bond_groups[key].append((i, j))

    return [pairs for pairs in bond_groups.values() if len(pairs) > 1]


def find_equivalent_angle_triples(mol) -> List[List[Tuple[int, int, int]]]:
    """Tier B: angle triples (i, c, j) grouped by chemical equivalence.

    For every heavy atom c with >= 2 neighbours we generate all
    unique (i, c, j) triples (i < j) and group them by the tuple
    (rank[c], sorted(rank[i], rank[j])).

    Returns
    -------
    list of lists of (i, c, j) tuples
        Only multi-member classes (>= 2 triples) are returned.
    """
    ranks = _canonical_ranks(mol)
    if ranks is None:
        return []

    triple_groups: Dict[Tuple[int, Tuple[int, int]],
                        List[Tuple[int, int, int]]] = defaultdict(list)
    try:
        atoms = list(mol.GetAtoms())
    except Exception:
        return []

    for c_atom in atoms:
        c = c_atom.GetIdx()
        nbrs = [n.GetIdx() for n in c_atom.GetNeighbors()]
        if len(nbrs) < 2:
            continue
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                i, j = nbrs[a], nbrs[b]
                pair_key = tuple(sorted([ranks[i], ranks[j]]))
                key = (ranks[c], pair_key)
                triple_groups[key].append((i, c, j))

    return [triples for triples in triple_groups.values() if len(triples) > 1]


# ---------------------------------------------------------------------------
# Tier A - Hungarian donor -> ideal-slot assignment
# ---------------------------------------------------------------------------

def hungarian_donor_assignment_cost_matrix(
    current_units: np.ndarray,
    ideal_units: np.ndarray,
    donor_types: List[str],
    slot_types_preferred: Optional[List[str]] = None,
) -> np.ndarray:
    """Build cost matrix for Hungarian donor-to-slot assignment.

    cost[r, c] = - cos_sim(current_units[r], ideal_units[c])
                 + type_mismatch_penalty(donor_types[r], slot_types_preferred[c])

    Lower cost = better match. Cosine similarity is in [-1, 1] so the
    geometric term is in [-1, 1]; the mismatch penalty (default 0.0,
    0.25 if a preferred slot type is given and disagrees with the donor
    element symbol) is small enough not to override geometry but enough
    to break ties.
    """
    cur = np.asarray(current_units, dtype=float)
    ide = np.asarray(ideal_units, dtype=float)
    if cur.ndim != 2 or ide.ndim != 2 or cur.shape[1] != 3 or ide.shape[1] != 3:
        raise ValueError("current_units / ideal_units must be (n, 3) arrays")

    # Cosine similarity (rows already unit-norm by contract, but stay safe).
    cur_norm = np.linalg.norm(cur, axis=1, keepdims=True)
    ide_norm = np.linalg.norm(ide, axis=1, keepdims=True)
    cur_safe = cur / np.clip(cur_norm, 1e-9, None)
    ide_safe = ide / np.clip(ide_norm, 1e-9, None)
    cos_sim = cur_safe @ ide_safe.T
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    cost = -cos_sim

    if slot_types_preferred is not None and len(slot_types_preferred) == ide.shape[0]:
        penalty = 0.25
        for r, dtype in enumerate(donor_types):
            for c, stype in enumerate(slot_types_preferred):
                if stype and dtype and stype != dtype:
                    cost[r, c] += penalty
    return cost


def hungarian_assign_donors_to_slots(
    coords: np.ndarray,
    mol,
    metal_idx: int,
) -> Dict[int, np.ndarray]:
    """Tier A: Hungarian assignment of metal's donors to ideal polyhedron slots.

    Steps:
      1. Identify donors (heavy neighbours of metal that are themselves
         not metals).
      2. CN = number of donors; classify polyhedron from (CN, donor types).
      3. Get ideal donor unit vectors from :mod:`delfin._polyhedron_targets`.
      4. Compute current donor unit vectors (M -> donor direction).
      5. Build cost matrix (negative cosine similarity).
      6. ``scipy.optimize.linear_sum_assignment`` -> donor -> slot mapping.
      7. Target = metal_pos + ideal_unit * d_ideal(M, donor).

    Returns
    -------
    dict
        ``{donor_atom_idx: target_3d_position (np.ndarray shape (3,))}``.
        Empty dict on any failure (missing sister module, no donors,
        degenerate geometry, ...).
    """
    # Lazy imports - graceful fallback if anything is missing.
    try:
        from delfin._polyhedron_targets import (
            classify_geometry_from_cn_donors,
            get_ideal_donor_vectors,
        )
    except ImportError:
        return {}

    try:
        from delfin.smiles_converter import _METAL_SET, _get_ml_bond_length
    except ImportError:
        return {}

    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError:
        return {}

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        return {}

    try:
        m_atom = mol.GetAtomWithIdx(int(metal_idx))
    except Exception:
        return {}

    donors_neighbors = [
        n for n in m_atom.GetNeighbors()
        if n.GetSymbol() not in _METAL_SET
    ]
    if not donors_neighbors:
        return {}

    donor_indices = [n.GetIdx() for n in donors_neighbors]
    donor_types = [mol.GetAtomWithIdx(d).GetSymbol() for d in donor_indices]
    cn = len(donor_indices)

    # Classify polyhedron geometry from CN + donor element types.
    try:
        geom_result = classify_geometry_from_cn_donors(cn, donor_types)
    except Exception:
        return {}
    if isinstance(geom_result, tuple):
        geometry = geom_result[0]
    else:
        geometry = geom_result

    try:
        ideal_vectors = get_ideal_donor_vectors(cn, geometry)
    except Exception:
        return {}
    ideal_vectors = np.asarray(ideal_vectors, dtype=float)
    if ideal_vectors.ndim != 2 or ideal_vectors.shape != (cn, 3):
        return {}

    # Current donor unit vectors (M -> donor).
    metal_pos = coords[int(metal_idx)]
    current_units = np.zeros((cn, 3), dtype=float)
    for r, d_idx in enumerate(donor_indices):
        diff = coords[d_idx] - metal_pos
        norm = float(np.linalg.norm(diff))
        if norm < 1e-6:
            # Degenerate (donor sits on metal) - point along +x as filler.
            current_units[r] = np.array([1.0, 0.0, 0.0])
        else:
            current_units[r] = diff / norm

    # Cost matrix + Hungarian.
    try:
        cost = hungarian_donor_assignment_cost_matrix(
            current_units, ideal_vectors, donor_types,
        )
        row_ind, col_ind = linear_sum_assignment(cost)
    except Exception:
        return {}

    # Build target positions.
    metal_sym = m_atom.GetSymbol()
    targets: Dict[int, np.ndarray] = {}
    for r, c in zip(row_ind, col_ind):
        donor_idx = donor_indices[int(r)]
        donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
        try:
            d_ideal = float(_get_ml_bond_length(metal_sym, donor_sym))
        except Exception:
            d_ideal = 2.0  # safe fallback (~typical M-L bond)
        targets[donor_idx] = metal_pos + ideal_vectors[int(c)] * d_ideal

    return targets


# ---------------------------------------------------------------------------
# Tier D - Global point-group detection
# ---------------------------------------------------------------------------

# Candidates in approximate order of *symmetry order* (highest first).
# We stop at the first group whose ALL operations validate against the
# 3D structure within tolerance.
_POINT_GROUP_CANDIDATES: Tuple[str, ...] = (
    "Oh", "Td",
    "D6h", "D5h", "D4h", "D3h", "D2h",
    "D3", "D2",
    "C6v", "C5v", "C4v", "C3v", "C2v",
    "C3", "C2",
    "S6", "S4",
    "Cs", "Ci",
)


def _get_point_group_operations(pg_name: str) -> Optional[List[np.ndarray]]:
    """Fetch operations for a point group (lazy import, None on failure)."""
    try:
        from delfin import _point_group_ops as pg_mod  # type: ignore
    except ImportError:
        return None
    getter = getattr(pg_mod, "get_operations", None)
    if getter is None:
        return None
    try:
        ops = getter(pg_name)
    except (KeyError, ValueError):
        return None
    except Exception:
        return None
    if ops is None:
        return None
    out: List[np.ndarray] = []
    for op in ops:
        arr = np.asarray(op, dtype=float)
        if arr.shape != (3, 3):
            return None
        out.append(arr)
    return out


def detect_global_point_group(
    mol,
    coords: np.ndarray,
    tolerance: float = 0.3,
) -> Tuple[str, List[np.ndarray], Dict[int, Dict[int, int]]]:
    """Tier D: detect actual molecular point group from 3D structure.

    Algorithm
    ---------
    1. Quick check: if every Morgan canonical rank is unique, there can
       be no graph automorphism and therefore no point-group symmetry
       above ``C1`` -> return early.
    2. Iterate candidate groups in descending symmetry order. For each:
         a. Fetch operations.
         b. For every operation, find a self-consistent permutation
            of atoms that maps within ``tolerance``. Atom partners must
            share the same canonical rank.
         c. If every operation succeeds -> accept this group.
    3. Fallback ``C1``.

    Returns
    -------
    (group_name, operations, atom_perms)
        ``operations`` is a list of 3x3 numpy arrays (identity first).
        ``atom_perms`` is ``{op_index: {atom_idx: partner_atom_idx}}``.
        For ``C1`` both lists/dicts are empty.

    Notes
    -----
    * Tolerance default 0.3 A is chemically reasonable for slightly
      distorted but recognisable symmetric structures.
    * Complexity per candidate: O(n_atoms^2 * n_ops). For molecules with
      >150 atoms we short-circuit to C1 (Tier D is only meaningful for
      small, near-symmetric cores anyway).
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        return "C1", [], {}

    n_atoms = coords.shape[0]
    if n_atoms == 0:
        return "C1", [], {}

    # Hard cap on size - Tier D is for cores, not full huge ligands.
    if n_atoms > 150:
        return "C1", [], {}

    ranks = _canonical_ranks(mol)
    if ranks is None or len(ranks) != n_atoms:
        return "C1", [], {}

    # Step 1: all-unique ranks => no possible non-trivial automorphism.
    if len(set(ranks)) == len(ranks):
        return "C1", [], {}

    # Center coordinates at centroid (point groups assume origin = COM-ish).
    centroid = np.mean(coords, axis=0)
    coords_c = coords - centroid

    # Pre-index atoms by canonical rank for cheap partner lookup.
    rank_to_atoms: Dict[int, List[int]] = defaultdict(list)
    for i, r in enumerate(ranks):
        rank_to_atoms[r].append(i)

    tol2 = float(tolerance) * float(tolerance)

    for pg_name in _POINT_GROUP_CANDIDATES:
        ops = _get_point_group_operations(pg_name)
        if ops is None or len(ops) == 0:
            continue

        all_valid = True
        atom_perms: Dict[int, Dict[int, int]] = {}

        for op_idx, op_matrix in enumerate(ops):
            # Identity passes trivially (and avoids degenerate mapping work).
            if np.allclose(op_matrix, np.eye(3), atol=1e-8):
                atom_perms[op_idx] = {i: i for i in range(n_atoms)}
                continue

            mapping: Dict[int, int] = {}
            used_partners: Set[int] = set()
            op_valid = True

            for i in range(n_atoms):
                transformed = op_matrix @ coords_c[i]
                candidates = rank_to_atoms.get(ranks[i], ())
                best_j = -1
                best_d2 = float("inf")
                for j in candidates:
                    if j in used_partners:
                        continue
                    delta = transformed - coords_c[j]
                    d2 = float(delta @ delta)
                    if d2 < best_d2:
                        best_d2 = d2
                        best_j = j
                if best_j < 0 or best_d2 > tol2:
                    op_valid = False
                    break
                mapping[i] = best_j
                used_partners.add(best_j)

            if not op_valid:
                all_valid = False
                break
            atom_perms[op_idx] = mapping

        if all_valid:
            return pg_name, ops, atom_perms

    return "C1", [], {}


# ---------------------------------------------------------------------------
# Self-test (small synthetic molecules)
# ---------------------------------------------------------------------------

def _selftest() -> None:  # pragma: no cover - manual smoke-test
    try:
        from rdkit import Chem
    except ImportError:
        print("[selftest] RDKit not available - skipped")
        return

    cases = {
        "H2O":     "O",
        "NH3":     "N",
        "CH4":     "C",
        "Benzene": "c1ccccc1",
        "CH3Cl":   "CCl",   # C3v-like methyl chloride (no stereo)
    }

    # Geometries are placeholders - Tier D will probably classify as C1
    # unless the embed produces near-perfect symmetry, but find_equivalent_*
    # functions work purely on connectivity.
    for name, smi in cases.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"[selftest] {name}: cannot parse SMILES")
            continue
        m = Chem.AddHs(m)
        try:
            from rdkit.Chem import AllChem
            res = AllChem.EmbedMolecule(m, randomSeed=42)
            if res != 0:
                print(f"[selftest] {name}: embed failed (rc={res})")
                continue
            AllChem.UFFOptimizeMolecule(m)
        except Exception as exc:
            print(f"[selftest] {name}: embed exception {exc}")
            continue

        conf = m.GetConformer()
        coords = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())],
            dtype=float,
        )

        eq_atoms = find_equivalent_atoms(m)
        eq_bonds = find_equivalent_bond_pairs(m)
        eq_angles = find_equivalent_angle_triples(m)
        pg, ops, perms = detect_global_point_group(m, coords, tolerance=0.3)
        print(
            f"[selftest] {name:8s} n_atoms={m.GetNumAtoms():3d} "
            f"eq_atom_classes={len(eq_atoms):2d} "
            f"eq_bond_classes={len(eq_bonds):2d} "
            f"eq_angle_classes={len(eq_angles):2d} "
            f"point_group={pg} (ops={len(ops)})"
        )


if __name__ == "__main__":  # pragma: no cover
    _selftest()
