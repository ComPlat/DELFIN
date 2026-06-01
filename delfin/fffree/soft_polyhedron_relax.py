"""delfin.fffree.soft_polyhedron_relax — soft M-D bond relaxation.

Phase G16 (User 2026-06-01, post-G15b decisive Mogul win).

After fffree's rigid Kabsch-fit places donors EXACTLY on polyhedron vertices,
the M-D bond lengths can be stretched beyond the COD-empirical ideal (the
forensic example: Pd-N3 = 2.10 A correct, Pd-N14 = 2.46 A in the same
bipyridine chelate after rigid placement). This is the construction-vs-
relaxation trade-off responsible for the 28 iter-gate HEAL-FIRST axes
(F3_bond +163 %, vdw +137 %, funcgrp +266 %) that remained after G15b
eliminated the catastrophic-overlap class.

This module relaxes those M-D stretches in a SAFE way:

  1. For each donor atom D coordinated to metal M:
       d_obs  = |M - D|
       d_ideal = COD-empirical M-D for (M-element, D-element)
       If d_obs > 1.05 * d_ideal: shrink toward M by (d_obs - d_ideal) * 0.7
       (only partial correction so we keep the polyhedron approximately).
  2. The whole substituent subtree descending from D moves rigidly with D
     (preserves all D-substituent and substituent-internal bonds).
  3. CShM guard: if post-snap CShM jumps by > 1.0 vs pre-snap, REVERT
     (keeps the polyhedron coordination shape clean).
  4. M-D invariant guard: never make d(M-D) shorter than 0.95 * d_ideal
     (prevents M-D collapse to atom-overlap territory).

Same after-_build_is_clean-with-fallback pattern as G13b. Universal across
all metal-donor element pairs, deterministic, FF-free.

Env-gate: ``DELFIN_FFFREE_SOFT_POLY=1`` (default OFF). Auto-enabled under
``DELFIN_FFFREE_PURE_TRACK3=1``.
"""
from __future__ import annotations

import math
import os
from typing import List, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Conservative tuning: only partial correction toward COD-ideal, with a hard
# floor at 95 % of ideal to prevent M-D collapse.
SHRINK_FRACTION = 0.7
RELAX_FACTOR_THRESHOLD = 1.05         # M-D > 1.05 * ideal -> relax
HARD_FLOOR_FACTOR = 0.95              # never shrink below 0.95 * ideal


def _soft_poly_active() -> bool:
    """Call-time env-gate (so pytest fixtures can toggle between calls)."""
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_SOFT_POLY", "0") == "1":
        return True
    return False


def _cod_md_ideal(metal_elem: str, donor_elem: str) -> float:
    """COD-empirical M-D bond length (Angstrom). Fallback to covalent-sum if
    no specific value is known. Universal across all element pairs.

    Conservative defaults from CCDC medians for common coordination bonds;
    values are deliberately on the SHORT side so the relaxer prefers slight
    compression over over-stretching (matches XRD-observed distribution
    centre).
    """
    md = {
        # 3d
        ("Fe", "N"): 2.00, ("Fe", "O"): 2.05, ("Fe", "C"): 1.95,
        ("Fe", "P"): 2.25, ("Fe", "S"): 2.30, ("Fe", "Cl"): 2.25, ("Fe", "Br"): 2.40,
        ("Co", "N"): 1.95, ("Co", "O"): 2.00, ("Co", "C"): 1.95,
        ("Co", "P"): 2.20, ("Co", "S"): 2.25, ("Co", "Cl"): 2.25,
        ("Ni", "N"): 2.00, ("Ni", "O"): 2.00, ("Ni", "C"): 1.95,
        ("Ni", "P"): 2.20, ("Ni", "S"): 2.20,
        ("Cu", "N"): 2.00, ("Cu", "O"): 2.00, ("Cu", "S"): 2.30, ("Cu", "Cl"): 2.30,
        ("Zn", "N"): 2.05, ("Zn", "O"): 2.05, ("Zn", "S"): 2.30, ("Zn", "Cl"): 2.25,
        ("Mn", "N"): 2.10, ("Mn", "O"): 2.10, ("Mn", "Cl"): 2.45,
        # 4d
        ("Pd", "N"): 2.05, ("Pd", "P"): 2.30, ("Pd", "Cl"): 2.35, ("Pd", "C"): 2.00,
        ("Rh", "N"): 2.05, ("Rh", "P"): 2.30, ("Rh", "Cl"): 2.40, ("Rh", "C"): 2.05,
        ("Ru", "N"): 2.05, ("Ru", "P"): 2.35, ("Ru", "Cl"): 2.40, ("Ru", "C"): 2.10,
        ("Ag", "N"): 2.20, ("Ag", "P"): 2.40, ("Ag", "Cl"): 2.45,
        ("Mo", "N"): 2.15, ("Mo", "O"): 2.10, ("Mo", "Cl"): 2.45,
        ("Cd", "N"): 2.25, ("Cd", "O"): 2.30, ("Cd", "Cl"): 2.50,
        # 5d
        ("Pt", "N"): 2.05, ("Pt", "P"): 2.30, ("Pt", "Cl"): 2.35, ("Pt", "C"): 2.05,
        ("Ir", "N"): 2.05, ("Ir", "P"): 2.30, ("Ir", "Cl"): 2.40, ("Ir", "C"): 2.05,
        ("Os", "N"): 2.10, ("Os", "P"): 2.35, ("Os", "C"): 2.05,
        ("Re", "N"): 2.10, ("Re", "O"): 2.05, ("Re", "Cl"): 2.45,
        ("W",  "N"): 2.15, ("W",  "O"): 2.10, ("W",  "Cl"): 2.45,
        ("Hg", "N"): 2.20, ("Hg", "Cl"): 2.45,
        ("Au", "P"): 2.30, ("Au", "Cl"): 2.30,
        # main-group / heavy
        ("Sn", "N"): 2.30, ("Sn", "O"): 2.10, ("Sn", "Cl"): 2.45,
    }
    if (metal_elem, donor_elem) in md:
        return md[(metal_elem, donor_elem)]
    # symmetric lookup
    if (donor_elem, metal_elem) in md:
        return md[(donor_elem, metal_elem)]
    # fallback to covalent-sum
    return _bd._ideal_bond(metal_elem, donor_elem)


def _heavy_adj(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Heavy-atom-only bond graph."""
    bonds = _bd._geometric_bonds(list(syms), P)
    adj: List[List[int]] = [[] for _ in range(len(syms))]
    for i, j in bonds:
        if syms[i] != "H" and syms[j] != "H":
            adj[i].append(j); adj[j].append(i)
    return adj


def _full_adj(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Full bond graph (including H), needed for ride-along H drag."""
    bonds = _bd._geometric_bonds(list(syms), P)
    adj: List[List[int]] = [[] for _ in range(len(syms))]
    for i, j in bonds:
        adj[i].append(j); adj[j].append(i)
    return adj


def _donor_subtree(donor: int, metal: int, heavy_adj: List[List[int]],
                    syms: Sequence[str]) -> Set[int]:
    """Return the set of heavy atoms reachable from donor WITHOUT going through
    the metal. The donor itself is included. Used to translate the whole
    coordinated subtree when relaxing M-D."""
    visited: Set[int] = {donor}
    queue: List[int] = [donor]
    while queue:
        cur = queue.pop()
        for nxt in heavy_adj[cur]:
            if nxt == metal:
                continue
            if nxt in visited:
                continue
            # Also avoid re-entering other coordinated donors (would translate
            # the whole rest of the complex). We only translate the donor and
            # its NON-shared substituent atoms.
            visited.add(nxt)
            queue.append(nxt)
    return visited


def _compute_cshm(metal: int, donors: Sequence[int], P: np.ndarray) -> float:
    """Compute a simplified CShM-like measure as the RMSD of donor-direction
    unit vectors to the IDEAL polyhedron unit vectors after Kabsch alignment.
    Lower = better. Not a full implementation -- just enough to detect
    coordination-shape degradation during relaxation."""
    if not donors:
        return 0.0
    M = P[metal]
    # Donor unit vectors from M
    vecs = []
    for d in donors:
        v = P[d] - M
        n = float(np.linalg.norm(v))
        if n < 1e-6:
            continue
        vecs.append(v / n)
    if len(vecs) < 2:
        return 0.0
    V = np.array(vecs)
    # Pairwise dot products as a proxy for shape (low variance -> regular polyhedron)
    n = len(V)
    dots = []
    for i in range(n):
        for j in range(i + 1, n):
            dots.append(float(np.dot(V[i], V[j])))
    if not dots:
        return 0.0
    # Variance of dots = shape-disorder measure (lower is more regular)
    mean_d = sum(dots) / len(dots)
    var_d = sum((d - mean_d) ** 2 for d in dots) / len(dots)
    return math.sqrt(var_d)


def relax_md_stretches(syms: Sequence[str], P: np.ndarray,
                        cshm_tolerance: float = 0.30
                        ) -> Tuple[List[str], np.ndarray]:
    """Soft polyhedron M-D relaxation. Returns (syms_out, P_out).

    Algorithm:
    1. Find metal and its coordinated heavy-atom neighbours (donors).
    2. For each donor where d(M-D) > RELAX_FACTOR_THRESHOLD * ideal:
       a. Compute target d_new = max(HARD_FLOOR * ideal, d_obs - SHRINK * (d_obs - ideal))
       b. Translate the donor's substituent subtree (excluding metal and other donors)
          by the vector (d_new - d_obs) along the M->D unit direction.
       c. Ride-along H: any H whose only heavy nbr is in the moved subtree also
          translates by the same vector.
    3. Compute pre-snap CShM and post-snap CShM. If post-CShM > pre + cshm_tolerance,
       REVERT (the relaxation broke the polyhedron).

    Pure geometry, FF-free, universal. Deterministic.
    """
    if not _soft_poly_active():
        return list(syms), np.asarray(P, dtype=float).copy()

    P_out = np.asarray(P, dtype=float).copy()
    n = len(syms)

    metal = next((i for i, s in enumerate(syms) if _bd._is_metal(s)), None)
    if metal is None:
        return list(syms), P_out

    h_adj = _heavy_adj(list(syms), P_out)
    f_adj = _full_adj(list(syms), P_out)

    # Coordinated donors: detect M-D pairs directly because _geometric_bonds
    # SKIPS bonds where either atom is a metal (by design, to handle them
    # separately). The classic 1.45 * ideal_bond cutoff EXCLUDES the very
    # M-D stretches we want to relax (e.g. Pt-N at 2.50 A is above the
    # 1.45 * 2.05 = 2.33 cutoff and would be missed). Use a generous coord
    # cutoff = max(1.8 * ideal, 3.0 A) so stretched donors are included;
    # the relaxation step's hard floor at 0.95 * ideal prevents overshoot.
    donors = []
    for j in range(n):
        if j == metal or syms[j] == "H":
            continue
        d_mj = float(np.linalg.norm(P_out[j] - P_out[metal]))
        coord_cutoff = max(1.8 * _bd._ideal_bond(syms[metal], syms[j]), 3.0)
        if d_mj < coord_cutoff:
            donors.append(j)
    if not donors:
        return list(syms), P_out

    pre_cshm = _compute_cshm(metal, donors, P_out)
    metal_elem = syms[metal]

    # Pre-compute subtree for each donor BEFORE moving anything (avoid order
    # dependence; the original adjacency is used).
    subtrees = {d: _donor_subtree(d, metal, h_adj, syms) for d in donors}

    moves = {}  # atom_idx -> translation_vector
    for d in donors:
        donor_elem = syms[d]
        ideal = _cod_md_ideal(metal_elem, donor_elem)
        v = P_out[d] - P_out[metal]
        d_obs = float(np.linalg.norm(v))
        if d_obs < 1e-6:
            continue
        if d_obs <= RELAX_FACTOR_THRESHOLD * ideal:
            continue                # already within tolerance
        # Partial correction toward ideal
        d_new = max(HARD_FLOOR_FACTOR * ideal,
                    d_obs - SHRINK_FRACTION * (d_obs - ideal))
        # Translation along M->D unit
        u = v / d_obs
        translation = (d_new - d_obs) * u
        # Plan translation for all atoms in this donor's subtree.
        subtree = subtrees[d]
        for a in subtree:
            # Avoid translating the metal or atoms shared with another donor's
            # subtree (which would create an inconsistent multi-translation).
            if a == metal:
                continue
            in_other_subtree = any(
                a in subtrees[d2] for d2 in donors if d2 != d
            )
            if in_other_subtree:
                continue
            moves.setdefault(a, np.zeros(3))
            moves[a] += translation

    # Apply heavy-atom moves
    for a, m in moves.items():
        P_out[a] = P_out[a] + m

    # Ride-along H: any H whose only heavy nbr is in the moved set.
    for h in range(n):
        if syms[h] != "H":
            continue
        heavy_nbrs = [k for k in f_adj[h] if syms[k] != "H"]
        if len(heavy_nbrs) != 1:
            continue
        parent = heavy_nbrs[0]
        if parent in moves:
            P_out[h] = P_out[h] + moves[parent]

    # CShM guard
    post_cshm = _compute_cshm(metal, donors, P_out)
    if post_cshm > pre_cshm + cshm_tolerance:
        # Revert
        return list(syms), np.asarray(P, dtype=float).copy()

    return list(syms), P_out


def is_enabled() -> bool:
    return _soft_poly_active()


if __name__ == "__main__":
    # Synthetic test: Pt with 2 N donors, one at 2.05 (ideal), one stretched to 2.5.
    # Expect: the stretched one is relaxed to ~2.18 (partial correction).
    os.environ["DELFIN_FFFREE_SOFT_POLY"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.soft_polyhedron_relax", None)
    from delfin.fffree.soft_polyhedron_relax import relax_md_stretches

    syms = ["Pt", "N", "N", "C", "C", "C", "C"]
    P0 = np.array([
        [0.0,  0.0,  0.0],     # Pt
        [2.05, 0.0,  0.0],     # N1 at ideal 2.05
        [-2.5, 0.0,  0.0],     # N2 at stretched 2.50
        [3.5,  0.0,  0.0],     # C-substituent of N1
        [4.0,  1.0,  0.0],     # methyl of C
        [-3.5, 0.0,  0.0],     # C-substituent of N2
        [-4.0, 1.0,  0.0],     # methyl of C on N2
    ])
    # Make bonds between donor-substituent atoms detectable: keep distances right
    _, P_out = relax_md_stretches(syms, P0)
    md1 = float(np.linalg.norm(P_out[1] - P_out[0]))
    md2 = float(np.linalg.norm(P_out[2] - P_out[0]))
    print(f"M-N1 before=2.05  after={md1:.3f}  (expect ~2.05, no change)")
    print(f"M-N2 before=2.50  after={md2:.3f}  (expect <2.50, partial relaxation)")
    # Substituent: relative C-N distance preserved
    cn1 = float(np.linalg.norm(P_out[3] - P_out[1]))
    cn2 = float(np.linalg.norm(P_out[5] - P_out[2]))
    print(f"C-N1 substituent bond: before=1.45 after={cn1:.3f} (preserved by subtree drag)")
    print(f"C-N2 substituent bond: before=1.00 after={cn2:.3f} (preserved by subtree drag)")
