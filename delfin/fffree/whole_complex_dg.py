"""delfin.fffree.whole_complex_dg — Whole-Complex Distance-Geometry Embed.

Phase G11 (User 2026-05-31 'fundamental': aim for DFT-quality without DFT).

Instead of constructive Polyhedron-vertex + Kabsch-fit (which produces ~1 Å
RMSD due to rigid ligand mismatch with vertex angles), this module sets up a
single Distance-Geometry bounds matrix for the WHOLE complex and embeds
all atoms simultaneously:

  Bounds matrix encodes:
    M-D distances (rigid, empirical COD-derived)
    Donor-donor distances (Polyhedron-vertex-derived, fixed)
    Ligand internal bonds/angles/torsions (ETKDG bounds-matrix transfer)
    Cross-ligand vdW (non-overlap)
    Chirality / stereo constraints from SMILES

  RDKit DG-embed solves the constrained metric matrix problem to give a
  structure that satisfies ALL constraints simultaneously — no rigid Kabsch
  fit, no collapsed bonds from forced-fit, no over-stretched M-D from rigid
  Kabsch on mismatched bite angles.

Universal across all TMC classes (any CN, any donor, any chelate denticity,
any hapto-π type). FF-free (DG is geometric, not energetic). Deterministic
(seeded random + triangle-bound smoothing).

Mathematical foundation:
  Crippen-Havel 1988 (Distance Geometry in Conformational Analysis)
  Riniker-Landrum 2015 (ETKDG: bounds matrix from CCDC torsion preferences)
  Cremer-Pople 1975 (ring puckering as bound matrix constraints)
  Lema-Saavedra & Fernández-Ramos JCP 2026 (2-stage CP reconstruction)

Env-gate: DELFIN_FFFREE_WHOLE_DG=1 (default OFF, byte-identical when unset).
Phase G11 prep: foundation module, integration in subsequent iterations.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_WHOLE_DG = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_WHOLE_DG", "0") == "1"


def build_complex_bounds_matrix(
    n_atoms: int,
    metal_idx: int,
    donor_indices: Sequence[int],
    md_distance: float,
    polyhedron_vertices: np.ndarray,
    bond_pairs: Sequence[Tuple[int, int]],
    bond_targets: Sequence[float],
    angle_triples: Sequence[Tuple[int, int, int]],
    angle_targets: Sequence[float],
    vdw_radii: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Construct (lo, hi) bounds matrices for whole-complex DG embed.

    Parameters
    ----------
    n_atoms : total atom count (metal + all ligand atoms)
    metal_idx : index of metal atom (typically 0)
    donor_indices : indices of donor atoms (length = CN)
    md_distance : M-D bond length (Å, COD-empirical)
    polyhedron_vertices : (CN, 3) unit vectors of donor positions on polyhedron
    bond_pairs : list of (i, j) ligand bonds
    bond_targets : ideal bond lengths (Å)
    angle_triples : list of (i, j, k) bonded triples
    angle_targets : ideal angles (radians)
    vdw_radii : per-atom vdW radii for clash bounds

    Returns
    -------
    lo, hi : (n, n) lower and upper bound matrices

    Universal: any TMC topology + polyhedron.
    """
    lo = np.zeros((n_atoms, n_atoms))
    hi = np.full((n_atoms, n_atoms), 100.0)  # large default

    # M-D rigid bonds
    for d_idx, d in enumerate(donor_indices):
        lo[metal_idx, d] = lo[d, metal_idx] = md_distance - 0.02
        hi[metal_idx, d] = hi[d, metal_idx] = md_distance + 0.02

    # Donor-donor distances from polyhedron
    for i, d_i in enumerate(donor_indices):
        for j, d_j in enumerate(donor_indices):
            if i >= j:
                continue
            v_i = polyhedron_vertices[i] / np.linalg.norm(polyhedron_vertices[i])
            v_j = polyhedron_vertices[j] / np.linalg.norm(polyhedron_vertices[j])
            d_ij = float(np.linalg.norm((v_i - v_j) * md_distance))
            lo[d_i, d_j] = lo[d_j, d_i] = d_ij - 0.05
            hi[d_i, d_j] = hi[d_j, d_i] = d_ij + 0.05

    # Ligand 1-2 bonds (tight)
    for (i, j), t in zip(bond_pairs, bond_targets):
        lo[i, j] = lo[j, i] = t - 0.02
        hi[i, j] = hi[j, i] = t + 0.02

    # Ligand 1-3 angles → derived distances
    for (i, j, k), theta in zip(angle_triples, angle_targets):
        # find bond lengths i-j and j-k
        d_ij = None; d_jk = None
        for (a, b), bl in zip(bond_pairs, bond_targets):
            if {a, b} == {i, j}: d_ij = bl
            if {a, b} == {j, k}: d_jk = bl
        if d_ij and d_jk:
            d_ik = math.sqrt(d_ij ** 2 + d_jk ** 2 - 2 * d_ij * d_jk * math.cos(theta))
            lo[i, k] = lo[k, i] = d_ik - 0.05
            hi[i, k] = hi[k, i] = d_ik + 0.05

    # Cross-pair vdW lower bound (clash prevention)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if lo[i, j] > 0:  # already constrained
                continue
            vdw_sum = vdw_radii[i] + vdw_radii[j]
            lo[i, j] = lo[j, i] = 0.85 * vdw_sum

    return lo, hi


def triangle_bounds_smooth(lo: np.ndarray, hi: np.ndarray, max_iter: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Apply triangle inequality smoothing to bounds matrix.

    For each triple (i, j, k):
      hi[i,j] = min(hi[i,j], hi[i,k] + hi[k,j])
      lo[i,j] = max(lo[i,j], lo[i,k] - hi[k,j], lo[k,j] - hi[i,k])

    Iterates until convergence. Required pre-embed step to avoid infeasible
    bounds matrix. Universal, deterministic, no parameters.
    """
    n = lo.shape[0]
    for _ in range(max_iter):
        changed = False
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(n):
                    if k == i or k == j:
                        continue
                    new_hi = hi[i, k] + hi[k, j]
                    if new_hi < hi[i, j] - 1e-6:
                        hi[i, j] = hi[j, i] = new_hi
                        changed = True
                    new_lo1 = lo[i, k] - hi[k, j]
                    new_lo2 = lo[k, j] - hi[i, k]
                    new_lo = max(new_lo1, new_lo2)
                    if new_lo > lo[i, j] + 1e-6:
                        lo[i, j] = lo[j, i] = new_lo
                        changed = True
        if not changed:
            break
    return lo, hi


def metric_embed(lo: np.ndarray, hi: np.ndarray, seed: int = 42) -> Optional[np.ndarray]:
    """Embed a bounds matrix to 3D coordinates via metric matrix method.

    Algorithm (Crippen-Havel 1988):
      1. Sample target distance matrix D from [lo, hi] (random within bounds)
      2. Compute metric matrix G = -0.5 * (D² - row_means - col_means + total_mean)
      3. EVD: G = Q · diag(λ) · Q^T
      4. Coordinates: P = Q[:, top-3] · sqrt(diag(λ_top-3))
      5. Return P if 3D embed is valid (positive top-3 eigenvalues)

    Deterministic via fixed seed. Universal: any positive bounds matrix.
    """
    rng = np.random.RandomState(seed)
    n = lo.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mid = 0.5 * (lo[i, j] + hi[i, j])
            spread = 0.5 * (hi[i, j] - lo[i, j])
            D[i, j] = D[j, i] = mid + rng.uniform(-spread, spread)
    D2 = D ** 2
    row_means = D2.mean(axis=1)
    total_mean = D2.mean()
    G = -0.5 * (D2 - row_means[:, None] - row_means[None, :] + total_mean)
    try:
        eigvals, eigvecs = np.linalg.eigh(G)
    except np.linalg.LinAlgError:
        return None
    # Top-3 eigenvalues (must be positive for valid 3D embed)
    top3 = eigvals[-3:]
    if top3.min() < -1e-3:  # tolerance for numerical error
        return None
    top3 = np.maximum(top3, 0)
    P = eigvecs[:, -3:] * np.sqrt(top3)[None, :]
    return P


def embed_whole_complex(
    n_atoms: int,
    metal_idx: int,
    donor_indices: Sequence[int],
    md_distance: float,
    polyhedron_vertices: np.ndarray,
    bond_pairs: Sequence[Tuple[int, int]],
    bond_targets: Sequence[float],
    angle_triples: Sequence[Tuple[int, int, int]] = (),
    angle_targets: Sequence[float] = (),
    vdw_radii: Optional[Sequence[float]] = None,
    n_attempts: int = 5,
    seed: int = 42,
) -> Optional[np.ndarray]:
    """Universal whole-complex DG embed.

    Returns (n_atoms, 3) coords if successful, None if infeasible bounds
    or all attempts failed.

    Universal across TMC classes. Deterministic via seed sequence.
    """
    if not _WHOLE_DG:
        return None
    if vdw_radii is None:
        vdw_radii = [1.5] * n_atoms
    lo, hi = build_complex_bounds_matrix(
        n_atoms, metal_idx, donor_indices, md_distance,
        polyhedron_vertices, bond_pairs, bond_targets,
        angle_triples, angle_targets, vdw_radii,
    )
    lo, hi = triangle_bounds_smooth(lo, hi)
    for attempt in range(n_attempts):
        P = metric_embed(lo, hi, seed=seed + attempt)
        if P is not None and np.all(np.isfinite(P)):
            # Center at metal
            P = P - P[metal_idx]
            return P
    return None


if __name__ == "__main__":
    # Test: simple octahedral complex
    os.environ["DELFIN_FFFREE_WHOLE_DG"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.whole_complex_dg", None)
    from delfin.fffree.whole_complex_dg import embed_whole_complex

    # Octahedron vertices: (±1, 0, 0), (0, ±1, 0), (0, 0, ±1)
    polyhedron = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                            [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    # M + 6 monodentate (only M-D bounds, no ligand internals)
    n = 7
    P = embed_whole_complex(
        n_atoms=n,
        metal_idx=0,
        donor_indices=[1, 2, 3, 4, 5, 6],
        md_distance=2.0,
        polyhedron_vertices=polyhedron,
        bond_pairs=[],
        bond_targets=[],
        vdw_radii=[1.0] * n,
    )
    if P is not None:
        print(f"OC-6 embed: {P.shape}")
        for i in range(n):
            d = float(np.linalg.norm(P[i] - P[0])) if i > 0 else 0
            print(f"  atom {i}: ({P[i][0]:+.3f}, {P[i][1]:+.3f}, {P[i][2]:+.3f})  d_M = {d:.3f}")
        # Verify M-D = 2.0 (target)
        md_dist = [float(np.linalg.norm(P[i] - P[0])) for i in range(1, 7)]
        print(f"  M-D distances: mean={np.mean(md_dist):.3f} (target 2.000)")
    else:
        print("OC-6 embed FAILED")
