"""donor_slide.py — Thomson-problem-based donor relaxation on M-D sphere.

Fundamental physics-principle implementation for poly_mean_dev_max reduction
WITHOUT templates, WITHOUT force-field, preserving M-D distance exactly.

Approach (Thomson 1904):
  N charges constrained to sphere repel each other → equilibrium = max-distance
  configuration = naturally polyhedron-like (3=triangle, 4=tetrahedron, 6=octahedron, etc.)

For each donor:
  - Compute pairwise repulsion force from other donors (1/r^2 falloff)
  - Project force to TANGENT plane of sphere (remove radial component)
  - Step in tangent direction
  - Re-project to sphere |M-D|=const (drift correction)

Plus optional chelate-bond constraint: chelate donor pairs have additional
spring-bond keeping their natural bite preserved.

Universal: any TMC with any CN
Deterministic: no RNG, fixed iteration
M-D invariant: preserved exactly by sphere constraint
No templates: physics principle (Thomson + Coulomb repulsion)
No force field: simple analytical 1/r^2, this is classical mechanics

Env-gate: DELFIN_FFFREE_DONOR_SLIDE=1 (default OFF, byte-identical when unset)
Configurable: DELFIN_FFFREE_DONOR_SLIDE_ITER (default 50), _STEP (default 0.05)
"""
from __future__ import annotations

import math
import os
from typing import List, Optional, Sequence, Set, Tuple

import numpy as np


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_DONOR_SLIDE", "0") == "1"


def _config() -> Tuple[int, float, float]:
    """Returns (max_iter, step_size, convergence_tol)."""
    try:
        n_iter = int(os.environ.get("DELFIN_FFFREE_DONOR_SLIDE_ITER", "50"))
    except Exception:
        n_iter = 50
    try:
        step = float(os.environ.get("DELFIN_FFFREE_DONOR_SLIDE_STEP", "0.05"))
    except Exception:
        step = 0.05
    try:
        tol = float(os.environ.get("DELFIN_FFFREE_DONOR_SLIDE_TOL", "1e-4"))
    except Exception:
        tol = 1e-4
    return n_iter, step, tol


def constrained_donor_slide(
    coords: np.ndarray,
    metal_idx: int,
    donor_indices: Sequence[int],
    chelate_pairs: Optional[Sequence[Tuple[int, int]]] = None,
) -> np.ndarray:
    """Sphere-constrained donor relaxation via Thomson-problem dynamics.

    Each donor sits on its M-D sphere |M-D|=const and feels Coulomb-like
    repulsion from other donors. Tangential gradient -> equilibrium ~ ideal
    polyhedron (octahedron/tetrahedron/etc) regardless of starting position.

    Returns updated coords (donor positions slid; non-donors unchanged).

    Parameters
    ----------
    coords : (n, 3) array - atomic positions
    metal_idx : int - metal atom index (must NOT be in donor_indices)
    donor_indices : sequence of int - donor atom indices
    chelate_pairs : optional sequence of (i, j) - chelate donor pairs whose
        natural bite distance should be preserved (additional spring).

    Returns
    -------
    coords_out : (n, 3) array - copy of coords with donor positions updated.

    Notes
    -----
    M-D distance invariant exactly (geometric projection to sphere each step).
    Universal across all TMC classes. Deterministic. No templates.
    """
    if not _enabled():
        return coords
    if metal_idx is None or len(donor_indices) < 2:
        return coords

    coords = np.asarray(coords, dtype=float).copy()
    M = coords[metal_idx].copy()
    donors = list(donor_indices)
    n_d = len(donors)

    md_lengths = np.array([float(np.linalg.norm(coords[d] - M)) for d in donors])
    if (md_lengths < 1e-6).any():
        return coords

    chelate_targets = {}
    if chelate_pairs:
        for i, j in chelate_pairs:
            if i in donors and j in donors:
                d_ij = float(np.linalg.norm(coords[i] - coords[j]))
                chelate_targets[(i, j)] = d_ij

    max_iter, step_size, tol = _config()

    for iteration in range(max_iter):
        prev = coords[donors].copy()
        donor_pos = np.array([coords[d] - M for d in donors])

        forces = np.zeros((n_d, 3), dtype=float)
        for i in range(n_d):
            for j in range(i + 1, n_d):
                r_ij = donor_pos[i] - donor_pos[j]
                d_ij = float(np.linalg.norm(r_ij))
                if d_ij < 1e-6:
                    continue
                f_mag = 1.0 / (d_ij * d_ij)
                f_ij = (r_ij / d_ij) * f_mag
                forces[i] += f_ij
                forces[j] -= f_ij

        for (a, b), d_target in chelate_targets.items():
            i = donors.index(a)
            j = donors.index(b)
            r_ij = donor_pos[i] - donor_pos[j]
            d_now = float(np.linalg.norm(r_ij))
            if d_now < 1e-6:
                continue
            delta = d_now - d_target
            f_spring = -delta * (r_ij / d_now)
            forces[i] += f_spring
            forces[j] -= f_spring

        for i in range(n_d):
            r_hat = donor_pos[i] / md_lengths[i]
            radial = float(forces[i] @ r_hat) * r_hat
            forces[i] = forces[i] - radial

        new_donor_pos = donor_pos + step_size * forces
        for i in range(n_d):
            new_r = new_donor_pos[i]
            new_r_norm = float(np.linalg.norm(new_r))
            if new_r_norm < 1e-6:
                continue
            new_donor_pos[i] = new_r * (md_lengths[i] / new_r_norm)

        for i, d in enumerate(donors):
            coords[d] = M + new_donor_pos[i]

        dpos = coords[donors] - prev
        max_disp = float(np.max(np.linalg.norm(dpos, axis=1)))
        if max_disp < tol:
            break

    return coords
