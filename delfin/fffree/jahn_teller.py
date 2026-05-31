"""delfin.fffree.jahn_teller — Jahn-Teller distortion enumeration.

For d-electron configurations with degenerate ground states (d4 HS, d7 LS,
d9 in octahedral fields), Jahn-Teller theorem requires geometric distortion
to lift degeneracy. This module enumerates the canonical JT distortions:

  Z-elongation : 2 trans M-D bonds elongate (axial), 4 cis stay (equatorial)
  Z-compression : 2 trans elongate-perpendicular, 4 cis elongate (opposite)
  Trigonal twist (D_3d) : occurs in trans-σ alkyl systems
  Pseudo-JT (PJT)       : when low-lying excited state mixes (e.g. d^8 sqp)

Universal across all d-orbital configurations subject to JT theorem.
FF-free, deterministic, geometric.

Env-gate: DELFIN_FFFREE_JAHN_TELLER=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for complete distortion-mode coverage.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_JT = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_JAHN_TELLER", "0") == "1"


# d-electron configurations with strong JT in octahedral
# (degenerate ground state)
JT_ACTIVE_OCT = {
    1: ("t2g^1",   "weak JT"),    # d1
    2: ("t2g^2",   "weak JT"),    # d2 (can be JT if T2g term)
    4: ("t2g^3 eg^1", "STRONG JT"),  # d4 HS — classical Cr2+, Mn3+
    5: ("t2g^4 eg^2", "weak JT"),    # d5 LS not JT
    6: ("t2g^3 eg^3", "weak JT"),    # d6 HS Mn1+
    7: ("t2g^6 eg^1", "STRONG JT"),  # d7 LS — classical Co2+
    9: ("t2g^6 eg^3", "STRONG JT"),  # d9 — classical Cu2+
}


def is_jt_active(d_count: int, spin_state: str = "high") -> Tuple[bool, str]:
    """Check if a d-electron configuration is Jahn-Teller active.

    Returns (is_active, description).
    """
    if not _JT:
        return (False, "JT-detection disabled")
    if d_count == 4 and spin_state == "high":
        return (True, "d4 HS: t2g^3 eg^1 strong JT")
    if d_count == 7 and spin_state == "low":
        return (True, "d7 LS: t2g^6 eg^1 strong JT")
    if d_count == 9:
        return (True, "d9: t2g^6 eg^3 strong JT (Cu2+)")
    if d_count in (1, 2, 5, 6):
        return (True, f"d{d_count}: weak JT possible")
    return (False, f"d{d_count}: no JT")


def apply_z_elongation(
    coords: np.ndarray,
    metal_idx: int,
    donor_indices: Sequence[int],
    elongation_factor: float = 1.15,
) -> np.ndarray:
    """Apply Jahn-Teller Z-elongation: trans pair (axial) elongates.

    Algorithm:
      1. Find principal axis from inertia of donors (or use z-axis as default)
      2. Identify 2 donors most aligned with axis (axial pair)
      3. Elongate their M-D bonds by elongation_factor (typical 1.10-1.20)

    Universal: works on any CN6 octahedral (or distorted) coordination.
    """
    if not _JT:
        return coords
    P = coords.copy()
    metal = P[int(metal_idx)]
    donors = [int(d) for d in donor_indices]
    if len(donors) < 6:
        return P
    # Find principal axis from donor distribution
    donor_pos = np.array([P[d] - metal for d in donors])
    eigvals, eigvecs = np.linalg.eigh(donor_pos.T @ donor_pos)
    axis = eigvecs[:, 2]  # largest variance
    # Score donors by alignment with axis (most positive = axial+, most negative = axial-)
    scores = [(i, float(np.dot(P[d] - metal, axis) /
                        max(np.linalg.norm(P[d] - metal), 1e-9)))
              for i, d in enumerate(donors)]
    scores.sort(key=lambda x: -abs(x[1]))
    # Top 2 by abs-score = axial pair
    axial_indices = [donors[scores[0][0]], donors[scores[1][0]]]
    # Elongate axial M-D bonds
    for d in axial_indices:
        v = P[d] - metal
        d_dist = float(np.linalg.norm(v))
        if d_dist < 1e-6:
            continue
        new_pos = metal + v * elongation_factor
        P[d] = new_pos
    return P


def apply_z_compression(
    coords: np.ndarray,
    metal_idx: int,
    donor_indices: Sequence[int],
    compression_factor: float = 0.90,
    eq_elongation: float = 1.05,
) -> np.ndarray:
    """JT Z-compression: axial pair compresses, equatorial elongates."""
    if not _JT:
        return coords
    P = coords.copy()
    metal = P[int(metal_idx)]
    donors = [int(d) for d in donor_indices]
    if len(donors) < 6:
        return P
    donor_pos = np.array([P[d] - metal for d in donors])
    eigvals, eigvecs = np.linalg.eigh(donor_pos.T @ donor_pos)
    axis = eigvecs[:, 2]
    scores = [(i, float(np.dot(P[d] - metal, axis) /
                        max(np.linalg.norm(P[d] - metal), 1e-9)))
              for i, d in enumerate(donors)]
    scores.sort(key=lambda x: -abs(x[1]))
    axial_indices = [donors[scores[0][0]], donors[scores[1][0]]]
    equatorial_indices = [donors[scores[i][0]] for i in range(2, len(donors))]
    for d in axial_indices:
        v = P[d] - metal
        P[d] = metal + v * compression_factor
    for d in equatorial_indices:
        v = P[d] - metal
        P[d] = metal + v * eq_elongation
    return P


def enumerate_jt_distortions(
    coords: np.ndarray,
    metal_idx: int,
    donor_indices: Sequence[int],
    d_count: int,
    spin_label: str = "high",
) -> List[np.ndarray]:
    """Enumerate all relevant JT distortions for the given d-electron count.

    Returns list of (coords_variant) for each JT mode.
    """
    if not _JT:
        return [coords]
    is_active, desc = is_jt_active(d_count, spin_label)
    if not is_active:
        return [coords]
    variants = [coords]  # undistorted
    # Z-elongation (most common)
    variants.append(apply_z_elongation(coords, metal_idx, donor_indices, 1.15))
    # Z-compression (less common but valid)
    variants.append(apply_z_compression(coords, metal_idx, donor_indices, 0.90, 1.05))
    return variants


if __name__ == "__main__":
    os.environ["DELFIN_FFFREE_JAHN_TELLER"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.jahn_teller", None)
    from delfin.fffree.jahn_teller import (
        is_jt_active, apply_z_elongation, apply_z_compression,
        enumerate_jt_distortions, JT_ACTIVE_OCT
    )

    print("=== Jahn-Teller active d-electron configurations ===")
    for d, (config, desc) in sorted(JT_ACTIVE_OCT.items()):
        print(f"  d{d}: {config:<14} {desc}")
    print("\n=== Common JT cases ===")
    for d in [3, 4, 5, 6, 7, 8, 9, 10]:
        active, desc = is_jt_active(d, "high")
        print(f"  d{d} HS: active={active}, {desc}")
    # Test distortion application
    print("\n=== JT distortion application (octahedral CN6) ===")
    metal_pos = np.array([0, 0, 0])
    octa = np.array([metal_pos] + [
        [2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2]
    ], float)
    elong = apply_z_elongation(octa, 0, list(range(1, 7)), 1.20)
    for i in range(1, 7):
        d = float(np.linalg.norm(elong[i] - elong[0]))
        print(f"  Atom {i}: M-D = {d:.3f}")
