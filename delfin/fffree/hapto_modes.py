"""delfin.fffree.hapto_modes — Universal η-hapto mode enumeration + placement.

For π-bound ligands (η²-η⁸), enumerate all hapto modes that a given fragment
can present to the metal. Some ligands have multiple binding modes:
  Cyclopentadienyl (Cp): η¹ (σ-bound), η³ (allyl-like), η⁵ (full)
  Indenyl: η¹, η³, η⁵
  Pentadienyl: η¹, η³, η⁵
  Benzene/arenes: η², η⁴, η⁶
  Cyclooctatetraene (COT): η², η⁴, η⁶, η⁸
  Cycloheptatrienyl: η¹, η³, η⁵, η⁷

Plus piano-stool placement: metal at axis perpendicular to ring centroid
at appropriate M-centroid distance.

Universal across all hapto-π ligands.
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_HAPTO_MODES=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_HAPTO_MODES = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_HAPTO_MODES", "0") == "1"


# Empirical M-centroid distances per (metal, hapto-mode) (Å, CCDC p50)
HAPTO_M_CENTROID_DIST = {
    # η⁵-Cp distances
    ("Fe", 5): 1.65,    ("Mn", 5): 1.78,  ("Co", 5): 1.67,
    ("Ni", 5): 1.78,    ("Ru", 5): 1.82,  ("Os", 5): 1.85,
    ("Rh", 5): 1.85,    ("Ir", 5): 1.88,  ("V", 5): 1.90,
    ("Cr", 5): 1.80,    ("Mo", 5): 1.96,  ("W", 5): 1.96,
    # η⁶-arene
    ("Cr", 6): 1.62,    ("Fe", 6): 1.55,  ("Ru", 6): 1.69,
    ("Os", 6): 1.73,    ("Co", 6): 1.61,  ("Mo", 6): 1.85,
    # η⁴-diene (e.g. Fe(η⁴-COT)2)
    ("Fe", 4): 1.70,    ("Ru", 4): 1.83,  ("Co", 4): 1.72,
    # η³-allyl
    ("Pd", 3): 2.00,    ("Pt", 3): 2.05,  ("Ni", 3): 1.95,
    # η²-alkene/olefin
    ("Pt", 2): 2.10,    ("Pd", 2): 2.10,  ("Cu", 2): 2.20,
    ("Ag", 2): 2.40,    ("Rh", 2): 2.20,  ("Ir", 2): 2.25,
    # η⁷-cycloheptatrienyl
    ("Cr", 7): 1.40,    ("Mo", 7): 1.60,
    # η⁸-COT
    ("U", 8): 1.95,     ("Th", 8): 2.00,  ("Ce", 8): 1.97,
}


def m_centroid_distance(metal: str, eta: int) -> float:
    """Universal M-centroid distance lookup with fallback."""
    if (metal, eta) in HAPTO_M_CENTROID_DIST:
        return HAPTO_M_CENTROID_DIST[(metal, eta)]
    # Fallback: scale by eta (more atoms = lower distance per atom contribution)
    base = 1.8  # typical
    return base * (5.0 / max(eta, 2)) ** 0.5


def possible_hapto_modes_for_fragment(n_ring_atoms: int) -> List[int]:
    """For a hapto-π ring with N atoms, return possible η-modes.

    5-ring (Cp): [1, 3, 5]
    6-ring (arene): [2, 4, 6]
    7-ring (cycloheptatrienyl): [1, 3, 5, 7]
    8-ring (COT): [2, 4, 6, 8]
    Open chain (pentadienyl etc.): [1, 3, 5]

    Universal: odd-N supports η¹/η³/η⁵/.../η^N; even-N supports η²/η⁴/.../η^N.
    """
    modes = []
    if n_ring_atoms < 2:
        return modes
    is_odd = n_ring_atoms % 2 == 1
    start = 1 if is_odd else 2
    for eta in range(start, n_ring_atoms + 1, 2):
        modes.append(eta)
    return modes


def place_metal_piano_stool(
    ring_atom_positions: np.ndarray,
    metal_symbol: str,
    eta: int,
) -> np.ndarray:
    """Universal piano-stool placement: metal perpendicular to ring centroid.

    Algorithm:
      1. Compute ring centroid (mean of ring atom positions)
      2. Compute ring normal via SVD (smallest eigenvalue direction)
      3. Place metal at centroid + normal × M-centroid distance
      4. Distance scaled per (metal, eta) empirical table

    Returns: (3,) metal position vector.
    """
    if not _HAPTO_MODES:
        # Default placement at simple offset
        centroid = ring_atom_positions.mean(axis=0)
        return centroid + np.array([0, 0, 1.8])
    centroid = ring_atom_positions.mean(axis=0)
    centered = ring_atom_positions - centroid
    _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    normal = Vt[-1]
    # Ensure consistent normal direction (positive z if possible)
    if normal[2] < 0:
        normal = -normal
    distance = m_centroid_distance(metal_symbol, eta)
    return centroid + normal * distance


def enumerate_hapto_modes(
    ring_atom_positions: np.ndarray,
    metal_symbol: str,
) -> List[Tuple[int, np.ndarray]]:
    """Enumerate all possible η-modes + metal placements for a hapto ring.

    Returns list of (eta_value, metal_position) tuples for each valid mode.
    """
    if not _HAPTO_MODES:
        return []
    n = len(ring_atom_positions)
    modes = possible_hapto_modes_for_fragment(n)
    out = []
    for eta in modes:
        m_pos = place_metal_piano_stool(ring_atom_positions, metal_symbol, eta)
        out.append((eta, m_pos))
    return out


if __name__ == "__main__":
    print("=== Universal hapto-mode enumeration ===")
    for n_ring in range(2, 9):
        modes = possible_hapto_modes_for_fragment(n_ring)
        print(f"  N_ring={n_ring}: possible η-modes = {modes}")

    print("\n=== M-centroid distance examples ===")
    examples = [("Fe", 5, "ferrocene"), ("Cr", 6, "Cr(η⁶-C6H6)2"),
                ("Ru", 5, "Ru(Cp)2"), ("U", 8, "U(COT)2"), ("Cr", 7, "Cr(η⁷-Tro)")]
    for m, eta, label in examples:
        d = m_centroid_distance(m, eta)
        print(f"  {label}: M-centroid = {d:.2f} Å")

    # Placement test on benzene ring
    os.environ["DELFIN_FFFREE_HAPTO_MODES"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.hapto_modes", None)
    from delfin.fffree.hapto_modes import (
        place_metal_piano_stool, enumerate_hapto_modes
    )
    angles = np.array([i * math.pi / 3 for i in range(6)])
    benzene = np.array([[math.cos(a) * 1.4, math.sin(a) * 1.4, 0] for a in angles])
    m_pos = place_metal_piano_stool(benzene, "Cr", 6)
    print(f"\nCr above benzene ring centroid: M-pos = {m_pos}")
    modes = enumerate_hapto_modes(benzene, "Cr")
    print(f"Enumerated η-modes for Cr-benzene: {[m[0] for m in modes]}")
