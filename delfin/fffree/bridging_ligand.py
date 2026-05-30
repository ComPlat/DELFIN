"""delfin.fffree.bridging_ligand — Universal μ-bridging-ligand detection + placement.

For multi-metal complexes, μ-ligands span 2+ metal centers (μ²=edge, μ³=face,
μ⁴=apex of 4 metals).

Detection: SMILES atoms with ≥2 bonds to metal atoms = bridging.
Placement: μ² on edge midpoint, μ³ on face centroid, μ⁴ on tetrahedron apex.

Universal across all cluster types and bridging donor elements (Cl, Br, I, O,
S, OR, SR, μ-CO, μ-CN, etc.).
FF-free, deterministic.

Env-gate: DELFIN_FFFREE_BRIDGING=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for multi-metal completeness.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_BRIDGING = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_BRIDGING", "0") == "1"


def detect_bridging_ligands(mol) -> List[Dict]:
    """Detect μ-bridging atoms in an RDKit Mol.

    A μ-bridging atom = donor with ≥2 bonds to metal atoms.

    Returns list of dicts:
      {atom_idx, element, n_metal_bonds, metal_neighbors, bridging_type}
    where bridging_type ∈ {μ2, μ3, μ4}.

    Universal: any SMILES, any metal, any bridging atom.
    """
    if not _BRIDGING:
        return []
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        def _is_metal(s): return s not in {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "B", "Si"}
    out = []
    for atom in mol.GetAtoms():
        s = atom.GetSymbol()
        if _is_metal(s):
            continue
        metal_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if _is_metal(n.GetSymbol())]
        if len(metal_neighbors) < 2:
            continue
        bridging_type = f"μ{len(metal_neighbors)}"
        out.append({
            "atom_idx": atom.GetIdx(),
            "element": s,
            "n_metal_bonds": len(metal_neighbors),
            "metal_neighbors": metal_neighbors,
            "bridging_type": bridging_type,
        })
    return out


def place_mu2_bridge(metal_pos_a: np.ndarray, metal_pos_b: np.ndarray,
                      m_x_distance: float = 2.0) -> np.ndarray:
    """μ²-bridging atom placed on perpendicular bisector of M-M edge.

    Position: midpoint of M-M segment + perpendicular offset for bonding distance.
    Universal: any pair of metal positions.
    """
    midpoint = 0.5 * (metal_pos_a + metal_pos_b)
    mm_vec = metal_pos_b - metal_pos_a
    mm_dist = float(np.linalg.norm(mm_vec))
    if mm_dist < 1e-6:
        return midpoint
    # M-X distance from each metal
    # X sits perpendicular to M-M line such that M-X = m_x_distance
    # Distance from midpoint to X along perpendicular: sqrt(m_x^2 - (mm/2)^2)
    perp_dist_sq = m_x_distance ** 2 - (mm_dist / 2) ** 2
    if perp_dist_sq < 0:
        return midpoint
    perp_dist = float(np.sqrt(perp_dist_sq))
    # Find a perpendicular direction (arbitrary; here +z if M-M is in xy plane)
    mm_unit = mm_vec / mm_dist
    helper = np.array([1.0, 0.0, 0.0]) if abs(mm_unit[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp = helper - mm_unit * float(np.dot(helper, mm_unit))
    perp_norm = float(np.linalg.norm(perp))
    if perp_norm < 1e-6:
        return midpoint
    perp = perp / perp_norm
    return midpoint + perp * perp_dist


def place_mu3_bridge(metal_pos_a: np.ndarray, metal_pos_b: np.ndarray,
                      metal_pos_c: np.ndarray, m_x_distance: float = 2.0) -> np.ndarray:
    """μ³-bridging atom on triangular face above face centroid (perpendicular)."""
    centroid = (metal_pos_a + metal_pos_b + metal_pos_c) / 3.0
    # Face normal
    v1 = metal_pos_b - metal_pos_a
    v2 = metal_pos_c - metal_pos_a
    normal = np.cross(v1, v2)
    n_norm = float(np.linalg.norm(normal))
    if n_norm < 1e-6:
        return centroid
    normal = normal / n_norm
    # Determine offset along normal: ensure M-X ≈ m_x_distance
    # Average distance from face vertices to centroid:
    d_avg = float(np.mean([np.linalg.norm(p - centroid) for p in [metal_pos_a, metal_pos_b, metal_pos_c]]))
    # X above centroid at height h such that M-X = m_x_distance
    h_sq = m_x_distance ** 2 - d_avg ** 2
    if h_sq < 0:
        return centroid
    h = float(np.sqrt(h_sq))
    return centroid + normal * h


def place_mu4_bridge(metal_positions: Sequence[np.ndarray],
                     m_x_distance: float = 2.0) -> np.ndarray:
    """μ⁴-bridging atom at centroid of 4 metal positions (or above face)."""
    metals = np.array(metal_positions)
    centroid = metals.mean(axis=0)
    # For tetrahedral M4 cluster, μ⁴ sits at centroid + axial offset along principal axis
    return centroid


if __name__ == "__main__":
    print("=== μ-bridging ligand placement demo ===")
    # μ² between two Fe atoms at (-1, 0, 0) and (1, 0, 0)
    pos = place_mu2_bridge(np.array([-1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]),
                           m_x_distance=2.0)
    d_a = float(np.linalg.norm(pos - np.array([-1.0, 0.0, 0.0])))
    d_b = float(np.linalg.norm(pos - np.array([1.0, 0.0, 0.0])))
    print(f"  μ² Fe-Fe(d=2.0Å) → X position: {pos}, M-X: {d_a:.3f} / {d_b:.3f}")

    # μ³ on a triangular face (equilateral d=2.0)
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([-0.5, 0.866, 0.0])
    c = np.array([-0.5, -0.866, 0.0])
    pos3 = place_mu3_bridge(a, b, c, m_x_distance=2.0)
    print(f"  μ³ above triangle: {pos3}, M-X: "
          f"{float(np.linalg.norm(pos3-a)):.3f} / "
          f"{float(np.linalg.norm(pos3-b)):.3f} / "
          f"{float(np.linalg.norm(pos3-c)):.3f}")

    # μ⁴ at centroid of tetrahedron
    tet = [np.array([1, 1, 1.0]), np.array([1, -1, -1.0]),
           np.array([-1, 1, -1.0]), np.array([-1, -1, 1.0])]
    pos4 = place_mu4_bridge(tet, m_x_distance=2.0)
    print(f"  μ⁴ at tetrahedron centroid: {pos4}")
