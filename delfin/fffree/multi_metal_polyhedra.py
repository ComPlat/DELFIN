"""delfin.fffree.multi_metal_polyhedra — Universal multi-metal cluster polyhedra.

FF-free, deterministic vertex math for multi-metal complexes:
  M_2 dimers (M-M bond + per-metal coordination sphere)
  M_3 triangles (Os3, Ru3 carbonyl-type clusters)
  M_4 tetrahedra ([Co4(CO)12], [Fe4S4] cubane)
  M_6 octahedra ([Re6S8] cluster, [Mo6Cl12])
  M_8 cubes (some Fe-S cubanes)

Each cluster polyhedron defines:
  vertex_positions : (n_metal, 3) array of M positions on the metal-metal polyhedron
  m_m_distance     : M-M bond length (COD-empirical per metal pair)
  faces            : list of M-index tuples for μ3-bridging ligand placement
  edges            : list of (i, j) for μ2-bridging ligand placement

Universal: any combination of N metal atoms forming a polyhedral cluster.
Fundamental: pure geometric construction from polyhedron math.
NO FF: vertex coordinates from regular-polyhedron formulas.

Env-gate: DELFIN_FFFREE_MULTI_METAL=1 (default OFF, byte-identical when unset).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_MULTI_METAL = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_MULTI_METAL", "0") == "1"

# Empirical M-M bond lengths (Å) for common metal pairs (COD-empirical p50)
_MM_BONDS = {
    ("Mn", "Mn"): 2.92, ("Re", "Re"): 3.02, ("Fe", "Fe"): 2.51,
    ("Co", "Co"): 2.49, ("Ru", "Ru"): 2.85, ("Os", "Os"): 2.88,
    ("Rh", "Rh"): 2.66, ("Ir", "Ir"): 2.70, ("Pd", "Pd"): 2.60,
    ("Pt", "Pt"): 2.62, ("Mo", "Mo"): 2.20,  # quadruple bond
    ("Cu", "Cu"): 2.46, ("Ag", "Ag"): 2.97, ("Au", "Au"): 2.78,
    ("Ni", "Ni"): 2.45, ("Cr", "Cr"): 2.49, ("W", "W"): 2.30,
}


def m_m_distance(M_a: str, M_b: str) -> float:
    """Universal M-M bond length. Symmetric lookup.

    Falls back to sum of covalent radii * 0.95 (typical metal contraction)
    if pair not in empirical table.
    """
    key1 = (M_a, M_b); key2 = (M_b, M_a)
    if key1 in _MM_BONDS:
        return _MM_BONDS[key1]
    if key2 in _MM_BONDS:
        return _MM_BONDS[key2]
    # Fallback: covalent radii sum × contraction factor
    _COV_R = {"Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
              "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
              "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47,
              "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
              "La": 2.07, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51,
              "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32}
    r_a = _COV_R.get(M_a, 1.5)
    r_b = _COV_R.get(M_b, 1.5)
    return 0.95 * (r_a + r_b)


def m2_dimer(metal_a: str, metal_b: str) -> Dict:
    """M_2 dimer: 2 metals along z-axis at distance d_MM/2 from center.

    Coordination spheres extend perpendicular to the M-M axis for each metal.
    """
    d = m_m_distance(metal_a, metal_b)
    return {
        "metals": [metal_a, metal_b],
        "vertex_positions": np.array([[0, 0, -d / 2], [0, 0, d / 2]]),
        "mm_axis": np.array([0, 0, 1]),
        "mm_distance": d,
        "edges": [(0, 1)],
        "faces": [],
    }


def m3_triangle(metals: List[str]) -> Dict:
    """M_3 equilateral triangle (e.g., Os3(CO)12, Ru3(CO)12).

    Universal D_3h symmetry placement.
    """
    if len(metals) != 3:
        raise ValueError("m3_triangle requires exactly 3 metals")
    d = max(m_m_distance(metals[i], metals[(i + 1) % 3]) for i in range(3))
    R = d / math.sqrt(3)
    angles = [0, 2 * math.pi / 3, 4 * math.pi / 3]
    pos = np.array([[R * math.cos(a), R * math.sin(a), 0] for a in angles])
    return {
        "metals": metals,
        "vertex_positions": pos,
        "mm_distance": d,
        "edges": [(0, 1), (1, 2), (2, 0)],
        "faces": [(0, 1, 2)],  # one triangular face for μ3-ligand
    }


def m4_tetrahedron(metals: List[str]) -> Dict:
    """M_4 tetrahedral cluster ([Co4(CO)12], [Fe4S4] cubane core).

    Universal T_d symmetry. 4 vertices, 6 edges (μ2 sites), 4 faces (μ3 sites).
    """
    if len(metals) != 4:
        raise ValueError("m4_tetrahedron requires exactly 4 metals")
    d = max(m_m_distance(metals[i], metals[j]) for i in range(4) for j in range(i + 1, 4))
    # Regular tetrahedron with edge length d
    # Vertices: (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1) scaled
    a = d / math.sqrt(8)  # edge = sqrt(8)·a for vertex coords above
    pos = a * np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
    ], float)
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    return {
        "metals": metals,
        "vertex_positions": pos,
        "mm_distance": d,
        "edges": edges,
        "faces": faces,
    }


def m6_octahedron(metals: List[str]) -> Dict:
    """M_6 octahedral cluster ([Re6S8(L)6], [Mo6Cl12]).

    Universal O_h symmetry. 6 vertices, 12 edges (μ2 sites), 8 triangular faces.
    """
    if len(metals) != 6:
        raise ValueError("m6_octahedron requires exactly 6 metals")
    d = max(m_m_distance(metals[i], metals[j]) for i in range(6) for j in range(i + 1, 6))
    a = d / math.sqrt(2)  # edge length = sqrt(2) for unit octahedron at coords ±1
    pos = a * np.array([
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
    ], float)
    # Edges (12): adjacent vertex pairs only (perpendicular axes, dot ≈ 0).
    # Antipodal pairs (dot ≈ -a²) excluded.
    edges = []
    for i in range(6):
        for j in range(i + 1, 6):
            if abs(np.dot(pos[i], pos[j])) < 0.5 * a * a:
                edges.append((i, j))
    # Faces (8 triangular)
    faces = [(0, 2, 4), (0, 2, 5), (0, 3, 4), (0, 3, 5),
             (1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5)]
    return {
        "metals": metals,
        "vertex_positions": pos,
        "mm_distance": d,
        "edges": edges,
        "faces": faces,
    }


def m4_square_planar(metals: List[str]) -> Dict:
    """M_4 square-planar cluster (D_4h, e.g. Pd_4 acetate squares).

    Universal D_4h symmetry. 4 vertices on a square, 4 edges, 1 face (square).
    """
    if len(metals) != 4:
        raise ValueError("m4_square_planar requires exactly 4 metals")
    d = max(m_m_distance(metals[i], metals[(i + 1) % 4]) for i in range(4))
    R = d / math.sqrt(2)
    angles = [i * math.pi / 2 for i in range(4)]
    pos = np.array([[R * math.cos(a), R * math.sin(a), 0] for a in angles])
    return {
        "metals": metals,
        "vertex_positions": pos,
        "mm_distance": d,
        "edges": [(0, 1), (1, 2), (2, 3), (3, 0)],
        "faces": [(0, 1, 2, 3)],
    }


CLUSTER_BUILDERS = {
    "M2_dimer": m2_dimer,
    "M3_triangle": m3_triangle,
    "M4_tetrahedron": m4_tetrahedron,
    "M4_square": m4_square_planar,
    "M6_octahedron": m6_octahedron,
}


def build_cluster(cluster_type: str, metals: List[str]) -> Dict:
    """Universal cluster builder dispatcher.

    cluster_type: one of CLUSTER_BUILDERS keys.
    metals: list of metal symbols, length must match the cluster size.
    """
    builder = CLUSTER_BUILDERS.get(cluster_type)
    if builder is None:
        raise ValueError(f"Unknown cluster: {cluster_type}")
    if cluster_type == "M2_dimer":
        return builder(metals[0], metals[1])
    return builder(metals)


if __name__ == "__main__":
    print("=== Universal multi-metal cluster polyhedra ===")
    for name, _ in CLUSTER_BUILDERS.items():
        if name == "M2_dimer":
            c = build_cluster(name, ["Fe", "Fe"])
        elif name == "M3_triangle":
            c = build_cluster(name, ["Os", "Os", "Os"])
        elif name == "M4_tetrahedron":
            c = build_cluster(name, ["Co", "Co", "Co", "Co"])
        elif name == "M4_square":
            c = build_cluster(name, ["Pd", "Pd", "Pd", "Pd"])
        elif name == "M6_octahedron":
            c = build_cluster(name, ["Re"] * 6)
        print(f"  {name}: {len(c['vertex_positions'])} metals, "
              f"d_MM={c['mm_distance']:.2f}Å, "
              f"{len(c['edges'])} edges (μ2 sites), "
              f"{len(c['faces'])} faces (μ3/μ4 sites)")
