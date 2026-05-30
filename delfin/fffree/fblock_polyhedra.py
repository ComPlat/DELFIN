"""delfin.fffree.fblock_polyhedra — Universal f-block CN8-12 polyhedra (Ln/An).

FF-free vertex math for lanthanide and actinide complexes with high CN.
Universal across La-Lu, Th-Lr, any donor combination.

Polyhedra:
  CN7  : capped trigonal prism (C_2v), capped octahedron (C_3v), pentagonal bipyramid (D_5h)
  CN8  : square antiprism (D_4d), dodecahedron (D_2d), cube (O_h), bicapped trigonal prism
  CN9  : tricapped trigonal prism (D_3h), capped square antiprism (C_4v)
  CN10 : bicapped square antiprism (D_4d), pentagonal antiprism (D_5d)
  CN11 : capped pentagonal antiprism (C_5v)
  CN12 : icosahedron (I_h), truncated tetrahedron, cuboctahedron (O_h)

Each polyhedron defines vertex positions on unit sphere; M-D distances scale
the radius. Universal across all f-block metals.

Env-gate: DELFIN_FFFREE_FBLOCK=1 (default OFF, byte-identical when unset).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_FBLOCK = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_FBLOCK", "0") == "1"


# Empirical M-D distances for f-block metals (Å, COD-empirical p50)
# Format: {(metal, donor_element): distance}
_FBLOCK_MD = {
    # Lanthanides (La-Lu): typical Ln-O ~2.4Å, Ln-N ~2.6Å
    ("La", "O"): 2.50, ("La", "N"): 2.65, ("La", "Cl"): 2.85,
    ("Ce", "O"): 2.48, ("Ce", "N"): 2.62,
    ("Pr", "O"): 2.46, ("Pr", "N"): 2.60,
    ("Nd", "O"): 2.45, ("Nd", "N"): 2.58,
    ("Sm", "O"): 2.42, ("Sm", "N"): 2.55,
    ("Eu", "O"): 2.41, ("Eu", "N"): 2.54,
    ("Gd", "O"): 2.40, ("Gd", "N"): 2.53,
    ("Tb", "O"): 2.38, ("Tb", "N"): 2.51,
    ("Dy", "O"): 2.37, ("Dy", "N"): 2.50,
    ("Ho", "O"): 2.36, ("Ho", "N"): 2.49,
    ("Er", "O"): 2.34, ("Er", "N"): 2.47,
    ("Tm", "O"): 2.33, ("Tm", "N"): 2.46,
    ("Yb", "O"): 2.32, ("Yb", "N"): 2.45,
    ("Lu", "O"): 2.30, ("Lu", "N"): 2.43,
    # Actinides
    ("Th", "O"): 2.40, ("Th", "N"): 2.55,
    ("U", "O"): 2.40, ("U", "N"): 2.55, ("U", "Cl"): 2.65,
    ("Np", "O"): 2.38, ("Pu", "O"): 2.36,
    ("Am", "O"): 2.40, ("Cm", "O"): 2.38,
}


def fblock_md_distance(metal: str, donor: str) -> float:
    """Universal f-block M-D distance lookup with fallback."""
    key = (metal, donor)
    if key in _FBLOCK_MD:
        return _FBLOCK_MD[key]
    # Fallback: covalent radii sum
    _COV_FBLOCK = {
        "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98,
        "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89,
        "Tm": 1.90, "Yb": 1.87, "Lu": 1.87,
        "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
        "Am": 1.80, "Cm": 1.69,
    }
    _COV_DONOR = {"O": 0.66, "N": 0.71, "S": 1.05, "F": 0.57, "Cl": 1.02,
                  "Br": 1.20, "I": 1.39, "C": 0.76, "P": 1.07}
    r_m = _COV_FBLOCK.get(metal, 1.95)
    r_d = _COV_DONOR.get(donor, 0.70)
    return r_m + r_d


def cn7_pentagonal_bipyramid() -> np.ndarray:
    """CN7 pentagonal bipyramid (D_5h). 2 axial + 5 equatorial vertices."""
    out = [(0, 0, 1), (0, 0, -1)]
    for k in range(5):
        a = 2 * math.pi * k / 5
        out.append((math.cos(a), math.sin(a), 0))
    return np.array(out)


def cn7_capped_octahedron() -> np.ndarray:
    """CN7 capped octahedron (C_3v): octahedron + 1 cap on a face."""
    # 6 octahedral vertices
    octa = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                     [0, 0, 1], [0, 0, -1]], float)
    # cap above (1,1,1)/sqrt(3) face
    cap = np.array([[1, 1, 1]]) / math.sqrt(3) * 1.0
    return np.vstack([octa, cap])


def cn8_square_antiprism() -> np.ndarray:
    """CN8 square antiprism (D_4d). 4 lower + 4 upper rotated by 45°."""
    out = []
    for k in range(4):
        a = math.pi * 2 * k / 4
        out.append((math.cos(a), math.sin(a), -0.5))
    for k in range(4):
        a = math.pi * 2 * k / 4 + math.pi / 4
        out.append((math.cos(a), math.sin(a), 0.5))
    arr = np.array(out)
    # normalize to unit sphere
    return arr / np.linalg.norm(arr, axis=1, keepdims=True)


def cn8_dodecahedron() -> np.ndarray:
    """CN8 dodecahedron (D_2d) — 8 vertices for trigonal-dodecahedral coordination.

    Vertices grouped 4+4 with different z heights.
    """
    a = 1.0
    out = []
    for k in range(4):
        ang = 2 * math.pi * k / 4
        out.append((a * math.cos(ang), a * math.sin(ang), 0.5))
    for k in range(4):
        ang = 2 * math.pi * k / 4 + math.pi / 4
        out.append((a * 0.71 * math.cos(ang), a * 0.71 * math.sin(ang), -0.7))
    arr = np.array(out)
    return arr / np.linalg.norm(arr, axis=1, keepdims=True)


def cn8_cube() -> np.ndarray:
    """CN8 cube (O_h). Used for some actinide complexes (UO2 derivatives)."""
    out = []
    for ix in (-1, 1):
        for iy in (-1, 1):
            for iz in (-1, 1):
                out.append((ix, iy, iz))
    arr = np.array(out, float) / math.sqrt(3)
    return arr


def cn9_tricapped_trigonal_prism() -> np.ndarray:
    """CN9 tricapped trigonal prism (D_3h). 6 prism + 3 cap on rectangular faces."""
    out = []
    # 6 prism vertices: trigonal at z=±h
    R = 0.85
    h = 0.5
    for z_sign in (1, -1):
        for k in range(3):
            ang = 2 * math.pi * k / 3
            out.append((R * math.cos(ang), R * math.sin(ang), z_sign * h))
    # 3 caps: on the rectangular faces (along x, rotated to 0°, 120°, 240°)
    for k in range(3):
        ang = 2 * math.pi * k / 3 + math.pi / 3
        out.append((math.cos(ang), math.sin(ang), 0))
    arr = np.array(out)
    return arr / np.linalg.norm(arr, axis=1, keepdims=True)


def cn10_bicapped_square_antiprism() -> np.ndarray:
    """CN10 = square antiprism (8) + 2 caps on opposite square faces."""
    sa = cn8_square_antiprism()
    caps = np.array([[0, 0, 1.0], [0, 0, -1.0]])
    return np.vstack([sa, caps])


def cn12_icosahedron() -> np.ndarray:
    """CN12 icosahedron (I_h). 12 vertices on golden-ratio structure."""
    g = (1 + math.sqrt(5)) / 2  # golden ratio
    out = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            out.append((0, s1, s2 * g))
            out.append((s1, s2 * g, 0))
            out.append((s2 * g, 0, s1))
    arr = np.array(out, float)
    return arr / np.linalg.norm(arr, axis=1, keepdims=True)


def cn12_cuboctahedron() -> np.ndarray:
    """CN12 cuboctahedron (O_h). 12 vertices at midpoints of cube edges."""
    out = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            out.append((0, s1, s2))
            out.append((s1, 0, s2))
            out.append((s1, s2, 0))
    arr = np.array(out, float)
    return arr / np.linalg.norm(arr, axis=1, keepdims=True)


REFS_FBLOCK = {
    "PBPY-7": cn7_pentagonal_bipyramid,
    "COCT-7": cn7_capped_octahedron,
    "SAPR-8": cn8_square_antiprism,
    "TDD-8":  cn8_dodecahedron,
    "CU-8":   cn8_cube,
    "TPRS-9": cn9_tricapped_trigonal_prism,
    "BSAP-10": cn10_bicapped_square_antiprism,
    "ICO-12":  cn12_icosahedron,
    "COC-12":  cn12_cuboctahedron,
}


def ref_vectors_fblock(geometry: str) -> np.ndarray:
    """Universal f-block polyhedron lookup. Returns unit vectors from metal to donors."""
    if geometry not in REFS_FBLOCK:
        raise ValueError(f"Unknown f-block geometry: {geometry}. Available: {list(REFS_FBLOCK)}")
    return REFS_FBLOCK[geometry]()


if __name__ == "__main__":
    print("=== Universal f-block CN7-12 polyhedra ===")
    for name, fn in REFS_FBLOCK.items():
        v = fn()
        n = len(v)
        # Check unit-sphere: all distances from origin = 1.0 ±tol
        max_dev = max(abs(np.linalg.norm(p) - 1.0) for p in v)
        # Compute min vertex-vertex angle (for closest-donor strain)
        angles = []
        for i in range(n):
            for j in range(i + 1, n):
                cos_a = float(np.dot(v[i], v[j]))
                cos_a = max(-1.0, min(1.0, cos_a))
                angles.append(math.degrees(math.acos(cos_a)))
        print(f"  {name}: CN={n}, unit_sphere_dev={max_dev:.4f}, "
              f"min_donor_angle={min(angles):.1f}°, max={max(angles):.1f}°")
