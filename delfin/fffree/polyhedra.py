"""Reference coordination polyhedra (unit vectors) + covalent radii for the
metal-FF-free builder (delfin.fffree).

Self-contained: ideal CN4/5/6 polyhedron vertex sets (tetrahedron, square planar,
trigonal bipyramid, square pyramid, octahedron, trigonal prism) and a covalent-radii
table.  Metal-donor distances use covalent-radii sums.
"""
from __future__ import annotations
import math
import numpy as np


def _norm_rows(V: np.ndarray) -> np.ndarray:
    return V / np.linalg.norm(V, axis=1, keepdims=True)


def _ref_polyhedra():
    R = {}
    t = 1 / math.sqrt(3)
    R[("CN4", "T-4 tetrahedron")] = np.array(
        [[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]])
    R[("CN4", "SP-4 square planar")] = np.array(
        [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
    R[("CN5", "TBP-5 trigonal bipyramid")] = np.array(
        [[0, 0, 1], [0, 0, -1], [1, 0, 0],
         [-0.5, math.sqrt(3) / 2, 0], [-0.5, -math.sqrt(3) / 2, 0]])
    R[("CN5", "SPY-5 square pyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [1, 1, 0.2], [1, -1, 0.2], [-1, 1, 0.2], [-1, -1, 0.2]], float))
    R[("CN6", "OC-6 octahedron")] = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    R[("CN6", "TPR-6 trigonal prism")] = _norm_rows(np.array(
        [[1, 0, 0.7], [-0.5, math.sqrt(3) / 2, 0.7], [-0.5, -math.sqrt(3) / 2, 0.7],
         [1, 0, -0.7], [-0.5, math.sqrt(3) / 2, -0.7], [-0.5, -math.sqrt(3) / 2, -0.7]]))
    return {k: _norm_rows(v) for k, v in R.items()}


REFS = _ref_polyhedra()

GEOM_BY_CN = {
    4: ["T-4 tetrahedron", "SP-4 square planar"],
    5: ["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"],
    6: ["OC-6 octahedron", "TPR-6 trigonal prism"],
}

# covalent radii (metal subset + donors); M-D = r(M) + r(D)
COV = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39, "Se": 1.20, "As": 1.19,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "La": 2.07,
}


def ref_vectors(geometry: str) -> np.ndarray:
    for (cn, shape), v in REFS.items():
        if shape == geometry:
            return v
    raise KeyError(geometry)


def md_distance(metal: str, donor: str) -> float:
    return COV.get(metal, 1.5) + COV.get(donor, 0.75)


def _kabsch_resid(P: np.ndarray, Q: np.ndarray) -> float:
    """Min mean-squared deviation aligning P onto Q by proper rotation (Kabsch,
    determinant-corrected to forbid reflection)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = P @ R.T - Q
    return float((diff * diff).sum() / len(P))


def cshm(observed_vecs, geometry: str) -> float:
    """Continuous shape measure (0 = ideal, larger = worse) of the observed
    donor unit-vectors against the ideal ``geometry`` polyhedron.  Scale-,
    rotation- and permutation-invariant: each vertex permutation of the ideal
    set is Kabsch-aligned and the smallest normalised residual is returned,
    mapped to S = 100·min_resid/obs_var.  Self-contained (no external dep);
    matches the offline calibration definition used to set the shape threshold.
    Deterministic.  Returns 0.0 for degenerate / size-mismatched input."""
    import itertools
    try:
        Q = ref_vectors(geometry)
    except KeyError:
        return 0.0
    P = np.asarray(observed_vecs, dtype=float)
    n = len(P)
    if n < 2 or n != len(Q):
        return 0.0

    def _unit_rms(a):
        rms = math.sqrt(float((a * a).sum()) / len(a))
        return a / rms if rms > 1e-9 else a

    P = _unit_rms(P)
    Q = _unit_rms(Q)
    obs_var = float((P * P).sum() / n)
    if obs_var < 1e-9:
        return 0.0
    best = float("inf")
    for perm in itertools.permutations(range(n)):
        r = _kabsch_resid(P, Q[list(perm)])
        if r < best:
            best = r
    return 100.0 * best / obs_var
