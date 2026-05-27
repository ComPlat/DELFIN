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
    # basal vertices in 90deg-cyclic order (45,135,225,315) so the proper-rotation
    # C4 in polya_isomer_count._spy_group (1->2->3->4) is a real geometric symmetry of
    # this vertex set -> chelate cis-edge enumeration & isomer dedup are consistent with
    # placement (the index space here IS the one assemble_from_config places into).
    R[("CN5", "SPY-5 square pyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [1, 1, 0.2], [-1, 1, 0.2], [-1, -1, 0.2], [1, -1, 0.2]], float))
    R[("CN6", "OC-6 octahedron")] = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    R[("CN6", "TPR-6 trigonal prism")] = _norm_rows(np.array(
        [[1, 0, 0.7], [-0.5, math.sqrt(3) / 2, 0.7], [-0.5, -math.sqrt(3) / 2, 0.7],
         [1, 0, -0.7], [-0.5, math.sqrt(3) / 2, -0.7], [-0.5, -math.sqrt(3) / 2, -0.7]]))
    # --- High-CN polyhedra (CN7-9). Vertex INDEX ORDER is chosen to match the
    # proper-rotation generators in polya_isomer_count (_pentagonal_bipyramid_group /
    # _square_antiprism_group / _tricapped_trigonal_prism_group) so isomer dedup is a
    # real geometric symmetry of this vertex set (same contract as SPY-5/TPR-6 above).
    # CN7 PB: idx 0,1 = axial (+z,-z); idx 2-6 = equatorial regular pentagon (0,72,...,288 deg).
    _pent = [[math.cos(2 * math.pi * k / 5), math.sin(2 * math.pi * k / 5), 0.0] for k in range(5)]
    R[("CN7", "PB-7 pentagonal bipyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [0, 0, -1]] + _pent, float))
    # CN8 square antiprism: idx 0-3 = top square (0,90,180,270 deg, +z); idx 4-7 = bottom
    # square (45,135,225,315 deg, -z) -- the 45deg stagger that defines the antiprism.
    _h8 = 0.62
    _top = [[math.cos(math.pi * k / 2), math.sin(math.pi * k / 2), _h8] for k in range(4)]
    _bot = [[math.cos(math.pi * (k + 0.5) / 2 + 0.0), math.sin(math.pi * (k + 0.5) / 2), -_h8] for k in range(4)]
    # bottom at 45,135,225,315: angle = 45 + 90k
    _bot = [[math.cos(math.radians(45 + 90 * k)), math.sin(math.radians(45 + 90 * k)), -_h8] for k in range(4)]
    R[("CN8", "SQAP-8 square antiprism")] = _norm_rows(np.array(_top + _bot, float))
    # CN9 tricapped trigonal prism: idx 0-2 = top triangle (0,120,240 deg, +z); idx 3-5 =
    # bottom triangle (eclipsed, -z); idx 6-8 = caps on the 3 rectangular faces (60,180,300 deg, z=0).
    _tri_t = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), 0.7] for k in range(3)]
    _tri_b = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), -0.7] for k in range(3)]
    _caps = [[1.3 * math.cos(math.radians(60 + 120 * k)), 1.3 * math.sin(math.radians(60 + 120 * k)), 0.0] for k in range(3)]
    R[("CN9", "TTP-9 tricapped trigonal prism")] = _norm_rows(np.array(_tri_t + _tri_b + _caps, float))
    return {k: _norm_rows(v) for k, v in R.items()}


REFS = _ref_polyhedra()

GEOM_BY_CN = {
    4: ["T-4 tetrahedron", "SP-4 square planar"],
    5: ["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"],
    6: ["OC-6 octahedron", "TPR-6 trigonal prism"],
    7: ["PB-7 pentagonal bipyramid"],
    8: ["SQAP-8 square antiprism"],
    9: ["TTP-9 tricapped trigonal prism"],
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
    """Min mean-squared deviation aligning P onto Q by a proper rotation
    (Kabsch, determinant-corrected to forbid reflection)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = P @ R.T - Q
    return float((diff * diff).sum() / len(P))


def cshm(observed_vecs, geometry: str) -> float:
    """Continuous shape measure (0 = ideal, larger = worse) of observed donor
    unit-vectors against the ideal ``geometry`` polyhedron — scale-, rotation-
    and permutation-invariant (Kabsch over all vertex permutations, S =
    100·min_resid/obs_var).  Self-contained, deterministic; returns 0.0 for
    degenerate / size-mismatched input.  Used by the self-gate to reject
    catastrophic-coordination-shape outlier builds (#39)."""
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
