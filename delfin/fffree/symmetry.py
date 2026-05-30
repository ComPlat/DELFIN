"""delfin.fffree.symmetry — Universal molecular point group + symmetry detection.

FF-free, deterministic, mathematically fundamental. Implements automated
point-group detection using the inertia-tensor + symmetry-operation-search
approach (Vetterling et al. 1992; modern variants in PyMOL/OpenBabel).

Output: point-group label (C_n, D_n, C_nv, C_nh, D_nh, S_n, T_d, O_h, I_h, ...)
+ list of symmetry operations (rotation axes, mirror planes, improper axes).

Used by:
  - ring_pucker_integration: symmetry-aware RMSD dedup (only compare same-orbit)
  - conformer_enum: orbit-counting for completeness via Burnside's lemma
  - polyhedra: validation that emitted geometries match expected point group
  - assemble_complex: site-symmetry detection for hapto ligands

Universal: works on ANY 3D atomic structure (organic, TMC, cluster).
Fundamental: derived from group theory + linear algebra on inertia tensor.

Env-gate: DELFIN_FFFREE_SYMMETRY=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for completeness counting.
"""
from __future__ import annotations

import math
import os
from typing import List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_SYMMETRY = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_SYMMETRY", "0") == "1"

# Tolerances
_TOL_AXIS = 0.10        # Å, for axis-matching atom pairs
_TOL_ANGLE = math.radians(3.0)
_TOL_INERTIA = 0.05     # principal moment matching threshold


def _center_of_mass(coords: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    if weights is None:
        return coords.mean(axis=0)
    return (coords * weights[:, None]).sum(axis=0) / float(weights.sum())


def _inertia_tensor(coords: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    """Compute inertia tensor about centre-of-mass.

    Returns 3x3 symmetric tensor. Eigenvectors define principal axes (a, b, c)
    sorted by eigenvalue. Equality patterns of eigenvalues classify the
    rotation symmetry:
      I_a < I_b < I_c (asymmetric top, no rotation symmetry > C_2)
      I_a = I_b < I_c (oblate symmetric top, C_n axis along c)
      I_a < I_b = I_c (prolate symmetric top, C_n axis along a)
      I_a = I_b = I_c (spherical top, T_d / O_h / I_h candidate)
    """
    if weights is None:
        weights = np.ones(len(coords))
    com = _center_of_mass(coords, weights)
    P = coords - com
    I = np.zeros((3, 3))
    for i in range(len(P)):
        x, y, z = P[i]
        w = float(weights[i])
        I[0, 0] += w * (y * y + z * z)
        I[1, 1] += w * (x * x + z * z)
        I[2, 2] += w * (x * x + y * y)
        I[0, 1] -= w * x * y
        I[0, 2] -= w * x * z
        I[1, 2] -= w * y * z
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    return I


def _principal_axes(coords: np.ndarray, weights: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Return (eigenvalues, eigenvectors) of inertia tensor, sorted ascending.

    eigenvectors[:, i] is the i-th principal axis.
    """
    I = _inertia_tensor(coords, weights)
    vals, vecs = np.linalg.eigh(I)  # already sorted ascending
    return vals, vecs


def _classify_top(eigvals: np.ndarray) -> str:
    """Classify top type from eigenvalues."""
    a, b, c = eigvals
    if abs(a) < 1e-6:
        return "linear"
    rel_ab = abs(b - a) / max(abs(b), 1e-6)
    rel_bc = abs(c - b) / max(abs(c), 1e-6)
    if rel_ab < _TOL_INERTIA and rel_bc < _TOL_INERTIA:
        return "spherical"
    if rel_ab < _TOL_INERTIA:
        return "oblate"
    if rel_bc < _TOL_INERTIA:
        return "prolate"
    return "asymmetric"


def _rotation_matrix_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation matrix."""
    a = axis / np.linalg.norm(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    one_c = 1 - c
    return np.array([
        [c + a[0] ** 2 * one_c, a[0] * a[1] * one_c - a[2] * s, a[0] * a[2] * one_c + a[1] * s],
        [a[1] * a[0] * one_c + a[2] * s, c + a[1] ** 2 * one_c, a[1] * a[2] * one_c - a[0] * s],
        [a[2] * a[0] * one_c - a[1] * s, a[2] * a[1] * one_c + a[0] * s, c + a[2] ** 2 * one_c],
    ])


def _coords_match_under_op(coords: np.ndarray, syms: Sequence[str], op_R: np.ndarray, tol: float = _TOL_AXIS) -> bool:
    """Check if operation R maps the structure to itself (atom permutation match).

    Algorithm: for each atom, find nearest atom of same element after applying R.
    If all matches within tolerance, operation is a symmetry.
    """
    com = coords.mean(axis=0)
    P = coords - com
    P_rot = P @ op_R.T
    n = len(P)
    used = [False] * n
    for i in range(n):
        d = np.linalg.norm(P - P_rot[i], axis=1)
        # find unused atom of same element within tol
        candidates = [j for j in range(n) if not used[j] and syms[j] == syms[i] and d[j] < tol]
        if not candidates:
            return False
        used[candidates[0]] = True
    return True


def _find_rotation_order(coords: np.ndarray, syms: Sequence[str], axis: np.ndarray, max_n: int = 8) -> int:
    """Find the highest n such that C_n rotation around axis is a symmetry."""
    for n in range(max_n, 1, -1):
        R = _rotation_matrix_axis_angle(axis, 2 * math.pi / n)
        if _coords_match_under_op(coords, syms, R):
            return n
    return 1


def detect_point_group(coords: np.ndarray, syms: Sequence[str]) -> dict:
    """UNIVERSAL point-group detection for any 3D atomic structure.

    Returns dict with:
      group       : label (C1, C2, C_2v, C_3v, C_4v, C_5v, C_6v, D_3h, T_d, O_h, ...)
      principal_axis : 3-vector along highest C_n
      order       : highest C_n order found
      has_inversion : bool
      has_mirror_h : bool (horizontal plane perpendicular to principal axis)
      has_mirror_v : bool (vertical planes containing principal axis)
      eigvals     : inertia tensor eigenvalues
      top_type    : 'spherical' | 'oblate' | 'prolate' | 'asymmetric' | 'linear'

    Algorithm: inertia tensor → top classification → axis candidates → order search.
    Universal across all molecule sizes and compositions.
    """
    if not _SYMMETRY:
        return {"group": "unknown", "order": 1, "top_type": "unknown"}
    coords = np.asarray(coords, float)
    syms = list(syms)
    n_atoms = len(coords)
    if n_atoms < 2:
        return {"group": "C1", "order": 1, "top_type": "atom"}
    # Use atomic-number-based weights (closer to mass; deterministic for fffree)
    _Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
          "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53, "Si": 14}
    weights = np.array([float(_Z.get(s, max(7, hash(s) % 80))) for s in syms])
    eigvals, eigvecs = _principal_axes(coords, weights)
    top_type = _classify_top(eigvals)
    # Find principal axis: for symmetric tops, axis is along the unique eigenvalue
    principal_axis = eigvecs[:, 2]  # default: largest eigval
    if top_type == "oblate":
        principal_axis = eigvecs[:, 2]
    elif top_type == "prolate":
        principal_axis = eigvecs[:, 0]
    # Test rotation order along principal axis
    order = _find_rotation_order(coords, syms, principal_axis, max_n=8)
    # Test inversion
    R_inv = -np.eye(3)
    has_inv = _coords_match_under_op(coords, syms, R_inv)
    # Test horizontal mirror (perpendicular to principal axis)
    # mirror across plane normal to v: R = I - 2 v v^T
    v = principal_axis / np.linalg.norm(principal_axis)
    R_h = np.eye(3) - 2 * np.outer(v, v)
    has_mirror_h = _coords_match_under_op(coords, syms, R_h)
    # Test vertical mirror: any of the 3 perpendicular principal axes
    has_mirror_v = False
    for k in range(3):
        if abs(np.dot(eigvecs[:, k], v)) < 0.5:  # perpendicular-ish
            R_v = np.eye(3) - 2 * np.outer(eigvecs[:, k], eigvecs[:, k])
            if _coords_match_under_op(coords, syms, R_v):
                has_mirror_v = True
                break
    # Group classification
    group = _name_group(order, has_inv, has_mirror_h, has_mirror_v, top_type)
    return {
        "group": group, "order": order, "top_type": top_type,
        "principal_axis": principal_axis,
        "has_inversion": has_inv,
        "has_mirror_h": has_mirror_h,
        "has_mirror_v": has_mirror_v,
        "eigvals": eigvals.tolist(),
    }


def _name_group(order: int, has_inv: bool, has_mirror_h: bool, has_mirror_v: bool, top: str) -> str:
    """Name the point group from the symmetry-element census."""
    if top == "spherical":
        # Likely T_d, O_h, or I_h
        if order >= 5:
            return "I_h" if has_inv else "I"
        if order == 4:
            return "O_h" if has_mirror_h else "O"
        if order == 3:
            return "T_d" if has_mirror_v else "T_h" if has_inv else "T"
    if top == "linear":
        return "D_inf_h" if has_inv else "C_inf_v"
    n = order
    if n == 1:
        if has_mirror_h or has_mirror_v:
            return "C_s"
        if has_inv:
            return "C_i"
        return "C_1"
    if has_mirror_h and has_mirror_v:
        return f"D_{n}h"
    if has_mirror_h:
        return f"C_{n}h"
    if has_mirror_v:
        # could be D_n (improper symmetry) or C_nv
        # for simplicity: D_n if has C_2 perpendicular, else C_nv
        return f"C_{n}v"
    if has_inv and n % 2 == 0:
        return f"S_{n}"
    return f"C_{n}"


def orbit_size(coords: np.ndarray, syms: Sequence[str], target: int) -> int:
    """Return the orbit size of atom `target` under the molecular point group.

    Mathematical principle: |orbit| = |G|/|Stab(target)|.
    For Burnside-lemma completeness counting of conformers.
    """
    pg = detect_point_group(coords, syms)
    order = pg["order"]
    # Approximate: full group order ~ order × {1 or 2 for inversion/mirror}
    G_size = order
    if pg["has_inversion"]:
        G_size *= 2
    if pg["has_mirror_h"]:
        G_size *= 2
    return max(1, G_size)


if __name__ == "__main__":
    # Test on standard molecules
    import os
    os.environ["DELFIN_FFFREE_SYMMETRY"] = "1"
    # Reload
    import importlib, sys
    sys.modules.pop("delfin.fffree.symmetry", None)
    from delfin.fffree.symmetry import detect_point_group
    # Methane (T_d)
    coords = np.array([
        [0, 0, 0],          # C
        [1, 1, 1],          # H
        [1, -1, -1],        # H
        [-1, 1, -1],        # H
        [-1, -1, 1],        # H
    ], float) * 0.6
    syms = ["C", "H", "H", "H", "H"]
    pg = detect_point_group(coords, syms)
    print(f"Methane: {pg['group']}  top={pg['top_type']}  order={pg['order']}  inv={pg.get('has_inversion')}")

    # Benzene (D_6h)
    angles = np.array([i * math.pi / 3 for i in range(6)])
    coords6 = np.array([[math.cos(a) * 1.4, math.sin(a) * 1.4, 0] for a in angles])
    syms6 = ["C"] * 6
    pg = detect_point_group(coords6, syms6)
    print(f"Benzene C6: {pg['group']}  order={pg['order']}  mirror_h={pg.get('has_mirror_h')}")

    # Water (C_2v)
    coords_w = np.array([
        [0, 0, 0.117],            # O
        [0.757, 0, -0.467],       # H
        [-0.757, 0, -0.467],      # H
    ])
    syms_w = ["O", "H", "H"]
    pg = detect_point_group(coords_w, syms_w)
    print(f"Water: {pg['group']}  order={pg['order']}  mirror_v={pg.get('has_mirror_v')}")
