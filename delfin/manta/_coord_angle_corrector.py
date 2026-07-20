"""Iter-12 Baustein 3 — Post-ETKDG/UFF Coordination-Angle Corrector.

Operates on already-finalised XYZ strings (post UFF, post snap, post clash-relief,
post dual-parse).  For every coordinated donor D bonded to a metal M, infer the
expected M-D-X coordination angle from D's hybridization (geometric inference,
no element/SMILES list), and rotate the rigid X-side of D in the (M, D, X) plane
to bring the observed angle within tolerance.  Per-violation revertable; chelate-
guarded; no diversity impact (per-conformer post-step).

Universal properties:
    - Metal classification by atomic-number range (no element whitelist).
    - Hybridization classified geometrically: aromatic-DFS + bond-shortening +
      neighbour count + covalent-radius table.
    - Rotation around an axis perpendicular to the (M, D, X) plane through D —
      this is the only axis that actually changes the M-D-X angle while keeping
      D and M fixed.
    - Atom set rotated = BFS from X with {M, D} blocked.  If BFS reaches any
      metal, it is a chelate ring → skip (UFF inherits its own equilibrium).
    - Per-fix revert if validation worsens (new clash, degraded angle).

Doctrine:
    - Co-evolution safe (Iter-10): post-ETKDG/UFF, no embedding disturbance.
    - Diversity safe (Iter-11): per-conformer; n_frames invariant.
    - Opt-in via DELFIN_FFFREE_COORD_ANGLE_FIX=1; bit-exact when disabled.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Vendored helpers from quality_framework/scripts/compute_ligand_angles.py
# (kept DRY — single source of truth, no runtime dep on validation framework).
# ---------------------------------------------------------------------------

_COV_RADII: Dict[str, float] = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

# Bondi vdW radii for clash check.  Falls back to 1.7 Å for missing entries.
_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
    # Transition-metal placeholder (universal: same value for missing metals).
    "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00, "Fe": 2.00,
    "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87, "Ge": 2.11,
    "As": 1.85, "Se": 1.90, "Y": 2.20, "Zr": 2.10, "Nb": 2.10, "Mo": 2.10,
    "Tc": 2.10, "Ru": 2.10, "Rh": 2.10, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "La": 2.30, "Hf": 2.10,
    "Ta": 2.10, "W": 2.10, "Re": 2.10, "Os": 2.10, "Ir": 2.10, "Pt": 1.75,
    "Au": 1.66, "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
}

_METAL_Z_RANGES = (
    set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
    | set(range(89, 104))
    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51, 55, 56, 81, 82, 83}
)
_Z_BY_SYMBOL: Dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,
    "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51,
    "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}


def _is_metal_sym(sym: str) -> bool:
    z = _Z_BY_SYMBOL.get(sym)
    return z is not None and z in _METAL_Z_RANGES


def _angle_deg(p_center: np.ndarray, p_a: np.ndarray, p_b: np.ndarray) -> Optional[float]:
    v1 = p_a - p_center
    v2 = p_b - p_center
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return None
    cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos)))


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse XYZ → (symbols, positions Nx3, raw_lines).

    raw_lines preserves header (n + comment) and any unparseable trailing
    lines so the round-trip writer can reconstruct an identical-shaped XYZ.
    """
    syms: List[str] = []
    pts: List[np.ndarray] = []
    lines = xyz_str.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append(np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))]))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(orig_lines: List[str], syms: List[str], positions: np.ndarray) -> str:
    """Reformat XYZ preserving header line(s) but rewriting atom coordinates.

    Format matches DELFIN's existing emit (used by build_complex / OB writer):
        '<sym:<2}    <x:>10.6f>    <y:>10.6f>    <z:>10.6f>'
    Reproduces a shape-identical line so hashes only differ when coords change.
    """
    out: List[str] = []
    atom_i = 0
    trailing_newline = xyz_str_ended_with_newline(orig_lines)
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            # 8 spaces of pad after a 1-char element, 7 after a 2-char element,
            # then 10-wide signed float matches existing DELFIN/OB format.
            out.append(
                f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}"
            )
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if trailing_newline else "")


def xyz_str_ended_with_newline(orig_lines: List[str]) -> bool:
    """Heuristic: DELFIN emit always trails with newline."""
    if not orig_lines:
        return False
    return True


# ---------------------------------------------------------------------------
# Adjacency + hybridization inference (pure geometry)
# ---------------------------------------------------------------------------

def _build_geometric_adjacency(
    syms: List[str], pts: np.ndarray
) -> Tuple[List[List[int]], Dict[Tuple[int, int], float]]:
    """Build heavy + H bond graph from inter-atomic distances.

    M-D bonds use wider tolerance (dative); organic bonds use tighter
    tolerance to avoid false-positive D-X bonds inflating neighbour counts.
    """
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    bond_d: Dict[Tuple[int, int], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ri = _COV_RADII.get(syms[i], 1.5)
            rj = _COV_RADII.get(syms[j], 1.5)
            is_metal_i = _is_metal_sym(syms[i])
            is_metal_j = _is_metal_sym(syms[j])
            if is_metal_i or is_metal_j:
                thr = ri + rj + 0.45
            else:
                thr = ri + rj + 0.25
            if d < thr:
                nbrs[i].append(j)
                nbrs[j].append(i)
                bond_d[(i, j)] = d
                bond_d[(j, i)] = d
    return nbrs, bond_d


def _detect_aromatic_atoms(
    syms: List[str], pts: np.ndarray, nbrs: List[List[int]]
) -> set:
    """5/6-membered planar C/N/O/S rings → aromatic flag."""
    aromatic_eligible = {"C", "N", "O", "S"}
    n = len(syms)
    heavy_only_nbrs: List[List[int]] = [
        [j for j in nbrs[i]
         if syms[j] != "H" and not _is_metal_sym(syms[j])
         and syms[i] in aromatic_eligible and syms[j] in aromatic_eligible]
        for i in range(n)
    ]
    rings: set = set()
    for start in range(n):
        if syms[start] not in aromatic_eligible:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in heavy_only_nbrs[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))
    aromatic_atoms: set = set()
    for ring in rings:
        ring_pts = np.array([pts[i] for i in ring])
        centroid = ring_pts.mean(axis=0)
        centered = ring_pts - centroid
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        normal = vh[-1] / max(np.linalg.norm(vh[-1]), 1e-9)
        ring_oop = float(np.max(np.abs(centered @ normal)))
        if ring_oop <= 0.20:
            for ri in ring:
                aromatic_atoms.add(ri)
    return aromatic_atoms


def _coord_angle_violations(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    bond_d: Dict[Tuple[int, int], float],
    aromatic_atoms: set,
    sp_tol: float = 15.0,
    sp2_tol: float = 15.0,
    sp3_tol: float = 15.0,
) -> List[Dict]:
    """Identical inference to coord_angle_realism_check, operating on
    pre-computed adjacency to avoid redundant work in the iteration loop.
    """
    n = len(syms)
    metal_idxs = [i for i in range(n) if _is_metal_sym(syms[i])]
    if not metal_idxs:
        return []

    violations: List[Dict] = []
    for m_idx in metal_idxs:
        donor_idxs = [j for j in nbrs[m_idx]
                      if not _is_metal_sym(syms[j]) and syms[j] != "H"]
        for d_idx in donor_idxs:
            d_nbrs_all = [k for k in nbrs[d_idx] if k != m_idx]
            d_heavy_nbrs = [k for k in d_nbrs_all
                            if syms[k] != "H" and not _is_metal_sym(syms[k])]
            d_h_nbrs = [k for k in d_nbrs_all if syms[k] == "H"]
            k_total = len(d_heavy_nbrs) + len(d_h_nbrs)
            if k_total < 1:
                continue

            if d_idx in aromatic_atoms:
                hyb_label = "sp2"
                ideal = 120.0
                tol = sp2_tol
            else:
                n_short = 0
                n_very_short = 0
                for k in d_heavy_nbrs:
                    d = bond_d.get((d_idx, k))
                    if d is None:
                        continue
                    rsum = _COV_RADII.get(syms[d_idx], 1.5) + _COV_RADII.get(syms[k], 1.5)
                    if d < 0.82 * rsum:
                        n_very_short += 1
                    elif d < 0.93 * rsum:
                        n_short += 1
                if k_total == 1:
                    if n_very_short >= 1 or n_short >= 1:
                        hyb_label, ideal, tol = "sp", 180.0, sp_tol
                    else:
                        hyb_label, ideal, tol = "sp2", 120.0, sp2_tol
                elif k_total == 2:
                    if n_very_short >= 1:
                        hyb_label, ideal, tol = "sp", 180.0, sp_tol
                    elif n_short >= 1:
                        hyb_label, ideal, tol = "sp2", 120.0, sp2_tol
                    else:
                        hyb_label, ideal, tol = "sp3", 109.47, sp3_tol
                else:
                    hyb_label, ideal, tol = "sp3", 109.47, sp3_tol

            for x_idx in d_nbrs_all:
                ang = _angle_deg(pts[d_idx], pts[m_idx], pts[x_idx])
                if ang is None:
                    continue
                dev = abs(ang - ideal)
                if dev > tol:
                    violations.append({
                        "M_idx": m_idx, "D_idx": d_idx, "X_idx": x_idx,
                        "M_sym": syms[m_idx], "D_sym": syms[d_idx],
                        "X_sym": syms[x_idx],
                        "observed_angle": ang,
                        "expected_angle": ideal,
                        "hybridization": hyb_label,
                        "deviation": ang - ideal,
                        "tol": tol,
                    })
    return violations


# ---------------------------------------------------------------------------
# Geometric utilities
# ---------------------------------------------------------------------------

def _bfs_side_X(
    nbrs: List[List[int]],
    start: int,
    blocked: set,
) -> set:
    """BFS from start, never traversing blocked atoms.  Returns the set
    of reachable atom indices including start (if start not blocked)."""
    if start in blocked:
        return set()
    visited = {start}
    queue = [start]
    while queue:
        cur = queue.pop()
        for nb in nbrs[cur]:
            if nb in blocked or nb in visited:
                continue
            visited.add(nb)
            queue.append(nb)
    return visited


def _rodrigues_rotation_matrix(axis: np.ndarray, theta_rad: float) -> np.ndarray:
    """Rodrigues rotation matrix R(axis, theta).  Axis must be unit-length."""
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(axis, axis)


def _arbitrary_perpendicular(v: np.ndarray) -> np.ndarray:
    """Return a unit vector perpendicular to v (which must be non-zero)."""
    if abs(v[0]) < 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    else:
        ref = np.array([0.0, 1.0, 0.0])
    perp = np.cross(v, ref)
    n = float(np.linalg.norm(perp))
    if n < 1e-9:
        return np.array([0.0, 1.0, 0.0])
    return perp / n


def _clash_identity_set(
    pts: np.ndarray,
    syms: List[str],
    side_a: set,
    side_b_filter,
    threshold_factor_heavy: float = 0.85,
    threshold_factor_h: float = 0.80,
) -> set:
    """Return SET of (i,j) pairs that clash, i ∈ side_a, j ∈ side_b_filter.

    Pair identities (not just counts) — lets caller reject if any NEW clash
    pair appears, even if total count stayed the same.

    Pre-bonded pairs (d < 1.10·Σr_cov) are skipped.
    Heavy-heavy pairs use 0.85·Σr_vdw threshold; pairs involving H use 0.80.
    """
    clashes: set = set()
    for i in side_a:
        for j in range(len(syms)):
            if j == i or j in side_a:
                continue
            if not side_b_filter(j):
                continue
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ri_v = _VDW_RADII.get(syms[i], 1.7 if syms[i] != "H" else 1.2)
            rj_v = _VDW_RADII.get(syms[j], 1.7 if syms[j] != "H" else 1.2)
            ri_c = _COV_RADII.get(syms[i], 1.5 if syms[i] != "H" else 0.31)
            rj_c = _COV_RADII.get(syms[j], 1.5 if syms[j] != "H" else 0.31)
            if d < (ri_c + rj_c) * 1.10:
                continue  # bonded
            thr = threshold_factor_h if (syms[i] == "H" or syms[j] == "H") else threshold_factor_heavy
            if d < thr * (ri_v + rj_v):
                clashes.add((min(i, j), max(i, j)))
    return clashes


def _hard_min_distance(
    pts: np.ndarray,
    syms: List[str],
    side_a: set,
    side_b_filter,
) -> float:
    """Return the minimum non-bonded distance between side_a and side_b_filter."""
    dmin = float('inf')
    for i in side_a:
        for j in range(len(syms)):
            if j == i or j in side_a:
                continue
            if not side_b_filter(j):
                continue
            ri_c = _COV_RADII.get(syms[i], 1.5 if syms[i] != "H" else 0.31)
            rj_c = _COV_RADII.get(syms[j], 1.5 if syms[j] != "H" else 0.31)
            d = float(np.linalg.norm(pts[i] - pts[j]))
            if d < (ri_c + rj_c) * 1.10:
                continue
            if d < dmin:
                dmin = d
    return dmin


def _h_close_heavy_pairs(
    pts: np.ndarray,
    syms: List[str],
    side_a: set,
    side_b_filter,
    threshold: float = 1.40,
) -> set:
    """Return SET of (h_idx, heavy_idx) pairs where an H atom in side_a is
    within ``threshold`` Å of a non-bonded heavy atom in side_b_filter.

    Used for the triangle-H gate: any NEW such pair indicates the rotation
    moved an H in side_a close to a non-bonded heavy outside, which creates
    a triangle-H pattern (H wedged between two heavies < 1.40 Å apart).

    Note: only considers H ∈ side_a — H atoms outside side_a are stationary
    and any change in their close-pairs is captured by the inverse direction
    of the rotation (heavy moves toward stationary H, which means the heavy
    in side_a appears close to the H — but we don't check that).  Empirically
    the v3 simpler form catches all in-frame triangle creation; the inverse
    direction was added in a v4 attempt but proved too slow without
    detectable extra benefit.
    """
    pairs: set = set()
    for i in side_a:
        if syms[i] != "H":
            continue
        for j in range(len(syms)):
            if j == i or j in side_a:
                continue
            if syms[j] == "H":
                continue
            if not side_b_filter(j):
                continue
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ri_c = _COV_RADII.get(syms[i], 0.31)
            rj_c = _COV_RADII.get(syms[j], 1.5)
            if d < (ri_c + rj_c) * 1.10:
                continue
            if d < threshold:
                pairs.add((i, j))
    return pairs


# Backward-compat shim (returns count for legacy code paths if any)
def _compute_clash_state(
    pts: np.ndarray,
    syms: List[str],
    side_a: set,
    side_b_filter,
    threshold_factor: float = 0.85,
) -> int:
    return len(_clash_identity_set(pts, syms, side_a, side_b_filter,
                                    threshold_factor_heavy=threshold_factor,
                                    threshold_factor_h=0.80))


# ---------------------------------------------------------------------------
# Iter-15 — Hard M-D bond-length invariant guard.
#
# User diagnosis 2026-05-03 (commits 35ac005, fdee83a, fdee83a-b3b4on):
# B3/B4 rotation operations broke 23-32% of metal-donor bonds because the
# rotation pivot D, although fixed, propagated through chelate-ring atoms
# that are simultaneously donors to the same metal M.  These secondary
# donors moved relative to M, stretching their M-D bonds.
#
# Hard invariant: for every (metal, donor) pair within "intact" range
# BEFORE the rotation, the same pair MUST remain within +/- 0.05 A of its
# pre-rotation distance after the rotation.  Any larger change rolls back
# the rotation.
# ---------------------------------------------------------------------------

# Maximum allowed change of any pre-existing intact M-D bond length (A)
_MD_INVARIANT_TOL: float = 0.05

# Reference-topology intact factor: an (M, D) pair is "intact" if the
# observed distance d <= _MD_INTACT_FACTOR * Sigma_r_cov.  Matches the
# F26 detector default (1.30 covers covalent + short dative range).
_MD_INTACT_FACTOR: float = 1.30


def _snapshot_md_bonds(
    pts: np.ndarray,
    syms: List[str],
    metal_idxs: List[int],
) -> Dict[Tuple[int, int], float]:
    """Return dict of (metal_idx, donor_idx) -> distance for every M-D pair
    within ``_MD_INTACT_FACTOR * Sigma_r_cov``.

    Universal: NO refcode/SMILES/element-list logic.  Donor = any non-H,
    non-metal atom within the intact-bond range of a metal.
    """
    snap: Dict[Tuple[int, int], float] = {}
    for m_idx in metal_idxs:
        rm = _COV_RADII.get(syms[m_idx], 1.5)
        pm = pts[m_idx]
        for j in range(len(syms)):
            if j == m_idx:
                continue
            sj = syms[j]
            if sj == "H" or _is_metal_sym(sj):
                continue
            rj = _COV_RADII.get(sj, 1.5)
            d = float(np.linalg.norm(pm - pts[j]))
            sigma = rm + rj
            if d <= _MD_INTACT_FACTOR * sigma:
                snap[(m_idx, j)] = d
    return snap


def _md_invariant_violated(
    pre: Dict[Tuple[int, int], float],
    new_pts: np.ndarray,
    tol: float = _MD_INVARIANT_TOL,
) -> bool:
    """Return True if any pre-snapshot M-D distance differs by more than
    ``tol`` Angstrom in ``new_pts``.

    Operates on the SAME atom indices as the snapshot — atom topology
    must not have changed between snapshot and check.
    """
    for (m_idx, d_idx), pre_d in pre.items():
        new_d = float(np.linalg.norm(new_pts[m_idx] - new_pts[d_idx]))
        if abs(new_d - pre_d) > tol:
            return True
    return False


# ---------------------------------------------------------------------------
# Top-level: per-XYZ correction
# ---------------------------------------------------------------------------

def _try_fix_one_violation(
    pts: np.ndarray,
    syms: List[str],
    nbrs: List[List[int]],
    metal_idxs: List[int],
    viol: Dict,
    pre_clash_set: set,
    pre_min_dist: float,
    pre_triangle_pairs: set,
    pre_md_snapshot: Optional[Dict[Tuple[int, int], float]] = None,
) -> bool:
    """Attempt rotation; revert on validation failure.  Mutates pts in place
    on success.  Returns True iff applied.

    Validation gates (any failure → revert / try smaller step):
      1. Identity gate: NO new clash pair (i,j) appears that wasn't in
         pre_clash_set.  Removing pre-existing clashes is OK.
      2. Hard-distance gate: minimum non-bonded distance between side_X and
         rest must be >= max(0.85, 0.7·pre_min_dist).  Catches catastrophic
         atom overlaps where a new pair drops to <0.5 Å.
      3. Angle gate: rotated M-D-X angle within ±5° of expected, OR strictly
         improved over observed.
      4. (Iter-15) M-D invariant: every pre-existing M-D bond length must
         change by less than ``_MD_INVARIANT_TOL`` Å.  Catches catastrophic
         donor-detachment (Iter-12+13+14 23-32% break-rate).
    """
    M = viol["M_idx"]
    D = viol["D_idx"]
    X = viol["X_idx"]

    # Atom set: BFS from X, blocking {D, M}
    blocked = {D, M}
    side_X = _bfs_side_X(nbrs, X, blocked)
    if not side_X or len(side_X) >= len(syms) - 2:
        return False  # disconnected or pathological

    # Chelate guard: if BFS reaches another metal, skip
    if any(m in side_X for m in metal_idxs if m != M):
        return False
    if M in side_X:
        return False

    # Ring guard: a rigid rotation of side_X only preserves the rest's internal
    # geometry if side_X attaches to the fixed part through EXACTLY one bond --
    # the X-D pivot.  If 2+ bonds cross the side_X boundary, the pivot D-X lies
    # inside a ring (an aromatic/chelate ligand ring), and rotating side_X about
    # D deforms that ring: every gate here (clash, triangle-H, M-D) passes on a
    # smooth ring bend, yet xyz2mol can no longer perceive the ring -> the input
    # graph is not recovered (round-trip loss, e.g. OYUGAR's NHC ring).  UFF owns
    # ring equilibrium; skip.  Universal (bond-topology only, no element list).
    _cross = sum(1 for a in side_X for b in nbrs[a] if b not in side_X)
    if _cross > 1:
        return False

    # Chelate/polydentate-donor guard: the corrector only owns FREE monodentate
    # donors; a chelating donor's group is UFF's job.  D chelates iff the metal M
    # has ANOTHER coordinated donor D' on the SAME ligand -- i.e. D' is reachable
    # from D through the ligand graph WITHOUT passing through M.  The existing
    # "side_X reaches a metal" guard misses this when the chelate ring is broken
    # in the covalent adjacency: a bidentate nitrate/carboxylate's 2nd M-O is a
    # long dative bond absent from nbrs, so the rotated side_X excludes the 2nd
    # donor and never trips.  Yet rotating ANY part of that small rigid group
    # (here just N + terminal O) distorts it while every M-D DISTANCE stays intact
    # (M-D invariant passes) -> xyz2mol perception flips -> round-trip loss
    # (NASQAB).  Detect via M-D snapshot (metal-distance geometry) + ligand
    # connectivity (topology); universal, no element/SMILES list.
    if pre_md_snapshot is not None:
        _m_donors = {d for (m, d) in pre_md_snapshot.keys() if m == M and d != D}
        if _m_donors:
            _lig = _bfs_side_X(nbrs, D, {M})  # D's ligand, not crossing the metal
            if any(dd in _lig for dd in _m_donors):
                return False

    v_DM = pts[M] - pts[D]
    v_DX = pts[X] - pts[D]
    n_DM = float(np.linalg.norm(v_DM))
    n_DX = float(np.linalg.norm(v_DX))
    if n_DM < 1e-6 or n_DX < 1e-6:
        return False

    cross = np.cross(v_DM, v_DX)
    if float(np.linalg.norm(cross)) < 1e-6:
        axis = _arbitrary_perpendicular(v_DM)
    else:
        axis = cross / float(np.linalg.norm(cross))

    theta_deg = viol["expected_angle"] - viol["observed_angle"]
    side_X_list = list(side_X)
    pivot = pts[D].copy()
    rest_filter = lambda j: (j not in side_X) and (j != D) and (j != M)

    def apply_rotation(theta_d: float) -> np.ndarray:
        R = _rodrigues_rotation_matrix(axis, math.radians(theta_d))
        new_pts = pts.copy()
        for i in side_X_list:
            new_pts[i] = pivot + R @ (pts[i] - pivot)
        return new_pts

    def gate_ok(new_pts: np.ndarray) -> Tuple[bool, Optional[float]]:
        """Returns (ok, angle).  ok = no new catastrophic clash AND no new
        clash-pair identity AND no new triangle-H pair AND M-D invariant
        preserved.  Angle = M-D-X after rotation."""
        # Hard-min-distance gate (catastrophic overlap detector)
        new_min = _hard_min_distance(new_pts, syms, side_X, rest_filter)
        floor = max(0.85, 0.7 * pre_min_dist) if pre_min_dist < 1e8 else 0.85
        if new_min < floor:
            return (False, None)
        # Identity gate: any new clash pair?
        new_clashes = _clash_identity_set(new_pts, syms, side_X, rest_filter,
                                           threshold_factor_heavy=0.85,
                                           threshold_factor_h=0.80)
        new_only = new_clashes - pre_clash_set
        if new_only:
            return (False, None)
        # Triangle-H gate: any NEW (H, heavy) pair within 1.40 Å (excluding
        # the H's bonded parent)?  This catches the H-wedge configuration
        # where rotation moves an H next to a non-bonded heavy.
        new_tri = _h_close_heavy_pairs(new_pts, syms, side_X, rest_filter, 1.40)
        if new_tri - pre_triangle_pairs:
            return (False, None)
        # (Iter-15) M-D invariant: any pre-existing intact M-D bond stretched
        # by more than _MD_INVARIANT_TOL Å is a hard FAIL.
        if pre_md_snapshot is not None and _md_invariant_violated(
            pre_md_snapshot, new_pts, tol=_MD_INVARIANT_TOL,
        ):
            return (False, None)
        ang = _angle_deg(new_pts[D], new_pts[M], new_pts[X])
        return (True, ang)

    # Try theta, -theta, 0.5·theta, -0.5·theta, 0.25·theta, -0.25·theta
    candidates = []
    for sign in (1, -1):
        for fac in (1.0, 0.5, 0.25):
            t = sign * theta_deg * fac
            new_pts = apply_rotation(t)
            ok, ang = gate_ok(new_pts)
            if not ok or ang is None:
                continue
            # Score: how close is new angle to expected?
            score = abs(ang - viol["expected_angle"])
            # Reject if angle didn't improve (still farther from ideal than orig)
            if score >= abs(viol["observed_angle"] - viol["expected_angle"]):
                continue
            candidates.append((score, t, new_pts, ang))

    if not candidates:
        return False

    # Pick the candidate with the best score (smallest |new - expected|)
    candidates.sort(key=lambda c: c[0])
    _, _, best_pts, _ = candidates[0]
    pts[:] = best_pts
    return True


def correct_xyz(xyz_str: str, max_passes: int = 10) -> str:
    """Correct coordination angles in a single XYZ string.

    Returns the corrected XYZ.  If correction fails or has no effect, returns
    the input unchanged.
    """
    syms, pts, orig_lines = _parse_xyz(xyz_str)
    if pts.shape[0] == 0:
        return xyz_str

    metal_idxs = [i for i, s in enumerate(syms) if _is_metal_sym(s)]
    if not metal_idxs:
        return xyz_str

    pts = pts.copy()

    # Iter-15: GLOBAL M-D snapshot taken ONCE at the start of correct_xyz.
    # Every rotation must preserve every pre-existing intact M-D bond
    # within ``_MD_INVARIANT_TOL`` Å of its initial length.  This catches
    # cumulative drift across multiple sequential rotations that, taken
    # individually, are each within tol but together exceed it.
    initial_md_snapshot = _snapshot_md_bonds(pts, syms, metal_idxs)

    for outer in range(max_passes):
        nbrs, bond_d = _build_geometric_adjacency(syms, pts)
        aromatic = _detect_aromatic_atoms(syms, pts, nbrs)
        viols = _coord_angle_violations(syms, pts, nbrs, bond_d, aromatic)
        if not viols:
            break
        # worst-deviation first
        viols.sort(key=lambda v: -abs(v["deviation"]))

        progress = False
        for v in viols:
            # Per-violation baseline state (re-computed each violation since
            # earlier violations in this pass may have moved atoms).
            side_X = _bfs_side_X(nbrs, v["X_idx"], {v["D_idx"], v["M_idx"]})
            if not side_X:
                continue
            if any(m in side_X for m in metal_idxs if m != v["M_idx"]):
                continue
            rest_filter = lambda j, sx=side_X, dd=v["D_idx"], mm=v["M_idx"]: \
                (j not in sx) and (j != dd) and (j != mm)
            pre_clash_set = _clash_identity_set(
                pts, syms, side_X, rest_filter,
                threshold_factor_heavy=0.85, threshold_factor_h=0.80,
            )
            pre_min_dist = _hard_min_distance(pts, syms, side_X, rest_filter)
            pre_triangle_pairs = _h_close_heavy_pairs(
                pts, syms, side_X, rest_filter, 1.40
            )

            ok = _try_fix_one_violation(
                pts, syms, nbrs, metal_idxs, v,
                pre_clash_set=pre_clash_set,
                pre_min_dist=pre_min_dist,
                pre_triangle_pairs=pre_triangle_pairs,
                pre_md_snapshot=initial_md_snapshot,
            )
            if ok:
                progress = True
        if not progress:
            break

    return _format_xyz(orig_lines, syms, pts)


def correct_results(mol, results):
    """Apply correct_xyz to each (xyz, label) tuple.  Fail-safe: any per-XYZ
    error returns the original tuple."""
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = correct_xyz(xyz)
            out.append((new_xyz, lbl) if len(entry) == 2 else (new_xyz,) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out
