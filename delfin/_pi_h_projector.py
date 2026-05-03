"""Iter-14 Baustein 4 — Post-ETKDG/UFF Rigid-π H Projector.

Operates on already-finalised XYZ strings (post UFF, post snap, post clash-relief,
post dual-parse, post Iter-12/13 Baustein 3 coord-angle correction).  For every
detected π-system (5/6-membered planar ring of C/N/O/S — aromatic / Cp / η^n
ligand), project ring-attached H atoms onto the ring plane along the outward
in-plane radial direction (parent → centroid).

Per User-Direktive ("π-System rotiert mit entsprechenden H atomen"): every
geometric transformation upstream that rotates/tilts/shifts a π-frame must
drag attached H atoms rigidly with the ring.  Iter-9 H1 already handles
ring-snap-time projection inside ``_snap_aromatic_rings_to_plane`` (heavy-atom
snap in the molecule conformer).  Iter-12/13 B3 already handles rotation of
the X-side around a coordinated donor (BFS includes attached H).  Other
upstream π-operations (sandwich rotation, hapto-shift, ring tilt, Cp
orientation) may leave attached H atoms in stale positions.  Baustein 4 is a
post-pass that re-projects them universally on the FINAL XYZ.

Doctrine:
    - Co-evolution safe (Iter-10): post-ETKDG/UFF, no embedding disturbance.
    - Diversity safe (Iter-11): per-conformer; n_frames invariant.
    - Universal (Master-Prinzip 2): geometric ring detection only — no
      refcode, SMILES substring, element-listing, or hardcoded atom-id logic.
    - Per-H rollback: if projection introduces a NEW non-bonded clash (with
      any heavy atom or H not in the same ring) below 0.85 · Σr_vdw threshold,
      revert that single H back to its original position.
    - Opt-in via ``DELFIN_BAUSTEIN4=1``; bit-exact when disabled.

Algorithm (per XYZ):
    1. Build geometric adjacency (same as B3, vendored helpers).
    2. Detect 5/6-rings of C/N/O/S that are planar (SVD max-OOP <= 0.10 Å).
    3. For each ring atom that has exactly one H neighbour:
       a. Compute H out-of-plane distance to ring's SVD plane.
       b. If |OOP| <= 0.10 Å, skip (already in plane).
       c. Compute outward radial unit vector ``r_hat = (parent - centroid)
          / ||parent - centroid||`` projected to ring plane (in-plane component).
       d. Place new H at ``parent + d_XH · r_hat`` where d_XH is the
          element-specific X-H ideal length (C: 1.08, N: 1.01, O: 0.96,
          S: 1.34, fallback 1.08).
       e. Validate: if the new H position introduces a clash pair (with any
          atom not in the same ring's heavy atoms or attached H atoms) at
          d < 0.80 · Σr_vdw, revert this H.
    4. Re-emit XYZ with shape-identical formatting.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Vendored helpers — kept consistent with delfin/_coord_angle_corrector.py
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

_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
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


# Element-specific X-H ideal bond lengths (Å).  Same constants as the
# Iter-9 H1 in-place projector inside _snap_aromatic_rings_to_plane.
_IDEAL_XH: Dict[str, float] = {"C": 1.08, "N": 1.01, "O": 0.96, "S": 1.34}
_IDEAL_XH_FALLBACK: float = 1.08

# Out-of-plane tolerance: H atoms within this distance of the ring plane
# are deemed in-plane and untouched.  Matches Iter-9 H1.
_H_OOP_TOL: float = 0.10

# Ring-planarity filter: only rings with max(OOP) <= this value are
# treated as π-systems.  Matches the ring_tol_ang in F20 detector.
_RING_PLANAR_TOL: float = 0.10


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
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

    Identical formatting to ``_coord_angle_corrector._format_xyz`` so the
    two helpers can run in sequence without diverging line shape.
    """
    out: List[str] = []
    atom_i = 0
    trailing_newline = bool(orig_lines)  # DELFIN emit always trails newline
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(
                f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}"
            )
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if trailing_newline else "")


def _build_geometric_adjacency(
    syms: List[str], pts: np.ndarray
) -> List[List[int]]:
    """Build heavy + H bond graph from inter-atomic distances.

    Same tolerance scheme as ``_coord_angle_corrector._build_geometric_adjacency``:
    M-D bonds use Σr_cov + 0.45 Å (dative); organic bonds Σr_cov + 0.25 Å.
    """
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
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
    return nbrs


def _detect_planar_rings(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
) -> List[Tuple[Tuple[int, ...], np.ndarray, np.ndarray]]:
    """Detect 5/6-membered rings of C/N/O/S that are planar within
    ``_RING_PLANAR_TOL``.

    Returns list of (ring_atom_idxs, centroid, plane_normal).  Rings
    containing a metal atom or a hydrogen are not considered (metals
    use their own coordination plane; H is never a ring atom).
    """
    aromatic_eligible = {"C", "N", "O", "S"}
    n = len(syms)
    heavy_only_nbrs: List[List[int]] = [
        [j for j in nbrs[i]
         if syms[j] != "H" and not _is_metal_sym(syms[j])
         and syms[i] in aromatic_eligible and syms[j] in aromatic_eligible]
        for i in range(n)
    ]
    rings_canon: set = set()
    for start in range(n):
        if syms[start] not in aromatic_eligible:
            continue
        # Depth-bounded DFS from each starting atom; collect canonical sorted tuples.
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in heavy_only_nbrs[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings_canon.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))

    out: List[Tuple[Tuple[int, ...], np.ndarray, np.ndarray]] = []
    for ring in rings_canon:
        ring_pts = np.array([pts[i] for i in ring])
        centroid = ring_pts.mean(axis=0)
        centered = ring_pts - centroid
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        normal_raw = vh[-1]
        nn = float(np.linalg.norm(normal_raw))
        if nn < 1e-9:
            continue
        normal = normal_raw / nn
        ring_oop = float(np.max(np.abs(centered @ normal)))
        if ring_oop > _RING_PLANAR_TOL:
            continue  # ring not planar enough to be considered π-system
        out.append((ring, centroid, normal))
    return out


def _h_introduces_clash(
    new_h_pos: np.ndarray,
    h_idx: int,
    parent_idx: int,
    ring_set: set,
    syms: List[str],
    pts: np.ndarray,
    threshold_factor: float = 0.80,
) -> bool:
    """Return True if placing H at ``new_h_pos`` creates a NEW clash pair
    (d < ``threshold_factor`` · Σr_vdw) with any atom that is not the
    parent and not in the same ring (heavy ring atoms or H attached to
    same ring atom indices).

    Pre-bonded check: skip pairs that would be a normal X-H bond at
    typical bond distance (0.85 · Σr_cov < d < 1.10 · Σr_cov).  Pairs
    closer than 0.85·Σr_cov indicate atomic overlap (NOT a bond) and
    MUST be treated as clashes — without this lower bound, a projected
    H landing on an oxygen at d=0.4 Å would be falsely classified as a
    short bond and the projection accepted (Iter-14 Fe2-bridged-OMe
    debug case: H projected onto OMe oxygen at 0.413 Å).
    """
    n = len(syms)
    rH_v = _VDW_RADII.get("H", 1.20)
    rH_c = _COV_RADII.get("H", 0.31)
    for j in range(n):
        if j == h_idx or j == parent_idx or j in ring_set:
            continue
        d = float(np.linalg.norm(new_h_pos - pts[j]))
        rj_v = _VDW_RADII.get(syms[j], 1.7 if syms[j] != "H" else 1.20)
        rj_c = _COV_RADII.get(syms[j], 1.5 if syms[j] != "H" else 0.31)
        # Bonded range: 0.85·Σr_cov ≤ d ≤ 1.10·Σr_cov (typical X-H bonds
        # are 0.95-1.05·Σr_cov).  Anything BELOW 0.85 is overlap, not a
        # bond — must be flagged as clash.
        sum_cov = rH_c + rj_c
        if 0.85 * sum_cov <= d <= 1.10 * sum_cov:
            continue  # bonded
        if d < threshold_factor * (rH_v + rj_v):
            return True
    return False


def project_ring_h_atoms(
    syms: List[str],
    pts: np.ndarray,
) -> int:
    """Project ring-attached H atoms onto their ring plane (mutates ``pts``
    in place).  Returns count of H atoms moved.

    For each detected planar 5/6-ring of C/N/O/S, every ring atom that
    has exactly one H neighbour gets that H projected to
    ``parent + d_XH · r_hat`` where ``r_hat`` is the in-plane outward
    unit vector from centroid to parent.  Per-H rollback if a new clash
    pair appears.
    """
    n = len(syms)
    if n == 0:
        return 0

    nbrs = _build_geometric_adjacency(syms, pts)
    rings = _detect_planar_rings(syms, pts, nbrs)
    if not rings:
        return 0

    moved = 0
    for ring, centroid, normal in rings:
        ring_set = set(ring)
        for atom_idx in ring:
            sym = syms[atom_idx]
            if sym not in _IDEAL_XH and sym != "C":
                # only project H on C/N/O/S parents
                if sym not in ("C", "N", "O", "S"):
                    continue
            # Find single H neighbour of this ring atom
            h_nbrs = [j for j in nbrs[atom_idx] if syms[j] == "H"]
            if len(h_nbrs) != 1:
                continue  # CH2/NH2 etc.: not a single-H sp2 ring atom
            h_idx = h_nbrs[0]
            # Skip if this H is bonded to more than one heavy atom (shouldn't
            # happen for H, but guard anyway).
            if len(nbrs[h_idx]) != 1:
                continue
            parent_pos = pts[atom_idx]
            h_pos = pts[h_idx]
            h_oop = float(abs(np.dot(h_pos - centroid, normal)))
            if h_oop <= _H_OOP_TOL:
                continue  # already in-plane
            # In-plane outward radial from centroid to parent
            outward = parent_pos - centroid
            # Remove out-of-plane component (project to plane)
            outward -= float(np.dot(outward, normal)) * normal
            n_out = float(np.linalg.norm(outward))
            if n_out < 1e-6:
                continue
            outward_unit = outward / n_out
            d_xh = _IDEAL_XH.get(sym, _IDEAL_XH_FALLBACK)
            new_h_pos = parent_pos + d_xh * outward_unit
            # Validate: no new clash pair
            if _h_introduces_clash(
                new_h_pos, h_idx, atom_idx, ring_set, syms, pts,
                threshold_factor=0.80,
            ):
                continue  # rollback this single H
            pts[h_idx] = new_h_pos
            moved += 1
    return moved


def correct_xyz(xyz_str: str, max_passes: int = 2) -> str:
    """Apply ring-attached H projection to a single XYZ string.

    Two passes maximum: the second pass catches H atoms whose parents
    moved by the first pass via shared ring membership (rare but possible
    for fused rings).  Returns the corrected XYZ; on any error returns
    the input unchanged.
    """
    try:
        syms, pts, orig_lines = _parse_xyz(xyz_str)
        if pts.shape[0] == 0:
            return xyz_str
        pts = pts.copy()
        any_moved = False
        for _ in range(max_passes):
            moved = project_ring_h_atoms(syms, pts)
            if moved == 0:
                break
            any_moved = True
        if not any_moved:
            return xyz_str
        return _format_xyz(orig_lines, syms, pts)
    except Exception:
        return xyz_str


def correct_results(mol, results):
    """Apply ``correct_xyz`` to each (xyz, label) tuple.  Fail-safe: any
    per-XYZ exception returns the original tuple unchanged.

    Signature mirrors ``_coord_angle_corrector.correct_results`` so the
    dispatch helper in ``smiles_converter.py`` can call either one with
    the same arguments.
    """
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = correct_xyz(xyz)
            out.append((new_xyz, lbl) if len(entry) == 2 else (new_xyz,) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out
