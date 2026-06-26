"""Targeted F19 fixer — sp3 H-atom tetrahedrality correction.

Detect sp3 heavy atoms (C, N, O, P, S, ...) whose attached H atoms
participate in (X-A-H) or (H-A-H) angles deviating from the ideal
tetrahedral 109.5° by more than ``tolerance_deg`` (default 10°), then
repair by rotating ONLY the offending H atoms around an in-place axis
through the parent A, preserving the A-H bond length exactly.

Architectural philosophy
------------------------
DIFFERENT from a universal post-processor (e.g. ``_pi_h_projector``):
this fixer is a targeted detector+fixer for a single chemistry pattern.
It fires only on detected violations and touches only H atoms — never a
heavy atom and never a metal.  Topology preservation is guaranteed by
construction: rotating H around the parent A in 3-D leaves every bond
involving non-H atoms unchanged, and the A-H bond length is preserved
(rotation around the parent A).

Doctrine
--------
- Topology-safe: metal atoms NEVER moved, non-H atoms NEVER moved.
- Bond-length preserving: every A-H stays at its original A-H distance.
- Per-H rollback: any H move that introduces a NEW non-bonded clash
  (heavy-H or H-H below 0.85·Σr_vdw) is reverted.
- Conservative tolerance (10° default) — only large violations fire.
- No env-flag inside module: the caller (pipeline wrapper) decides.
- Lazy imports for RDKit and scipy; numpy at module-level (project norm).
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants (vendored from delfin/_pi_h_projector.py for parity)
# ---------------------------------------------------------------------------

_IDEAL_TETRA_DEG: float = 109.471  # arccos(-1/3) in degrees
_IDEAL_XH_DEFAULT: float = 1.09  # used only when an embedded A-H is degenerate

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


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


# ---------------------------------------------------------------------------
# XYZ helpers (parity with _pi_h_projector / _coord_angle_corrector)
# ---------------------------------------------------------------------------


def _parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    syms: List[str] = []
    pts: List[np.ndarray] = []
    lines = xyz_str.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append(np.array([float(m.group(2)), float(m.group(3)),
                                 float(m.group(4))]))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(orig_lines: List[str], syms: List[str],
                positions: np.ndarray) -> str:
    out: List[str] = []
    atom_i = 0
    trailing_newline = bool(orig_lines)
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if trailing_newline else "")


# ---------------------------------------------------------------------------
# Geometric helpers
# ---------------------------------------------------------------------------


def _compute_angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle (a-b-c) in degrees, with b as the vertex."""
    v1 = a - b
    v2 = c - b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return float("nan")
    cosv = float(np.dot(v1, v2) / (n1 * n2))
    cosv = max(-1.0, min(1.0, cosv))
    return math.degrees(math.acos(cosv))


def _build_geometric_adjacency(syms: List[str],
                                pts: np.ndarray) -> List[List[int]]:
    """Distance-based adjacency, matching _pi_h_projector convention.

    M-D dative bonds use Σr_cov + 0.45 Å; organic bonds Σr_cov + 0.25 Å.
    """
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(pts[i] - pts[j]))
            ri = _COV_RADII.get(syms[i], 1.5)
            rj = _COV_RADII.get(syms[j], 1.5)
            metal = _is_metal_sym(syms[i]) or _is_metal_sym(syms[j])
            thr = ri + rj + (0.45 if metal else 0.25)
            if d < thr:
                nbrs[i].append(j)
                nbrs[j].append(i)
    return nbrs


def _h_introduces_clash(new_h_pos: np.ndarray, h_idx: int, parent_idx: int,
                        syms: List[str], pts: np.ndarray,
                        sibling_idxs: Optional[set] = None,
                        threshold_factor: float = 0.70,
                        ignore_h_to_h: bool = False) -> bool:
    """True iff placing H at ``new_h_pos`` creates a NEW clash with any
    atom that is neither the H itself nor its parent nor a sibling
    (sibling = another atom bonded to the same parent — geminal pair).

    Geminal H pairs in sp3 atoms legitimately come close to one another
    at tetrahedral geometry (e.g. CH4 H-H = 1.78 Å ≈ 0.74·Σr_vdw), so
    geminal sibling pairs are excluded from the clash check.  Real
    inter-molecular clashes (H landing inside another fragment) still
    trigger via the heavy-atom and non-sibling-H pathway.

    Bonded-range exclusion mirrors _pi_h_projector: a distance inside
    [0.85·Σr_cov, 1.10·Σr_cov] indicates a normal X-H bond and is not a
    clash.  Distances below 0.85·Σr_cov indicate atomic overlap.

    ``threshold_factor`` defaults to 0.70 (≈1.68 Å for H-H) so that
    geometries with tightly-packed but chemically correct H atoms are
    not flagged, while genuine atomic overlap is.
    """
    n = len(syms)
    rH_v = _VDW_RADII.get("H", 1.20)
    rH_c = _COV_RADII.get("H", 0.31)
    siblings = sibling_idxs if sibling_idxs is not None else set()
    for j in range(n):
        if j == h_idx or j == parent_idx or j in siblings:
            continue
        if ignore_h_to_h and syms[j] == "H":
            continue
        d = float(np.linalg.norm(new_h_pos - pts[j]))
        rj_v = _VDW_RADII.get(syms[j], 1.7 if syms[j] != "H" else 1.20)
        rj_c = _COV_RADII.get(syms[j], 1.5 if syms[j] != "H" else 0.31)
        sum_cov = rH_c + rj_c
        if 0.85 * sum_cov <= d <= 1.10 * sum_cov:
            continue  # bonded
        if d < threshold_factor * (rH_v + rj_v):
            return True
    return False


# ---------------------------------------------------------------------------
# RDKit-driven sp3 detection
# ---------------------------------------------------------------------------


def _sp3_heavy_atom_indices(mol) -> List[int]:
    """Return indices of non-H sp3 atoms in ``mol``.

    Uses RDKit hybridization annotation; falls back to empty list if the
    mol has none (caller skips silently).
    """
    from rdkit import Chem  # lazy
    sp3 = Chem.HybridizationType.SP3
    out: List[int] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        if atom.GetHybridization() == sp3:
            out.append(atom.GetIdx())
    return out


# ---------------------------------------------------------------------------
# FULL umbrella reconstruction (graph-driven, replaces per-angle rotation)
# ---------------------------------------------------------------------------
#
# Root-of-partial (legacy ``fix_sp3_h_tetrahedrality`` below): it rotated each
# H independently to make ONE (other-A-H) angle 109.5°, never reconstructing a
# joint tetrahedral umbrella, and it derived H-count from the *embedded*
# geometry — which in the unsanitized fallback is degenerate (AddHs overlaps
# sp3 H, so geometric adjacency over-counts H and the count-4 gate skips the
# centre entirely).  The functions below fix BOTH: they read connectivity from
# the RDKit ``mol`` bond graph (correct H-count per centre) and place EVERY H
# of a centre jointly on the ideal tetrahedral / pyramidal / bent sites left
# open by the heavy neighbours, at the original A-H bond length.
#
# Universal: element + bond-graph + hybridization only.  Topology-safe (only H
# move; A-H lengths preserved).  Deterministic (fixed tables, sorted iteration,
# no RNG/clock).  Never emits non-finite coordinates (NaN guards + per-H
# rollback to the original position).


_IDEAL_TETRA_RAD: float = math.acos(-1.0 / 3.0)  # 109.471° in radians


def _ortho_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Two deterministic unit vectors perpendicular to unit ``axis``."""
    helper = (np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9
              else np.array([0.0, 1.0, 0.0]))
    u = helper - float(np.dot(helper, axis)) * axis
    nu = float(np.linalg.norm(u))
    if nu < 1e-9:
        u = np.array([0.0, 1.0, 0.0])
        nu = 1.0
    u = u / nu
    v = np.cross(axis, u)
    nv = float(np.linalg.norm(v))
    if nv < 1e-9:
        return u, np.array([0.0, 0.0, 1.0])
    return u, v / nv


def _umbrella_dirs_one_heavy(axis_in: np.ndarray, n_h: int,
                             apex_rad: float) -> List[np.ndarray]:
    """``n_h`` outward unit directions making angle ``apex_rad`` with
    ``-axis_in`` (the single heavy bond), equispaced around ``axis_in``.

    ``axis_in`` = unit vector parent -> the single heavy neighbour.  The
    umbrella opens AWAY from that bond.  Phase 0 is deterministic.
    """
    a = axis_in
    u, v = _ortho_basis(a)
    cos_t = math.cos(apex_rad)
    sin_t = math.sin(apex_rad)
    dirs: List[np.ndarray] = []
    k_max = max(1, n_h)
    for k in range(n_h):
        phi = 2.0 * math.pi * k / k_max
        d = (-a) * cos_t + (math.cos(phi) * u + math.sin(phi) * v) * sin_t
        nd = float(np.linalg.norm(d))
        if nd > 1e-9:
            dirs.append(d / nd)
    return dirs


def _umbrella_dirs_two_heavy(a: np.ndarray, b: np.ndarray, n_h: int,
                             apex_rad: float) -> List[np.ndarray]:
    """``n_h`` outward unit directions for a centre with TWO heavy bonds
    ``a``,``b`` (unit, parent->neighbour).  H straddle the plane that
    bisects the two bonds; H-A-H opening = ``apex_rad``.
    """
    bis = a + b
    nbis = float(np.linalg.norm(bis))
    if nbis < 1e-9:
        # collinear heavy bonds — pick an arbitrary perpendicular plane
        u, w = _ortho_basis(a)
        bis_n = u
        out_axis = w
    else:
        bis_n = bis / nbis
        out_axis = np.cross(a, b)
        no = float(np.linalg.norm(out_axis))
        if no < 1e-9:
            u, w = _ortho_basis(bis_n)
            out_axis = w
        else:
            out_axis = out_axis / no
    half = apex_rad / 2.0
    cos_h = math.cos(half)
    sin_h = math.sin(half)
    d1 = -bis_n * cos_h + out_axis * sin_h
    d2 = -bis_n * cos_h - out_axis * sin_h
    out: List[np.ndarray] = []
    for d in (d1, d2):
        nd = float(np.linalg.norm(d))
        if nd > 1e-9:
            out.append(d / nd)
    return out[:n_h]


def _umbrella_dirs_three_heavy(a: np.ndarray, b: np.ndarray, c: np.ndarray
                               ) -> List[np.ndarray]:
    """Single outward direction for a centre with THREE heavy bonds
    (methine): opposite the sum of the three bond vectors."""
    s = a + b + c
    d = -s
    nd = float(np.linalg.norm(d))
    if nd < 1e-9:
        return [-a]
    return [d / nd]


# Apex angle (H-A-H / X-A-H opening) per heavy-atom element.  Tetrahedral for
# C/Si/N+/quaternary; mild pyramidal/bent narrowing for the lone-pair-bearing
# neutral heteroatoms; wider for the heavy/soft donors.  All in radians.
_APEX_DEG: Dict[str, float] = {
    "C": 109.471, "Si": 109.471, "Ge": 109.471, "Sn": 109.471,
    "N": 107.0, "P": 93.5, "As": 92.0,
    "O": 108.5, "S": 92.0, "Se": 91.0,
    "B": 109.471,
}
_APEX_DEFAULT_DEG: float = 109.471


def _apex_rad_for(sym: str, formal_charge: int, n_heavy: int, n_h: int
                  ) -> float:
    """Ideal apex angle (radians) for an sp3 centre.

    Charged centres with no lone pair (e.g. ammonium N+, oxonium O+) take
    the full tetrahedral angle; neutral lone-pair heteroatoms keep the mild
    VSEPR narrowing.  A fully-substituted centre (n_heavy + n_h == 4) is
    tetrahedral regardless.
    """
    if n_heavy + n_h >= 4:
        return _IDEAL_TETRA_RAD
    if sym == "N" and formal_charge > 0:
        return _IDEAL_TETRA_RAD
    if sym == "O" and formal_charge > 0:
        return _IDEAL_TETRA_RAD
    if sym == "S" and formal_charge > 0:
        return _IDEAL_TETRA_RAD
    return math.radians(_APEX_DEG.get(sym, _APEX_DEFAULT_DEG))


def _ideal_h_dirs_full(sym: str, formal_charge: int,
                       parent_pos: np.ndarray,
                       heavy_pos: List[np.ndarray],
                       n_h: int) -> List[np.ndarray]:
    """Return ``n_h`` outward unit H directions for an sp3 centre given the
    fixed heavy-neighbour positions.  Joint umbrella, not per-angle."""
    bonds: List[np.ndarray] = []
    for hp in heavy_pos:
        v = np.asarray(hp) - parent_pos
        nv = float(np.linalg.norm(v))
        if nv > 1e-9:
            bonds.append(v / nv)
    n_heavy = len(bonds)
    apex = _apex_rad_for(sym, formal_charge, n_heavy, n_h)
    if n_heavy == 0:
        tet = [np.array([1.0, 1.0, 1.0]), np.array([1.0, -1.0, -1.0]),
               np.array([-1.0, 1.0, -1.0]), np.array([-1.0, -1.0, 1.0])]
        return [v / math.sqrt(3.0) for v in tet][:n_h]
    if n_heavy == 1:
        return _umbrella_dirs_one_heavy(bonds[0], n_h, apex)
    if n_heavy == 2:
        return _umbrella_dirs_two_heavy(bonds[0], bonds[1], n_h, apex)
    return _umbrella_dirs_three_heavy(bonds[0], bonds[1], bonds[2])[:n_h]


def _greedy_match(h_vecs: List[np.ndarray],
                  target_dirs: List[np.ndarray]) -> List[int]:
    """Assign each existing H vector (parent->H) to the closest unused
    target direction (min angle).  Deterministic greedy; returns a list
    ``mapping`` with ``target_dirs[mapping[i]]`` chosen for H ``i``."""
    nt = len(target_dirs)
    used = [False] * nt
    mapping: List[int] = []
    for hv in h_vecs:
        nh = float(np.linalg.norm(hv))
        best_k = -1
        best_cos = -2.0
        for k in range(nt):
            if used[k]:
                continue
            nd = float(np.linalg.norm(target_dirs[k]))
            if nh < 1e-9 or nd < 1e-9:
                cos = -1.0
            else:
                cos = float(np.dot(hv, target_dirs[k]) / (nh * nd))
            if cos > best_cos:
                best_cos = cos
                best_k = k
        if best_k < 0:
            best_k = next((k for k in range(nt) if not used[k]), 0) \
                if nt else -1
        if 0 <= best_k < nt:
            used[best_k] = True
        mapping.append(best_k)
    return mapping


def _rotate_dirs_about_axis(dirs: List[np.ndarray], axis: np.ndarray,
                            angle: float) -> List[np.ndarray]:
    """Rodrigues rotation of each unit dir about unit ``axis`` by ``angle``."""
    a = axis
    c = math.cos(angle)
    s = math.sin(angle)
    out: List[np.ndarray] = []
    for d in dirs:
        out.append(d * c + np.cross(a, d) * s
                   + a * float(np.dot(a, d)) * (1.0 - c))
    return out


def _count_dir_clashes(parent_pos: np.ndarray, dirs: List[np.ndarray],
                       h_idxs: List[int], d_xh: float, atom_idx: int,
                       syms: List[str], pts: np.ndarray,
                       siblings: set) -> int:
    """Count how many of the placed H (parent + d_xh·dir) introduce a clash."""
    c = 0
    for d in dirs:
        cand = parent_pos + d_xh * d
        # Ignore H-H during phase probing: other H atoms are still at their
        # stale (pre-rebuild) positions, so an H-H "clash" is uninformative.
        if _h_introduces_clash(cand, -1, atom_idx, syms, pts,
                               sibling_idxs=siblings,
                               threshold_factor=0.70,
                               ignore_h_to_h=True):
            c += 1
    return c


def _best_h_dirs_full(sym: str, formal_charge: int, parent_pos: np.ndarray,
                      heavy_pos: List[np.ndarray], h_idxs: List[int],
                      atom_idx: int, syms: List[str], pts: np.ndarray,
                      siblings: set) -> List[np.ndarray]:
    """Ideal H directions for a centre, with rigid phase optimisation for the
    single-heavy (methyl/aminyl) umbrella so a crowded centre rotates its H
    umbrella to the least-clashing phase instead of being rolled back.

    For n_heavy != 1 the geometry is phase-fixed (methylene plane / methine
    apex / free tetrahedron), so the base template is returned directly.
    """
    base = _ideal_h_dirs_full(sym, formal_charge, parent_pos, heavy_pos,
                              len(h_idxs))
    # Phase optimisation only helps the single-heavy umbrella with >=2 H
    # (CH3 / NH3+ / SiH3 / PH2-on-single-C etc.); other geometries are rigid.
    if len(heavy_pos) != 1 or len(h_idxs) < 2 or not base:
        return base
    axis = np.asarray(heavy_pos[0]) - parent_pos
    na = float(np.linalg.norm(axis))
    if na < 1e-9:
        return base
    axis = axis / na
    # Representative A-H length for the clash probe (mean of existing H).
    d_xh = 0.0
    for h in h_idxs:
        d_xh += float(np.linalg.norm(pts[h] - parent_pos))
    d_xh = (d_xh / len(h_idxs)) if h_idxs else _IDEAL_XH_DEFAULT
    if d_xh < 1e-6:
        d_xh = _IDEAL_XH_DEFAULT
    n_probes = 12
    period = 2.0 * math.pi / 3.0  # 3-fold symmetric umbrella
    best_dirs = base
    best_clash = _count_dir_clashes(parent_pos, base, h_idxs, d_xh,
                                    atom_idx, syms, pts, siblings)
    best_dev = _phase_dev(base, h_idxs, parent_pos, pts)
    for k in range(1, n_probes):
        phase = period * (k / n_probes)
        cand = _rotate_dirs_about_axis(base, axis, phase)
        clash = _count_dir_clashes(parent_pos, cand, h_idxs, d_xh,
                                   atom_idx, syms, pts, siblings)
        dev = _phase_dev(cand, h_idxs, parent_pos, pts)
        if (clash, dev) < (best_clash, best_dev):
            best_clash = clash
            best_dev = dev
            best_dirs = cand
    return best_dirs


def _phase_dev(dirs: List[np.ndarray], h_idxs: List[int],
               parent_pos: np.ndarray, pts: np.ndarray) -> float:
    """Sum of angles between existing H vectors and their greedy-matched
    target dirs — used as a deterministic tiebreak so a clash-free umbrella
    prefers the phase closest to the input (minimal rotation)."""
    h_vecs = [pts[h] - parent_pos for h in h_idxs]
    mapping = _greedy_match(h_vecs, dirs)
    tot = 0.0
    for li, hv in enumerate(h_vecs):
        k = mapping[li]
        if k < 0 or k >= len(dirs):
            continue
        nh = float(np.linalg.norm(hv))
        nd = float(np.linalg.norm(dirs[k]))
        if nh < 1e-9 or nd < 1e-9:
            continue
        cv = max(-1.0, min(1.0, float(np.dot(hv, dirs[k]) / (nh * nd))))
        tot += math.degrees(math.acos(cv))
    return tot


def fix_sp3_h_tetrahedrality_full(xyz: str, mol,
                                  tolerance_deg: float = 10.0
                                  ) -> Tuple[str, dict]:
    """FULL graph-driven sp3-H umbrella reconstruction.

    For every sp3 heavy centre in ``mol`` (C/N/O/P/S/Si/B/..., including
    protonated / charged centres), rebuild ALL of its attached H atoms onto
    the ideal tetrahedral / pyramidal / bent sites implied by the FIXED
    heavy-neighbour frame, at the original per-H A-H bond length.

    Connectivity comes from the RDKit ``mol`` bond graph (so a methylene C
    is known to carry exactly two H even when the embedded geometry overlaps
    them).  Atom ordering must match the XYZ — guaranteed by the fallback
    wiring, where both the XYZ and ``mol`` come from the same post-AddHs mol.

    Guarantees: only H move; A-H length preserved (within FP epsilon); no
    non-finite output (NaN guard + per-centre rollback); deterministic.

    Returns ``(new_xyz, report)`` with the same report keys as the legacy
    fixer plus ``"centres_rebuilt"``.
    """
    report: dict = {
        "sp3_atoms_checked": 0,
        "violations_detected": 0,
        "h_atoms_moved": 0,
        "centres_rebuilt": 0,
        "max_angle_correction_deg": 0.0,
        "topology_preserved": True,
    }
    if mol is None:
        return xyz, report
    try:
        syms, pts, orig_lines = _parse_xyz(xyz)
        if pts.shape[0] == 0:
            return xyz, report
        pts = pts.copy()
        n_atoms = pts.shape[0]
        n_mol = mol.GetNumAtoms()

        sp3_heavy = _sp3_heavy_atom_indices(mol)
        report["sp3_atoms_checked"] = len(sp3_heavy)
        if not sp3_heavy:
            return xyz, report

        # ---- Pass 1: rebuild every sp3 centre's H umbrella jointly. ----
        # Heavy-atom clashes block the centre (real overlap with a non-H);
        # H-H is IGNORED here because the OTHER centres' H are still at their
        # stale positions — an H-H "clash" against a stale H is uninformative.
        # Pass 2 re-validates H-H once every centre has been relocated.
        max_corr = 0.0
        # Track which H atoms were moved + their parent + old position, for
        # the pass-2 re-validation/rollback.
        moved_h: Dict[int, Tuple[int, np.ndarray]] = {}
        centre_h: Dict[int, List[int]] = {}
        for atom_idx in sorted(sp3_heavy):
            if atom_idx >= n_mol or atom_idx >= n_atoms:
                continue
            if _is_metal_sym(syms[atom_idx]):
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            nbr_idxs = [nb.GetIdx() for nb in atom.GetNeighbors()]
            # Mol-graph H neighbours (true H-count); skip metal-bonded H.
            h_idxs: List[int] = []
            for j in nbr_idxs:
                if j >= n_atoms:
                    continue
                if syms[j] != "H":
                    continue
                h_idxs.append(j)
            if not h_idxs:
                continue
            heavy_idxs = [j for j in nbr_idxs
                          if j < n_atoms and syms[j] != "H"
                          and not _is_metal_sym(syms[j])]
            # Defensive: a sane sp3 centre has heavy+H <= 4 (C/N/O/F/B) /
            # <= 6 for hypervalent; skip pathological perception artefacts.
            cap = 4 if syms[atom_idx] in ("C", "N", "O", "F", "B") else 6
            if len(heavy_idxs) + len(h_idxs) > cap:
                continue

            parent_pos = pts[atom_idx]
            heavy_pos = [pts[j] for j in heavy_idxs]
            siblings = set(nbr_idxs) | {atom_idx}
            target_dirs = _best_h_dirs_full(
                syms[atom_idx], int(atom.GetFormalCharge()),
                parent_pos, heavy_pos, h_idxs, atom_idx, syms, pts,
                siblings,
            )
            if len(target_dirs) < len(h_idxs):
                continue  # no clean template — leave centre untouched

            h_vecs = [pts[h] - parent_pos for h in h_idxs]
            mapping = _greedy_match(h_vecs, target_dirs)

            old_pos = {h: pts[h].copy() for h in h_idxs}
            new_pos: Dict[int, np.ndarray] = {}
            this_corr = 0.0
            bad = False
            for li, h in enumerate(h_idxs):
                k = mapping[li]
                if k < 0 or k >= len(target_dirs):
                    bad = True
                    break
                d_ah = float(np.linalg.norm(old_pos[h] - parent_pos))
                if d_ah < 1e-6:
                    d_ah = _IDEAL_XH_DEFAULT
                cand = parent_pos + d_ah * target_dirs[k]
                if not np.all(np.isfinite(cand)):
                    bad = True
                    break
                new_pos[h] = cand
                cur = old_pos[h] - parent_pos
                nc = float(np.linalg.norm(cur))
                nt = float(np.linalg.norm(target_dirs[k]))
                if nc > 1e-9 and nt > 1e-9:
                    cv = max(-1.0, min(1.0,
                                       float(np.dot(cur, target_dirs[k])
                                             / (nc * nt))))
                    this_corr = max(this_corr, math.degrees(math.acos(cv)))
            if bad:
                continue
            # Skip if every H is already within tolerance (preserve geometry
            # / avoid needless re-emit) — joint check, byte-stable on good in.
            if this_corr < tolerance_deg:
                continue

            # Apply jointly; per-H rollback only on a genuine HEAVY OVERLAP
            # (factor 0.55 ≈ true atomic interpenetration, not mere vdW
            # proximity — a correct umbrella in a crowded centre legitimately
            # places H in soft vdW contact with a neighbour heavy atom; the
            # downstream VSEPR / rotamer relief resolves the soft contact).
            # H-H deferred to pass-2.  M / heavy atoms never move.
            placed: List[int] = []
            for h in h_idxs:
                if _h_introduces_clash(new_pos[h], h, atom_idx, syms, pts,
                                       sibling_idxs=siblings,
                                       threshold_factor=0.55,
                                       ignore_h_to_h=True):
                    continue
                pts[h] = new_pos[h]
                moved_h[h] = (atom_idx, old_pos[h])
                placed.append(h)
            if placed:
                centre_h[atom_idx] = placed
                if this_corr > max_corr:
                    max_corr = this_corr

        # ---- Pass 2: roll back any H that landed in a GENUINE atomic
        # overlap (heavy OR H, factor 0.55 = real interpenetration) now that
        # all centres are relocated.  Soft vdW contacts are left for the
        # downstream VSEPR / rotamer relief.  Deterministic order. ----
        for h in sorted(moved_h.keys()):
            parent_idx, old = moved_h[h]
            sib = {parent_idx, h}
            try:
                pa = mol.GetAtomWithIdx(parent_idx)
                sib |= {nb.GetIdx() for nb in pa.GetNeighbors()}
            except Exception:
                pass
            if _h_introduces_clash(pts[h], h, parent_idx, syms, pts,
                                   sibling_idxs=sib,
                                   threshold_factor=0.55,
                                   ignore_h_to_h=False):
                pts[h] = old
                if parent_idx in centre_h and h in centre_h[parent_idx]:
                    centre_h[parent_idx].remove(h)

        moved = sum(len(v) for v in centre_h.values())
        centres = sum(1 for v in centre_h.values() if v)

        report["h_atoms_moved"] = moved
        report["centres_rebuilt"] = centres
        report["violations_detected"] = centres
        report["max_angle_correction_deg"] = float(max_corr)
        if moved == 0:
            return xyz, report
        if not np.all(np.isfinite(pts)):
            # Never emit non-finite — fall back to input.
            return xyz, report
        return _format_xyz(orig_lines, syms, pts), report
    except Exception:
        return xyz, report


# ---------------------------------------------------------------------------
# Public detector
# ---------------------------------------------------------------------------


def detect_sp3_h_violations(xyz: str, mol,
                            tolerance_deg: float = 10.0) -> List[dict]:
    """Detect sp3 heavy atoms whose H-involving angles deviate from
    109.5° by more than ``tolerance_deg``.

    Returns a list of dicts, one per (atom_idx, h_idx, other_idx) triple
    that violates the tolerance, each with::

        {"atom_idx": int,
         "h_idx": int,
         "other_idx": int,
         "current_angle_deg": float,
         "ideal": 109.5}

    Implementation note: the angle of interest is (other-atom-H) with the
    sp3 atom as the vertex.  Both (H-A-H) and (X-A-H) are checked; the
    distinction is captured by ``other_idx`` (whether the other neighbour
    is itself an H or a heavy atom).
    """
    syms, pts, _ = _parse_xyz(xyz)
    if pts.shape[0] == 0 or mol is None:
        return []

    sp3_heavy = _sp3_heavy_atom_indices(mol)
    if not sp3_heavy:
        return []

    violations: List[dict] = []
    n_mol = mol.GetNumAtoms()
    for atom_idx in sp3_heavy:
        if atom_idx >= n_mol or atom_idx >= len(syms):
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        h_nbrs = [j for j in nbrs
                  if j < len(syms) and syms[j] == "H"]
        if not h_nbrs:
            continue
        for h_idx in h_nbrs:
            for other_idx in nbrs:
                if other_idx == h_idx or other_idx >= len(syms):
                    continue
                ang = _compute_angle_deg(pts[other_idx], pts[atom_idx],
                                          pts[h_idx])
                if math.isnan(ang):
                    continue
                if abs(ang - _IDEAL_TETRA_DEG) > tolerance_deg:
                    violations.append({
                        "atom_idx": atom_idx,
                        "h_idx": h_idx,
                        "other_idx": other_idx,
                        "current_angle_deg": ang,
                        "ideal": 109.5,
                    })
    return violations


# ---------------------------------------------------------------------------
# Fixer core: rotate one H so the (other-A-H) angle becomes tetrahedral
# ---------------------------------------------------------------------------


def _rotate_h_to_target_angle(pts: np.ndarray, h_idx: int, atom_idx: int,
                               other_idx: int, target_deg: float,
                               fallback_dir: Optional[np.ndarray] = None
                               ) -> Optional[np.ndarray]:
    """Compute a new H position so the (other-A-H) angle becomes
    ``target_deg``, by rotating the (A->H) vector toward/away from the
    (A->other) direction inside the plane spanned by those two vectors.

    Preserves the A-H distance exactly.  Returns ``None`` if the geometry
    is fully degenerate (zero vectors).  If A-H is collinear with
    A-other, ``fallback_dir`` is used as the in-plane perpendicular
    hint (e.g. the direction toward another neighbour of A); when not
    provided, an arbitrary orthogonal axis is constructed.
    """
    pivot = pts[atom_idx]
    v_h = pts[h_idx] - pivot
    v_other = pts[other_idx] - pivot
    d_ah = float(np.linalg.norm(v_h))
    d_ao = float(np.linalg.norm(v_other))
    if d_ah < 1e-6 or d_ao < 1e-6:
        return None
    u_other = v_other / d_ao
    # Decompose v_h into components parallel and perpendicular to u_other
    parallel = float(np.dot(v_h, u_other)) * u_other
    perp = v_h - parallel
    n_perp = float(np.linalg.norm(perp))
    if n_perp < 1e-6:
        # Degenerate: A-H is exactly collinear with A-other.  Use the
        # fallback direction (typically another A-neighbour) projected
        # to the plane perpendicular to u_other; if unavailable, build
        # an arbitrary orthonormal axis.
        if fallback_dir is not None:
            fd = np.asarray(fallback_dir, dtype=float)
            fd_perp = fd - float(np.dot(fd, u_other)) * u_other
            n_fd = float(np.linalg.norm(fd_perp))
            if n_fd > 1e-6:
                u_perp = fd_perp / n_fd
            else:
                u_perp = None  # fallback collinear too; build arbitrary
        else:
            u_perp = None
        if u_perp is None:
            # Build an arbitrary unit vector orthogonal to u_other
            seed = np.array([1.0, 0.0, 0.0]) \
                if abs(u_other[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            seed -= float(np.dot(seed, u_other)) * u_other
            u_perp = seed / float(np.linalg.norm(seed))
    else:
        u_perp = perp / n_perp
    # Build new v_h with given (other-A-H) angle and same A-H length
    theta = math.radians(target_deg)
    new_v_h = d_ah * (math.cos(theta) * u_other + math.sin(theta) * u_perp)
    return pivot + new_v_h


# ---------------------------------------------------------------------------
# Public fixer
# ---------------------------------------------------------------------------


def fix_sp3_h_tetrahedrality(xyz: str, mol, tolerance_deg: float = 10.0,
                              max_violations_per_atom: int = 4
                              ) -> Tuple[str, dict]:
    """Detect and repair sp3 H-tetrahedrality violations on a single XYZ.

    Topology guarantees:
        * Metal atoms are never moved.
        * Non-H atoms are never moved.
        * Only H positions are adjusted.
        * Each A-H bond length is preserved exactly.

    Algorithm
    ---------
    1. Detect violations via :func:`detect_sp3_h_violations`.
    2. Group violations by ``(atom_idx, h_idx)``; only the (other-A-H)
       triple with the largest deviation drives the move for that H.
    3. Cap moves per parent at ``max_violations_per_atom``.
    4. For each chosen (atom_idx, h_idx, other_idx), rotate the H within
       the plane (A-H, A-other) so the new (other-A-H) angle = 109.5°.
    5. Per-H rollback: if the rotation introduces a NEW non-bonded clash
       (heavy-H or H-H below 0.85·Σr_vdw), revert that H.

    Returns
    -------
    (new_xyz, report) — ``report`` keys:

        - ``"sp3_atoms_checked"``: int
        - ``"violations_detected"``: int
        - ``"h_atoms_moved"``: int
        - ``"max_angle_correction_deg"``: float (0.0 if no move)
        - ``"topology_preserved"``: bool (always True by construction)
    """
    report: dict = {
        "sp3_atoms_checked": 0,
        "violations_detected": 0,
        "h_atoms_moved": 0,
        "max_angle_correction_deg": 0.0,
        "topology_preserved": True,
    }
    if mol is None:
        return xyz, report
    try:
        syms, pts, orig_lines = _parse_xyz(xyz)
        if pts.shape[0] == 0:
            return xyz, report
        pts = pts.copy()

        sp3_heavy = _sp3_heavy_atom_indices(mol)
        report["sp3_atoms_checked"] = len(sp3_heavy)
        if not sp3_heavy:
            return xyz, report

        violations = detect_sp3_h_violations(xyz, mol,
                                              tolerance_deg=tolerance_deg)
        report["violations_detected"] = len(violations)
        if not violations:
            return xyz, report

        # For each (atom_idx, h_idx), keep only the violation with the
        # largest deviation — that drives the corrective rotation.
        worst_per_h: Dict[Tuple[int, int], dict] = {}
        for v in violations:
            key = (v["atom_idx"], v["h_idx"])
            dev = abs(v["current_angle_deg"] - _IDEAL_TETRA_DEG)
            if key not in worst_per_h:
                worst_per_h[key] = v
            else:
                cur = worst_per_h[key]
                cur_dev = abs(cur["current_angle_deg"] - _IDEAL_TETRA_DEG)
                if dev > cur_dev:
                    worst_per_h[key] = v

        # Apply per-parent cap: keep at most ``max_violations_per_atom``
        # H moves per sp3 atom, prioritising the largest deviations.
        by_parent: Dict[int, List[dict]] = {}
        for v in worst_per_h.values():
            by_parent.setdefault(v["atom_idx"], []).append(v)
        ordered_moves: List[dict] = []
        for parent, vs in by_parent.items():
            vs_sorted = sorted(
                vs,
                key=lambda x: -abs(x["current_angle_deg"]
                                    - _IDEAL_TETRA_DEG),
            )
            ordered_moves.extend(vs_sorted[:max_violations_per_atom])

        moved = 0
        max_corr = 0.0
        # Pre-compute neighbour map for sibling/fallback lookups
        nbr_map: Dict[int, List[int]] = {}
        for parent in by_parent.keys():
            try:
                a = mol.GetAtomWithIdx(parent)
                nbr_map[parent] = [n.GetIdx() for n in a.GetNeighbors()]
            except Exception:
                nbr_map[parent] = []
        for v in ordered_moves:
            atom_idx = v["atom_idx"]
            h_idx = v["h_idx"]
            other_idx = v["other_idx"]
            # Skip if pivot atom is a metal (defense in depth — sp3 RDKit
            # hybridisation should not flag metals, but guard).
            if atom_idx >= len(syms) or h_idx >= len(syms):
                continue
            if _is_metal_sym(syms[atom_idx]):
                continue
            if syms[h_idx] != "H":
                continue
            old_h = pts[h_idx].copy()
            # Build fallback in-plane direction from another A-neighbour
            # (anything that is neither the H being moved nor the pivot
            # nor the ``other`` reference).
            fb_dir: Optional[np.ndarray] = None
            for cand in nbr_map.get(atom_idx, []):
                if cand in (h_idx, other_idx) or cand >= len(syms):
                    continue
                vec = pts[cand] - pts[atom_idx]
                nv = float(np.linalg.norm(vec))
                if nv > 1e-6:
                    fb_dir = vec / nv
                    break
            new_h = _rotate_h_to_target_angle(
                pts, h_idx, atom_idx, other_idx, _IDEAL_TETRA_DEG,
                fallback_dir=fb_dir,
            )
            if new_h is None:
                continue
            # Bond-length sanity (rotation preserves length within FP
            # epsilon; guard anyway against arithmetic drift).
            d_old = float(np.linalg.norm(old_h - pts[atom_idx]))
            d_new = float(np.linalg.norm(new_h - pts[atom_idx]))
            if d_old > 1e-6 and abs(d_new - d_old) > 1e-3:
                continue
            # Speculatively place, check clash, rollback if needed.  Siblings
            # (other neighbours of the same parent — including geminal H)
            # are excluded so legitimate tetrahedral H-H pairs don't trip.
            siblings = set(nbr_map.get(atom_idx, []))
            pts[h_idx] = new_h
            if _h_introduces_clash(new_h, h_idx, atom_idx, syms, pts,
                                   sibling_idxs=siblings,
                                   threshold_factor=0.70):
                pts[h_idx] = old_h
                continue
            moved += 1
            corr = abs(v["current_angle_deg"] - _IDEAL_TETRA_DEG)
            if corr > max_corr:
                max_corr = corr

        report["h_atoms_moved"] = moved
        report["max_angle_correction_deg"] = float(max_corr)
        if moved == 0:
            return xyz, report
        return _format_xyz(orig_lines, syms, pts), report
    except Exception:
        # Fail-safe: never break a caller pipeline on this targeted fix.
        return xyz, report


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


def _selftest() -> None:
    """Three sanity checks (callable as ``python -m delfin._fix_sp3_h_tetrahedrality``):

    1. Synthetic methane with one H pushed to 180° from another H — fixer
       should restore tetrahedral geometry within tolerance.
    2. Real-ish CH3-OH (methanol) with one H bent — fixer should bend it
       back; oxygen stays put.
    3. Pristine tetrahedral methane — no-op (moves=0).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def _xyz_from_coords(syms_list, coords_arr):
        lines = [f"{len(syms_list)}", "selftest"]
        for s, c in zip(syms_list, coords_arr):
            lines.append(f"{s:4s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}")
        return "\n".join(lines) + "\n"

    # --- Case 1: methane, one H pushed to be collinear-ish with another ---
    mol_ch4 = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol_ch4, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_ch4)
    conf = mol_ch4.GetConformer()
    syms = [a.GetSymbol() for a in mol_ch4.GetAtoms()]
    coords = np.array([list(conf.GetAtomPosition(i))
                       for i in range(mol_ch4.GetNumAtoms())])
    # Identify C and its Hs
    c_idx = next(i for i, s in enumerate(syms) if s == "C")
    h_idxs = [i for i, s in enumerate(syms) if s == "H"]
    # Push one H so its angle to another H is ~170° (forced bad)
    c = coords[c_idx]
    h0, h1 = h_idxs[0], h_idxs[1]
    # Place h0 at opposite ray of h1 (≈180°)
    d_ah = float(np.linalg.norm(coords[h0] - c))
    dir_h1 = (coords[h1] - c)
    dir_h1 /= np.linalg.norm(dir_h1)
    coords[h0] = c - dir_h1 * d_ah  # angle h1-c-h0 ≈ 180°
    bad_xyz = _xyz_from_coords(syms, coords)
    new_xyz, rep = fix_sp3_h_tetrahedrality(bad_xyz, mol_ch4,
                                             tolerance_deg=10.0)
    print("CASE 1 (forced 180° H-C-H in CH4):")
    print(f"  violations={rep['violations_detected']} "
          f"moved={rep['h_atoms_moved']} "
          f"max_corr={rep['max_angle_correction_deg']:.2f}°")
    # Verify post-fix geometry
    syms2, pts2, _ = _parse_xyz(new_xyz)
    ang_after = _compute_angle_deg(pts2[h1], pts2[c_idx], pts2[h0])
    print(f"  H-C-H angle after fix: {ang_after:.2f}° (target 109.47°)")
    # Verify C-H bond lengths
    for h in h_idxs:
        d = float(np.linalg.norm(pts2[h] - pts2[c_idx]))
        print(f"  C-H[{h}] = {d:.4f} Å")
    print(f"  topology_preserved={rep['topology_preserved']}")

    # --- Case 2: methanol, force a bend, fixer corrects, O does not move ---
    mol_meoh = Chem.AddHs(Chem.MolFromSmiles("CO"))
    AllChem.EmbedMolecule(mol_meoh, randomSeed=7)
    AllChem.UFFOptimizeMolecule(mol_meoh)
    conf2 = mol_meoh.GetConformer()
    syms_m = [a.GetSymbol() for a in mol_meoh.GetAtoms()]
    coords_m = np.array([list(conf2.GetAtomPosition(i))
                         for i in range(mol_meoh.GetNumAtoms())])
    c_m = next(i for i, s in enumerate(syms_m) if s == "C")
    o_m = next(i for i, s in enumerate(syms_m) if s == "O")
    hs_on_c = [i for i, a in enumerate(mol_meoh.GetAtoms())
               if a.GetSymbol() == "H"
               and any(n.GetIdx() == c_m for n in a.GetNeighbors())]
    # Force one CH3 hydrogen to make a near-180° O-C-H angle
    h_bend = hs_on_c[0]
    d_ah_m = float(np.linalg.norm(coords_m[h_bend] - coords_m[c_m]))
    dir_oc = (coords_m[o_m] - coords_m[c_m])
    dir_oc /= np.linalg.norm(dir_oc)
    coords_m[h_bend] = coords_m[c_m] - dir_oc * d_ah_m  # ~180° O-C-H
    o_pos_pre = coords_m[o_m].copy()
    bad_xyz_m = _xyz_from_coords(syms_m, coords_m)
    new_xyz_m, rep_m = fix_sp3_h_tetrahedrality(bad_xyz_m, mol_meoh,
                                                  tolerance_deg=10.0)
    print("\nCASE 2 (CH3-OH with one H bent to 180° O-C-H):")
    print(f"  violations={rep_m['violations_detected']} "
          f"moved={rep_m['h_atoms_moved']} "
          f"max_corr={rep_m['max_angle_correction_deg']:.2f}°")
    syms3, pts3, _ = _parse_xyz(new_xyz_m)
    ang2 = _compute_angle_deg(pts3[o_m], pts3[c_m], pts3[h_bend])
    print(f"  O-C-H after fix: {ang2:.2f}° (target 109.47°)")
    d_o_drift = float(np.linalg.norm(pts3[o_m] - o_pos_pre))
    print(f"  O drift: {d_o_drift:.6f} Å (must be 0)")

    # --- Case 3: pristine methane — no-op ---
    mol_ch4b = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol_ch4b, randomSeed=1)
    AllChem.UFFOptimizeMolecule(mol_ch4b)
    conf3 = mol_ch4b.GetConformer()
    syms_p = [a.GetSymbol() for a in mol_ch4b.GetAtoms()]
    coords_p = np.array([list(conf3.GetAtomPosition(i))
                         for i in range(mol_ch4b.GetNumAtoms())])
    good_xyz = _xyz_from_coords(syms_p, coords_p)
    _, rep_g = fix_sp3_h_tetrahedrality(good_xyz, mol_ch4b,
                                          tolerance_deg=10.0)
    print("\nCASE 3 (pristine methane):")
    print(f"  violations={rep_g['violations_detected']} "
          f"moved={rep_g['h_atoms_moved']}")
    print("  (expect 0 moves)")


if __name__ == "__main__":
    _selftest()
