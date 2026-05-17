"""Welle-5b A — Universal VSEPR-correct H-placement post-pass.

End-of-pipeline pass that snaps every X-H bond direction to the
VSEPR-correct local geometry of its heavy parent.  Operates on the
finalised XYZ string (post UFF, post snap, post clash relief, post
B3 coord-angle correction, post B4 ring-H projection, post hapto
adjustment) — same wire-in position as Baustein 4.

Scope (heavy parent X with attached hydrogens):

    sp      X≡Y-H              linear            (1 H)
    sp²     X=CH-Y, =CH₂       trigonal planar   (1-2 H in plane)
    sp²     Ar-C-H, Ar-N-H     trigonal planar   (1 H in ring plane)
    sp²     X=NH, R-CO-OH      planar bent       (1 H in plane)
    sp³     X-CH₃              tetrahedral       (3 H umbrella, 109.5°)
    sp³     X-CH₂-Y            tetrahedral       (2 H, 109.5°)
    sp³     X-CH(-Y)-Z         tetrahedral       (1 H, 109.5°)
    sp³     X-NH₂              pyramidal         (2 H, ~107°)
    sp³     X-NH-Y             pyramidal         (1 H, ~107°)
    sp³     X-OH               bent              (1 H, ~109°)
    sp³     X-SH               bent              (1 H, ~92°)
    sp³     X-PH₂              pyramidal         (2 H, ~93°)
    sp³     X-BH₂              trigonal          (2 H, sp²-like 120°)
    sp³     X-SiH₃             tetrahedral       (3 H umbrella)

X-H ideal bond length table (Å):

    C-H 1.06 (sp), 1.08 (sp²), 1.09 (sp³)
    N-H 1.01, O-H 0.96, S-H 1.34, P-H 1.42
    B-H 1.19, Si-H 1.48, As-H 1.52, Se-H 1.46

Universal — element + bond-graph + bond-order + hybridization only.
NO SMILES literals, NO refcodes, NO named-ligand patterns.

Doctrine:
    - Default-OFF (env-flag opt-in).  Bit-exact when disabled.
    - Per-H rollback on (a) M-D invariant break, (b) new heavy clash.
    - Aromatic ring H is delegated to Baustein 4 / Iter-9 H1 — skipped
      here to avoid double-projection.
    - Metal-hydride (M-H) bonds are preserved as-is — bond-length
      reference table is in T3.1 (`DELFIN_MH_TABLE_FALLBACK`) and a
      VSEPR rotation around a metal centre is out of scope.
"""
from __future__ import annotations

import os
import math
import re
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Constants — kept consistent with delfin/_pi_h_projector.py.
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


# ---------------------------------------------------------------------------
# X-H ideal bond-length table — heavy-element × hybridization.
# ---------------------------------------------------------------------------

# Hybridization-keyed lookup: { sym: { 'sp': l, 'sp2': l, 'sp3': l } }
_IDEAL_XH: Dict[str, Dict[str, float]] = {
    "C":  {"sp": 1.06, "sp2": 1.08, "sp3": 1.09},
    "N":  {"sp": 1.00, "sp2": 1.00, "sp3": 1.01},
    "O":  {"sp": 0.96, "sp2": 0.96, "sp3": 0.96},
    "S":  {"sp": 1.34, "sp2": 1.34, "sp3": 1.34},
    "P":  {"sp": 1.42, "sp2": 1.42, "sp3": 1.42},
    "B":  {"sp": 1.19, "sp2": 1.19, "sp3": 1.19},
    "Si": {"sp": 1.48, "sp2": 1.48, "sp3": 1.48},
    "As": {"sp": 1.52, "sp2": 1.52, "sp3": 1.52},
    "Se": {"sp": 1.46, "sp2": 1.46, "sp3": 1.46},
}
_IDEAL_XH_FALLBACK: float = 1.09


# VSEPR rollback / detection tolerances.
_MD_INVARIANT_TOL: float = 0.05   # M-D break threshold (Å)
_MD_INTACT_FACTOR: float = 1.30   # Σr_cov scale for M-D snapshot
_CLASH_VDW_FACTOR: float = 0.80   # new-clash threshold

# H atoms whose current placement is within this angle (degrees) of the
# VSEPR-correct direction are left untouched.  Below this we re-place.
# Generous so we never re-emit micro-corrections that change the bit
# pattern when the input already looks fine.
_ANGLE_TOL_DEG: float = 5.0


# ---------------------------------------------------------------------------
# XYZ parsing (shape-preserving) — copy of the well-tested pattern from
# _pi_h_projector.py so a hand-written file with arbitrary header lines
# survives a round-trip with identical formatting on the non-touched lines.
# ---------------------------------------------------------------------------

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
            pts.append(np.array(
                [float(m.group(2)), float(m.group(3)), float(m.group(4))]
            ))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(
    orig_lines: List[str],
    syms: List[str],
    positions: np.ndarray,
) -> str:
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
# Bond-graph (geometric).  Identical scheme to _pi_h_projector.py.
# ---------------------------------------------------------------------------


def _build_geometric_adjacency(
    syms: List[str], pts: np.ndarray
) -> List[List[int]]:
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


# ---------------------------------------------------------------------------
# M-D invariant snapshot / check.
# ---------------------------------------------------------------------------


def _snapshot_md_bonds(
    pts: np.ndarray,
    syms: List[str],
) -> Dict[Tuple[int, int], float]:
    snap: Dict[Tuple[int, int], float] = {}
    n = len(syms)
    for i in range(n):
        if not _is_metal_sym(syms[i]):
            continue
        ri = _COV_RADII.get(syms[i], 1.5)
        for j in range(n):
            if j == i:
                continue
            sj = syms[j]
            if _is_metal_sym(sj):
                continue
            rj = _COV_RADII.get(sj, 1.5)
            d = float(np.linalg.norm(pts[i] - pts[j]))
            sigma = ri + rj
            if d <= _MD_INTACT_FACTOR * sigma:
                snap[(i, j)] = d
    return snap


def _md_invariant_violated(
    pre: Dict[Tuple[int, int], float],
    new_pts: np.ndarray,
    tol: float = _MD_INVARIANT_TOL,
) -> bool:
    for (m_idx, d_idx), pre_d in pre.items():
        new_d = float(np.linalg.norm(new_pts[m_idx] - new_pts[d_idx]))
        if abs(new_d - pre_d) > tol:
            return True
    return False


# ---------------------------------------------------------------------------
# Aromatic-ring detection — H on these atoms is owned by Baustein 4 and
# skipped here to avoid double-projection.
# ---------------------------------------------------------------------------


def _planar_ring_atom_set(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    ring_planar_tol: float = 0.10,
) -> set:
    """Return set of atom indices that belong to a planar 5/6-ring of
    C/N/O/S.  Conservative: any H on these atoms is considered owned by
    Baustein 4's projector and skipped here.
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

    planar_atoms: set = set()
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
        if ring_oop > ring_planar_tol:
            continue
        planar_atoms.update(ring)
    return planar_atoms


# ---------------------------------------------------------------------------
# Hybridization derivation — purely from the heavy bond-graph.
# ---------------------------------------------------------------------------


def _heavy_neighbors(
    idx: int, syms: List[str], nbrs: List[List[int]]
) -> List[int]:
    return [j for j in nbrs[idx]
            if syms[j] != "H" and not _is_metal_sym(syms[j])]


def _h_neighbors(
    idx: int, syms: List[str], nbrs: List[List[int]]
) -> List[int]:
    return [j for j in nbrs[idx] if syms[j] == "H"]


def _metal_neighbors(
    idx: int, syms: List[str], nbrs: List[List[int]]
) -> List[int]:
    return [j for j in nbrs[idx] if _is_metal_sym(syms[j])]


def _derive_hybridization(
    idx: int,
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    planar_ring_atoms: set,
) -> str:
    """Best-effort hybridization label for a heavy parent atom.

    Returns one of ``"sp"``, ``"sp2"``, ``"sp3"``.  Rules:
      - Member of a planar aromatic 5/6-ring → ``sp2``.
      - Total steric number 2 (one heavy nbr, no H, plus implicit lone-pair
        bias for C/N) → ``sp`` (alkyne / nitrile / cumulene C).
      - Heavy-neighbour count exactly 2 with one short bond (<1.30 Å for
        C-C, <1.25 Å for C-N) → ``sp2`` (vinyl, imine).
      - Heavy-neighbour count == 3 → ``sp2`` (planar by VSEPR; e.g.
        carbonyl C, =CH₂ terminal C with the implicit double bond not
        shown).
      - Else ``sp3``.

    sp / sp2 thresholds are conservative.  Bond-length triple-bond test
    uses pure C-C cov-radius sum scaled to 1.20 Å (C≡C ~1.20-1.22 Å).
    """
    sym = syms[idx]
    if idx in planar_ring_atoms:
        return "sp2"
    heavy_n = _heavy_neighbors(idx, syms, nbrs)
    n_heavy = len(heavy_n)

    # CN-2 sp1 detection (alkyne / nitrile / allene-end)
    if n_heavy == 1:
        j = heavy_n[0]
        d = float(np.linalg.norm(pts[idx] - pts[j]))
        # heavy partner has >=2 heavy nbrs and short distance to idx → triple/cumulene
        # C≡C ~1.20, C≡N ~1.16, C=O terminal ~1.21 → only triple-bond test
        # uses d < 1.27 Å for the C-X partner.
        if d < 1.27:
            # exclude carbonyl/imine: those have *two* heavy partners typically
            return "sp"
    # CN-2 with both heavy nbrs and a short bond to one of them → sp2 (vinyl, imine).
    if n_heavy == 2:
        for j in heavy_n:
            d = float(np.linalg.norm(pts[idx] - pts[j]))
            if (sym, syms[j]) in (("C", "C"), ("C", "N"), ("N", "C"), ("N", "N"),
                                  ("C", "O"), ("O", "C")):
                if d < 1.40:
                    return "sp2"
    # CN-3 heavy → trigonal planar (default sp2)
    if n_heavy == 3:
        return "sp2"
    return "sp3"


# ---------------------------------------------------------------------------
# Clash check.
# ---------------------------------------------------------------------------


def _h_introduces_clash(
    new_h_pos: np.ndarray,
    h_idx: int,
    parent_idx: int,
    skip_idxs: set,
    syms: List[str],
    pts: np.ndarray,
    threshold_factor: float = _CLASH_VDW_FACTOR,
    ignore_h_to_h: bool = False,
) -> bool:
    """Heavy clash check.  If ``ignore_h_to_h`` is True, H/H pairs are
    excluded — used during pass-1 to avoid blocking moves whose only
    "clash" is with another stale-position H atom that pass-2 will also
    relocate.  Pass-2 re-enables the check (the H positions are then
    final).
    """
    n = len(syms)
    rH_v = _VDW_RADII.get("H", 1.20)
    rH_c = _COV_RADII.get("H", 0.31)
    for j in range(n):
        if j == h_idx or j == parent_idx or j in skip_idxs:
            continue
        if ignore_h_to_h and syms[j] == "H":
            continue
        d = float(np.linalg.norm(new_h_pos - pts[j]))
        rj_v = _VDW_RADII.get(syms[j], 1.7 if syms[j] != "H" else 1.20)
        rj_c = _COV_RADII.get(syms[j], 1.5 if syms[j] != "H" else 0.31)
        sum_cov = rH_c + rj_c
        # Pre-bonded range: 0.85·Σr_cov < d < 1.10·Σr_cov is a normal bond.
        if 0.85 * sum_cov <= d <= 1.10 * sum_cov:
            continue
        if d < threshold_factor * (rH_v + rj_v):
            return True
    return False


# ---------------------------------------------------------------------------
# Direction generators per VSEPR template.
# ---------------------------------------------------------------------------


def _normalize(v: np.ndarray) -> Optional[np.ndarray]:
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        return None
    return v / n


def _methyl_umbrella_dirs(
    axis_in: np.ndarray, phase_rad: float, n_h: int = 3,
) -> List[np.ndarray]:
    """Return n_h methyl umbrella directions at sp3 angle 109.47° from
    ``-axis_in``, rotated around ``axis_in`` by ``phase_rad``.

    ``axis_in`` is the unit vector parent -> heavy nbr (single C-X bond).
    The umbrella opens AWAY from this axis (so H makes 109.47° with C-X).
    """
    a = axis_in
    u, v = _orthogonal_basis(a)
    cos_t = 1.0 / 3.0
    sin_t = math.sqrt(1.0 - cos_t * cos_t)
    dirs: List[np.ndarray] = []
    for k in range(3):
        ang = 2.0 * math.pi * k / 3.0 + phase_rad
        d = (-a) * cos_t + (math.cos(ang) * u + math.sin(ang) * v) * sin_t
        dn = _normalize(d)
        if dn is not None:
            dirs.append(dn)
    return dirs[:n_h]


def _orthogonal_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors perpendicular to ``axis`` (axis must be
    unit-length).  Deterministic — uses standard-basis trick.
    """
    if abs(axis[0]) < 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    else:
        helper = np.array([0.0, 1.0, 0.0])
    u = helper - np.dot(helper, axis) * axis
    un = float(np.linalg.norm(u))
    if un < 1e-9:
        u = np.array([0.0, 1.0, 0.0])
        un = 1.0
    u = u / un
    v = np.cross(axis, u)
    return u, v


def _vsepr_directions(
    parent_pos: np.ndarray,
    heavy_nbr_positions: List[np.ndarray],
    hybridization: str,
    n_h: int,
) -> List[np.ndarray]:
    """Return ``n_h`` unit vectors in the local VSEPR-correct directions
    away from ``parent_pos`` given the heavy-neighbour positions and the
    target hybridization.

    All vectors point outward from ``parent_pos``.  Caller multiplies by
    the X-H ideal bond length and adds ``parent_pos``.
    """
    n_heavy = len(heavy_nbr_positions)
    # Vectors parent -> heavy_nbr (unit).
    bonds_in: List[np.ndarray] = []
    for hp in heavy_nbr_positions:
        v = _normalize(np.asarray(hp) - parent_pos)
        if v is not None:
            bonds_in.append(v)

    # ---------------- sp (linear) --------------------------------------
    if hybridization == "sp":
        if not bonds_in:
            return [np.array([1.0, 0.0, 0.0])][:n_h]
        # H opposite of single heavy nbr
        return [-bonds_in[0]][:n_h]

    # ---------------- sp2 (trigonal planar 120°) -----------------------
    if hybridization == "sp2":
        if n_heavy == 0:
            # free planar centre — pick three vectors in xy-plane.
            base = []
            for k in range(3):
                ang = 2.0 * math.pi * k / 3.0
                base.append(np.array([math.cos(ang), math.sin(ang), 0.0]))
            return base[:n_h]
        if n_heavy == 1:
            a = bonds_in[0]
            # plane defined arbitrarily — pick orthogonal in-plane axis
            u, _ = _orthogonal_basis(a)
            # The two sp2 vectors at 120° from -a in the (a, u) plane.
            d1 = -math.cos(math.pi / 3.0) * a + math.sin(math.pi / 3.0) * u
            d2 = -math.cos(math.pi / 3.0) * a - math.sin(math.pi / 3.0) * u
            out = [_normalize(d1), _normalize(d2)]
            return [d for d in out if d is not None][:n_h]
        if n_heavy >= 2:
            # Plane defined by parent + 2 heavy nbrs.  Third direction is
            # the in-plane vector that bisects the gap between bonds_in.
            a, b = bonds_in[0], bonds_in[1]
            bis = a + b
            d = _normalize(-bis)
            if d is None:
                # collinear → fall back to perpendicular
                d = _normalize(np.cross(a, np.array([1.0, 0.0, 0.0])))
                if d is None:
                    d = np.array([0.0, 0.0, 1.0])
            return [d][:n_h]
        return []

    # ---------------- sp3 (tetrahedral, 109.5°) ------------------------
    if n_heavy == 0:
        # free centre — pick the canonical tetrahedron.
        tet = [
            np.array([+1.0, +1.0, +1.0]),
            np.array([+1.0, -1.0, -1.0]),
            np.array([-1.0, +1.0, -1.0]),
            np.array([-1.0, -1.0, +1.0]),
        ]
        return [v / math.sqrt(3.0) for v in tet][:n_h]
    if n_heavy == 1:
        # Methyl umbrella around the parent->heavy axis (phase 0).
        return _methyl_umbrella_dirs(bonds_in[0], phase_rad=0.0, n_h=n_h)
    if n_heavy == 2:
        # Methylene: H atoms straddle the plane defined by the two heavy
        # bonds.  Bisector = (a + b) normalised.  In-plane "out" is
        # -bisector.  Out-of-plane axis = cross(a, b).  H₁,H₂ at ±cos(α)
        # along out-of-plane, with α/2 = (180 - 109.5)/2 = 35.26°.
        a, b = bonds_in[0], bonds_in[1]
        bis_raw = a + b
        bis_n = _normalize(bis_raw)
        if bis_n is None:
            # collinear heavy — fall back to one perpendicular
            u, _ = _orthogonal_basis(a)
            bis_n = u
        out_axis = _normalize(np.cross(a, b))
        if out_axis is None:
            # collinear — pick arbitrary perpendicular
            u, v = _orthogonal_basis(bis_n)
            out_axis = v
        # 109.5° H-X-H → half angle = 54.75° from the bisector plane,
        # i.e. out-of-plane component = sin(54.75°) ≈ 0.8165.
        half = math.radians(109.47 / 2.0)
        cos_h = math.cos(half)  # in-plane (toward -bisector)
        sin_h = math.sin(half)  # out-of-plane
        d1 = -bis_n * cos_h + out_axis * sin_h
        d2 = -bis_n * cos_h - out_axis * sin_h
        return [_normalize(d1), _normalize(d2)][:n_h]
    if n_heavy >= 3:
        # Methine — H opposite the (a + b + c) vector.
        s = sum(bonds_in[:3])
        d = _normalize(-s)
        if d is None:
            d = bonds_in[0] * -1.0
        return [d][:n_h]
    return []


# ---------------------------------------------------------------------------
# Main entry.
# ---------------------------------------------------------------------------


def _angle_between(u: np.ndarray, v: np.ndarray) -> float:
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu < 1e-9 or nv < 1e-9:
        return 0.0
    c = float(np.dot(u, v) / (nu * nv))
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def _greedy_assign(
    h_positions: List[np.ndarray],
    target_dirs: List[np.ndarray],
    parent_pos: np.ndarray,
) -> List[int]:
    """Greedy matching: assign each existing H to the closest target dir
    (lowest angle).  Returns a list ``mapping`` such that
    ``target_dirs[mapping[i]]`` is the chosen target for existing H ``i``.

    Greedy minimum-angle matching is deterministic for any inputs and is
    cheap (≤ 3 H per typical parent).
    """
    n_t = len(target_dirs)
    used: List[bool] = [False] * n_t
    mapping: List[int] = []
    for hp in h_positions:
        v = hp - parent_pos
        best_k = -1
        best_ang = 1e9
        for k in range(n_t):
            if used[k]:
                continue
            ang = _angle_between(v, target_dirs[k])
            if ang < best_ang:
                best_ang = ang
                best_k = k
        if best_k < 0:
            mapping.append(0 if n_t > 0 else -1)
        else:
            mapping.append(best_k)
            used[best_k] = True
    return mapping


def _collect_parent_plans(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    planar_atoms: set,
) -> List[Tuple[int, List[int], List[int], List[np.ndarray], float]]:
    """Return per-parent snap plans: list of tuples
    ``(parent_idx, h_idxs, heavy_idxs, target_dirs, d_xh)``.

    A plan is only emitted when:
      - parent is heavy non-metal not in a planar aromatic ring
      - parent has ≥ 1 H neighbour, none bonded to a metal
      - heavy + H count ≤ 4 for C/N/O/F (defensive against perception
        artifacts)
      - a clean VSEPR template (≥ n_h directions) is available.
    """
    plans: List[Tuple[int, List[int], List[int], List[np.ndarray], float]] = []
    n = len(syms)
    for parent_idx in range(n):
        sym = syms[parent_idx]
        if sym == "H" or _is_metal_sym(sym):
            continue
        if parent_idx in planar_atoms:
            continue
        h_idxs = _h_neighbors(parent_idx, syms, nbrs)
        if not h_idxs:
            continue
        if any(_metal_neighbors(h, syms, nbrs) for h in h_idxs):
            continue
        heavy_n = _heavy_neighbors(parent_idx, syms, nbrs)
        if len(heavy_n) + len(h_idxs) > 4 and sym in ("C", "N", "O", "F"):
            continue
        hybridization = _derive_hybridization(
            parent_idx, syms, pts, nbrs, planar_atoms,
        )
        target_dirs = _vsepr_directions(
            pts[parent_idx],
            [pts[j] for j in heavy_n],
            hybridization,
            len(h_idxs),
        )
        if len(target_dirs) < len(h_idxs):
            continue
        d_xh = _IDEAL_XH.get(sym, {}).get(hybridization, _IDEAL_XH_FALLBACK)
        plans.append((parent_idx, h_idxs, heavy_n, target_dirs, d_xh))
    return plans


def _snap_xyz(
    syms: List[str],
    pts: np.ndarray,
) -> int:
    """Mutate ``pts`` in-place: snap every H bonded to a non-aromatic-ring
    heavy parent to the VSEPR direction implied by the heavy bond-graph
    and target hybridization.  Returns the count of moved H atoms.

    Two-pass strategy:
      pass-1: snap all H to VSEPR-correct directions, ignoring H-H
              clashes (input H atoms may also be misplaced — a "clash"
              against a stale H is not informative).  Heavy-atom clashes
              and M-D invariant breaks ARE checked per-H.
      pass-2: re-evaluate after all H are at their target positions.
              For each H, if a heavy clash now exists (or the move did
              NOT happen in pass-1 due to heavy clash), the H stays at
              its current position.
    """
    n = len(syms)
    if n == 0:
        return 0
    nbrs = _build_geometric_adjacency(syms, pts)
    planar_atoms = _planar_ring_atom_set(syms, pts, nbrs)
    pre_md = _snapshot_md_bonds(pts, syms)

    plans = _collect_parent_plans(syms, pts, nbrs, planar_atoms)
    if not plans:
        return 0

    # Pass-1 — per-parent placement:
    #   * sp3 methyl (CN-1 + 3H): rigid umbrella with phase optimised to
    #     minimise heavy clashes (12 candidate phases spanning 0..2π/3).
    #   * sp3 methylene (CN-2 + 2H): both H placed at the two computed
    #     directions; per-H heavy-clash check.
    #   * all other VSEPR cases: greedy H→target assignment.
    # skip-set per H: parent ∪ heavy_n ∪ H-siblings.  H-H clashes
    # ignored (stale-H artifacts).  Heavy clashes blocked.  M-D invariant
    # blocked with per-H rollback.
    moved = 0
    for parent_idx, h_idxs, heavy_n, target_dirs, d_xh in plans:
        skip_set = {parent_idx} | set(h_idxs) | set(heavy_n)
        sym = syms[parent_idx]
        n_h = len(h_idxs)

        # Methyl umbrella special case: rotate as rigid unit around C-X axis.
        is_methyl = (n_h == 3 and len(heavy_n) == 1
                     and sym in ("C", "Si", "N", "P", "B"))
        if is_methyl:
            a = _normalize(pts[heavy_n[0]] - pts[parent_idx])
            if a is None:
                continue
            best_phase = 0.0
            best_clash = 1 << 30
            best_dev = 1e9
            # 12 candidate phases over one third-symmetric period.
            n_probes = 12
            for k in range(n_probes):
                phase = (2.0 * math.pi / 3.0) * (k / n_probes)
                dirs = _methyl_umbrella_dirs(a, phase_rad=phase, n_h=3)
                if len(dirs) < 3:
                    continue
                # Greedy: count heavy clashes for this phase
                hp = [pts[parent_idx] + d_xh * d for d in dirs]
                clash_count = 0
                for h_idx, new_pos in zip(h_idxs, hp):
                    if _h_introduces_clash(
                        new_pos, h_idx, parent_idx, skip_set, syms, pts,
                        ignore_h_to_h=True,
                    ):
                        clash_count += 1
                # Tiebreak: prefer phase closest to existing H positions
                # (minimise rotation), so default-OFF behaviour is bit-exact
                # for already-good methyls when env-flag is OFF.
                mapping_tmp = _greedy_assign(
                    [pts[h] for h in h_idxs], dirs, pts[parent_idx],
                )
                dev = 0.0
                for li, h_idx in enumerate(h_idxs):
                    v = pts[h_idx] - pts[parent_idx]
                    dev += _angle_between(v, dirs[mapping_tmp[li]])
                if (clash_count, dev) < (best_clash, best_dev):
                    best_clash = clash_count
                    best_dev = dev
                    best_phase = phase
            best_dirs = _methyl_umbrella_dirs(a, phase_rad=best_phase, n_h=3)
            if len(best_dirs) < 3:
                continue
            mapping = _greedy_assign(
                [pts[h] for h in h_idxs], best_dirs, pts[parent_idx],
            )
            # Skip if every existing H is already within tolerance of its
            # assigned target — preserves bit-exactness for happy inputs.
            already_correct = True
            for li, h_idx in enumerate(h_idxs):
                v = pts[h_idx] - pts[parent_idx]
                if _angle_between(v, best_dirs[mapping[li]]) >= _ANGLE_TOL_DEG:
                    already_correct = False
                    break
            if already_correct:
                continue
            # Apply 3 H moves; rollback if M-D invariant breaks.
            old_h_pos = {h: pts[h].copy() for h in h_idxs}
            for li, h_idx in enumerate(h_idxs):
                tdir = best_dirs[mapping[li]]
                if _h_introduces_clash(
                    pts[parent_idx] + d_xh * tdir, h_idx, parent_idx,
                    skip_set, syms, pts, ignore_h_to_h=True,
                ):
                    # Skip this individual H; keep others in umbrella
                    continue
                pts[h_idx] = pts[parent_idx] + d_xh * tdir
                moved += 1
            if _md_invariant_violated(pre_md, pts, tol=_MD_INVARIANT_TOL):
                for h, p in old_h_pos.items():
                    pts[h] = p
                # subtract from moved count
                moved -= n_h
            continue

        # Generic per-H mapping (sp2 / sp3 CN>=2 / pyramidal / bent / ...).
        mapping = _greedy_assign(
            [pts[h] for h in h_idxs], target_dirs, pts[parent_idx],
        )
        for local_i, h_idx in enumerate(h_idxs):
            tgt_k = mapping[local_i]
            if tgt_k < 0 or tgt_k >= len(target_dirs):
                continue
            tdir = target_dirs[tgt_k]
            if tdir is None:
                continue
            cur_vec = pts[h_idx] - pts[parent_idx]
            cur_ang = _angle_between(cur_vec, tdir)
            if cur_ang < _ANGLE_TOL_DEG:
                continue
            new_pos = pts[parent_idx] + d_xh * tdir
            if _h_introduces_clash(
                new_pos, h_idx, parent_idx, skip_set, syms, pts,
                ignore_h_to_h=True,
            ):
                continue
            old_pos = pts[h_idx].copy()
            pts[h_idx] = new_pos
            if _md_invariant_violated(pre_md, pts, tol=_MD_INVARIANT_TOL):
                pts[h_idx] = old_pos
                continue
            moved += 1
    return moved


# ---------------------------------------------------------------------------
# Welle-5f D — Inter-substituent rotamer rescue (PMe3 / NMe3 / tBu pattern).
#
# Pattern: heavy sp3 atom X with >=2 alkyl substituents Ck where each Ck has
# >=2 hydrogens (e.g. PMe3, NMe3, CMe4, SiMe3).  After pass-1 places each
# methyl umbrella in isolation (heavy-clash-only optimisation), the three
# CH3 groups can still be in an "eclipsed" rotamer relative to each other
# at the shared parent.  UFF doesn't always reach the staggered minimum.
#
# Strategy:
#   1. Identify shared parents (sp3 non-metal heavy atom with >=2 alkyl
#      neighbours, each having >=2 H).
#   2. For each substituent in deterministic order, rotate it rigidly
#      around the parent-substituent bond axis to minimise the H-H clash
#      score against H atoms on the OTHER substituents of the SAME shared
#      parent.
#   3. Sweep until no substituent changes phase (max N_SWEEPS).
#
# Default-OFF — opt-in via ``DELFIN_5F_D_ALKYL_ROTAMER=1``.  Bit-exact
# when disabled (early-return before any geometry touch).  Per-substituent
# rollback on M-D invariant break or net heavy-clash worsening.
# ---------------------------------------------------------------------------

_ROTAMER_N_PROBES: int = 12
_ROTAMER_MAX_SWEEPS: int = 3
_ROTAMER_HH_CLASH_THR: float = 2.4  # A - count H-H pairs below this
_ROTAMER_HH_WEIGHT_THR: float = 2.0  # A - pairs below this weighted 4x


def _rotate_points_around_axis(
    points: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float,
) -> np.ndarray:
    """Rotate ``points`` (N x 3) around ``axis`` (unit) through ``origin``
    by ``angle`` radians.  Rodrigues formula.  Deterministic numpy.
    """
    a = axis
    c = math.cos(angle)
    s = math.sin(angle)
    one_c = 1.0 - c
    K = np.array([
        [c + a[0] * a[0] * one_c,
         a[0] * a[1] * one_c - a[2] * s,
         a[0] * a[2] * one_c + a[1] * s],
        [a[1] * a[0] * one_c + a[2] * s,
         c + a[1] * a[1] * one_c,
         a[1] * a[2] * one_c - a[0] * s],
        [a[2] * a[0] * one_c - a[1] * s,
         a[2] * a[1] * one_c + a[0] * s,
         c + a[2] * a[2] * one_c],
    ])
    shifted = points - origin
    rotated = shifted @ K.T
    return rotated + origin


def _substituent_h_indices(
    sub_c: int, parent: int, syms: List[str], nbrs: List[List[int]],
) -> List[int]:
    """Return H atom indices attached directly to ``sub_c`` (one bond),
    not crossing back through ``parent``.  For a CH3 group this is the
    3 methyl hydrogens.  Conservative: first shell only.
    """
    return [j for j in nbrs[sub_c]
            if syms[j] == "H" and j != parent]


def _hh_clash_score_between(
    moving_h: List[int],
    fixed_h: List[int],
    pts: np.ndarray,
) -> float:
    """Sum-of-clash score: 4x weight for pairs <2.0 A, 1x weight for
    pairs in [2.0, 2.4) A.  Lower is better.  Deterministic.
    """
    score = 0.0
    for hi in moving_h:
        for hj in fixed_h:
            d = float(np.linalg.norm(pts[hi] - pts[hj]))
            if d < _ROTAMER_HH_WEIGHT_THR:
                score += 4.0 * (_ROTAMER_HH_WEIGHT_THR - d)
            elif d < _ROTAMER_HH_CLASH_THR:
                score += (_ROTAMER_HH_CLASH_THR - d)
    return score


def _heavy_clash_worsens(
    moving_idxs: List[int],
    pts_new: np.ndarray,
    pts_old: np.ndarray,
    syms: List[str],
    skip: set,
    threshold_factor: float = _CLASH_VDW_FACTOR,
) -> bool:
    """Return True if rotating ``moving_idxs`` introduces NEW heavy-atom
    clashes (H to non-H) that didn't exist before.  Pure-H clashes are
    not counted here (handled separately by the H-H scorer).
    """
    n = len(syms)
    rH_v = _VDW_RADII.get("H", 1.20)
    for hi in moving_idxs:
        if syms[hi] != "H":
            continue
        for j in range(n):
            if j == hi or j in skip or syms[j] == "H":
                continue
            rj_v = _VDW_RADII.get(syms[j], 1.70)
            rj_c = _COV_RADII.get(syms[j], 1.50)
            sum_cov = _COV_RADII.get("H", 0.31) + rj_c
            d_new = float(np.linalg.norm(pts_new[hi] - pts_new[j]))
            if 0.85 * sum_cov <= d_new <= 1.10 * sum_cov:
                continue
            if d_new < threshold_factor * (rH_v + rj_v):
                d_old = float(np.linalg.norm(pts_old[hi] - pts_old[j]))
                if d_new < d_old - 1e-6:
                    return True
    return False


def _collect_shared_parents(
    syms: List[str], pts: np.ndarray, nbrs: List[List[int]],
    planar_atoms: set,
) -> List[Tuple[int, List[int]]]:
    """Return list of (parent_idx, [sub_c_idxs]) for every sp3 non-metal
    heavy atom that has >=2 alkyl substituents (each with >=2 H atoms
    directly attached).  Deterministic iteration order (ascending idx).
    """
    out: List[Tuple[int, List[int]]] = []
    for parent in range(len(syms)):
        sym = syms[parent]
        if sym == "H" or _is_metal_sym(sym):
            continue
        if parent in planar_atoms:
            continue
        if sym not in ("P", "N", "C", "B", "Si", "As", "Sb", "S", "O"):
            continue
        heavy_n = _heavy_neighbors(parent, syms, nbrs)
        alkyl_subs: List[int] = []
        for c in heavy_n:
            if syms[c] not in ("C", "N", "Si"):
                continue
            if c in planar_atoms:
                continue
            h_on_c = _h_neighbors(c, syms, nbrs)
            if len(h_on_c) < 2:
                continue
            alkyl_subs.append(c)
        if len(alkyl_subs) >= 2:
            out.append((parent, sorted(alkyl_subs)))
    return out


def _snap_inter_substituent_rotamers(
    syms: List[str], pts: np.ndarray,
) -> int:
    """Rotate each alkyl substituent rigidly around its parent-axis to
    minimise inter-substituent H-H clashes.  Mutates ``pts`` in-place.
    Returns number of substituents rotated.  Deterministic.
    """
    n = len(syms)
    if n == 0:
        return 0
    nbrs = _build_geometric_adjacency(syms, pts)
    planar_atoms = _planar_ring_atom_set(syms, pts, nbrs)
    pre_md = _snapshot_md_bonds(pts, syms)
    shared = _collect_shared_parents(syms, pts, nbrs, planar_atoms)
    if not shared:
        return 0

    n_rotated = 0
    for parent, sub_cs in shared:
        sub_h_map: Dict[int, List[int]] = {}
        for c in sub_cs:
            sub_h_map[c] = _substituent_h_indices(c, parent, syms, nbrs)
        if not any(sub_h_map.values()):
            continue

        for sweep in range(_ROTAMER_MAX_SWEEPS):
            any_change = False
            for c in sub_cs:
                moving_h = sub_h_map[c]
                if not moving_h:
                    continue
                axis = _normalize(pts[c] - pts[parent])
                if axis is None:
                    continue
                fixed_h: List[int] = []
                for c2 in sub_cs:
                    if c2 == c:
                        continue
                    fixed_h.extend(sub_h_map[c2])
                if not fixed_h:
                    continue
                skip = {parent, c} | set(sub_cs)
                for hl in sub_h_map.values():
                    skip.update(hl)
                old_pts = pts.copy()
                base_score = _hh_clash_score_between(moving_h, fixed_h, pts)
                best_angle = 0.0
                best_score = base_score
                period = (2.0 * math.pi / 3.0) if len(moving_h) == 3 else 2.0 * math.pi
                for k in range(1, _ROTAMER_N_PROBES):
                    angle = period * (k / _ROTAMER_N_PROBES)
                    moving_pts_arr = np.array([old_pts[h] for h in moving_h])
                    rotated = _rotate_points_around_axis(
                        moving_pts_arr, pts[parent], axis, angle,
                    )
                    trial_pts = pts.copy()
                    for li, h_idx in enumerate(moving_h):
                        trial_pts[h_idx] = rotated[li]
                    if _heavy_clash_worsens(
                        moving_h, trial_pts, old_pts, syms, skip,
                    ):
                        continue
                    if _md_invariant_violated(pre_md, trial_pts):
                        continue
                    score = _hh_clash_score_between(moving_h, fixed_h, trial_pts)
                    if score < best_score - 1e-9:
                        best_score = score
                        best_angle = angle
                if best_angle != 0.0:
                    moving_pts_arr = np.array([pts[h] for h in moving_h])
                    rotated = _rotate_points_around_axis(
                        moving_pts_arr, pts[parent], axis, best_angle,
                    )
                    for li, h_idx in enumerate(moving_h):
                        pts[h_idx] = rotated[li]
                    n_rotated += 1
                    any_change = True
            if not any_change:
                break
    return n_rotated


def correct_xyz(xyz_str: str) -> str:
    """Apply VSEPR H snap to a single XYZ.  Failsafe: any exception
    returns the input unchanged.

    Welle-5f D: after the per-parent VSEPR snap, an optional second pass
    rotates each alkyl substituent rigidly around its parent-axis to
    minimise inter-substituent H-H clashes (PMe3 / NMe3 / tBu pattern).
    Opt-in via ``DELFIN_5F_D_ALKYL_ROTAMER=1`` (default OFF, bit-exact
    when disabled).
    """
    try:
        syms, pts, orig_lines = _parse_xyz(xyz_str)
        if pts.shape[0] == 0:
            return xyz_str
        pts = pts.copy()
        moved = _snap_xyz(syms, pts)
        rot_moved = 0
        try:
            if os.environ.get("DELFIN_5F_D_ALKYL_ROTAMER", "0") not in ("", "0"):
                rot_moved = _snap_inter_substituent_rotamers(syms, pts)
        except Exception:
            rot_moved = 0
        if moved == 0 and rot_moved == 0:
            return xyz_str
        return _format_xyz(orig_lines, syms, pts)
    except Exception:
        return xyz_str


def snap_h_to_vsepr(
    xyz_str: str,
    mol=None,
    *,
    mode: str = "end_of_pipeline",
) -> str:
    """Public entry: VSEPR-correct snap of every X-H bond direction.

    ``mol`` is currently unused — hybridization is derived from the
    bond-graph geometry, not from RDKit.  The parameter is accepted so
    the dispatch helper in ``smiles_converter.py`` can pass the parent
    molecule in case a future revision switches to RDKit-derived
    hybridization (``GetHybridization()``).

    ``mode`` is reserved for future use ("pre_uff", "post_b4", etc).
    """
    return correct_xyz(xyz_str)


def correct_results(mol, results):
    """Apply VSEPR snap to each (xyz, label) entry.  Failsafe per-entry.

    Signature mirrors ``_pi_h_projector.correct_results`` so the
    dispatch helper in ``smiles_converter.py`` can call it with the
    same arguments.
    """
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = snap_h_to_vsepr(xyz, mol)
            out.append(
                (new_xyz, lbl) if len(entry) == 2
                else (new_xyz,) + tuple(entry[1:])
            )
        except Exception:
            out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Welle-5f-E: Aliphatic alkyl-donor → metal axis rescue.
#
# D-XALKEA Zr CN3 pattern: ETKDG/UFF leaves the C atom of an M-CH2-X
# alkyl donor at d_MC ~ 2.40 A while the real sigma-bond is ~2.20 A.
# UFF cannot rotate the CH2 substituent into the sp3-tetrahedral axis
# that points at the metal.  This post-pass snaps such broken-but-near
# alkyl donors back to the ideal M-C bond length along the existing
# M-C direction, rigidly carrying the entire substituent fragment.
#
# Universal: element + bond-graph only.  No SMILES literals / refcodes.
# Default-OFF env-flag.  Bit-exact when disabled.  Per-violation
# rollback on M-D invariant break or new heavy clash.
# ---------------------------------------------------------------------------


# Per-(metal,C) contraction factor on the covalent-radius sum that gives
# an empirically reasonable M-C(sp3) sigma-bond length.  Calibrated
# against CSD-typical ranges: Zr-C ~2.20 A (cov sum 2.51 -> 0.877),
# Ti-C ~2.10 (sum 2.36 -> 0.890), Hf-C ~2.22 (sum 2.51 -> 0.885),
# Mo-C ~2.20 (sum 2.30 -> 0.957), Fe-C ~2.05 (sum 2.08 -> 0.986).
# A single factor of 0.88 gives <0.05 A error for early-TM cases that
# trigger this pattern; later/3d metals are within ~0.10 A.  We do not
# need exact bond length here -- only a sane scalar in the
# stretched-vs-correct discriminator band 0.85..1.15.
_ALKYL_DONOR_MC_FACTOR: float = 0.88

# Trigger band: only act on candidates whose current d_MC is in
# (_ALKYL_DONOR_LO * ideal, _ALKYL_DONOR_HI * ideal).  Below LO =
# already bonded, leave alone (avoid bit-flip on good frames).  Above
# HI = not a donor (truly non-bonded), leave to upstream isomer
# enumerator.
_ALKYL_DONOR_LO: float = 1.15
_ALKYL_DONOR_HI: float = 1.45


def _ideal_mc_length(metal_sym: str, c_sym: str = "C") -> Optional[float]:
    rm = _COV_RADII.get(metal_sym)
    rc = _COV_RADII.get(c_sym)
    if rm is None or rc is None:
        return None
    return _ALKYL_DONOR_MC_FACTOR * (rm + rc)


def _collect_alkyl_donor_candidates(
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    planar_ring_atoms: set,
) -> List[Tuple[int, int, float]]:
    """Return list of (c_idx, m_idx, ideal_mc) tuples.

    Criteria (universal, element + bond-graph only):
      * C atom (NOT in any planar aromatic ring).
      * Exactly ONE heavy non-metal neighbour X (M-CH2-X / M-CHR-X /
        M-CR2-X tail of an alkyl donor).  Aromatic / vinyl carbons are
        excluded (planar_ring_atoms gate + heavy-neighbour count gate).
      * Geometrically close-but-stretched to a metal M (single M-C
        pair in the band _ALKYL_DONOR_LO * ideal < d < _ALKYL_DONOR_HI
        * ideal).  Picks the SINGLE closest metal in that band.
    """
    candidates: List[Tuple[int, int, float]] = []
    n = len(syms)
    for c_idx in range(n):
        if syms[c_idx] != "C":
            continue
        if c_idx in planar_ring_atoms:
            continue
        heavy_n = _heavy_neighbors(c_idx, syms, nbrs)
        if len(heavy_n) != 1:
            continue
        best_m = -1
        best_d = 1e9
        best_ideal = 0.0
        for m_idx in range(n):
            if m_idx == c_idx:
                continue
            if not _is_metal_sym(syms[m_idx]):
                continue
            ideal = _ideal_mc_length(syms[m_idx], "C")
            if ideal is None:
                continue
            d = float(np.linalg.norm(pts[c_idx] - pts[m_idx]))
            if _ALKYL_DONOR_LO * ideal < d < _ALKYL_DONOR_HI * ideal:
                if d < best_d:
                    best_d = d
                    best_m = m_idx
                    best_ideal = ideal
        if best_m >= 0:
            candidates.append((c_idx, best_m, best_ideal))
    return candidates


def _collect_rigid_fragment(
    seed_idx: int,
    blocked: int,
    syms: List[str],
    nbrs: List[List[int]],
) -> List[int]:
    """BFS from ``seed_idx`` over non-metal atoms; never cross ``blocked``.

    Returns the substituent atom indices to translate as a rigid unit
    (seed + H neighbours + the entire X chain attached via the heavy
    neighbour, excluding the metal direction).
    """
    visited = {seed_idx}
    stack = [seed_idx]
    out: List[int] = [seed_idx]
    while stack:
        cur = stack.pop()
        for nx in nbrs[cur]:
            if nx in visited or nx == blocked:
                continue
            if _is_metal_sym(syms[nx]):
                # never carry a metal as a fragment atom
                continue
            visited.add(nx)
            out.append(nx)
            stack.append(nx)
    return out


def _snap_alkyl_axis_xyz(
    syms: List[str],
    pts: np.ndarray,
) -> int:
    """In-place snap of stretched M-CHx-X candidates to ideal M-C length.

    Per-candidate strategy (translate-only, no rotation in v1):
      1. Compute fragment = {C, H-neighbours, X-side substituent tree}.
      2. Translate the entire fragment along the (C->M) direction by
         delta = (d_now - ideal).  Positive delta = pull toward metal.
      3. Snapshot M-D bonds AFTER move; if any pre-existing intact M-D
         bond (other than the new one we just made) broke beyond
         _MD_INVARIANT_TOL, ROLLBACK the fragment translation.
      4. New-clash check: if any heavy-heavy distance in the moved
         fragment is now below _CLASH_VDW_FACTOR * (rA + rB) against
         a non-fragment heavy atom (excluding the metal M, since the
         new M-C bond is the intended outcome), ROLLBACK.
    Returns count of candidates that were actually snapped.
    """
    n = len(syms)
    if n == 0:
        return 0
    nbrs = _build_geometric_adjacency(syms, pts)
    planar_atoms = _planar_ring_atom_set(syms, pts, nbrs)
    cands = _collect_alkyl_donor_candidates(syms, pts, nbrs, planar_atoms)
    if not cands:
        return 0
    pre_md = _snapshot_md_bonds(pts, syms)
    moved = 0
    for c_idx, m_idx, ideal in cands:
        v_cm = pts[m_idx] - pts[c_idx]
        d_now = float(np.linalg.norm(v_cm))
        if d_now < 1e-6:
            continue
        u_cm = v_cm / d_now
        delta = d_now - ideal  # positive -> shrink the M-C distance
        if abs(delta) < 1e-3:
            continue
        frag = _collect_rigid_fragment(c_idx, m_idx, syms, nbrs)
        old_pos = {i: pts[i].copy() for i in frag}
        shift = delta * u_cm
        for i in frag:
            pts[i] = pts[i] + shift
        # Drop the (m_idx, c_idx) pair from the M-D snapshot before the
        # check: changing this bond is the intended outcome, not a
        # violation.  Check every OTHER intact M-D bond for breakage.
        pre_md_other = {
            k: v for k, v in pre_md.items() if k != (m_idx, c_idx)
        }
        if _md_invariant_violated(
            pre_md_other, pts, tol=_MD_INVARIANT_TOL
        ):
            for i, p in old_pos.items():
                pts[i] = p
            continue
        # New heavy-clash check.  Skip H atoms (handled by H-VSEPR
        # pass).  Skip the metal we just bonded to.  Skip self.
        frag_set = set(frag)
        clash = False
        for i in frag:
            if syms[i] == "H":
                continue
            ri_v = _VDW_RADII.get(syms[i], 1.7)
            for j in range(n):
                if j in frag_set:
                    continue
                if j == m_idx:
                    continue
                if syms[j] == "H":
                    continue
                rj_v = _VDW_RADII.get(syms[j], 1.7)
                d = float(np.linalg.norm(pts[i] - pts[j]))
                if d < _CLASH_VDW_FACTOR * (ri_v + rj_v):
                    clash = True
                    break
            if clash:
                break
        if clash:
            for i, p in old_pos.items():
                pts[i] = p
            continue
        moved += 1
    return moved


def snap_alkyl_donor_to_metal_axis(
    xyz_str: str,
    mol=None,
    *,
    env_flag: str = "DELFIN_5F_E_ALKYL_DONOR_AXIS",
) -> str:
    """Welle-5f-E: rescue stretched M-CHx-X alkyl donors.

    For each C atom with exactly one heavy non-metal neighbour whose
    distance to the closest metal is in the band
    (_ALKYL_DONOR_LO, _ALKYL_DONOR_HI) * ideal_M-C, rigidly translate
    the C + H neighbours + X-side fragment along the (C->M) direction
    so that d_MC matches the ideal M-C(sp3) sigma length.

    Designed for the D-XALKEA Zr CN3 pattern: ETKDG/UFF leaves the
    alkyl C at ~2.40 A while the real sigma-bond is ~2.20 A and only
    one frame out of four happens to be axially correct.  After this
    pass, all qualifying frames are snapped to sigma length.

    Failsafe: any exception returns the input unchanged.
    Bit-exact when no candidate is in the trigger band.

    The ``env_flag`` parameter is informational only -- gating happens
    in the smiles_converter dispatch helper, not here.
    """
    try:
        syms, pts, orig_lines = _parse_xyz(xyz_str)
        if pts.shape[0] == 0:
            return xyz_str
        pts = pts.copy()
        moved = _snap_alkyl_axis_xyz(syms, pts)
        if moved == 0:
            return xyz_str
        return _format_xyz(orig_lines, syms, pts)
    except Exception:
        return xyz_str


def correct_results_alkyl_donor_axis(mol, results):
    """Apply alkyl-donor-axis snap to each (xyz, label) entry.

    Per-entry failsafe.  Signature mirrors ``correct_results`` so the
    smiles_converter dispatch helper can call it with the same
    arguments.
    """
    out = []
    for entry in results:
        try:
            xyz, lbl = entry[0], entry[1]
            new_xyz = snap_alkyl_donor_to_metal_axis(xyz, mol)
            out.append(
                (new_xyz, lbl) if len(entry) == 2
                else (new_xyz,) + tuple(entry[1:])
            )
        except Exception:
            out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Welle-5f-C — Post-UFF M-H atom-overlap rescue (D-QOPVEZ Fe-H 0.4 A pattern)
# ---------------------------------------------------------------------------
#
# Symptom: UFF (or UFF-soft donor pin) collapses one or more H atoms onto a
# metal centre — frames showing M-H below 1.3 A (catastrophic atom-overlap
# territory; real M-H ~1.5-1.85 A; see feedback_fe_h_atom_overlap_post_uff).
#
# Strategy:
#   1. Detect every (M, H) pair with d_MH < ``min_dist``.
#   2. Push the H radially outward along the M->H axis to the ideal M-H
#      bond length looked up in ``_M_H_IDEAL`` (per metal symbol).
#   3. Verify the rescue does NOT introduce a new heavy-atom clash
#      (vdW * ``_CLASH_VDW_FACTOR`` rule) and does NOT break any other
#      M-D bond (M-D invariant snapshot).
#   4. If a rescue is geometrically impossible (cannot find a clean push
#      direction within ``_RESCUE_N_TRIES`` perturbations of the radial
#      axis), the *frame* is discarded — caller receives ``None``.
#
# Universal — only depends on element symbols + bond-graph geometry.
# Env-flag default-OFF (``DELFIN_5F_C_MH_OVERLAP_RESCUE``).

# Per-metal ideal M-H bond length (Angstrom).  Subset of the master table
# in ``smiles_converter._METAL_HYDRIDE_BOND_LENGTHS``; duplicated here so
# this helper has no import-cycle on ``smiles_converter``.
_M_H_IDEAL: Dict[str, float] = {
    # 3d TM
    "Sc": 1.85, "Ti": 1.75, "V":  1.70, "Cr": 1.65, "Mn": 1.62,
    "Fe": 1.55, "Co": 1.53, "Ni": 1.50, "Cu": 1.55, "Zn": 1.60,
    # 4d TM
    "Y":  1.95, "Zr": 1.85, "Nb": 1.75, "Mo": 1.70, "Tc": 1.68,
    "Ru": 1.60, "Rh": 1.58, "Pd": 1.55, "Ag": 1.65, "Cd": 1.70,
    # 5d TM + La
    "La": 2.05, "Hf": 1.85, "Ta": 1.75, "W":  1.70, "Re": 1.68,
    "Os": 1.62, "Ir": 1.58, "Pt": 1.55, "Au": 1.62, "Hg": 1.70,
    # Main-group metals
    "Li": 1.60, "Na": 1.85, "K":  2.20, "Be": 1.35, "Mg": 1.70,
    "Ca": 2.00, "Sr": 2.15, "Ba": 2.30, "Al": 1.55, "Ga": 1.55,
    "In": 1.70, "Tl": 1.85, "Sn": 1.70, "Pb": 1.80, "Bi": 1.85,
}
_M_H_IDEAL_FALLBACK: float = 1.60  # generic mid-row TM

_RESCUE_DEFAULT_MIN_DIST: float = 1.30  # Å (rescue threshold)
_RESCUE_N_TRIES: int = 6                # axis perturbation attempts


def _ideal_m_h(metal_sym: str) -> float:
    return _M_H_IDEAL.get(metal_sym, _M_H_IDEAL_FALLBACK)


def _heavy_clash_after_move(
    new_pos: np.ndarray,
    moved_h_idx: int,
    syms: List[str],
    pts: np.ndarray,
    nbrs: List[List[int]],
    vdw_factor: float = _CLASH_VDW_FACTOR,
) -> bool:
    """Return True if placing ``moved_h_idx`` at ``new_pos`` would clash
    with any heavy atom that is not its parent (covalent neighbour)."""
    parent_set = set(nbrs[moved_h_idx])
    r_h = _VDW_RADII.get("H", 1.20)
    for j in range(len(syms)):
        if j == moved_h_idx or j in parent_set:
            continue
        sj = syms[j]
        if sj == "H":
            # H-H clash check is delegated to downstream H detectors;
            # here we only police heavy-atom overlap created by the rescue.
            continue
        r_j = _VDW_RADII.get(sj, 1.70)
        thr = vdw_factor * (r_h + r_j)
        d = float(np.linalg.norm(new_pos - pts[j]))
        if d < thr:
            return True
    return False


def rescue_m_h_atom_overlap(
    xyz_str: str,
    mol=None,
    *,
    min_dist: float = _RESCUE_DEFAULT_MIN_DIST,
) -> Optional[str]:
    """Post-UFF rescue for catastrophically short M-H bonds.

    For every (metal M, hydrogen H) pair where ``d(M, H) < min_dist``
    push H radially outward to the ideal M-H bond length from
    ``_M_H_IDEAL``.  If the rescue introduces a new heavy-atom clash
    or breaks any other M-D bond, retry with a perturbed direction up
    to ``_RESCUE_N_TRIES`` times.  If all attempts fail the frame is
    discarded (returns ``None``).

    Returns
    -------
    Optional[str]
        The rescued XYZ string (DELFIN format preserved), the original
        string when nothing needed rescuing, or ``None`` if at least
        one M-H overlap was unrescuable — caller must drop the frame.
    """
    try:
        syms, pts, orig_lines = _parse_xyz(xyz_str)
    except Exception:
        return xyz_str
    if pts.shape[0] == 0:
        return xyz_str

    # Quick scan: any (M, H) below threshold?
    metal_idxs = [i for i, s in enumerate(syms) if _is_metal_sym(s)]
    if not metal_idxs:
        return xyz_str
    h_idxs = [i for i, s in enumerate(syms) if s == "H"]
    if not h_idxs:
        return xyz_str

    overlap_pairs: List[Tuple[int, int, float]] = []
    for m in metal_idxs:
        for h in h_idxs:
            d = float(np.linalg.norm(pts[m] - pts[h]))
            if d < min_dist:
                overlap_pairs.append((m, h, d))
    if not overlap_pairs:
        return xyz_str

    # Snapshot M-D invariant on the ORIGINAL geometry; any rescue must
    # preserve every other M-D bond within ``_MD_INVARIANT_TOL``.
    md_pre = _snapshot_md_bonds(pts, syms)
    nbrs = _build_geometric_adjacency(syms, pts)

    new_pts = pts.copy()
    for m_idx, h_idx, _d in overlap_pairs:
        ideal = _ideal_m_h(syms[m_idx])
        m_pos = new_pts[m_idx]
        h_pos = new_pts[h_idx]

        # Primary radial direction (M -> H).  If degenerate (zero
        # vector — H sits exactly on M), seed with a stable axis.
        direction = h_pos - m_pos
        nd = float(np.linalg.norm(direction))
        if nd < 1e-6:
            direction = np.array([1.0, 0.0, 0.0])
        else:
            direction = direction / nd

        rescued = False
        # Deterministic perturbation list: identity + 6 orthogonal tilts.
        u, v = _orthogonal_basis(direction)
        candidates = [direction]
        for k in range(1, _RESCUE_N_TRIES):
            ang = (math.pi / 6.0) * k    # 30 deg, 60 deg, ...
            tilt = (
                math.cos(ang) * direction
                + math.sin(ang) * (u if k % 2 == 0 else v)
            )
            n_t = float(np.linalg.norm(tilt))
            if n_t > 1e-9:
                candidates.append(tilt / n_t)

        for axis in candidates:
            cand_pos = m_pos + axis * ideal
            # Try-move and validate.
            trial = new_pts.copy()
            trial[h_idx] = cand_pos
            if _heavy_clash_after_move(
                cand_pos, h_idx, syms, trial, nbrs
            ):
                continue
            if _md_invariant_violated(md_pre, trial):
                continue
            new_pts = trial
            rescued = True
            break

        if not rescued:
            # Cannot rescue this pair safely — discard whole frame.
            return None

    return _format_xyz(orig_lines, syms, new_pts)


def rescue_results(mol, results):
    """Apply M-H overlap rescue to each (xyz, label) entry.

    Frames whose rescue is geometrically impossible are dropped — they
    would otherwise propagate catastrophic atom-overlap into downstream
    scoring.  Failsafe per-entry: any exception keeps the original entry.
    Env-flag-gated by the caller (``DELFIN_5F_C_MH_OVERLAP_RESCUE``).
    """
    out = []
    for entry in results:
        try:
            xyz = entry[0]
            new_xyz = rescue_m_h_atom_overlap(xyz, mol)
            if new_xyz is None:
                # Drop frame entirely — rescue impossible.
                continue
            out.append(
                (new_xyz,) + tuple(entry[1:])
                if len(entry) > 1 else (new_xyz,)
            )
        except Exception:
            out.append(entry)
    return out
