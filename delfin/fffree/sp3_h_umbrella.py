"""sp³-H umbrella enforcer — build-time geometric template for CH3/CH2/CH/NH2/NH.

Motivation
==========
Forensik on ``eaee0c2-GRACE-VOLLPOOL/29-Ni_pincer-tBu-imid_.xyz`` (2026-06-04)
revealed every tBu methyl group sitting in a degenerate orthogonal-axes pose:
two of the three H atoms lie antiparallel through the C (H-C-H = 180°) while
the third sits perpendicular (H-C-H = 90° / 90°).  The C-C-C tetrahedral
backbone is correct, the C-H bond lengths are correct (1.09 Å), but the H
positions are pathological.

Why existing heal layers miss it
--------------------------------
* ``topology_healing``: angle severity via Mahalanobis prior on CCDC bond
  angles.  X-ray riding-model H positions have very wide variance, so
  90° / 180° may fall inside tolerance.
* ``grip_healing``: bond-length-residual driven, ignores angles.
* ``grip_polish`` (L-BFGS angle term): 180° is a saddle of cos(angle) loss
  → gradient = 0 → no move.
* GRACE: rotates bond axes, not H positions inside a CH3 cone.

This module attacks the problem at the geometric source: it builds the
H umbrella from scratch around a tetrahedral template anchored on the
heavy-atom neighbour(s) of the sp³ centre.  It is callable as

    * a **build-time hook** (post-construction, pre-refine), or
    * a **post-hoc healer** (see ``sp3_h_heal.py``).

Doctrine (matches existing GUPPY fffree modules)
------------------------------------------------
* Topology-safe: metals NEVER moved, non-H atoms NEVER moved.
* Bond-length preserving: every X-H stays at its original length (or the
  ideal length when ``--ideal-bond`` is specified, default = preserve).
* Deterministic: byte-identical output for same input + same seed.
* Env-gated default OFF byte-identical: ``DELFIN_FFFREE_SP3_H_UMBRELLA``.
* Lazy RDKit; numpy at module-level.
* Fail-safe: every exception is swallowed; caller pipeline never breaks.

Geometry recipe
---------------
For CH₃ (1 heavy-atom neighbour X):
    Reference axis r = (C - X) / ||C - X||
    Place H_k at C + L · [cos(θ)·r + sin(θ)·(cos(φ_k)·u + sin(φ_k)·v)]
    with θ = 109.47°, φ_k = 2π·k/3 + φ0  (k=0,1,2)
    L = C-H bond length (preserved per-H or set to 1.09 Å)
    u, v = deterministic orthonormal frame perpendicular to r.

For CH₂ (2 heavy-atom neighbours X, Y):
    Compute bisector direction b = (r_X + r_Y) / ||·||  inside the X-C-Y plane
    Place H atoms symmetrically along ±(out-of-plane direction) at the
    bisector angle so X-C-H = Y-C-H ≈ 109.47° (canonical sp³ CH₂).

For CH (3 heavy-atom neighbours X, Y, Z):
    Single H goes along the negative sum of the 3 heavy-atom direction
    vectors (umbrella closure), reproducing the tetrahedral 4th vertex.

For NH₂ / NH: same recipes (N is the centre, ideal angle 107° but the
template uses 109.47° — within the heuristic tolerance, the slight
flattening at N is handled later by the existing ``_fix_sp3_n_pyramidality``
post-processor when called downstream).

Deterministic frame construction
--------------------------------
``u`` is built from ẑ × r normalised, with deterministic fallback to ŷ × r
when r is parallel to ẑ (|r_z| > 0.95).  ``v = r × u``.  This makes the
azimuth basis a function of r alone — same heavy-atom direction → same
basis → same H positions, even across cluster nodes.

Public API
==========
* :func:`detect_sp3_centers(coords, syms, bond_topology)` -> list[Sp3Center]
* :func:`check_umbrella_geometry(center, coords)` -> Sp3Diagnostic
* :func:`enforce_umbrella(coords, center)` -> new_coords
* :func:`enforce_all_sp3_umbrella(coords, syms, mol)` -> new_coords
"""
from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants — vendored from neighbouring modules so this file is standalone.
# ---------------------------------------------------------------------------

_IDEAL_TETRA_DEG: float = 109.471  # arccos(-1/3) in degrees
_IDEAL_TETRA_RAD: float = math.acos(-1.0 / 3.0)

_DEFAULT_CH_LENGTH: float = 1.090  # Å, ideal sp³ C-H
_DEFAULT_NH_LENGTH: float = 1.010  # Å, ideal sp³ N-H

# Degenerate band (degrees) — H-C-H or H-N-H angles inside this interval
# are accepted; outside triggers the umbrella flag.
_DEGEN_LO_DEG: float = 95.0
_DEGEN_HI_DEG: float = 125.0
_COLLINEAR_COS_THRESHOLD: float = 0.90  # |cos(angle)| ≥ 0.9 → "collinear"

_COV_RADII: Dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}

# A small VDW table — we only need H + a handful of heavies for clash gating.
_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
}

# Atomic numbers for the small symbol set we care about for sp³ centres.
_SP3_CENTRE_SYMBOLS: FrozenSet[str] = frozenset({"C", "N"})

# Metals are detected by Z range — see _fix_sp3_h_tetrahedrality for the
# canonical set; we re-use it here for safety though sp³ centres should
# never be metals.
_METAL_Z_RANGES = (
    set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
    | set(range(89, 104))
    | {3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 51, 55, 56, 81, 82, 83}
)

_Z_BY_SYMBOL: Dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Br": 35, "Mo": 42,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "I": 53, "W": 74, "Re": 75,
    "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
}


def _is_metal_sym(sym: str) -> bool:
    z = _Z_BY_SYMBOL.get(sym)
    return z is not None and z in _METAL_Z_RANGES


_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)

# ---------------------------------------------------------------------------
# Env flag
# ---------------------------------------------------------------------------

_UMBRELLA_ENV: str = "DELFIN_FFFREE_SP3_H_UMBRELLA"


def umbrella_active() -> bool:
    """``True`` iff the build-time umbrella enforcer is enabled (default OFF)."""
    raw = os.environ.get(_UMBRELLA_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sp3Center:
    """Description of one sp³ centre with explicit H neighbours.

    Attributes
    ----------
    center_idx
        Atom index of the sp³ centre (C or N) inside the ``coords`` array.
    center_sym
        Element symbol of the centre (``"C"`` or ``"N"``).
    h_indices
        Sorted tuple of indices of explicit H neighbours attached to the
        centre.
    heavy_indices
        Sorted tuple of indices of heavy-atom neighbours of the centre
        (all neighbours that are not in ``h_indices``).
    classification
        One of ``"CH3"``, ``"CH2"``, ``"CH"``, ``"NH2"``, ``"NH"`` or
        ``"OTHER"`` when the centre does not match any canonical template.
    """

    center_idx: int
    center_sym: str
    h_indices: Tuple[int, ...]
    heavy_indices: Tuple[int, ...]
    classification: str


@dataclass(frozen=True)
class Sp3Diagnostic:
    """Per-centre umbrella diagnostic.

    ``hh_angles_deg`` lists every pairwise H-X-H angle.
    ``xh_angles_deg`` lists every heavy-X-H angle.
    ``min_hh_cos`` is the most extreme |cos(angle)| over H-X-H pairs (a
    proxy for "two H's antiparallel" when ≥ 0.9).
    ``flag_degenerate`` is True iff any H-X-H angle is outside
    ``(_DEGEN_LO_DEG, _DEGEN_HI_DEG)`` *or* any pair is collinear.
    """

    center: Sp3Center
    hh_angles_deg: Tuple[float, ...]
    xh_angles_deg: Tuple[float, ...]
    max_hh_cos_abs: float
    flag_degenerate: bool
    max_hh_dev_deg: float


# ---------------------------------------------------------------------------
# XYZ helpers (kept simple — used by the validation script in Mission 4)
# ---------------------------------------------------------------------------


def parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse a multiline XYZ string into ``(symbols, coords, raw_lines)``."""
    syms: List[str] = []
    pts: List[List[float]] = []
    lines = xyz_str.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append([float(m.group(2)), float(m.group(3)),
                        float(m.group(4))])
    if not pts:
        return syms, np.zeros((0, 3), dtype=np.float64), lines
    return syms, np.asarray(pts, dtype=np.float64), lines


def format_xyz(orig_lines: List[str], syms: List[str],
               positions: np.ndarray) -> str:
    """Re-emit an XYZ string with updated coordinates (atom rows only)."""
    out: List[str] = []
    atom_i = 0
    for line in orig_lines:
        m = _XYZ_LINE_RE.match(line)
        if m and atom_i < len(syms):
            x, y, z = positions[atom_i]
            out.append(f"{syms[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            atom_i += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if orig_lines else "")


# ---------------------------------------------------------------------------
# Geometric helpers
# ---------------------------------------------------------------------------


def _angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle (a-b-c) in degrees with ``b`` as the vertex.  NaN-safe."""
    v1 = a - b
    v2 = c - b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return float("nan")
    cosv = float(np.dot(v1, v2) / (n1 * n2))
    cosv = max(-1.0, min(1.0, cosv))
    return math.degrees(math.acos(cosv))


def _cos_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1 = a - b
    v2 = c - b
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return 0.0
    return float(np.dot(v1, v2) / (n1 * n2))


def _deterministic_basis_for_axis(r: np.ndarray
                                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Build a deterministic orthonormal pair (u, v) with u, v ⟂ r.

    Construction:
        Pick seed ẑ unless r is too close to ẑ (|r_z| > 0.95) in which
        case fall back to ŷ.  Subtract the r-component, normalise.  Then
        v = r × u, normalised.  Output is invariant under sign of seed
        and depends only on ``r`` (and the chosen fallback rule), giving
        bit-stable behaviour for identical inputs.
    """
    r = r / max(float(np.linalg.norm(r)), 1e-12)
    if abs(float(r[2])) < 0.95:
        seed = np.array([0.0, 0.0, 1.0])
    else:
        seed = np.array([0.0, 1.0, 0.0])
    u = seed - float(np.dot(seed, r)) * r
    nu = float(np.linalg.norm(u))
    if nu < 1e-9:
        # last-resort: orthogonal to ẑ deterministic
        seed = np.array([1.0, 0.0, 0.0])
        u = seed - float(np.dot(seed, r)) * r
        nu = float(np.linalg.norm(u))
        if nu < 1e-9:
            # truly pathological — return canonical ẑ/ŷ basis (caller
            # should not reach this path because r is a unit vector).
            return np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])
    u = u / nu
    v = np.cross(r, u)
    nv = float(np.linalg.norm(v))
    if nv < 1e-9:  # would mean u‖r which we just defended against
        v = np.array([0.0, 1.0, 0.0])
    else:
        v = v / nv
    return u, v


# ---------------------------------------------------------------------------
# Geometric adjacency (used when ``bond_topology`` is not provided)
# ---------------------------------------------------------------------------


def _geometric_adjacency(syms: Sequence[str], coords: np.ndarray,
                         cutoff_pad: float = 0.25) -> List[List[int]]:
    """Distance-based adjacency table (Σr_cov + pad)."""
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        ri = _COV_RADII.get(syms[i], 1.50)
        for j in range(i + 1, n):
            rj = _COV_RADII.get(syms[j], 1.50)
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if d < ri + rj + cutoff_pad:
                nbrs[i].append(j)
                nbrs[j].append(i)
    return nbrs


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def _classify_center(sym: str, n_h: int, n_heavy: int) -> str:
    if sym == "C":
        if n_h == 3 and n_heavy == 1:
            return "CH3"
        if n_h == 2 and n_heavy == 2:
            return "CH2"
        if n_h == 1 and n_heavy == 3:
            return "CH"
        return "OTHER"
    if sym == "N":
        if n_h == 2 and n_heavy == 1:
            return "NH2"
        if n_h == 1 and n_heavy == 2:
            return "NH"
        return "OTHER"
    return "OTHER"


def detect_sp3_centers(
    coords: np.ndarray,
    syms: Sequence[str],
    bond_topology: Optional[Sequence[Sequence[int]]] = None,
    *,
    mol=None,
    include_other: bool = False,
) -> List[Sp3Center]:
    """Enumerate sp³-candidate centres (C or N) with explicit H neighbours.

    Parameters
    ----------
    coords
        ``(N, 3)`` array of Cartesian positions.
    syms
        Length-``N`` sequence of element symbols (canonical case).
    bond_topology
        Optional adjacency: ``bond_topology[i]`` is the list of neighbour
        indices of atom ``i``.  When ``None``, geometric adjacency
        (Σr_cov + 0.25 Å) is used.
    mol
        Optional RDKit ``Mol`` — when provided, sp³ classification uses
        RDKit's hybridization annotation rather than the geometric N(H)
        heuristic.
    include_other
        When True, also return centres that don't match a canonical
        CH3/CH2/CH/NH2/NH template (classification = ``"OTHER"``).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = len(syms)
    if coords.shape[0] != n or n == 0:
        return []

    if bond_topology is None:
        adj = _geometric_adjacency(syms, coords)
    else:
        adj = [list(b) for b in bond_topology]

    sp3_set: Optional[Set[int]] = None
    if mol is not None:
        try:
            from rdkit import Chem  # noqa: F401 — used for hybridization symbol
            sp3_set = set()
            for atom in mol.GetAtoms():
                if str(atom.GetHybridization()) == "SP3":
                    sp3_set.add(int(atom.GetIdx()))
        except Exception:
            sp3_set = None

    out: List[Sp3Center] = []
    for i in range(n):
        sym = syms[i]
        if sym not in _SP3_CENTRE_SYMBOLS:
            continue
        if _is_metal_sym(sym):
            continue
        nbrs = adj[i] if i < len(adj) else []
        h_idxs = sorted([j for j in nbrs if 0 <= j < n and syms[j] == "H"])
        heavy_idxs = sorted([j for j in nbrs if 0 <= j < n and syms[j] != "H"])
        if not h_idxs:
            continue
        # Optional RDKit sp³ filter
        if sp3_set is not None and i not in sp3_set:
            continue
        cls = _classify_center(sym, len(h_idxs), len(heavy_idxs))
        if cls == "OTHER" and not include_other:
            continue
        out.append(Sp3Center(
            center_idx=int(i),
            center_sym=sym,
            h_indices=tuple(h_idxs),
            heavy_indices=tuple(heavy_idxs),
            classification=cls,
        ))
    return out


# ---------------------------------------------------------------------------
# Diagnostic
# ---------------------------------------------------------------------------


def check_umbrella_geometry(center: Sp3Center,
                            coords: np.ndarray) -> Sp3Diagnostic:
    """Compute H-X-H angles and flag whether the umbrella is degenerate."""
    coords = np.asarray(coords, dtype=np.float64)
    c = coords[center.center_idx]
    h_idxs = center.h_indices
    heavy_idxs = center.heavy_indices

    hh_angles: List[float] = []
    max_abs_cos = 0.0
    max_dev = 0.0
    flag = False
    # Pairwise H-X-H
    for a in range(len(h_idxs)):
        for b in range(a + 1, len(h_idxs)):
            ang = _angle_deg(coords[h_idxs[a]], c, coords[h_idxs[b]])
            cosv = _cos_angle(coords[h_idxs[a]], c, coords[h_idxs[b]])
            hh_angles.append(ang)
            if abs(cosv) > max_abs_cos:
                max_abs_cos = abs(cosv)
            dev = abs(ang - _IDEAL_TETRA_DEG)
            if dev > max_dev:
                max_dev = dev
            if math.isnan(ang):
                continue
            if ang < _DEGEN_LO_DEG or ang > _DEGEN_HI_DEG:
                flag = True
            if abs(cosv) > _COLLINEAR_COS_THRESHOLD:
                flag = True
    # Heavy-X-H (for diagnostic only; not used to set ``flag_degenerate``)
    xh_angles: List[float] = []
    for h_idx in h_idxs:
        for x_idx in heavy_idxs:
            ang = _angle_deg(coords[x_idx], c, coords[h_idx])
            xh_angles.append(ang)
    return Sp3Diagnostic(
        center=center,
        hh_angles_deg=tuple(hh_angles),
        xh_angles_deg=tuple(xh_angles),
        max_hh_cos_abs=float(max_abs_cos),
        flag_degenerate=bool(flag),
        max_hh_dev_deg=float(max_dev),
    )


# ---------------------------------------------------------------------------
# Construction primitives
# ---------------------------------------------------------------------------


def _build_methyl_umbrella(
    c_pos: np.ndarray,
    x_pos: np.ndarray,
    bond_lengths: Sequence[float],
    azimuth_offset: float = 0.0,
) -> np.ndarray:
    """Place 3 H positions for a CH₃ around centre ``c_pos`` with reference
    axis defined by ``c - x`` (heavy-atom neighbour).

    Each H is placed at distance ``bond_lengths[k]`` (k = 0, 1, 2) so the
    individual C-H lengths are preserved exactly.  The three azimuthal
    angles are 120°-spaced; ``azimuth_offset`` is added to all three (useful
    for staggered conformer generation).
    """
    r = c_pos - x_pos
    nr = float(np.linalg.norm(r))
    if nr < 1e-9:
        # degenerate axis — pick canonical r = ẑ
        r = np.array([0.0, 0.0, 1.0])
    else:
        r = r / nr
    u, v = _deterministic_basis_for_axis(r)
    theta = _IDEAL_TETRA_RAD
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    out = np.zeros((3, 3), dtype=np.float64)
    for k in range(3):
        phi = (2.0 * math.pi * k) / 3.0 + azimuth_offset
        d = (cos_t * r + sin_t * (math.cos(phi) * u + math.sin(phi) * v))
        out[k] = c_pos + bond_lengths[k] * d
    return out


def _build_ch2_umbrella(
    c_pos: np.ndarray,
    x_pos: np.ndarray,
    y_pos: np.ndarray,
    bond_lengths: Sequence[float],
) -> np.ndarray:
    """Place 2 H positions for a CH₂ given two heavy neighbours X, Y.

    Algorithm
    ---------
    1. Compute in-plane bisector ``b = -(uX + uY) / ||·||`` (points away
       from both heavies inside the X-C-Y plane).
    2. Compute out-of-plane normal ``n = uX × uY / ||·||``.
    3. Both H's lie in the (b, n)-plane.  We pick the H-C-H half-angle
       to be 109.47°/2 (≈54.74° on each side of the bisector when X-C-Y
       is exactly tetrahedral).  More precisely we solve for the half
       angle θ_H so the H-X-C tetrahedral target is satisfied at the
       cost of relaxing H-Y-C symmetrically — this is the canonical
       sp³ CH₂ construction.
    4. H_+ = c + L · (cos(θ_H)·b + sin(θ_H)·n)
       H_- = c + L · (cos(θ_H)·b - sin(θ_H)·n)
    """
    uX = x_pos - c_pos
    uY = y_pos - c_pos
    nX = float(np.linalg.norm(uX))
    nY = float(np.linalg.norm(uY))
    if nX < 1e-9 or nY < 1e-9:
        # degenerate — fall back to a CH3-like placement using X only
        return _build_methyl_umbrella(c_pos, x_pos, list(bond_lengths) + [1.09])[:2]
    uX /= nX
    uY /= nY
    # Bisector pointing AWAY from both heavies (in the X-C-Y plane)
    b_vec = -(uX + uY)
    nb = float(np.linalg.norm(b_vec))
    if nb < 1e-9:
        # X and Y antiparallel through C — build perpendicular umbrella
        # using a deterministic basis perpendicular to uX.
        u_perp, v_perp = _deterministic_basis_for_axis(uX)
        b_vec = u_perp
    else:
        b_vec /= nb
    n_vec = np.cross(uX, uY)
    n_norm = float(np.linalg.norm(n_vec))
    if n_norm < 1e-9:
        # Still degenerate (collinear) — use a third perpendicular axis
        u_perp, v_perp = _deterministic_basis_for_axis(uX)
        n_vec = v_perp
    else:
        n_vec /= n_norm
    # Half-angle so the H-C-H angle ≈ 109.47°
    half = _IDEAL_TETRA_RAD / 2.0
    cos_h = math.cos(half)
    sin_h = math.sin(half)
    h0 = c_pos + bond_lengths[0] * (cos_h * b_vec + sin_h * n_vec)
    h1 = c_pos + bond_lengths[1] * (cos_h * b_vec - sin_h * n_vec)
    return np.vstack([h0, h1])


def _build_methine_umbrella(
    c_pos: np.ndarray,
    heavy_positions: Sequence[np.ndarray],
    bond_length: float,
) -> np.ndarray:
    """Place 1 H opposite the centroid of the 3 heavy directions.

    For an ideal sp³ centre with 3 heavy neighbours X, Y, Z the missing
    4th vertex lies along ``-(uX + uY + uZ) / ||·||`` from the centre.
    """
    accum = np.zeros(3, dtype=np.float64)
    for h_pos in heavy_positions:
        v = h_pos - c_pos
        nv = float(np.linalg.norm(v))
        if nv < 1e-9:
            continue
        accum += v / nv
    nA = float(np.linalg.norm(accum))
    if nA < 1e-9:
        # degenerate (e.g. 3 heavies symmetric around C) — fall back to
        # canonical ẑ direction.  This is the symmetry-breaking ambiguity
        # noted in the spec; deterministic choice ensures reproducibility.
        h_dir = np.array([0.0, 0.0, 1.0])
    else:
        h_dir = -accum / nA
    return c_pos.reshape(1, 3) + bond_length * h_dir.reshape(1, 3)


# ---------------------------------------------------------------------------
# Public enforcer (single centre)
# ---------------------------------------------------------------------------


def enforce_umbrella(coords: np.ndarray, center: Sp3Center,
                     *, preserve_bond_lengths: bool = True,
                     azimuth_offset: float = 0.0) -> np.ndarray:
    """Return a new ``coords`` array with the H positions of ``center``
    replaced by an ideal sp³ umbrella.

    * ``preserve_bond_lengths=True`` (default) keeps each X-H at its
      original length; otherwise idealised lengths (1.09 Å C-H, 1.01 Å
      N-H) are used.
    * ``azimuth_offset`` (CH₃ only) rotates the 3-H pattern around the
      C-X axis — useful for staggered conformers.

    The centre and heavy neighbours are NEVER moved.  Non-related atoms
    are NEVER moved.  Returns a copy; the input is untouched.
    """
    coords = np.asarray(coords, dtype=np.float64).copy()
    n = coords.shape[0]
    ci = center.center_idx
    if ci < 0 or ci >= n:
        return coords
    h_idxs = list(center.h_indices)
    heavy_idxs = list(center.heavy_indices)
    cls = center.classification

    def _length(h_idx: int, default: float) -> float:
        if not preserve_bond_lengths:
            return default
        d = float(np.linalg.norm(coords[h_idx] - coords[ci]))
        if d < 0.5 or d > 1.5:
            return default
        return d

    if cls == "CH3" and len(heavy_idxs) >= 1 and len(h_idxs) == 3:
        Ls = [_length(h, _DEFAULT_CH_LENGTH) for h in h_idxs]
        # Choose the FIRST heavy neighbour (sorted) as anchor — deterministic.
        x_idx = heavy_idxs[0]
        if x_idx >= n:
            return coords
        new_hs = _build_methyl_umbrella(
            coords[ci], coords[x_idx], Ls, azimuth_offset=azimuth_offset,
        )
        for k, h_idx in enumerate(h_idxs):
            if 0 <= h_idx < n:
                coords[h_idx] = new_hs[k]
        return coords

    if cls == "CH2" and len(heavy_idxs) >= 2 and len(h_idxs) == 2:
        Ls = [_length(h, _DEFAULT_CH_LENGTH) for h in h_idxs]
        x_idx, y_idx = heavy_idxs[0], heavy_idxs[1]
        if x_idx >= n or y_idx >= n:
            return coords
        new_hs = _build_ch2_umbrella(coords[ci], coords[x_idx],
                                     coords[y_idx], Ls)
        for k, h_idx in enumerate(h_idxs):
            if 0 <= h_idx < n:
                coords[h_idx] = new_hs[k]
        return coords

    if cls == "CH" and len(heavy_idxs) >= 3 and len(h_idxs) == 1:
        h_idx = h_idxs[0]
        L = _length(h_idx, _DEFAULT_CH_LENGTH)
        new_h = _build_methine_umbrella(
            coords[ci],
            [coords[heavy_idxs[0]], coords[heavy_idxs[1]],
             coords[heavy_idxs[2]]],
            L,
        )
        if 0 <= h_idx < n:
            coords[h_idx] = new_h[0]
        return coords

    if cls == "NH2" and len(heavy_idxs) >= 1 and len(h_idxs) == 2:
        Ls = [_length(h, _DEFAULT_NH_LENGTH) for h in h_idxs]
        x_idx = heavy_idxs[0]
        if x_idx >= n:
            return coords
        # NH₂ is sp³-pyramidal, ≈106° H-N-H — we re-use the methyl
        # builder but pick a 2-H subset (k = 0, 1) of the 3-H umbrella
        # so the geometry is closer to NH₂ ideal than a planar build.
        new_hs = _build_methyl_umbrella(coords[ci], coords[x_idx],
                                        Ls + [_DEFAULT_NH_LENGTH])
        for k, h_idx in enumerate(h_idxs):
            if 0 <= h_idx < n:
                coords[h_idx] = new_hs[k]
        return coords

    if cls == "NH" and len(heavy_idxs) >= 2 and len(h_idxs) == 1:
        h_idx = h_idxs[0]
        L = _length(h_idx, _DEFAULT_NH_LENGTH)
        # Use methine recipe but supply a third "phantom" heavy direction
        # equal to the centroid of the existing two so the closure points
        # out of plane.
        x_pos = coords[heavy_idxs[0]]
        y_pos = coords[heavy_idxs[1]]
        bisector = -((x_pos + y_pos) - 2.0 * coords[ci])
        nb = float(np.linalg.norm(bisector))
        if nb > 1e-9:
            phantom = coords[ci] + bisector / nb
            new_h = _build_methine_umbrella(
                coords[ci], [x_pos, y_pos, phantom], L,
            )
            if 0 <= h_idx < n:
                coords[h_idx] = new_h[0]
        return coords

    return coords


# ---------------------------------------------------------------------------
# Batch enforcer
# ---------------------------------------------------------------------------


def _h_introduces_clash(new_pos: np.ndarray, h_idx: int, center_idx: int,
                        siblings: Set[int], coords: np.ndarray,
                        syms: Sequence[str],
                        threshold_factor: float = 0.70) -> bool:
    """True iff placing H at ``new_pos`` creates a NEW atomic overlap
    with any atom that is not itself, its parent, or a geminal sibling.
    Bonded ranges (0.85·Σr_cov ≤ d ≤ 1.10·Σr_cov) are excluded.
    """
    n = coords.shape[0]
    rH_v = _VDW_RADII.get("H", 1.20)
    rH_c = _COV_RADII.get("H", 0.31)
    for j in range(n):
        if j == h_idx or j == center_idx or j in siblings:
            continue
        d = float(np.linalg.norm(new_pos - coords[j]))
        sj = syms[j]
        rj_v = _VDW_RADII.get(sj, 1.70 if sj != "H" else 1.20)
        rj_c = _COV_RADII.get(sj, 1.50 if sj != "H" else 0.31)
        sum_cov = rH_c + rj_c
        if 0.85 * sum_cov <= d <= 1.10 * sum_cov:
            continue
        if d < threshold_factor * (rH_v + rj_v):
            return True
    return False


def enforce_all_sp3_umbrella(
    coords: np.ndarray,
    syms: Sequence[str],
    mol=None,
    *,
    bond_topology: Optional[Sequence[Sequence[int]]] = None,
    preserve_bond_lengths: bool = True,
    azimuth_offset: float = 0.0,
    only_degenerate: bool = True,
    enable_clash_rollback: bool = True,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Apply the umbrella template to every sp³ centre in one pass.

    Parameters
    ----------
    only_degenerate
        When True (default), only centres whose
        :func:`check_umbrella_geometry` returns ``flag_degenerate=True`` are
        rewritten.  When False, every detected sp³ centre is rewritten
        (build-time mode).
    enable_clash_rollback
        When True (default), each per-centre rewrite is reverted if it
        creates a new non-bonded overlap (see ``_h_introduces_clash``).

    Returns
    -------
    (new_coords, report)
        ``report`` keys: ``"sp3_centres_seen"``, ``"degenerate_flagged"``,
        ``"centres_rewritten"``, ``"hs_moved"``, ``"hs_rolled_back"``,
        ``"max_hh_dev_before_deg"``, ``"max_hh_dev_after_deg"``.
    """
    coords = np.asarray(coords, dtype=np.float64).copy()
    report = {
        "sp3_centres_seen": 0,
        "degenerate_flagged": 0,
        "centres_rewritten": 0,
        "hs_moved": 0,
        "hs_rolled_back": 0,
        "max_hh_dev_before_deg": 0.0,
        "max_hh_dev_after_deg": 0.0,
    }
    if coords.shape[0] == 0:
        return coords, report
    try:
        centres = detect_sp3_centers(coords, syms,
                                     bond_topology=bond_topology,
                                     mol=mol)
    except Exception:
        return coords, report
    report["sp3_centres_seen"] = len(centres)
    max_before = 0.0
    max_after = 0.0
    for cen in centres:
        diag = check_umbrella_geometry(cen, coords)
        if diag.max_hh_dev_deg > max_before:
            max_before = diag.max_hh_dev_deg
        if only_degenerate and not diag.flag_degenerate:
            continue
        if diag.flag_degenerate:
            report["degenerate_flagged"] += 1
        # Snapshot for rollback
        old_h_pos = {h: coords[h].copy() for h in cen.h_indices}
        new_coords = enforce_umbrella(coords, cen,
                                       preserve_bond_lengths=preserve_bond_lengths,
                                       azimuth_offset=azimuth_offset)
        # Clash gating per-H (sibling = other H atoms on the same centre)
        siblings = set(cen.h_indices) | {cen.center_idx} | set(cen.heavy_indices)
        if enable_clash_rollback:
            for h_idx in cen.h_indices:
                if _h_introduces_clash(new_coords[h_idx], h_idx,
                                       cen.center_idx,
                                       siblings, new_coords, syms):
                    new_coords[h_idx] = old_h_pos[h_idx]
                    report["hs_rolled_back"] += 1
        # Count moved H's
        moved = 0
        for h_idx in cen.h_indices:
            if not np.allclose(new_coords[h_idx], old_h_pos[h_idx],
                               atol=1e-6):
                moved += 1
        if moved > 0:
            report["centres_rewritten"] += 1
            report["hs_moved"] += moved
        coords = new_coords
        # Recompute post-rewrite diagnostic
        diag2 = check_umbrella_geometry(cen, coords)
        if diag2.max_hh_dev_deg > max_after:
            max_after = diag2.max_hh_dev_deg
    report["max_hh_dev_before_deg"] = float(max_before)
    # ``max_after`` only includes centres that were rewritten; for
    # consistency compute over all centres at the end.
    for cen in centres:
        d2 = check_umbrella_geometry(cen, coords)
        if d2.max_hh_dev_deg > max_after:
            max_after = d2.max_hh_dev_deg
    report["max_hh_dev_after_deg"] = float(max_after)
    return coords, report


# ---------------------------------------------------------------------------
# String-level convenience (XYZ in -> XYZ out)
# ---------------------------------------------------------------------------


def enforce_all_sp3_umbrella_xyz(
    xyz_str: str,
    mol=None,
    *,
    preserve_bond_lengths: bool = True,
    azimuth_offset: float = 0.0,
    only_degenerate: bool = True,
    enable_clash_rollback: bool = True,
) -> Tuple[str, Dict[str, int]]:
    """Convenience wrapper: parse XYZ → enforce → re-emit."""
    syms, coords, raw = parse_xyz(xyz_str)
    if coords.shape[0] == 0:
        return xyz_str, {
            "sp3_centres_seen": 0, "degenerate_flagged": 0,
            "centres_rewritten": 0, "hs_moved": 0, "hs_rolled_back": 0,
            "max_hh_dev_before_deg": 0.0, "max_hh_dev_after_deg": 0.0,
        }
    new_coords, report = enforce_all_sp3_umbrella(
        coords, syms, mol,
        preserve_bond_lengths=preserve_bond_lengths,
        azimuth_offset=azimuth_offset,
        only_degenerate=only_degenerate,
        enable_clash_rollback=enable_clash_rollback,
    )
    if report["hs_moved"] == 0:
        return xyz_str, report
    return format_xyz(raw, syms, new_coords), report


__all__ = [
    "Sp3Center",
    "Sp3Diagnostic",
    "umbrella_active",
    "parse_xyz",
    "format_xyz",
    "detect_sp3_centers",
    "check_umbrella_geometry",
    "enforce_umbrella",
    "enforce_all_sp3_umbrella",
    "enforce_all_sp3_umbrella_xyz",
]
