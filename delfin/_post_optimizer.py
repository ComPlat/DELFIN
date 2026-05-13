"""Baustein 5 — PBD-based Post-UFF Geometry Corrector.

Position-Based Dynamics (PBD) style relaxation operating on already-finalised
XYZ strings (post UFF, post snap, post clash-relief, post dual-parse).  The
optimizer combines four elements:

    1. Phase A — One-shot topology repair.  Catastrophic M-D bond breaks
       (|M-D| > 1.30 · ideal) are healed by rigid-body translation of the
       offending donor's BFS fragment back toward the metal; if translation
       introduces clashes, a rotation about the M-axis is tried instead.

    2. Phase A.5 — Symmetry projection.  Per metal centre, ideal donor unit
       vectors are pulled from ``_TOPO_GEOMETRY_VECTORS`` (smiles_converter)
       for the matched (CN, geometry) and a Hungarian assignment matches
       current donors to ideal slots.  Each donor's rigid fragment is then
       translated a fraction ``step_size`` toward its ideal position.

    3. Phase B — Iterative Gauss-Seidel constraint resolution.  Up to
       ``max_iter`` sweeps over bonds, optional angles, and clash pairs.
       Each constraint is resolved as a mass-weighted projection that
       conserves the centre of mass of the constrained pair (or rigid
       rotation around an axis for angle / clash cases).

    4. Phase C — Hard topology gate.  After every commit, ``_passes_topology``
       verifies that all M-D distances stay within [0.85, 1.10] of ideal
       and that no non-bonded heavy-pair has collapsed inside 0.85·Σr_cov.
       Failed steps are rolled back atomically; if the total violation
       budget rises by > 5 %, the step is reverted and ``step_size`` halved.

Doctrine
--------
- Universal (no SMILES/refcode/element-list shortcuts).  Donor classification
  is purely structural: bond connectivity to a metal atom.
- Hard topology preservation.  A commit may only land if it preserves the
  graph: no M-D break, no new spurious heavy-heavy bond.
- Mass-weighted moves keep heavy fragments still and light atoms mobile,
  matching the rigid-body intuition of UFF post-relaxation.
- Standalone module.  Imports from delfin.smiles_converter for the master
  constants (``_TOPO_GEOMETRY_VECTORS``, ``_get_ml_bond_length``,
  ``_METAL_SET``, ``_COVALENT_RADII``).  Optional helpers from
  ``delfin._atom_weights`` and ``delfin._vdw_radii`` are imported lazily;
  inline fallbacks are used when those modules are not yet available.
- Always returns a valid XYZ string.  On any internal exception the original
  XYZ is returned unchanged with a populated report dict.

Entry point
-----------
``post_optimize_geometry(xyz, mol, class_label=..., **params) -> (xyz, report)``
"""
from __future__ import annotations

import math
import re
import traceback
from collections import deque
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment  # type: ignore
    from scipy.spatial.transform import Rotation as _Rotation  # type: ignore
    _SCIPY_OK = True
except Exception:  # pragma: no cover - SciPy expected in DELFIN env
    _SCIPY_OK = False

# ---------------------------------------------------------------------------
# Constants from delfin.smiles_converter (single source of truth).
# Fallbacks are conservative — only used in standalone testing of this file.
# ---------------------------------------------------------------------------
try:
    from delfin.smiles_converter import (  # type: ignore
        _TOPO_GEOMETRY_VECTORS,
        _get_ml_bond_length,
        _METAL_SET,
        _COVALENT_RADII,
    )
except Exception:  # pragma: no cover
    _TOPO_GEOMETRY_VECTORS = {}
    _METAL_SET = set()
    _COVALENT_RADII = {}

    def _get_ml_bond_length(metal_symbol: str, donor_symbol: str) -> float:  # type: ignore
        return 2.0


# ---------------------------------------------------------------------------
# Mass / vdW helpers — try the dedicated stub modules first; fall back inline.
# ---------------------------------------------------------------------------
try:
    from delfin._atom_weights import get_atom_weight as _ext_get_atom_weight  # type: ignore
    _HAVE_EXT_WEIGHTS = True
except Exception:
    _HAVE_EXT_WEIGHTS = False

try:
    from delfin._vdw_radii import get_vdw_sum as _ext_get_vdw_sum, VDW_RADII as _EXT_VDW  # type: ignore
    _HAVE_EXT_VDW = True
except Exception:
    _HAVE_EXT_VDW = False


# Inline atomic mass table (a.u.). Used when delfin._atom_weights is not present.
_ATOMIC_MASS: Dict[str, float] = {
    "H": 1.008, "He": 4.003, "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990,
    "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906,
    "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.906, "Pd": 106.42,
    "Ag": 107.868, "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.76,
    "Te": 127.60, "I": 126.904, "Cs": 132.905, "Ba": 137.327, "La": 138.905,
    "Ce": 140.116, "Pr": 140.908, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.925, "Dy": 162.500, "Ho": 164.930,
    "Er": 167.259, "Tm": 168.934, "Yb": 173.045, "Lu": 174.967, "Hf": 178.49,
    "Ta": 180.948, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
    "Pt": 195.084, "Au": 196.967, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.980, "Ac": 227.0, "Th": 232.038, "Pa": 231.036, "U": 238.029,
    "Np": 237.0, "Pu": 244.0,
}

# Inline Bondi vdW radii (Å). Used when delfin._vdw_radii is not present.
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


def get_atom_weight(mol, idx: int) -> float:
    """Atomic mass in a.u. (fallback to inline table)."""
    if _HAVE_EXT_WEIGHTS:
        try:
            return float(_ext_get_atom_weight(mol, idx))
        except Exception:
            pass
    try:
        sym = mol.GetAtomWithIdx(idx).GetSymbol()
    except Exception:
        return 12.0
    return float(_ATOMIC_MASS.get(sym, 12.0))


def get_vdw_sum(sym_i: str, sym_j: str) -> float:
    """Sum of Bondi vdW radii (fallback 1.70 Å per atom)."""
    if _HAVE_EXT_VDW:
        try:
            return float(_ext_get_vdw_sum(sym_i, sym_j))
        except Exception:
            pass
    ri = _VDW_RADII.get(sym_i, 1.70)
    rj = _VDW_RADII.get(sym_j, 1.70)
    return float(ri + rj)


def _cov_radius(sym: str) -> float:
    """Covalent radius (Å). Falls back to 1.50 for unknown symbols."""
    r = _COVALENT_RADII.get(sym)
    if r is None:
        r = 0.76 if sym == "C" else 1.50
    return float(r)


def _is_metal(sym: str) -> bool:
    return sym in _METAL_SET if _METAL_SET else False


# ---------------------------------------------------------------------------
# XYZ I/O — shape-preserving (matches _coord_angle_corrector format).
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
            pts.append(np.array([float(m.group(2)), float(m.group(3)),
                                 float(m.group(4))]))
    if not pts:
        return syms, np.zeros((0, 3)), lines
    return syms, np.vstack(pts), lines


def _format_xyz(orig_lines: List[str], syms: List[str],
                positions: np.ndarray) -> str:
    out: List[str] = []
    atom_i = 0
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
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Graph + fragment helpers
# ---------------------------------------------------------------------------

def _mol_neighbors(mol) -> List[List[int]]:
    """Adjacency list from RDKit Mol — bond-based, not distance-based."""
    n = mol.GetNumAtoms()
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        nbrs[i].append(j)
        nbrs[j].append(i)
    return nbrs


def _metal_indices(mol) -> List[int]:
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() in _METAL_SET]


def _donor_indices(mol, m_idx: int) -> List[int]:
    out: List[int] = []
    atom = mol.GetAtomWithIdx(m_idx)
    for nb in atom.GetNeighbors():
        if nb.GetSymbol() in _METAL_SET:
            continue
        if nb.GetSymbol() == "H":
            continue
        out.append(nb.GetIdx())
    return out


def _bfs_fragment(nbrs: List[List[int]], start: int,
                  blocked: Set[int]) -> Set[int]:
    """BFS reachable atoms from start, never traversing ``blocked``."""
    if start in blocked:
        return set()
    visited = {start}
    queue = deque([start])
    while queue:
        cur = queue.popleft()
        for nb in nbrs[cur]:
            if nb in blocked or nb in visited:
                continue
            visited.add(nb)
            queue.append(nb)
    return visited


def _decompose_into_fragments(mol, metals: Sequence[int]) -> List[Dict]:
    """For every donor D, compute the ligand fragment behind D (BFS from D
    with all metals blocked, but allowing traversal into any chelate co-donors
    naturally — i.e. only metals are blocked).

    Returns list of dicts with keys: ``atoms`` (set[int]), ``anchor_M`` (int),
    ``anchor_D`` (int).
    """
    nbrs = _mol_neighbors(mol)
    blocked = set(metals)
    frags: List[Dict] = []
    seen_anchors: Set[Tuple[int, int]] = set()
    for m_idx in metals:
        for d_idx in _donor_indices(mol, m_idx):
            key = (m_idx, d_idx)
            if key in seen_anchors:
                continue
            seen_anchors.add(key)
            atoms = _bfs_fragment(nbrs, d_idx, blocked)
            frags.append({
                "atoms": atoms,
                "anchor_M": m_idx,
                "anchor_D": d_idx,
            })
    return frags


# ---------------------------------------------------------------------------
# Rigid-body transformations (in-place on a returned copy)
# ---------------------------------------------------------------------------

def _fragment_translate(coords: np.ndarray, frag_atoms: Set[int],
                        translation: np.ndarray) -> np.ndarray:
    out = coords.copy()
    idx = np.fromiter(frag_atoms, dtype=int)
    if idx.size == 0:
        return out
    out[idx] += translation
    return out


def _rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues' rotation formula — works without SciPy too."""
    a = axis / max(float(np.linalg.norm(axis)), 1e-12)
    if _SCIPY_OK:
        return _Rotation.from_rotvec(a * angle_rad).as_matrix()
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    ux, uy, uz = a
    return np.array([
        [c + ux*ux*(1-c),     ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s],
        [uy*ux*(1-c) + uz*s,  c + uy*uy*(1-c),    uy*uz*(1-c) - ux*s],
        [uz*ux*(1-c) - uy*s,  uz*uy*(1-c) + ux*s, c + uz*uz*(1-c)],
    ])


def _fragment_rotate(coords: np.ndarray, frag_atoms: Set[int],
                     axis_point: np.ndarray, axis_dir: np.ndarray,
                     angle_rad: float) -> np.ndarray:
    out = coords.copy()
    if not frag_atoms or abs(angle_rad) < 1e-9:
        return out
    R = _rotation_matrix(axis_dir, angle_rad)
    idx = np.fromiter(frag_atoms, dtype=int)
    rel = out[idx] - axis_point
    out[idx] = (rel @ R.T) + axis_point
    return out


# ---------------------------------------------------------------------------
# Constraint primitives
# ---------------------------------------------------------------------------

def _vec_angle_deg(p_center: np.ndarray, p_a: np.ndarray,
                   p_b: np.ndarray) -> Optional[float]:
    v1 = p_a - p_center
    v2 = p_b - p_center
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return None
    cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(math.degrees(math.acos(cos)))


def _non_bonded_heavy_pairs(mol, max_pairs: int = 4000) -> List[Tuple[int, int]]:
    """Enumerate non-bonded heavy-heavy index pairs (i<j).  Capped for very
    large molecules to keep the topology gate O(N)."""
    n = mol.GetNumAtoms()
    bonded: Set[Tuple[int, int]] = set()
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonded.add((min(i, j), max(i, j)))
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pairs: List[Tuple[int, int]] = []
    for a_i in range(len(heavy)):
        for a_j in range(a_i + 1, len(heavy)):
            i = heavy[a_i]
            j = heavy[a_j]
            if (i, j) in bonded:
                continue
            pairs.append((i, j))
            if len(pairs) >= max_pairs:
                return pairs
    return pairs


def _md_pairs(mol, metals: Sequence[int]) -> List[Tuple[int, int, float]]:
    """List of (m_idx, d_idx, ideal_length)."""
    out: List[Tuple[int, int, float]] = []
    for m in metals:
        m_sym = mol.GetAtomWithIdx(m).GetSymbol()
        for d in _donor_indices(mol, m):
            d_sym = mol.GetAtomWithIdx(d).GetSymbol()
            try:
                d_ideal = float(_get_ml_bond_length(m_sym, d_sym))
            except Exception:
                d_ideal = 2.0
            out.append((m, d, d_ideal))
    return out


def _passes_topology(coords: np.ndarray, mol, metals: Sequence[int],
                     md_pairs: Sequence[Tuple[int, int, float]],
                     nb_pairs: Sequence[Tuple[int, int]]) -> bool:
    """Hard topology gate — see module docstring.

    v2 (post wave8-b5 forensik): tightened M-D window from [0.85, 1.10] →
    [0.93, 1.07] after observed Ir-N compression to 1.71Å (0.83×ideal) on
    02-Ir(ppy)2(acac). 15% compression chemically catastrophic.
    """
    # M-D distance window — tighter to prevent compression cascade.
    for (m, d, d_ideal) in md_pairs:
        d_cur = float(np.linalg.norm(coords[m] - coords[d]))
        if d_cur < 0.93 * d_ideal or d_cur > 1.07 * d_ideal:
            return False
    # Non-bonded heavy-heavy collapse → new spurious bond.
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    for (i, j) in nb_pairs:
        # Skip M-(heavy) non-bond pairs: M already covered above and counter-ion
        # / outer-sphere fragments are allowed close approach.  Skip pairs that
        # share a metal neighbour (chelate close-contacts).
        if syms[i] in _METAL_SET or syms[j] in _METAL_SET:
            continue
        r_sum = _cov_radius(syms[i]) + _cov_radius(syms[j])
        d_cur = float(np.linalg.norm(coords[i] - coords[j]))
        if d_cur < 0.85 * r_sum:
            return False
    return True


def _total_violation(coords: np.ndarray, mol, metals: Sequence[int],
                     md_pairs: Sequence[Tuple[int, int, float]],
                     nb_pairs: Sequence[Tuple[int, int]],
                     clash_factor: float) -> float:
    """Sum of squared violations across M-D bonds and non-bonded clashes."""
    s = 0.0
    for (m, d, d_ideal) in md_pairs:
        d_cur = float(np.linalg.norm(coords[m] - coords[d]))
        s += (d_cur - d_ideal) ** 2
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    for (i, j) in nb_pairs:
        d_cur = float(np.linalg.norm(coords[i] - coords[j]))
        thr = get_vdw_sum(syms[i], syms[j]) * clash_factor
        if d_cur < thr:
            s += (thr - d_cur) ** 2
    return float(s)


def _max_md_deviation(coords: np.ndarray,
                      md_pairs: Sequence[Tuple[int, int, float]]) -> float:
    worst = 0.0
    for (m, d, d_ideal) in md_pairs:
        d_cur = float(np.linalg.norm(coords[m] - coords[d]))
        worst = max(worst, abs(d_cur - d_ideal))
    return worst


def _count_clashes(coords: np.ndarray, mol,
                   nb_pairs: Sequence[Tuple[int, int]],
                   clash_factor: float) -> int:
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    n = 0
    for (i, j) in nb_pairs:
        d_cur = float(np.linalg.norm(coords[i] - coords[j]))
        thr = get_vdw_sum(syms[i], syms[j]) * clash_factor
        if d_cur < thr:
            n += 1
    return n


# ---------------------------------------------------------------------------
# Phase A — Topology repair (catastrophic M-D break)
# ---------------------------------------------------------------------------

def _find_catastrophic_md_breaks(
    coords: np.ndarray,
    mol,
    md_pairs: Sequence[Tuple[int, int, float]],
    break_factor: float = 1.30,
) -> List[Tuple[int, int, float, float]]:
    """Return [(m, d, d_cur, d_ideal), ...] for M-D pairs that are >1.30·ideal."""
    out: List[Tuple[int, int, float, float]] = []
    for (m, d, d_ideal) in md_pairs:
        d_cur = float(np.linalg.norm(coords[m] - coords[d]))
        if d_cur > break_factor * d_ideal:
            out.append((m, d, d_cur, d_ideal))
    return out


def _repair_md_break_translation(
    coords: np.ndarray,
    mol,
    nbrs: List[List[int]],
    metals: Sequence[int],
    m: int, d: int, d_ideal: float,
    nb_pairs: Sequence[Tuple[int, int]],
    md_pairs: Sequence[Tuple[int, int, float]],
    clash_factor: float,
) -> Optional[np.ndarray]:
    """Rigid-translate the fragment behind D so that |M-D| == d_ideal.  Reject
    if topology gate fails afterwards."""
    blocked = set(metals)
    frag = _bfs_fragment(nbrs, d, blocked)
    if not frag:
        return None
    cur = coords[d] - coords[m]
    cur_d = float(np.linalg.norm(cur))
    if cur_d < 1e-9:
        return None
    direction = cur / cur_d
    target = coords[m] + direction * d_ideal
    translation = target - coords[d]
    new_coords = _fragment_translate(coords, frag, translation)
    if _passes_topology(new_coords, mol, metals, md_pairs, nb_pairs):
        return new_coords
    return None


def _repair_md_break_rotation(
    coords: np.ndarray,
    mol,
    nbrs: List[List[int]],
    metals: Sequence[int],
    m: int, d: int, d_ideal: float,
    nb_pairs: Sequence[Tuple[int, int]],
    md_pairs: Sequence[Tuple[int, int, float]],
    clash_factor: float,
    n_angles: int = 6,
) -> Optional[np.ndarray]:
    """If pure translation would clash, try rotating the fragment around an
    axis through M perpendicular to the M→D direction (sweep a few angles)."""
    blocked = set(metals)
    frag = _bfs_fragment(nbrs, d, blocked)
    if not frag:
        return None
    direction = coords[d] - coords[m]
    cur_d = float(np.linalg.norm(direction))
    if cur_d < 1e-9:
        return None
    # Build a perpendicular axis (any vector orthogonal to direction).
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(direction, ref) / cur_d)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    axis = np.cross(direction, ref)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    best = None
    best_score = math.inf
    for k in range(1, n_angles + 1):
        for sign in (+1.0, -1.0):
            ang = sign * (k * math.pi / (n_angles + 1) / 2.0)
            tried = _fragment_rotate(coords, frag, coords[m], axis, ang)
            # After rotation, also pull along new M→D direction toward d_ideal.
            new_dir = tried[d] - coords[m]
            new_len = float(np.linalg.norm(new_dir))
            if new_len < 1e-9:
                continue
            unit = new_dir / new_len
            target = coords[m] + unit * d_ideal
            tried = _fragment_translate(tried, frag, target - tried[d])
            if not _passes_topology(tried, mol, metals, md_pairs, nb_pairs):
                continue
            score = _total_violation(tried, mol, metals, md_pairs, nb_pairs,
                                     clash_factor)
            if score < best_score:
                best_score = score
                best = tried
    return best


# ---------------------------------------------------------------------------
# Phase A.5 — Symmetry projection to ideal polyhedron
# ---------------------------------------------------------------------------

# Map CN → ordered list of geometry codes to try (preferred first).  These
# names match the keys of ``_TOPO_GEOMETRY_VECTORS``.
_CN_TO_GEOMS: Dict[int, List[str]] = {
    2: ["LIN"],
    3: ["TP", "TS"],
    4: ["SQ", "TH", "SS"],
    5: ["TBP", "SP"],
    6: ["OH", "TPR"],
    7: ["PBP", "COH"],
    8: ["SAP", "DD"],
    9: ["TTP"],
    10: ["BCSAP", "PAP"],
    11: ["CPAP"],
    12: ["ICOS", "CUBO", "HBP"],
}


def _best_polyhedron_for_metal(
    coords: np.ndarray, m_idx: int, donor_idxs: Sequence[int]
) -> Optional[Tuple[str, List[Tuple[float, float, float]]]]:
    """Pick the (geometry, vectors) that best matches the current donor layout
    using the Hungarian-assignment cost on unit vectors."""
    cn = len(donor_idxs)
    geom_keys = _CN_TO_GEOMS.get(cn, [])
    if not geom_keys or not _SCIPY_OK:
        return None
    # Current donor unit vectors (M → D).
    cur_units = []
    for d in donor_idxs:
        v = coords[d] - coords[m_idx]
        n = float(np.linalg.norm(v))
        if n < 1e-9:
            return None
        cur_units.append(v / n)
    cur_arr = np.vstack(cur_units)

    best: Optional[Tuple[str, List[Tuple[float, float, float]], float]] = None
    for key in geom_keys:
        vecs = _TOPO_GEOMETRY_VECTORS.get(key)
        if vecs is None or len(vecs) != cn:
            continue
        ideal = np.array(vecs, dtype=float)
        norms = np.linalg.norm(ideal, axis=1, keepdims=True)
        ideal_units = ideal / np.maximum(norms, 1e-12)
        cost = -(cur_arr @ ideal_units.T)
        try:
            row, col = linear_sum_assignment(cost)
        except Exception:
            continue
        total = float(cost[row, col].sum())  # more negative = better
        if best is None or total < best[2]:
            best = (key, [tuple(v) for v in vecs], total)
    if best is None:
        return None
    return best[0], best[1]


def _project_donors_to_polyhedron(
    coords: np.ndarray,
    mol,
    nbrs: List[List[int]],
    metals: Sequence[int],
    md_pairs: Sequence[Tuple[int, int, float]],
    nb_pairs: Sequence[Tuple[int, int]],
    step_size: float,
    clash_factor: float,
) -> Tuple[np.ndarray, int]:
    """For each metal, pick the matching ideal polyhedron, then translate each
    donor's BFS fragment toward its ideal position by ``step_size`` of the
    delta.  Rolls each fragment-move back if it breaks topology.

    Returns (new_coords, n_projections_committed).
    """
    if not _SCIPY_OK:
        return coords, 0
    out = coords.copy()
    n_done = 0
    blocked = set(metals)
    for m in metals:
        donor_idxs = _donor_indices(mol, m)
        if len(donor_idxs) < 2:
            continue
        pick = _best_polyhedron_for_metal(out, m, donor_idxs)
        if pick is None:
            continue
        _, ideal_vecs = pick

        ideal = np.array(ideal_vecs, dtype=float)
        ideal_units = ideal / np.maximum(
            np.linalg.norm(ideal, axis=1, keepdims=True), 1e-12)
        cur_units = []
        cur_lens = []
        for d in donor_idxs:
            v = out[d] - out[m]
            n = float(np.linalg.norm(v))
            cur_lens.append(n)
            cur_units.append(v / max(n, 1e-12))
        cur_arr = np.vstack(cur_units)
        cost = -(cur_arr @ ideal_units.T)
        try:
            row, col = linear_sum_assignment(cost)
        except Exception:
            continue
        for r, c in zip(row, col):
            d = donor_idxs[r]
            target_unit = ideal_units[c]
            target_pos = out[m] + target_unit * cur_lens[r]
            delta = target_pos - out[d]
            shift = delta * float(step_size)
            frag = _bfs_fragment(nbrs, d, blocked)
            if not frag:
                continue
            tried = _fragment_translate(out, frag, shift)
            if _passes_topology(tried, mol, metals, md_pairs, nb_pairs):
                out = tried
                n_done += 1
    return out, n_done


# ---------------------------------------------------------------------------
# Phase B — Iterative Gauss-Seidel sweeps
# ---------------------------------------------------------------------------

def _stage_bonds(coords: np.ndarray, mol,
                 md_pairs: Sequence[Tuple[int, int, float]],
                 metals: Sequence[int], nb_pairs: Sequence[Tuple[int, int]],
                 step_size: float, bond_tol: float) -> Tuple[np.ndarray, bool, int]:
    """Mass-weighted bond-distance projection over all M-D pairs (and any
    further covalent bonds that are far from ideal).  Returns (coords, moved,
    n_repairs)."""
    out = coords.copy()
    moved = False
    repairs = 0
    # Resolve in order of largest violation first.
    ranked = []
    for (m, d, d_ideal) in md_pairs:
        d_cur = float(np.linalg.norm(out[m] - out[d]))
        err = d_cur - d_ideal
        if abs(err) > bond_tol * d_ideal:
            ranked.append((abs(err), m, d, d_ideal, d_cur, err))
    ranked.sort(reverse=True)
    for (_, m, d, d_ideal, d_cur, err) in ranked:
        w_i = get_atom_weight(mol, m)
        w_j = get_atom_weight(mol, d)
        if w_i + w_j < 1e-9:
            continue
        a_i = w_i / (w_i + w_j)
        a_j = w_j / (w_i + w_j)
        if d_cur < 1e-9:
            continue
        unit = (out[d] - out[m]) / d_cur
        shift = err * step_size
        new = out.copy()
        new[m] = out[m] + unit * shift * a_j
        new[d] = out[d] - unit * shift * a_i
        if _passes_topology(new, mol, metals, md_pairs, nb_pairs):
            out = new
            moved = True
            repairs += 1
    return out, moved, repairs


def _stage_angles(coords: np.ndarray, mol, nbrs: List[List[int]],
                  metals: Sequence[int],
                  md_pairs: Sequence[Tuple[int, int, float]],
                  nb_pairs: Sequence[Tuple[int, int]],
                  step_size: float, angle_tol: float
                  ) -> Tuple[np.ndarray, bool, int]:
    """Rotate non-metal donor X around (M→D) axis through D to bring M-D-X
    closer to a hybridization-dependent ideal.  Conservative: only acts on
    light X atoms (single neighbour, e.g. H) or terminal heavy atoms to avoid
    fighting the angle-corrector module."""
    out = coords.copy()
    moved = False
    repairs = 0
    for (m, d, _) in md_pairs:
        d_atom = mol.GetAtomWithIdx(d)
        d_nbrs = [n for n in nbrs[d] if n != m]
        if not d_nbrs:
            continue
        # Hybridization heuristic via RDKit when available.
        try:
            hyb_label = str(d_atom.GetHybridization())
        except Exception:
            hyb_label = ""
        if "SP3" in hyb_label:
            ideal = 109.47
        elif "SP2" in hyb_label or d_atom.GetIsAromatic():
            ideal = 120.0
        elif "SP" == hyb_label or "SP " in hyb_label or "SP," in hyb_label:
            ideal = 180.0
        else:
            ideal = 120.0
        for x in d_nbrs:
            ang = _vec_angle_deg(out[d], out[m], out[x])
            if ang is None:
                continue
            dev = ang - ideal
            if abs(dev) < angle_tol:
                continue
            # Only act if x's subtree is small (terminal-ish) — keeps the
            # operation safe and complementary to B3.
            blocked = set(metals) | {d}
            subtree = _bfs_fragment(nbrs, x, blocked)
            if len(subtree) > 6:
                continue
            v_md = out[m] - out[d]
            v_dx = out[x] - out[d]
            axis = np.cross(v_md, v_dx)
            if float(np.linalg.norm(axis)) < 1e-6:
                continue
            angle_rad = -math.radians(dev) * step_size
            new = _fragment_rotate(out, subtree, out[d], axis, angle_rad)
            if _passes_topology(new, mol, metals, md_pairs, nb_pairs):
                out = new
                moved = True
                repairs += 1
    return out, moved, repairs


def _stage_clashes(coords: np.ndarray, mol, nbrs: List[List[int]],
                   metals: Sequence[int],
                   md_pairs: Sequence[Tuple[int, int, float]],
                   nb_pairs: Sequence[Tuple[int, int]],
                   step_size: float, clash_factor: float
                   ) -> Tuple[np.ndarray, bool, int]:
    """Resolve non-bonded clashes using a 3-strategy ladder.

    Strategy 1: rotate one of the offending atoms' fragment around its M-D
    axis (good for inter-ligand clash between two donor-side fragments).
    Strategy 2: mass-weighted pair-push along the pair direction.
    Strategy 3: small perpendicular fragment translation.
    """
    out = coords.copy()
    moved = False
    repairs = 0
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    blocked = set(metals)

    # Index atom → enclosing fragment (anchor M, D, atoms set).
    atom_to_frag: Dict[int, Tuple[int, int, Set[int]]] = {}
    for m in metals:
        for d in _donor_indices(mol, m):
            frag = _bfs_fragment(nbrs, d, blocked)
            for a in frag:
                atom_to_frag.setdefault(a, (m, d, frag))

    # Rank clashes by overlap depth.
    ranked = []
    for (i, j) in nb_pairs:
        d_cur = float(np.linalg.norm(out[i] - out[j]))
        thr = get_vdw_sum(syms[i], syms[j]) * clash_factor
        if d_cur < thr:
            ranked.append((thr - d_cur, i, j, d_cur, thr))
    ranked.sort(reverse=True)

    for (overlap, i, j, d_cur, thr) in ranked:
        if d_cur < 1e-9:
            continue
        unit = (out[j] - out[i]) / d_cur
        committed = False

        # Strategy 1 — fragment rotation around M-D axis if both atoms live
        # in different ligand fragments under the same metal.
        fi = atom_to_frag.get(i)
        fj = atom_to_frag.get(j)
        if fi is not None and fj is not None and fi[0] == fj[0] and fi[1] != fj[1]:
            for victim in (fi, fj):
                m_a, d_a, atoms_a = victim
                axis_dir = out[d_a] - out[m_a]
                axis_n = float(np.linalg.norm(axis_dir))
                if axis_n < 1e-9:
                    continue
                axis_dir = axis_dir / axis_n
                # Try a few small rotation magnitudes.
                for ang_deg in (4.0, -4.0, 8.0, -8.0):
                    ang_rad = math.radians(ang_deg) * step_size
                    new = _fragment_rotate(out, atoms_a, out[m_a],
                                           axis_dir, ang_rad)
                    if not _passes_topology(new, mol, metals, md_pairs, nb_pairs):
                        continue
                    new_d = float(np.linalg.norm(new[i] - new[j]))
                    if new_d > d_cur + 1e-3:
                        out = new
                        moved = True
                        repairs += 1
                        committed = True
                        break
                if committed:
                    break
        if committed:
            continue

        # Strategy 2 — mass-weighted pair push apart (purely along unit).
        w_i = get_atom_weight(mol, i)
        w_j = get_atom_weight(mol, j)
        if w_i + w_j < 1e-9:
            continue
        a_i = w_i / (w_i + w_j)
        a_j = w_j / (w_i + w_j)
        push = overlap * step_size
        new = out.copy()
        new[i] = out[i] - unit * push * a_j
        new[j] = out[j] + unit * push * a_i
        if _passes_topology(new, mol, metals, md_pairs, nb_pairs):
            out = new
            moved = True
            repairs += 1
            continue

        # Strategy 3 — perpendicular fragment shift of the lighter side.
        target = j if w_j <= w_i else i
        f = atom_to_frag.get(target)
        if f is None:
            continue
        _, _, frag_atoms = f
        # Perpendicular vector to the pair axis, in the plane with z.
        ref = np.array([0.0, 0.0, 1.0])
        perp = np.cross(unit, ref)
        if float(np.linalg.norm(perp)) < 1e-6:
            perp = np.cross(unit, np.array([1.0, 0.0, 0.0]))
        perp = perp / max(float(np.linalg.norm(perp)), 1e-12)
        shift = perp * (overlap * step_size * 0.5)
        new = _fragment_translate(out, frag_atoms, shift)
        if _passes_topology(new, mol, metals, md_pairs, nb_pairs):
            out = new
            moved = True
            repairs += 1
    return out, moved, repairs


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def post_optimize_geometry(
    xyz: str,
    mol,
    class_label: str = "sigma",
    max_iter: int = 10,
    step_size: float = 0.1,
    bond_tol: float = 0.15,
    angle_tol: float = 10.0,
    clash_factor: float = 0.85,
    enable_symmetry: bool = False,
    enable_angles: bool = True,
    md_drift_max: float = 0.05,
) -> Tuple[str, Dict]:
    """Apply PBD post-optimization to ``xyz`` using bonding info from ``mol``.

    Parameters
    ----------
    xyz : str
        Input XYZ-formatted string (with N + comment header).
    mol : rdkit.Chem.Mol
        Molecule with a single conformer matching ``xyz`` atom-for-atom.
    class_label : str
        Coordination class hint (``sigma``, ``hapto``, ``multi_sigma``,
        ``multi_hapto``, ``no_metal``).  Currently only used to disable
        symmetry projection for ``hapto`` / ``no_metal``.
    max_iter : int
        Cap on Gauss-Seidel sweeps.
    step_size : float
        Relaxation factor for each projection (0..1).
    bond_tol : float
        Relative tolerance for M-D bond deviation (fraction of ideal length).
    angle_tol : float
        Absolute tolerance (degrees) before an M-D-X angle is corrected.
    clash_factor : float
        Multiplier on Σr_vdw used as the clash threshold.
    enable_symmetry : bool
        Run Phase A.5 polyhedron projection.
    enable_angles : bool
        Run Stage 2 angle corrections inside the iterative loop.

    Returns
    -------
    (new_xyz, report)
        ``new_xyz`` is a shape-identical XYZ string.  ``report`` keys:
        ``iterations``, ``repairs``, ``final_md_max_dev``,
        ``final_clash_count``, ``converged``, ``topology_preserved``.
    """
    report: Dict = {
        "iterations": 0,
        "repairs": 0,
        "final_md_max_dev": 0.0,
        "final_clash_count": 0,
        "converged": False,
        "topology_preserved": True,
    }
    if mol is None:
        return xyz, report
    try:
        syms, coords, orig_lines = _parse_xyz(xyz)
        if coords.shape[0] != mol.GetNumAtoms():
            # Atom count mismatch — fall back gracefully.
            return xyz, report

        metals = _metal_indices(mol)
        nbrs = _mol_neighbors(mol)
        md_pairs = _md_pairs(mol, metals)
        nb_pairs = _non_bonded_heavy_pairs(mol)

        # If no metals, no constraints to enforce — return unchanged.
        if not metals:
            report["converged"] = True
            return xyz, report

        # Verify initial topology — if it is already broken, do not commit
        # any moves that would propagate failure.  We still run the repair
        # phase, which can restore the gate.
        initial_topology_ok = _passes_topology(
            coords, mol, metals, md_pairs, nb_pairs)

        total_repairs = 0

        # ---- Phase A: catastrophic M-D break repair ------------------------
        breaks = _find_catastrophic_md_breaks(coords, mol, md_pairs)
        for (m, d, d_cur, d_ideal) in breaks:
            new = _repair_md_break_translation(
                coords, mol, nbrs, metals, m, d, d_ideal,
                nb_pairs, md_pairs, clash_factor)
            if new is None:
                new = _repair_md_break_rotation(
                    coords, mol, nbrs, metals, m, d, d_ideal,
                    nb_pairs, md_pairs, clash_factor)
            if new is not None:
                coords = new
                total_repairs += 1

        # ---- Phase A.5: symmetry projection --------------------------------
        if enable_symmetry and class_label not in ("hapto", "multi_hapto",
                                                    "no_metal"):
            coords, n_sym = _project_donors_to_polyhedron(
                coords, mol, nbrs, metals, md_pairs, nb_pairs,
                step_size=step_size, clash_factor=clash_factor)
            total_repairs += n_sym

        # ---- Phase B: Gauss-Seidel sweeps ----------------------------------
        cur_step = float(step_size)
        prev_total = _total_violation(
            coords, mol, metals, md_pairs, nb_pairs, clash_factor)

        converged = False
        last_iter = 0
        for it in range(max_iter):
            last_iter = it + 1
            moved_any = False
            snapshot = coords.copy()

            new, m1, r1 = _stage_bonds(
                coords, mol, md_pairs, metals, nb_pairs,
                cur_step, bond_tol)
            coords = new
            moved_any = moved_any or m1
            total_repairs += r1

            if enable_angles:
                new, m2, r2 = _stage_angles(
                    coords, mol, nbrs, metals, md_pairs, nb_pairs,
                    cur_step, angle_tol)
                coords = new
                moved_any = moved_any or m2
                total_repairs += r2

            new, m3, r3 = _stage_clashes(
                coords, mol, nbrs, metals, md_pairs, nb_pairs,
                cur_step, clash_factor)
            coords = new
            moved_any = moved_any or m3
            total_repairs += r3

            cur_total = _total_violation(
                coords, mol, metals, md_pairs, nb_pairs, clash_factor)
            if prev_total > 1e-9 and cur_total > prev_total * 1.05:
                # Rolled back this sweep — halve step and restart from snapshot.
                coords = snapshot
                cur_step *= 0.5
                if cur_step < 0.02:
                    break
                continue
            prev_total = cur_total

            if not moved_any:
                converged = True
                break

        # ---- Finalise ------------------------------------------------------
        topology_ok = _passes_topology(coords, mol, metals, md_pairs, nb_pairs)
        if not topology_ok and initial_topology_ok:
            # Worst-case: we broke the gate while the input was fine.  Fall
            # back to the input coordinates so the caller never gets a worse
            # structure than what they passed in.
            _, coords, _ = _parse_xyz(xyz)
            topology_ok = True
            converged = False

        # ---- v2 quality preservation gate (post wave8-b5 forensik) ---------
        # Reject result if ANY M-D bond drifted more than md_drift_max from
        # input, UNLESS the input itself was broken (Phase A catastrophic
        # repair is supposed to drift large to fix breaks). Triggers only
        # when input was already topology-OK and we drifted bonds anyway.
        if topology_ok and md_pairs and initial_topology_ok:
            _, coords_input, _ = _parse_xyz(xyz)
            drift_max = 0.0
            for (m_idx, d_idx, _d_ideal) in md_pairs:
                d_before = float(np.linalg.norm(coords_input[m_idx] - coords_input[d_idx]))
                d_after = float(np.linalg.norm(coords[m_idx] - coords[d_idx]))
                drift_max = max(drift_max, abs(d_after - d_before))
            if drift_max > md_drift_max:
                # Optimization caused excessive M-D drift on already-good
                # input — fall back to input (Bundle-1-equivalent behavior).
                coords = coords_input
                converged = False
                report["fallback_md_drift"] = drift_max

        report["iterations"] = last_iter
        report["repairs"] = total_repairs
        report["final_md_max_dev"] = _max_md_deviation(coords, md_pairs)
        report["final_clash_count"] = _count_clashes(
            coords, mol, nb_pairs, clash_factor)
        report["final_violations"] = float(_total_violation(
            coords, mol, metals, md_pairs, nb_pairs, clash_factor))
        report["converged"] = bool(converged)
        report["topology_preserved"] = bool(topology_ok)
        new_xyz = _format_xyz(orig_lines, syms, coords)
        return new_xyz, report
    except Exception as exc:  # pragma: no cover - defensive
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc(limit=4)
        return xyz, report


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _selftest() -> None:  # pragma: no cover
    """Run on a handful of synthetic structures.  Requires RDKit."""
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
    except Exception as exc:
        print("[selftest] RDKit not available:", exc)
        return

    samples = [
        ("acetate Pt(II)", "[Pt]([NH3])([NH3])(Cl)Cl"),
        ("ferrocene-like", "[Fe]"),  # no donors — should be no-op
        ("Cu pyridine", "[Cu](C1=CC=NC=C1)(C2=CC=NC=C2)(Cl)Cl"),
    ]
    for label, smi in samples:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            print(f"[{label}] mol parse failed for {smi!r}")
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        except Exception:
            pass
        if mol.GetNumConformers() == 0:
            print(f"[{label}] embed failed; skipping")
            continue
        conf = mol.GetConformer(0)
        lines = [str(mol.GetNumAtoms()), label]
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            lines.append(
                f"{mol.GetAtomWithIdx(i).GetSymbol():4s} {p.x:12.6f} "
                f"{p.y:12.6f} {p.z:12.6f}"
            )
        xyz = "\n".join(lines) + "\n"
        new_xyz, report = post_optimize_geometry(
            xyz, mol, class_label="sigma", max_iter=10, step_size=0.3)
        print(f"[{label}] report = {report}")


if __name__ == "__main__":  # pragma: no cover
    _selftest()
