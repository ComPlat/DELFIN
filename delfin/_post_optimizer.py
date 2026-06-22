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
import os
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


def _compute_h_neighbors(mol) -> List[List[int]]:
    """Per-atom list of bonded H atom indices.

    Used as a rigid-H tracking aid for single-atom moves in
    :func:`_stage_bonds` (M, D shifts) and :func:`_stage_clashes` Strategy 2
    (clash pair push). Dragging bonded H atoms with their heavy parent
    by the same delta preserves C-H / N-H / O-H bond lengths that the
    heavy-atom-only topology gate does not check.

    Returns a list of length ``mol.GetNumAtoms()``; entries indexed by an H
    atom are always empty so the helper is safe to use on any atom index.
    """
    n = mol.GetNumAtoms()
    out: List[List[int]] = [[] for _ in range(n)]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        idx = atom.GetIdx()
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "H":
                out[idx].append(nb.GetIdx())
    return out


def _has_bridging_bh_hydride(mol) -> bool:
    """Detect a μ²-B-H-M bridging hydride (Patch P-BH-NOGO precheck).

    Returns True iff the molecule contains at least one hydrogen atom that is
    simultaneously bonded to a boron atom (atomic number 5) AND to a metal
    centre (symbol in :data:`_METAL_SET`). This is the topological signature
    of bridging-hydride borohydrides (e.g. Zr[(μ-H)₂BR₂]₄, scorpionate
    tris(pyrazolyl)borate-Au, Re-cyclooctyl-borate) where rigid-H tracking
    breaks the bridge: the H sits between B and M, and dragging it with one
    parent during a clash push catastrophically detaches it from the other.

    The check is purely graph-based — no SMILES strings, no element-list
    shortcuts beyond the universal ``_METAL_SET`` and the boron atomic-number
    comparison. Returns False for:

    * Molecules without boron.
    * Borohydrides isolated from a metal (B-H bonds present, but the H is not
      simultaneously bonded to any metal).
    * Boron compounds where every B-H neighbour is not also M-bonded.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule whose bond topology will be inspected. Must include explicit
        H atoms with the relevant B-H and M-H bonds present.

    Returns
    -------
    bool
        True if a μ²-B-H-M bridge is present; False otherwise.
    """
    if mol is None:
        return False
    try:
        # Pre-collect metal indices once.
        metal_idxs: Set[int] = {
            a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET
        }
        if not metal_idxs:
            return False
        for atom in mol.GetAtoms():
            # Atomic number 5 == Boron. Use atomic number (not symbol) so
            # detection is robust to symbol-case quirks in input files.
            if atom.GetAtomicNum() != 5:
                continue
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() != "H":
                    continue
                # H is bonded to this B. Check if it is also bonded to any
                # metal — that is the μ²-bridge signature.
                for nb2 in nb.GetNeighbors():
                    if nb2.GetIdx() in metal_idxs:
                        return True
    except Exception:
        # Any RDKit oddity must not crash the optimizer; fall back to "no
        # bridge detected" so default behaviour is unchanged.
        return False
    return False


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


def _ring_aware_subtree(
    mol,
    nbrs: List[List[int]],
    x: int,
    blocked: Set[int],
) -> Set[int]:
    """BFS subtree from ``x`` expanded to ring-consistent membership.

    Behaviour:

    - Compute the standard :func:`_bfs_fragment` subtree from ``x`` with
      ``blocked``.
    - If ``x`` belongs to one or more SSSR rings (RDKit ``GetRingInfo``),
      union the subtree with **every atom of every such ring** plus each
      ring atom's bonded H children.  Atoms already in ``blocked`` are
      excluded from the union so a ring closing through the donor ``d``
      does not pull ``d`` into the rotated set.

    Rationale (Patch P-RING-AWARE):  the legacy size-6 BFS subtree gate in
    :func:`_stage_angles` rotates only a fraction of an aromatic/aliphatic
    ring when the ring is small enough to fit under the cap, leaving the
    remaining ring atoms and their bonded H's stationary.  The partial
    rotation stretches the ring's heavy backbone *and* the C-H bonds on
    the unrotated ring atoms; in hapto-CN7 systems with η⁵-Cp or η⁶-arene
    substituent rings the result is one or more ORPHAN_H (nearest-heavy
    > 1.40 Å).  Including the full ring + ring-bonded H's keeps the ring
    rigid under the rotation so C-H bonds and ring-internal heavy bonds
    are preserved.

    Pure RDKit graph-based — no SMILES patterns, no refcode shortcuts.
    Safe on molecules without rings (returns the plain BFS subtree) and
    on RDKit oddities (any exception falls back to the plain subtree).
    """
    base = _bfs_fragment(nbrs, x, blocked)
    if not base or mol is None:
        return base
    try:
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings() if ring_info is not None else ()
    except Exception:
        return base
    if not atom_rings:
        return base
    expanded: Set[int] = set(base)
    for ring in atom_rings:
        if x not in ring:
            continue
        for ridx in ring:
            if ridx in blocked:
                continue
            expanded.add(ridx)
            # H children of the ring atom (bonded H only, via nbrs).
            for nb in nbrs[ridx]:
                if nb in blocked:
                    continue
                try:
                    if mol.GetAtomWithIdx(nb).GetSymbol() == "H":
                        expanded.add(nb)
                except Exception:
                    continue
    return expanded


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


# ---------------------------------------------------------------------------
# Smart Phase A.5 — Jahn-Teller / multiple-bond / cyclometal aware
# ---------------------------------------------------------------------------

# Common Jahn-Teller (E_g term ground-state) d^n configurations that exhibit
# strong tetragonal distortion in octahedral coordination. Cu(II) d^9 is the
# canonical "always-distorted" case; Mn(III) d^4 high-spin and Cr(II) d^4
# high-spin also show pronounced JT elongation. Conservative list — we only
# protect cases where the distortion is essentially guaranteed.
#
# Format: (element symbol, formal charge) → label
_JAHN_TELLER_CASES: Dict[Tuple[str, int], str] = {
    ("Cu", 2): "jahn_teller_cu2",   # d^9, the textbook JT system
    ("Mn", 3): "jahn_teller_mn3",   # d^4 HS, eg^1 → axial elongation
    ("Cr", 2): "jahn_teller_cr2",   # d^4 HS
    ("Ag", 2): "jahn_teller_ag2",   # d^9 (rare but unambiguous)
    ("Ni", 3): "jahn_teller_ni3",   # d^7 LS, t2g^6 eg^1
}

# Metal-oxo / metal-imido bond-length threshold (Å). Bonds shorter than this
# between a metal and a terminal O / N donor strongly imply M=X multiple bond
# character (typical M=O is 1.6-1.8 Å for early transition metals; M-O single
# bond is 1.9-2.2 Å). 1.85 Å is a conservative discriminator.
_METAL_OXO_BOND_MAX: float = 1.85

# Cyclometalation ring-size that we protect. Aryl C-H activation forms an
# essentially-rigid 5-membered M-C-C-C-X ring (X = donor heteroatom). We do
# NOT try to regularise such rings to standard polyhedra — they are part of
# the rigid backbone, not a free coordination sphere.
_CYCLOMETAL_RING_SIZE: int = 5


def _classify_distortion_mandatory(
    mol,
    coords: np.ndarray,
    metal_idx: int,
    donor_idxs: Sequence[int],
) -> Optional[str]:
    """Detect whether a metal centre carries a chemically-mandatory distortion.

    Returns one of the canonical distortion labels — or ``None`` when the
    metal coordination is "free" to be projected to a standard polyhedron.
    Pure graph + geometry detection: no SMILES strings, no refcode shortcuts.

    Detected patterns:

    * ``"jahn_teller_cu2"`` etc. — metal+formal-charge combination drawn from
      :data:`_JAHN_TELLER_CASES`. Octahedral-or-higher CN required (CN >= 5
      to cover 4+1, 4+2 and full Oh). The label is suffixed by the metal so
      callers can act per-element if needed.
    * ``"metal_oxo"`` — at least one terminal O or N donor lies within
      :data:`_METAL_OXO_BOND_MAX` of the metal (signature of M=O / M=NR
      multiple bond, e.g. vanadyl V=O, oxo-tungsten W=O, dioxo-osmium
      cis-OsO2). Terminal = the donor has no further heavy-atom neighbours
      except the metal (or it is bonded by a formally double bond).
    * ``"cyclometal_5ring"`` — at least one donor participates in a
      5-membered ring containing the metal AND an sp2/aromatic carbon
      bonded to the metal (the classic cyclometalation signature, e.g.
      Ir(ppy)3 phenyl-pyridyl ring closure).
    * ``"mer_tridentate"`` — three of the donors are pair-wise nearly-linear
      (one pair at ~180°, two pairs at ~90°) which is the meridional
      signature of a planar-tridentate ligand (e.g. terpyridine,
      pincer-NCN). Conservative cosine threshold (-0.90) so genuine
      planar-tridentates are caught but bent fac-arrangements are not.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule with bond topology consistent with ``coords``.
    coords : np.ndarray
        ``(n, 3)`` Cartesian coordinates indexed by atom idx.
    metal_idx : int
        Atom index of the central metal.
    donor_idxs : Sequence[int]
        Atom indices of donor atoms directly bonded to ``metal_idx``.

    Returns
    -------
    Optional[str]
        Distortion label if any pattern matches; ``None`` otherwise. On any
        RDKit exception, returns ``None`` (fail-safe → standard projection).
    """
    if mol is None or len(donor_idxs) == 0:
        return None
    try:
        m_atom = mol.GetAtomWithIdx(int(metal_idx))
    except Exception:
        return None
    m_sym = m_atom.GetSymbol()
    try:
        m_charge = int(m_atom.GetFormalCharge())
    except Exception:
        m_charge = 0
    cn = len(donor_idxs)

    # ----- Pattern 1: Jahn-Teller-distorted d^n configuration -----
    # Trigger when (metal, formal_charge) is in the JT list AND CN >= 5
    # (the 4+1/4+2 elongated octahedron is the most chemically frequent JT
    # case; CN < 5 cannot show axial elongation by construction). JT is the
    # most specific test (metal+charge match) and is therefore checked first.
    #
    # Welle-3 T2.2 extension (env-gated, default OFF): when
    # DELFIN_JT_CN4_DISTORTION=1, the gate is relaxed to CN >= 4 so
    # Cu(II)/Mn(III)/Cr(II)/Ag(II)/Ni(III) at CN=4 are also marked
    # distortion-mandatory.  Square-planar is the natural d9/d4-HS CN=4
    # polyhedron and must not be projected to tetrahedral by the
    # standard polyhedron picker.
    import os as _os_pp
    jt_label = _JAHN_TELLER_CASES.get((m_sym, m_charge))
    _jt_min_cn = 4 if (
        _os_pp.environ.get("DELFIN_JT_CN4_DISTORTION", "0")
        in ("1", "true", "True")
    ) else 5
    if jt_label is not None and cn >= _jt_min_cn:
        return jt_label

    # ----- Pattern 2: Cyclometalation 5-ring (M-C-...-X-M closure) -----
    # Checked before metal_oxo because a cyclometal chelate often forces the
    # heteroatom donor close to the metal (~1.95-2.05 Å) and the M-C may
    # also be short (~2.0 Å). Connectivity (ring closure) is the load-
    # bearing chemical signal; we look for a 5-membered ring containing the
    # metal where M is bonded to BOTH an aromatic / sp2 carbon AND another
    # heteroatom donor in the same ring.
    try:
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        if not rings:
            # RingInfo not initialised yet (common for RWMol built from
            # scratch without sanitisation). FastFindRings is the safe
            # non-sanitising option that populates the ring perception.
            try:
                from rdkit.Chem import FastFindRings  # type: ignore
                FastFindRings(mol)
                rings = mol.GetRingInfo().AtomRings()
            except Exception:
                rings = ()
    except Exception:
        rings = ()
    metal_int = int(metal_idx)
    donor_set = set(int(x) for x in donor_idxs)
    for ring in rings:
        if metal_int not in ring:
            continue
        if len(ring) != _CYCLOMETAL_RING_SIZE:
            continue
        donors_in_ring = donor_set.intersection(ring)
        if len(donors_in_ring) < 2:
            continue
        has_aryl_c = False
        has_hetero = False
        for d in donors_in_ring:
            try:
                d_atom = mol.GetAtomWithIdx(int(d))
            except Exception:
                continue
            d_sym = d_atom.GetSymbol()
            if d_sym == "C":
                try:
                    if d_atom.GetIsAromatic() or "SP2" in str(
                        d_atom.GetHybridization()
                    ):
                        has_aryl_c = True
                except Exception:
                    has_aryl_c = True   # be permissive on RDKit hiccups
                # Welle-3 T2.2 extension (env-gated, default OFF):
                # graph-only fallback when RDKit aromatic/hybridisation
                # perception was skipped (RWMol dative-bond-first builders).
                # If the carbon donor has >= 1 carbon ring-neighbour in
                # the same 5-ring (M-C-C-...-X-M closure), treat it as
                # cyclometal aryl.
                if (
                    not has_aryl_c
                    and _os_pp.environ.get(
                        "DELFIN_CYCLOMETAL_GRAPH_FALLBACK", "0"
                    ) in ("1", "true", "True")
                ):
                    try:
                        ring_set = set(int(x) for x in ring)
                        c_ring_nbrs = [
                            nb for nb in d_atom.GetNeighbors()
                            if int(nb.GetIdx()) in ring_set
                            and nb.GetSymbol() == "C"
                        ]
                        if len(c_ring_nbrs) >= 1:
                            has_aryl_c = True
                    except Exception:
                        pass
            elif d_sym in ("N", "O", "P", "S"):
                has_hetero = True
        if has_aryl_c and has_hetero:
            return "cyclometal_5ring"

    # ----- Pattern 3: Metal-oxo / metal-imido multiple bond -----
    # Any donor with element O or N at distance < 1.85 Å from the metal that
    # is also "terminal" (only the metal as heavy neighbour) or formally
    # double-bonded to the metal. Bond order check is RDKit-tolerant — many
    # M=O bonds are parsed as SINGLE due to dative-bond conventions, so we
    # fall back to a geometric criterion. We additionally require the donor
    # to be terminal in the heavy-atom subgraph (only the metal as heavy
    # neighbour) OR the bond is formally double-bonded — this excludes
    # short chelated donors from cyclometalation rings (already caught by
    # Pattern 2 above when applicable, otherwise harmless to skip metal_oxo).
    for d_idx in donor_idxs:
        try:
            d_atom = mol.GetAtomWithIdx(int(d_idx))
        except Exception:
            continue
        d_sym = d_atom.GetSymbol()
        if d_sym not in ("O", "N"):
            continue
        d_cur = float(np.linalg.norm(coords[int(metal_idx)] - coords[int(d_idx)]))
        # Heavy-atom neighbour count (excluding the metal). Terminal donors
        # have zero such neighbours; chelated ones have one or more.
        try:
            heavy_nbr_n = sum(
                1 for nb in d_atom.GetNeighbors()
                if nb.GetIdx() != metal_int and nb.GetAtomicNum() > 1
            )
        except Exception:
            heavy_nbr_n = -1
        # Bond-order criterion (catches explicit DOUBLE bonds).
        bond_is_double = False
        try:
            bond = mol.GetBondBetweenAtoms(int(metal_idx), int(d_idx))
            if bond is not None and float(bond.GetBondTypeAsDouble()) >= 1.5:
                bond_is_double = True
        except Exception:
            bond_is_double = False
        if bond_is_double:
            return "metal_oxo"
        # Geometric short-bond criterion: only fire when donor is terminal
        # (no other heavy neighbours) to avoid false-positives on chelated
        # heteroatoms that are pulled close by ring strain.
        if d_cur < _METAL_OXO_BOND_MAX and heavy_nbr_n == 0:
            return "metal_oxo"

    # ----- Pattern 4: Meridional planar tridentate -----
    # Three donors that are roughly co-linear-and-perpendicular: pairwise
    # cosines should be (-1, 0, 0) in the ideal mer-tridentate arrangement.
    # Restricted to CN == 3 (an isolated planar-tridentate ligand without
    # spectators) because in higher CN the pattern can be matched by any
    # three orthogonal donors of an octahedron — which we DO want projected.
    if cn == 3:
        units: List[np.ndarray] = []
        for d_idx in donor_idxs:
            v = coords[int(d_idx)] - coords[int(metal_idx)]
            n = float(np.linalg.norm(v))
            if n < 1e-9:
                continue
            units.append(v / n)
        n_units = len(units)
        for ii in range(n_units):
            for jj in range(ii + 1, n_units):
                cij = float(np.dot(units[ii], units[jj]))
                if cij > -0.90:        # need near-trans pair
                    continue
                for kk in range(n_units):
                    if kk in (ii, jj):
                        continue
                    cik = float(np.dot(units[ii], units[kk]))
                    cjk = float(np.dot(units[jj], units[kk]))
                    # k should be near-perpendicular to BOTH i and j
                    if abs(cik) < 0.30 and abs(cjk) < 0.30:
                        return "mer_tridentate"
    # CN >= 4: ALSO detect mer when the three near-collinear donors are
    # graph-connected via a single contiguous chelate chain (i.e., share a
    # ligand path that does not pass through the metal). This catches
    # tridentate-N3 + 3 spectator donors at CN=6 without false-positives
    # on free monodentate octahedral sets.
    if cn >= 4:
        mer_label = _try_detect_chelate_mer_tridentate(
            mol, coords, metal_int, list(donor_idxs)
        )
        if mer_label is not None:
            return mer_label

    return None


def _try_detect_chelate_mer_tridentate(
    mol,
    coords: np.ndarray,
    metal_idx: int,
    donor_idxs: List[int],
) -> Optional[str]:
    """Detect mer-tridentate at CN >= 4 using graph + geometry.

    A planar-tridentate ligand has three donors that are all reachable
    from each other through bonds in the same ligand fragment (without
    crossing the metal). Combined with the geometric mer-cosine pattern
    (i-M-j ~180°, k ~90° to both), this is the universal signature of
    pincer / terpy / mer-N3 chelates.

    Returns ``"mer_tridentate"`` if the pattern is found, else ``None``.
    """
    try:
        nbrs = _mol_neighbors(mol)
    except Exception:
        return None
    # BFS each donor through the ligand (blocking the metal). Group donors
    # by their reachable atom-set hash to find donors that share a chelate.
    blocked = {int(metal_idx)}
    reach_by_donor: Dict[int, frozenset] = {}
    for d in donor_idxs:
        try:
            reach = frozenset(_bfs_fragment(nbrs, int(d), blocked))
        except Exception:
            continue
        reach_by_donor[int(d)] = reach
    # Find triples that are all in the same connected (ligand) fragment.
    units: Dict[int, np.ndarray] = {}
    for d in donor_idxs:
        v = coords[int(d)] - coords[int(metal_idx)]
        n = float(np.linalg.norm(v))
        if n < 1e-9:
            continue
        units[int(d)] = v / n
    donor_list = [d for d in donor_idxs if d in units and d in reach_by_donor]
    nd = len(donor_list)
    for i in range(nd):
        di = donor_list[i]
        ri = reach_by_donor[di]
        for j in range(i + 1, nd):
            dj = donor_list[j]
            # Same fragment iff dj is reachable from di through ligand bonds.
            if dj not in ri:
                continue
            cij = float(np.dot(units[di], units[dj]))
            if cij > -0.90:        # need near-trans pair
                continue
            for k in range(nd):
                if k in (i, j):
                    continue
                dk = donor_list[k]
                # k must also belong to the same chelate fragment.
                if dk not in ri:
                    continue
                cik = float(np.dot(units[di], units[dk]))
                cjk = float(np.dot(units[dj], units[dk]))
                if abs(cik) < 0.30 and abs(cjk) < 0.30:
                    return "mer_tridentate"
    return None


def _smart_phase_a5_should_skip(label: Optional[str]) -> bool:
    """Whether a distortion label means Phase A.5 must be SKIPPED entirely.

    Skip when the chemically-correct geometry is non-symmetric (JT,
    cyclometal 5-ring) — the standard polyhedron would actively damage it.
    For ``"metal_oxo"`` we also skip because the short multiple bond changes
    the M-D radius and we do not yet have an asymmetric-distance projection.
    ``"mer_tridentate"`` is preserved by skipping too (we'd need an
    isomer-aware target which is not yet generated automatically).

    Returns False for None (no distortion → run standard projection) and
    True for any non-None label.
    """
    return label is not None


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
    smart_mode: bool = False,
    report: Optional[Dict] = None,
) -> Tuple[np.ndarray, int]:
    """For each metal, pick the matching ideal polyhedron, then translate each
    donor's BFS fragment toward its ideal position by ``step_size`` of the
    delta.  Rolls each fragment-move back if it breaks topology.

    Parameters
    ----------
    smart_mode : bool
        When True, each metal is first classified via
        :func:`_classify_distortion_mandatory`; if a distortion-mandatory
        pattern is detected (Jahn-Teller, M=O, cyclometal 5-ring,
        mer-tridentate), the standard symmetric projection is skipped for
        that metal so the chemically-correct distortion is preserved. The
        skip decision is recorded in ``report["smart_a5_skipped_metals"]``
        as a list of ``(metal_idx, label)`` tuples.
    report : Optional[Dict]
        Mutable report dictionary; populated with per-metal skip metadata
        when ``smart_mode`` is True. Passing None disables reporting.

    Returns
    -------
    (new_coords, n_projections_committed)
    """
    if not _SCIPY_OK:
        return coords, 0
    out = coords.copy()
    n_done = 0
    blocked = set(metals)
    skipped: List[Tuple[int, str]] = []
    for m in metals:
        donor_idxs = _donor_indices(mol, m)
        if len(donor_idxs) < 2:
            continue
        if smart_mode:
            distortion_label = _classify_distortion_mandatory(
                mol, out, m, donor_idxs
            )
            if _smart_phase_a5_should_skip(distortion_label):
                skipped.append((int(m), str(distortion_label)))
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
    if smart_mode and report is not None and skipped:
        report["smart_a5_skipped_metals"] = skipped
    return out, n_done


# ---------------------------------------------------------------------------
# Phase B — Iterative Gauss-Seidel sweeps
# ---------------------------------------------------------------------------

def _stage_bonds(coords: np.ndarray, mol,
                 md_pairs: Sequence[Tuple[int, int, float]],
                 metals: Sequence[int], nb_pairs: Sequence[Tuple[int, int]],
                 step_size: float, bond_tol: float,
                 h_nbrs: Optional[List[List[int]]] = None,
                 ) -> Tuple[np.ndarray, bool, int]:
    """Mass-weighted bond-distance projection over all M-D pairs (and any
    further covalent bonds that are far from ideal).  Returns (coords, moved,
    n_repairs).

    If ``h_nbrs`` is provided, bonded H atoms of M and D are dragged with the
    same delta as their heavy parent (rigid-H tracking, Patch P-H-TRACK).
    """
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
        delta_m = unit * shift * a_j
        delta_d = -unit * shift * a_i
        new[m] = out[m] + delta_m
        new[d] = out[d] + delta_d
        if h_nbrs is not None:
            for h_idx in h_nbrs[m]:
                new[h_idx] = out[h_idx] + delta_m
            for h_idx in h_nbrs[d]:
                new[h_idx] = out[h_idx] + delta_d
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
    fighting the angle-corrector module.

    Patch P-RING-AWARE (env-flag ``DELFIN_B5_STAGE2_RING_AWARE``, default 0):
    when set to ``1``, and ``x`` lies on an SSSR ring, the rotated subtree is
    expanded to cover **every atom of every ring containing ``x``** plus the
    bonded H children of those ring atoms (see :func:`_ring_aware_subtree`).
    The terminal-ish size cap is bumped from 6 to 12 to accommodate a typical
    aromatic ring + ortho-H bundle (6 ring C + 5 H = 11).  When the env-flag
    is unset / 0, the function is bit-exact to its pre-patch behaviour.

    Rationale: the legacy size-6 BFS subtree rotates only a fraction of a
    small ring, leaving the remainder stationary; the partial rotation
    stretches ring-internal C-C bonds and C-H bonds on the unrotated ring
    atoms, producing ORPHAN_H (nearest heavy > 1.40 Å) cascades in
    hapto-CN7 η-arene / η⁵-Cp complexes (per Agent 8 forensik, +249
    ORPHAN_H over 28 files in rigid-H pool).  Ring-aware expansion keeps
    the ring rigid under the rotation.
    """
    out = coords.copy()
    moved = False
    repairs = 0
    # Read the env-flag once per call; failure → default OFF.
    try:
        _ring_env_on = int(os.environ.get(
            "DELFIN_B5_STAGE2_RING_AWARE", "0")) != 0
    except Exception:
        _ring_env_on = False
    # Size cap: classic 6 for the BFS-only path, 12 for the ring-aware
    # expansion (covers benzene + ortho-H = 11).
    _size_cap = 12 if _ring_env_on else 6
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
            # Ring-aware subtree expansion only when env-flag is set AND
            # x is a ring atom; otherwise classic BFS preserves bit-exact
            # pre-patch behaviour.
            x_in_ring = False
            if _ring_env_on:
                try:
                    x_in_ring = bool(mol.GetAtomWithIdx(x).IsInRing())
                except Exception:
                    x_in_ring = False
            if _ring_env_on and x_in_ring:
                subtree = _ring_aware_subtree(mol, nbrs, x, blocked)
            else:
                subtree = _bfs_fragment(nbrs, x, blocked)
            if len(subtree) > _size_cap:
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
                   step_size: float, clash_factor: float,
                   h_nbrs: Optional[List[List[int]]] = None,
                   ) -> Tuple[np.ndarray, bool, int]:
    """Resolve non-bonded clashes using a 3-strategy ladder.

    Strategy 1: rotate one of the offending atoms' fragment around its M-D
    axis (good for inter-ligand clash between two donor-side fragments).
    Strategy 2: mass-weighted pair-push along the pair direction.
    Strategy 3: small perpendicular fragment translation.

    If ``h_nbrs`` is provided, Strategy 2 drags bonded H atoms of the
    pushed-apart pair (i, j) with the same delta as their heavy parent
    (rigid-H tracking, Patch P-H-TRACK). Strategies 1 & 3 already operate
    on BFS fragments which include H descendants and need no change.
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
        delta_i = -unit * push * a_j
        delta_j = unit * push * a_i
        new[i] = out[i] + delta_i
        new[j] = out[j] + delta_j
        if h_nbrs is not None:
            for h_idx in h_nbrs[i]:
                new[h_idx] = out[h_idx] + delta_i
            for h_idx in h_nbrs[j]:
                new[h_idx] = out[h_idx] + delta_j
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
    rigid_h: bool = False,
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
    rigid_h : bool
        If True, bonded H atoms are dragged with their heavy parent in
        single-atom moves (Stage 1 bond shift on M/D, Stage 3 Strategy 2
        clash push on i/j). Preserves C-H / N-H / O-H bond lengths that the
        heavy-atom-only topology gate does not check. Default False for
        bit-exact pre-patch behaviour.

        Patch P-BH-NOGO (2026-05-13): if the molecule contains a bridging
        μ²-B-H-M hydride (RDKit graph: H bonded to both a B atom and a
        metal in ``_METAL_SET``), rigid_h is forced to False for that frame
        regardless of the caller's request. Dragging a bridging H with its
        boron parent during a clash push breaks the B-H-M bridge geometry
        (CANFAY Zr-CN11 borohydride: 4/9 B-H bonds catastrophically broken,
        1.16 → 0.95 Å). The override is universal (graph-based, no SMILES /
        refcode patterns) and only triggers when boron + metal coexist in
        the same molecule.

    Notes
    -----
    Phase A.5 polyhedron projection (``enable_symmetry``) is governed by an
    additional global env-flag ``DELFIN_B5_PHASE_A5_ENABLE`` (default 0).
    When the env-flag is unset/0 the parameter is ignored and Phase A.5 is
    forced OFF, regardless of caller, protecting downstream code from future
    regressions if the default is flipped. Set ``DELFIN_B5_PHASE_A5_ENABLE=1``
    to honour the ``enable_symmetry`` argument. Forensik 2026-05-13 showed
    Phase A.5 over-idealises Cu(II) Jahn-Teller distortions, M=O multiple
    bonds and cyclometallation 5-rings (~5 % of σ-class structures).

    A second env-flag ``DELFIN_B5_PHASE_A5_SMART`` (default 0) enables the
    distortion-aware "smart" mode added 2026-05-13. When set to 1, Phase A.5
    is unconditionally enabled but each metal is first classified via
    :func:`_classify_distortion_mandatory`. Metals carrying a chemically-
    mandatory distortion (Jahn-Teller Cu(II)/Mn(III)/Cr(II)/Ag(II)/Ni(III);
    metal-oxo M=O / M=NR multiple bond; cyclometal 5-ring closure;
    meridional planar tridentate) are SKIPPED so their geometry is
    preserved. Standard metals continue to receive standard polyhedron
    projection. The ``report`` dict gains ``smart_a5_skipped_metals`` —
    a list of ``(metal_idx, label)`` tuples — when any metal is skipped.
    Setting ``DELFIN_B5_PHASE_A5_SMART=0`` (default) preserves the legacy
    P-A5-ENV behaviour bit-exact.

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

        # ---- Patch P-BH-NOGO: bridging μ²-B-H-M hydride precheck -----------
        # If the molecule has a B-H-M bridge, rigid-H dragging detaches the
        # bridging hydrogen and breaks the bond geometry (see CANFAY Zr-CN11
        # borohydride: 4/9 B-H bonds collapsed to ~0.95 Å). Force rigid_h
        # OFF for the entire post-optimize call when such a bridge exists.
        rigid_h_effective: bool = bool(rigid_h)
        if rigid_h_effective and _has_bridging_bh_hydride(mol):
            rigid_h_effective = False
            report["bh_bridge_nogo"] = True
        h_nbrs: Optional[List[List[int]]] = (
            _compute_h_neighbors(mol) if rigid_h_effective else None
        )

        # ---- Patch P-A5-ENV: Phase A.5 polyhedron-projection env-gate ------
        # Phase A.5 over-idealises Cu(II) Jahn-Teller, M=O double bonds and
        # cyclometallation 5-rings (forensik 2026-05-13). Default-OFF env
        # override prevents accidental re-enable from upstream defaults.
        # Set DELFIN_B5_PHASE_A5_ENABLE=1 to honour the caller's parameter.
        try:
            _a5_env_raw = os.environ.get("DELFIN_B5_PHASE_A5_ENABLE", "0")
            _a5_env_on = int(_a5_env_raw) != 0
        except Exception:
            _a5_env_on = False
        enable_symmetry_effective: bool = (
            bool(enable_symmetry) if _a5_env_on else False
        )

        # ---- Patch P-A5-SMART: Jahn-Teller / multiple-bond / cyclometal-aware
        # smart-projection mode. When DELFIN_B5_PHASE_A5_SMART=1, Phase A.5
        # is unconditionally enabled (regardless of P-A5-ENV) AND each metal
        # is classified by :func:`_classify_distortion_mandatory` before
        # projection. Metals carrying a chemically-mandatory distortion are
        # SKIPPED so their geometry is preserved (Cu(II) 4+2 elongation,
        # M=O short bond, cyclometalation 5-ring, mer-tridentate plane).
        # Standard metals continue to receive standard polyhedron projection.
        #
        # When the env-flag is unset / 0, behaviour is identical to the
        # legacy P-A5-ENV gate (Phase A.5 disabled unless P-A5-ENV=1).
        try:
            _a5_smart_raw = os.environ.get("DELFIN_B5_PHASE_A5_SMART", "0")
            _a5_smart_on = int(_a5_smart_raw) != 0
        except Exception:
            _a5_smart_on = False
        if _a5_smart_on:
            # Smart mode overrides the legacy enable gate — Phase A.5 is
            # always considered for the chemistry-class allowlist below,
            # but the smart projector itself skips distortion-mandatory
            # metals so the override is safe.
            enable_symmetry_effective = True

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
        # Gated by ``enable_symmetry_effective`` (the env-aware override).
        # ``_a5_smart_on`` activates Jahn-Teller / M=O / cyclometal /
        # mer-tridentate aware skipping inside the projector.
        if enable_symmetry_effective and class_label not in (
                "hapto", "multi_hapto", "no_metal"):
            coords, n_sym = _project_donors_to_polyhedron(
                coords, mol, nbrs, metals, md_pairs, nb_pairs,
                step_size=step_size, clash_factor=clash_factor,
                smart_mode=_a5_smart_on, report=report)
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
                cur_step, bond_tol, h_nbrs=h_nbrs)
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
                cur_step, clash_factor, h_nbrs=h_nbrs)
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
