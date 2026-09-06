"""Baustein 6 — Variational L-BFGS-B Post-Refiner with 4-Tier Symmetry Awareness.

Complements Baustein 5 (PBD, ``delfin.manta._post_optimizer``).  Where B5 handles
catastrophic moves + hard topology repair via constraint projection, B6
performs smooth simultaneous balance of all forces (bonds, angles, clashes,
hard topology, plus four tiers of symmetry pressure) via gradient-based
minimization of an 8-term energy functional ``U_total``.

Pipeline
--------
1. Pre-compute symmetry info (Tier A / B / C / D — one-shot, expensive).
2. Build the ``U_total`` + ``grad_U_total`` closures via ``delfin.manta._energy_terms``.
3. Run ``scipy.optimize.minimize(method="L-BFGS-B", jac=True)``.
4. Validate result against the Baustein 5 topology hard-gate; on failure
   the input XYZ is returned unchanged.
5. Return refined XYZ string + a structured report dict.

Doctrine
--------
- Universal.  No SMILES/refcode/element-list shortcuts.  Symmetry detection
  uses Morgan ranks + SMARTS archetypes + connectivity automorphisms.
- Topology preservation is a hard gate (M-D ∈ [0.85, 1.10] × ideal, no new
  spurious heavy-heavy bonds inside 0.85 · Σr_cov).
- Always returns a valid XYZ string.  On any failure (missing sister module,
  scipy exception, topology check fail) the original XYZ is returned with a
  populated report dict and ``fallback_used=True``.
- Lazy imports for RDKit / SciPy / sister B6 modules.  Importing this module
  is cheap and never raises even if helpers are missing.

Entry point
-----------
``variational_refine(xyz, mol, class_label=..., **params) -> (xyz, report)``
"""
from __future__ import annotations

import math
import os
import re
import traceback
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Patch P-H-TRACK (B6 port) — environment flag helper.
# ---------------------------------------------------------------------------

def _delfin_env_int(name: str, default: int) -> int:
    """Read ``int`` value from environment variable ``name`` with safe fallback."""
    try:
        return int(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Class-conditional hyperparameters (see BAUSTEIN6_MASTERPLAN section 6).
# ---------------------------------------------------------------------------

_CLASS_HYPER: Dict[str, Dict[str, Any]] = {
    "sigma": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0, "k_torsion": 50.0,
        "k_topology": 10000.0,
        "k_A": 100.0, "k_B": 50.0, "k_C": 80.0, "k_D": 30.0,
        "enable_D": True,
    },
    "hapto": {
        "k_bond": 800.0, "k_angle": 80.0, "k_clash": 400.0, "k_torsion": 30.0,
        "k_topology": 10000.0,
        "k_A": 30.0, "k_B": 30.0, "k_C": 50.0, "k_D": 20.0,
        "enable_D": True,
    },
    "multi_sigma": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0, "k_torsion": 50.0,
        "k_topology": 10000.0,
        "k_A": 100.0, "k_B": 60.0, "k_C": 100.0, "k_D": 40.0,
        "enable_D": True,
    },
    "multi_hapto": {
        "k_bond": 600.0, "k_angle": 50.0, "k_clash": 300.0, "k_torsion": 30.0,
        "k_topology": 10000.0,
        "k_A": 20.0, "k_B": 20.0, "k_C": 30.0, "k_D": 10.0,
        "enable_D": False,
    },
    "no_metal": {
        "k_bond": 1000.0, "k_angle": 100.0, "k_clash": 500.0, "k_torsion": 50.0,
        "k_topology": 0.0,
        "k_A": 0.0, "k_B": 100.0, "k_C": 150.0, "k_D": 50.0,
        "enable_D": True,
    },
}


# ---------------------------------------------------------------------------
# Constants — lazy fallbacks for _METAL_SET / topology helpers.
# ---------------------------------------------------------------------------

def _load_metal_set() -> set:
    """Lazy import of the canonical metal set."""
    try:
        from delfin.smiles_converter import _METAL_SET  # type: ignore
        return set(_METAL_SET)
    except Exception:
        # Minimal d/f-block fallback so module-load never fails.
        return {
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
            "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
            "Er", "Tm", "Yb", "Lu", "Th", "U",
        }


# ---------------------------------------------------------------------------
# XYZ I/O — DELFIN-style (preserves header lines).
# ---------------------------------------------------------------------------

_XYZ_LINE_RE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
    r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)


def _parse_xyz_to_array(xyz: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Parse a DELFIN-style XYZ string.

    Returns
    -------
    symbols : list[str]
        Element symbol per atom (length N).
    coords : ndarray, shape (N, 3)
        Cartesian coordinates in Å.
    orig_lines : list[str]
        Raw input lines (header preserved for round-trip rewrite).
    """
    syms: List[str] = []
    pts: List[List[float]] = []
    lines = xyz.splitlines()
    for line in lines:
        m = _XYZ_LINE_RE.match(line)
        if m:
            syms.append(m.group(1))
            pts.append([float(m.group(2)), float(m.group(3)), float(m.group(4))])
    if not pts:
        return syms, np.zeros((0, 3), dtype=float), lines
    return syms, np.asarray(pts, dtype=float), lines


def _array_to_xyz(coords: np.ndarray, mol, orig_lines: Optional[List[str]] = None,
                  symbols: Optional[List[str]] = None) -> str:
    """Render coords back into a DELFIN-style XYZ string.

    If ``orig_lines`` are supplied, header / blank lines are preserved verbatim
    and only the atom rows are rewritten in place.  Otherwise a minimal
    two-line header (atom count + blank comment) is emitted.
    """
    if symbols is None:
        try:
            symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        except Exception:
            symbols = ["X"] * len(coords)

    if orig_lines is not None:
        out: List[str] = []
        atom_i = 0
        for line in orig_lines:
            m = _XYZ_LINE_RE.match(line)
            if m and atom_i < len(symbols):
                x, y, z = coords[atom_i]
                out.append(f"{symbols[atom_i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
                atom_i += 1
            else:
                out.append(line)
        return "\n".join(out) + "\n"

    # Plain emission (used when no original header is available).
    out = [f"{len(coords)}", ""]
    for sym, (x, y, z) in zip(symbols, coords):
        out.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Topology check — reuse Baustein 5 if available, otherwise inline minimal.
# ---------------------------------------------------------------------------

def _topology_violation(coords: np.ndarray, mol, metal_set: set):
    """Quantify the hard-gate violation instead of collapsing it to a bool.

    WHY (2026-07-30).  :func:`_topology_check` is ABSOLUTE: it asks "is this frame inside
    the band", never "is this frame better than the one I was given".  Measured on the
    champion path, ``topo_ok_input=False`` on every frame -- the gate was unpassable
    BEFORE the optimiser ran, so a minimisation that took U from 609439 to 646 in 138
    iterations was discarded as "topology not preserved".  An absolute gate over a
    population that mostly starts out of band is structurally dead code.

    Returns ``(n_md_out, worst_md_excess, n_collapse, worst_collapse_deficit)``:
      * ``n_md_out``     — M-D pairs outside the gate's window
      * ``worst_md_excess`` — largest RELATIVE excursion past the window edge (0 inside)
      * ``n_collapse``   — non-bonded heavy pairs below the collapse line
      * ``worst_collapse_deficit`` — largest relative shortfall (0 if none)

    A vector, not a scalar, so no axis can be traded against another.  Note the
    reference length cancels between input and output, which makes the RELATIVE
    comparison robust to a wrong ideal M-D -- the very quantity last night's CSD
    measurement showed we do not have right yet.

    Returns ``None`` when the canonical Baustein-5 helpers are unavailable (caller then
    keeps the absolute behaviour rather than inventing a second definition).
    """
    try:
        from delfin.manta._post_optimizer import (  # type: ignore
            _metal_indices, _md_pairs, _non_bonded_heavy_pairs, _cov_radius,
            _METAL_SET, _MD_LO, _MD_HI, _COLLAPSE_FRAC,
        )
    except Exception:
        return None
    try:
        metals = _metal_indices(mol)
        md = _md_pairs(mol, metals)
        nb = _non_bonded_heavy_pairs(mol)
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
    except Exception:
        return None

    n_md_out = 0
    worst_md = 0.0
    for (m, d, d_ideal) in md:
        if not d_ideal:
            continue
        r = float(np.linalg.norm(coords[m] - coords[d])) / float(d_ideal)
        # Same window as the absolute gate -- ONE definition, read from one place.
        if r < _MD_LO:
            n_md_out += 1
            worst_md = max(worst_md, _MD_LO - r)
        elif r > _MD_HI:
            n_md_out += 1
            worst_md = max(worst_md, r - _MD_HI)

    n_col = 0
    worst_col = 0.0
    for (i, j) in nb:
        if syms[i] in _METAL_SET or syms[j] in _METAL_SET:
            continue
        r_sum = _cov_radius(syms[i]) + _cov_radius(syms[j])
        if not r_sum:
            continue
        r = float(np.linalg.norm(coords[i] - coords[j])) / float(r_sum)
        if r < _COLLAPSE_FRAC:
            n_col += 1
            worst_col = max(worst_col, _COLLAPSE_FRAC - r)

    return (int(n_md_out), float(worst_md), int(n_col), float(worst_col))


def _topology_not_worse(before, after) -> bool:
    """True iff ``after`` is no worse than ``before`` on EVERY axis of the violation
    vector.  Componentwise, so a smaller M-D excursion can never buy a new clash."""
    if before is None or after is None:
        return False
    return all(a <= b + 1e-12 for a, b in zip(after, before))


def _coordination_sphere_indices(mol, metal_set: set) -> List[int]:
    """Metal atoms and their bonded donors — the atoms the SEATING places.

    WHY (2026-07-30).  M-D lengths come from the measured CN-resolved table and the
    chelate bite from ``bite_close``; both are SET, not searched.  Once a coordinate is
    set to its measured value there is nothing left to improve, so any displacement the
    optimiser applies to it is damage by definition.  Freezing is therefore strictly
    better than penalising it afterwards: it costs no term, no weight, no barrier.

    It also explains the two anti-correlated terms.  ``U_A`` and ``U_topology`` both score
    0.94 -- they penalise real crystals MORE than our own frames, i.e. they assert our
    build is more realistic than reality.  They exist mainly to stop the optimiser from
    destroying what the seating placed; with the sphere frozen that job is gone.
    """
    out: List[int] = []
    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in metal_set:
                continue
            out.append(atom.GetIdx())
            for nb in atom.GetNeighbors():
                out.append(nb.GetIdx())
    except Exception:
        return []
    return sorted(set(out))


# Normal-distribution conversion: p90 - p10 = 2 * 1.2816 * sigma.  Imported from the
# module that owns the band readers so the gate and U_bond's inverse-variance weight can
# never disagree on what a sigma is; the literal is only a standalone-import fallback.
try:
    from delfin.manta._energy_terms import _P10_P90_TO_SIGMA  # type: ignore
except Exception:
    _P10_P90_TO_SIGMA = 2.5631


def _realism_deviation(coords: np.ndarray, mol):
    """``(R_bond, R_angle)`` — how far the frame lies OUTSIDE its own measured bands.

    Per item: the excursion past the nearer band edge (exactly 0 inside the band),
    expressed in units of that bin's own measured sigma.  Summed per class.

    The sigma normalisation is not cosmetic.  It is what makes Å and degrees
    commensurable, and it IS the weighting: a tightly measured bin (the five-ring
    chelate bite is 3.2 deg wide) scores a given deviation harder than a loose one.
    Measured, not guessed -- and crucially it couples to the WIDTH OF THE BAND, a
    property of chemistry, never to the SIZE OF THE ERROR.  Coupling weight to how
    broken a site is would hand a torn ligand enormous forces and blow it apart.

    M-D and D-M-D are deliberately absent: those are handled by freezing the
    coordination sphere, which is stronger than gating it after the fact.

    Returns ``None`` per component when no table is loaded for it, so an unset table
    makes the gate inert on that axis rather than inventing a reference.
    """
    try:
        from delfin.manta._energy_terms import (  # type: ignore
            _measured_bond_band, _measured_angle_band, _enumerate_bonds,
            _enumerate_angles, _bond_band_table, _angle_band_table,
        )
        from delfin.manta._energy_terms import _smiles_converter_metals  # type: ignore
    except Exception:
        return (None, None)

    coords = np.asarray(coords, dtype=np.float64)
    try:
        metals = _smiles_converter_metals()
    except Exception:
        metals = set()

    def _excess(x, band):
        p10, _p50, p90 = band
        sigma = (p90 - p10) / _P10_P90_TO_SIGMA
        if not (sigma > 1.0e-9):
            return None            # degenerate bin — no scale, no statement
        if x < p10:
            return (p10 - x) / sigma
        if x > p90:
            return (x - p90) / sigma
        return 0.0

    r_bond = None
    try:
        if _bond_band_table():
            acc = 0.0
            for (i, j, _order) in _enumerate_bonds(mol):
                if (mol.GetAtomWithIdx(i).GetSymbol() in metals
                        or mol.GetAtomWithIdx(j).GetSymbol() in metals):
                    continue        # M-D: frozen, not gated
                band = _measured_bond_band(mol, i, j)
                if band is None:
                    continue
                e = _excess(float(np.linalg.norm(coords[i] - coords[j])), band)
                if e is not None:
                    acc += e
            r_bond = float(acc)
    except Exception:
        r_bond = None

    r_angle = None
    try:
        if _angle_band_table():
            acc = 0.0
            for (i, j, k) in _enumerate_angles(mol):
                if mol.GetAtomWithIdx(j).GetSymbol() in metals:
                    continue        # D-M-D: frozen, not gated
                band = _measured_angle_band(mol, i, j, k)
                if band is None:
                    continue
                v1 = coords[i] - coords[j]
                v2 = coords[k] - coords[j]
                n1 = float(np.linalg.norm(v1))
                n2 = float(np.linalg.norm(v2))
                if n1 < 1.0e-8 or n2 < 1.0e-8:
                    continue
                cs = float(np.dot(v1, v2) / (n1 * n2))
                cs = max(-1.0 + 1.0e-12, min(1.0 - 1.0e-12, cs))
                e = _excess(math.degrees(math.acos(cs)), band)
                if e is not None:
                    acc += e
            r_angle = float(acc)
    except Exception:
        r_angle = None

    return (r_bond, r_angle)


def _realism_not_worse(before, after) -> bool:
    """True iff ``after`` is no worse than ``before`` on EVERY class.

    Componentwise, never a sum.  ``mean_delta = 0.585`` is exactly what a mean that
    trades good against bad looks like: a summed criterion would let twenty bond lengths
    improve by a permille and pay for it with one destroyed angle.  A component with no
    table on either side is skipped rather than counted as a pass on evidence we do not
    have.
    """
    for b, a in zip(before, after):
        if b is None or a is None:
            continue
        if a > b + 1.0e-9 + 1.0e-6 * abs(b):
            return False
    return True


def _topology_check(coords: np.ndarray, mol, metal_set: set) -> bool:
    """Conservative topology hard-gate.

    Tries the canonical ``delfin.manta._post_optimizer._passes_topology`` first; if
    that module is unavailable, falls back to an inline check that mirrors its
    semantics (M-D window 0.85-1.10 × ideal, non-bonded heavy-pair collapse
    below 0.85 · Σr_cov).
    """
    try:
        from delfin.manta._post_optimizer import (  # type: ignore
            _passes_topology, _metal_indices, _md_pairs, _non_bonded_heavy_pairs,
        )
        metals = _metal_indices(mol)
        md = _md_pairs(mol, metals)
        nb = _non_bonded_heavy_pairs(mol)
        return bool(_passes_topology(coords, mol, metals, md, nb))
    except Exception:
        pass

    # Inline fallback — only M-D window, conservative.
    try:
        metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in metal_set]
    except Exception:
        return True
    if not metals:
        return True

    try:
        from delfin.smiles_converter import _get_ml_bond_length  # type: ignore
    except Exception:
        def _get_ml_bond_length(_a: str, _b: str) -> float:  # type: ignore
            return 2.0

    for m in metals:
        m_atom = mol.GetAtomWithIdx(m)
        m_sym = m_atom.GetSymbol()
        for nb_atom in m_atom.GetNeighbors():
            if nb_atom.GetSymbol() in metal_set or nb_atom.GetSymbol() == "H":
                continue
            d_idx = nb_atom.GetIdx()
            d_sym = nb_atom.GetSymbol()
            try:
                d_ideal = float(_get_ml_bond_length(m_sym, d_sym))
            except Exception:
                d_ideal = 2.0
            d_cur = float(np.linalg.norm(coords[m] - coords[d_idx]))
            if d_cur < 0.85 * d_ideal or d_cur > 1.10 * d_ideal:
                return False
    return True


# ---------------------------------------------------------------------------
# Patch P-H-TRACK (B6 port) — rigid-H DoF reduction helpers.
#
# Mirrors :func:`delfin.manta._post_optimizer._compute_h_neighbors` and adds the
# extra plumbing required for L-BFGS-B coordinate substitution.
#
# Approach
# --------
# B5 corrects in discrete per-atom moves where dragging bonded H is trivial
# (apply the same delta to every H child).  B6 minimises a smooth U_total via
# L-BFGS-B with 3N independent DoFs — H atoms with their own DoFs lag behind
# their heavy parent under angle/clash/topology pressure, stretching C-H /
# N-H / O-H bonds beyond the heavy-atom-only topology gate's notice.
#
# The fix is a coordinate substitution: treat eligible H as rigid offsets
# from their unique heavy parent.  L-BFGS-B sees a reduced coord array
# (heavy atoms + non-rigid H only).  At every evaluation:
#
#   1.  Expand reduced -> full coords by placing each rigid H at
#       ``parent_pos + cached_offset[h]``.
#   2.  Call the full U_total to get (E, grad_full).
#   3.  Fold each rigid-H gradient row into its parent's row
#       (chain rule: H = parent + const => dE/dparent gains dE/dH).
#   4.  Strip the rigid-H rows and return the reduced gradient.
#
# An H atom is eligible iff it has exactly one heavy neighbour (the typical
# terminal C-H / N-H / O-H pattern).  Bridging H, lone H, and any H whose
# only neighbour is itself H are kept as free DoFs.
# ---------------------------------------------------------------------------


def _compute_h_neighbors(mol) -> List[List[int]]:
    """Per-atom list of bonded H atom indices.

    Mirrors :func:`delfin.manta._post_optimizer._compute_h_neighbors`.  Returns a
    list of length ``mol.GetNumAtoms()``; entries indexed by an H atom are
    always empty so the helper is safe to use on any atom index.
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


def _compute_rigid_h_map(mol) -> Tuple[List[int], List[int]]:
    """Identify H atoms eligible to be rigidly tied to a unique heavy parent.

    Returns
    -------
    rigid_h_indices : list[int]
        Sorted indices of H atoms tracked rigidly.
    h_parent : list[int]
        Per-atom parent index (length ``mol.GetNumAtoms()``).  Non-rigid H
        entries (and every non-H atom) hold ``-1``.

    An H is eligible iff:
      * Symbol is ``"H"``.
      * Exactly one neighbour and that neighbour is not H.
    """
    n = mol.GetNumAtoms()
    rigid: List[int] = []
    parent: List[int] = [-1] * n
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) != 1:
            continue
        nb = nbrs[0]
        if nb.GetSymbol() == "H":
            continue
        h_idx = int(atom.GetIdx())
        parent[h_idx] = int(nb.GetIdx())
        rigid.append(h_idx)
    rigid.sort()
    return rigid, parent


def _build_reduce_expand(
    coords: np.ndarray,
    rigid_h_indices: List[int],
    h_parent: List[int],
    frozen_indices: Optional[List[int]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pre-compute coordinate-substitution book-keeping for L-BFGS-B.

    Parameters
    ----------
    coords : (N, 3)
        Initial full Cartesian coordinates (Å).
    rigid_h_indices : list[int]
        H atoms tracked rigidly (output of :func:`_compute_rigid_h_map`).
    h_parent : list[int]
        Per-atom parent index (output of :func:`_compute_rigid_h_map`).
    frozen_indices : list[int], optional
        Atoms held at their INPUT position (output of
        :func:`_coordination_sphere_indices`).  A third category next to "free"
        and "rigid H": rigid H is *reconstructed* from a parent, frozen is simply
        never moved.

    Returns
    -------
    free_indices : (M,) int array
        Atom indices remaining as free L-BFGS-B DoFs (heavy + non-rigid H,
        minus anything frozen).
    free_pos_of_atom : (N,) int array
        Reverse map: ``free_pos_of_atom[atom_idx]`` = row in the reduced
        array for that atom, or ``-1`` for rigid H and for frozen atoms.
    h_offsets : (H, 3) float array
        Cached parent->H offset vectors (same order as ``rigid_h_indices``).
    rigid_h_arr : (H,) int array
        ``np.asarray(rigid_h_indices)`` for fast indexing.
    frozen_arr : (F,) int array
        Frozen atom indices.
    frozen_coords : (F, 3) float array
        Their input positions, written back verbatim at every expansion.
    """
    n = coords.shape[0]
    rigid_set = set(rigid_h_indices)
    frozen_set = set(frozen_indices or ())
    # A frozen atom is never also a rigid-H DoF; frozen wins (it is the stronger
    # statement, and the H offset would be a no-op anyway).
    rigid_h_indices = [h for h in rigid_h_indices if h not in frozen_set]
    rigid_set = set(rigid_h_indices)
    free = [i for i in range(n) if i not in rigid_set and i not in frozen_set]
    free_indices = np.asarray(free, dtype=np.int64)
    free_pos = np.full(n, -1, dtype=np.int64)
    for pos, i in enumerate(free):
        free_pos[i] = pos
    rigid_h_arr = np.asarray(rigid_h_indices, dtype=np.int64)
    if rigid_h_arr.size:
        parents = np.asarray([h_parent[h] for h in rigid_h_indices], dtype=np.int64)
        h_offsets = coords[rigid_h_arr] - coords[parents]
    else:
        h_offsets = np.zeros((0, 3), dtype=np.float64)
    frozen_arr = np.asarray(sorted(frozen_set), dtype=np.int64)
    frozen_coords = (coords[frozen_arr].copy() if frozen_arr.size
                     else np.zeros((0, 3), dtype=np.float64))
    return (free_indices, free_pos, h_offsets, rigid_h_arr,
            frozen_arr, frozen_coords)


def _expand_full(
    reduced: np.ndarray,
    free_indices: np.ndarray,
    rigid_h_arr: np.ndarray,
    h_offsets: np.ndarray,
    h_parent: List[int],
    n_atoms: int,
    frozen_arr: Optional[np.ndarray] = None,
    frozen_coords: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Reconstruct full ``(N, 3)`` coords from the reduced DoF array.

    Order matters: free rows, then frozen rows, then rigid H -- a rigid H may hang
    off a FROZEN parent (an N-H on an amine donor), so both parent categories have
    to be in place before the offsets are applied.
    """
    full = np.empty((n_atoms, 3), dtype=np.float64)
    full[free_indices] = reduced
    if frozen_arr is not None and frozen_arr.size:
        full[frozen_arr] = frozen_coords
    if rigid_h_arr.size:
        parents = np.asarray([h_parent[h] for h in rigid_h_arr.tolist()],
                             dtype=np.int64)
        full[rigid_h_arr] = full[parents] + h_offsets
    return full


def _reduce_grad(
    grad_full: np.ndarray,
    free_indices: np.ndarray,
    free_pos: np.ndarray,
    rigid_h_arr: np.ndarray,
    h_parent: List[int],
) -> np.ndarray:
    """Fold rigid-H gradient rows into the parent rows (chain rule), then
    slice down to the reduced ``(M, 3)`` shape used by L-BFGS-B.
    """
    grad_red = grad_full[free_indices].copy()
    if rigid_h_arr.size:
        for h_idx in rigid_h_arr.tolist():
            p = h_parent[h_idx]
            pos = int(free_pos[p])
            if pos < 0:
                # Parent is FROZEN -> this H is frozen with it.  Its gradient row has
                # no free DoF to fold into and must be dropped, not written to row -1
                # (which would silently corrupt the last free atom's gradient).
                continue
            grad_red[pos] += grad_full[h_idx]
    return grad_red


# ---------------------------------------------------------------------------
# Minimal fallback energy term (U_bond + U_clash + harmonic M-D).
# Used only when delfin.manta._energy_terms is not yet available.
# ---------------------------------------------------------------------------

def _fallback_U_total(coords_2d: np.ndarray, mol, sym_info: Dict[str, Any],
                      params: Dict[str, Any]) -> Tuple[float, np.ndarray]:
    """Bare-minimum U_total used when the full _energy_terms module is missing.

    Includes:
      - Harmonic U_bond over all RDKit bonds (ideal = covalent-radius sum).
      - Harmonic donor-target pull for the Tier-A targets supplied in
        ``sym_info["donor_targets"]``.
      - One-sided quadratic clash penalty on non-bonded heavy pairs.

    Returns (U_total, grad shape (N, 3)).
    """
    n = coords_2d.shape[0]
    grad = np.zeros_like(coords_2d)
    U = 0.0

    # Element symbols (cheap, RDKit-backed).
    try:
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
    except Exception:
        syms = ["C"] * n

    # Covalent radii fallback table.
    cov = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
           "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
           "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
           "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Pt": 1.36}

    def _r(sym: str) -> float:
        return cov.get(sym, 0.80 if sym != "H" else 0.31)

    k_b = float(params.get("k_bond", 1000.0))
    try:
        bonds = list(mol.GetBonds())
    except Exception:
        bonds = []
    for b in bonds:
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        d_ideal = _r(syms[i]) + _r(syms[j])
        diff = coords_2d[i] - coords_2d[j]
        d = float(np.linalg.norm(diff))
        if d < 1e-9:
            continue
        delta = d - d_ideal
        U += k_b * delta * delta
        g = (2.0 * k_b * delta) * (diff / d)
        grad[i] += g
        grad[j] -= g

    # Tier-A donor-target harmonic.
    k_A = float(params.get("k_A", 0.0))
    targets = sym_info.get("donor_targets", {})
    if k_A > 0 and targets:
        for d_idx, tgt in targets.items():
            if d_idx < 0 or d_idx >= n:
                continue
            diff = coords_2d[d_idx] - np.asarray(tgt, dtype=float)
            U += k_A * float(np.dot(diff, diff))
            grad[d_idx] += 2.0 * k_A * diff

    # Light clash term — only sample a bounded subset of non-bonded pairs.
    k_c = float(params.get("k_clash", 0.0))
    if k_c > 0 and n <= 400:
        bonded = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in bonds}
        for i in range(n):
            if syms[i] == "H":
                continue
            for j in range(i + 1, n):
                if syms[j] == "H":
                    continue
                if (i, j) in bonded:
                    continue
                thr = (_r(syms[i]) + _r(syms[j])) * 1.5  # vdW-ish proxy
                diff = coords_2d[i] - coords_2d[j]
                d = float(np.linalg.norm(diff))
                if d < 1e-9 or d >= thr:
                    continue
                pen = thr - d
                U += k_c * pen * pen
                g = (-2.0 * k_c * pen) * (diff / d)
                grad[i] += g
                grad[j] -= g

    return float(U), grad


# ---------------------------------------------------------------------------
# Symmetry pre-compute — calls sister modules when available, returns minimal
# sym_info struct when they are not.
# ---------------------------------------------------------------------------

def _precompute_symmetry(mol, coords: np.ndarray, class_label: str,
                         params: Dict[str, Any], metal_set: set,
                         enable_global_pg: bool) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Build the symmetry-info struct that ``U_total`` consumes.

    Returns ``(sym_info, meta)`` where ``meta`` carries diagnostic counts for
    the final report.  Every key in ``sym_info`` falls back to a safe empty
    container if the corresponding helper is missing.
    """
    sym_info: Dict[str, Any] = {
        "donor_targets": {},
        "equiv_atoms": [],
        "equiv_bond_pairs": [],
        "equiv_angle_triples": [],
        "fragments": [],
        "global_pg": "C1",
        "global_ops": [],
        "atom_perms": {},
    }
    meta: Dict[str, Any] = {
        "global_pg": "C1",
        "fragments_detected": 0,
        "equiv_classes": 0,
    }

    # ----- Tier A: per-metal donor targets via Hungarian assignment -----
    try:
        from delfin.manta._symmetry_detection import (  # type: ignore
            hungarian_assign_donors_to_slots,
        )
        metals = [a.GetIdx() for a in mol.GetAtoms()
                  if a.GetSymbol() in metal_set]
        targets: Dict[int, np.ndarray] = {}
        for m_idx in metals:
            try:
                per_metal = hungarian_assign_donors_to_slots(coords, mol, m_idx)
                if per_metal:
                    targets.update(per_metal)
            except Exception:
                continue
        sym_info["donor_targets"] = targets
    except Exception:
        pass

    # ----- Tier B: equivalent atoms / bond-pairs / angle-triples -----
    try:
        from delfin.manta._symmetry_detection import (  # type: ignore
            find_equivalent_atoms,
            find_equivalent_bond_pairs,
            find_equivalent_angle_triples,
        )
        equiv_atoms = find_equivalent_atoms(mol) or []
        sym_info["equiv_atoms"] = equiv_atoms
        sym_info["equiv_bond_pairs"] = find_equivalent_bond_pairs(mol) or []
        sym_info["equiv_angle_triples"] = find_equivalent_angle_triples(mol) or []
        meta["equiv_classes"] = len(equiv_atoms)
    except Exception:
        pass

    # ----- Tier C and Tier D: REMOVED 2026-07-30 -----
    # Both were measured against 107 real crystals and both report EXACTLY ZERO force, on
    # our own frames as well as on the crystals.  A term whose gradient vanishes everywhere
    # cannot change any outcome, so computing it is pure cost -- and this precompute is the
    # expensive half of the refiner by its own docstring ("one-shot, expensive").
    #
    # The reason is chemistry rather than a fixable bug.  Every complex inspected came out
    # pg=C1: real coordination complexes have no global point group, so Tier D has no
    # operation to apply.  And a ligand inside a complex is distorted enough that the
    # fragment-orbit guard (residual <= 0.35 A) rejects its automorphisms, so Tier C has no
    # fragment to enforce.  That verdict includes the automorphism rewrite made earlier the
    # same day, which replaced a hand-written SMARTS table of 14 named groups: the table was
    # the wrong way to find fragments, and the right way finds none that matter here.
    #
    # Tier A (coordination sphere) and Tier B (Morgan equivalence, ratio 1.83) stay: both
    # measure non-zero and B discriminates.
    #
    # detect_fragment_orbits() and U_C_fragment/U_D_global are left in place, exercised by
    # the module self-tests, so an ablation A/B can be read against the old logs and the
    # idea can be revived if a future class of system actually carries the symmetry.
    #
    # DELFIN_FFREE_TIER_CD=1 restores the precompute, so the ABLATION ARM can run the old
    # behaviour end to end.  Without it the sym_info keys stay empty and U_total skips both
    # terms anyway -- but a removal that cannot be re-run cannot be confirmed, and this
    # project confirms removals rather than assuming them.
    if os.environ.get("DELFIN_FFREE_TIER_CD", "0") == "1":
        try:
            from delfin.manta._fragment_archetypes import (  # type: ignore
                detect_fragment_orbits,
            )
            _tol = float(os.environ.get("DELFIN_FFREE_ORBIT_RMS_TOL", "0.35"))
            frags = detect_fragment_orbits(mol, coords, rms_tol=_tol) or []
            sym_info["fragments"] = frags
            meta["fragments_detected"] = len(frags)
        except Exception:
            pass
        if enable_global_pg and bool(params.get("enable_D", False)):
            try:
                from delfin.manta._symmetry_detection import (  # type: ignore
                    detect_global_point_group,
                )
                pg, ops, perms = detect_global_point_group(mol, coords)
                sym_info["global_pg"] = pg or "C1"
                sym_info["global_ops"] = list(ops) if ops is not None else []
                sym_info["atom_perms"] = perms or {}
                meta["global_pg"] = sym_info["global_pg"]
            except Exception:
                pass

    return sym_info, meta


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

def variational_refine(
    xyz: str,
    mol,
    class_label: str = "sigma",
    max_iter: int = 200,
    ftol: float = 1e-6,
    enable_global_pg: bool = True,
    rigid_h: Optional[bool] = None,
) -> Tuple[str, Dict[str, Any]]:
    """L-BFGS-B variational post-refinement with 4-tier symmetry awareness.

    Parameters
    ----------
    xyz : str
        DELFIN-style XYZ string (post Baustein 5).
    mol : rdkit.Chem.Mol
        Topology source.  Must carry a 3D conformer aligned with ``xyz``.
    class_label : {"sigma", "hapto", "multi_sigma", "multi_hapto", "no_metal"}
        Selects the hyperparameter preset.  Unknown labels fall back to
        ``"sigma"``.
    max_iter : int, default 200
        L-BFGS-B iteration cap.
    ftol : float, default 1e-6
        Objective relative tolerance passed to L-BFGS-B.
    enable_global_pg : bool, default True
        If ``False``, Tier D (global molecular point group) is skipped
        regardless of the class default.
    rigid_h : bool or None, default None
        Patch P-H-TRACK (B6 port).  When ``True``, terminal H atoms are
        tied rigidly to their unique heavy parent via coordinate
        substitution so X-H bond lengths cannot drift under angle /
        clash / topology pressure.  When ``False``, every atom keeps its
        own 3 DoFs (bit-exact pre-patch behaviour).  When ``None``
        (default), the value is read from the environment variable
        ``DELFIN_B6_RIGID_H`` (``0`` ⇒ False, anything else ⇒ True),
        falling back to ``False`` when unset.  Default OFF mirrors the
        B5 :func:`delfin.manta._post_optimizer.post_optimize_geometry` patch.

    Returns
    -------
    refined_xyz : str
        Optimized geometry (or the original ``xyz`` if anything failed /
        topology was broken).
    report : dict
        Diagnostic record.  Always contains the keys ``iterations``,
        ``converged``, ``energy_initial``, ``energy_final``,
        ``topology_preserved``, ``fallback_used``, ``global_pg``,
        ``fragments_detected``, ``equiv_classes``, ``rigid_h``,
        ``rigid_h_count``, and (on failure) ``error``.
    """
    # Resolve the rigid_h flag: explicit arg wins, else env, else False.
    if rigid_h is None:
        rigid_h_flag = bool(_delfin_env_int("DELFIN_B6_RIGID_H", 0))
    else:
        rigid_h_flag = bool(rigid_h)
    # Skeleton report (mandatory keys present even on early return).
    report: Dict[str, Any] = {
        "iterations": 0,
        "converged": False,
        "energy_initial": float("nan"),
        "energy_final": float("nan"),
        "topology_preserved": True,
        "fallback_used": False,
        "global_pg": "C1",
        "fragments_detected": 0,
        "equiv_classes": 0,
        "rigid_h": rigid_h_flag,
        "rigid_h_count": 0,
    }

    # ----- Step 1: parse XYZ -----
    try:
        symbols, coords, orig_lines = _parse_xyz_to_array(xyz)
    except Exception as exc:
        report.update({"error": f"xyz parse: {exc}", "fallback_used": True})
        return xyz, report

    n_atoms = coords.shape[0]
    if n_atoms < 2:
        report.update({"error": "fewer than 2 atoms parsed",
                       "fallback_used": True})
        return xyz, report

    # Optional cross-check: mol atom count should match.
    try:
        if mol is not None and mol.GetNumAtoms() != n_atoms:
            report.update({"error": "mol/xyz atom-count mismatch",
                           "fallback_used": True})
            return xyz, report
    except Exception:
        pass

    # ----- Step 2: hyperparameters -----
    params = dict(_CLASS_HYPER.get(class_label, _CLASS_HYPER["sigma"]))

    # ----- Step 3: SciPy availability gate -----
    try:
        from scipy.optimize import minimize  # type: ignore
    except Exception as exc:
        report.update({"error": f"scipy unavailable: {exc}",
                       "fallback_used": True})
        return xyz, report

    metal_set = _load_metal_set()

    # ----- Step 3.5: gate the INPUT too (diagnosis, 2026-07-30) -----
    # Step 9 applies an ABSOLUTE gate ("topology preserved").  If the frame handed to us
    # ALREADY violates it, no optimiser can ever pass -- and the report would still read
    # "topology not preserved", blaming the functional for damage it did not do.  34 of 34
    # frames were rejected with exactly that message, so the two cases must be told apart:
    # ``topo_ok_input=False`` means the gate was unpassable from the start (the gate has to
    # become RELATIVE), ``True`` means the minimisation really did walk out of the band (the
    # barrier's [lo_frac, hi_frac] and the gate's [0.93, 1.07] disagree).  Diagnostic only.
    try:
        report["topo_ok_input"] = bool(_topology_check(coords, mol, metal_set))
    except Exception:
        report["topo_ok_input"] = None

    # ----- Step 4: symmetry pre-compute -----
    try:
        sym_info, sym_meta = _precompute_symmetry(
            mol, coords, class_label, params, metal_set, enable_global_pg
        )
        report["global_pg"] = sym_meta["global_pg"]
        report["fragments_detected"] = sym_meta["fragments_detected"]
        report["equiv_classes"] = sym_meta["equiv_classes"]
    except Exception as exc:
        report.update({"error": f"symmetry precompute: {exc}",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 5: pick U_total implementation -----
    try:
        from delfin.manta._energy_terms import U_total as _U_total  # type: ignore
        _have_full_U = True
    except Exception:
        _U_total = None  # type: ignore
        _have_full_U = False

    def _eval(coords_2d: np.ndarray) -> Tuple[float, np.ndarray]:
        if _have_full_U and _U_total is not None:
            try:
                U, g = _U_total(coords_2d, mol, sym_info, params)
                g_arr = np.asarray(g, dtype=float).reshape(coords_2d.shape)
                return float(U), g_arr
            except Exception as _u_exc:
                # SILENT DEMOTION (found 2026-07-30).  This bare except is how the 8-term
                # functional could stop running without anyone noticing: detect_fragments()
                # returns FragmentMatch NAMEDTUPLES while U_C_fragment calls frag.get(...),
                # so ANY molecule matching an archetype -- benzene, pyridine, methyl,
                # carboxylate -- raised AttributeError here and the refiner quietly
                # minimised the reduced fallback instead.  Record the first reason so a
                # demotion is visible in the report rather than inferred from a mood.
                if "u_total_fallback" not in report:
                    report["u_total_fallback"] = f"{type(_u_exc).__name__}: {_u_exc}"
                return _fallback_U_total(coords_2d, mol, sym_info, params)
        return _fallback_U_total(coords_2d, mol, sym_info, params)

    # ----- Step 5.5: build rigid-H DoF reduction map -----
    rigid_h_arr: np.ndarray = np.zeros((0,), dtype=np.int64)
    h_offsets: np.ndarray = np.zeros((0, 3), dtype=np.float64)
    free_indices: np.ndarray = np.arange(n_atoms, dtype=np.int64)
    free_pos: np.ndarray = np.arange(n_atoms, dtype=np.int64)
    h_parent: List[int] = [-1] * n_atoms
    frozen_arr: np.ndarray = np.zeros((0,), dtype=np.int64)
    frozen_coords: np.ndarray = np.zeros((0, 3), dtype=np.float64)

    # FREEZE THE COORDINATION SPHERE (2026-07-30, default OFF).  See
    # _coordination_sphere_indices: what the seating SET must not be re-searched.
    _freeze_req = bool(_delfin_env_int("DELFIN_FFREE_FREEZE_SPHERE", 0))
    _frozen_list: List[int] = []
    if _freeze_req:
        _frozen_list = _coordination_sphere_indices(mol, metal_set)
        # Leaving nothing to optimise would hand L-BFGS-B an empty x0.  A complex whose
        # every atom is metal-or-donor has no ligand interior to relax, so the freeze
        # simply does not apply -- report it instead of degrading to a zero-DoF run.
        if len(_frozen_list) >= n_atoms:
            report["freeze_sphere_skipped"] = "no free atoms left"
            _frozen_list = []
        report["frozen_count"] = len(_frozen_list)

    if rigid_h_flag or _frozen_list:
        try:
            if rigid_h_flag:
                rigid_h_indices, h_parent = _compute_rigid_h_map(mol)
            else:
                rigid_h_indices, h_parent = [], [-1] * n_atoms
            (free_indices, free_pos, h_offsets, rigid_h_arr,
             frozen_arr, frozen_coords) = _build_reduce_expand(
                coords, rigid_h_indices, h_parent, _frozen_list
            )
            report["rigid_h_count"] = int(rigid_h_arr.size)
        except Exception as exc:
            # If anything in the substitution prep fails, fall back to plain
            # full-DoF minimisation (still bit-exact when flag is False at
            # the caller's request).
            report["rigid_h"] = False
            rigid_h_flag = False
            rigid_h_arr = np.zeros((0,), dtype=np.int64)
            h_offsets = np.zeros((0, 3), dtype=np.float64)
            free_indices = np.arange(n_atoms, dtype=np.int64)
            free_pos = np.arange(n_atoms, dtype=np.int64)
            h_parent = [-1] * n_atoms
            frozen_arr = np.zeros((0,), dtype=np.int64)
            frozen_coords = np.zeros((0, 3), dtype=np.float64)
            report.setdefault("rigid_h_error", str(exc))

    n_free = int(free_indices.size)
    # One switch for both reductions: rigid H OR a frozen sphere puts us on the
    # reduced path.  Written once so the two can never disagree.
    _reduced = bool((rigid_h_flag and rigid_h_arr.size > 0) or frozen_arr.size > 0)

    def objective(x_flat: np.ndarray) -> Tuple[float, np.ndarray]:
        if _reduced:
            reduced = x_flat.reshape(n_free, 3)
            full = _expand_full(reduced, free_indices, rigid_h_arr,
                                h_offsets, h_parent, n_atoms,
                                frozen_arr, frozen_coords)
            U, grad_full = _eval(full)
            grad_red = _reduce_grad(grad_full, free_indices, free_pos,
                                    rigid_h_arr, h_parent)
            return float(U), np.asarray(grad_red, dtype=float).reshape(-1)
        coords_2d = x_flat.reshape(n_atoms, 3)
        U, grad = _eval(coords_2d)
        return float(U), np.asarray(grad, dtype=float).reshape(-1)

    # ----- Step 6: initial energy (use FULL coords for the diagnostic) -----
    try:
        U_initial, _ = _eval(coords)
        report["energy_initial"] = float(U_initial)
    except Exception as exc:
        report.update({"error": f"initial energy: {exc}",
                       "fallback_used": True})
        return xyz, report

    if not np.isfinite(U_initial):
        report.update({"error": "non-finite initial energy",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 7: L-BFGS-B run -----
    if _reduced:
        x0 = coords[free_indices].flatten()
    else:
        x0 = coords.flatten()
    try:
        result = minimize(
            fun=objective,
            x0=x0,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": int(max_iter), "ftol": float(ftol),
                     "gtol": 1e-5},
        )
        if _reduced:
            reduced_x = np.asarray(result.x, dtype=float).reshape(n_free, 3)
            new_coords = _expand_full(reduced_x, free_indices, rigid_h_arr,
                                      h_offsets, h_parent, n_atoms,
                                      frozen_arr, frozen_coords)
        else:
            new_coords = np.asarray(result.x, dtype=float).reshape(n_atoms, 3)
        U_final = float(result.fun)
        converged = bool(result.success)
        iterations = int(getattr(result, "nit", 0))
    except Exception as exc:
        report.update({
            "error": f"L-BFGS-B: {exc}",
            "fallback_used": True,
            "traceback": traceback.format_exc(limit=3),
        })
        return xyz, report

    report["iterations"] = iterations
    report["converged"] = converged
    report["energy_final"] = U_final

    # ----- Step 8: numerical sanity (NaN/Inf guard) -----
    if not np.all(np.isfinite(new_coords)):
        report.update({"error": "non-finite coords after L-BFGS-B",
                       "fallback_used": True,
                       "topology_preserved": False})
        return xyz, report

    # Reject if the optimizer somehow blew the energy up.
    if not np.isfinite(U_final) or U_final > U_initial + 1.0 + abs(U_initial):
        report.update({"error": "energy increased / non-finite",
                       "fallback_used": True})
        return xyz, report

    # ----- Step 9: topology gate — absolute first, then NEVER-WORSE -----
    try:
        topo_ok = _topology_check(new_coords, mol, metal_set)
    except Exception:
        topo_ok = True  # B5 unavailable → don't penalize

    # RELATIVE FALLBACK (2026-07-30).  The absolute gate asks "is this frame in band",
    # which our own frames mostly are not: topo_ok_input was False on every measured
    # frame, so a minimisation that took U from 609439 to 646 was thrown away for damage
    # it did not cause.  When the INPUT already violates, the only honest question is
    # whether the output is WORSE -- componentwise over (M-D count, worst M-D excursion,
    # collapse count, worst collapse deficit), so no axis can be traded for another.
    # Accepting only non-worse frames is never-worse BY CONSTRUCTION in the refiner's own
    # terms, and the relative comparison is immune to a wrong ideal M-D length (it cancels).
    # Default ON: B6 as a whole is default-OFF, so the champion is untouched either way,
    # and a default-OFF lever inside a default-OFF module would be unmeasurable dead code.
    if not topo_ok and _delfin_env_int("DELFIN_B6_TOPO_RELATIVE", 1):
        _v_in = _topology_violation(coords, mol, metal_set)
        _v_out = _topology_violation(new_coords, mol, metal_set)
        report["topo_viol_input"] = _v_in
        report["topo_viol_output"] = _v_out
        # Only rescue frames the absolute gate could never have passed.  A frame that
        # STARTED in band and left it is real damage and stays rejected.
        if (_v_in is not None and any(_v_in)
                and _topology_not_worse(_v_in, _v_out)):
            topo_ok = True
            report["topo_relative_pass"] = True

    if not topo_ok:
        report.update({"error": "topology not preserved",
                       "fallback_used": True,
                       "topology_preserved": False})
        return xyz, report

    report["topology_preserved"] = True

    # ----- Step 9b: REALISM gate — a second never-worse, on a second axis -----
    # Topology only asks "is the frame still connected and un-collapsed".  It says nothing
    # about whether the geometry got more or less like a real crystal, and the measurement
    # says that is exactly where the functional is ambivalent: it repairs the coarse cases
    # and degrades the fine ones (mean_delta = 0.585).  A mean cannot express that, so the
    # comparison is componentwise over the measured classes, in units of each bin's own
    # sigma.  Reference-free: the bands ARE the reference, no crystal is consulted.
    if _delfin_env_int("DELFIN_FFREE_REALISM_GATE", 0):
        try:
            _r_in = _realism_deviation(coords, mol)
            _r_out = _realism_deviation(new_coords, mol)
        except Exception as exc:
            # A gate that cannot measure must not reject: falling through keeps the
            # pre-gate behaviour instead of discarding work on a failure of the ruler.
            _r_in = _r_out = (None, None)
            report.setdefault("realism_error", str(exc))
        report["realism_in"] = _r_in
        report["realism_out"] = _r_out
        if not _realism_not_worse(_r_in, _r_out):
            report.update({"error": "realism regressed",
                           "fallback_used": True})
            return xyz, report
        report["realism_gate_pass"] = True

    # ----- Step 10: write output -----
    try:
        refined_xyz = _array_to_xyz(new_coords, mol, orig_lines, symbols)
    except Exception as exc:
        report.update({"error": f"xyz write: {exc}",
                       "fallback_used": True})
        return xyz, report

    return refined_xyz, report


# ---------------------------------------------------------------------------
# Self-test — minimal end-to-end smoke (synthetic, no DELFIN-data dependency).
# ---------------------------------------------------------------------------

def _synthetic_xyz_water() -> Tuple[str, Any]:
    """Mildly distorted water molecule (a triangle with bond lengths perturbed)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    # Perturb each H slightly so U_bond has work to do.
    for i, dx in enumerate([(0.05, 0.0, 0.0), (0.0, -0.08, 0.04)]):
        if i + 1 < mol.GetNumAtoms():
            p = conf.GetAtomPosition(i + 1)
            conf.SetAtomPosition(i + 1, (p.x + dx[0], p.y + dx[1], p.z + dx[2]))
    lines = [str(mol.GetNumAtoms()), "synthetic water"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _synthetic_xyz_methane() -> Tuple[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    conf = mol.GetConformer()
    # Slightly squash one C-H.
    p = conf.GetAtomPosition(1)
    conf.SetAtomPosition(1, (p.x * 0.85, p.y * 0.85, p.z * 0.85))
    lines = [str(mol.GetNumAtoms()), "synthetic methane"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _synthetic_xyz_acetone() -> Tuple[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None

    mol = Chem.MolFromSmiles("CC(C)=O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=11)
    AllChem.UFFOptimizeMolecule(mol, maxIters=50)
    conf = mol.GetConformer()
    # Perturb a hydrogen position.
    p = conf.GetAtomPosition(mol.GetNumAtoms() - 1)
    conf.SetAtomPosition(mol.GetNumAtoms() - 1, (p.x + 0.12, p.y - 0.08, p.z + 0.06))
    lines = [str(mol.GetNumAtoms()), "synthetic acetone"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _synthetic_xyz_metal() -> Tuple[str, Any]:
    """Diphenylzinc, one ring carbon displaced.

    A metal with two donors AND a real ligand interior -- the minimum shape that can
    tell "sphere frozen" apart from "nothing moved at all".
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return "", None
    try:
        mol = Chem.MolFromSmiles("c1ccccc1[Zn]c1ccccc1")
        if mol is None:
            return "", None
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=7) != 0:
            return "", None
    except Exception:
        return "", None
    conf = mol.GetConformer()
    # Displace one carbon well away from where the ring wants it.
    for a in mol.GetAtoms():
        if a.GetSymbol() == "C":
            p = conf.GetAtomPosition(a.GetIdx())
            conf.SetAtomPosition(a.GetIdx(), (p.x + 0.25, p.y - 0.18, p.z + 0.11))
            break
    lines = [str(mol.GetNumAtoms()), "synthetic diphenylzinc"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines) + "\n", mol


def _xyz_coords(xyz: str) -> np.ndarray:
    rows = []
    for ln in xyz.splitlines()[2:]:
        p = ln.split()
        if len(p) >= 4:
            rows.append([float(p[1]), float(p[2]), float(p[3])])
    return np.asarray(rows, dtype=np.float64)


def _self_test_freeze_and_gate() -> None:
    """Prove the two new levers actually DO something when armed.

    A lever that is byte-identical when off and ALSO byte-identical when on is dead code
    that looks alive -- the exact failure mode that hid four separate holes in this
    functional.  Off-identity is checked by diffing this file's output against HEAD;
    these cases check the other half.
    """
    import tempfile

    # --- 1. freeze: the coordination sphere must not move, the ligand must ---
    xyz, mol = _synthetic_xyz_metal()
    if not xyz or mol is None:
        print("[freeze  ] SKIP — RDKit could not build the metal case")
    else:
        c_in = _xyz_coords(xyz)
        sphere = _coordination_sphere_indices(mol, _load_metal_set())
        os.environ["DELFIN_FFREE_FREEZE_SPHERE"] = "1"
        try:
            new_xyz, rep = variational_refine(xyz, mol, class_label="sigma_coord",
                                              max_iter=80, ftol=1e-5,
                                              enable_global_pg=False)
        finally:
            os.environ.pop("DELFIN_FFREE_FREEZE_SPHERE", None)
        c_out = _xyz_coords(new_xyz)
        if c_out.shape != c_in.shape or not sphere:
            print(f"[freeze  ] SKIP — sphere={len(sphere)} shape={c_out.shape}")
        else:
            moved = np.linalg.norm(c_out - c_in, axis=1)
            free = [i for i in range(c_in.shape[0]) if i not in set(sphere)]
            d_sphere = float(moved[sphere].max())
            d_free = float(moved[free].max()) if free else 0.0
            ok = (d_sphere == 0.0) and (d_free > 1.0e-6)
            print(f"[freeze  ] sphere={len(sphere)} frozen_max_move={d_sphere:.3e} "
                  f"free_max_move={d_free:.4f} fallback={rep.get('fallback_used')} "
                  f"OK={ok}")

    # --- 2. _realism_not_worse: the pure decision logic ---
    cases = [
        ((1.0, 2.0), (1.0, 2.0), True,  "identical"),
        ((1.0, 2.0), (0.5, 1.0), True,  "both better"),
        ((1.0, 2.0), (0.5, 9.0), False, "bond better, angle much worse"),
        ((1.0, 2.0), (1.0, 2.1), False, "angle slightly worse"),
        ((None, 2.0), (None, 1.0), True, "missing table skipped"),
        ((None, 2.0), (None, 3.0), False, "missing table does not excuse the other"),
    ]
    bad = 0
    for before, after, want, why in cases:
        got = _realism_not_worse(before, after)
        if got != want:
            bad += 1
            print(f"[gatelogic] FAIL {why}: {before} -> {after} got={got} want={want}")
    print(f"[gatelogic] {len(cases) - bad}/{len(cases)} cases OK")

    # --- 3. _realism_deviation against a synthetic band table ---
    xyz, mol = _synthetic_xyz_acetone()
    if not xyz or mol is None:
        print("[bands   ] SKIP — RDKit not available")
        return
    coords = _xyz_coords(xyz)
    r_none = _realism_deviation(coords, mol)
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        wide = fh.name
        fh.write("L4\tC|C|1.0\t9999\t0.10\t1.50\t9.00\n")
        fh.write("L4\tC|H|1.0\t9999\t0.10\t1.09\t9.00\n")
        fh.write("L4\tC|O|2.0\t9999\t0.10\t1.21\t9.00\n")
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        narrow = fh.name
        fh.write("L4\tC|C|1.0\t9999\t4.90\t5.00\t5.10\n")
        fh.write("L4\tC|H|1.0\t9999\t4.90\t5.00\t5.10\n")
        fh.write("L4\tC|O|2.0\t9999\t4.90\t5.00\t5.10\n")

    os.environ["DELFIN_FFREE_BOND_BANDS"] = wide
    r_wide = _realism_deviation(coords, mol)
    os.environ["DELFIN_FFREE_BOND_BANDS"] = narrow
    r_narrow = _realism_deviation(coords, mol)
    os.environ.pop("DELFIN_FFREE_BOND_BANDS", None)

    ok = (r_none[0] is None                      # no table -> inert, not zero
          and r_wide[0] == 0.0                   # everything inside -> exactly 0
          and r_narrow[0] is not None and r_narrow[0] > 0.0)
    print(f"[bands   ] no_table={r_none[0]} wide={r_wide[0]} "
          f"narrow={r_narrow[0]:.1f} OK={ok}")

    # --- 4. end-to-end: a band the INPUT satisfies and the optimiser must leave ---
    # p10/p90 are pinned to the input's own min/max C-H, so R_in = 0 by construction and
    # any regularisation of the perturbed H pushes R_out above it -> the gate must reject
    # and hand back the INPUT xyz unchanged.
    ch = [float(np.linalg.norm(coords[b.GetBeginAtomIdx()] - coords[b.GetEndAtomIdx()]))
          for b in mol.GetBonds()
          if {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()} == {"C", "H"}]
    if not ch:
        print("[e2e     ] SKIP — no C-H bonds")
        return
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        pinned = fh.name
        fh.write(f"L4\tC|H|1.0\t9999\t{min(ch):.6f}\t"
                 f"{0.5 * (min(ch) + max(ch)):.6f}\t{max(ch):.6f}\n")
    os.environ["DELFIN_FFREE_BOND_BANDS"] = pinned
    os.environ["DELFIN_FFREE_REALISM_GATE"] = "1"
    try:
        gated_xyz, rep_g = variational_refine(xyz, mol, class_label="no_metal",
                                              max_iter=100, ftol=1e-5,
                                              enable_global_pg=False)
    finally:
        os.environ.pop("DELFIN_FFREE_REALISM_GATE", None)
        os.environ.pop("DELFIN_FFREE_BOND_BANDS", None)
    r_in = rep_g.get("realism_in", (None, None))
    r_out = rep_g.get("realism_out", (None, None))
    rejected = rep_g.get("error") == "realism regressed"
    unchanged = (gated_xyz == xyz)
    consistent = (rejected != bool(rep_g.get("realism_gate_pass", False)))
    print(f"[e2e     ] R_in={r_in[0]} R_out={r_out[0]} rejected={rejected} "
          f"input_returned={unchanged} consistent={consistent} "
          f"OK={consistent and (not rejected or unchanged)}")


def _self_test() -> None:
    """Run three synthetic cases and print a one-line summary each."""
    cases = [
        ("water", "no_metal", _synthetic_xyz_water),
        ("methane", "no_metal", _synthetic_xyz_methane),
        ("acetone", "no_metal", _synthetic_xyz_acetone),
    ]
    for name, klass, builder in cases:
        xyz, mol = builder()
        if not xyz or mol is None:
            print(f"[{name}] SKIP — RDKit not available")
            continue
        new_xyz, report = variational_refine(
            xyz, mol, class_label=klass, max_iter=100, ftol=1e-5,
            enable_global_pg=False,
        )
        ok = (report.get("topology_preserved", False)
              and not report.get("fallback_used", True))
        e0 = report.get("energy_initial", float("nan"))
        e1 = report.get("energy_final", float("nan"))
        delta = (e0 - e1) if (np.isfinite(e0) and np.isfinite(e1)) else float("nan")
        print(
            f"[{name:8s} class={klass:9s}] iters={report['iterations']:3d} "
            f"converged={report['converged']!s:5s} "
            f"E0={e0:10.4f} -> E1={e1:10.4f} ΔE={delta:+.4f} "
            f"topo_ok={report['topology_preserved']} "
            f"fallback={report['fallback_used']} OK={ok}"
        )


if __name__ == "__main__":  # pragma: no cover
    _self_test()
    if os.environ.get("DELFIN_B6_SELFTEST_GATE", "0") == "1":
        _self_test_freeze_and_gate()
