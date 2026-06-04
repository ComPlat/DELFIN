"""delfin.fffree.hapto_honest_construction — Build-time hapto-honest
construction corrector.

Construction-side partner to ``agent_workspace/quality_framework/scripts/
hapto_only_cshm.py`` (the metric).  When a built complex contains hapto
units (η²-η⁸ of contiguous donor carbons), this module enforces three
honest geometry constraints that the metric scores:

  1. **M-centroid distance** equal to the empirical CCDC p50 ideal value
     (from :mod:`delfin.fffree.hapto_modes`).  Set per (metal, η).

  2. **Perpendicularity** of M to the ring plane: the M-centroid axis is
     aligned with the ring normal.

  3. **Ring planarity preserved**: the ring atoms (and any attached
     ring-H atoms) are moved as a RIGID BLOCK so internal C-C distances
     and ring-H positions stay exact.  No internal ring deformation.

Rigid-block movement keeps the σ-only sub-polyhedron of the complex
INVARIANT (only the hapto-cluster atoms shift; σ donors are not touched).

The implementation runs after the standard #131 Construction-Fix hooks
and BEFORE the final ``refine.refine()`` clash-relief pass.  Default-OFF:
when none of the env flags below is set, the function returns the input
coordinates byte-identically (``P.copy()`` only — no float arithmetic).

Env flags
---------
* ``DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION=1`` — per-fix activation.
* ``DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1``       — master activation.

Universality contract
---------------------
* Graph + geometry only (no SMILES patterns).
* Deterministic (sorted iteration, no RNG).
* Metal + σ donors are NEVER moved (σ-only metric is byte-invariant).
* M-D distance for σ donors is exactly invariant (they are not in any
  hapto cluster).

Public API
----------
``apply_hapto_honest(symbols, coords, metal_idx, donor_idxs=None)`` ->
``(P_new, n_units_fixed)``.

* ``donor_idxs`` is accepted for parity with sibling corrector modules
  but is NOT used: hapto detection is geometry-only (we look at carbons
  in the M-shell).  When the env flag is unset, returns
  ``(coords.copy(), 0)``.
"""
from __future__ import annotations

import itertools
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np


__all__ = [
    "apply_hapto_honest",
    "detect_hapto_units",
    "honest_active",
]

# ---------------------------------------------------------------------------
# Constants — mirror the metric module so behaviour stays coherent.
# ---------------------------------------------------------------------------
_DEFAULT_MD_CUT: float = 2.80
_HAPTO_CC_MAX:   float = 1.70
_HAPTO_CENTROID_CUT: float = 2.6

#: Tight tolerance for the rigid-block validator (Å).  After moving the
#: hapto cluster + attached H atoms, all *intra-block* pairwise distances
#: must reproduce within this tolerance — guarantees rigidity.
_RIGID_TOL: float = 1e-6

#: Per-(metal, η) ideal M-centroid distances (Å), duplicated from
#: :mod:`delfin.fffree.hapto_modes` so this module never imports it
#: (keeps the env-OFF byte-identical guarantee airtight: no Python-level
#: side effects from a transitive import).
_HAPTO_M_CENTROID_IDEAL: Dict[Tuple[str, int], float] = {
    ("Fe", 5): 1.65, ("Mn", 5): 1.78, ("Co", 5): 1.67,
    ("Ni", 5): 1.78, ("Ru", 5): 1.82, ("Os", 5): 1.85,
    ("Rh", 5): 1.85, ("Ir", 5): 1.88, ("V",  5): 1.90,
    ("Cr", 5): 1.80, ("Mo", 5): 1.96, ("W",  5): 1.96,
    ("Cr", 6): 1.62, ("Fe", 6): 1.55, ("Ru", 6): 1.69,
    ("Os", 6): 1.73, ("Co", 6): 1.61, ("Mo", 6): 1.85,
    ("Fe", 4): 1.70, ("Ru", 4): 1.83, ("Co", 4): 1.72,
    ("Pd", 3): 2.00, ("Pt", 3): 2.05, ("Ni", 3): 1.95,
    ("Pt", 2): 2.10, ("Pd", 2): 2.10, ("Cu", 2): 2.20,
    ("Ag", 2): 2.40, ("Rh", 2): 2.20, ("Ir", 2): 2.25,
    ("Cr", 7): 1.40, ("Mo", 7): 1.60,
    ("U",  8): 1.95, ("Th", 8): 2.00, ("Ce", 8): 1.97,
}

_METAL_SYMBOLS: Set[str] = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Th", "U",
    "Al", "Ga", "In", "Sn", "Pb", "Sb", "Bi",
}

# Covalent radii (Å) — used to attach H atoms to ring carbons (rigid drag).
_COV_R: Dict[str, float] = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "B": 0.84, "Si": 1.11, "Se": 1.20, "Te": 1.38, "As": 1.19,
}
_DEFAULT_COV: float = 0.76


def _ideal_m_centroid_distance(metal: str, eta: int) -> float:
    key = (metal, int(eta))
    if key in _HAPTO_M_CENTROID_IDEAL:
        return _HAPTO_M_CENTROID_IDEAL[key]
    return 1.8 * (5.0 / max(int(eta), 2)) ** 0.5


# ---------------------------------------------------------------------------
# Env flag.
# ---------------------------------------------------------------------------
def honest_active() -> bool:
    """Return ``True`` when the hapto-honest construction is enabled.

    Either the per-fix env var or the global ``CONSTRUCTION_FIX_ALL``
    master switch activates the corrector.  When neither is set the
    corrector is byte-identically off.
    """
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_HAPTO_HONEST_CONSTRUCTION", "0") == "1"


def _post_correction_clash_gate_active() -> bool:
    """Surgical-fix env-flag (2026-06-04, b00f9a0 forensik).

    When set, after each rigid-block hapto correction the global heavy-heavy
    / X-H collapse count is recomputed via ``build_time_clash_gate`` and the
    correction is reverted when the count INCREASES vs the pre-correction
    state.  This catches the failure mode where snapping a ring centroid to
    its ideal M-centroid distance pushes ring atoms into nearby ligand
    atoms.  Default OFF -> byte-identical to b00f9a0.
    """
    return os.environ.get(
        "DELFIN_FFFREE_HAPTO_HONEST_CLASH_GATE", "0"
    ) == "1"


# ---------------------------------------------------------------------------
# Detection — graph + geometry only.
# ---------------------------------------------------------------------------
@dataclass
class HaptoUnit:
    indices: List[int]       # sorted carbon indices comprising the hapto unit
    eta: int                 # = len(indices)
    centroid: np.ndarray     # (3,) geometric centroid
    normal: np.ndarray       # (3,) unit ring normal (oriented towards metal)
    centroid_dist: float     # current |M - centroid|
    attached_h: List[int]    # H atom indices bonded to ring carbons
    metal_symbol: str        # for ideal-distance lookup


def _is_metal(sym: str) -> bool:
    return sym in _METAL_SYMBOLS


def _md_shell_carbons(symbols: Sequence[str], coords: np.ndarray,
                       metal: int, md_cut: float) -> List[int]:
    """Return sorted donor-C indices within ``md_cut`` of the metal."""
    if metal < 0 or metal >= len(symbols):
        return []
    M = coords[metal]
    out: List[int] = []
    for i, s in enumerate(symbols):
        if i == metal or s != "C":
            continue
        if float(np.linalg.norm(coords[i] - M)) <= md_cut:
            out.append(i)
    return sorted(out)


def _attached_hydrogens(symbols: Sequence[str], coords: np.ndarray,
                        ring_indices: Sequence[int]) -> List[int]:
    """Return sorted H atom indices bonded to ring carbons (the rigid-drag
    set).  Bond = within covalent-sum × 1.30 Å, with H ≤ 1.3 Å.
    """
    out: Set[int] = set()
    ring_set = set(int(i) for i in ring_indices)
    for i, s in enumerate(symbols):
        if s != "H" or i in ring_set:
            continue
        h_pos = coords[i]
        for c in ring_indices:
            d = float(np.linalg.norm(coords[c] - h_pos))
            cov = _COV_R.get("H", _DEFAULT_COV) + _COV_R.get("C", _DEFAULT_COV)
            if d <= cov * 1.30:
                out.add(int(i))
                break
    return sorted(out)


def _ring_normal_towards(centred: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Compute the ring normal via SVD and orient it so that
    ``dot(normal, axis) > 0`` (towards the metal)."""
    try:
        _U, _S, Vt = np.linalg.svd(centred, full_matrices=False)
        normal = Vt[-1]
    except np.linalg.LinAlgError:
        normal = np.array([0.0, 0.0, 1.0])
    nrm = float(np.linalg.norm(normal))
    if nrm < 1e-9:
        normal = np.array([0.0, 0.0, 1.0])
    else:
        normal = normal / nrm
    if float(np.dot(normal, axis)) < 0.0:
        normal = -normal
    return normal


def detect_hapto_units(symbols: Sequence[str], coords: np.ndarray,
                       metal_idx: int,
                       md_cut: float = _DEFAULT_MD_CUT) -> List[HaptoUnit]:
    """Identify hapto units (η²-η⁸) around ``metal_idx`` from geometry.

    A unit is a connected component (over C-C ≤ ``HAPTO_CC_MAX``) of
    donor carbons in the M-shell with size ≥ 2 AND mean atom-to-centroid
    distance ≤ ``HAPTO_CENTROID_CUT``.  Returns a deterministically
    ordered list (sorted by min atom index).
    """
    if metal_idx < 0 or metal_idx >= len(symbols):
        return []
    M = coords[metal_idx]
    donor_C = _md_shell_carbons(symbols, coords, metal_idx, md_cut)
    if len(donor_C) < 2:
        return []
    # Build C-C adjacency by geometric proximity.
    adj: Dict[int, Set[int]] = {i: set() for i in donor_C}
    for a, b in itertools.combinations(donor_C, 2):
        if float(np.linalg.norm(coords[a] - coords[b])) <= _HAPTO_CC_MAX:
            adj[a].add(b)
            adj[b].add(a)
    seen: Set[int] = set()
    out: List[HaptoUnit] = []
    metal_sym = str(symbols[metal_idx])
    for start in sorted(donor_C):
        if start in seen:
            continue
        stack = [start]
        comp: List[int] = []
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            comp.append(cur)
            for nb in sorted(adj[cur]):
                if nb not in seen:
                    stack.append(nb)
        comp.sort()
        if len(comp) < 2 or len(comp) > 8:
            continue
        ring_xyz = coords[comp]
        cen = ring_xyz.mean(axis=0)
        mean_r = float(np.mean(np.linalg.norm(ring_xyz - cen, axis=1)))
        if mean_r > _HAPTO_CENTROID_CUT:
            continue
        axis = cen - M
        normal = _ring_normal_towards(ring_xyz - cen, axis)
        attached_h = _attached_hydrogens(symbols, coords, comp)
        out.append(HaptoUnit(
            indices=comp,
            eta=len(comp),
            centroid=cen,
            normal=normal,
            centroid_dist=float(np.linalg.norm(axis)),
            attached_h=attached_h,
            metal_symbol=metal_sym,
        ))
    out.sort(key=lambda u: u.indices[0])
    return out


# ---------------------------------------------------------------------------
# Construction logic — rigid-block transform per hapto unit.
# ---------------------------------------------------------------------------
def _rotation_matrix_from_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return the rotation matrix that aligns unit vector ``a`` onto unit
    vector ``b``.  Uses the Rodrigues formula; falls back to a 180°
    rotation about a perpendicular axis when the vectors are anti-parallel.
    """
    aa = a / max(float(np.linalg.norm(a)), 1e-12)
    bb = b / max(float(np.linalg.norm(b)), 1e-12)
    c = float(np.dot(aa, bb))
    if c > 1.0 - 1e-12:
        return np.eye(3)
    if c < -1.0 + 1e-12:
        # 180° rotation about any axis perpendicular to aa.
        ref = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(aa, ref))) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        axis = np.cross(aa, ref)
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        # Rodrigues for 180°: R = 2 outer(axis, axis) - I.
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    v = np.cross(aa, bb)
    s = float(np.linalg.norm(v))
    K = np.array([
        [0.0, -v[2],  v[1]],
        [v[2],  0.0, -v[0]],
        [-v[1], v[0],  0.0],
    ])
    return np.eye(3) + K + K @ K * ((1.0 - c) / (s * s))


def _apply_unit_correction(coords: np.ndarray,
                           metal_pos: np.ndarray,
                           unit: HaptoUnit) -> np.ndarray:
    """Return new coords where ``unit.indices`` + ``unit.attached_h`` are
    rigid-block transformed so that:

      * the ring centroid sits at distance ``ideal(metal, η)`` from M, AND
      * the M-centroid axis is parallel to the ring normal.

    The metal AND all atoms NOT in the rigid block are byte-invariant.
    """
    ideal_d = _ideal_m_centroid_distance(unit.metal_symbol, unit.eta)
    cen = unit.centroid
    normal = unit.normal
    # Target centroid: along ring-normal from the metal at the ideal
    # distance.  We slide along normal (not along M->centroid) so the
    # FINAL geometry has the ring axis along the M-centroid line.
    # Equivalently we move the rigid block so its centroid lands at
    # M + ideal_d * normal_oriented_towards_centroid.  Then we rotate
    # the block so the ring normal aligns onto (target_centroid - M).
    target_axis = normal  # already oriented towards metal-positive side
    target_centroid = metal_pos + target_axis * ideal_d
    # Translation: rigidly shift the block centroid from `cen` to
    # `target_centroid`.
    block_idx = list(unit.indices) + list(unit.attached_h)
    P_new = coords.copy()
    delta = target_centroid - cen
    for k in block_idx:
        P_new[k] = coords[k] + delta
    # After translation the ring normal hasn't changed direction (rigid).
    # The M-centroid vector now is target_centroid - metal_pos == ideal_d*normal,
    # so by construction the M-centroid axis is ALREADY aligned with the
    # ring normal.  No rotation is needed in the ideal case.
    return P_new


def _validate_rigid_block(original: np.ndarray, updated: np.ndarray,
                          block_idx: Sequence[int],
                          tol: float = _RIGID_TOL) -> bool:
    """Verify intra-block pairwise distances are preserved within ``tol``.
    Returns ``True`` when the block is rigid, ``False`` otherwise.
    """
    if len(block_idx) < 2:
        return True
    for a, b in itertools.combinations(block_idx, 2):
        d_old = float(np.linalg.norm(original[a] - original[b]))
        d_new = float(np.linalg.norm(updated[a] - updated[b]))
        if abs(d_old - d_new) > tol:
            return False
    return True


def _validate_sigma_invariance(original: np.ndarray, updated: np.ndarray,
                                metal: int, block_idx: Sequence[int],
                                tol: float = 1e-9) -> bool:
    """Verify that no atom OUTSIDE the rigid block (and not the metal) is
    moved.  Defence in depth — we promise σ donors are untouched."""
    block = set(int(i) for i in block_idx)
    for i in range(len(original)):
        if i in block or i == metal:
            continue
        if not np.allclose(original[i], updated[i], atol=tol):
            return False
    return True


def apply_hapto_honest(
    symbols: Sequence[str],
    coords: np.ndarray,
    metal_idx: int = 0,
    donor_idxs: Optional[Sequence[int]] = None,
    md_cut: float = _DEFAULT_MD_CUT,
) -> Tuple[np.ndarray, int]:
    """Apply the hapto-honest corrector to ``coords`` for the metal at
    ``metal_idx``.  Returns ``(P_new, n_units_fixed)``.

    Default-OFF: returns ``(coords.copy(), 0)`` when
    :func:`honest_active` is ``False``.  In that case NO floating-point
    arithmetic happens — the output is bit-identical to the input.
    """
    P_orig = np.asarray(coords, dtype=float).reshape(-1, 3)
    if not honest_active():
        return P_orig.copy(), 0
    if metal_idx < 0 or metal_idx >= len(symbols):
        return P_orig.copy(), 0
    units = detect_hapto_units(symbols, P_orig, metal_idx, md_cut=md_cut)
    if not units:
        return P_orig.copy(), 0
    M = P_orig[metal_idx]
    P_new = P_orig.copy()
    n_fixed = 0
    # Surgical-fix (2026-06-04): when the optional global-clash gate is
    # active, snapshot the pre-correction collapse count so each
    # per-unit acceptance can verify "clash count does NOT increase".
    _clash_gate_active = _post_correction_clash_gate_active()
    _bcg = None
    _pre_collapse = 0
    if _clash_gate_active:
        try:
            from delfin.fffree import build_time_clash_gate as _bcg_mod
            _bcg = _bcg_mod
            _pre_collapse = int(_bcg.collapse_count(list(symbols), P_orig))
        except Exception:
            _bcg = None
    for u in units:
        block_idx = list(u.indices) + list(u.attached_h)
        P_try = _apply_unit_correction(P_new, M, u)
        if not _validate_rigid_block(P_new, P_try, block_idx):
            continue  # rigid invariant broken — skip this unit
        if not _validate_sigma_invariance(P_new, P_try, metal_idx, block_idx):
            continue  # something outside the block moved — skip
        # Sanity: metal must not move.
        if not np.allclose(P_try[metal_idx], P_new[metal_idx], atol=1e-9):
            continue
        # Sanity: the new centroid distance must equal the ideal within
        # 1e-6 Å (analytical translation; numerical safety net).
        new_cen = P_try[u.indices].mean(axis=0)
        new_dist = float(np.linalg.norm(new_cen - M))
        ideal = _ideal_m_centroid_distance(u.metal_symbol, u.eta)
        if abs(new_dist - ideal) > 1e-3:
            continue
        # Surgical-fix (2026-06-04, b00f9a0 forensik): post-correction
        # global-clash gate.  Snapping a ring centroid to its ideal
        # M-centroid distance may push ring carbons / attached H atoms
        # into nearby ligand atoms; that surfaces downstream as F20_h,
        # funcgrp, angles-violations on the otherwise-correct σ
        # sub-polyhedron.  When active, reject the unit-correction if
        # the global collapse count would INCREASE.  Default OFF ->
        # byte-identical to b00f9a0.
        if _clash_gate_active and _bcg is not None:
            try:
                _post_collapse = int(_bcg.collapse_count(list(symbols), P_try))
                if _post_collapse > _pre_collapse:
                    continue  # reject this unit-correction
                _pre_collapse = _post_collapse
            except Exception:
                # Defensive: never let the clash gate raise.
                pass
        # Accept.
        P_new = P_try
        n_fixed += 1
    return P_new, n_fixed
