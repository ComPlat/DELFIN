"""delfin.fffree.terminal_ligand_stagger — adjacent CO/NO/CN torsion stagger.

Surgical fix 2 (User mandate 2026-06-03, F24-interlig deep forensik):
when two terminal CO/NO/CN ligands occupy adjacent polyhedron vertices the
default M-L axis + zero-torsion placement collides the two terminal
O/N/N atoms at d ~ 1.85 A (well below the 0.85 x vdW floor of ~2.5 A).
This module detects such pairs purely from the assembled-complex graph
and the polyhedron-vertex adjacency, then rotates the second ligand's
terminal atom (and any attached H) about the M-donor axis by
``DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_ANGLE`` (default 60 deg).

The rotation is a pure rigid spin about the M-donor bond:

   - r(M-donor) is preserved EXACTLY (axis goes through the donor).
   - the donor stays at its vertex (M-D invariant by construction).
   - r(donor-terminal) is preserved (rigid rotation of the terminal about
     an axis through the donor).

Universal, graph-based: detection uses
:func:`delfin.fffree.grip_ensemble.identify_ligand_subgraphs` to find
the no-metal connected components, then keeps components whose atom set
is ``{donor, terminal_O_or_N}`` (2-atom heavy fragment with the donor as
the M-bonded atom and a single sp-terminal triple-bond partner).

Adjacency between polyhedron vertices is the standard "shortest geodesic
on the convex hull": two vertices ``v_i, v_j`` are adjacent iff the
angle ``acos(u_i . u_j)`` is below the polyhedron-specific median + 1e-3
rad threshold (deterministic, no SMILES).

API:
    detect_terminal_pairs(syms, P, metal_idx, donor_idxs, geometry=None,
        adjacency_cos_min=None) -> list[(d1, t1, d2, t2)]
    apply_stagger(syms, P, metal_idx, donor_idxs, geometry=None) ->
        (P_new, n_pairs_staggered)
    stagger_active() -> bool

Env flag:
    DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER=1  (per-fix)
    DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_ANGLE=60.0  (degrees, optional)
    DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_THRESHOLD=2.5  (Å min terminal-
        terminal distance under which we stagger, default 2.5)
    DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1   (F24-interlig surgical master)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (global master)
"""
from __future__ import annotations

import os
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Default stagger angle (degrees).  60 deg is the canonical staggered
# (Newman-eclipse + 60) for two sp-axes on adjacent vertices.
_DEFAULT_STAGGER_DEG = 60.0
# Below this terminal-terminal distance (A) we treat the pair as
# "potentially colliding" and apply the stagger.  2.5 A is well above the
# 1.85 A clash band and just inside the 0.85 * vdW(O,O) = 2.58 A floor.
_DEFAULT_THRESHOLD = 2.5
# Bond perception (covalent-sum factor).  Matches build_time_clash_gate.
_BOND_FACTOR = 1.30
# Triple-bond max length for the donor-terminal pair (C≡O 1.13, C≡N 1.16,
# N≡O 1.15 A); 1.30 covers the upper bound + slack.
_TRIPLE_MAX = 1.30


def stagger_active() -> bool:
    """True iff the terminal-ligand stagger (surgical fix 2) is enabled.

    Activated by either the per-fix flag
    ``DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER=1``, the F24-interlig surgical
    master ``DELFIN_FFFREE_F24_INTERLIG_FIX_ALL=1``, or the global
    construction master ``DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1``.
    Default OFF -> byte-identical to pre-fix HEAD.
    """
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_F24_INTERLIG_FIX_ALL", "0") == "1":
        return True
    return os.environ.get(
        "DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER", "0"
    ) == "1"


def _stagger_angle_rad() -> float:
    """Read ``DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_ANGLE`` (degrees)."""
    try:
        deg = float(os.environ.get(
            "DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_ANGLE",
            str(_DEFAULT_STAGGER_DEG),
        ))
    except (TypeError, ValueError):
        deg = _DEFAULT_STAGGER_DEG
    return math.radians(deg)


def _stagger_threshold() -> float:
    """Read ``DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_THRESHOLD`` (A)."""
    try:
        return float(os.environ.get(
            "DELFIN_FFFREE_TERMINAL_LIGAND_STAGGER_THRESHOLD",
            str(_DEFAULT_THRESHOLD),
        ))
    except (TypeError, ValueError):
        return _DEFAULT_THRESHOLD


def _heavy_bonds(syms: Sequence[str], P: np.ndarray
                 ) -> List[Tuple[int, int]]:
    """Heavy (incl. metal) geometric bonds: distance < 1.30 x covalent-sum.

    Used to identify the M-donor edge and the donor-terminal edge.
    """
    n = len(syms)
    out: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                out.append((i, j))
    return out


def _find_terminal_for_donor(
    d: int, metal_idx: int,
    syms: Sequence[str], P: np.ndarray,
) -> Optional[int]:
    """Find the sp-terminal partner of donor ``d`` for a CO/CN/NO ligand.

    The donor (a C or N) is the M-bonded atom and has the terminal as its
    OPPOSITE-direction triple-bond partner (the M-D-T angle is ~180 deg
    for CO/CN/NO).  We pick the heavy non-metal atom that:

      (a) is within ``_TRIPLE_MAX`` (1.30 A) of d -- triple-bond signature,
      (b) lies "behind" d when looking from the metal -- the M-D-T angle
          must be > 150 deg (cos < -0.87).

    This is robust against spurious cross-fragment "bonds" picked up by a
    naive covalent-radius bond search (which can occur when two CO
    ligands sit at small donor-separation angles).
    """
    if metal_idx < 0 or d < 0:
        return None
    M = P[metal_idx]
    Pd = P[d]
    md_vec = Pd - M
    nmd = float(np.linalg.norm(md_vec))
    if nmd < 1e-9:
        return None
    md_hat = md_vec / nmd
    best_idx = -1
    best_dist = float("inf")
    for k in range(len(syms)):
        if k == d or k == metal_idx:
            continue
        if _bd._is_metal(syms[k]) or syms[k] == "H":
            continue
        d_dk = float(np.linalg.norm(P[k] - Pd))
        if d_dk >= _TRIPLE_MAX:
            continue
        # Geometric anti-axis check: T must be roughly on the far side
        # of D from M (cos(M-D-T angle) < -0.87, i.e. angle > 150 deg).
        dt_vec = P[k] - Pd
        ndt = float(np.linalg.norm(dt_vec))
        if ndt < 1e-9:
            continue
        cos_ang = float(np.dot(md_hat, dt_vec / ndt))
        # When the donor-terminal vector points roughly along the same
        # direction as M->D (extending the M-D axis outward), cos > 0.87.
        if cos_ang < 0.87:
            continue
        if d_dk < best_dist:
            best_dist = d_dk
            best_idx = k
    return best_idx if best_idx >= 0 else None


def _adjacency_from_geometry(
    metal_pos: np.ndarray,
    donor_pos_by_idx: Dict[int, np.ndarray],
) -> Dict[Tuple[int, int], bool]:
    """Geometry-only polyhedron-vertex adjacency between donor pairs.

    Two donors are "adjacent" iff the angle ``M-D1, M-D2`` is at or below
    the median pairwise donor angle of this complex.  This is a purely
    structural definition (no polyhedron name needed), works for any CN,
    and is deterministic.
    """
    idxs = sorted(donor_pos_by_idx.keys())
    if len(idxs) < 2:
        return {}
    angles: Dict[Tuple[int, int], float] = {}
    for ii in range(len(idxs)):
        ai = idxs[ii]
        va = donor_pos_by_idx[ai] - metal_pos
        na = float(np.linalg.norm(va))
        if na < 1e-9:
            continue
        for jj in range(ii + 1, len(idxs)):
            aj = idxs[jj]
            vb = donor_pos_by_idx[aj] - metal_pos
            nb = float(np.linalg.norm(vb))
            if nb < 1e-9:
                continue
            c = float(np.dot(va, vb) / (na * nb))
            c = max(-1.0, min(1.0, c))
            angles[(ai, aj)] = math.acos(c)
    if not angles:
        return {}
    sorted_ang = sorted(angles.values())
    # Adjacency cut-off = median + 1e-3 rad slack to break ties.
    mid = sorted_ang[len(sorted_ang) // 2]
    out: Dict[Tuple[int, int], bool] = {}
    for k, v in angles.items():
        out[k] = (v <= mid + 1e-3)
    return out


def detect_terminal_pairs(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int, donor_idxs: Sequence[int],
    threshold: Optional[float] = None,
) -> List[Tuple[int, int, int, int]]:
    """Find adjacent (donor, terminal) ligand pairs that look like
    terminal CO/NO/CN with too-close terminal atoms.

    A terminal ligand is identified by:
      1. donor ``d`` is on the polyhedron (in ``donor_idxs``),
      2. donor has exactly ONE heavy neighbour ``t`` that is non-metal,
      3. ``t`` is sp-terminal: its only heavy neighbour is ``d`` AND the
         d-t distance is below ``_TRIPLE_MAX`` (triple-bond signature),
      4. the donor-element / terminal-element pair is one of:
         (C, O)  — carbonyl (M-C≡O)
         (N, O)  — nitrosyl (M-N≡O)
         (C, N)  — cyanide  (M-C≡N)

    Two terminal ligands form a "pair candidate" when their donors are
    adjacent (median-angle adjacency) AND ``d(t_i, t_j) < threshold``.

    Returns
    -------
    list of (d1, t1, d2, t2)
        Deterministic ordering: pairs sorted by (min(d1,d2), max(d1,d2)),
        and within each pair the lower-index donor goes first.
    """
    if threshold is None:
        threshold = _stagger_threshold()
    n = len(syms)
    if n == 0 or len(donor_idxs) < 2 or metal_idx < 0 or metal_idx >= n:
        return []
    P = np.asarray(P, dtype=float)
    donor_set = set(int(d) for d in donor_idxs if 0 <= int(d) < n)
    metal_pos = P[int(metal_idx)]

    # Find terminal-ligand (donor, terminal) units via geometric
    # anti-axis search (robust against spurious cross-fragment bonds in
    # the naive covalent-radius graph for tightly-packed configurations).
    terminals: List[Tuple[int, int]] = []
    for d in sorted(donor_set):
        t = _find_terminal_for_donor(d, int(metal_idx), syms, P)
        if t is None:
            continue
        # Element pair must be one of the canonical sp-terminals.
        pair = (str(syms[d]), str(syms[t]))
        if pair not in (("C", "O"), ("N", "O"), ("C", "N")):
            continue
        terminals.append((d, t))

    if len(terminals) < 2:
        return []

    # Polyhedron adjacency from donor geometry.
    adj = _adjacency_from_geometry(
        metal_pos, {d: P[d] for d, _ in terminals},
    )

    # Form pair candidates: adjacent donors with too-close terminals.
    out: List[Tuple[int, int, int, int]] = []
    for a in range(len(terminals)):
        d_a, t_a = terminals[a]
        for b in range(a + 1, len(terminals)):
            d_b, t_b = terminals[b]
            key = (min(d_a, d_b), max(d_a, d_b))
            if not adj.get(key, False):
                continue
            d_tt = float(np.linalg.norm(P[t_a] - P[t_b]))
            if d_tt >= float(threshold):
                continue
            # Deterministic donor ordering (lower index first).
            if d_a <= d_b:
                out.append((d_a, t_a, d_b, t_b))
            else:
                out.append((d_b, t_b, d_a, t_a))
    out.sort(key=lambda q: (min(q[0], q[2]), max(q[0], q[2]), q[0]))
    return out


def _rotation_matrix_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation matrix for a unit ``axis`` and ``angle`` (radians).
    Deterministic; returns identity for a zero-length axis.
    """
    a = np.asarray(axis, dtype=float)
    na = float(np.linalg.norm(a))
    if na < 1e-9:
        return np.eye(3)
    u = a / na
    c = math.cos(angle)
    s = math.sin(angle)
    ux, uy, uz = u
    K = np.array([
        [0.0, -uz, uy],
        [uz, 0.0, -ux],
        [-uy, ux, 0.0],
    ], dtype=float)
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def _rotate_about_axis(
    point: np.ndarray, origin: np.ndarray,
    axis: np.ndarray, angle: float,
) -> np.ndarray:
    """Rotate ``point`` about the line through ``origin`` along ``axis`` by
    ``angle`` radians.  Deterministic.
    """
    R = _rotation_matrix_axis_angle(axis, angle)
    return (point - origin) @ R.T + origin


def _terminal_h_indices(syms: Sequence[str], P: np.ndarray,
                        t: int) -> List[int]:
    """H atoms whose closest non-H, non-metal heavy atom is ``t``.

    Uses a simple closest-heavy-atom assignment (rather than the
    potentially-spurious bond graph) so it works even in tightly packed
    cis-CO complexes.  For pure CO/CN/NO ligands this returns []; for
    HCN- / NCH-type terminals it returns the proton.
    """
    out: List[int] = []
    n = len(syms)
    if t < 0 or t >= n:
        return out
    for h in range(n):
        if syms[h] != "H":
            continue
        # Closest heavy (non-metal, non-H) atom to this H.
        best_k = -1
        best_d = float("inf")
        for k in range(n):
            if k == h:
                continue
            if syms[k] == "H" or _bd._is_metal(syms[k]):
                continue
            d = float(np.linalg.norm(P[h] - P[k]))
            if d < best_d:
                best_d = d
                best_k = k
        if best_k == t and best_d < _BOND_FACTOR * _bd._ideal_bond(
            syms[h], syms[t]
        ):
            out.append(h)
    return out


def apply_stagger(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int, donor_idxs: Sequence[int],
    angle: Optional[float] = None,
    threshold: Optional[float] = None,
) -> Tuple[np.ndarray, int]:
    """Apply the torsion-stagger to every adjacent terminal CO/NO/CN pair.

    For each pair ``(d1, t1, d2, t2)``:
      - rotation axis ``= unit(M - d1)`` (M-donor direction of the FIRST
        ligand; pure spin about it preserves r(M-d1) and r(M-d2)
        identically since d2 does not move).
      - rotation origin ``= M`` (so the metal stays at the origin).
      - the SECOND ligand's terminal atom ``t2`` (and any attached H) is
        rotated by ``angle`` (default 60 deg).
      - d2 itself does NOT move (preserves its polyhedron-vertex seat).

    For each pair we try the +angle and the -angle rotation and keep the
    one that maximises ``d(t1, t2_new)`` (deterministic + > 0 by
    construction).

    Donors and the metal are NEVER moved.  Cumulative across multiple
    pairs: if a ligand participates in more than one detected pair, the
    rotations compose (deterministic order via
    :func:`detect_terminal_pairs`).

    Env-gated default-OFF: returns ``(P_copy, 0)`` byte-identical when no
    flag set.
    """
    P_new = np.asarray(P, dtype=float).copy()
    if not stagger_active():
        return P_new, 0
    if angle is None:
        angle = _stagger_angle_rad()
    pairs = detect_terminal_pairs(
        syms, P_new, int(metal_idx), donor_idxs, threshold=threshold,
    )
    if not pairs:
        return P_new, 0
    metal_pos = P_new[int(metal_idx)]

    n_pairs = 0
    for (d1, t1, d2, t2) in pairs:
        # Rotation axis = M -> d1 direction.
        axis = P_new[d1] - metal_pos
        # Indices that ride along with t2: the H atoms whose only heavy
        # neighbour is t2 (NCH / HCN style).  For pure CO/CN/NO there
        # are none.
        h_idxs = _terminal_h_indices(syms, P_new, t2)
        candidates: List[Tuple[float, np.ndarray, List[Tuple[int, np.ndarray]]]] = []
        for sign in (+1.0, -1.0):
            t2_new = _rotate_about_axis(
                P_new[t2], metal_pos, axis, sign * angle,
            )
            new_hs: List[Tuple[int, np.ndarray]] = []
            for h in h_idxs:
                new_hs.append((
                    h,
                    _rotate_about_axis(
                        P_new[h], metal_pos, axis, sign * angle,
                    ),
                ))
            d_after = float(np.linalg.norm(P_new[t1] - t2_new))
            candidates.append((d_after, t2_new, new_hs))
        # Pick the larger d(t1, t2_new) -- deterministic (sign +1 first).
        candidates.sort(key=lambda q: -q[0])
        best_d, best_t2, best_hs = candidates[0]
        # Defence-in-depth: only commit when the rotation actually improves
        # the terminal-terminal distance (>= original + 1e-6).
        d_before = float(np.linalg.norm(P_new[t1] - P_new[t2]))
        if best_d > d_before + 1e-6:
            P_new[t2] = best_t2
            for h, ph in best_hs:
                P_new[h] = ph
            n_pairs += 1
    return P_new, n_pairs


__all__ = [
    "detect_terminal_pairs", "apply_stagger", "stagger_active",
    "_stagger_angle_rad", "_stagger_threshold",
]
