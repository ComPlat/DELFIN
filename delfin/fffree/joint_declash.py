"""delfin.fffree.joint_declash — JOINT global INTER-LIGAND heavy-heavy declash.

The deepest FF-free recall lever for the "class-B" bulk: large / multi-ligand
complexes whose COORDINATION CORE is already IDEAL (donor-angle-RMSD ~0-8°, sane
M-D) but which fall to the distorted legacy-UFF fallback purely because the
ligand BODIES clash — the self-gate ``_build_is_clean`` rejects on genuine
INTER-LIGAND heavy-heavy overlaps (~1.9-2.2 A) even though the first coordination
shell is perfect.

Why a dedicated pass on top of #308 ``torsion_relax``
----------------------------------------------------
#308 minimises a SINGLE global clash SUM that is H-inclusive and counts BOTH
intra- and inter-ligand contacts equally.  On a crowded class-B complex the H-H
and intra-ligand terms DOMINATE that sum, so the coordinate descent spends its
moves relieving cheap H-H / intra clashes and the never-worse-on-MIN floor is set
by an H-H pair — leaving the load-bearing HEAVY-HEAVY inter-ligand overlap (the
ONLY thing the self-gate actually rejects on) under-optimised.  This module makes
the objective the thing the gate measures: a GLOBAL INTER-LIGAND HEAVY-HEAVY
clash sum (H terms kept only as a light secondary tie-breaker), so every move is
spent opening the contacts that block the gate.

Degrees of freedom (identical kinematics to #308 -> identical safety property)
------------------------------------------------------------------------------
For every ligand, jointly:
  * its WHOLE-BODY RIGID ROTATION about the M-donor axis (anchor = metal, pivot =
    donor; the donor lies ON the axis so the M-D distance is exactly preserved),
  * its INTERNAL rotatable-bond torsions (rigid distal sub-tree about each
    non-ring single-bond axis).
Both are pure rigid rotations of a sub-tree about a bond axis: they change ONLY
the dihedral / orientation, never a bond length or a bond angle.  The metal +
ALL donor atoms are FROZEN, so the coordination polyhedron is invariant by
construction (proven in the test-suite: bond-length + bond-angle RMSD before/
after == 0).

Method
------
Internal-coordinate coordinate descent over the joint DOF set (whole-ligand M-D
rotations FIRST — they move the most mass and relieve the inter-ligand overlap
most directly — then internal torsions), each DOF swept on a fixed deterministic
angular grid to its objective-minimising angle.  Accept-if-better (strictly
decreasing), never-worse-on-MIN floor on the inter-ligand heavy minimum, hard
M-D invariant guard (<= md_tol, default 0.05 A; any violation rolls the whole
relaxation back to the input).  Deterministic: fixed DOF order, fixed grid, no
RNG.  Never returns non-finite; on any exception the input frame is returned.

Integration
-----------
Env-gated ``DELFIN_FFFREE_JOINT_DECLASH`` (default ``"0"`` => byte-identical: the
pass is never invoked AND the coupled decompose gate-lift stays at 8).  Runs as a
post-placement pass on each assembled frame, AFTER #308 (if also on) and BEFORE
the self-gate, so a declashed class-B build now PASSES ``_build_is_clean``.

License: open-source / pure geometry only (Bondi-style vdW radii reused from the
FF-free refiner).  No CSD/CCDC data.
"""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd
from delfin.fffree import torsion_relax as _TR
from delfin.fffree.refine import _vdw

# Geometric clash factor — identical to the self-gate / #308 / the spec's f.
_CLASH_F = 0.75

# H-H and X-H contacts kept as a LIGHT secondary tie-breaker: the self-gate
# rejects on HEAVY-HEAVY only, so heavy-heavy must dominate the objective.
_H_WEIGHT = 0.05

# Default coordinate-descent controls (env-overridable; bounded + deterministic).
_DEF_GRID = 24          # angular grid steps per DOF
_DEF_PASSES = 8         # max coordinate-descent passes over all DOFs
_DEF_MD_TOL = 0.05      # A, hard M-D invariant guard
_DEF_MAX_DOFS = 96      # cap DOFs optimised (bulkiest first; whole-ligand spins kept)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_JOINT_DECLASH", "0") == "1"


# ---------------------------------------------------------------------------
# Per-atom ligand membership (so we can score INTER-ligand contacts only)
# ---------------------------------------------------------------------------


def _ligand_of_atom(n: int, syms: Sequence[str],
                    bond_pairs: Optional[Sequence[Tuple[int, int]]],
                    P: np.ndarray) -> np.ndarray:
    """Per-atom ligand id (0..k-1) for the assembled complex; the metal gets -1.

    The ligand graph is the bond graph with the metal-donor coordination bonds
    REMOVED — each connected component on the non-metal atoms is one ligand.
    Using the threaded ``bond_pairs`` (true connectivity) is strongly preferred so
    two interpenetrating ligands at a fortuitous bonding distance are not fused
    into one component (which would hide their inter-ligand clash from the
    objective).  Falls back to geometric perception when ``bond_pairs`` is None.
    """
    adj, _bonds = _TR._adjacency(syms, P, bond_pairs)
    metals = {i for i in range(n) if _bd._is_metal(syms[i])}
    lig = np.full(n, -1, dtype=int)
    cur = 0
    for start in range(n):
        if start in metals or lig[start] != -1:
            continue
        # BFS over non-metal atoms only (coordination bonds to the metal are not
        # traversed -> ligands stay separate components).
        stack = [start]
        lig[start] = cur
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b in metals or lig[b] != -1:
                    continue
                lig[b] = cur
                stack.append(b)
        cur += 1
    return lig


# ---------------------------------------------------------------------------
# Inter-ligand clash objective (heavy-heavy dominant, H light tie-breaker)
# ---------------------------------------------------------------------------


def _inter_mask(syms: Sequence[str], lig: np.ndarray,
                excl: Sequence[Set[int]]) -> Tuple[np.ndarray, np.ndarray]:
    """Two upper-triangular (n,n) boolean masks of counted pairs:

      * ``heavy``: i<j, DIFFERENT ligands, neither H, neither metal, not 1-2/1-3
        (the self-gate-relevant inter-ligand HEAVY-HEAVY contacts), and
      * ``light``: i<j, DIFFERENT ligands, at least one H, not 1-2/1-3 (kept only
        as a small secondary tie-breaker so the descent does not introduce gross
        H clashes while opening heavy ones).

    Intra-ligand contacts are NOT counted (this is the *inter-ligand* declash;
    intra-ligand strain is #308's / the refiner's job and is geometry-fixed by the
    rigid sub-tree kinematics anyway).  Pairs sharing the metal (i.e. a ligand
    atom vs the metal) are excluded — the metal is frozen and on the axis.
    """
    n = len(syms)
    heavy = np.zeros((n, n), dtype=bool)
    light = np.zeros((n, n), dtype=bool)
    is_h = np.array([s == "H" for s in syms])
    is_m = np.array([_bd._is_metal(s) for s in syms])
    iu, ju = np.triu_indices(n, k=1)
    for k in range(len(iu)):
        i, j = int(iu[k]), int(ju[k])
        if is_m[i] or is_m[j]:
            continue
        if lig[i] == lig[j]:                      # same ligand (or both -1) -> intra
            continue
        if lig[i] < 0 or lig[j] < 0:              # safety: unassigned -> skip
            continue
        if j in excl[i]:                          # 1-2/1-3 (cannot happen inter, but safe)
            continue
        if is_h[i] or is_h[j]:
            light[i, j] = True
        else:
            heavy[i, j] = True
    return heavy, light


def _objective(P: np.ndarray, heavy: np.ndarray, light: np.ndarray,
               rsum: np.ndarray) -> Tuple[float, float]:
    """Return ``(L, heavy_dmin)``.

    ``L = sum_{heavy pairs} over^2 + _H_WEIGHT * sum_{light pairs} over^2`` where
    ``over = max(0, f*(vdw_i+vdw_j) - d)``.  ``heavy_dmin`` is the minimum
    inter-ligand HEAVY-HEAVY distance (the quantity the self-gate gates on).
    """
    diff = P[:, None, :] - P[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    over_h = np.where(heavy, rsum - dist, 0.0)
    over_h = np.where(over_h > 0.0, over_h, 0.0)
    over_l = np.where(light, rsum - dist, 0.0)
    over_l = np.where(over_l > 0.0, over_l, 0.0)
    loss = float((over_h * over_h).sum()) + _H_WEIGHT * float((over_l * over_l).sum())
    hd = dist[heavy]
    hdmin = float(hd.min()) if hd.size else float("inf")
    return loss, hdmin


# ---------------------------------------------------------------------------
# Core: joint inter-ligand coordinate-descent declash
# ---------------------------------------------------------------------------


def declash(syms: Sequence[str], P, frozen: Iterable[int],
            grid: int = _DEF_GRID, passes: int = _DEF_PASSES,
            md_tol: float = _DEF_MD_TOL, max_dofs: int = _DEF_MAX_DOFS,
            bond_pairs: Optional[Sequence[Tuple[int, int]]] = None) -> np.ndarray:
    """Joint global INTER-LIGAND heavy-heavy declash of an assembled complex.

    ``frozen``: indices that must not move (metal + all donor atoms).
    ``bond_pairs`` (optional, strongly preferred): the TRUE connectivity threaded
    from the builder (see :func:`torsion_relax._adjacency`) — both for robust DOF
    identification and for correct ligand-membership partition.

    DOFs: each ligand's whole-body rotation about its M-D axis PLUS its internal
    rotatable-bond torsions (reusing #308's :func:`identify_dofs`, which already
    treats the M-D bond as a rotatable axis with the metal on the anchor side).
    The whole-ligand M-D spins are floated FIRST (they move the most mass and
    relieve the inter-ligand overlap most directly), then internal torsions.

    Returns the relaxed coordinate array.  Accept-if-better, never-worse-on
    heavy-MIN, hard M-D invariant; deterministic; never raises (input on failure).
    """
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return np.array(P, dtype=float)
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return P0
    frozen_set = set(int(x) for x in frozen)

    dofs = _TR.identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                             bond_pairs=bond_pairs)
    if not dofs:
        return P0

    # Reorder: whole-ligand M-D spins (anchor == metal) first, then internal
    # torsions; both already sorted bulkiest-first within #308's identify_dofs,
    # so this is a stable partition preserving determinism.
    metals = {i for i in range(n) if _bd._is_metal(syms[i])}
    md_spins = [d for d in dofs if d["anchor"] in metals]
    internal = [d for d in dofs if d["anchor"] not in metals]
    ordered = md_spins + internal

    adj, _bonds = _TR._adjacency(syms, P0, bond_pairs)
    excl = _TR._excl_1_2_3(adj, n)
    lig = _ligand_of_atom(n, syms, bond_pairs, P0)
    heavy, light = _inter_mask(syms, lig, excl)
    if not heavy.any() and not light.any():
        return P0                                     # nothing inter-ligand to declash

    base_md = _TR._md_pairs(syms, P0)
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[:, None] + vdw[None, :])

    Pcur = P0.copy()
    best_loss, base_hdmin = _objective(Pcur, heavy, light, rsum)
    if best_loss <= 1e-12:
        return Pcur
    # never-worse-on-heavy-MIN floor: a move is rejected if it would push the worst
    # inter-ligand HEAVY contact below the input frame's worst (a SUM objective
    # could otherwise relieve several mild clashes by crushing one comfortable
    # heavy pair -> lower L but a WORSE self-gate-relevant minimum).
    hdmin_floor = base_hdmin - 1e-6

    angles = [2.0 * math.pi * k / float(grid) for k in range(grid)]
    for _ in range(max(1, passes)):
        improved = False
        for dof in ordered:
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            rot = dof["rotating"]
            origin = Pcur[anchor]
            axis = Pcur[pivot] - Pcur[anchor]
            if float(np.linalg.norm(axis)) < 1e-9:
                continue
            local_best_loss = best_loss
            local_best_P = None
            for ang in angles:
                if abs(ang) < 1e-12:
                    continue
                trial = _TR._rotate_subtree(Pcur, origin, axis, ang, rot)
                if not np.all(np.isfinite(trial)):
                    continue
                if not _TR._md_ok(base_md, trial, md_tol):
                    continue
                tl, thd = _objective(trial, heavy, light, rsum)
                if thd < hdmin_floor:
                    continue                          # would worsen worst heavy contact
                if tl < local_best_loss - 1e-12:
                    local_best_loss = tl
                    local_best_P = trial
            if local_best_P is not None:
                Pcur = local_best_P
                best_loss = local_best_loss
                improved = True
            if best_loss <= 1e-12:
                break
        if not improved or best_loss <= 1e-12:
            break

    # final hard guard: roll the whole declash back to the input on any M-D
    # violation / non-finite / worse-than-input loss or heavy-min (never-worse).
    if not np.all(np.isfinite(Pcur)) or not _TR._md_ok(base_md, Pcur, md_tol):
        return P0
    fin_loss, fin_hdmin = _objective(Pcur, heavy, light, rsum)
    if fin_loss > best_loss + 1e-9 or fin_hdmin < hdmin_floor:
        return P0
    return Pcur


def declash_if_enabled(syms: Sequence[str], P, frozen: Iterable[int],
                       bond_pairs: Optional[Sequence[Tuple[int, int]]] = None):
    """Wire-in entry: when ``DELFIN_FFFREE_JOINT_DECLASH=1`` run :func:`declash`,
    else return ``P`` unchanged (byte-identical default-OFF).  ``bond_pairs``
    (optional) is the true connectivity from the builder.  Never raises."""
    if not _enabled():
        return P
    try:
        grid = _TR._env_int("DELFIN_FFFREE_JOINT_GRID", _DEF_GRID, 4, 72)
        passes = _TR._env_int("DELFIN_FFFREE_JOINT_PASSES", _DEF_PASSES, 1, 32)
        md_tol = _TR._env_float("DELFIN_FFFREE_JOINT_MD_TOL", _DEF_MD_TOL, 0.0, 2.0)
        max_dofs = _TR._env_int("DELFIN_FFFREE_JOINT_MAX_DOFS", _DEF_MAX_DOFS, 1, 256)
        return declash(syms, P, frozen, grid=grid, passes=passes,
                       md_tol=md_tol, max_dofs=max_dofs, bond_pairs=bond_pairs)
    except Exception:
        return P
