"""delfin.manta.torsion_relax — FF-free whole-complex TORSION-SPACE clash relax.

When rigid M-D-axis rotation (#307 inter-ligand frame selection) is not enough —
when LIGAND-INTERNAL rotation (P-C, then C-C, ...) is needed — this module
optimizes ALL rotatable single bonds of the ASSEMBLED complex SIMULTANEOUSLY
(one joint dihedral vector) against a purely geometric clash objective.  NO force
field, no atom types, no charges — only van-der-Waals overlap penalties.

Why it cannot distort the internal geometry (the key safety property the
whole-DG-embed #82/#100 violated):  a torsion = rigid rotation of the distal
sub-tree about a bond axis.  It changes ONLY the dihedral angle; every bond
length and every bond angle is preserved EXACTLY (internal-coordinate
kinematics).  Bonds and angles are not optimization variables, so the optimizer
cannot produce the distortion that loose 1-2/1-3 DG bounds did.

DOF
---
Every rotatable single bond (heavy-heavy, NOT in a ring), INCLUDING the
M-donor bond (= whole-ligand rotation about the M-D axis, the continuous
generalisation of #307's discrete spin).  For each rotatable bond, the side that
does NOT contain the metal + the donor atoms is the moving distal sub-tree.

Frozen
------
The coordination core — metal (index 0) and all donor atoms — stays fixed
(donor on its vertex, M-D distance + polyhedron preserved).  Rings are NEVER
torsion-rotated (that would break the ring); ring flex stays with Cremer-Pople
elsewhere.

Objective
---------
``L(theta) = sum_{i<j, non-bonded, not 1-2/1-3} max(0, f*(vdw_i+vdw_j) - d_ij)^2``
H-inclusive (the AVECUA-critical clash was H-H 0.73 A), inter- AND intra-ligand,
``f = 0.75`` (the same factor + the same ``_vdw`` table the rest of the FF-free
builder's ``_clash_count`` uses).

Method
------
Internal-coordinate coordinate descent: cycle over every rotatable bond in a
fixed deterministic order; rotate each distal sub-tree about its axis to the
L-minimizing angle found on a fixed deterministic grid; iterate to convergence /
max passes.  Accept-if-better (a move is taken only if L strictly decreases —
never worse).  M-D-invariant guard (<= 0.05 A) else the whole relaxation rolls
back to the input frame.  Deterministic: fixed bond order, fixed grid, no RNG.

Integration
-----------
Env-gated ``DELFIN_FFFREE_TORSION_RELAX`` (default ``"0"`` = byte-identical, the
relax is never invoked).  When ``"1"`` it runs as a post-placement step on each
assembled frame in the ensemble / heteroleptic paths of ``assemble_complex.py``
(after #307 frame selection, around ``refine``).  Never returns non-finite; on
any exception the un-relaxed frame is returned unchanged.

License: open-source / pure geometry only (Bondi-style vdW radii reused from the
FF-free refiner).  No CSD/CCDC data.
"""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd
from delfin.manta.refine import _vdw

# Geometric clash factor — identical to the builder's _clash_count / the spec's f.
_CLASH_F = 0.75

# Default coordinate-descent controls (env-overridable; all bounded + deterministic).
_DEF_GRID = 18          # angular grid steps per rotatable bond (12-24 per spec)
_DEF_PASSES = 6         # max coordinate-descent passes over all bonds
_DEF_MD_TOL = 0.05      # A, hard M-D invariant guard
_DEF_MAX_DOFS = 64      # cap rotatable bonds optimized (bulkiest first)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_TORSION_RELAX", "0") == "1"


def _env_int(name: str, default: int, lo: int, hi: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, int(raw)))
    except (TypeError, ValueError):
        return default


def _env_float(name: str, default: float, lo: float, hi: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, float(raw)))
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Geometric graph (no SMILES, no Open Babel — matches the FF-free philosophy)
# ---------------------------------------------------------------------------


# A heavy-heavy contact this much shorter than its ideal bond is an OVERLAP, not a
# bond — exactly the gross-overlap floor the FF-free self-gate uses (_build_is_clean
# rejects heavy-heavy < 0.60*ideal).  Perceiving such an overlap as a "bond" would put
# the two interpenetrating atoms into each other's 1-2/1-3 exclusion set, so the clash
# objective would IGNORE the very overlap it must relieve (masking bug).  Keeping it a
# non-bond lets the loss penalise it and lets a torsion try to open it up.
_BOND_LO_FRAC = 0.60


def _adjacency(syms: Sequence[str], P: np.ndarray,
               bond_pairs: Optional[Sequence[Tuple[int, int]]] = None
               ) -> Tuple[List[List[int]], List[Tuple[int, int]]]:
    """Bond graph of the assembled complex: heavy-heavy + X-H + M-donor.

    ``bond_pairs`` (preferred): the TRUE connectivity threaded from the builder's
    ligand mols (atom order preserved through assembly), plus the M-donor pairs.
    Using it avoids the fundamental failure of pure distance perception on a
    crowded complex — two atoms from DIFFERENT ligands at a fortuitous bonding
    distance (e.g. interpenetrating PMe3 methyls ~1.5 A) would otherwise be read
    as a covalent bond, which (a) puts the clashing pair into each other's 1-2/1-3
    exclusion set so the loss IGNORES it, and (b) makes the M-D coordination bonds
    look like ring bonds (cross-ligand cycle) so whole-ligand M-D rotation is lost.

    When ``bond_pairs`` is None we fall back to geometric perception
    (``_bd._geometric_bonds`` + auto M-D), with one safety: a heavy-heavy pair
    closer than ``_BOND_LO_FRAC*ideal`` is an OVERLAP, not a bond, and is dropped
    (so the loss penalises it).  X-H pairs keep the full window (X-ray-short X-H
    ~0.6-0.95 A are real).
    """
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    bonds: List[Tuple[int, int]] = []
    seen: Set[Tuple[int, int]] = set()

    def _add(i: int, j: int) -> None:
        if i == j or i < 0 or j < 0 or i >= n or j >= n:
            return
        key = (min(i, j), max(i, j))
        if key in seen:
            return
        seen.add(key)
        bonds.append(key)
        adj[i].append(j)
        adj[j].append(i)

    if bond_pairs is not None:
        for i, j in bond_pairs:
            _add(int(i), int(j))
    else:
        # organic bonds (heavy-heavy + X-H), metal pairs excluded by the helper;
        # reject heavy-heavy OVERLAPS so they stay clashes, not 1-2 neighbours.
        for i, j in _bd._geometric_bonds(syms, P):
            if syms[i] != "H" and syms[j] != "H":
                d = float(np.linalg.norm(P[i] - P[j]))
                if d < _BOND_LO_FRAC * _bd._ideal_bond(syms[i], syms[j]):
                    continue
            _add(i, j)
    # coordination bonds: metal -> first-shell heavy atoms (added in both modes so a
    # threaded bond list need not pre-include them; idempotent via the seen-set).
    for i in range(n):
        if not _bd._is_metal(syms[i]):
            continue
        for j in range(n):
            if i == j or syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _bd._ideal_bond(syms[i], syms[j]):
                _add(i, j)
    return adj, bonds


def _ring_bonds(adj: Sequence[Sequence[int]],
                bonds: Sequence[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    """Return the set of bonds (i<j) that lie on at least one cycle.

    A bond (i,j) is in a ring iff j is still reachable from i after the bond is
    removed.  Deterministic BFS; cheap for molecular sizes.  Coordination bonds
    are handled the same way, so a chelate's M-D bonds (which close a metallacycle
    through the backbone) are correctly recognised as ring bonds and excluded as
    rotatable axes — exactly what we want (the chelate core must stay rigid).
    """
    bondset = {(min(a, b), max(a, b)) for a, b in bonds}
    ring: Set[Tuple[int, int]] = set()
    for (i, j) in bondset:
        # BFS from i to j without using the direct i-j edge
        stack = [i]
        seen = {i}
        reached = False
        while stack:
            cur = stack.pop()
            for nb in adj[cur]:
                if cur == i and nb == j:
                    continue          # forbid the direct edge in one direction
                if cur == j and nb == i:
                    continue          # forbid the direct edge in the other
                if nb in seen:
                    continue
                if nb == j:
                    reached = True
                    stack = []
                    break
                seen.add(nb)
                stack.append(nb)
        if reached:
            ring.add((i, j))
    return ring


def _distal_subtree(adj: Sequence[Sequence[int]], anchor: int, pivot: int
                    ) -> List[int]:
    """Atoms on the *pivot* side of the bond anchor-pivot (BFS from pivot, not
    crossing back through anchor).  Returns pivot + everything distal to it."""
    visited = {anchor, pivot}
    out = [pivot]
    stack = [pivot]
    while stack:
        cur = stack.pop()
        for nb in adj[cur]:
            if nb in visited:
                continue
            visited.add(nb)
            out.append(nb)
            stack.append(nb)
    return out


# ---------------------------------------------------------------------------
# Rotatable-DOF identification
# ---------------------------------------------------------------------------


def identify_dofs(syms: Sequence[str], P: np.ndarray, frozen: Set[int],
                  max_dofs: int = _DEF_MAX_DOFS,
                  bond_pairs: Optional[Sequence[Tuple[int, int]]] = None
                  ) -> List[dict]:
    """Identify rotatable single bonds of the assembled complex.

    A bond (a,b) is a rotatable axis when:
      * it is NOT in a ring (rings stay rigid),
      * the side NOT containing the frozen core (metal + donors) is non-empty
        and contains at least one atom that actually moves under rotation
        (heavy atom OR >= 3 H),
      * after fixing the orientation so the frozen-core side is the ANCHOR side,
        the anchor end has at least one heavy reference neighbour (else the
        rotation axis is cylindrically degenerate, e.g. a terminal C-H).

    The M-donor bond qualifies (whole-ligand rotation about the M-D axis); the
    metal is always on the anchor side so the metal + the donor stay put.

    Returns a list of dicts: ``{"anchor", "pivot", "rotating", "score"}`` sorted
    bulkiest-rotating-group first, then by (anchor, pivot) for determinism.
    """
    n = len(syms)
    if n == 0:
        return []
    adj, bonds = _adjacency(syms, P, bond_pairs)
    ring = _ring_bonds(adj, bonds)
    frozen = set(frozen)
    metals = {i for i in range(n) if _bd._is_metal(syms[i])}

    dofs: List[dict] = []
    for (a, b) in bonds:
        if (a, b) in ring:
            continue
        # The two rigid halves split by the bond.  ANCHOR must be the half that
        # carries the metal (so the metal stays put); the other half (PIVOT side)
        # rotates about the bond axis.
        side_a = set(_distal_subtree(adj, b, a))   # a-side (b removed = anchor cut)
        side_b = set(_distal_subtree(adj, a, b))   # b-side
        a_has_metal = bool(side_a & metals)
        b_has_metal = bool(side_b & metals)
        if a_has_metal and b_has_metal:
            continue                    # bond is INSIDE the metal core -> never rotate
        if a_has_metal:
            anchor, pivot = a, b
        elif b_has_metal:
            anchor, pivot = b, a
        else:
            # No metal on either side (safety branch; a connected complex always has
            # the metal on one side): rotate the smaller half deterministically.
            if len(side_a) <= len(side_b):
                anchor, pivot = b, a
            else:
                anchor, pivot = a, b

        rotating = _distal_subtree(adj, anchor, pivot)
        rot_set = set(rotating)
        # A frozen atom (metal / donor) in the rotating half is only acceptable if it
        # is the PIVOT itself: the pivot lies ON the rotation axis, so it does not
        # move.  This is exactly the M-donor whole-ligand-rotation case (anchor=metal,
        # pivot=donor; donor frozen but on-axis -> the M-D distance is preserved).
        # Any OTHER frozen atom in the rotating half WOULD move -> reject.
        if (rot_set - {pivot}) & frozen:
            continue
        # movement filter: a heavy atom moves, or it is a >=3 H rotor (methyl)
        rot_heavy = sum(1 for idx in rotating if syms[idx] != "H")
        rot_h = sum(1 for idx in rotating if syms[idx] == "H")
        rot_heavy_distal = rot_heavy - 1            # exclude the pivot itself
        if rot_heavy_distal <= 0 and rot_h < 3:
            continue
        # anchor needs a heavy reference neighbour other than the pivot, else the
        # rotation is cylindrically degenerate (no preferred dihedral reference).
        anc_heavy_ref = sum(1 for nb in adj[anchor]
                            if nb != pivot and syms[nb] != "H")
        if anc_heavy_ref < 1:
            continue
        score = rot_heavy * 100 + len(rotating)
        dofs.append({
            "anchor": int(anchor),
            "pivot": int(pivot),
            "rotating": [int(x) for x in rotating],
            "score": int(score),
        })

    dofs.sort(key=lambda d: (-d["score"], d["anchor"], d["pivot"]))
    if len(dofs) > max_dofs:
        dofs = dofs[:max_dofs]
    return dofs


# ---------------------------------------------------------------------------
# Rigid rotation (Rodrigues) of a sub-tree about a bond axis
# ---------------------------------------------------------------------------


def _rotate_subtree(P: np.ndarray, origin: np.ndarray, axis: np.ndarray,
                    angle: float, idx: Sequence[int]) -> np.ndarray:
    """Return a copy of P with atoms ``idx`` rotated by ``angle`` (rad) about the
    line through ``origin`` along unit ``axis``.  Vectorised Rodrigues."""
    out = P.copy()
    nrm = float(np.linalg.norm(axis))
    if nrm < 1e-12 or not idx:
        return out
    u = axis / nrm
    c = math.cos(angle)
    s = math.sin(angle)
    sub = out[idx] - origin                       # (k,3)
    dot = sub @ u                                  # (k,)
    cross = np.cross(np.broadcast_to(u, sub.shape), sub)
    rot = sub * c + cross * s + np.outer(dot * (1.0 - c), u)
    out[idx] = rot + origin
    return out


# ---------------------------------------------------------------------------
# Clash objective L (pure geometry, H-inclusive)
# ---------------------------------------------------------------------------


def _excl_1_2_3(adj: Sequence[Sequence[int]], n: int) -> List[Set[int]]:
    """1-2 + 1-3 exclusion sets from the bond graph (same closure the detector +
    refiner use): excl[i] = adj[i] ∪ (⋃_{j∈adj[i]} adj[j]) − {i}."""
    excl: List[Set[int]] = []
    for i in range(n):
        s: Set[int] = set(adj[i])
        for j in adj[i]:
            s.update(adj[j])
        s.discard(i)
        excl.append(s)
    return excl


def _nonbonded_mask(syms: Sequence[str],
                    excl: Sequence[Set[int]]) -> np.ndarray:
    """Boolean (n,n) upper-triangular mask of pairs that count toward the clash
    objective: i<j and j not in excl[i] (not 1-2/1-3).  Built once per relax."""
    n = len(syms)
    m = np.zeros((n, n), dtype=bool)
    iu, ju = np.triu_indices(n, k=1)
    keep = np.fromiter((ju[k] not in excl[iu[k]] for k in range(len(iu))),
                       dtype=bool, count=len(iu))
    m[iu[keep], ju[keep]] = True
    return m


def loss_and_min(syms: Sequence[str], P: np.ndarray, mask: np.ndarray,
                 rsum: np.ndarray) -> Tuple[float, float]:
    """Return (L, d_min) where ``L`` is the clash objective and ``d_min`` is the
    minimum non-bonded (counted) pair distance.

    ``L = sum_{counted pairs} max(0, f*(vdw_i+vdw_j) - d)^2`` — H-inclusive, inter-
    AND intra-ligand, factor f=0.75.  ``mask`` (from :func:`_nonbonded_mask`) and
    ``rsum`` (= f*(vdw_i+vdw_j), precomputed) make this a pure vectorised pass."""
    diff = P[:, None, :] - P[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    over = np.where(mask, rsum - dist, 0.0)
    over = np.where(over > 0.0, over, 0.0)
    loss = float((over * over).sum())
    md = dist[mask]
    dmin = float(md.min()) if md.size else float("inf")
    return loss, dmin


def clash_loss(syms: Sequence[str], P: np.ndarray,
               excl: Sequence[Set[int]]) -> float:
    """L = sum over non-bonded (not 1-2/1-3) pairs of max(0, f*(vdw_i+vdw_j)-d)^2.
    H-inclusive, inter- AND intra-ligand, f=0.75.  Convenience wrapper used by
    tests/diagnostics; the hot loop uses :func:`loss_and_min` with cached arrays."""
    n = len(syms)
    if n < 2:
        return 0.0
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[:, None] + vdw[None, :])
    mask = _nonbonded_mask(syms, excl)
    return loss_and_min(syms, P, mask, rsum)[0]


# ---------------------------------------------------------------------------
# M-D invariant guard
# ---------------------------------------------------------------------------


def _md_pairs(syms: Sequence[str], P: np.ndarray) -> List[Tuple[int, int, float]]:
    """(metal, donor, distance) for every first-shell M-heavy coordination bond."""
    pairs: List[Tuple[int, int, float]] = []
    n = len(syms)
    for i in range(n):
        if not _bd._is_metal(syms[i]):
            continue
        for j in range(n):
            if i == j or syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _bd._ideal_bond(syms[i], syms[j]):
                pairs.append((i, j, d))
    return pairs


def _md_ok(base_pairs: Sequence[Tuple[int, int, float]], P: np.ndarray,
           tol: float) -> bool:
    for i, j, d0 in base_pairs:
        if abs(float(np.linalg.norm(P[i] - P[j])) - d0) > tol:
            return False
    return True


# ---------------------------------------------------------------------------
# Core: coordinate-descent torsion relaxation
# ---------------------------------------------------------------------------


def relax(syms: Sequence[str], P, frozen: Iterable[int],
          grid: int = _DEF_GRID, passes: int = _DEF_PASSES,
          md_tol: float = _DEF_MD_TOL,
          max_dofs: int = _DEF_MAX_DOFS,
          bond_pairs: Optional[Sequence[Tuple[int, int]]] = None) -> np.ndarray:
    """Joint torsion-space clash relaxation of an assembled complex.

    ``frozen``: indices that must not move (metal + donors).  ``bond_pairs``
    (optional): the true connectivity threaded from the builder (see
    :func:`_adjacency`); strongly preferred on crowded complexes.  Returns the
    relaxed coordinate array (float, same shape).  Accept-if-better, never-worse;
    the M-D invariant is enforced (>tol -> the offending move is rejected, and a
    final guard rolls the whole result back to the input on any violation).
    Deterministic; never raises (falls back to the input on any failure)."""
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return np.array(P, dtype=float)
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return P0
    frozen_set = set(int(x) for x in frozen)
    dofs = identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                         bond_pairs=bond_pairs)
    if not dofs:
        return P0

    adj, _bonds = _adjacency(syms, P0, bond_pairs)
    excl = _excl_1_2_3(adj, n)
    base_md = _md_pairs(syms, P0)
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[:, None] + vdw[None, :])
    mask = _nonbonded_mask(syms, excl)

    Pcur = P0.copy()
    best_loss, base_dmin = loss_and_min(syms, Pcur, mask, rsum)
    if best_loss <= 1e-12:
        return Pcur
    # never-worse-on-min floor: a move is rejected if it pushes the worst counted
    # contact below the input frame's worst contact.  The loss is a SUM, so without
    # this floor the optimiser could relieve one sub-threshold clash by squeezing an
    # already-comfortable pair (e.g. AVECUA: relieve a 1.7 A pair while pulling a
    # 2.48 A pair down to 1.8) -> a lower total L but a WORSE best-of-ensemble MIN.
    # Allow a hair of slack so finite-grid round-off cannot block a real improvement.
    dmin_floor = base_dmin - 1e-6

    angles = [2.0 * math.pi * k / float(grid) for k in range(grid)]
    for _ in range(max(1, passes)):
        improved = False
        for dof in dofs:
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
                    continue                       # identity = current state
                trial = _rotate_subtree(Pcur, origin, axis, ang, rot)
                if not np.all(np.isfinite(trial)):
                    continue
                if not _md_ok(base_md, trial, md_tol):
                    continue
                tl, tmin = loss_and_min(syms, trial, mask, rsum)
                if tmin < dmin_floor:
                    continue                       # would worsen the worst contact
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

    # final hard guard: if anything broke the M-D invariant or went non-finite, or
    # the worst contact ended up below the input's worst contact, roll the whole
    # relaxation back to the input frame (never-worse contract).
    if not np.all(np.isfinite(Pcur)) or not _md_ok(base_md, Pcur, md_tol):
        return P0
    fin_loss, fin_dmin = loss_and_min(syms, Pcur, mask, rsum)
    if fin_loss > best_loss + 1e-9 or fin_dmin < dmin_floor:
        return P0
    return Pcur


def bonds_from_blocks(metal_idx: int,
                      blocks: Sequence[Tuple[int, Sequence[Tuple[int, int]], int]]
                      ) -> List[Tuple[int, int]]:
    """Build a global bond list from per-ligand blocks.

    ``blocks``: one ``(offset, local_bonds, donor_local)`` per ligand, where
    ``offset`` is the ligand's first global atom index, ``local_bonds`` are its
    intra-ligand (i,j) pairs in LOCAL indices (atom order preserved through
    assembly), and ``donor_local`` is the donor's local index (a M-donor bond is
    emitted to ``metal_idx``; <0 to skip).  This is the TRUE connectivity that
    makes the relaxer robust against crowded-complex distance mis-perception."""
    out: List[Tuple[int, int]] = []
    for offset, local_bonds, donor_local in blocks:
        for (li, lj) in local_bonds:
            out.append((offset + int(li), offset + int(lj)))
        if donor_local is not None and int(donor_local) >= 0:
            out.append((metal_idx, offset + int(donor_local)))
    return out


def relax_if_enabled(syms: Sequence[str], P, frozen: Iterable[int],
                     bond_pairs: Optional[Sequence[Tuple[int, int]]] = None):
    """Wire-in entry: when ``DELFIN_FFFREE_TORSION_RELAX=1`` run :func:`relax`,
    else return ``P`` unchanged (byte-identical default-OFF).  ``bond_pairs``
    (optional) is the true connectivity from the builder.  Never raises."""
    if not _enabled():
        return P
    try:
        grid = _env_int("DELFIN_FFFREE_TORSION_GRID", _DEF_GRID, 4, 72)
        passes = _env_int("DELFIN_FFFREE_TORSION_PASSES", _DEF_PASSES, 1, 32)
        md_tol = _env_float("DELFIN_FFFREE_TORSION_MD_TOL", _DEF_MD_TOL, 0.0, 2.0)
        max_dofs = _env_int("DELFIN_FFFREE_TORSION_MAX_DOFS", _DEF_MAX_DOFS, 1, 256)
        return relax(syms, P, frozen, grid=grid, passes=passes,
                     md_tol=md_tol, max_dofs=max_dofs, bond_pairs=bond_pairs)
    except Exception:
        return P
