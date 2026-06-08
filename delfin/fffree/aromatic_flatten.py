"""delfin.fffree.aromatic_flatten — universal aromatic-ring planarity gate.

The Mogul-primary DG embed honours every pairwise (lo, hi) bound it is given,
but the embedder leaves a soft-minimum residual on weakly-constrained pairs.
For aromatic ring systems (phenyl, pyridyl, NHC, ...) this residual shows up
as out-of-plane buckling: phenyl rings come out with max-deviation of
0.15-0.30 Å against the best-fit plane, where the CCDC distribution mean is
0.02-0.05 Å.

This module provides a UNIVERSAL post-embed projection that:

  1. Walks every aromatic atom set the RDKit SSSR + ``GetIsAromatic()``
     perception reports.
  2. Fits the best plane through the heavy ring atoms via SVD.
  3. If the max-deviation exceeds a threshold (default 0.10 Å), projects each
     ring atom (and any rigidly-attached H) onto that plane.

Universality contract:

  * No SMILES patterns, no per-class branches.
  * Only RDKit graph perception (``mol.GetAtomWithIdx(i).GetIsAromatic()``
    plus the SSSR ring set) is used to define the ring set.  Both are
    Hückel-graph-only and apply equally to phenyl, pyridyl, pyrazolyl,
    imidazolyl, NHC carbene rings, furan, thiophene, ...
  * The metal centre and any atom that sits inside the M-shell are FROZEN
    (zero displacement) so the M-D invariant is preserved.

Determinism contract:

  * Deterministic ring iteration order: sorted by (min atom idx, ...).
  * Pure SVD + linear algebra; no RNG, no I/O.
  * Idempotent: a second call on already-planar rings is a no-op.

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations

import os
from typing import List, Optional, Sequence, Set, Tuple

import numpy as np


# Default planarity threshold: max perpendicular deviation per ring atom
# beyond which we project the ring.  0.10 Å sits comfortably above the
# CCDC distribution mean (0.02-0.05 Å) but well below the soft-DG residual
# we observe on phenyl rings without this gate (0.15-0.30 Å).
_DEFAULT_MAX_DEV: float = 0.10

# Below this absolute per-atom out-of-plane offset, we never bother
# moving an atom — that is at the noise floor of the embedder.
_NOOP_FLOOR: float = 0.005


def _env_on(flag: str, default: str = "0") -> bool:
    return os.environ.get(flag, default).strip() == "1"


def aromatic_planarity_gate_enabled() -> bool:
    """True iff the post-embed aromatic-planarity gate is active.

    Default OFF.  When ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` the wrapper in
    ``assemble_via_mogul`` flips the default to ON for the Mogul-primary
    construction path (set ``DELFIN_FFFREE_AROMATIC_PLANARITY_GATE=0``
    explicitly to compare bytes against the pre-gate behaviour).
    """
    return _env_on("DELFIN_FFFREE_AROMATIC_PLANARITY_GATE", "0")


def fused_aromatic_plane_enabled() -> bool:
    """True iff the joint-plane projection for fused aromatic systems is active.

    Per-ring planarity (``flatten_aromatic_rings``) flattens each ring in
    isolation, but adjacent fused rings can still pivot against each other
    around their shared edge — phenanthroline (3 fused 6-rings, 14 atoms)
    comes out with plane-RMS 0.62 Å in the wild even when each individual
    ring lands inside the per-ring tolerance.  The fused-system pass groups
    rings that share ≥ 2 atoms and projects ALL atoms of the joint system
    onto ONE common SVD best-fit plane.

    Default OFF (byte-identity).  Auto-flipped to ON in
    ``assemble_via_mogul`` when ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` and the
    operator did not explicitly set this flag.  Universal: pure RDKit
    aromaticity + ring-overlap perception, no SMILES patterns.
    """
    return _env_on("DELFIN_FFFREE_FUSED_AROMATIC_PLANE", "0")


def _max_dev_from_plane(P_ring: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """Return ``(max_dev, centroid, normal)`` for a ring point cloud.

    ``max_dev`` is the largest absolute perpendicular deviation of any ring
    atom from the SVD best-fit plane.
    """
    centroid = P_ring.mean(axis=0)
    M = P_ring - centroid
    try:
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        return float("inf"), centroid, np.array([0.0, 0.0, 1.0])
    normal = Vt[-1]
    nn = float(np.linalg.norm(normal))
    if nn < 1e-9:
        return float("inf"), centroid, np.array([0.0, 0.0, 1.0])
    normal = normal / nn
    offsets = M @ normal
    return float(np.max(np.abs(offsets))), centroid, normal


def measure_ring_planarity(P: np.ndarray, ring_atom_idxs: Sequence[int]) -> float:
    """Public helper: max perpendicular deviation of ring atoms from their plane.

    Returns +inf on degenerate input (< 3 points or collinear).
    """
    pts = np.asarray(P, dtype=float)[list(ring_atom_idxs)]
    if pts.shape[0] < 3:
        return float("inf")
    max_dev, _, _ = _max_dev_from_plane(pts)
    return max_dev


def detect_aromatic_rings(mol) -> List[Tuple[int, ...]]:
    """Return all aromatic rings in ``mol`` as sorted tuples of atom indices.

    A ring is "aromatic" iff EVERY atom in it has ``GetIsAromatic() == True``
    AND the ring has at least 3 atoms.  Universal: pure RDKit Hückel
    perception, no SMILES, no element ``if``-branches.

    Deterministic order: sorted lex by atom-idx tuple.
    """
    out: List[Tuple[int, ...]] = []
    try:
        rings = mol.GetRingInfo().AtomRings()
    except Exception:
        return out
    for ring in rings:
        ring_ints = sorted(int(a) for a in ring)
        if len(ring_ints) < 3:
            continue
        is_arom = True
        for a in ring_ints:
            try:
                if not mol.GetAtomWithIdx(a).GetIsAromatic():
                    is_arom = False
                    break
            except Exception:
                is_arom = False
                break
        if is_arom:
            out.append(tuple(ring_ints))
    out.sort()
    return out


def _ride_along_hs(
    mol,
    ring_atoms: Sequence[int],
) -> List[Tuple[int, int]]:
    """Return ``(H_idx, parent_idx)`` pairs for H atoms attached to ring atoms.

    Only direct H neighbours of a ring atom are considered (Cipso-H, Nipso-H).
    Hydrogens of substituent groups (e.g. methyl Hs on a tolyl ring) are NOT
    included — they ride with the substituent's sp3 carbon, not with the ring.
    """
    pairs: List[Tuple[int, int]] = []
    ring_set: Set[int] = set(int(a) for a in ring_atoms)
    for parent in sorted(ring_set):
        try:
            patom = mol.GetAtomWithIdx(int(parent))
        except Exception:
            continue
        for nb in patom.GetNeighbors():
            try:
                if nb.GetSymbol() == "H":
                    pairs.append((int(nb.GetIdx()), int(parent)))
            except Exception:
                continue
    pairs.sort()
    return pairs


def flatten_aromatic_rings(
    mol,
    P: np.ndarray,
    *,
    max_dev_threshold: float = _DEFAULT_MAX_DEV,
    frozen_idxs: Optional[Sequence[int]] = None,
) -> Tuple[np.ndarray, List[dict]]:
    """Project every aromatic ring in ``mol`` onto its SVD best-fit plane.

    Parameters
    ----------
    mol : RDKit Mol
        Source molecule with full bond + aromaticity perception applied.
    P : (N, 3) ndarray
        Cartesian coordinates.
    max_dev_threshold : float
        Skip rings whose max deviation is already at or below this value.
    frozen_idxs : sequence of int, optional
        Atom indices that MUST NOT move (typically the metal + every
        coordinating donor; this preserves the M-D invariant).  Frozen ring
        atoms anchor the plane in place (the plane is still fit through
        their positions); only NON-frozen ring atoms are translated.
        Hydrogens attached to frozen ring atoms also stay put.

    Returns
    -------
    P_out : (N, 3) ndarray
        Updated coordinates (copy; input ``P`` is not mutated).
    diagnostics : list of dict
        One entry per ring processed (whether flattened or skipped),
        containing ``{ring, max_dev_before, max_dev_after, flattened}``.

    Notes
    -----
    Universal: depends only on RDKit aromaticity perception + SVD.  No
    SMILES patterns, no per-class branches, no element-specific tables.
    """
    P_out = np.asarray(P, dtype=float).copy()
    frozen: Set[int] = set(int(i) for i in (frozen_idxs or ()))
    diagnostics: List[dict] = []

    rings = detect_aromatic_rings(mol)
    for ring in rings:
        pts = P_out[list(ring)]
        if pts.shape[0] < 3:
            continue
        max_dev_before, centroid, normal = _max_dev_from_plane(pts)
        if not np.isfinite(max_dev_before):
            diagnostics.append({
                "ring": tuple(ring),
                "max_dev_before": float("inf"),
                "max_dev_after": float("inf"),
                "flattened": False,
            })
            continue
        if max_dev_before <= max_dev_threshold:
            diagnostics.append({
                "ring": tuple(ring),
                "max_dev_before": max_dev_before,
                "max_dev_after": max_dev_before,
                "flattened": False,
            })
            continue

        # Identify frozen ring atoms (donor / metal atoms that sit IN this
        # ring) — these atoms cannot move, so the plane must be ANCHORED on
        # them.  Otherwise the projection translates non-frozen atoms onto
        # a plane that the frozen atom doesn't sit on, and the post max-dev
        # equals the frozen atom's offset (the bug user-eye SIYMEU exposes:
        # NHC carbene Cipso pinned 0.14 Å off the SVD plane).
        frozen_in_ring = [k for k, idx in enumerate(ring) if int(idx) in frozen]
        if frozen_in_ring:
            frozen_pts = pts[frozen_in_ring]
            anchor = frozen_pts.mean(axis=0)
            # Re-fit the plane centred on the frozen-atom centroid so the
            # plane passes through the anchor by construction.  The normal
            # is the smallest singular vector of (pts - anchor).
            M_anchor = pts - anchor
            try:
                _, _, Vt_anch = np.linalg.svd(M_anchor, full_matrices=False)
                n_anch = Vt_anch[-1]
                nn = float(np.linalg.norm(n_anch))
                if nn > 1e-9:
                    normal = n_anch / nn
                    centroid = anchor
            except np.linalg.LinAlgError:
                pass  # fall through with original plane

        # Per-atom offsets (signed) against the (possibly anchored) plane
        offsets = (pts - centroid) @ normal
        # Build move map: ring atoms get -offset * normal; frozen atoms get
        # zero displacement (plane is anchored on them, so their post-offset
        # is < 0.005 Å by construction for the anchored case; for the
        # un-frozen case the plane is the SVD best-fit so the centroid lies
        # on it).
        moves = {}
        for k, idx in enumerate(ring):
            if int(idx) in frozen:
                moves[int(idx)] = np.zeros(3, dtype=float)
                continue
            if abs(offsets[k]) < _NOOP_FLOOR:
                moves[int(idx)] = np.zeros(3, dtype=float)
                continue
            moves[int(idx)] = -float(offsets[k]) * normal
        # Apply ring-atom moves
        for idx, mv in moves.items():
            P_out[int(idx)] = P_out[int(idx)] + mv
        # Ride-along Hs (only Cipso-H / Nipso-H)
        for h_idx, parent in _ride_along_hs(mol, ring):
            mv = moves.get(int(parent))
            if mv is None:
                continue
            P_out[int(h_idx)] = P_out[int(h_idx)] + mv

        # Recompute max dev after projection
        pts_after = P_out[list(ring)]
        max_dev_after, _, _ = _max_dev_from_plane(pts_after)
        diagnostics.append({
            "ring": tuple(ring),
            "max_dev_before": max_dev_before,
            "max_dev_after": max_dev_after,
            "flattened": True,
        })

    return P_out, diagnostics


# ---------------------------------------------------------------------------
# Fused aromatic system handling (phenanthroline / naphthalene / indole, ...)
# ---------------------------------------------------------------------------
#
# Per-ring planarity (above) is NECESSARY but not SUFFICIENT for fused
# poly-aromatic systems.  Two fused 6-rings can each pass an isolated
# per-ring planarity check while pivoting against each other around their
# shared edge — phenanthroline (3 fused 6-rings, 14 atoms) is seen in the
# wild with plane-RMS 0.62 Å and max-deviation 1.04 Å.
#
# This fused-system pass groups rings that share ≥ 2 atoms (one shared
# edge) via Union-Find on the SSSR ring set, fits ONE common plane
# through the joint heavy-atom cloud, and projects every atom of the
# joint system onto that single plane.  It is the universal geometric
# expression of "extended π-system flatness": RDKit's Hückel perception
# decides which atoms are π-coupled, the SVD decides where the common
# π-plane sits.
#
# Universality contract (mirrors the per-ring pass):
#
#   * No SMILES patterns, no per-class branches.
#   * Detection: ring/ring adjacency via shared-atom-count ≥ 2 from the
#     SSSR ring set.  Equally captures naphthalene, anthracene, pyrene,
#     phenanthrene, phenanthroline, indole, quinoline, carbazole, ...
#   * Frozen atoms (metal + every coordinating donor) ANCHOR the plane:
#     when any frozen atom sits inside the fused system, the joint plane
#     is re-fit so it passes through the frozen-atom centroid; frozen
#     atoms then receive zero displacement.
#   * Determinism: rings sorted lex by min atom-idx, fused-system index
#     assigned by lex-order of the lowest ring inside it.


def group_fused_aromatic_rings(
    rings: Sequence[Sequence[int]],
) -> List[List[Tuple[int, ...]]]:
    """Group aromatic rings into fused-systems.

    Two rings are considered "fused" iff they share ≥ 2 atom indices
    (one shared edge).  Returns a list of fused-system groups; each
    group is a list of ring tuples.  Single-ring (isolated) groups are
    also returned, so the caller may decide to skip them (the per-ring
    pass handles those).

    Deterministic order: groups sorted by the lex-lowest ring inside;
    rings inside a group also sorted lex.
    """
    R = [tuple(sorted(int(a) for a in r)) for r in rings]
    n = len(R)
    if n == 0:
        return []
    parent = list(range(n))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    sets = [set(r) for r in R]
    for i in range(n):
        for j in range(i + 1, n):
            if len(sets[i] & sets[j]) >= 2:
                _union(i, j)

    buckets: dict = {}
    for i in range(n):
        root = _find(i)
        buckets.setdefault(root, []).append(R[i])

    groups: List[List[Tuple[int, ...]]] = []
    for root in sorted(buckets.keys()):
        grp = sorted(buckets[root])
        groups.append(grp)
    groups.sort(key=lambda g: g[0] if g else tuple())
    return groups


def detect_fused_aromatic_systems(mol) -> List[Tuple[int, ...]]:
    """Return joint atom-sets of fused aromatic systems with ≥ 2 rings.

    Each returned tuple is the SORTED union of atom indices belonging
    to all rings of a single fused-system that contains MORE THAN one
    aromatic ring.  Isolated (single-ring) aromatic rings are NOT
    returned here — they are handled by ``flatten_aromatic_rings``.

    Deterministic: union-tuples sorted by their first atom-idx.
    """
    rings = detect_aromatic_rings(mol)
    groups = group_fused_aromatic_rings(rings)
    out: List[Tuple[int, ...]] = []
    for grp in groups:
        if len(grp) < 2:
            continue
        union: Set[int] = set()
        for r in grp:
            union.update(int(a) for a in r)
        out.append(tuple(sorted(union)))
    out.sort()
    return out


def measure_fused_planarity(
    P: np.ndarray,
    atom_idxs: Sequence[int],
) -> float:
    """Public helper: max perpendicular deviation of fused-system atoms.

    Returns +inf for degenerate input (< 3 points or collinear).
    """
    pts = np.asarray(P, dtype=float)[list(atom_idxs)]
    if pts.shape[0] < 3:
        return float("inf")
    max_dev, _, _ = _max_dev_from_plane(pts)
    return max_dev


def flatten_fused_aromatic_systems(
    mol,
    P: np.ndarray,
    *,
    max_dev_threshold: float = _DEFAULT_MAX_DEV,
    frozen_idxs: Optional[Sequence[int]] = None,
) -> Tuple[np.ndarray, List[dict]]:
    """Project every fused aromatic system in ``mol`` onto its joint plane.

    A "fused system" is a connected set of aromatic rings where every
    adjacent pair shares ≥ 2 atoms (one shared edge).  All heavy atoms
    of the joint system are fit to ONE common SVD best-fit plane; each
    atom is then shifted perpendicular to that plane by its signed
    offset.  Frozen atoms (typically metal + coordinating donors) sit
    on the plane by ANCHORING: the plane is re-fit so it passes through
    the frozen-atom centroid before atoms are moved.

    Parameters
    ----------
    mol : RDKit Mol
        Source molecule with aromaticity perception applied.
    P : (N, 3) ndarray
        Cartesian coordinates.
    max_dev_threshold : float
        Skip fused-systems whose joint max-deviation is already ≤ this
        value.  Default 0.10 Å (same scale as the per-ring pass).
    frozen_idxs : sequence of int, optional
        Atom indices that MUST NOT move (typically metal + every donor).
        Frozen atoms inside the joint system anchor the plane (SVD is
        still done across all joint atoms, but the plane is re-centred
        on the frozen-atom centroid); frozen atoms receive zero
        displacement.  H atoms attached to frozen ring atoms also
        stay put.

    Returns
    -------
    P_out : (N, 3) ndarray
        Updated coordinates (copy).
    diagnostics : list of dict
        One entry per fused-system processed:
        ``{atoms, n_rings, max_dev_before, max_dev_after, flattened}``.

    Universal: depends only on RDKit aromaticity + ring overlap + SVD.
    No SMILES patterns, no element-specific tables.  Single-ring
    aromatics are SKIPPED here (the per-ring pass owns them).
    """
    P_out = np.asarray(P, dtype=float).copy()
    frozen: Set[int] = set(int(i) for i in (frozen_idxs or ()))
    diagnostics: List[dict] = []

    rings = detect_aromatic_rings(mol)
    groups = group_fused_aromatic_rings(rings)

    for grp in groups:
        if len(grp) < 2:
            # Single-ring system — leave to flatten_aromatic_rings.
            continue
        joint: List[int] = sorted({int(a) for r in grp for a in r})
        if len(joint) < 3:
            continue

        pts = P_out[joint]
        max_dev_before, centroid, normal = _max_dev_from_plane(pts)
        if not np.isfinite(max_dev_before):
            diagnostics.append({
                "atoms": tuple(joint),
                "n_rings": len(grp),
                "max_dev_before": float("inf"),
                "max_dev_after": float("inf"),
                "flattened": False,
            })
            continue
        if max_dev_before <= max_dev_threshold:
            diagnostics.append({
                "atoms": tuple(joint),
                "n_rings": len(grp),
                "max_dev_before": max_dev_before,
                "max_dev_after": max_dev_before,
                "flattened": False,
            })
            continue

        # Identify frozen atoms inside the joint set and anchor the plane
        # on their centroid (mirrors the per-ring pass — NHC carbene Cipso
        # / donor N atoms cannot move, so the plane must pass through
        # them by construction).
        frozen_in_joint = [k for k, idx in enumerate(joint) if int(idx) in frozen]
        if frozen_in_joint:
            frozen_pts = pts[frozen_in_joint]
            anchor = frozen_pts.mean(axis=0)
            M_anchor = pts - anchor
            try:
                _, _, Vt_anch = np.linalg.svd(M_anchor, full_matrices=False)
                n_anch = Vt_anch[-1]
                nn = float(np.linalg.norm(n_anch))
                if nn > 1e-9:
                    normal = n_anch / nn
                    centroid = anchor
            except np.linalg.LinAlgError:
                pass

        offsets = (pts - centroid) @ normal
        moves: dict = {}
        for k, idx in enumerate(joint):
            if int(idx) in frozen:
                moves[int(idx)] = np.zeros(3, dtype=float)
                continue
            if abs(offsets[k]) < _NOOP_FLOOR:
                moves[int(idx)] = np.zeros(3, dtype=float)
                continue
            moves[int(idx)] = -float(offsets[k]) * normal
        for idx, mv in moves.items():
            P_out[int(idx)] = P_out[int(idx)] + mv
        for h_idx, parent in _ride_along_hs(mol, joint):
            mv = moves.get(int(parent))
            if mv is None:
                continue
            P_out[int(h_idx)] = P_out[int(h_idx)] + mv

        pts_after = P_out[joint]
        max_dev_after, _, _ = _max_dev_from_plane(pts_after)
        diagnostics.append({
            "atoms": tuple(joint),
            "n_rings": len(grp),
            "max_dev_before": max_dev_before,
            "max_dev_after": max_dev_after,
            "flattened": True,
        })

    return P_out, diagnostics


__all__ = [
    "aromatic_planarity_gate_enabled",
    "detect_aromatic_rings",
    "detect_fused_aromatic_systems",
    "flatten_aromatic_rings",
    "flatten_fused_aromatic_systems",
    "fused_aromatic_plane_enabled",
    "group_fused_aromatic_rings",
    "measure_fused_planarity",
    "measure_ring_planarity",
]
