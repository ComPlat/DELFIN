"""delfin.fffree.conformer_dedup — Kabsch-RMSD pairwise matrix + Butina
clustering + severity-best per cluster.  The post-GRIP "essential set"
filter that removes the ~67% near-duplicate redundancy documented in
``project_grip_rmsd_dedup_essential_conformers_2026_06_02``.

Algorithm (Kabsch-RMSD + Butina clustering):

  1. Compute pairwise Kabsch-RMSD matrix (O(N^2 * M) where N = #frames,
     M = atoms).  Both structures centered, then optimal rotation via
     SVD (no reflection).
  2. Butina clustering:
       a) Sort frames by NUMBER OF NEIGHBOURS (within threshold) DESC.
       b) Greedy: take the highest-neighbour frame as cluster centroid,
          assign all its neighbours, mark assigned.
       c) Repeat until every frame has a cluster.
     Deterministic tiebreak: smaller frame index wins (input order is
     itself deterministic via Pólya x CP x rotamer ordering).
  3. Per cluster, keep the member with MINIMUM severity.  If severity
     is missing / inf, fall back to lowest frame index in the cluster.
  4. Return the kept frames sorted by severity ascending.

Universal: works on any (N, 3) coordinate set with a matching symbol
list.  Determinism: identical inputs -> identical output ordering.

Env-gates (all default OFF):
  - ``DELFIN_FFFREE_RMSD_DEDUP=1``           -> enable in emit pipeline
  - ``DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD``    -> Å (default 0.15, was 0.30
                                                  before subagent #129 follow-up;
                                                  tighter = more clusters, less
                                                  collapse of structurally
                                                  distinct conformers).
  - ``DELFIN_FFFREE_RMSD_DEDUP_MAX_KEEP``     -> cap (default 50)
  - ``DELFIN_FFFREE_SYMMETRIC_REP``           -> 0/1, default 0.  When 1,
                                                  the cluster representative
                                                  is chosen by a combined
                                                  point-group-order +
                                                  severity score (see
                                                  :func:`select_symmetric_representative`).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Set, Tuple

import numpy as np


# ------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------

DEFAULT_RMSD_THRESHOLD = 0.15   # Å — empirical near-duplicate cutoff
                                # (subagent #129 follow-up tightened from 0.30
                                # so structurally distinct rotamer pairs no
                                # longer collapse into the same cluster).
DEFAULT_MAX_KEEP = 50           # cap on returned frames per SMILES
DEFAULT_USE_HEAVY_ONLY = True   # ignore H atoms in RMSD (less noise)


# ------------------------------------------------------------------
# Env-tunable helpers
# ------------------------------------------------------------------


def _env_threshold(default: float = DEFAULT_RMSD_THRESHOLD) -> float:
    try:
        return float(os.environ.get(
            "DELFIN_FFFREE_RMSD_DEDUP_THRESHOLD",
            str(default),
        ))
    except (TypeError, ValueError):
        return default


def _env_max_keep(default: int = DEFAULT_MAX_KEEP) -> int:
    try:
        return max(1, int(os.environ.get(
            "DELFIN_FFFREE_RMSD_DEDUP_MAX_KEEP",
            str(default),
        )))
    except (TypeError, ValueError):
        return default


def _env_use_heavy_only(default: bool = DEFAULT_USE_HEAVY_ONLY) -> bool:
    v = os.environ.get("DELFIN_FFFREE_RMSD_DEDUP_HEAVY_ONLY", "1" if default else "0")
    return str(v).strip() not in ("0", "false", "False", "off")


def _env_pre_cluster_emit(default: bool = False) -> bool:
    """Pre-cluster-emit toggle (subagent #129 follow-up insight #3).

    When True (env ``DELFIN_FFFREE_PRE_CLUSTER_EMIT=1``), the caller is
    expected to cluster PRE-polished candidates first, then keep originals
    from busy clusters.  The flag is read by callers (converter_backend.py)
    and does not change :func:`dedup_by_rmsd` behaviour directly.
    """
    v = os.environ.get("DELFIN_FFFREE_PRE_CLUSTER_EMIT", "1" if default else "0")
    return str(v).strip() in ("1", "true", "True", "on")


def dedup_by_rmsd_preserve_originals(
    frames: Sequence,
    n_originals: int,
    threshold: Optional[float] = None,
    max_keep: Optional[int] = None,
    syms: Optional[Sequence[str]] = None,
    symmetric_rep: Optional[bool] = None,
) -> List:
    """RMSD-cluster ``frames``, then keep ALL originals (the first
    ``n_originals`` entries) plus the per-cluster best of every cluster
    that contains at least one non-original member.

    Subagent #129 follow-up: this implements the "RMSD-dedup BEFORE
    additive emission -- cluster first, then keep originals from busy
    clusters" insight.  Compared to :func:`dedup_by_rmsd`, originals are
    NEVER dropped even when an additive rotamer/pucker frame has lower
    severity in the same cluster -- they are still in the output.

    Returns a list of length ``<= max_keep``.  Originals preserve their
    input order; additives come after, sorted by severity.
    """
    th = float(threshold) if threshold is not None else _env_threshold()
    cap = int(max_keep) if max_keep is not None else _env_max_keep()
    if not frames:
        return []
    # Normalise input.
    norm: List[Tuple[str, np.ndarray, float]] = []
    is_frame = isinstance(frames[0], Frame) if len(frames) else False
    for f in frames:
        if isinstance(f, Frame):
            norm.append((str(f.label), np.asarray(f.P, dtype=float),
                          float(f.severity)))
        elif isinstance(f, tuple) and len(f) == 3:
            norm.append((str(f[0]), np.asarray(f[1], dtype=float),
                          float(f[2])))
        elif isinstance(f, tuple) and len(f) == 2:
            norm.append((str(f[0]), np.asarray(f[1], dtype=float), 0.0))
        elif isinstance(f, tuple) and len(f) >= 3:
            norm.append((str(f[0]), np.asarray(f[1], dtype=float),
                          float(f[2])))
        else:
            continue
    n = len(norm)
    if n == 0:
        return []
    n_orig = max(0, min(int(n_originals), n))
    coords_list = [t[1] for t in norm]
    sev_list = [t[2] for t in norm]
    M = pairwise_rmsd_matrix(coords_list, syms=syms)
    clusters = butina_cluster(M, th)
    if symmetric_rep is None:
        symmetric_rep = _env_symmetric_rep()

    keepers: List[int] = []
    seen: Set = set()  # type: ignore[name-defined]
    # 1) All originals first (in input order) -- ensures the as-built
    #    geometry never disappears even when a rotamer wins by score.
    for i in range(n_orig):
        if i in seen:
            continue
        keepers.append(int(i))
        seen.add(int(i))
    # 2) For every cluster that has a non-original member NOT already
    #    represented by an original we kept, pick its representative.
    orig_set = set(range(n_orig))
    for c in clusters:
        cluster_origs = [int(idx) for idx in c if int(idx) in orig_set]
        if cluster_origs:
            continue  # this cluster already has an original kept
        if symmetric_rep and syms is not None:
            best = select_symmetric_representative(
                c, coords_list, sev_list, list(syms),
            )
        else:
            best = min(c, key=lambda idx: (float(norm[idx][2]), int(idx)))
        if int(best) not in seen:
            keepers.append(int(best))
            seen.add(int(best))
    keepers = keepers[:cap]

    out = []
    from typing import cast as _cast  # local import to keep top-level clean
    for idx in keepers:
        src = frames[idx]
        if is_frame and isinstance(src, Frame):
            out.append(src)
        elif isinstance(src, tuple) and len(src) >= 3:
            out.append(src)
        elif isinstance(src, tuple) and len(src) == 2:
            out.append(src)
        else:
            out.append((norm[idx][0], norm[idx][1], norm[idx][2]))
    return out


def _env_symmetric_rep(default: bool = False) -> bool:
    """Cluster-representative strategy switch (subagent #129 follow-up).

    When True (env ``DELFIN_FFFREE_SYMMETRIC_REP=1``), the per-cluster
    representative is chosen by the combined point-group order + severity
    score implemented in :func:`select_symmetric_representative`.  When
    False (default, byte-identical to HEAD), we keep the historical
    minimum-severity-then-input-order rule.
    """
    v = os.environ.get("DELFIN_FFFREE_SYMMETRIC_REP", "1" if default else "0")
    return str(v).strip() in ("1", "true", "True", "on")


# ------------------------------------------------------------------
# Kabsch RMSD
# ------------------------------------------------------------------


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Return the optimal-rotation RMSD between two (N, 3) coordinate
    arrays via the Kabsch SVD algorithm (proper rotation only -- no
    reflection).

    Both arrays MUST have the same row count; otherwise returns +inf.
    Empty arrays return 0.0.

    Determinism: pure numpy operations; identical inputs -> identical
    output up to floating-point precision.
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape:
        return float("inf")
    n = P.shape[0]
    if n == 0:
        return 0.0
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    try:
        U, _, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        # Fall back to no-rotation RMSD.
        diff = Pc - Qc
        return float(np.sqrt(np.sum(diff * diff) / n))
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    if d == 0.0:
        d = 1.0
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    rotated = Pc @ R.T
    diff = rotated - Qc
    return float(np.sqrt(np.sum(diff * diff) / n))


def pairwise_rmsd_matrix(
    coords_list: Sequence[np.ndarray],
    syms: Optional[Sequence[str]] = None,
    heavy_only: Optional[bool] = None,
) -> np.ndarray:
    """Return the (N, N) Kabsch-RMSD distance matrix for ``coords_list``.

    When ``heavy_only`` is True (default from env) and ``syms`` is
    provided, H atoms are dropped before computing RMSD -- a common
    convention for conformer comparison.

    The matrix is symmetric with zero diagonal.  Order matches the
    input list.  Robust against mismatched array shapes (returns +inf
    for the offending pair).
    """
    if heavy_only is None:
        heavy_only = _env_use_heavy_only()
    n = len(coords_list)
    M = np.zeros((n, n), dtype=float)
    if n == 0:
        return M
    # Pre-extract the heavy-atom slice once per frame.
    if heavy_only and syms is not None:
        heavy_mask = np.array([s != "H" for s in syms])
        sliced = [np.asarray(P, dtype=float)[heavy_mask] for P in coords_list]
    else:
        sliced = [np.asarray(P, dtype=float) for P in coords_list]
    for i in range(n):
        for j in range(i + 1, n):
            r = kabsch_rmsd(sliced[i], sliced[j])
            M[i, j] = r
            M[j, i] = r
    return M


# ------------------------------------------------------------------
# Butina clustering
# ------------------------------------------------------------------


def butina_cluster(
    distances: np.ndarray,
    threshold: float,
) -> List[List[int]]:
    """Butina clustering on a symmetric distance matrix.

    Steps:

    1. For each frame, count neighbours at distance < threshold.
    2. Sort frames by neighbour count DESC; tiebreak by lower frame index.
    3. Greedy: pop the highest-neighbour frame as a new cluster centroid,
       assign all its unassigned neighbours to that cluster, mark
       everything assigned, repeat.

    Returns ``[cluster_0_indices, cluster_1_indices, ...]`` with each
    cluster list sorted ascending by frame index.

    Deterministic: stable np.argsort + ascending tiebreak.
    """
    n = int(distances.shape[0])
    if n == 0:
        return []
    if n == 1:
        return [[0]]
    # Boolean neighbour matrix (exclude diagonal).
    neigh_mat = (distances < float(threshold)) & ~np.eye(n, dtype=bool)
    neigh_count = neigh_mat.sum(axis=1)
    # Build deterministic sort key: (-count, index) so largest count first,
    # smallest index breaks ties.
    order = sorted(range(n), key=lambda i: (-int(neigh_count[i]), int(i)))
    assigned = [False] * n
    clusters: List[List[int]] = []
    for centroid in order:
        if assigned[centroid]:
            continue
        cluster = [int(centroid)]
        assigned[centroid] = True
        # Walk the neighbour set in ascending order for determinism.
        for j in range(n):
            if j == centroid or assigned[j]:
                continue
            if neigh_mat[centroid, j]:
                cluster.append(int(j))
                assigned[j] = True
        cluster.sort()
        clusters.append(cluster)
    return clusters


# ------------------------------------------------------------------
# Frame container
# ------------------------------------------------------------------


@dataclass
class Frame:
    """Lightweight container for a frame to be deduped."""

    label: str
    P: np.ndarray
    severity: float = float("inf")
    extra: dict = field(default_factory=dict)


# ------------------------------------------------------------------
# Symmetry-aware cluster representative (subagent #129 follow-up)
# ------------------------------------------------------------------


def _try_spglib_order(coords: np.ndarray, elements: Sequence[str], tol: float = 0.1) -> Optional[int]:
    """Return the point-group order via spglib if available, else ``None``.

    Wraps ``spglib.get_symmetry_dataset`` on a single-molecule "crystal" with
    a large cell (sufficient vacuum so the algorithm doesn't see periodic
    images).  Returns ``None`` if spglib is missing or the call fails.
    """
    try:
        import spglib  # type: ignore
    except Exception:
        return None
    try:
        P = np.asarray(coords, dtype=float)
        if P.size == 0:
            return None
        # Build a large cubic cell around the molecule (20 Å padding).
        lo = P.min(axis=0) - 10.0
        hi = P.max(axis=0) + 10.0
        cell = np.diag(hi - lo)
        frac = (P - lo) / np.maximum(hi - lo, 1e-9)
        # Map symbols -> integer Z-like labels for spglib.
        sym_to_id: dict = {}
        ids: List[int] = []
        for s in elements:
            if s not in sym_to_id:
                sym_to_id[s] = len(sym_to_id) + 1
            ids.append(sym_to_id[s])
        ds = spglib.get_symmetry_dataset((cell, frac, ids), symprec=float(tol))
        if ds is None:
            return None
        # spglib returns rotations under PBC -> #rotations is the group order
        # (it correctly returns 1 for C1, 2 for Ci/Cs, etc. with vacuum cell).
        n_ops = len(ds.get("rotations", []))
        return max(1, int(n_ops))
    except Exception:
        return None


def _custom_point_group_order(coords: np.ndarray, elements: Sequence[str], tol: float = 0.1) -> int:
    """Lightweight point-group order detector (Cn / σ / i) without spglib.

    Algorithm:
      1. Center the molecule at the centroid (weighted by 1 per heavy atom).
      2. Test inversion ``i``: P -> -P; if every -P[i] maps to a same-symbol
         P[j] within ``tol``, the order doubles.
      3. Test mirror planes σ along the principal-inertia axes: reflect P
         through each plane (x=0, y=0, z=0) and check the same mapping.
      4. Test rotation Cn (n=2..6) around each principal-inertia axis:
         apply the rotation, check the mapping.
      5. Return the largest detected order (lower bound; misses S2n axes).

    This is a fallback for environments without spglib; it returns 1 for C1
    molecules and rises with detectable symmetry.  Deterministic.
    """
    P = np.asarray(coords, dtype=float)
    n = P.shape[0]
    if n == 0:
        return 1
    P = P - P.mean(axis=0)
    syms = [str(s) for s in elements]

    def _maps_to_same(P_t: np.ndarray) -> bool:
        used = [False] * n
        for i in range(n):
            best_j = -1
            for j in range(n):
                if used[j] or syms[j] != syms[i]:
                    continue
                if float(np.linalg.norm(P_t[i] - P[j])) < float(tol):
                    best_j = j
                    break
            if best_j < 0:
                return False
            used[best_j] = True
        return True

    order = 1

    # 1) Inversion.
    if _maps_to_same(-P):
        order = max(order, 2)

    # 2) Mirror planes along Cartesian axes.
    for axis in range(3):
        P_t = P.copy()
        P_t[:, axis] = -P_t[:, axis]
        if _maps_to_same(P_t):
            order = max(order, 2)

    # 3) Cn rotations around Cartesian axes.  We try the canonical axes only
    #    (fast + universal).  Principal-inertia alignment would catch more,
    #    but the canonical sweep covers most relevant cases when the caller
    #    has built ideal-polyhedra-oriented complexes.
    for axis_idx in range(3):
        axis_vec = np.zeros(3)
        axis_vec[axis_idx] = 1.0
        for n_fold in (2, 3, 4, 5, 6):
            theta = 2.0 * np.pi / n_fold
            c, s = np.cos(theta), np.sin(theta)
            C = 1.0 - c
            x, y, z = float(axis_vec[0]), float(axis_vec[1]), float(axis_vec[2])
            R = np.array([
                [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
                [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
                [z * x * C - y * s, z * y * C + x * s, c + z * z * C   ],
            ], dtype=float)
            if _maps_to_same(P @ R.T):
                order = max(order, n_fold)
    return int(order)


def point_group_order(coords: np.ndarray, elements: Sequence[str],
                      tol: float = 0.1) -> int:
    """Return the point-group order of ``(coords, elements)``.

    Tries :func:`_try_spglib_order` first (when spglib is installed) and
    falls back to :func:`_custom_point_group_order` (canonical-axis Cn /
    σ / i detection).  Always returns at least 1.

    Universal: depends only on coordinates + element labels.  Deterministic:
    no randomness.

    Parameters
    ----------
    coords : (N, 3) ndarray
        Cartesian coordinates.
    elements : sequence of str
        Element symbols aligned with ``coords``.
    tol : float, default 0.1
        Position-deviation tolerance for equivalence.
    """
    if coords is None or len(coords) == 0:
        return 1
    o = _try_spglib_order(coords, elements, tol=tol)
    if isinstance(o, int) and o >= 1:
        return o
    return _custom_point_group_order(coords, elements, tol=tol)


def select_symmetric_representative(
    cluster_indices: Sequence[int],
    coords_list: Sequence[np.ndarray],
    severities: Sequence[float],
    symbols: Sequence[str],
    tol: float = 0.1,
) -> int:
    """Pick the cluster index that maximises ``pg_order + 1/(1+severity)``.

    Implements the user-anchored "symmetry priority" rule from
    ``project_grip_symmetry_priority_rotamers_2026_06_02``:

      score = pg_order * 1.0 + 1.0 / (1.0 + severity)

    Higher point-group symmetry wins; ties broken by lower severity, then
    by lower input index (deterministic).

    The function takes the FULL ``coords_list`` / ``severities`` arrays and
    a cluster-index list so the caller can re-use them across clusters.
    """
    best_idx = int(cluster_indices[0])
    best_score = -float("inf")
    for idx in cluster_indices:
        i = int(idx)
        try:
            pg = point_group_order(coords_list[i], symbols, tol=tol)
        except Exception:
            pg = 1
        sev = float(severities[i])
        score = float(pg) * 1.0 + 1.0 / (1.0 + max(sev, 0.0))
        # Tiebreak: prefer lower severity, then lower input index.
        if (score > best_score
            or (abs(score - best_score) < 1e-12
                and (sev, i) < (float(severities[best_idx]), int(best_idx)))):
            best_score = score
            best_idx = i
    return int(best_idx)


# ------------------------------------------------------------------
# Public API: dedup_by_rmsd
# ------------------------------------------------------------------


def dedup_by_rmsd(
    frames: Sequence,
    threshold: Optional[float] = None,
    max_keep: Optional[int] = None,
    syms: Optional[Sequence[str]] = None,
    symmetric_rep: Optional[bool] = None,
) -> List:
    """Cluster ``frames`` by Kabsch-RMSD, keep the minimum-severity member
    per cluster, and return the survivors sorted by severity ascending.

    ``frames`` can be:

    - a sequence of :class:`Frame` instances
    - a sequence of 3-tuples ``(label, P, severity)``
    - a sequence of 2-tuples ``(label, P)`` -- severity assumed 0.0

    Parameters
    ----------
    threshold : float, optional
        Å cutoff for cluster membership; default from env (0.3).
    max_keep : int, optional
        Cap on returned frames; default from env (50).
    syms : sequence of element symbols, optional
        When provided and heavy-only mode is on (default), H atoms are
        dropped from the RMSD computation.  All frames are assumed to
        share the same atom ordering.

    Returns
    -------
    list
        Same element type as the input.  Length <= ``min(len(frames),
        max_keep)``.  Sorted by ``severity`` ascending; ties broken by
        the original input order.

    Determinism: pairwise RMSD is deterministic; Butina order is
    deterministic; severity tiebreak by input order is deterministic.
    """
    th = float(threshold) if threshold is not None else _env_threshold()
    cap = int(max_keep) if max_keep is not None else _env_max_keep()
    if not frames:
        return []
    # Normalize the input into (label, P, severity) triples + remember the
    # original element type so we can rebuild the output in the same shape.
    norm: List[Tuple[str, np.ndarray, float]] = []
    is_frame = isinstance(frames[0], Frame) if len(frames) else False
    for f in frames:
        if isinstance(f, Frame):
            norm.append((str(f.label), np.asarray(f.P, dtype=float),
                          float(f.severity)))
        elif isinstance(f, tuple) and len(f) == 3:
            norm.append((str(f[0]), np.asarray(f[1], dtype=float),
                          float(f[2])))
        elif isinstance(f, tuple) and len(f) == 2:
            norm.append((str(f[0]), np.asarray(f[1], dtype=float), 0.0))
        elif isinstance(f, tuple) and len(f) >= 3:
            # Allow longer tuples (extra fields stored only for round-trip).
            norm.append((str(f[0]), np.asarray(f[1], dtype=float),
                          float(f[2])))
        else:
            # Unknown shape -> skip defensively (preserves "never crash").
            continue
    n = len(norm)
    if n == 0:
        return []
    coords_list = [t[1] for t in norm]
    sev_list = [t[2] for t in norm]
    M = pairwise_rmsd_matrix(coords_list, syms=syms)
    clusters = butina_cluster(M, th)
    # Resolve the representative-selection strategy.  Per-call kwarg wins,
    # then env-flag.  Default OFF -> byte-identical to HEAD.
    if symmetric_rep is None:
        symmetric_rep = _env_symmetric_rep()
    keepers: List[int] = []
    if symmetric_rep and syms is not None:
        for c in clusters:
            best = select_symmetric_representative(
                c, coords_list, sev_list, list(syms),
            )
            keepers.append(int(best))
    else:
        for c in clusters:
            best = min(c, key=lambda idx: (float(norm[idx][2]), int(idx)))
            keepers.append(int(best))
    # Sort kept frames by severity (ties: original input order).
    keepers.sort(key=lambda idx: (float(norm[idx][2]), int(idx)))
    keepers = keepers[:cap]
    # Rebuild output in the original input shape.
    out = []
    for idx in keepers:
        src = frames[idx]
        if is_frame and isinstance(src, Frame):
            out.append(src)
        elif isinstance(src, tuple) and len(src) >= 3:
            out.append(src)
        elif isinstance(src, tuple) and len(src) == 2:
            out.append(src)
        else:
            # Fallback: build a triple
            out.append((norm[idx][0], norm[idx][1], norm[idx][2]))
    return out


# ------------------------------------------------------------------
# Standalone smoke check
# ------------------------------------------------------------------


if __name__ == "__main__":
    # Two near-duplicate frames + one distinct frame -> dedup should
    # yield exactly 2 clusters.
    P0 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=float)
    P1 = P0 + np.array([0.05, 0.0, 0.0])         # near-duplicate
    P2 = np.array([[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 3, 0]], dtype=float)
    frames = [("A", P0, 1.0), ("B", P1, 2.0), ("C", P2, 1.5)]
    kept = dedup_by_rmsd(frames, threshold=0.3)
    print(f"input: {len(frames)} frames -> kept: {len(kept)} clusters")
    for lab, P, sev in kept:
        print(f"  label={lab} sev={sev}")
