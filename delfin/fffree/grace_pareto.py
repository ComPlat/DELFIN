"""delfin.fffree.grace_pareto — Pareto-frontier selection + adaptive
ensemble size for GRACE.

Single-score severity ranking biases the GRACE survivor set toward one
metric (the linear combination ``severity + 10·clash + cshm``) and may
discard Pareto-optimal conformers that are slightly worse on one axis
but strictly better on another.  This module replaces the top-K
severity sort with a deterministic Pareto-frontier selection over four
independent axes:

    1. ``mogul_sev``       — Mogul/CCDC Mahalanobis severity (minimise)
    2. ``interlig_clash``  — inter-ligand clash count (minimise)
    3. ``cshm``            — donor-polyhedron continuous shape measure
                             (minimise; ``inf`` → 1e6 sentinel)
    4. ``arom_plan``       — aromatic-ring out-of-plane maximum (minimise;
                             0 = perfectly planar)

The Pareto-frontier set is *all* candidates that no other candidate
strictly dominates on every axis.  Dominated candidates are ranked by
the legacy single-score (severity-weighted) and appended to fill the
ensemble up to ``max_keep``.

Determinism contract
====================

* All tie-breaks use the deterministic key
  ``(score, isomer_id, ring_id, rotamer_id, label)`` matching the
  ``_rmsd_dedup`` legacy invariant.
* The frontier is computed via an O(n²) sweep with strict ``<``
  comparisons; outputs are lex-sorted.
* No RNG, no time-dependent fall-backs.

Default-OFF byte-identity
=========================

When ``DELFIN_FFFREE_GRACE_PARETO=1`` is NOT set, :func:`pareto_active`
returns ``False`` and the consumer (``grace_ensemble._rmsd_dedup``)
keeps the legacy severity-sort code path bit-for-bit.

When ``DELFIN_FFFREE_GRACE_ADAPTIVE=1`` is NOT set,
:func:`adaptive_active` returns ``False`` and the consumer iterates
the full deterministic loop without an early-stopping check.

The two flags are independent: Pareto selection can run with the legacy
fixed-N iteration budget, and adaptive iteration can run with the
legacy severity-sort.

Public API
==========

* :func:`pareto_active`                 — env gate (default OFF)
* :func:`adaptive_active`               — env gate (default OFF)
* :func:`resolve_max_ensemble`          — env-driven hard cap
* :func:`resolve_adaptive_patience`     — plateau iteration count
* :func:`compute_arom_plan_score`       — aromatic out-of-plane scalar
* :func:`compute_pareto_frontier`       — frontier + rank for dominated
* :func:`pareto_select`                 — frontier + top-K dominated
* :func:`adaptive_plateau_check`        — stop criterion for streaming
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "ENV_PARETO_ENABLE",
    "ENV_ADAPTIVE_ENABLE",
    "ENV_MAX_ENSEMBLE",
    "ENV_ADAPTIVE_PATIENCE",
    "DEFAULT_MAX_ENSEMBLE",
    "DEFAULT_ADAPTIVE_PATIENCE",
    "PARETO_AXES",
    "pareto_active",
    "adaptive_active",
    "resolve_max_ensemble",
    "resolve_adaptive_patience",
    "compute_arom_plan_score",
    "candidate_axes",
    "dominates",
    "compute_pareto_frontier",
    "pareto_select",
    "adaptive_plateau_check",
    "AdaptiveState",
]

_LOG = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Env-gate constants
# ---------------------------------------------------------------------------

ENV_PARETO_ENABLE = "DELFIN_FFFREE_GRACE_PARETO"
ENV_ADAPTIVE_ENABLE = "DELFIN_FFFREE_GRACE_ADAPTIVE"
ENV_MAX_ENSEMBLE = "DELFIN_FFFREE_GRACE_MAX_ENSEMBLE"
ENV_ADAPTIVE_PATIENCE = "DELFIN_FFFREE_GRACE_ADAPTIVE_PATIENCE"

DEFAULT_MAX_ENSEMBLE: int = 50
DEFAULT_ADAPTIVE_PATIENCE: int = 3

# Lex-stable axis order (deterministic frontier traversal).
PARETO_AXES: Tuple[str, ...] = (
    "mogul_sev",
    "interlig_clash",
    "cshm",
    "arom_plan",
)


# ---------------------------------------------------------------------------
# Env parsers
# ---------------------------------------------------------------------------
def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name, "").strip().lower()
    if not raw:
        return bool(default)
    if raw in ("1", "true", "yes", "on"):
        return True
    if raw in ("0", "false", "no", "off"):
        return False
    return bool(default)


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return int(default)
    try:
        return int(raw)
    except (TypeError, ValueError):
        return int(default)


def pareto_active() -> bool:
    """``True`` iff :data:`ENV_PARETO_ENABLE` is set.  Default OFF."""
    return _env_bool(ENV_PARETO_ENABLE, default=False)


def adaptive_active() -> bool:
    """``True`` iff :data:`ENV_ADAPTIVE_ENABLE` is set.  Default OFF."""
    return _env_bool(ENV_ADAPTIVE_ENABLE, default=False)


def resolve_max_ensemble(arg: Optional[int] = None) -> int:
    """Hard upper cap on the per-isomer ensemble size."""
    if arg is not None and int(arg) > 0:
        return int(arg)
    return max(1, _env_int(ENV_MAX_ENSEMBLE, DEFAULT_MAX_ENSEMBLE))


def resolve_adaptive_patience(arg: Optional[int] = None) -> int:
    """Number of consecutive plateau iterations before adaptive stop."""
    if arg is not None and int(arg) > 0:
        return int(arg)
    return max(1, _env_int(ENV_ADAPTIVE_PATIENCE, DEFAULT_ADAPTIVE_PATIENCE))


# ---------------------------------------------------------------------------
# Auxiliary axis: aromatic ring planarity
# ---------------------------------------------------------------------------
def compute_arom_plan_score(P: np.ndarray, mol) -> float:
    """Max out-of-plane deviation of aromatic ring atoms (Å).

    For every aromatic ring (as flagged by RDKit), fits a least-squares
    plane through the ring atoms and records the maximum perpendicular
    deviation of any ring atom to that plane.  Returns the maximum over
    all aromatic rings; 0.0 when no aromatic rings are present.  Pure /
    deterministic / no RNG.

    On any RDKit / numpy failure returns 0.0 (defensive: do not
    penalise a candidate for a metric the pipeline could not score).
    """
    if mol is None:
        return 0.0
    try:
        ri = mol.GetRingInfo()
        rings = ri.AtomRings()
    except Exception:
        return 0.0
    if not rings:
        return 0.0
    try:
        P_arr = np.asarray(P, dtype=np.float64)
    except Exception:
        return 0.0

    worst = 0.0
    n_atoms_total = P_arr.shape[0]
    for ring in rings:
        try:
            # Aromatic gate: at least one bond in the ring is aromatic.
            arom = False
            for i, a_idx in enumerate(ring):
                j_idx = ring[(i + 1) % len(ring)]
                b = mol.GetBondBetweenAtoms(int(a_idx), int(j_idx))
                if b is not None and b.GetIsAromatic():
                    arom = True
                    break
            if not arom:
                continue

            idx = [int(a) for a in ring if 0 <= int(a) < n_atoms_total]
            if len(idx) < 4:
                continue

            pts = P_arr[idx]
            c = pts.mean(axis=0)
            X = pts - c
            # Plane normal = smallest singular vector of X.
            _, _, Vt = np.linalg.svd(X, full_matrices=False)
            n = Vt[-1]
            n_norm = float(np.linalg.norm(n))
            if n_norm < 1e-12:
                continue
            n = n / n_norm
            dev = float(np.max(np.abs(X @ n)))
            if dev > worst:
                worst = dev
        except Exception:
            continue
    return float(worst)


# ---------------------------------------------------------------------------
# Axis vector extraction
# ---------------------------------------------------------------------------
def candidate_axes(
    candidate,
    mol=None,
    axes: Sequence[str] = PARETO_AXES,
) -> Tuple[float, ...]:
    """Pure axis-vector extraction from a GraceCandidate.

    ``inf`` cshm values are mapped to 1e6 so the Pareto comparison is
    finite-valued and well-ordered.  ``arom_plan`` is computed inline
    from the candidate's coords + the assembled ``mol``; when ``mol``
    is None the axis collapses to 0.0 (neutral).
    """
    vec: List[float] = []
    for axis in axes:
        if axis == "mogul_sev":
            v = float(getattr(candidate, "severity", float("inf")))
            if not np.isfinite(v):
                v = 1e9
            vec.append(v)
        elif axis == "interlig_clash":
            v = int(getattr(candidate, "clash_count", 0))
            vec.append(float(v))
        elif axis == "cshm":
            v = float(getattr(candidate, "cshm", float("inf")))
            if not np.isfinite(v):
                v = 1e6
            vec.append(v)
        elif axis == "arom_plan":
            try:
                v = float(compute_arom_plan_score(
                    np.asarray(getattr(candidate, "P"), dtype=np.float64),
                    mol,
                ))
            except Exception:
                v = 0.0
            vec.append(v)
        else:
            # Unknown axis: neutral.
            vec.append(0.0)
    return tuple(vec)


# ---------------------------------------------------------------------------
# Dominance
# ---------------------------------------------------------------------------
def dominates(a: Sequence[float], b: Sequence[float], tol: float = 1e-12) -> bool:
    """``True`` iff ``a`` Pareto-dominates ``b`` (all axes minimised).

    Convention: ``a`` dominates ``b`` when ``a[i] <= b[i] + tol`` for
    every axis AND ``a[i] < b[i] - tol`` for at least one axis.
    """
    if len(a) != len(b):
        raise ValueError("dominates: vector length mismatch")
    weak_better = True
    strict_any = False
    for ai, bi in zip(a, b):
        if ai > bi + tol:
            weak_better = False
            break
        if ai < bi - tol:
            strict_any = True
    return weak_better and strict_any


# ---------------------------------------------------------------------------
# Pareto frontier
# ---------------------------------------------------------------------------
@dataclass
class _CandSlot:
    cand: Any
    vec: Tuple[float, ...]
    idx: int
    score: float
    sort_key: Tuple[Any, ...]
    dominated: bool = False


def _sort_key(c: Any) -> Tuple[Any, ...]:
    return (
        float(getattr(c, "score", float("inf"))),
        int(getattr(c, "isomer_id", 0)),
        int(getattr(c, "ring_id", 0)),
        int(getattr(c, "rotamer_id", 0)),
        str(getattr(c, "label", "")),
    )


def compute_pareto_frontier(
    candidates: Sequence[Any],
    mol=None,
    axes: Sequence[str] = PARETO_AXES,
) -> Tuple[List[Any], List[Tuple[Any, int]]]:
    """Compute the Pareto frontier over the given axes.

    Returns ``(dominant, dominated_with_rank)`` where:

    * ``dominant`` — lex-sorted list of frontier candidates.
    * ``dominated_with_rank`` — list of ``(candidate, dominator_count)``
      pairs for non-frontier candidates, sorted by
      ``(dominator_count, _sort_key(c))`` so the user can promote the
      "least dominated" first.

    Determinism: input order is irrelevant; outputs are lex-sorted.
    """
    if not candidates:
        return [], []

    slots: List[_CandSlot] = []
    for i, c in enumerate(candidates):
        vec = candidate_axes(c, mol=mol, axes=axes)
        slots.append(_CandSlot(
            cand=c,
            vec=vec,
            idx=i,
            score=float(getattr(c, "score", float("inf"))),
            sort_key=_sort_key(c),
        ))

    # O(n²) dominance sweep with dominator counts.
    dom_count = [0] * len(slots)
    for i, si in enumerate(slots):
        for j, sj in enumerate(slots):
            if i == j:
                continue
            if dominates(sj.vec, si.vec):
                dom_count[i] += 1
        if dom_count[i] > 0:
            si.dominated = True

    dominant = [s for s in slots if not s.dominated]
    dominated = [(slots[i], dom_count[i]) for i in range(len(slots))
                 if slots[i].dominated]

    dominant.sort(key=lambda s: s.sort_key)
    dominated.sort(key=lambda t: (int(t[1]), t[0].sort_key))

    dominant_cands = [s.cand for s in dominant]
    dominated_cands = [(s.cand, int(c)) for (s, c) in dominated]
    return dominant_cands, dominated_cands


def pareto_select(
    candidates: Sequence[Any],
    max_keep: int,
    mol=None,
    axes: Sequence[str] = PARETO_AXES,
) -> List[Any]:
    """Pareto-frontier ensemble; fill to ``max_keep`` with dominated tail.

    Algorithm:

      1. ``dominant, dominated = compute_pareto_frontier(...)``
      2. If ``len(dominant) >= max_keep``: return dominant lex-sorted
         then truncated.
      3. Otherwise append the dominated candidates in
         dominator-count-then-sort-key order until ``max_keep`` is hit.

    The returned ensemble always contains the lex-smallest-key candidate
    of the frontier (matches legacy ``best`` semantics).  Empty input
    returns empty output.
    """
    if not candidates or max_keep <= 0:
        return []

    dominant, dominated = compute_pareto_frontier(
        candidates, mol=mol, axes=axes
    )

    if len(dominant) >= max_keep:
        return list(dominant[: int(max_keep)])

    out: List[Any] = list(dominant)
    for cand, _rank in dominated:
        if len(out) >= int(max_keep):
            break
        out.append(cand)
    return out


# ---------------------------------------------------------------------------
# Adaptive plateau detection
# ---------------------------------------------------------------------------
@dataclass
class AdaptiveState:
    """Streaming state for adaptive ensemble growth.

    The caller calls :func:`update` once per emitted candidate; the
    state tracks (a) the running Pareto frontier and (b) the count of
    consecutive iterations without a new frontier point.  When the
    plateau counter exceeds :attr:`patience` the caller should stop.
    """

    patience: int = DEFAULT_ADAPTIVE_PATIENCE
    max_ensemble: int = DEFAULT_MAX_ENSEMBLE
    mol: Any = None
    axes: Tuple[str, ...] = PARETO_AXES

    frontier_vecs: List[Tuple[float, ...]] = None  # type: ignore[assignment]
    n_seen: int = 0
    plateau_run: int = 0

    def __post_init__(self) -> None:
        if self.frontier_vecs is None:
            self.frontier_vecs = []

    def update(self, candidate: Any) -> bool:
        """Ingest one candidate; return ``True`` iff frontier grew.

        Side-effects:
          * increments ``n_seen``;
          * updates ``frontier_vecs`` (with pruning of newly-dominated);
          * resets / increments ``plateau_run``.
        """
        self.n_seen += 1
        vec = candidate_axes(candidate, mol=self.mol, axes=self.axes)

        # Is the new point dominated by anyone in the current frontier?
        for fv in self.frontier_vecs:
            if dominates(fv, vec):
                self.plateau_run += 1
                return False

        # Prune frontier members now dominated by the new point.
        kept: List[Tuple[float, ...]] = []
        for fv in self.frontier_vecs:
            if dominates(vec, fv):
                continue
            kept.append(fv)
        kept.append(vec)
        self.frontier_vecs = kept
        self.plateau_run = 0
        return True

    def should_stop(self) -> bool:
        """Adaptive stop: plateau reached OR max_ensemble hit."""
        if self.n_seen >= int(self.max_ensemble):
            return True
        if self.plateau_run >= int(self.patience):
            return True
        return False


def adaptive_plateau_check(
    candidates_seen: Sequence[Any],
    patience: int = DEFAULT_ADAPTIVE_PATIENCE,
    mol=None,
    axes: Sequence[str] = PARETO_AXES,
) -> bool:
    """Convenience wrapper: ``True`` iff the last ``patience`` candidates
    contributed no new frontier point.

    Pure / deterministic / does NOT mutate input.
    """
    if len(candidates_seen) < int(patience):
        return False
    state = AdaptiveState(
        patience=int(patience),
        max_ensemble=10 ** 9,  # effectively unbounded for the check
        mol=mol,
        axes=tuple(axes),
    )
    for c in candidates_seen:
        state.update(c)
    return state.plateau_run >= int(patience)
