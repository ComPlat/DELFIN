"""delfin.fffree.build_time_clash_gate — Per-candidate collapse / clash gate.

Construction-side fix #3 (User mandate 2026-06-03): structqual_xh_collapse
+154% and heavy_collapse +/- 0 still fire in voll-pool because the
construction phase emits candidates that already contain X-H or heavy-heavy
collapses (ETKDG numerics + degenerate metallacycle embeds), and the existing
`_has_collapsed_heavy_bonds` guard is env-gated default-OFF.

This module exposes a UNIVERSAL build-time check that mirrors the structqual
criteria EXACTLY (so a candidate that passes the gate cannot fail structqual):

  1. SUPERPOSITION   any-pair distance < SUPERPOS_TOL (0.5 A).
  2. HEAVY_COLLAPSE  any heavy-heavy bonded pair < HH_COLLAPSE_RATIO (0.70) x
                     element-pair ideal length.
  3. XH_COLLAPSE     any X-H bonded pair < XH_COLLAPSE_TOL (0.70 A).

Geometry-only, deterministic, no SMILES.  Env-gated default-OFF byte-identical.

API:
    has_collapse(syms, P) -> bool
    collapse_count(syms, P) -> int       # total violation count for ranking
    collapse_details(syms, P) -> list[dict]   # per-violation records
    apply_build_gate(syms, P, candidates_iter) -> chosen_P | None

Env flag:
    DELFIN_FFFREE_BUILD_CLASH_GATE=1    (per-fix activation)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (master activation)
"""
from __future__ import annotations

import os
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Calibrated thresholds: identical to structure_quality.py so the build-gate
# fires on EXACTLY the cases structqual would fail.
SUPERPOS_TOL = 0.5
HH_COLLAPSE_RATIO = 0.70
XH_COLLAPSE_TOL = 0.70
_BOND_FACTOR = 1.30


def _flag_active() -> bool:
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_BUILD_CLASH_GATE", "0") == "1"


def _min_pair_distance(P: np.ndarray) -> float:
    """Minimum atom-atom distance in the structure (any pair)."""
    n = len(P)
    if n < 2:
        return float("inf")
    mn = float("inf")
    for i in range(n):
        if i + 1 >= n:
            break
        d = np.linalg.norm(P[i + 1:] - P[i], axis=1)
        if d.size:
            md = float(d.min())
            if md < mn:
                mn = md
    return mn


def _geometric_bonds(syms: Sequence[str], P: np.ndarray
                     ) -> List[Tuple[int, int]]:
    """All geometric bonds (including X-H) via covalent-radii sum cutoff.

    Excludes metal bonds (matches structure_quality behaviour).
    """
    n = len(syms)
    out: List[Tuple[int, int]] = []
    for i in range(n):
        if _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                out.append((i, j))
    return out


def has_collapse(syms: Sequence[str], P: np.ndarray) -> bool:
    """True iff ANY of the 3 structqual collapse criteria fires.

    Env-gated default-OFF: returns False when no flag set (no rejection ->
    byte-identical behaviour).
    """
    if not _flag_active():
        return False
    if len(syms) == 0:
        return False
    P = np.asarray(P, dtype=float)
    # 1. superposition
    if _min_pair_distance(P) < SUPERPOS_TOL:
        return True
    # 2 + 3 via geometric bonds
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            # X-H bond
            if d < XH_COLLAPSE_TOL:
                return True
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                return True
    return False


def collapse_count(syms: Sequence[str], P: np.ndarray) -> int:
    """Total number of collapse-violation atom-pairs (used for ranking).

    Always counts (no env flag) -- this is read-only.  Callers gate on
    `has_collapse` and use `collapse_count` to pick the candidate with the
    FEWEST violations when none is fully clean.
    """
    if len(syms) == 0:
        return 0
    P = np.asarray(P, dtype=float)
    n_pairs = 0
    # 1. superposition
    if _min_pair_distance(P) < SUPERPOS_TOL:
        # count all overlapping pairs
        n = len(P)
        for i in range(n):
            for j in range(i + 1, n):
                if float(np.linalg.norm(P[i] - P[j])) < SUPERPOS_TOL:
                    n_pairs += 1
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            if d < XH_COLLAPSE_TOL:
                n_pairs += 1
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                n_pairs += 1
    return n_pairs


def collapse_details(syms: Sequence[str], P: np.ndarray) -> List[dict]:
    """Per-violation records (used for tests/debug)."""
    if len(syms) == 0:
        return []
    P = np.asarray(P, dtype=float)
    out: List[dict] = []
    # 1. superposition: only any-pair distances under tol
    n = len(P)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < SUPERPOS_TOL:
                out.append({"mode": "superposition", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4)})
    for i, j in _geometric_bonds(syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        if syms[i] == "H" or syms[j] == "H":
            if d < XH_COLLAPSE_TOL:
                out.append({"mode": "xh_collapse", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4)})
        else:
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d < HH_COLLAPSE_RATIO * ideal:
                out.append({"mode": "heavy_collapse", "i": i, "j": j,
                            "a": syms[i], "b": syms[j], "d": round(d, 4),
                            "ratio": round(d / max(ideal, 1e-9), 4)})
    return out


def select_clean_candidate(
    syms_per_candidate: Sequence[Sequence[str]],
    P_per_candidate: Sequence[np.ndarray],
) -> Optional[int]:
    """Pick the index of the first FULLY-CLEAN candidate (collapse_count==0),
    falling back to the candidate with the FEWEST collapse violations.

    Returns None only when there are no candidates.

    Env-gated default-OFF: when no flag, returns 0 (first candidate)
    unconditionally -> byte-identical to picking the first candidate.
    """
    if not syms_per_candidate or not P_per_candidate:
        return None
    if not _flag_active():
        return 0
    best_idx = 0
    best_count = collapse_count(syms_per_candidate[0], P_per_candidate[0])
    if best_count == 0:
        return 0
    for k in range(1, len(syms_per_candidate)):
        c = collapse_count(syms_per_candidate[k], P_per_candidate[k])
        if c < best_count:
            best_count = c
            best_idx = k
            if best_count == 0:
                return best_idx
    return best_idx


__all__ = [
    "has_collapse", "collapse_count", "collapse_details",
    "select_clean_candidate", "_flag_active",
    "SUPERPOS_TOL", "HH_COLLAPSE_RATIO", "XH_COLLAPSE_TOL",
]
