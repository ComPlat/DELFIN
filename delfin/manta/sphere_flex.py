"""delfin.manta.sphere_flex — soft coordination-sphere clash relaxation (NO FF, NO QM).

Why this exists (data-driven, 2026-06-29 clash forensik)
--------------------------------------------------------
The rigid ideal-polyhedron placement seats every donor EXACTLY on its vertex at
the ideal M-D distance and then ``refine`` relieves clashes with the metal + ALL
donors FROZEN.  Measured against real CCDC crystals, the FF-free build has a
residual inter-ligand heavy-heavy clash gap: ~48% of structures carry a worst
overlap > 0.5 A vs ~17% for real crystals (mean 0.48 vs 0.28 A) — MOSTLY MILD
(only ~6% gross > 1.0 A).  A whole-body M-D-axis rotation CANNOT relieve these
(empirically de-risked: rotating a crowded ligand body just spins the clash
around — the overlap is intrinsic to the rigid sphere, not the body orientation).

What relieves a mild intrinsic clash is letting the coordination sphere BREATHE a
little — exactly what UFF does (and why legacy sometimes looked cleaner) — but
UFF collapses the metal topology.  This pass lets the donors deviate by a SMALL,
HARD-BOUNDED amount under a harmonic restoring spring toward their ideal vertex
position AND ideal M-D distance, so the sphere flexes just enough to open the
inter-ligand contacts while staying within a few hundredths of an Angstrom of the
ideal polyhedron.  Pure geometry (Bondi vdW repulsion + harmonic restraints),
deterministic, never-worse-gated, hard-bounded; the metal stays fixed at origin.

Env-gated ``DELFIN_FFFREE_SPHERE_FLEX`` (default ``"0"`` => byte-identical: never
invoked).  Runs as a post-placement pass on each assembled frame, after refine /
torsion_relax / joint_declash and before the self-gate.  Bounds + tolerances are
env-overridable but clamped.  Never raises (returns the input frame on any error).
"""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin.manta._bond_decollapse as _bd
from delfin.manta.refine import _vdw

# inter-ligand heavy-heavy clash factor — matches the self-gate / detector intent
# (detector flags < 0.65*vdwsum; we relax toward >= 0.70 to clear it with margin).
_CLASH_F = 0.70
_DEF_MAX_DISP = 0.30     # A, hard cap on per-donor displacement from its ideal vertex
_DEF_K = 1.5             # harmonic restoring strength (toward ideal vertex + ideal M-D)
_DEF_PASSES = 60
_DEF_DAMP = 0.5
_DEF_MD_TOL = 0.30       # A, max tolerated M-D drift (also bounded by max_disp)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_SPHERE_FLEX", "0") == "1"


def _env_float(name, default, lo, hi):
    try:
        v = float(os.environ.get(name, default))
    except Exception:
        return default
    return max(lo, min(hi, v))


def _env_int(name, default, lo, hi):
    try:
        v = int(os.environ.get(name, default))
    except Exception:
        return default
    return max(lo, min(hi, v))


def _ligand_of_atom(syms: Sequence[str], P: np.ndarray,
                    bond_pairs: Optional[Sequence[Tuple[int, int]]]) -> np.ndarray:
    """Per-atom ligand id; metal = -1.  Components of the bond graph with the
    metal-donor coordination bonds removed (so ligands stay separate)."""
    n = len(syms)
    if bond_pairs is not None:
        adj: List[List[int]] = [[] for _ in range(n)]
        for i, j in bond_pairs:
            adj[i].append(j); adj[j].append(i)
    else:
        b = _bd._geometric_bonds(syms, P)
        adj = [[] for _ in range(n)]
        for i, j in b:
            adj[i].append(j); adj[j].append(i)
    metals = {i for i in range(n) if _bd._is_metal(syms[i])}
    lig = np.full(n, -1, dtype=int)
    cur = 0
    for start in range(n):
        if start in metals or lig[start] != -1:
            continue
        stack = [start]; lig[start] = cur
        while stack:
            a = stack.pop()
            for nb in adj[a]:
                if nb in metals or lig[nb] != -1:
                    continue
                lig[nb] = cur; stack.append(nb)
        cur += 1
    return lig


def _inter_heavy_pairs(syms: Sequence[str], lig: np.ndarray) -> np.ndarray:
    """(K,2) array of inter-ligand heavy-heavy atom-index pairs (neither H, neither
    metal, different ligand)."""
    n = len(syms)
    heavy = [i for i in range(n)
             if syms[i] != "H" and not _bd._is_metal(syms[i]) and lig[i] >= 0]
    pairs = []
    for a in range(len(heavy)):
        for b in range(a + 1, len(heavy)):
            i, j = heavy[a], heavy[b]
            if lig[i] != lig[j]:
                pairs.append((i, j))
    return np.array(pairs, dtype=int) if pairs else np.zeros((0, 2), dtype=int)


def _clash_loss(P: np.ndarray, pairs: np.ndarray, rsum: np.ndarray) -> Tuple[float, int]:
    """(sum of squared overlaps, n_clashing) over the inter-ligand heavy pairs."""
    if pairs.shape[0] == 0:
        return 0.0, 0
    d = np.linalg.norm(P[pairs[:, 0]] - P[pairs[:, 1]], axis=1)
    over = rsum - d
    over = np.where(over > 0.0, over, 0.0)
    return float((over * over).sum()), int((over > 0.0).sum())


def sphere_flex(syms: Sequence[str], P, fixed: Iterable[int],
                bond_pairs: Optional[Sequence[Tuple[int, int]]] = None,
                max_disp: float = _DEF_MAX_DISP, k: float = _DEF_K,
                passes: int = _DEF_PASSES, damp: float = _DEF_DAMP,
                md_tol: float = _DEF_MD_TOL) -> np.ndarray:
    """Soft-sphere inter-ligand clash relaxation.  ``fixed`` = metal + donor indices
    (metal stays frozen; donors are soft-restrained, not frozen).  Returns relaxed
    P; never-worse on the inter-ligand heavy clash count, hard-bounded donor drift,
    deterministic.  Never raises."""
    try:
        P = np.array(P, dtype=float).copy()
    except Exception:
        return np.array(P, dtype=float)
    n = len(syms)
    if n < 3 or P.shape != (n, 3) or not np.all(np.isfinite(P)):
        return P
    fixed = set(int(x) for x in fixed)
    metals = [i for i in range(n) if _bd._is_metal(syms[i])]
    if len(metals) != 1:
        return P
    m = metals[0]
    donors = sorted(d for d in fixed if d != m)
    if not donors:
        return P

    lig = _ligand_of_atom(syms, P, bond_pairs)
    pairs = _inter_heavy_pairs(syms, lig)
    if pairs.shape[0] == 0:
        return P
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[pairs[:, 0]] + vdw[pairs[:, 1]])

    P0 = P.copy()                                  # ideal donor positions (= vertices)
    md0 = {d: float(np.linalg.norm(P[d] - P[m])) for d in donors}
    base_loss, base_nclash = _clash_loss(P, pairs, rsum)
    if base_nclash == 0:
        return P                                   # nothing to do
    donor_set = set(donors)

    cur_damp = damp
    best = P.copy()
    best_loss, best_nclash = base_loss, base_nclash

    # geometric clash-repulsion + harmonic donor restraint coordinate descent.
    from delfin.manta.refine import _bonds_adj, _violations
    for _ in range(max(1, passes)):
        _, bonded, adj = _bonds_adj(syms, P)
        _loss_unused, moves = _violations(syms, P, bonded, adj)
        disp = np.zeros((n, 3)); cnt = np.zeros(n)
        for idx, dvec in moves:
            if idx == m:                           # metal stays put
                continue
            disp[idx] += dvec; cnt[idx] += 1
        for i in range(n):
            if cnt[i] > 0:
                disp[i] /= cnt[i]
        # harmonic restraint on donors: pull back toward ideal vertex AND ideal M-D
        for d in donors:
            disp[d] += k * (P0[d] - P[d])          # spring to ideal vertex position
            v = P[d] - P[m]; nv = float(np.linalg.norm(v))
            if nv > 1e-9:
                disp[d] += k * (md0[d] - nv) * (v / nv)   # spring to ideal M-D length
        trial = P + cur_damp * disp
        # hard bound: clamp each donor to within max_disp of its ideal vertex
        for d in donors:
            v = trial[d] - P0[d]; nv = float(np.linalg.norm(v))
            if nv > max_disp:
                trial[d] = P0[d] + v / nv * max_disp
        if not np.all(np.isfinite(trial)):
            cur_damp *= 0.6
            if cur_damp < 0.03:
                break
            continue
        tloss, tnclash = _clash_loss(trial, pairs, rsum)
        # accept if the inter-ligand heavy clash strictly improves (fewer clashes,
        # or same count but lower overlap energy)
        better = (tnclash < best_nclash) or (tnclash == best_nclash and tloss < best_loss - 1e-9)
        if better:
            P = trial; best = trial.copy(); best_loss, best_nclash = tloss, tnclash
            if best_nclash == 0:
                break
        else:
            cur_damp *= 0.6
            if cur_damp < 0.03:
                break

    # final hard guards: donor drift / M-D drift within bound, never-worse on clash
    for d in donors:
        if float(np.linalg.norm(best[d] - P0[d])) > max_disp + 1e-6:
            return P0
        if abs(float(np.linalg.norm(best[d] - best[m])) - md0[d]) > md_tol + 1e-6:
            return P0
    if not np.all(np.isfinite(best)):
        return P0
    fin_loss, fin_nclash = _clash_loss(best, pairs, rsum)
    if fin_nclash > base_nclash or (fin_nclash == base_nclash and fin_loss > base_loss + 1e-9):
        return P0                                  # never-worse
    return best


def flex_if_enabled(syms: Sequence[str], P, fixed: Iterable[int],
                    bond_pairs: Optional[Sequence[Tuple[int, int]]] = None):
    """Wire-in entry: when ``DELFIN_FFFREE_SPHERE_FLEX=1`` run :func:`sphere_flex`,
    else return ``P`` unchanged (byte-identical default-OFF).  Never raises."""
    if not _enabled():
        return P
    try:
        max_disp = _env_float("DELFIN_FFFREE_FLEX_MAX_DISP", _DEF_MAX_DISP, 0.0, 1.0)
        k = _env_float("DELFIN_FFFREE_FLEX_K", _DEF_K, 0.1, 10.0)
        passes = _env_int("DELFIN_FFFREE_FLEX_PASSES", _DEF_PASSES, 1, 200)
        md_tol = _env_float("DELFIN_FFFREE_FLEX_MD_TOL", _DEF_MD_TOL, 0.0, 2.0)
        return sphere_flex(syms, P, fixed, bond_pairs=bond_pairs,
                           max_disp=max_disp, k=k, passes=passes, md_tol=md_tol)
    except Exception:
        return P
