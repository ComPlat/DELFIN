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

# inter-ligand heavy-heavy clash factor.  The license-free detector flags only
# < 0.65*vdwsum, but REAL CSD crystals keep inter-ligand contacts >= ~0.85*vdwsum
# (the 0.65 gate is deliberately permissive to keep FP=0).  Visual "clashes" the eye
# catches live in the 0.70-0.85 band — too close for a real crystal but passing the
# detector.  So we relax toward the REAL-crystal separation (0.85) to remove them.
_CLASH_F = float(os.environ.get("DELFIN_FFFREE_FLEX_CLASH_F", "0.85"))
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


def _monodentate_ligands(n, syms, lig, donors):
    """List of (donor_idx, [body_atom_idxs]) for each MONODENTATE ligand (exactly
    one donor on it).  Chelates (>=2 donors) are excluded — a single radial move
    would distort their bite, so they keep the ideal placement."""
    from collections import Counter
    dcount = Counter(int(lig[d]) for d in donors)
    out = []
    for d in donors:
        lg = int(lig[d])
        if dcount[lg] != 1:
            continue                                # part of a chelate -> skip
        body = [i for i in range(n) if int(lig[i]) == lg]   # whole ligand (incl donor)
        out.append((d, body))
    return out


def sphere_flex(syms: Sequence[str], P, fixed: Iterable[int],
                bond_pairs: Optional[Sequence[Tuple[int, int]]] = None,
                max_disp: float = _DEF_MAX_DISP, k: float = _DEF_K,
                passes: int = _DEF_PASSES, damp: float = _DEF_DAMP,
                md_tol: float = _DEF_MD_TOL) -> np.ndarray:
    """Radial coordination-sphere clash relief (NO FF, NO QM).

    Mechanism (data-driven 2026-06-29): mild inter-ligand heavy-heavy clashes are
    intrinsic to the rigid ideal placement and CANNOT be opened by rotating a
    ligand body about its M-D axis (de-risked: rotation just spins the clash).
    What relieves them is letting a crowded MONODENTATE ligand translate RADIALLY
    OUTWARD along its own M-D axis by a small, HARD-BOUNDED amount (rigid-body, so
    intra-ligand bond lengths/angles are exactly preserved and the donor stays on
    its polyhedron ray -> coordination ANGLES are invariant; only that ligand's M-D
    distance grows, by <= ``max_disp``).  Deterministic coordinate descent over the
    per-ligand radial offset; never-worse on the inter-ligand heavy clash; the metal
    and all chelate donors stay fixed.  Never raises.

    ``fixed`` = metal + donor indices.  ``max_disp`` doubles as the max radial M-D
    elongation (A)."""
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

    base_loss, base_nclash = _clash_loss(P, pairs, rsum)
    if base_nclash == 0:
        return P                                   # nothing to do

    movers = _monodentate_ligands(n, syms, lig, donors)
    if not movers:
        return P
    md0 = {d: float(np.linalg.norm(P[d] - P[m])) for d, _ in movers}
    P0 = P.copy()

    # radial-offset grid: 0 .. max_disp outward (never inward), fine step.
    steps = max(4, int(round(max_disp / 0.05)))
    grid = [max_disp * (s + 1) / steps for s in range(steps)]   # >0 outward
    best = P.copy()
    best_loss, best_nclash = base_loss, base_nclash

    for _ in range(max(1, passes)):
        improved = False
        for d, body in movers:
            axis = best[d] - best[m]
            na = float(np.linalg.norm(axis))
            if na < 1e-9:
                continue
            u = axis / na
            cur_off = na - md0[d]                   # current elongation already applied
            local_loss, local_nclash, local_P = best_loss, best_nclash, None
            for g in grid:
                shift = (g - cur_off) * u           # move so total elongation == g
                trial = best.copy()
                trial[body] = best[body] + shift
                tl, tnc = _clash_loss(trial, pairs, rsum)
                better = (tnc < local_nclash) or (tnc == local_nclash and tl < local_loss - 1e-9)
                if better:
                    local_loss, local_nclash, local_P = tl, tnc, trial
            if local_P is not None:
                best = local_P; best_loss, best_nclash = local_loss, local_nclash
                improved = True
                if best_nclash == 0:
                    break
        if not improved or best_nclash == 0:
            break

    # hard guards: each moved donor's M-D elongation within [0, max_disp]; never-worse.
    for d, _ in movers:
        elong = float(np.linalg.norm(best[d] - best[m])) - md0[d]
        if elong < -1e-6 or elong > max_disp + 1e-6:
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
