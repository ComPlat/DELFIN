"""delfin.fffree.ligand_relax — COD-empirical ligand-internal relaxer (#38).

THE GRUNDPFEILER (division of labor): fffree fixes the metal coordination sphere
BY CONSTRUCTION (enumeration + ideal polyhedra, metal-FF-free) and relaxes ONLY the
ligand internals with a COD-empirical organic potential, metal + donor atoms FROZEN.
This is the inverse of UFF (which relaxes everything with a poor metal force field).

The existing `refine.py` minimises a discrete DEFECT COUNT (clashes/collapses).  This
module adds the CONTINUOUS COD loss that the placement alone cannot supply and that is
exactly where UFF wins today (F23 funcgrp internal +261%, F24 amide +56%, lig_realistic
-34% on the clean 75ce258 baseline):

  loss = w_bond·Σ(d_ij - r0)²            bonds  -> COD-empirical r0 (bd._ideal_bond)
       + w_ang ·Σ(d_ik - d0(θ0))²        angles -> Urey-Bradley 1-3 spring at θ0 by
                                                   steric number (109.47/120/180)
       + w_plan·Σ oop²                   sp2-like residual planarity (amide/aromatic/
                                                   carbonyl) -> drives residual pucker
                                                   of ALREADY-near-planar 3-coord centres
                                                   to 0 (does NOT flatten genuine sp3 —
                                                   guards against the #35 sp3-N over-
                                                   pyramidalisation regression)
       + w_clash·Σ(tgt - d)²             non-bonded clash penalty (keeps the relaxation
                                                   from trading internal geometry for a
                                                   new steric overlap)

Deterministic coordinate descent (no RNG, fixed damping schedule, accept-if-loss-
decreases gate).  Metal + donor atoms frozen.  Hard ‖M-D‖ rollback, per-step clamp,
finite guard.  Geometry-only + (where available) graph hybridisation -> universal,
atom-order-independent.  Disabled by default; enable via DELFIN_FFFREE_LIGAND_RELAX=1.
"""
from __future__ import annotations
from typing import List, Set, Tuple
import numpy as np
import delfin._bond_decollapse as bd

# Bondi vdW radii (match the vdw_clash detector exactly via bd._VDW).
_VDW = bd._VDW
_VDW_D = bd._VDW_DEFAULT
_CLASH_FACTOR = 0.85            # clash if d < factor·(vdW_i+vdW_j)  (matches detector)

# loss weights (calibration knobs — tuned against the F23/F24/F3/vdw detectors on smoke)
_W_BOND = 1.0
_W_ANG = 0.6
_W_PLAN = 1.2
_W_CLASH = 2.0

# steric-number -> ideal bond angle (deg).  2=linear, 3=trigonal, 4=tetrahedral.
_THETA0 = {2: 180.0, 3: 120.0, 4: 109.47}

_MAX_STEP = 0.30               # per-atom per-pass clamp (Å) — iter26 NaN lesson
_PLANAR_OOP_TOL = 0.50         # only planarise 3-coord centres already within this (Å)


def _vdw(s: str) -> float:
    return _VDW.get(s, _VDW_D)


def _bonds_adj(syms, P):
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    adj: List[List[int]] = [[] for _ in syms]
    for i, j in b:
        adj[i].append(j); adj[j].append(i)
    return bonded, adj


def _loss_and_moves(syms, P, bonded, adj):
    """Continuous COD loss (scalar) + per-atom displacement moves (steepest-descent
    style, accumulated like refine.py).  Pure function, deterministic."""
    n = len(syms)
    loss = 0.0
    moves: List[Tuple[int, np.ndarray]] = []

    # --- bond term: pull each bond toward COD-empirical r0 ---
    for i, j in bonded:
        d = float(np.linalg.norm(P[i] - P[j]))
        if d <= 1e-6:
            continue
        r0 = bd._ideal_bond(syms[i], syms[j])
        dev = d - r0
        loss += _W_BOND * dev * dev
        u = (P[j] - P[i]) / d
        corr = 0.5 * _W_BOND * (-dev)        # +stretch if short, -compress if long
        moves.append((j, +corr * u)); moves.append((i, -corr * u))

    # --- angle term: Urey-Bradley 1-3 spring toward θ0(steric number) ---
    for j in range(n):
        if bd._is_metal(syms[j]):
            continue
        nbrs = adj[j]
        if len(nbrs) < 2:
            continue
        theta0 = _THETA0.get(len(nbrs))
        if theta0 is None:
            continue
        ct = np.cos(np.radians(theta0))
        for a_i in range(len(nbrs)):
            for b_i in range(a_i + 1, len(nbrs)):
                a, b = nbrs[a_i], nbrs[b_i]
                rja = float(np.linalg.norm(P[a] - P[j]))
                rjb = float(np.linalg.norm(P[b] - P[j]))
                d_ab = float(np.linalg.norm(P[a] - P[b]))
                if rja <= 1e-6 or rjb <= 1e-6 or d_ab <= 1e-6:
                    continue
                d0 = np.sqrt(max(rja * rja + rjb * rjb - 2 * rja * rjb * ct, 1e-6))
                dev = d_ab - d0
                loss += _W_ANG * dev * dev
                u = (P[a] - P[b]) / d_ab
                corr = 0.5 * _W_ANG * (-dev)
                moves.append((a, +corr * u)); moves.append((b, -corr * u))

    # --- planarity term: residual de-pucker of ALREADY-near-planar 3-coord centres ---
    # (sp2-like: amide N, carbonyl/aromatic C).  Genuinely pyramidal sp3 (large oop) is
    # left untouched -> guards against #35 sp3-N over-pyramidalisation.
    for j in range(n):
        if syms[j] == "H" or bd._is_metal(syms[j]):
            continue
        nbrs = [k for k in adj[j] if not bd._is_metal(syms[k])]
        if len(nbrs) != 3:
            continue
        a, b, c = (P[k] for k in nbrs)
        nrm = np.cross(b - a, c - a)
        ln = float(np.linalg.norm(nrm))
        if ln <= 1e-6:
            continue
        nrm = nrm / ln
        oop = float(np.dot(P[j] - a, nrm))      # signed out-of-plane distance of centre
        if abs(oop) > _PLANAR_OOP_TOL:           # genuine pyramid -> leave alone
            continue
        loss += _W_PLAN * oop * oop
        corr = 0.5 * _W_PLAN * (-oop)
        moves.append((j, corr * nrm))

    # --- clash penalty: non-bonded heavy/H pairs inside the vdW clash radius ---
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= 1e-6:
                continue
            tgt = _CLASH_FACTOR * (_vdw(syms[i]) + _vdw(syms[j]))
            if d < tgt:
                dev = tgt - d
                loss += _W_CLASH * dev * dev
                u = (P[i] - P[j]) / d
                push = 0.5 * _W_CLASH * dev
                moves.append((i, +push * u)); moves.append((j, -push * u))

    return loss, moves


def relax(syms, P, fixed_idx: Set[int], max_passes: int = 120, damp: float = 0.5):
    """COD-loss coordinate descent; metal + donor atoms (fixed_idx) FROZEN.
    Deterministic. Accept-if-loss-decreases with damping; per-step clamp + finite +
    ‖M-D‖ guard.  Returns the relaxed P, or the original P if nothing improved."""
    P = np.array(P, float).copy()
    if P.size == 0 or not np.all(np.isfinite(P)):
        return P
    fixed = set(int(i) for i in fixed_idx)
    n = len(syms)
    md_snap = bd._md_snapshot(syms, P)
    bonded, adj = _bonds_adj(syms, P)
    best_loss, _ = _loss_and_moves(syms, P, bonded, adj)
    if best_loss <= 1e-9:
        return P
    for _ in range(max_passes):
        bonded, adj = _bonds_adj(syms, P)
        loss, moves = _loss_and_moves(syms, P, bonded, adj)
        disp = np.zeros((n, 3)); cnt = np.zeros(n)
        for idx, dvec in moves:
            if idx in fixed:
                continue
            disp[idx] += dvec; cnt[idx] += 1
        for i in range(n):
            if cnt[i] > 0:
                disp[i] /= cnt[i]
        step = damp * disp
        # per-atom clamp (iter26 NaN lesson) — no atom moves more than _MAX_STEP/pass
        norms = np.linalg.norm(step, axis=1)
        big = norms > _MAX_STEP
        if np.any(big):
            step[big] *= (_MAX_STEP / norms[big])[:, None]
        trial = P + step
        if not np.all(np.isfinite(trial)):          # finite guard
            break
        if bd._md_broken(md_snap, trial):           # ‖M-D‖ ±tol hard rollback
            damp *= 0.6
            if damp < 0.03:
                break
            continue
        tb, ta = _bonds_adj(syms, trial)
        tloss, _ = _loss_and_moves(syms, trial, tb, ta)
        if tloss < best_loss - 1e-9:                # accept-if-better gate
            P = trial; best_loss = tloss
        else:
            damp *= 0.6
            if damp < 0.03:
                break
    return P


if __name__ == "__main__":
    # smoke self-test: a deliberately puckered planar amide + a stretched bond.
    # relax must (a) reduce loss, (b) keep metal+donor frozen, (c) be deterministic,
    # (d) not break ‖M-D‖.
    syms = ["Fe", "N", "C", "O", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],     # Fe (metal, frozen)
        [2.0, 0.0, 0.0],     # N  donor (frozen)
        [3.0, 0.6, 0.4],     # C  (puckered out of N-O-H plane)
        [3.2, 1.6, 0.0],     # O
        [2.6, -0.7, 0.0],    # H
    ], float)
    fixed = {0, 1}
    out1 = relax(syms, P, fixed)
    out2 = relax(syms, P, fixed)
    b0, a0 = _bonds_adj(syms, P)
    l0, _ = _loss_and_moves(syms, P, b0, a0)
    b1, a1 = _bonds_adj(syms, out1)
    l1, _ = _loss_and_moves(syms, out1, b1, a1)
    print(f"loss {l0:.4f} -> {l1:.4f}  (reduced: {l1 < l0})")
    print(f"metal frozen: {np.allclose(out1[0], P[0])}  donor frozen: {np.allclose(out1[1], P[1])}")
    print(f"deterministic: {np.allclose(out1, out2)}")
    print(f"M-D intact: {not bd._md_broken(bd._md_snapshot(syms, P), out1)}")
    print(f"finite: {np.all(np.isfinite(out1))}")
