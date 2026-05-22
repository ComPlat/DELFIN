"""delfin.fffree.refine — COD-loss refiner (NO force field).

Deterministic coordinate descent that directly minimises the structural-defect
loss the battery measures (clashes / H-overcoordination / bond collapse), moving
only the ligand periphery with the metal + donor atoms FROZEN.  Per-pass
accept-if-better gate with rollback (the _bond_decollapse pattern, generalised).

This is the mechanism that makes the geometric metal-FF-free placement clean
WITHOUT a force field on the metal: the placement gives correct topology +
M-D geometry; the refiner removes the residual clashes that rigid placement
leaves, which is exactly where UFF-relaxed structures win.
"""
from __future__ import annotations
from typing import List, Set, Tuple
import numpy as np
import delfin._bond_decollapse as bd

_VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10, "Se": 1.90,
        "As": 1.85, "Te": 2.06}
_VDW_D = 1.7
_CLASH = 0.70          # non-bonded pair < this * (vdW_i+vdW_j) = clash
_COLLAPSE = 0.70       # bond < this * ideal = collapsed


def _vdw(s):
    return _VDW.get(s, _VDW_D)


def _bonds_adj(syms, P):
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    adj: List[List[int]] = [[] for _ in syms]
    for i, j in b:
        adj[i].append(j); adj[j].append(i)
    return b, bonded, adj


def _violations(syms, P, bonded, adj):
    """Return (loss, moves) where moves = list of (atom_idx, displacement) that
    would relieve a violation.  loss = #clashes + #collapses + #h_overcoord."""
    n = len(syms)
    loss = 0
    moves: List[Tuple[int, np.ndarray]] = []
    # clashes (non-bonded, incl. H)
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1e-6:
                continue
            tgt = _CLASH * (_vdw(syms[i]) + _vdw(syms[j]))
            if d < tgt:
                loss += 1
                u = (P[i] - P[j]) / d
                push = 0.5 * (tgt - d)
                moves.append((i, +push * u)); moves.append((j, -push * u))
    # collapsed bonds -> stretch to ideal
    for i, j in bonded:
        d = float(np.linalg.norm(P[i] - P[j]))
        ideal = bd._ideal_bond(syms[i], syms[j])
        if d < _COLLAPSE * ideal and d > 1e-6:
            loss += 1
            u = (P[j] - P[i]) / d
            stretch = 0.5 * (ideal - d)
            moves.append((j, +stretch * u)); moves.append((i, -stretch * u))
    # H over-coordination -> push H off the non-parent heavy
    for h in range(n):
        if syms[h] != "H":
            continue
        heavies = [k for k in range(n) if k != h and syms[k] != "H"
                   and not bd._is_metal(syms[k])
                   and np.linalg.norm(P[h] - P[k]) < 1.3 * bd._ideal_bond(syms[h], syms[k])]
        if len(heavies) >= 2:
            loss += 1
            parent = min(heavies, key=lambda k: np.linalg.norm(P[h] - P[k]))
            for k in heavies:
                if k == parent:
                    continue
                d = float(np.linalg.norm(P[h] - P[k]))
                if d > 1e-6:
                    moves.append((h, 0.4 * (P[h] - P[k]) / d))
    return loss, moves


def refine(syms, P, fixed_idx: Set[int], max_passes: int = 40, damp: float = 0.6):
    """Coordinate descent minimising the defect loss; metal+donors frozen.
    Deterministic. Returns refined P (or original if no improvement)."""
    P = np.array(P, float).copy()
    fixed = set(fixed_idx)
    n = len(syms)
    _, bonded, adj = _bonds_adj(syms, P)
    best_loss, _ = _violations(syms, P, bonded, adj)
    if best_loss == 0:
        return P
    for _ in range(max_passes):
        _, bonded, adj = _bonds_adj(syms, P)
        loss, moves = _violations(syms, P, bonded, adj)
        if loss == 0:
            break
        disp = np.zeros((n, 3))
        cnt = np.zeros(n)
        for idx, d in moves:
            if idx in fixed:
                continue
            disp[idx] += d; cnt[idx] += 1
        for i in range(n):
            if cnt[i] > 0:
                disp[i] /= cnt[i]
        trial = P + damp * disp
        _, tb, ta = _bonds_adj(syms, trial)
        tloss, _ = _violations(syms, trial, tb, ta)
        if tloss < best_loss:          # accept-if-better gate
            P = trial; best_loss = tloss
        else:
            damp *= 0.6                # smaller step
            if damp < 0.05:
                break
    return P
