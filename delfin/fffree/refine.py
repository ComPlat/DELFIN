"""delfin.fffree.refine — geometry refiner (NO force field).

Deterministic coordinate descent that directly minimises a structural-defect
count (close-contact clashes / hydrogen over-coordination / collapsed bonds),
moving only the ligand periphery with the metal + donor atoms FROZEN.  Each pass
is accepted only if the defect count decreases, otherwise the step is halved and
rolled back.

This makes the rigid geometric placement clean without a force field on the
metal: the placement already gives correct connectivity and metal-donor
distances; this step removes the residual close contacts that rigid placement
leaves behind.
"""
from __future__ import annotations
from typing import List, Set, Tuple
import os
import numpy as np
import delfin._bond_decollapse as bd
from delfin.fffree import cod_ideals as _CODI

# Track 3 pure FF-free mode (User 2026-05-30): when PURE_TRACK3=1, MMFF is removed
# from the fffree pipeline (Phases 1+2). To compensate AND BEAT the FF baseline, the
# pure-Track-3 mode auto-enables the FF-free defect correctors below (COD-empirical
# bond targets + rigid-H drag) as part of the integrated FF-free correction stack.
# When PURE_TRACK3 is off, individual flags remain independent (byte-identical
# behaviour preserved). Race-strategy: FF-free must win, not just match.
_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"

# Iter 29 (construction-driver): refine heavy-heavy bonds toward COD-empirical lengths
# (real crystal p50) instead of the generic covalent-sum _ideal_bond, which is
# systematically too long for aromatic/conjugated ligands (C-C 1.52 vs COD 1.40) and
# was the dominant F3_bond / fragfit gap of FF-free vs UFF.  Env-gated, default OFF
# (byte-identical when unset).  CCDC-ready: cod_ideals is swappable for a CSD source.
# Auto-enabled under PURE_TRACK3 to close the bondlen frac_outlier gap vs MMFF.
_USE_COD_BONDS = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_COD_BONDS", "0") == "1"

# Iter 30 (User 2026-05-28, RUJSIY eye-validation): when a heavy atom (C/N/O of a ring)
# is moved by the refiner, the H atoms bonded to it MUST move along (rigid-H drag),
# else the X–H bond stretches from 1.08 Å to ~1.3 Å as the heavy drifts away — user
# observed "in den frames bewegen sich nur die C und N atome im ring aber die H bleiben
# an ihrer stelle und dadurch entstehen unrealitische strukturen".  Env-gated default
# OFF (byte-identical when unset); the propagation copies the heavy's accumulated
# displacement onto its bonded H atoms before the step is applied.
# Auto-enabled under PURE_TRACK3 to close the H-anomaly gap vs MMFF.
_RIGID_H_DRAG = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_RIGID_H_DRAG", "0") == "1"


def _bond_ideal(a, b, use_cod=None):
    """Bond ideal length. use_cod=None uses module-level _USE_COD_BONDS;
    pass use_cod=False to force covalent-sum (used in Phase 3b hapto gate)."""
    if use_cod is None:
        use_cod = _USE_COD_BONDS
    if use_cod:
        v = _CODI.cod_ideal_bond(a, b)
        if v is not None:
            return v
    return bd._ideal_bond(a, b)

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


def _violations(syms, P, bonded, adj, use_cod=None):
    """Return (loss, moves) where moves = list of (atom_idx, displacement) that
    would relieve a violation.  loss = #clashes + #collapses + #h_overcoord.

    Phase 3b: use_cod parameter forwards to _bond_ideal for class-conditional gating
    (hapto systems pass use_cod=False to skip COD_BONDS auto-enable)."""
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
    # collapsed bonds -> stretch to ideal; distorted heavy-heavy bonds -> toward ideal
    for i, j in bonded:
        d = float(np.linalg.norm(P[i] - P[j]))
        if d <= 1e-6:
            continue
        ideal = _bond_ideal(syms[i], syms[j], use_cod=use_cod)
        if d < _COLLAPSE * ideal:
            loss += 1
            u = (P[j] - P[i]) / d
            stretch = 0.5 * (ideal - d)
            moves.append((j, +stretch * u)); moves.append((i, -stretch * u))
        elif syms[i] != "H" and syms[j] != "H":
            dev = (d - ideal) / ideal
            if dev < -0.25 or dev > 0.085:        # bond_distort band
                loss += 1
                u = (P[j] - P[i]) / d
                corr = 0.4 * (ideal - d)           # +stretch if short, -compress if long
                moves.append((j, +corr * u)); moves.append((i, -corr * u))
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


def _has_hapto_signature(syms, P) -> bool:
    """Heuristic hapto detection: metal with ≥3 bonded C atoms clustered at
    ~2.0-2.5 Å (Cp/allyl/π-system signature). Returns True if any metal carries a
    hapto motif. Phase 3b uses this to gate Phase 3a auto-enable: hapto systems
    skip COD_BONDS + RIGID_H_DRAG to avoid refine-loop oscillation on Cp rings
    (Sn-Cp case: 3/4 multi-hapto SMILES timed out under Phase 3a auto-enable).

    Pure geometry, no FF, universal across all CN.
    """
    n = len(syms)
    for m in range(n):
        if not bd._is_metal(syms[m]):
            continue
        c_neighbors = [k for k in range(n)
                       if k != m and syms[k] == "C"
                       and float(np.linalg.norm(P[k] - P[m])) < 2.7]
        if len(c_neighbors) < 3:
            continue
        # check cluster: spread of M-C distances small (hapto = ~equidistant)
        d_to_m = sorted(float(np.linalg.norm(P[k] - P[m])) for k in c_neighbors)
        if d_to_m[2] - d_to_m[0] < 0.4:
            return True
    return False


def refine(syms, P, fixed_idx: Set[int], max_passes: int = 80, damp: float = 0.6):
    """Coordinate descent minimising the defect loss; metal+donors frozen.
    Deterministic. Returns refined P (or original if no improvement).

    Phase 3b (User 2026-05-30): when PURE_TRACK3 + hapto signature detected,
    disable COD_BONDS + RIGID_H_DRAG locally to prevent refine-loop oscillation
    on Cp-type π systems (preserves the +1460% multi-hapto coverage breakthrough
    of Phase 1+2 while keeping the +4-axis Phase 3a quality win on σ systems).
    """
    P = np.array(P, float).copy()
    fixed = set(fixed_idx)
    n = len(syms)
    # Phase 3b: class-conditional gate for auto-enabled Phase 3a flags.
    # Hapto systems (Cp/allyl/π) skip COD_BONDS + RIGID_H_DRAG to avoid refine-loop
    # oscillation that caused 3/4 Sn multi-hapto smoke timeouts under Phase 3a.
    _hapto = _PURE_TRACK3 and _has_hapto_signature(syms, P)
    _use_rigid_h = _RIGID_H_DRAG and not _hapto
    _use_cod = _USE_COD_BONDS and not _hapto
    _, bonded, adj = _bonds_adj(syms, P)
    best_loss, _ = _violations(syms, P, bonded, adj, use_cod=_use_cod)
    if best_loss == 0:
        return P
    for _ in range(max_passes):
        _, bonded, adj = _bonds_adj(syms, P)
        loss, moves = _violations(syms, P, bonded, adj, use_cod=_use_cod)
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
        # Iter 30: rigid-H drag — propagate heavy-atom displacement to bonded H atoms
        # so X-H bonds preserve length instead of stretching as the heavy drifts.
        # H atoms must not already have an own displacement (else they get double-pushed
        # — overwrite only when cnt[h]==0).
        if _use_rigid_h:
            # bonded H map: for each H, identify its closest heavy (parent) using the
            # current bonded graph that produced this pass (`adj`)
            for hi in range(n):
                if syms[hi] != "H" or hi in fixed or cnt[hi] > 0:
                    continue
                parents = [k for k in adj[hi] if syms[k] != "H" and not bd._is_metal(syms[k])]
                if not parents:
                    continue
                par = parents[0] if len(parents) == 1 else min(
                    parents, key=lambda k: float(np.linalg.norm(P[hi] - P[k])))
                if cnt[par] > 0:
                    disp[hi] = disp[par]      # rigid drag: H rides with its heavy parent
        trial = P + damp * disp
        _, tb, ta = _bonds_adj(syms, trial)
        tloss, _ = _violations(syms, trial, tb, ta, use_cod=_use_cod)
        # Phase G10 (User 2026-05-31): under PURE_TRACK3, accept moves even if
        # loss equals current — allows escape from saddle points where collapsed
        # bonds can't unwind because stretching one creates a temporary defect
        # elsewhere. Fundamental: gradient descent shouldn't reject equal-loss
        # moves on a complex landscape. Universal across all molecule classes.
        if _PURE_TRACK3:
            if tloss <= best_loss:
                P = trial; best_loss = tloss
            else:
                damp *= 0.6
                if damp < 0.05:
                    break
        elif tloss < best_loss:          # accept-if-better gate (legacy)
            P = trial; best_loss = tloss
        else:
            damp *= 0.6                # smaller step
            if damp < 0.05:
                break
    return P
