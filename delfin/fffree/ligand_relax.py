"""delfin.fffree.ligand_relax — COD-empirical ligand-internal relaxer (#38, v2).

THE GRUNDPFEILER (division of labor): fffree fixes the metal coordination sphere
BY CONSTRUCTION (enumeration + ideal polyhedra, metal-FF-free) and relaxes ONLY the
ligand internals with a COD-empirical organic potential, metal + donor atoms FROZEN.
This is the inverse of UFF (which relaxes everything with a poor metal force field).

v1 failed its smoke gate HARD (net -6, 14 severe: nitrate 0->100% broken, bondlen
outliers 0.36->15.3%, methyl broken, hanom +409%, coord_geom +100% despite "frozen").
v2 fixes the five root causes (see feedback_relax38_v1_failed_smoke):

  (1) MULTI-AXIS FIREWALL accept-gate (like the #32 decollapse): a pass is accepted
      ONLY if the continuous loss decreases AND no firewall axis (bond-distort, vdW
      clash, sp2 planarity, H-anomaly) rises above the INPUT baseline AND ‖M-D‖ holds.
      The scalar loss alone is insufficient — it can fall while a group breaks locally.
  (2) H RIGID WITH PARENT — H atoms are NOT free DOF; each H is translated rigidly with
      its parent heavy atom (preserves C-H length + H-C-H angles -> methyl intact).
  (3) RIGID GROUPS frozen — oxoanions / carboxylates (a C/N/S/P centre with >=2 bonded
      O, or >=3 for S/P) are held rigid (frozen) so the relaxer cannot distort them.
  (4) FROZEN SET expanded to donor + donor's heavy neighbours -> the constructed
      coordination polyhedron is inviolable (fixes the coord_geom leak).
  (5) Strict NEVER-WORSE gate -> the relaxed build is never worse than its input on any
      firewall axis, so wiring it can only help or no-op (self-gate philosophy #30).

Continuous loss on the remaining free heavy atoms: bond->r0 (COD covalent ideals),
Urey-Bradley angle->theta0(steric number), conservative sp2 residual de-puckering
(does NOT flatten genuine sp3), clash penalty.  Deterministic coordinate descent
(no RNG, fixed damping schedule, per-step clamp, finite guard).  Geometry-only +
graph-free -> universal, atom-order-independent.  Disabled by default; enable via
DELFIN_FFFREE_LIGAND_RELAX=1.
"""
from __future__ import annotations
from typing import Dict, List, Set, Tuple
import numpy as np
import delfin._bond_decollapse as bd

_VDW = bd._VDW
_VDW_D = bd._VDW_DEFAULT
_CLASH_FACTOR = 0.85

_W_BOND = 1.0
_W_ANG = 0.5
_W_PLAN = 1.0
_W_CLASH = 2.0

_THETA0 = {2: 180.0, 3: 120.0, 4: 109.47}
_MAX_STEP = 0.25
_PLANAR_OOP_TOL = 0.50


def _vdw(s: str) -> float:
    return _VDW.get(s, _VDW_D)


def _bonds_adj(syms, P):
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    adj: List[List[int]] = [[] for _ in syms]
    for i, j in b:
        adj[i].append(j); adj[j].append(i)
    return b, bonded, adj


def _h_parent_map(syms, P, adj) -> Dict[int, int]:
    """Each H -> its nearest bonded heavy atom (parent).  H is dragged rigidly with
    the parent (preserves local H geometry), never optimised as a free DOF."""
    hp: Dict[int, int] = {}
    for h in range(len(syms)):
        if syms[h] != "H":
            continue
        heavies = [k for k in adj[h] if syms[k] != "H"]
        if not heavies:
            heavies = [k for k in range(len(syms)) if k != h and syms[k] != "H"]
        if heavies:
            hp[h] = min(heavies, key=lambda k: float(np.linalg.norm(P[h] - P[k])))
    return hp


def _rigid_group_atoms(syms, P, adj) -> Set[int]:
    """Atoms of rigid functional groups that the per-atom loss would otherwise break:
    oxoanions / carboxylates -- a C/N centre with >=2 bonded O, or S/P with >=3 bonded
    O.  Freeze the centre + its O (and the O's H, e.g. -COOH) so they move only via the
    frozen coordination, never internally.  Geometry-only."""
    rigid: Set[int] = set()
    for c in range(len(syms)):
        sc = syms[c]
        if sc not in ("C", "N", "S", "P"):
            continue
        oxy = [k for k in adj[c] if syms[k] == "O"]
        need = 3 if sc in ("S", "P") else 2
        if len(oxy) >= need:
            rigid.add(c)
            for o in oxy:
                rigid.add(o)
                for k in adj[o]:           # -OH hydrogens / terminal
                    if syms[k] == "H":
                        rigid.add(k)
    return rigid


def _loss_and_moves(syms, P, bonded, adj, movable):
    """Continuous COD loss (scalar) + per-atom displacement moves for MOVABLE heavy
    atoms only.  Pure function, deterministic."""
    n = len(syms)
    loss = 0.0
    moves: List[Tuple[int, np.ndarray]] = []

    for i, j in bonded:
        if syms[i] == "H" or syms[j] == "H":
            continue
        d = float(np.linalg.norm(P[i] - P[j]))
        if d <= 1e-6:
            continue
        r0 = bd._ideal_bond(syms[i], syms[j])
        dev = d - r0
        loss += _W_BOND * dev * dev
        u = (P[j] - P[i]) / d
        corr = 0.5 * _W_BOND * (-dev)
        moves.append((j, +corr * u)); moves.append((i, -corr * u))

    for j in range(n):
        if bd._is_metal(syms[j]) or syms[j] == "H":
            continue
        nbrs = [k for k in adj[j] if syms[k] != "H"]
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
        oop = float(np.dot(P[j] - a, nrm))
        if abs(oop) > _PLANAR_OOP_TOL:
            continue
        loss += _W_PLAN * oop * oop
        moves.append((j, 0.5 * _W_PLAN * (-oop) * nrm))

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

    # keep only moves on movable heavy atoms
    moves = [(idx, d) for idx, d in moves if idx in movable]
    return loss, moves


def _count_h_bad(syms, P, bonds) -> int:
    """H-anomaly proxy: H whose bond to its parent is outside [0.7,1.4]*ideal, or H
    bonded to >=2 heavy atoms (over-coordination).  Mirrors the hanom explosion v1 hit."""
    n = len(syms)
    adjH: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adjH[i].append(j); adjH[j].append(i)
    cnt = 0
    for h in range(n):
        if syms[h] != "H":
            continue
        heavies = [k for k in adjH[h] if syms[k] != "H" and not bd._is_metal(syms[k])]
        if len(heavies) >= 2:
            cnt += 1; continue
        for k in adjH[h]:
            if syms[k] == "H":
                continue
            d = float(np.linalg.norm(P[h] - P[k]))
            ideal = bd._ideal_bond(syms[h], syms[k])
            if ideal > 1e-6 and (d < 0.7 * ideal or d > 1.4 * ideal):
                cnt += 1; break
    return cnt


def _firewall(syms, P):
    """Cheap geometry-only proxies mirroring the metrics v1 broke.  A pass that raises
    ANY of these above the input baseline is rejected (never-worse guarantee)."""
    bonds, bonded, _ = _bonds_adj(syms, P)
    excl = bd._neighbor_exclusion(len(syms), bonds)
    return {
        "distort": bd._count_bond_distort(syms, P, bonds),
        "vdw": bd._count_vdw_clashes(syms, P, excl),
        "hplanar": bd._count_h_planar_viol(syms, P),
        "hbad": _count_h_bad(syms, P, bonds),
    }


def relax(syms, P, fixed_idx: Set[int], max_passes: int = 120, damp: float = 0.4):
    """COD-loss coordinate descent on free ligand-internal heavy atoms.  Frozen:
    metal + donors + donor heavy-neighbours + rigid oxoanion/carboxylate groups.  H is
    dragged rigidly with its parent.  Multi-axis firewall + ‖M-D‖ never-worse gate.
    Deterministic.  Returns relaxed P, or the original if nothing strictly improved."""
    P = np.array(P, float).copy()
    if P.size == 0 or not np.all(np.isfinite(P)):
        return P
    n = len(syms)
    syms = list(syms)
    bonds, bonded, adj = _bonds_adj(syms, P)

    # frozen set: metal + donors (caller) + donor heavy-neighbours + rigid groups
    frozen = set(int(i) for i in fixed_idx)
    for d_ in list(frozen):
        if d_ < n:
            for k in adj[d_]:
                if syms[k] != "H":
                    frozen.add(k)
    frozen |= _rigid_group_atoms(syms, P, adj)

    hp = _h_parent_map(syms, P, adj)
    h_offset = {h: P[h] - P[par] for h, par in hp.items()}

    # movable = heavy, non-metal, not frozen, not an H
    movable = set(i for i in range(n)
                  if syms[i] != "H" and not bd._is_metal(syms[i]) and i not in frozen)
    if not movable:
        return P

    fw0 = _firewall(syms, P)
    md_snap = bd._md_snapshot(syms, P)
    best_loss, _ = _loss_and_moves(syms, P, bonded, adj, movable)
    if best_loss <= 1e-9:
        return P

    for _ in range(max_passes):
        bonds, bonded, adj = _bonds_adj(syms, P)
        loss, moves = _loss_and_moves(syms, P, bonded, adj, movable)
        disp = np.zeros((n, 3)); cnt = np.zeros(n)
        for idx, dvec in moves:
            disp[idx] += dvec; cnt[idx] += 1
        for i in range(n):
            if cnt[i] > 0:
                disp[i] /= cnt[i]
        step = damp * disp
        norms = np.linalg.norm(step, axis=1)
        big = norms > _MAX_STEP
        if np.any(big):
            step[big] *= (_MAX_STEP / norms[big])[:, None]
        trial = P + step
        # drag H rigidly with parent
        for h, par in hp.items():
            trial[h] = trial[par] + h_offset[h]
        if not np.all(np.isfinite(trial)):
            break
        # never-worse multi-axis firewall + M-D
        if bd._md_broken(md_snap, trial):
            damp *= 0.6
            if damp < 0.03:
                break
            continue
        tb, tbonded, ta = _bonds_adj(syms, trial)
        tloss, _ = _loss_and_moves(syms, trial, tbonded, ta, movable)
        fw = _firewall(syms, trial)
        worse = any(fw[k] > fw0[k] for k in fw0)
        if tloss < best_loss - 1e-9 and not worse:
            P = trial; best_loss = tloss
        else:
            damp *= 0.6
            if damp < 0.03:
                break
    return P


if __name__ == "__main__":
    # smoke self-test: puckered amide-C + stretched bond; rigid nitrate must stay intact.
    syms = ["Fe", "N", "C", "O", "H", "N", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.6, 0.4], [3.2, 1.6, 0.0],
        [2.6, -0.7, 0.0],
        [0.0, 3.0, 0.0], [1.2, 3.0, 0.0], [-0.6, 4.0, 0.0], [-0.6, 2.0, 0.0],   # nitrate NO3
    ], float)
    fixed = {0, 1}
    no3_before = P[5:9].copy()
    out1 = relax(syms, P, fixed); out2 = relax(syms, P, fixed)
    b0, bo0, a0 = _bonds_adj(syms, P); l0, _ = _loss_and_moves(syms, P, bo0, a0,
        set(range(len(syms))))
    b1, bo1, a1 = _bonds_adj(syms, out1); l1, _ = _loss_and_moves(syms, out1, bo1, a1,
        set(range(len(syms))))
    fw0 = _firewall(syms, P); fw1 = _firewall(syms, out1)
    print(f"loss {l0:.4f} -> {l1:.4f}")
    print(f"firewall in={fw0} out={fw1}  never-worse={all(fw1[k]<=fw0[k] for k in fw0)}")
    print(f"metal+donor frozen: {np.allclose(out1[0],P[0]) and np.allclose(out1[1],P[1])}")
    print(f"nitrate intact: {np.allclose(out1[5:9], no3_before)}")
    print(f"deterministic: {np.allclose(out1,out2)}  finite: {np.all(np.isfinite(out1))}")
    print(f"M-D intact: {not bd._md_broken(bd._md_snapshot(syms,P), out1)}")
