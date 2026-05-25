"""delfin.fffree.ligand_relax — ligand-internal relaxer (#38, v3: torsional/rigid-body).

THE GRUNDPFEILER (division of labor): fffree fixes the metal coordination sphere BY
CONSTRUCTION (enumeration + ideal polyhedra, metal-FF-free) and relaxes ONLY the ligand
internals, metal + donor atoms FROZEN.  This is the inverse of UFF.

v1 (force, no firewall) and v2 (force + every guard) both FAILED: a force-based relaxer
with independently-movable atoms cannot preserve a rigid multi-atom group's ANGLES
(nitrate O-N-O bent 18 deg despite bond-restraint + perpendicular projection + hard 5%
bond guard).  v3 is the architecturally-correct universal mechanism: RIGID-BODY /
TORSIONAL.

  - Rigid clusters are detected GEOMETRICALLY (element-agnostic, no element lists): atoms
    joined by non-rotatable bonds (ring bonds + short = multiple/conjugated bonds) form a
    rigid body.  Single, non-ring bonds are the only rotatable degrees of freedom.
  - Relaxation rotates the distal sub-tree about each rotatable bond (Rodrigues), which
    preserves EVERY bond length AND EVERY bond angle within rigid units BY CONSTRUCTION —
    nitrate, carboxylate, amide, aromatic ring stay exact for ANY element.  Only torsions
    change.  H rides with its rigid fragment automatically.
  - A rotation whose distal sub-tree contains a frozen atom (metal/donor) is skipped, so
    the constructed coordination is never disturbed (chelate second donors stay put).
  - Deterministic coordinate (torsion) descent: fixed bond order, fixed angle scan, greedy
    accept-best; a candidate is accepted only if the clash loss strictly drops AND the
    multi-axis firewall (clash/distort/planarity/H-anomaly/bad-angle) does not rise above
    the input baseline AND ‖M-D‖ holds.  Never-worse self-gate.

Targets the CONFORMATIONAL axes (vdw_clashes, lig_realistic, hanom).  It does NOT change
intra-group angles/planarity (F23/F24) -- those are placement-level and need a separate
construction lever.  Geometry-only, universal, deterministic.  Disabled by default; enable
via DELFIN_FFFREE_LIGAND_RELAX=1.
"""
from __future__ import annotations
from typing import Dict, List, Set, Tuple
import numpy as np
import delfin._bond_decollapse as bd

_VDW = bd._VDW
_VDW_D = bd._VDW_DEFAULT
_CLASH_FACTOR = 0.85
_MULTIPLE_BOND = 0.92          # bond shorter than this*ideal = multiple/conjugated -> rigid
_ANGLE_SCAN = (60.0, 120.0, 180.0, 240.0, 300.0)   # staggered torsion states
_MAX_SWEEPS = 2
_MAX_ROTATABLE = 12         # cap per-structure cost (huge flexible ligands -> skip extras)


def _vdw(s: str) -> float:
    return _VDW.get(s, _VDW_D)


def _bonds_adj(syms, P):
    b = bd._geometric_bonds(syms, P)
    bonded = {(min(i, j), max(i, j)) for i, j in b}
    adj: List[List[int]] = [[] for _ in syms]
    for i, j in b:
        adj[i].append(j); adj[j].append(i)
    return b, bonded, adj


def _ring_bonds(bonded, adj) -> Set[Tuple[int, int]]:
    """Bonds that lie in a cycle (removing the bond leaves its endpoints connected).
    Universal, geometry-only.  Ligands are small -> cheap."""
    ring: Set[Tuple[int, int]] = set()
    for (a, b) in bonded:
        seen = {a}; stack = [a]; reach = False
        while stack and not reach:
            x = stack.pop()
            for y in adj[x]:
                if (min(x, y), max(x, y)) == (a, b):
                    continue                       # the removed bond
                if y == b:
                    reach = True; break
                if y not in seen:
                    seen.add(y); stack.append(y)
        if reach:
            ring.add((a, b))
    return ring


def _bfs_dist(roots, adj, n) -> List[int]:
    """Hop distance from the nearest root atom (metal/donor) to every atom."""
    INF = 10 ** 9
    dist = [INF] * n
    from collections import deque
    q = deque()
    for r in roots:
        if 0 <= r < n:
            dist[r] = 0; q.append(r)
    while q:
        x = q.popleft()
        for y in adj[x]:
            if dist[y] == INF:
                dist[y] = dist[x] + 1; q.append(y)
    return dist


def _subtree(child, parent, adj) -> Set[int]:
    """Atoms reachable from `child` without crossing back through `parent` — the group
    that rotates about the (parent,child) bond."""
    seen = {parent, child}; stack = [child]; tree = {child}
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); tree.add(y); stack.append(y)
    return tree


def _rotate(P, idxs, pt, axis, theta_deg):
    """Rodrigues rotation of points `idxs` about the line (pt, axis) by theta_deg."""
    k = axis / (np.linalg.norm(axis) + 1e-12)
    th = np.radians(theta_deg); c, s = np.cos(th), np.sin(th)
    Q = P.copy()
    for i in idxs:
        v = P[i] - pt
        Q[i] = pt + v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1 - c)
    return Q


def _clash_loss(syms, P, bonded) -> float:
    """Sum of squared van-der-Waals overlaps of non-bonded pairs (heavy + H).  The
    conformational quantity torsions can relieve without touching local geometry."""
    n = len(syms); loss = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= 1e-6:
                continue
            tgt = _CLASH_FACTOR * (_vdw(syms[i]) + _vdw(syms[j]))
            if d < tgt:
                loss += (tgt - d) ** 2
    return loss


def _count_h_bad(syms, P, bonds) -> int:
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


def _metal_shell(syms, P) -> int:
    """Occupancy of every metal's first coordination shell (heavy atoms within bonding
    range).  A torsion must NEVER swing a backbone atom INTO the coordination sphere —
    that is what perturbs the constructed polyhedron (coord_geom).  Coordination is
    fffree's sacred strength; this guard keeps it inviolable.  Geometry-only, universal."""
    n = len(syms); shell = 0
    for m in range(n):
        if not bd._is_metal(syms[m]):
            continue
        for j in range(n):
            if j == m or syms[j] == "H":
                continue
            # mirror the coord_geom detector's M-D detection EXACTLY (min(1.65*ideal,
            # 2.90)) so any atom swinging into the polyhedron-detection sphere is caught
            # -> coord_geom cannot regress (build-against-the-metric).
            cut = min(1.65 * bd._ideal_bond(syms[m], syms[j]), 2.90)
            if float(np.linalg.norm(P[j] - P[m])) < cut:
                shell += 1
    return shell


def _firewall(syms, P):
    """Never-worse multi-axis guard mirroring the metrics earlier versions broke."""
    bonds, bonded, adj = _bonds_adj(syms, P)
    excl = bd._neighbor_exclusion(len(syms), bonds)
    return {
        "distort": bd._count_bond_distort(syms, P, bonds),
        "vdw": bd._count_vdw_clashes(syms, P, excl),
        "hplanar": bd._count_h_planar_viol(syms, P),
        "hbad": _count_h_bad(syms, P, bonds),
        "badang": bd._count_bad_angles(syms, P, adj),
        "mshell": _metal_shell(syms, P),    # coordination-sphere intrusion guard (v4)
    }


def relax(syms, P, fixed_idx: Set[int], max_sweeps: int = _MAX_SWEEPS):
    """Torsional/rigid-body relaxer.  Rotates distal sub-trees about single, non-ring,
    non-multiple bonds to minimise van-der-Waals clashes, preserving all bond lengths and
    angles within rigid fragments.  Metal + donors frozen; sub-trees containing a frozen
    atom are not rotated.  Multi-axis never-worse firewall + ‖M-D‖ gate.  Deterministic."""
    P = np.array(P, float).copy()
    if P.size == 0 or not np.all(np.isfinite(P)):
        return P
    n = len(syms); syms = list(syms)
    bonds, bonded, adj = _bonds_adj(syms, P)
    frozen = set(int(i) for i in fixed_idx)
    metals = [i for i in range(n) if bd._is_metal(syms[i])]
    # FREEZE every atom coord_geom counts as a coordinating donor (its cutoff
    # min(1.65*ideal, 2.90)).  The caller's donor set / a geometric fixed set can be
    # TIGHTER than this, leaving a borderline donor (1.45x-1.65x) movable -> rotating it
    # perturbs the constructed polyhedron (coord_geom).  Matching the detector's donor
    # radius makes the coordination inviolable in BOTH smoke and production. Universal.
    for m in metals:
        for j in range(n):
            if j == m or syms[j] == "H":
                continue
            if float(np.linalg.norm(P[j] - P[m])) < min(1.65 * bd._ideal_bond(syms[m], syms[j]), 2.90):
                frozen.add(j)

    ringb = _ring_bonds(bonded, adj)
    roots = metals if metals else ([min(frozen)] if frozen else [0])
    dist = _bfs_dist(roots, adj, n)

    # rotatable bonds: single, non-ring, non-multiple, heavy-heavy, non-metal.  Oriented
    # parent(closer to root)->child(farther); the child sub-tree rotates.
    rotatable: List[Tuple[int, int, frozenset]] = []
    for (a, b) in sorted(bonded):
        if (a, b) in ringb or syms[a] == "H" or syms[b] == "H":
            continue
        if bd._is_metal(syms[a]) or bd._is_metal(syms[b]):
            continue
        d = float(np.linalg.norm(P[a] - P[b]))
        ideal = bd._ideal_bond(syms[a], syms[b])
        if ideal > 1e-6 and d < _MULTIPLE_BOND * ideal:
            continue                                   # multiple/conjugated -> rigid
        parent, child = (a, b) if dist[a] <= dist[b] else (b, a)
        sub = _subtree(child, parent, adj)
        if sub & frozen or parent in sub:              # would move a frozen donor -> skip
            continue
        if len(sub) < 2:                               # nothing meaningful to rotate
            continue
        rotatable.append((parent, child, frozenset(sub)))
    if not rotatable:
        return P
    if len(rotatable) > _MAX_ROTATABLE:                # cap cost: largest sub-trees first
        rotatable.sort(key=lambda t: (-len(t[2]), t[0], t[1]))
        rotatable = rotatable[:_MAX_ROTATABLE]

    fw0 = _firewall(syms, P)
    md_snap = bd._md_snapshot(syms, P)
    best_clash = _clash_loss(syms, P, bonded)
    if best_clash <= 1e-9:
        return P

    # Speed: a rigid rotation PRESERVES bonds, so `bonded` is reused for clash scoring
    # (no per-angle bond rebuild).  The angle scan uses only the cheap clash loss; the
    # expensive multi-axis firewall + M-D check run ONCE per bond on the BEST candidate.
    for _ in range(max_sweeps):
        improved = False
        for parent, child, sub in rotatable:
            axis = P[child] - P[parent]
            if float(np.linalg.norm(axis)) < 1e-6:
                continue
            idxs = list(sub)
            best_theta, best_l = None, best_clash
            for theta in _ANGLE_SCAN:
                trial = _rotate(P, idxs, P[parent], axis, theta)
                cl = _clash_loss(syms, trial, bonded)     # bonds preserved -> reuse
                if cl < best_l - 1e-9:
                    best_theta, best_l = theta, cl
            if best_theta is None:
                continue
            trial = _rotate(P, idxs, P[parent], axis, best_theta)
            if not np.all(np.isfinite(trial)) or bd._md_broken(md_snap, trial):
                continue
            fw = _firewall(syms, trial)                   # ONCE on best candidate only
            if any(fw[k] > fw0[k] for k in fw0):
                continue
            P = trial; best_clash = best_l; improved = True
        if not improved:
            break
    return P


if __name__ == "__main__":
    # self-test: a clashing biphenyl-like pair across a rotatable C-C, plus a rigid
    # nitrate.  relax must relieve the clash by torsion while keeping nitrate EXACT and
    # all bonds/angles within fragments unchanged.
    syms = ["Fe", "O", "C", "C", "C", "O", "N", "O", "O", "O"]
    P = np.array([
        [0.0, 0.0, 0.0],            # Fe metal
        [1.9, 0.0, 0.0],            # O donor (frozen)
        [3.0, 0.3, 0.0],            # C (bonded to O donor)
        [4.2, -0.3, 0.0],           # C  (rotatable C2-C3)
        [5.3, 0.4, 0.1],            # C
        [3.1, 1.6, 0.2],            # O on C2 (will clash with C4 region)
        [0.0, 4.0, 0.0],            # nitrate N
        [1.2, 4.0, 0.0], [-0.6, 5.0, 0.0], [-0.6, 3.0, 0.0],
    ], float)
    fixed = {0, 1}
    def no3(X):
        N = X[6]; Os = X[7:10]
        bl = [float(np.linalg.norm(o - N)) for o in Os]
        ang = []
        for a in range(3):
            for b in range(a + 1, 3):
                va, vb = Os[a] - N, Os[b] - N
                ang.append(float(np.degrees(np.arccos(np.clip(
                    np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)), -1, 1)))))
        return np.array(bl), np.array(ang)
    bl0, ang0 = no3(P)
    _, bo0, _ = _bonds_adj(syms, P); cl0 = _clash_loss(syms, P, bo0)
    o1 = relax(syms, P, fixed); o2 = relax(syms, P, fixed)
    _, bo1, _ = _bonds_adj(syms, o1); cl1 = _clash_loss(syms, o1, bo1)
    bl1, ang1 = no3(o1)
    fw0 = _firewall(syms, P); fw1 = _firewall(syms, o1)
    print(f"clash loss {cl0:.4f} -> {cl1:.4f}")
    print(f"firewall in={fw0} out={fw1} never-worse={all(fw1[k]<=fw0[k] for k in fw0)}")
    print(f"metal+donor frozen: {np.allclose(o1[0],P[0]) and np.allclose(o1[1],P[1])}")
    print(f"nitrate EXACT: bonds dmax={np.abs(bl1-bl0).max():.5f}A "
          f"angles dmax={np.abs(ang1-ang0).max():.3f}deg")
    print(f"deterministic: {np.allclose(o1,o2)}  finite: {np.all(np.isfinite(o1))}")
    print(f"M-D intact: {not bd._md_broken(bd._md_snapshot(syms,P), o1)}")
