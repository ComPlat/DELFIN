"""Iter-25 (2026-05-20) — final bond-decollapse corrector.

Calibration breakthrough (project_calibration_breakthrough_2026_05_20): ~79% of
hapto structures emit COLLAPSED ligands — heavy-heavy pairs at 0.24-1.2 Å
(atoms overlapping / superimposed phenyls / substituents fused onto their
parent C).  The collapse originates in several builder paths (scaffold
substituent placement, BFS ring building, post-UFF on unparametrized metals);
rather than chase every origin, this is a FINAL guaranteed pass that runs on
the emitted XYZ after all builders/UFF/bandages and pushes any collapsed
non-metal geometry back to physical bond lengths.

Mechanism (per frame): a short spring relaxation on the heavy-atom graph —
  * bonded heavy-heavy pairs are pulled toward their element-pair ideal length
  * non-bonded heavy pairs closer than a floor are pushed apart (resolves
    superimposed phenyl rings)
  * metals AND every atom within M-D bonding distance are FROZEN, so the
    coordination sphere / η-ring / M-D invariant is never disturbed
H atoms are dragged rigidly with their heavy parent.

Per-frame rollback: only keep the relaxed frame if it strictly REDUCES the
collapsed-bond count and introduces no new metal-donor break.  Bit-exact when
no collapsed bond is present (early return).  XRD-validated target: drop the
collapse rate from ~79%.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

_METALS = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg La Ce Lu Sn Pb Ge Sb Bi".split())
_COV = {'C':0.76,'N':0.71,'O':0.66,'H':0.31,'S':1.05,'Cl':1.02,'P':1.07,'F':0.57,
        'Br':1.20,'I':1.39,'B':0.84,'Si':1.11,'Se':1.20,'As':1.19,'Te':1.38}
_MD_TOL = 0.05          # M-D invariant tolerance (Å)
_MD_FACTOR = 1.30       # M-D bond detection factor
_NONBOND_FLOOR = 0.78   # non-bonded heavy pair < this·Σcov → push apart


def _ideal_bond(a: str, b: str) -> float:
    return _COV.get(a, 0.9) + _COV.get(b, 0.9)


def _is_metal(s: str) -> bool:
    return s in _METALS


def _parse(xyz: str):
    lines = xyz.splitlines()
    syms: List[str] = []
    pts: List[List[float]] = []
    keep: List[str] = []
    for ln in lines:
        p = ln.split()
        if len(p) == 4:
            try:
                xyz_v = [float(p[1]), float(p[2]), float(p[3])]
            except ValueError:
                keep.append(ln)
                continue
            syms.append(p[0]); pts.append(xyz_v)
        else:
            keep.append(ln)
    return syms, np.array(pts, dtype=float), lines


def _count_collapsed(syms, P, bonds) -> int:
    n = 0
    for i, j in bonds:
        if _is_metal(syms[i]) or _is_metal(syms[j]):
            continue
        if float(np.linalg.norm(P[i] - P[j])) < 0.82 * _ideal_bond(syms[i], syms[j]):
            n += 1
    return n


def _md_snapshot(syms, P) -> Dict[Tuple[int, int], float]:
    snap = {}
    n = len(syms)
    for i in range(n):
        if not _is_metal(syms[i]):
            continue
        for j in range(n):
            if i == j or _is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= _MD_FACTOR * _ideal_bond(syms[i], syms[j]):
                snap[(i, j)] = d
    return snap


def _md_broken(snap, P) -> bool:
    for (i, j), d0 in snap.items():
        if abs(float(np.linalg.norm(P[i] - P[j])) - d0) > _MD_TOL:
            return True
    return False


def _geometric_bonds(syms, P) -> List[Tuple[int, int]]:
    """Heavy-heavy + X-H pairs within 1.3·Σcov are treated as bonds.  Atom-order
    independent (the emitted XYZ may not match the mol's atom indexing), and
    collapse-robust (a 0.24 Å pair is still < 1.3·Σcov so it is captured)."""
    n = len(syms)
    bonds: List[Tuple[int, int]] = []
    for i in range(n):
        if _is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _ideal_bond(syms[i], syms[j]):
                bonds.append((i, j))
    return bonds


def correct_xyz(mol, xyz: str) -> str:
    """Return a decollapsed copy of ``xyz`` (or the original if no gain).

    ``mol`` is accepted for signature parity but NOT used — connectivity is
    derived geometrically so the corrector is robust to mol/XYZ atom-order
    drift (the actual failure mode that made a mol-bond version detect 0
    collapses on a clearly-collapsed OKOKUU frame).
    """
    try:
        syms, P, _ = _parse(xyz)
    except Exception:
        return xyz
    n = len(syms)
    if n < 3:
        return xyz

    bonds = _geometric_bonds(syms, P)
    n_collapsed = _count_collapsed(syms, P, bonds)
    if n_collapsed == 0:
        return xyz  # bit-exact early return

    # Freeze metals + coordination sphere (atoms within M-D distance).
    frozen = set(i for i in range(n) if _is_metal(syms[i]))
    md = _md_snapshot(syms, P)
    for (mi, dj) in md:
        frozen.add(dj)
    md_snap = dict(md)

    # H → heavy parent (for rigid drag)
    h_parent: Dict[int, int] = {}
    heavy_nbr: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i, j in bonds:
        if syms[i] == 'H' and syms[j] != 'H':
            h_parent[i] = j
        elif syms[j] == 'H' and syms[i] != 'H':
            h_parent[j] = i
        if syms[i] != 'H' and syms[j] != 'H':
            heavy_nbr[i].append(j); heavy_nbr[j].append(i)

    work = P.copy()
    bonded_set = set()
    for i, j in bonds:
        bonded_set.add((min(i, j), max(i, j)))

    for _pass in range(250):
        forces = np.zeros_like(work)
        moved = False
        # bonded springs → ideal length
        for i, j in bonds:
            if syms[i] == 'H' or syms[j] == 'H':
                continue
            diff = work[j] - work[i]
            d = float(np.linalg.norm(diff))
            tgt = _ideal_bond(syms[i], syms[j])
            if d < 1e-6:
                diff = np.random.default_rng(_pass).standard_normal(3); d = float(np.linalg.norm(diff))
            if abs(d - tgt) > 0.02:
                f = 0.5 * (d - tgt) * (diff / d)
                forces[i] += f; forces[j] -= f; moved = True
        # non-bonded repulsion (resolve superimposed atoms)
        for a in range(n):
            if syms[a] == 'H':
                continue
            for b in range(a + 1, n):
                if syms[b] == 'H':
                    continue
                if (a, b) in bonded_set:
                    continue
                diff = work[b] - work[a]
                d = float(np.linalg.norm(diff))
                floor = _NONBOND_FLOOR * _ideal_bond(syms[a], syms[b])
                if d < floor and d > 1e-6:
                    f = 0.5 * (d - floor) * (diff / d)
                    forces[a] += f; forces[b] -= f; moved = True
        if not moved:
            break
        # apply, but never move frozen atoms
        for i in range(n):
            if i in frozen or syms[i] == 'H':
                continue
            work[i] = work[i] + forces[i]
    # drag H rigidly to follow their parent's displacement
    for h, par in h_parent.items():
        work[h] = work[h] + (work[par] - P[par])

    # acceptance: collapse strictly reduced + M-D not broken
    if _count_collapsed(syms, work, bonds) >= n_collapsed:
        return xyz
    if _md_broken(md_snap, work):
        return xyz

    # rebuild xyz preserving original non-atom lines
    out_lines = xyz.splitlines()
    ai = 0
    res = []
    for ln in out_lines:
        p = ln.split()
        if len(p) == 4:
            try:
                float(p[1])
            except ValueError:
                res.append(ln); continue
            x, y, z = work[ai]
            res.append(f"{p[0]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            ai += 1
        else:
            res.append(ln)
    return "\n".join(res) + ("\n" if xyz.endswith("\n") else "")


def correct_results(mol, results):
    """Apply correct_xyz to every (xyz, label) frame in results."""
    if not results:
        return results
    out = []
    for item in results:
        try:
            xyz, label = item
        except Exception:
            out.append(item); continue
        try:
            out.append((correct_xyz(mol, xyz), label))
        except Exception:
            out.append((xyz, label))
    return out
