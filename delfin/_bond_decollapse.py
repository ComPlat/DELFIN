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

import os
from typing import Dict, List, Tuple

import numpy as np

_METALS = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg La Ce Lu Sn Pb Ge Sb Bi".split())
_COV = {'C':0.76,'N':0.71,'O':0.66,'H':0.31,'S':1.05,'Cl':1.02,'P':1.07,'F':0.57,
        'Br':1.20,'I':1.39,'B':0.84,'Si':1.11,'Se':1.20,'As':1.19,'Te':1.38}
_MD_TOL = 0.05          # M-D invariant tolerance (Å)
_MD_FACTOR = 1.30       # M-D bond detection factor
_NONBOND_FLOOR = 0.78   # non-bonded heavy pair < this·Σcov → push apart

# Iter-26 (Task #32): vdw/angle-awareness so the decollapse pass cannot trade a
# collapse fix for vdw-clash / H-planarity / angle regressions (the Iter-25
# concentrated trade vs iter23: vdw +20%, F3_angle +33%, F20 +3.6pp).
# Bondi vdW radii — verbatim from quality_framework/scripts/nonbonded_vdw_check.py
# so the in-corrector clash proxy matches the detector bit-for-bit on radii.
_VDW = {
    "H": 1.20, "He": 1.40, "Li": 1.81, "Be": 1.52, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.15, "Ti": 2.11, "V": 2.07, "Cr": 2.06,
    "Mn": 2.05, "Fe": 2.04, "Co": 2.00, "Ni": 1.97, "Cu": 1.96, "Zn": 2.01,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18, "Mo": 2.17,
    "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43, "Hf": 2.23, "Ta": 2.22, "W": 2.18,
    "Re": 2.16, "Os": 2.16, "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23,
    "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
}
_VDW_DEFAULT = 1.7      # detector default for missing element
_VDW_FACTOR = 0.85      # clash if d < factor·(vdW_a + vdW_b)  (matches detector)
_ANGLE_MIN_DEG = 70.0   # optional angle proxy: heavy-heavy-heavy angle < this = bad
_F20_OOP_TOL = 0.20     # ring-H out-of-plane tolerance (Å)   (matches detector)
_F20_RING_TOL = 0.10    # ring planarity tolerance (Å)        (matches detector)
_F20_AROMATIC = {"C", "N", "O", "S"}


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


def _neighbor_exclusion(n: int, bonds: List[Tuple[int, int]]) -> List[set]:
    """1-2 + 1-3 neighbor exclusion sets derived purely from the geometric
    bonds (no SMILES graph).  Mirrors the detector's 1-3 closure
    (nonbonded_vdw_check._smiles_to_atom_pairs): excl[i] = adj[i] ∪
    (⋃_{j∈adj[i]} adj[j]) − {i}.  Used by both the vdw proxy and the
    pair-aware repulsion floor so true steric clashes are treated
    identically to the detector while 1-3 contacts are left alone."""
    adj: List[set] = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    excl: List[set] = []
    for i in range(n):
        s = set(adj[i])
        for j in adj[i]:
            s.update(adj[j])
        s.discard(i)
        excl.append(s)
    return excl


def _count_vdw_clashes(syms, P, excl) -> int:
    """Geometry-only reproduction of nonbonded_vdw_check.vdw_clash_report:
    non-bonded (not 1-2/1-3) heavy pair is a clash if d < 0.85·(vdW_a+vdW_b).
    Heavy-only, metal-metal pairs skipped."""
    n = len(syms)
    cnt = 0
    for i in range(n):
        if syms[i] == 'H':
            continue
        mi = _is_metal(syms[i])
        ri = _VDW.get(syms[i], _VDW_DEFAULT)
        ex = excl[i]
        for j in range(i + 1, n):
            if syms[j] == 'H' or j in ex:
                continue
            if mi and _is_metal(syms[j]):
                continue
            min_d = _VDW_FACTOR * (ri + _VDW.get(syms[j], _VDW_DEFAULT))
            if float(np.linalg.norm(P[i] - P[j])) < min_d:
                cnt += 1
    return cnt


def _count_h_planar_viol(syms, P) -> int:
    """F20 proxy — port of compute_ligand_angles.sp2_h_planarity_check onto
    (syms, P).  Find 5/6-rings of aromatic-eligible heavy atoms (C/N/O/S),
    confirm ring planar (≤ _F20_RING_TOL), count ring-bonded H whose OOP
    distance from the ring plane exceeds _F20_OOP_TOL.  Self-contained
    distance-bond adjacency (includes H) so it matches the detector."""
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        ri = _COV.get(syms[i], 0.76)
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * (ri + _COV.get(syms[j], 0.76)):
                nbrs[i].append(j)
                nbrs[j].append(i)
    heavy_only: List[List[int]] = [
        [j for j in nbrs[i] if syms[j] != 'H' and syms[i] in _F20_AROMATIC
         and syms[j] in _F20_AROMATIC]
        for i in range(n)
    ]
    rings: set = set()
    for start in range(n):
        if syms[start] not in _F20_AROMATIC:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in heavy_only[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))
    cnt = 0
    for ring in rings:
        ring_pts = np.array([P[i] for i in ring])
        centroid = ring_pts.mean(axis=0)
        centered = ring_pts - centroid
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        normal = vh[-1] / max(float(np.linalg.norm(vh[-1])), 1e-9)
        if float(np.max(np.abs(centered @ normal))) > _F20_RING_TOL:
            continue
        for ri_idx in ring:
            for hi in nbrs[ri_idx]:
                if syms[hi] != 'H':
                    continue
                if float(abs(np.dot(P[hi] - centroid, normal))) > _F20_OOP_TOL:
                    cnt += 1
    return cnt


def _count_bad_angles(syms, P, heavy_nbr) -> int:
    """Optional angle proxy — count heavy-heavy-heavy bond angles below
    _ANGLE_MIN_DEG (physically implausible for any hybridization, exactly
    what collapse + repulsion-overshoot produce).  Monotone proxy for F3,
    not the exact metric (F3 needs RDKit hybridization)."""
    cnt = 0
    for c in range(len(syms)):
        if _is_metal(syms[c]) or syms[c] == 'H':
            continue
        nb = [k for k in heavy_nbr[c] if syms[k] != 'H']
        for a_i in range(len(nb)):
            va = P[nb[a_i]] - P[c]
            na = float(np.linalg.norm(va))
            if na < 1e-6:
                continue
            for b_i in range(a_i + 1, len(nb)):
                vb = P[nb[b_i]] - P[c]
                nbn = float(np.linalg.norm(vb))
                if nbn < 1e-6:
                    continue
                cosang = float(np.dot(va, vb) / (na * nbn))
                cosang = max(-1.0, min(1.0, cosang))
                if np.degrees(np.arccos(cosang)) < _ANGLE_MIN_DEG:
                    cnt += 1
    return cnt


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

    # Iter-26: 1-2/1-3 exclusion + pre-relaxation baselines for the vdw/F20/
    # angle firewall (reject any relaxed frame that worsens a measured axis).
    excl = _neighbor_exclusion(n, bonds)
    anglegate = os.environ.get("DELFIN_BOND_DECOLLAPSE_ANGLEGATE", "0") == "1"
    clashes0 = _count_vdw_clashes(syms, P, excl)
    h_planar0 = _count_h_planar_viol(syms, P)
    angle0 = _count_bad_angles(syms, P, heavy_nbr) if anglegate else 0

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
                # Iter-26: for TRUE non-bonded pairs (not 1-3 neighbors) lift the
                # floor to the vdw-clash threshold so the relaxation stops
                # settling atoms into the clash band; 1-3 contacts keep the
                # gentler covalent floor (else rings blow apart).
                if b not in excl[a] and not (_is_metal(syms[a]) and _is_metal(syms[b])):
                    floor = max(floor, _VDW_FACTOR * (_VDW.get(syms[a], _VDW_DEFAULT)
                                                      + _VDW.get(syms[b], _VDW_DEFAULT)))
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

    # acceptance (Iter-26 multi-axis firewall): keep the relaxed frame ONLY if
    # collapse strictly drops AND no measured axis worsens AND M-D intact.
    # Cheapest checks first so most rejections short-circuit before the SVD.
    if _count_collapsed(syms, work, bonds) >= n_collapsed:
        return xyz
    if _count_vdw_clashes(syms, work, excl) > clashes0:          # vdw (mandatory)
        return xyz
    if _count_h_planar_viol(syms, work) > h_planar0:            # F20 firewall
        return xyz
    if anglegate and _count_bad_angles(syms, work, heavy_nbr) > angle0:  # opt-in
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
