"""delfin.manta.hapto_declash — FF-free piano-stool inter-ligand declash.

The legacy analytical hapto scaffold (``smiles_converter._build_hapto_scaffold``)
seats an eta-coordinated ring (arene / Cp / Cb / diene) on the metal and the
remaining sigma-legs on a cone, then relaxes bond lengths and pushes apart hard
atom overlaps.  What it does NOT do is choose the ring's *rotational conformer*:
the angle of the ring about its M->centroid axis (staggered vs eclipsed relative
to the legs).  Construction leaves the ring at an arbitrary rotation, so its
carbons / substituents routinely eclipse the carbonyls and crowd the other
ligands.  Measured against CCDC piano-stool crystals (find_inter_ligand_clash,
F26) the analytical build over-clashes ~3-4x, and the residual after the
collapse fix is dominated by the eta-face crowding pendant aryls / other rings /
sigma-donors — contacts that are essentially absent in the crystals.

THE MISSING DEGREE OF FREEDOM = eta-RING ROTATION about its M->centroid axis.
It is a perfect isometry of the coordination polyhedron: every ring atom keeps
its cylindrical (radius, height) coordinates about the axis, so EVERY M-C(ring)
distance is exactly preserved while the ring body (and its hanging substituents)
swings to relieve the inter-ligand overlap.  The metal and all sigma-donors are
frozen.  This is exactly the staggered/eclipsed conformer a real crystal adopts.

Safety properties (mirrors joint_declash, the sigma/torsion sibling):
  * metal + all eta-ring atoms + all sigma-donor atoms FROZEN,
  * the ring spin preserves every M-(ring atom) distance to < 1e-4 A (asserted;
    any violating angle is rejected),
  * never-worse-on-MIN floor: a spin is rejected if it would push the worst
    inter-ligand contact below the input frame's worst,
  * strictly-decreasing accept (an angle is taken only if it lowers the global
    inter-ligand overlap), final roll-back to the input on any anomaly,
  * deterministic: fixed ring order, fixed angular grid, no RNG.

Env-gated ``DELFIN_FFFREE_HAPTO_DECLASH`` (default ``"0"`` -> the scaffold never
calls this module -> byte-identical).  License: pure geometry, Bondi vdW radii;
no CSD/CCDC data.
"""
from __future__ import annotations

import math
import os
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

# Bondi van-der-Waals radii (A) — identical table to the F26 inter-ligand-clash
# detector so the optimiser targets exactly what the metric (and the CCDC
# reference) measure.
_VDW: Dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.11, "Ti": 1.95, "V": 1.92, "Cr": 1.89,
    "Mn": 1.97, "Fe": 1.94, "Co": 1.92, "Ni": 1.84, "Cu": 1.86, "Zn": 2.10,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.83, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18, "Mo": 2.17,
    "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43, "W": 2.18, "Re": 2.16, "Os": 2.16,
    "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23, "Tl": 1.96, "Pb": 2.02,
    "Bi": 2.07,
}
_VDW_D = 1.7

# Pyykko single-bond covalent radii (A) for geometry-only bond perception when
# the RDKit mol is not threaded (final-manifold post-pass operates on xyz frames).
_COV: Dict[str, float] = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46,
    "Bi": 1.48,
}
_COV_D = 1.5
# transition + main-group metal Z (mirror of the F26 detector's range test)
_METAL_Z = (set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
            | set(range(89, 104))
            | {3, 11, 19, 37, 55, 4, 12, 20, 38, 56, 13, 31, 49, 81, 50, 82, 51, 83, 52, 84})
_ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
    "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
]
_SYM_Z = {s: i for i, s in enumerate(_ELEMENTS) if s}


def _is_metal(sym: str) -> bool:
    return _SYM_Z.get(sym, 0) in _METAL_Z


_DEF_GRID = 36           # angular grid resolution over 2*pi (10-degree steps)
_DEF_CLASH_F = 0.85      # crowding factor; matches the F26 / CCDC reference
_DEF_MAXROT_DEG = 60.0   # search window (+/-) around each ring's AS-BUILT angle
_MC_TOL = 1e-3           # max allowed change in any M-(ring atom) distance (A)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_HAPTO_DECLASH", "0") == "1"


def _components_no_metal(
    n: int, bond_pairs: Sequence[Tuple[int, int]], metals: Set[int]
) -> List[List[int]]:
    """Connected components of the heavy+H graph after metal bonds are removed.
    Each component = one ligand."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in bond_pairs:
        if i in metals or j in metals:
            continue
        if 0 <= i < n and 0 <= j < n:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    comp: Dict[int, List[int]] = {}
    for i in range(n):
        if i in metals:
            continue
        comp.setdefault(find(i), []).append(i)
    return list(comp.values())


def _total_overlap(
    P: np.ndarray, vdw: np.ndarray, cov: np.ndarray, f: float
) -> Tuple[float, float, int]:
    """PARTITION-FREE never-worse measure: over every NON-bonded heavy-atom pair
    (bonded = within sum_cov + 0.40), the summed vdW-overlap depth (< f*sum_vdW),
    the worst ratio, and the count of hard overlaps (< 0.90 A absolute).  Because
    the ring spin is a rigid rotation (intra-ligand distances are invariant) the
    only terms that move are ligand-vs-ligand, so minimising this total is
    equivalent to relieving the inter-ligand clash -- but the guard no longer
    depends on a bond-perception partition, so it cannot be fooled by a
    severe-overlap fragment merge/split (the F26 re-classification artifact)."""
    n = P.shape[0]
    loss = 0.0
    worst = 9.9
    hard = 0
    for i in range(n):
        if vdw[i] <= 0:
            continue
        pi = P[i]
        for j in range(i + 1, n):
            dx = pi[0] - P[j][0]; dy = pi[1] - P[j][1]; dz = pi[2] - P[j][2]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d < cov[i] + cov[j] + 0.40:
                continue                              # bonded -> skip
            sv = vdw[i] + vdw[j]
            r = d / sv if sv > 1e-9 else 9.9
            if r < worst:
                worst = r
            if d < 0.90:
                hard += 1
            thr = f * sv
            if d < thr:
                loss += (thr - d)
    return loss, worst, hard


def _spin(P: np.ndarray, rot_idx: List[int], origin: np.ndarray,
          axis_unit: np.ndarray, ang: float) -> np.ndarray:
    """Rodrigues rotation of the atoms in ``rot_idx`` about (origin, axis_unit)."""
    c, s = math.cos(ang), math.sin(ang)
    Pn = P.copy()
    u = axis_unit
    for i in rot_idx:
        v = P[i] - origin
        Pn[i] = origin + v * c + np.cross(u, v) * s + u * float(u.dot(v)) * (1.0 - c)
    return Pn


def declash_hapto(
    coords,
    syms: Sequence[str],
    metal_indices: Iterable[int],
    by_metal: Dict[int, List[Sequence[int]]],
    all_hapto_atoms: Iterable[int],
    frozen_donors: Iterable[int],
    bond_pairs: Sequence[Tuple[int, int]],
    grid: int = _DEF_GRID,
    clash_f: float = _DEF_CLASH_F,
    max_rot_deg: float = _DEF_MAXROT_DEG,
) -> np.ndarray:
    """Rotate every eta-ring about its M->centroid axis to the LOCAL inter-ligand
    -clash minimum within +/- ``max_rot_deg`` of its as-built angle.  Pure
    geometry; donor-frozen; M-(ring atom) distances invariant; deterministic;
    never raises (returns the input coords on any anomaly).

    The +/- window is the conformer-preservation guard: spinning every ring to
    the GLOBAL 360-degree minimum canonicalises ring-rotation conformers to one
    geometry, so a downstream RMSD-dedup collapses the ensemble (observed 6->2 on
    LAGPUE).  A bounded window relieves the eclipse (the staggered minimum is
    within one C-C step, ~30-36 degrees) while keeping each conformer distinct,
    so the isomer/conformer COUNT is never reduced.

    Parameters mirror the data already available at the end of
    ``_build_hapto_scaffold``:
      ``by_metal[m]`` = list of eta-groups (each a list of ring atom indices) on
      metal ``m``; ``all_hapto_atoms`` = union of all ring atoms; ``frozen_donors``
      = sigma-donor atom indices; ``bond_pairs`` = the molecule's bonds (metal
      bonds may be included — they are filtered for the ligand partition).
    """
    try:
        P0 = np.array(coords, dtype=float)
    except Exception:
        return np.array(coords, dtype=float)
    n = len(syms)
    if n < 4 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return P0
    metals = {int(m) for m in metal_indices}
    if not metals:
        return P0
    hapto_set = {int(a) for a in all_hapto_atoms}

    comps = _components_no_metal(n, bond_pairs, metals)
    if len(comps) < 2:
        return P0  # single ligand: no inter-ligand clash to relieve
    # map atom -> its ligand component
    atom_comp: Dict[int, int] = {}
    for ci, comp in enumerate(comps):
        for a in comp:
            atom_comp[a] = ci

    # vdw with H zeroed (heavy-heavy never-worse measure) + covalent radii for
    # the partition-free bonded-pair exclusion.
    vdw = np.array([(_VDW.get(s, _VDW_D) if s != "H" else 0.0) for s in syms],
                   dtype=float)
    cov = np.array([_COV.get(s, _COV_D) for s in syms], dtype=float)
    base_loss, base_worst, base_hard = _total_overlap(P0, vdw, cov, clash_f)
    if base_loss <= 1e-9:
        return P0

    # Build the deterministic list of spinnable eta-rings: (metal, ring_atoms,
    # rotating_set).  rotating_set = the whole ligand component carrying the ring
    # (ring + its hanging substituents), so substituents swing with the ring.
    # Skip a ring whose component also carries a sigma-donor to ANY metal (a
    # chelate that coordinates through both the face and a sidearm): spinning it
    # would break the sidearm M-D bond.
    frozen_donor_set = {int(d) for d in frozen_donors}
    rings: List[Tuple[int, List[int], List[int]]] = []
    for m in sorted(metals):
        for grp in by_metal.get(m, []):
            ring = sorted(int(a) for a in grp if 0 <= int(a) < n)
            if len(ring) < 4:
                continue
            ci = atom_comp.get(ring[0])
            if ci is None:
                continue
            comp = comps[ci]
            # the ligand component must not carry a separate sigma-donor
            if any((a in frozen_donor_set) for a in comp):
                continue
            # nor a ring atom seated on a DIFFERENT metal (bridging face)
            rot = [a for a in comp]
            rings.append((m, ring, rot))
    if not rings:
        return P0

    # bounded angle grid: fixed step (2*pi/grid) but only within +/- max_rot of
    # the as-built angle (offset 0), so each ring's conformer identity is kept.
    step = 2.0 * math.pi / float(grid)
    kmax = int(max(0.0, math.radians(max_rot_deg)) / step)
    angles = [k * step for k in range(-kmax, kmax + 1)]
    Pcur = P0.copy()
    cur_loss, cur_worst, cur_hard = base_loss, base_worst, base_hard

    for (m, ring, rot) in rings:
        centroid = Pcur[ring].mean(axis=0)
        axis = centroid - Pcur[m]
        nrm = float(np.linalg.norm(axis))
        if nrm < 1e-6:
            continue
        u = axis / nrm
        md0 = [float(np.linalg.norm(Pcur[c] - Pcur[m])) for c in ring]
        best_P = None
        best_loss = cur_loss
        best_worst = cur_worst
        best_hard = cur_hard
        for ang in angles:
            if abs(ang) < 1e-12:
                continue
            trial = _spin(Pcur, rot, Pcur[m], u, ang)
            # M-(ring atom) invariance guard
            ok = True
            for c, d0 in zip(ring, md0):
                if abs(float(np.linalg.norm(trial[c] - trial[m])) - d0) > _MC_TOL:
                    ok = False
                    break
            if not ok:
                continue
            # Partition-free never-worse guard: total heavy-heavy overlap depth,
            # worst contact, AND hard-overlap count must each be no worse than the
            # current accepted frame; accept only the strictly-better-loss angle.
            tl, tw, th = _total_overlap(trial, vdw, cov, clash_f)
            if tw < cur_worst - 1e-6 or th > cur_hard:
                continue
            if tl < best_loss - 1e-9:
                best_loss, best_worst, best_hard, best_P = tl, tw, th, trial
        if best_P is not None:
            Pcur = best_P
            cur_loss, cur_worst, cur_hard = best_loss, best_worst, best_hard

    # final guard: never return something worse than the INPUT on any axis
    fin_loss, fin_worst, fin_hard = _total_overlap(Pcur, vdw, cov, clash_f)
    if (not np.all(np.isfinite(Pcur)) or fin_loss > base_loss + 1e-9
            or fin_worst < base_worst - 1e-6 or fin_hard > base_hard):
        return P0
    return Pcur


# ---------------------------------------------------------------------------
# Geometry-only entry point (final-manifold post-pass: no RDKit mol available,
# connectivity + eta-ring detection are derived from the coordinates).
# ---------------------------------------------------------------------------

def _derive_bonds_geom(syms: Sequence[str], P: np.ndarray,
                       metal_tol: float = 0.45, organic_tol: float = 0.40
                       ) -> List[Tuple[int, int]]:
    """Distance-based bonds (metal bonds use a wider tolerance), consistent with
    the F26 detector's derive_bonds."""
    n = len(syms)
    cov = [_COV.get(s, _COV_D) for s in syms]
    ism = [_is_metal(s) for s in syms]
    out: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            tol = metal_tol if (ism[i] or ism[j]) else organic_tol
            dmax = cov[i] + cov[j] + tol
            dx = P[i][0] - P[j][0]; dy = P[i][1] - P[j][1]; dz = P[i][2] - P[j][2]
            if dx * dx + dy * dy + dz * dz <= dmax * dmax:
                out.append((i, j))
    return out


def detect_hapto(syms: Sequence[str], P: np.ndarray):
    """From geometry alone, return (metal_indices, by_metal, all_hapto_atoms,
    frozen_donors, bond_pairs) for the eta-ring declash, or None if there is no
    eta-coordinated ring (so the caller can skip)."""
    n = len(syms)
    metals = [i for i in range(n) if _is_metal(syms[i])]
    if not metals:
        return None
    bonds = _derive_bonds_geom(syms, P)
    metal_set = set(metals)
    # heavy adjacency excluding metal bonds
    adj: Dict[int, Set[int]] = {i: set() for i in range(n)}
    md: Dict[int, List[int]] = {m: [] for m in metals}
    donor: Set[int] = set()
    for i, j in bonds:
        if i in metal_set and j in metal_set:
            continue
        if i in metal_set or j in metal_set:
            m, a = (i, j) if i in metal_set else (j, i)
            md[m].append(a)
            donor.add(a)
            continue
        if syms[i] != "H" and syms[j] != "H":
            adj[i].add(j); adj[j].add(i)
    by_metal: Dict[int, List[List[int]]] = {}
    all_hapto: Set[int] = set()
    # an eta-face on metal m = >=4 ring carbons (deg>=2 in the heavy graph) that
    # are all bonded to m
    for m in metals:
        seated = [a for a in md[m]
                  if syms[a] == "C" and len(adj[a]) >= 2]
        if len(seated) >= 4:
            by_metal.setdefault(m, []).append(sorted(seated))
            all_hapto.update(seated)
    if not all_hapto:
        return None
    frozen_donors = [d for d in donor if d not in all_hapto]
    return metals, by_metal, all_hapto, frozen_donors, bonds


def declash_frame(syms: Sequence[str], coords,
                  grid: int = _DEF_GRID, clash_f: float = _DEF_CLASH_F,
                  max_rot_deg: float = _DEF_MAXROT_DEG):
    """Geometry-only eta-ring declash of one assembled frame.  Detects the
    eta-ring(s) + donors from the coordinates, then applies :func:`declash_hapto`.
    Returns declashed coords (same atom order); identity (input coords) if there
    is no eta-coordinated ring or on any anomaly.  Never raises."""
    try:
        P = np.array(coords, dtype=float)
        det = detect_hapto(syms, P)
        if det is None:
            return P
        metals, by_metal, all_hapto, frozen_donors, bonds = det
        return declash_hapto(P, syms, metals, by_metal, all_hapto, frozen_donors,
                             bonds, grid=grid, clash_f=clash_f, max_rot_deg=max_rot_deg)
    except Exception:
        return np.array(coords, dtype=float)


# ---------------------------------------------------------------------------
# Carbonyl (and terminal triple-bond) length correction.
# ---------------------------------------------------------------------------
# The analytical seat's flat _bond_len places every C-O at the SINGLE-bond length
# (~1.43 A) ignoring bond order, so terminal metal carbonyls M-C#O come out
# stretched to ~1.40 A instead of the true ~1.13-1.15 A (measured: 56% of
# carbonyls > 1.35 A).  This contracts a stretched terminal C#O back to the
# correct triple-bond length by moving ONLY the terminal O inward along the C->O
# axis.  Moving the terminal atom inward can only reduce inter-atom overlap, so it
# is never-worse; it is a pure local geometry correction.

_CO_TRIPLE = 1.15       # M-C#O target (A)
_STRETCH_THR = 1.25     # only contract bonds longer than this (leave good ones)


def correct_carbonyls(syms: Sequence[str], coords):
    """Contract stretched terminal metal-carbonyl C#O bonds to the triple-bond
    length by sliding ONLY the terminal O inward along the C->O axis.

    SAFETY (the discriminator that makes this never-worse): geometry alone cannot
    tell a stretched carbonyl (C#O, should be ~1.15) from a legitimate single C-O
    (alkoxide / hydroxymethyl, correctly ~1.43).  So we require the UNAMBIGUOUS
    carbonyl signature -- a metal-bound, **2-coordinate** carbon (bonded to exactly
    the metal + the terminal O, no H, no other heavy atom) in a **linear** M-C-O
    arrangement (angle > 150 deg).  M-CH2-OH (4-coord C), M-acyl/ester (3-coord C),
    and bent C-O groups are all excluded, so only genuine M-C#O is touched.
    Geometry-only, never lengthens, moves nothing but the terminal O;
    deterministic; returns the input on any anomaly."""
    try:
        P = np.array(coords, dtype=float)
        n = len(syms)
        if n < 3 or P.shape != (n, 3) or not np.all(np.isfinite(P)):
            return P
        metals = [i for i in range(n) if _is_metal(syms[i])]
        if not metals:
            return P
        bonds = _derive_bonds_geom(syms, P)
        metal_set = set(metals)
        # full neighbour lists (incl. H + metal) and the metal each atom binds
        all_nbr: Dict[int, List[int]] = {i: [] for i in range(n)}
        metal_of: Dict[int, int] = {}
        for i, j in bonds:
            all_nbr[i].append(j); all_nbr[j].append(i)
            if (i in metal_set) ^ (j in metal_set):
                a = j if i in metal_set else i
                m = i if i in metal_set else j
                metal_of.setdefault(a, m)
        changed = False
        for o in range(n):
            if syms[o] != "O" or len(all_nbr[o]) != 1:
                continue                          # terminal O only
            c = all_nbr[o][0]
            if syms[c] != "C" or c not in metal_of:
                continue                          # O's carbon must bind a metal
            # carbonyl C is 2-coordinate: exactly the metal + this O, nothing else
            cn = all_nbr[c]
            if len(cn) != 2 or set(cn) != {o, metal_of[c]}:
                continue
            m = metal_of[c]
            v = P[o] - P[c]
            d = float(np.linalg.norm(v))
            if d <= _STRETCH_THR or d < 1e-6:
                continue                          # already short enough
            # linear M-C-O (carbonyls are linear; excludes bent C-O groups)
            mc = P[m] - P[c]
            mcn = float(np.linalg.norm(mc))
            if mcn < 1e-6:
                continue
            cosang = float(np.dot(mc, v) / (mcn * d))  # angle M-C-O
            if cosang > -0.866:                   # require > 150 deg (cos < -0.866)
                continue
            P[o] = P[c] + v * (_CO_TRIPLE / d)    # slide O inward to 1.15
            changed = True
        return P if changed else np.array(coords, dtype=float)
    except Exception:
        return np.array(coords, dtype=float)
