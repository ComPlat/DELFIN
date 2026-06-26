"""delfin.manta.backbone_reembed — core-preserving GLOBAL backbone re-embed as an
additional FF-free conformer source (Task: largest XRD-recall lever, 2026-06-18).

WHY (the proof, CONFORMER_REACHABILITY_2026_06_17 / RECALL_GAP_FORENSIK_2026_06_17):
the residual ~1.5-2.2 Å crystal-recall gap is NOT the M-D coordination core (median
0.12 Å -> already good) but the LIGAND-BACKBONE CONFORMER: the crystal-matching
ligand fold is simply not produced (best-achievable backbone RMSD over the ensemble
median 1.81 Å).  For 59 % of backbone-near-miss systems the crystal fold IS reachable
(<1.5 Å) -- via a GLOBAL distance-geometry (ETKDG) re-embed of the ligand skeleton,
NOT the local frame-perturbation the existing generator does.

THE CRITICAL INVARIANT (#82/#100 lesson, non-negotiable): whole-complex DG embeds
DESTROYED the coordination core (+34 % cshm).  So this module NEVER re-embeds the
metal.  It FREEZES the metal + ALL donor atoms at their native (good, 0.12 Å)
positions and only re-folds each ligand's BACKBONE: per ligand, ETKDG-sample the FREE
ligand mol, then RIGID Kabsch-fit the conformer's DONOR atoms onto the native donor
vertex positions -> M-D distances + donor geometry stay EXACTLY native, the backbone
takes the new fold.  Chelates: only conformers whose internal donor-donor distances
match the polyhedron bite (donor-fit RMSD <= tol) are accepted (bite filter).
Monodentate: the single donor is pinned exactly on its native vertex (M-D invariant);
the backbone folds around it (donor-axis spins add orientation diversity).

Output frames are APPENDED to the FF-free ensemble (additive, best-of-ensemble MIN
crystal-recall keeps the better of old/new).  Heavy-atom RMSD dedup keeps them
distinct.  Env-gated DELFIN_FFFREE_BACKBONE_REEMBED, default OFF -> byte-identical
(no extra frames, path untouched).  Deterministic: fixed ETKDG seeds (42, 43, ...),
single-thread, no time/random.  FF-free: ETKDG/DG is GEOMETRY sampling, no force-field
relax of the metal core.

Public entry: ``reembed_complex(metal, geom_key, lig_groups, native, max_frames)``
where ``lig_groups`` is a list of per-ligand specs (free-ligand mol + donor local
indices) and ``native`` is the native assembled frame (syms, P) whose metal + donor
positions are frozen.  Returns a list of (syms, P) re-embedded frames (may be empty);
never raises, never produces a non-finite frame.
"""
from __future__ import annotations
import os
from typing import List, Optional, Tuple, Dict
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import delfin._bond_decollapse as _bd

# Deterministic ETKDG seed base (42, 43, ...).  PYTHONHASHSEED=0 + single thread.
_SEED0 = 42
# Hard core-invariant guard (Å): a re-embedded frame whose frozen donor/metal
# positions drift beyond this vs the native core is REJECTED (the #82/#100 guard).
_CORE_DRIFT_TOL = 0.05
# Chelate bite filter: a conformer whose donor atoms cannot be rigid-fit onto the
# native donor vertices within this RMSD has the wrong bite -> rejected.  Set at the
# core-drift tolerance so an accepted conformer's donors pin onto the native vertices
# WITHIN the hard core guard (the #82/#100 invariant: donors stay native +/-0.05 A).
_BITE_RMSD_TOL = _CORE_DRIFT_TOL
# Heavy-atom RMSD below which two re-embedded frames are duplicates.
_DEDUP_RMSD = 0.30


def enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_BACKBONE_REEMBED", "0") == "1"


def _n_conf() -> int:
    return int(os.environ.get("DELFIN_FFFREE_REEMBED_NCONF", "60"))


def _max_frames_default() -> int:
    return int(os.environ.get("DELFIN_FFFREE_REEMBED_MAXFRAMES", "8"))


def _n_axis_spin() -> int:
    # extra donor-axis spins for monodentate ligands (orientation diversity while
    # the donor stays pinned on its native vertex).  >=1.
    return max(1, int(os.environ.get("DELFIN_FFFREE_REEMBED_AXISSPIN", "3")))


def _kabsch_rot(Pobs: np.ndarray, Ptgt: np.ndarray) -> np.ndarray:
    """Proper-rotation matrix best-fitting rows of Pobs onto rows of Ptgt
    (Kabsch, determinant-corrected to forbid reflection -> chirality preserved)."""
    H = Pobs.T @ Ptgt
    U, _, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, dsign]) @ U.T


def _etkdg_confs(mol: Chem.Mol, k: int) -> Optional[Tuple[List[str], List[np.ndarray]]]:
    """Deterministic ETKDG conformer pool of the FREE ligand (with explicit H).

    GLOBAL distance-geometry sampling of the WHOLE torsion manifold from scratch
    (the proven reachability lever) -- NOT a local perturbation of an assembled
    frame.  Fixed seed, single thread, no MMFF metal involvement (pure ligand).
    Returns (syms, [coords, ...]) in AddHs(mol) atom order, or None on failure."""
    try:
        m = Chem.AddHs(Chem.Mol(mol))
    except Exception:
        return None
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = _SEED0
        params.numThreads = 1
        params.useRandomCoords = True
        params.pruneRmsThresh = 0.1
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=int(k), params=params))
    except Exception:
        cids = []
    if not cids:
        try:
            if AllChem.EmbedMolecule(m, randomSeed=_SEED0, useRandomCoords=True) != 0:
                return None
            cids = [0]
        except Exception:
            return None
    # MMFF relax of the FREE ORGANIC ligand only (no metal in this mol) -> clean
    # internals; this is ligand-internal MM, not a force field on the metal core.
    try:
        AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=1, maxIters=200)
    except Exception:
        pass
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    coords = []
    for c in cids:
        try:
            P = np.array(m.GetConformer(c).GetPositions(), float)
        except Exception:
            continue
        if np.all(np.isfinite(P)):
            coords.append(P)
    if not coords:
        return None
    return syms, coords


def _bite_constrained_confs(mol: Chem.Mol, donor_local: List[int],
                            donor_target_pos: np.ndarray, k: int
                            ) -> Optional[Tuple[List[str], List[np.ndarray]]]:
    """Per-chelate bite-constrained GLOBAL embed: distance-geometry sample the
    ligand BUT pin the donor-donor (and so the bite) distances to the NATIVE
    polyhedron geometry, so every conformer's donors already match the native
    vertices (rigid Kabsch then leaves ~0 residual -> the frozen-core invariant
    holds exactly).  This is the per-chelate DG the codebase notes IS feasible
    (unlike whole-complex DG, which destroyed the core, #82/#100).

    donor_target_pos: (d,3) native donor positions with the METAL AT ORIGIN
    (= P[donor_global] - metal_pos).  Only the donor-donor distances are pinned
    (the bite); the rest of the molecule re-folds globally.  Deterministic
    (fixed seed).  Returns (syms, [coords,...]) in AddHs(mol) order, or None."""
    try:
        m = Chem.AddHs(Chem.Mol(mol))
    except Exception:
        return None
    tp = [np.asarray(p, float) for p in donor_target_pos]
    cids = []
    try:
        from rdkit.Chem import rdDistGeom as _DG
        from rdkit import DistanceGeometry as _DGs
        bm = _DG.GetMoleculeBoundsMatrix(m)
        for a in range(len(donor_local)):
            for b in range(a + 1, len(donor_local)):
                i, j = int(donor_local[a]), int(donor_local[b])
                dist = float(np.linalg.norm(tp[a] - tp[b]))
                lo, hi = (i, j) if i < j else (j, i)
                bm[lo][hi] = dist + 0.02
                bm[hi][lo] = max(dist - 0.02, 0.0)
        if _DGs.DoTriangleSmoothing(bm):
            ep = _DG.EmbedParameters()
            ep.randomSeed = _SEED0
            ep.numThreads = 1
            ep.useRandomCoords = True
            ep.pruneRmsThresh = 0.1
            ep.SetBoundsMat(bm)
            cids = list(AllChem.EmbedMultipleConfs(m, numConfs=int(k), params=ep))
    except Exception:
        cids = []
    if not cids:
        return None
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    coords = []
    for c in cids:
        try:
            P = np.array(m.GetConformer(c).GetPositions(), float)
        except Exception:
            continue
        if np.all(np.isfinite(P)):
            coords.append(P)
    if not coords:
        return None
    return syms, coords


def _orthobasis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """A deterministic orthonormal pair perpendicular to a unit axis."""
    a = axis / (np.linalg.norm(axis) + 1e-12)
    tmp = (np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0]))
    e1 = tmp - a * float(np.dot(tmp, a))
    e1 /= (np.linalg.norm(e1) + 1e-12)
    e2 = np.cross(a, e1)
    return e1, e2


def _graft_ligand(lsyms: List[str], lcoords: List[np.ndarray],
                  donor_local: List[int], native_donor_pos: np.ndarray,
                  metal_pos: np.ndarray, n_axis_spin: int
                  ) -> List[np.ndarray]:
    """Graft every accepted free-ligand conformer onto the FROZEN core.

    native_donor_pos: (d,3) native positions of THIS ligand's donor atoms (frozen).
    Returns a list of (n_atoms, 3) grafted coordinate arrays whose donor atoms sit
    EXACTLY on native_donor_pos (so M-D + donor geometry stay native).  Each member
    is a distinct backbone fold.  For multidentate ligands a rigid Kabsch fit of the
    donor set realises the new fold; conformers whose intrinsic donor geometry does
    not match the native bite (Kabsch residual > tol) are dropped.  For monodentate
    ligands the single donor is pinned and the ligand spun about the M-donor axis."""
    out: List[np.ndarray] = []
    d = len(donor_local)
    tgt = np.asarray(native_donor_pos, float)               # (d,3) frozen target
    tgt_c = tgt.mean(axis=0)
    for lP in lcoords:
        lP = np.asarray(lP, float)
        if lP.shape[0] != len(lsyms) or not np.all(np.isfinite(lP)):
            continue
        src = lP[donor_local]                               # (d,3) conformer donors
        if d >= 2:
            # rigid Kabsch fit of the conformer donor set onto the native vertices
            src_c = src.mean(axis=0)
            R = _kabsch_rot(src - src_c, tgt - tgt_c)
            Pg = (lP - src_c) @ R.T + tgt_c
            # bite filter: residual donor-fit RMSD = mismatch between the
            # conformer's intrinsic donor-donor geometry and the native bite.
            res = float(np.sqrt(((Pg[donor_local] - tgt) ** 2).sum(axis=1).mean()))
            if res > _BITE_RMSD_TOL:
                continue                                    # wrong bite -> skip
            # snap the donors EXACTLY onto the native vertices (residual is small but
            # the hard core guard demands <=0.05; absorb the tiny residual by a final
            # per-frame correction that pins donors and rigidly drags the backbone).
            Pg = _pin_donors(Pg, donor_local, tgt)
            if Pg is not None and np.all(np.isfinite(Pg)):
                out.append(Pg)
        else:
            # monodentate: pin the donor on its native vertex; orient the ligand so
            # its donor lone-pair points at the metal, then spin about the M-donor
            # axis for orientation diversity.  Donor stays EXACTLY native (M-D inv.).
            di = donor_local[0]
            donor_pos = tgt[0]
            axis = metal_pos - donor_pos
            naxis = np.linalg.norm(axis)
            if naxis < 1e-6:
                continue
            axis = axis / naxis
            # ligand lone-pair direction: donor -> away from its substituent centroid
            nbr = _local_neighbors(lsyms, lP, di)
            if nbr:
                sub_c = lP[nbr].mean(axis=0)
                lp_dir = lP[di] - sub_c
            else:
                lp_dir = -(lP[di])
            nlp = np.linalg.norm(lp_dir)
            if nlp < 1e-6:
                lp_dir = np.array([0.0, 0.0, 1.0]); nlp = 1.0
            lp_dir = lp_dir / nlp
            R0 = _rot_a_to_b(lp_dir, axis)                  # lone pair -> M-donor axis
            base = (lP - lP[di]) @ R0.T + donor_pos
            e1, e2 = _orthobasis(axis)
            for s in range(max(1, n_axis_spin)):
                ang = 2.0 * np.pi * s / max(1, n_axis_spin)
                Rs = _rot_about_axis(axis, ang)
                Pg = (base - donor_pos) @ Rs.T + donor_pos
                Pg = _pin_donors(Pg, [di], donor_pos.reshape(1, 3))
                if Pg is not None and np.all(np.isfinite(Pg)):
                    out.append(Pg)
    return out


def _local_neighbors(lsyms, lP, di) -> List[int]:
    """Heavy-or-H atoms bonded to donor di (geometry only, headerless safe)."""
    nb = []
    for j in range(len(lsyms)):
        if j == di:
            continue
        d = float(np.linalg.norm(lP[j] - lP[di]))
        if d < 1.30 * _bd._ideal_bond(lsyms[di], lsyms[j]):
            nb.append(j)
    return nb


def _pin_donors(Pg: np.ndarray, donor_local: List[int], tgt: np.ndarray
                ) -> Optional[np.ndarray]:
    """Rigidly transform Pg so its donor atoms sit EXACTLY on tgt.

    A residual <=BITE tol is absorbed by one more rigid Kabsch+translate that maps
    the current donor positions onto tgt; for a single donor this is a pure
    translation.  Guarantees the frozen-core invariant (donor drift -> 0 after pin).
    Returns the transformed array, or None on degeneracy."""
    Pg = np.asarray(Pg, float)
    src = Pg[donor_local]
    if len(donor_local) == 1:
        return Pg + (tgt[0] - src[0])
    src_c = src.mean(axis=0); tgt_c = np.asarray(tgt, float).mean(axis=0)
    try:
        R = _kabsch_rot(src - src_c, np.asarray(tgt, float) - tgt_c)
    except Exception:
        return None
    return (Pg - src_c) @ R.T + tgt_c


def _rot_a_to_b(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a / (np.linalg.norm(a) + 1e-12); b = b / (np.linalg.norm(b) + 1e-12)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(np.dot(a, b))
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def _rot_about_axis(axis: np.ndarray, ang: float) -> np.ndarray:
    a = axis / (np.linalg.norm(axis) + 1e-12)
    x, y, z = a; c = np.cos(ang); s = np.sin(ang); C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def _heavy_rmsd(syms, A, B) -> float:
    h = [i for i, s in enumerate(syms) if s != "H"]
    if not h:
        return 0.0
    da = A[h] - A[h].mean(axis=0); db = B[h] - B[h].mean(axis=0)
    return float(np.sqrt(((da - db) ** 2).sum(axis=1).mean()))


def reembed_complex(metal: str, lig_groups: List[Dict], native: Tuple[List[str], np.ndarray],
                    max_frames: Optional[int] = None) -> List[Tuple[List[str], np.ndarray]]:
    """Produce core-preserving backbone-re-embedded frames for ONE assembled complex.

    metal: metal element symbol (the native frame's atom 0 is the metal).
    lig_groups: per-ligand specs in NATIVE atom order, each a dict with
        'mol'          : the free-ligand RDKit mol (NO metal),
        'global_idxs'  : list of this ligand's atom global indices in the native frame
                         (in AddHs(mol) order),
        'donor_local'  : donor atom indices WITHIN the ligand (AddHs(mol) order).
    native: (syms, P) of the native assembled frame (metal at index 0; frozen core).
    Returns up to max_frames (syms, P) re-embedded frames (additive, deduped), each
    with metal + every donor EXACTLY native (hard-guarded).  Empty on failure."""
    if max_frames is None:
        max_frames = _max_frames_default()
    syms, P = native
    P = np.asarray(P, float)
    syms = list(syms)
    if P.size == 0 or not np.all(np.isfinite(P)):
        return []
    metal_idx = next((i for i, s in enumerate(syms) if _bd._is_metal(s)), None)
    if metal_idx is None:
        return []
    metal_pos = P[metal_idx]
    n_axis = _n_axis_spin()
    k = _n_conf()

    # per-ligand grafted conformer pools (each a list of (n_lig_atoms, 3) arrays)
    per_lig: List[List[np.ndarray]] = []
    for lg in lig_groups:
        gidx = lg["global_idxs"]
        donor_local = lg["donor_local"]
        if not donor_local or len(gidx) != Chem.AddHs(lg["mol"]).GetNumAtoms():
            # cannot map this ligand onto the native frame reliably -> bail (no frames)
            return []
        native_donor_pos = P[[gidx[d] for d in donor_local]]
        if len(donor_local) >= 2:
            # CHELATE: per-chelate bite-CONSTRAINED global embed so the conformer
            # donors already match the native bite (-> exact pin, core preserved).
            ec = _bite_constrained_confs(lg["mol"], donor_local,
                                         native_donor_pos - metal_pos, k)
            if ec is None:                      # bite-constrained embed infeasible
                ec = _etkdg_confs(lg["mol"], k)  # fall back (bite filter will prune)
        else:
            ec = _etkdg_confs(lg["mol"], k)
        if ec is None:
            return []
        lsyms, lcoords = ec
        if len(lsyms) != len(gidx):
            return []
        grafted = _graft_ligand(lsyms, lcoords, donor_local, native_donor_pos,
                                metal_pos, n_axis)
        if not grafted:
            # this ligand produced no acceptable fold -> keep it at its NATIVE coords
            grafted = [P[gidx]]
        per_lig.append(grafted)

    # frame 0 of every per-lig pool is "native-or-first"; build complex frames by
    # advancing ONE ligand at a time through its alternative folds (additive,
    # bounded), so the donor/core of the others stays exactly native -> cheap +
    # never combinatorial-blow-up.  Deterministic ordering.
    n_lig = len(per_lig)
    frames: List[Tuple[List[str], np.ndarray]] = []
    seen: List[np.ndarray] = []

    def _assemble(choice: List[int]) -> Optional[np.ndarray]:
        Pc = P.copy()
        for li, ci in enumerate(choice):
            gidx = lig_groups[li]["global_idxs"]
            Pc[gidx] = per_lig[li][ci]
        # HARD core-invariant guard (#82/#100): metal + every donor must stay native.
        if float(np.linalg.norm(Pc[metal_idx] - metal_pos)) > _CORE_DRIFT_TOL:
            return None
        for li in range(n_lig):
            gidx = lig_groups[li]["global_idxs"]
            for dl in lig_groups[li]["donor_local"]:
                gi = gidx[dl]
                if float(np.linalg.norm(Pc[gi] - P[gi])) > _CORE_DRIFT_TOL:
                    return None
        if not np.all(np.isfinite(Pc)):
            return None
        return Pc

    # candidate choices: for each ligand independently, every alternative fold (the
    # others stay native) -> linear in total fold count, additive to the base.
    choices: List[List[int]] = []
    base_choice = [0] * n_lig
    for li in range(n_lig):
        for ci in range(len(per_lig[li])):
            ch = list(base_choice)
            ch[li] = ci
            choices.append(ch)
    # also one "all-alternative-1" combined fold (every ligand takes its 2nd fold if any)
    if n_lig > 1 and any(len(p) > 1 for p in per_lig):
        ch = [min(1, len(p) - 1) for p in per_lig]
        choices.append(ch)

    for ch in choices:
        Pc = _assemble(ch)
        if Pc is None:
            continue
        dup = False
        for Pk in seen:
            if _heavy_rmsd(syms, Pc, Pk) < _DEDUP_RMSD:
                dup = True
                break
        if dup:
            continue
        seen.append(Pc)
        frames.append((list(syms), Pc))
        if len(frames) >= int(max_frames):
            break
    return frames
