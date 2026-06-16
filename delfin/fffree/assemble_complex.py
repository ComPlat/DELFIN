#!/usr/bin/env python3
"""assemble_complex.py — assemble a full 3D TMC from a metal + ligand(s) +
geometry, by orienting each ligand so its donor sits on the placed polyhedron
vertex with its lone pair pointing at the metal.  Metal-FF-free: the sphere is
geometric (metal_sphere_builder); ligands are rigidly oriented; (constrained MMFF
relax of the organic periphery with the core frozen = next step).

This ties the pieces together: enumeration (which donors / which vertex) +
sphere placement (where) + this (orient + merge) -> a DFT-startable structure.

Prototype: homoleptic monodentate [M(L)n].  Deterministic.
"""
from __future__ import annotations
import os
import itertools
from typing import List, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from delfin.fffree import metal_sphere_builder as MSB
import delfin._bond_decollapse as _bd

SEED = 42


def _rot_align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix mapping unit vector a -> unit vector b (Rodrigues)."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(np.dot(a, b))
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def _ligand_3d(smiles: str):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=SEED)
    AllChem.MMFFOptimizeMolecule(m)
    P = m.GetConformer().GetPositions()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return syms, np.array(P, float), m


def _kabsch_rot(Pobs: np.ndarray, Ptgt: np.ndarray) -> np.ndarray:
    """Proper-rotation matrix best-fitting rows of Pobs onto rows of Ptgt
    (Kabsch, determinant-corrected to forbid reflection)."""
    H = Pobs.T @ Ptgt
    U, _, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, dsign]) @ U.T


def _embed_metallacycle(lmol, donor_idxs, metal_sym, k=6, donor_target_pos=None):
    """Embed a chelating ligand TOGETHER with a placeholder metal bonded to its
    donor atoms, so the metallacycle ring forms with correct ring geometry.  A
    free-ligand embed yields the (energetically preferred) extended/anti conformer
    whose backbone, once the donor-donor vector is rigid-fit onto the metal, buckles
    INTO the coordination shell; embedding the actual M-N-...-N ring instead places
    the backbone on the ring's far side (carbons ~2.5-2.9 A from M, not ~1.7).

    Distance-geometry only needs the M-donor bond lengths (covalent-radii sums), so
    the real metal symbol is used as the placeholder; NO MMFF (no metal params ->
    would distort M-D).  Returns (lsyms, [coords,...], metal_local_idx) where lsyms +
    coords EXCLUDE the placeholder metal, are RECENTERED so the metal sits at the
    origin (donor positions are then the M->donor vectors), and match AddHs(lmol)
    atom order (donor indices preserved), or None on failure.  Deterministic."""
    try:
        rw = Chem.RWMol(lmol)
        mi = rw.AddAtom(Chem.Atom(metal_sym))
        for di in donor_idxs:
            rw.AddBond(int(di), mi, Chem.BondType.DATIVE)   # donor->metal: preserves donor valence/H
        m = rw.GetMol()
        try:
            Chem.SanitizeMol(m, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES))
        except Exception:
            pass
        mh = Chem.AddHs(m)
        # σ-tail fix (env DELFIN_FFFREE_CHELATE_BITE, default OFF): pin the donor-donor
        # (and M-donor) distances to the IDEAL polyhedron geometry of the assigned
        # vertices so the metallacycle embeds with the correct bite angle (e.g.
        # octahedral cis = 90°) instead of the ligand's free natural bite -> donors land
        # on the exact vertices -> clean coordination shape (CShM->0).  Per-chelate DG
        # bounds (feasible; the WHOLE-complex DG was not).  Validated: bite 80-101° ->
        # 87-94°.  Falls back to the unconstrained embed if infeasible/fails.
        cids = []
        if (donor_target_pos is not None and len(donor_target_pos) == len(donor_idxs)
                and os.environ.get("DELFIN_FFFREE_CHELATE_BITE", "0") == "1"):
            try:
                from rdkit.Chem import rdDistGeom as _DG
                from rdkit import DistanceGeometry as _DGs
                nA = mh.GetNumAtoms()
                bm = _DG.GetMoleculeBoundsMatrix(mh)
                dset = {int(d) for d in donor_idxs}
                tp = [np.asarray(p, float) for p in donor_target_pos]

                def _setb(i, j, dist, tol=0.05):
                    lo, hi = (i, j) if i < j else (j, i)
                    bm[lo][hi] = float(dist + tol)
                    bm[hi][lo] = float(max(dist - tol, 0.0))
                for i in range(nA):                       # reset metal->non-donor permissive
                    if i != mi and i not in dset:
                        bm[min(i, mi)][max(i, mi)] = 100.0
                        bm[max(i, mi)][min(i, mi)] = 1.2
                for a, da in enumerate(donor_idxs):       # pin M-D + donor-donor to ideal vertices
                    _setb(mi, int(da), float(np.linalg.norm(tp[a])))
                    for b in range(a + 1, len(donor_idxs)):
                        _setb(int(da), int(donor_idxs[b]), float(np.linalg.norm(tp[a] - tp[b])))
                if _DGs.DoTriangleSmoothing(bm):
                    ep = _DG.EmbedParameters()
                    ep.randomSeed = SEED
                    ep.useRandomCoords = True
                    ep.SetBoundsMat(bm)
                    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, params=ep))
            except Exception:
                cids = []
        if not cids:                                      # default / fallback: unconstrained embed
            cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, randomSeed=SEED,
                                                   numThreads=1, useRandomCoords=False))
        if not cids:
            if AllChem.EmbedMolecule(mh, randomSeed=SEED, useRandomCoords=True) != 0:
                return None
            cids = [0]
        keep = [i for i in range(mh.GetNumAtoms()) if i != mi]      # drop placeholder metal
        lsyms = [mh.GetAtomWithIdx(i).GetSymbol() for i in keep]
        coords = []
        for cid in cids:
            P = np.array(mh.GetConformer(cid).GetPositions(), float)
            coords.append(P[keep] - P[mi])                          # metal -> origin
        return lsyms, coords
    except Exception:
        return None


def _has_collapsed_heavy_bonds(syms, P, factor=0.70):
    """True if any heavy-heavy non-metal bonded pair sits below `factor` × the
    covalent-sum ideal — catches the YILNUF-class oxalate-embed failure where the
    DG metallacycle places C-C and C-O backbone bonds below the chemistry-possible
    threshold and the self-gate later rejects the whole build to legacy.

    Iter-32e (User 2026-05-28): a post-orient sanity check so a bad embed can fall
    back to the rigid-fit path BEFORE poisoning the assembled coords.  Universal,
    geometry-only, no SMILES knowledge.  Env-gated default OFF: byte-identical
    when DELFIN_FFFREE_CHELATE_REJECT_COLLAPSED unset.
    """
    if os.environ.get("DELFIN_FFFREE_CHELATE_REJECT_COLLAPSED", "0") != "1":
        return False
    n = len(syms)
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            ideal = _bd._ideal_bond(syms[i], syms[j])
            # only check pairs that ARE bonded (not random non-bonded close contacts)
            if d > 1.30 * ideal:
                continue
            if d < factor * ideal:
                return True
    return False


def _canonical_arm_order(lg, dent):
    """Return the chelate's donor-local indices reordered into a CANONICAL,
    instance-independent arm order so that the enumerator's arm index `i`
    refers to the SAME physical donor (by element) for every instance of a
    ligand type.

    Why this matters: ``decompose`` builds ``donor_local_idxs`` / ``donor_elems``
    in raw SMILES atom order, so two chemically-identical asymmetric chelates
    (e.g. an N,O-glycinate) can land with OPPOSITE arm ordering (one [O,N], the
    other [N,O]).  ``enumerate_chelate_configs`` labels both as the SAME ligand
    ``type`` and enumerates arm permutations assuming arm index `a` maps to a
    fixed element.  If the assembly seats arm `a` -> ``donor_local_idxs[a]``
    (raw order), the two instances interpret arm indices oppositely, so the
    enumerated configs no longer biject onto the distinct element-stereoisomers:
    some collapse to duplicates and others (the homo-trans ones) are never built
    -> the ~42% coordination-isomer 3D-collapse.

    Canonical order = sort donor arms by ``(element, original-local-index)`` so
    the i-th arm is deterministic and element-consistent across instances —
    matching the honest coverage detector's element-sorted arm convention.
    Returns the reordered ``dons_d`` (length ``dent``); deterministic.  The
    legacy raw order is restored byte-identically with
    DELFIN_LEGACY_CHELATE_SEAT=1."""
    dons = list(lg["donor_local_idxs"])[:dent]
    if os.environ.get("DELFIN_LEGACY_CHELATE_SEAT", "0") == "1":
        return dons
    elems = lg.get("donor_elems") or [None] * len(dons)
    elems = list(elems)[:dent]
    try:
        order = sorted(range(dent),
                       key=lambda i: (str(elems[i]) if i < len(elems) else "",
                                      int(dons[i])))
        return [dons[i] for i in order]
    except Exception:
        return dons


def _orient_chelate_to_vertices(lP, donor_idxs, targets, asym=True):
    """Rotate a metal-centered chelate conformer (from _embed_metallacycle) so its
    donors seat onto the target vertex directions, then per-donor rescale to the
    ideal M-donor distance.  The ring geometry (backbone clears the metal) is
    preserved as a rigid body.  Works for any denticity.

    CONFIG-FAITHFUL seating (default for ASYMMETRIC chelates; DELFIN_LEGACY_CHELATE_SEAT
    unset): donor_idxs[i] (= the config's arm i) is seated on targets[i] in FIXED
        correspondence (donor i -> target i), so the enumerated arm->vertex
        assignment is realised exactly.  Without this, a Kabsch best-permutation
        collapses every asymmetric-chelate config that differs only in arm seating
        onto the SAME geometry (a duplicate) and never realises the homo-trans
        configs (N-trans-N / O-trans-O / all-trans) -> ~42% of coordination
        isomers are lost as 3D duplicates.  A single rigid Kabsch rotation onto the
        fixed correspondence preserves the chelate's internal geometry (bite angle,
        backbone) — it only chooses the orientation, never permutes/distorts.

    SYMMETRIC chelates (``asym=False``, e.g. an all-S thioether crown or
        ethylenediamine): ALL arm->vertex permutations are the SAME stereoisomer,
        so the fixed correspondence gives no coverage benefit but can ADD ring
        strain on a rigid (macrocyclic / tridentate) backbone.  These keep the
        best-fit (lowest-residual permutation) seating — strictly better geometry,
        coverage-neutral.  (Gate evidence: bis-kappa3 all-S crowns HUKMEF/AQADIF
        went 0 -> 11 isolated-atom faults under forced seating; symmetric-scoping
        removes that with no coverage loss.)

    LEGACY seating (DELFIN_LEGACY_CHELATE_SEAT=1): always the best-fit permutation,
        byte-identical to the pre-fix behaviour for the ON-vs-OFF gate / escape hatch."""
    dvecs = [np.asarray(lP[d], float) for d in donor_idxs]
    nrm = [float(np.linalg.norm(v)) for v in dvecs]
    if any(x < 1e-6 for x in nrm):
        return None
    u = np.array([v / x for v, x in zip(dvecs, nrm)])
    Vt = [np.asarray(T, float) / np.linalg.norm(T) for T in targets]
    tgt_md = [float(np.linalg.norm(T)) for T in targets]
    legacy = (os.environ.get("DELFIN_LEGACY_CHELATE_SEAT", "0") == "1"
              or not asym)        # symmetric chelate -> best-fit (no isomer to lose)
    if legacy:
        best = None
        for perm in itertools.permutations(range(len(targets))):
            Varr = np.array([Vt[p] for p in perm])
            R = _kabsch_rot(u, Varr)
            resid = float(np.sum((u @ R.T - Varr) ** 2))
            if best is None or resid < best[0]:
                best = (resid, R, perm)
        R = best[1]; perm = best[2]
    else:
        # config-faithful: fixed correspondence donor i -> target i (identity perm),
        # one rigid proper-rotation Kabsch fit (no permutation search).  Falls back
        # to the legacy best-fit on any numerical failure (never crashes).
        try:
            perm = tuple(range(len(targets)))
            Varr = np.array([Vt[p] for p in perm])
            R = _kabsch_rot(u, Varr)
        except Exception:
            best = None
            for p_ in itertools.permutations(range(len(targets))):
                Varr = np.array([Vt[p] for p in p_])
                R_ = _kabsch_rot(u, Varr)
                resid = float(np.sum((u @ R_.T - Varr) ** 2))
                if best is None or resid < best[0]:
                    best = (resid, R_, p_)
            R = best[1]; perm = best[2]
    Q = lP @ R.T
    # Per-donor RADIAL placement at the exact ideal M-donor distance (NOT a uniform
    # scale, which preserved the ETKDG embed's M-D asymmetry -> over-contracted donors,
    # the FEKZON CCDC defect).  Keep each donor's Kabsch-rotated DIRECTION (so the embed's
    # natural bite angle is preserved) and set only its radius to md.  The constrained
    # relax (donors fixed here) then pulls the backbone into consistency.
    for i, di in enumerate(donor_idxs):
        r = float(np.linalg.norm(Q[di]))
        if r > 1e-6:
            Q[di] = Q[di] / r * tgt_md[perm[i]]
    return Q


def _donor_and_lp(syms, P, mol, donor_idx: int) -> np.ndarray:
    """Lone-pair direction at the donor = away from the centroid of its neighbours."""
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(donor_idx).GetNeighbors()]
    if not nbrs:
        return np.array([1.0, 0, 0])
    v = np.zeros(3)
    for n in nbrs:
        u = P[n] - P[donor_idx]; v += u / np.linalg.norm(u)
    lp = -v
    nn = np.linalg.norm(lp)
    return lp / nn if nn > 1e-6 else np.array([1.0, 0, 0])


def _axis_rot(axis: np.ndarray, theta: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis); c = np.cos(theta); s = np.sin(theta)
    x, y, z = a
    return np.array([
        [c + x*x*(1-c),   x*y*(1-c)-z*s, x*z*(1-c)+y*s],
        [y*x*(1-c)+z*s, c + y*y*(1-c),   y*z*(1-c)-x*s],
        [z*x*(1-c)-y*s, z*y*(1-c)+x*s, c + z*z*(1-c)]])


def _subtree(mol, start, blocked):
    """Atoms reachable from ``start`` over bonds without crossing ``blocked``."""
    seen = {start}; stack = [start]
    while stack:
        a = stack.pop()
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            j = nb.GetIdx()
            if j == blocked or j in seen:
                continue
            seen.add(j); stack.append(j)
    return seen


def _vsepr_reconstruct(lsyms, lP, lmol, di):
    """Re-pyramidalise the donor's LOCAL geometry to ideal VSEPR with one
    coordination vacancy for the metal, rigidly dragging each substituent's
    subtree so substituents point AWAY from the metal.

    Fixes donors placed with their free-ligand geometry: e.g. a planar sp2
    carbanion (–CH2– with Si+H+H) whose two H end up on the M–D bond axis, or
    any donor whose H/substituents point at the metal — the dominant source of
    the coordination-angle / H-anomaly deficit.  Returns ``(modified_lP,
    vacancy_direction)``; the caller aligns the vacancy at the metal.

    Falls back to the plain lone-pair direction (no change) for ring,
    hypervalent (>=4 substituents), single-substituent (already linear) or
    degenerate donors.  Universal, geometry-only.  Disable via
    DELFIN_FFFREE_DONOR_VSEPR=0."""
    if os.environ.get("DELFIN_FFFREE_DONOR_VSEPR", "1") == "0":
        return lP, _donor_and_lp(lsyms, lP, lmol, di)
    atom = lmol.GetAtomWithIdx(di)
    nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
    k = len(nbrs)
    if k == 0 or k >= 4 or atom.IsInRing():
        return lP, _donor_and_lp(lsyms, lP, lmol, di)
    lP = np.array(lP, float).copy()
    d = lP[di]
    u = []
    for ni in nbrs:
        w = lP[ni] - d; nw = np.linalg.norm(w)
        if nw < 1e-6:
            return lP, _donor_and_lp(lsyms, lP, lmol, di)
        u.append(w / nw)
    u = np.array(u)
    if k == 1:
        return lP, -u[0]                       # linear: metal opposite, no move
    # substituent subtrees must be disjoint (else a ring not through the donor)
    subs = [_subtree(lmol, nbrs[i], di) for i in range(k)]
    seen = set()
    for s in subs:
        if seen & s:
            return lP, _donor_and_lp(lsyms, lP, lmol, di)
        seen |= s
    # ideal angle of each substituent from the metal vacancy, by hybridisation
    hyb = str(atom.GetHybridization())
    theta = {"SP": np.radians(180.0), "SP2": np.radians(120.0)}.get(hyb, np.radians(109.47))
    # vacancy axis a = where the metal goes; prefer the lone-pair sum, fall back
    # to the substituent-plane normal when the donor is planar (sum ~ 0).
    a = -u.sum(axis=0); na = np.linalg.norm(a)
    if na < 0.20 and k >= 2:
        a = np.cross(u[0], u[1]); na = np.linalg.norm(a)
    if na < 1e-6:
        return lP, _donor_and_lp(lsyms, lP, lmol, di)
    a = a / na
    tmp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = tmp - a * float(np.dot(tmp, a)); e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    moved = lP.copy()
    for i, ni in enumerate(nbrs):
        ip = u[i] - a * float(np.dot(u[i], a))     # azimuthal component
        nip = np.linalg.norm(ip)
        if nip < 1e-6:
            phi = i * (2 * np.pi / k)               # on-axis: spread evenly
            azim = np.cos(phi) * e1 + np.sin(phi) * e2
        else:
            azim = ip / nip
        target = np.cos(theta) * a + np.sin(theta) * azim
        target = target / np.linalg.norm(target)
        R = _rot_align(u[i], target)               # rotate this subtree about donor
        for j in subs[i]:
            moved[j] = (lP[j] - d) @ R.T + d
    return moved, a


def assemble_monodentate(metal: str, ligand_smiles: str, donor_idx: int,
                         geometry: str, thetas=None) -> Tuple[List[str], np.ndarray]:
    """Assemble; thetas[i] = rotation of ligand i about its own M-D axis (free DOF
    used for clash relief — geometric, metal-FF-free)."""
    ref = MSB._ref_vectors(geometry)
    n = len(ref)
    lsyms, lP, lmol = _ligand_3d(ligand_smiles)
    donor_elem = lsyms[donor_idx]
    md = MSB.md_distance(metal, donor_elem)
    if thetas is None:
        thetas = [0.0] * n
    out_syms = [metal]
    blocks = [np.zeros((1, 3))]
    lp = _donor_and_lp(lsyms, lP, lmol, donor_idx)
    for i in range(n):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        vertex = Vunit * md
        R = _rot_align(lp, -Vunit)
        Q = (lP - lP[donor_idx]) @ R.T            # donor at origin, lp along -Vunit
        Q = (Q) @ _axis_rot(-Vunit, thetas[i]).T  # spin about the M-D axis
        Q = Q + vertex
        out_syms += lsyms
        blocks.append(Q)
    return out_syms, np.vstack(blocks)


def _min_interligand(syms, P, n_lig, lig_natoms):
    """Min heavy-heavy distance between atoms of DIFFERENT ligands (clash proxy)."""
    blocks = []  # (start,end) per ligand in P (atom 0 is metal)
    s = 1
    for _ in range(n_lig):
        blocks.append((s, s + lig_natoms)); s += lig_natoms
    mind = 1e9
    for a in range(n_lig):
        for b in range(a + 1, n_lig):
            for i in range(*blocks[a]):
                if syms[i] == "H": continue
                for j in range(*blocks[b]):
                    if syms[j] == "H": continue
                    d = float(np.linalg.norm(P[i] - P[j]))
                    if d < mind: mind = d
    return mind


def clash_relief(metal, ligand_smiles, donor_idx, geometry, grid=12, passes=3):
    """Coordinate-descent over per-ligand M-D-axis rotations to maximize the min
    inter-ligand heavy distance. Deterministic (fixed grid + order)."""
    ref = MSB._ref_vectors(geometry); n = len(ref)
    lsyms, _, _ = _ligand_3d(ligand_smiles)
    lig_natoms = len(lsyms)
    angles = [2 * np.pi * k / grid for k in range(grid)]
    thetas = [0.0] * n
    for _ in range(passes):
        improved = False
        for i in range(n):
            best_t, best_score = thetas[i], -1.0
            for t in angles:
                trial = list(thetas); trial[i] = t
                syms, P = assemble_monodentate(metal, ligand_smiles, donor_idx, geometry, trial)
                sc = _min_interligand(syms, P, n, lig_natoms)
                if sc > best_score:
                    best_score, best_t = sc, t
            if abs(best_t - thetas[i]) > 1e-9:
                thetas[i] = best_t; improved = True
        if not improved:
            break
    return assemble_monodentate(metal, ligand_smiles, donor_idx, geometry, thetas), thetas


def _bite_aware_targets(lP, d1, d2, T1, T2):
    """Contract the ideal vertex targets T1,T2 to the chelate's NATURAL bite
    (donor-donor distance from the ligand conformer), keeping each M-D distance
    and the vertex-pair bisector + plane.  A chelate's natural bite angle is a
    real structural feature (e.g. ethylenediamine ~78 deg, not the ideal 90 deg
    cis-edge): forcing donors onto the exact ideal vertices over-stretches the
    donor-donor distance, so the rigid ring buckles its backbone INWARD toward
    the metal -> shape-outlier / over-coordination.  Placing the donors at the
    natural bite keeps the ring relaxed and the backbone outside the shell, while
    the coordination stays realistic.  Only CONTRACTS (tight chelates); wide
    chelates keep the ideal vertices.  Universal, geometry-only, deterministic."""
    b_nat = float(np.linalg.norm(lP[d1] - lP[d2]))
    d_vert = float(np.linalg.norm(T1 - T2))
    if not (1e-6 < b_nat < d_vert):
        return T1, T2                          # wide/degenerate -> ideal vertices
    r1 = float(np.linalg.norm(T1)); r2 = float(np.linalg.norm(T2))
    if r1 < 1e-6 or r2 < 1e-6:
        return T1, T2
    u1, u2 = T1 / r1, T2 / r2
    bis = u1 + u2
    nb = np.linalg.norm(bis)
    pdir = u1 - u2
    pdir = pdir - (pdir @ (bis / nb)) * (bis / nb) if nb > 1e-9 else pdir
    npd = np.linalg.norm(pdir)
    if nb < 1e-9 or npd < 1e-9:
        return T1, T2                          # collinear donors -> can't contract in-plane
    bis /= nb; pdir /= npd
    # angle between the two donors that yields donor-donor distance == b_nat
    cos_t = (r1 * r1 + r2 * r2 - b_nat * b_nat) / (2.0 * r1 * r2)
    theta = float(np.arccos(np.clip(cos_t, -1.0, 1.0)))
    h = theta / 2.0
    T1n = r1 * (np.cos(h) * bis + np.sin(h) * pdir)
    T2n = r2 * (np.cos(h) * bis - np.sin(h) * pdir)
    return T1n, T2n


def _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2):
    """Rigid-fit a chelating ligand's donors onto targets T1,T2; return placed
    coords. (donor-donor vector -> target vector, midpoints matched, backbone
    rotated away from metal at origin.)

    Targets are first contracted to the chelate's natural bite (see
    _bite_aware_targets) so a tight ring is not over-stretched onto the ideal
    vertices.  The axial sweep then maximises the MINIMUM metal-distance of the
    non-donor heavy atoms (not the centroid distance): donors sit ON the rotation
    axis so only the backbone moves, and pushing its CLOSEST atom as far from the
    metal as possible keeps backbone atoms out of the coordination shell (the
    self-gate's cshm picks the cn closest atoms; a single intruder -> shape
    outlier / over-coordination).  Universal across all chelate geometries,
    deterministic."""
    T1, T2 = _bite_aware_targets(lP, d1, d2, T1, T2)
    Q = lP.copy()
    R1 = _rot_align(Q[d2] - Q[d1], T2 - T1)
    Q = (Q - Q[d1]) @ R1.T + Q[d1]
    Q = Q + (0.5 * (T1 + T2) - 0.5 * (Q[d1] + Q[d2]))
    axis = T2 - T1
    mid = 0.5 * (T1 + T2)
    body_idx = [i for i in range(len(lsyms))
                if i not in (d1, d2) and lsyms[i] != "H"]
    if not body_idx:
        return Q                              # diatomic chelate: no backbone to rotate
    best, bestQ = -1e9, Q
    for k in range(36):
        Qr = (Q - mid) @ _axis_rot(axis, 2 * np.pi * k / 36).T + mid
        score = min(float(np.linalg.norm(Qr[i])) for i in body_idx)   # closest backbone atom to metal
        if score > best:
            best, bestQ = score, Qr
    return bestQ


def assemble_multichelate(metal: str, geometry: str, chelate_specs):
    """chelate_specs: list of (smiles, [d1,d2], [v1,v2]).  Place several chelating
    ligands on their assigned cis vertex-pairs (e.g. tris-chelate [M(en)3])."""
    ref = MSB._ref_vectors(geometry)
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    for smi, dons, verts in chelate_specs:
        lsyms, lP, lmol = _ligand_3d(smi)
        d1, d2 = dons
        T1 = ref[verts[0]] / np.linalg.norm(ref[verts[0]]) * MSB.md_distance(metal, lsyms[d1])
        T2 = ref[verts[1]] / np.linalg.norm(ref[verts[1]]) * MSB.md_distance(metal, lsyms[d2])
        Q = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2)
        out_syms += lsyms; blocks.append(Q)
    return out_syms, np.vstack(blocks)


def assemble_chelate(metal: str, ligand_smiles: str, donor_indices: List[int],
                     vertex_indices: List[int], geometry: str):
    """Place a chelating ligand: rigid-fit its donor atoms onto the assigned
    polyhedron vertices (align donor-donor vector to vertex-vertex vector, match
    midpoints), then rotate about that axis so the backbone points away from the
    metal.  Bidentate; donors land near (not exactly on) vertices if the
    ligand bite != vertex spacing (a constrained relax closes that)."""
    ref = MSB._ref_vectors(geometry)
    lsyms, lP, lmol = _ligand_3d(ligand_smiles)
    d1, d2 = donor_indices[0], donor_indices[1]
    md1 = MSB.md_distance(metal, lsyms[d1]); md2 = MSB.md_distance(metal, lsyms[d2])
    T1 = ref[vertex_indices[0]] / np.linalg.norm(ref[vertex_indices[0]]) * md1
    T2 = ref[vertex_indices[1]] / np.linalg.norm(ref[vertex_indices[1]]) * md2
    # rigid-fit + backbone-away-from-metal axial sweep (single source of truth)
    bestQ = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2)
    syms = [metal] + lsyms
    P = np.vstack([np.zeros((1, 3)), bestQ])
    md_act = (float(np.linalg.norm(bestQ[d1])), float(np.linalg.norm(bestQ[d2])))
    return syms, P, md_act


def _ligand_3d_from_mol(frag_mol):
    """Embed a ligand fragment mol (heavy-atom indices preserved under AddHs)."""
    m = Chem.AddHs(frag_mol)
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return syms, m.GetConformer().GetPositions(), m


def _ligand_confs_from_mol(frag_mol, k=10):
    """UNIVERSAL multi-conformer generation for a ligand (deterministic): K diverse
    ETKDG conformers (fixed seed, single-thread) + MMFF.  Returns (syms, [coords],
    mol).  Used to pick the clash-minimal conformer per ligand at placement — a
    fundamental Layer-3 mechanism applied to every ligand, not a per-case patch."""
    m = Chem.AddHs(frag_mol)
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                           numThreads=1))
    if not cids:
        if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
            return None
        cids = [0]
    try:
        AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=1)
    except Exception:
        pass
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return syms, [np.array(m.GetConformer(c).GetPositions(), float) for c in cids], m


def _clash_count(Q, existing, syms_Q, syms_ex):
    """# heavy/H pairs between block Q and existing atoms closer than 0.7*(vdW sum)."""
    if len(existing) == 0:
        return 0
    from delfin.fffree.refine import _vdw
    c = 0
    for a in range(len(Q)):
        for b in range(len(existing)):
            d = float(np.linalg.norm(Q[a] - existing[b]))
            if d < 0.70 * (_vdw(syms_Q[a]) + _vdw(syms_ex[b])):
                c += 1
    return c


def assemble_heteroleptic_from_mols(metal: str, geometry: str, vertex_specs,
                                    refine: bool = True):
    """vertex_specs[i] = (frag_mol, donor_local_idx).  Takes ligand MOLS directly
    (preserves donor index).  refine=True runs the geometry refiner (deterministic
    coordinate descent, metal+donors frozen) to remove rigid-placement clashes."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    placed_P = [np.zeros(3)]; placed_syms = [metal]
    fixed = {0}                       # metal frozen
    pos = 1
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag)
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        md = MSB.md_distance(metal, lsyms[di])
        vertex = Vunit * md
        # UNIVERSAL: pick the conformer whose placement clashes least with the
        # metal + already-placed ligands (defect-count-driven conformer selection).
        best_Q, best_clash = None, 1e18
        for lP in coords_list:
            if len(lsyms) == 1:
                Q = vertex.reshape(1, 3)
            else:
                lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                R = _rot_align(lp, -Vunit)
                Q = (lPv - lPv[di]) @ R.T + vertex
            cl = _clash_count(Q, np.array(placed_P), lsyms, placed_syms)
            if cl < best_clash:
                best_clash, best_Q = cl, Q
            if cl == 0:
                break
        out_syms += lsyms; blocks.append(best_Q)
        for row in best_Q:
            placed_P.append(row)
        placed_syms += lsyms
        if len(lsyms) == 1:
            fixed.add(pos); pos += 1
        else:
            fixed.add(pos + di); pos += len(lsyms)
    P = np.vstack(blocks)
    if refine:
        try:
            from delfin.fffree.refine import refine as _refine
            P = _refine(out_syms, P, fixed)
        except Exception:
            pass
    return out_syms, P


def _complex_rmsd(syms, Pa, Pb):
    """Heavy-atom RMSD between two SAME-topology complex frames (identity
    correspondence; both built from the same atom ordering).  Translation-only
    aligned on the metal+heavy centroid (the core is already rigid/identical, so a
    full Kabsch is unnecessary and would mask genuine conformer differences)."""
    heavy = [k for k, s in enumerate(syms) if s != "H"]
    if not heavy:
        heavy = list(range(len(syms)))
    A = Pa[heavy]; B = Pb[heavy]
    A = A - A.mean(axis=0); B = B - B.mean(axis=0)
    R = _kabsch_rot(A, B)
    return float(np.sqrt(((A @ R.T - B) ** 2).sum(axis=1).mean()))


def assemble_heteroleptic_ensemble(metal: str, geometry: str, vertex_specs,
                                   n_frames: int = 6, rmsd_dedup: float = 0.5,
                                   per_lig_confs: int = 6, refine: bool = True):
    """Ensemble variant of ``assemble_heteroleptic_from_mols``: instead of keeping
    only the single clash-minimal conformer per ligand, emit a small RMSD-deduped
    ENSEMBLE of full-complex frames that vary ligand internal conformation +
    substituent/co-ligand orientation while keeping the ideal coordination core
    (vertex directions + M-D distances) FIXED.

    This reuses the SAME FF-free Layer-3 machinery as the single path —
    ``_ligand_confs_from_mol`` (deterministic ETKDG conformer pool + MMFF),
    ``_vsepr_reconstruct``/``_rot_align`` orientation, ``_clash_count`` scoring, and
    the geometric ``refine`` — only it RETAINS several distinct low-clash conformer
    *combinations* across ligands rather than collapsing to one.  Best-of-ensemble
    crystal-recall then gets a comparable conformer spray to the legacy multi-frame
    path.  Deterministic (fixed SEED, single-thread), graph-only, never non-finite.

    Returns a list of (syms, P) frames (>=1) or None on failure.  The first frame is
    byte-identical to ``assemble_heteroleptic_from_mols`` (the clash-minimal pick).
    """
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    # Per-vertex candidate placements: for each ligand build its conformer pool and
    # orient every conformer onto the vertex; keep the distinct low-clash candidates
    # (clash vs the metal only -> placement-order-independent + diverse).  An
    # all-monatomic ligand set has a single rigid placement (-> 1 frame).
    out_syms = [metal]
    fixed = {0}
    per_vertex_cands = []     # per vertex: list of (Q, internal_clash) sorted best-first
    vertex_lsyms = []
    pos = 1
    metal_P = np.zeros((1, 3))
    metal_sym = [metal]
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag, k=max(per_lig_confs, 1))
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        md = MSB.md_distance(metal, lsyms[di])
        vertex = Vunit * md
        cand = []                                   # (Q, clash_vs_metal)
        seen_local = []
        for lP in coords_list:
            if len(lsyms) == 1:
                Q = vertex.reshape(1, 3)
            else:
                lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                R = _rot_align(lp, -Vunit)
                Q = (lPv - lPv[di]) @ R.T + vertex
            if not np.all(np.isfinite(Q)):
                continue
            cl = _clash_count(Q, metal_P, lsyms, metal_sym)
            # dedup conformers of THIS ligand by intra-ligand RMSD (identity corr.)
            dup = False
            for Qs in seen_local:
                if Qs.shape == Q.shape:
                    c0 = Q - Q.mean(axis=0); c1 = Qs - Qs.mean(axis=0)
                    if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                        dup = True
                        break
            if dup:
                continue
            seen_local.append(Q)
            cand.append((Q, cl))
        if not cand:
            return None
        # Axial-spin rotamers (substituent/co-ligand orientation): for a near-rigid
        # donor (e.g. =S/=P thiourea/phosphine CN2 ligands) the dominant pose DOF vs
        # the crystal is the ROTATION of the whole ligand about the M-donor axis, which
        # the ETKDG internal-conformer pool does not sample.  Spin each candidate block
        # about its M-D axis (= Vunit through the donor) at deterministic increments,
        # keeping distinct low-clash variants.  Pure geometry (the linear core/M-D stay
        # fixed -> only orientation changes).  Opt-out via DELFIN_FFFREE_CN2_NOSPIN=1.
        if len(lsyms) > 1 and os.environ.get("DELFIN_FFFREE_CN2_NOSPIN", "0") != "1":
            n_spin = int(os.environ.get("DELFIN_FFFREE_CN2_SPINS", "6"))
            spun = []
            for Q0, _cl0 in list(cand):
                for s in range(1, max(n_spin, 1)):
                    ang = 2.0 * np.pi * s / max(n_spin, 1)
                    Rs = _axis_rot(Vunit, ang)
                    Qs = (Q0 - vertex) @ Rs.T + vertex   # spin about the donor vertex
                    if not np.all(np.isfinite(Qs)):
                        continue
                    cls = _clash_count(Qs, metal_P, lsyms, metal_sym)
                    dup = False
                    for Qe, _ in (cand + spun):
                        if Qe.shape == Qs.shape:
                            c0 = Qs - Qs.mean(axis=0); c1 = Qe - Qe.mean(axis=0)
                            if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                                dup = True
                                break
                    if not dup:
                        spun.append((Qs, cls))
            cand += spun
        cand.sort(key=lambda t: t[1])               # low clash-vs-metal first
        per_vertex_cands.append(cand)
        vertex_lsyms.append(lsyms)
        out_syms += lsyms
        if len(lsyms) == 1:
            fixed.add(pos); pos += 1
        else:
            fixed.add(pos + di); pos += len(lsyms)

    # Enumerate full-complex frames over the Cartesian product of per-vertex
    # conformer choices, greedily ordered (sum of candidate ranks -> the best
    # combinations first), score each by total inter-ligand clash, RMSD-dedup at the
    # complex level, keep up to n_frames.  Capped product keeps it bounded+fast.
    import itertools as _it
    rank_lists = [list(range(len(c))) for c in per_vertex_cands]
    combos = list(_it.product(*rank_lists))
    # order combos by total rank (frame 0 = all best = the single-path pick), then
    # lexicographically -> deterministic.
    combos.sort(key=lambda cb: (sum(cb), cb))
    MAX_EVAL = 64                                    # bound the build work
    frames = []                                      # (syms, P) kept (deduped)
    for cb in combos[:MAX_EVAL]:
        blocks = [np.zeros((1, 3))]
        placed_P = [np.zeros(3)]; placed_syms = [metal]
        ok = True
        for vi, ci in enumerate(cb):
            Q = per_vertex_cands[vi][ci][0]
            blocks.append(Q)
            for row in Q:
                placed_P.append(row)
            placed_syms += vertex_lsyms[vi]
        P = np.vstack(blocks)
        if not np.all(np.isfinite(P)):
            continue
        if refine:
            try:
                from delfin.fffree.refine import refine as _refine
                P = _refine(out_syms, P, fixed)
            except Exception:
                pass
        if not np.all(np.isfinite(P)):
            continue
        # complex-level RMSD dedup vs already-kept frames
        dup = False
        for _, Pk in frames:
            if Pk.shape == P.shape and _complex_rmsd(out_syms, P, Pk) < rmsd_dedup:
                dup = True
                break
        if dup:
            continue
        frames.append((list(out_syms), P))
        if len(frames) >= n_frames:
            break
    if not frames:
        return None
    return frames


def _constrained_uff_relax(ligmol, fixed_idx, max_its=500):
    """METAL-FREE constrained UFF relax (the _hapto_rigid_v2 pattern): the ligands-
    only mol (no metal, no M-D bonds) is UFF-minimized with the donor atoms FIXED
    at their placed vertex positions and inter-fragment vdW ON, so ligand
    peripheries relax + inter-ligand clashes/H-overcoord resolve WITHOUT UFF having
    to type the metal (which it can't reliably do for 4d/5d)."""
    try:
        ff = AllChem.UFFGetMoleculeForceField(ligmol, ignoreInterfragInteractions=False)
        if ff is None:
            return False
        ff.Initialize()                       # MUST precede AddFixedPoint, else the fixed
        for i in fixed_idx:                   # points are cleared -> donors drift (M-D break)
            ff.AddFixedPoint(int(i))
        ff.Minimize(maxIts=max_its)
        return True
    except Exception:
        return False


def build_and_relax(metal: str, geometry: str, vertex_specs, relax: bool = True):
    """Heteroleptic monodentate assembly + metal-free constrained relax (donors
    pinned at vertices).  vertex_specs[i] = (frag_mol, donor_local_idx).
    Returns (syms, P) = [metal] + relaxed ligand atoms.  Deterministic."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    lig = Chem.RWMol()              # ligands-only mol (no metal) for the FF
    conf_xyz = []
    fixed = []
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        emb = _ligand_3d_from_mol(frag)
        if emb is None:
            return None
        lsyms, lP, lmol = emb
        md = MSB.md_distance(metal, lsyms[di])
        vertex = Vunit * md
        if len(lsyms) == 1:
            Q = vertex.reshape(1, 3)
        else:
            lp = _donor_and_lp(lsyms, lP, lmol, di)
            R = _rot_align(lp, -Vunit)
            Q = (lP - lP[di]) @ R.T + vertex
        offset = lig.GetNumAtoms()
        for a in lmol.GetAtoms():
            lig.AddAtom(Chem.Atom(a.GetAtomicNum()))
        for b in lmol.GetBonds():
            lig.AddBond(b.GetBeginAtomIdx() + offset, b.GetEndAtomIdx() + offset,
                        b.GetBondType())
        fixed.append(offset + di)          # donor pinned at its vertex
        for row in Q:
            conf_xyz.append(np.asarray(row, float))
    if lig.GetNumAtoms() == 0:
        return None
    conf = Chem.Conformer(lig.GetNumAtoms())
    for k, xyz in enumerate(conf_xyz):
        conf.SetAtomPosition(k, [float(xyz[0]), float(xyz[1]), float(xyz[2])])
    lig.AddConformer(conf, assignId=True)
    if relax:
        try:
            Chem.SanitizeMol(lig, catchErrors=True)
        except Exception:
            pass
        _constrained_uff_relax(lig, fixed)
    LP = lig.GetConformer().GetPositions()
    syms = [metal] + [a.GetSymbol() for a in lig.GetAtoms()]
    P = np.vstack([np.zeros((1, 3)), np.array(LP, float)])
    return syms, P


def assemble_from_config(metal, geometry, config, ligands, refine=True):
    """Build a 3D complex from a chelate-isomer config (vertex -> (ligand_idx,
    arm_idx)) and the decomposed ligand list.  Chelating ligands are Kabsch-fit
    onto their two assigned vertices; monodentate ligands are oriented onto their
    vertex.  Multi-conformer selection per ligand + constrained refine.  Returns
    (syms, P) = [metal] + ligand atoms, or None on failure."""
    ref = MSB._ref_vectors(geometry)
    # group config by ligand instance: lig_idx -> [(vertex, arm), ...]
    by_lig = {}
    for v, (li, arm) in config.items():
        by_lig.setdefault(li, []).append((v, arm))
    out_syms = [metal]; placed = [np.zeros(3)]; placed_syms = [metal]
    fixed = {0}; pos = 1
    relax_frags = []          # (AddHs(lg.mol), ligands-only offset) for the internal relax
    for li, va in by_lig.items():
        lg = ligands[li]
        dons = lg["donor_local_idxs"]
        lig_offset = pos - 1                       # start index in the ligands-only frame
        # Chelates: embed the metallacycle (ligand + dative-bonded metal) so the
        # ring forms with correct geometry (backbone clears the metal) instead of
        # rigid-fitting a non-chelating free-ligand conformer.  Fall back to the
        # free-ligand path if the metallacycle embed fails.
        ring_confs = None
        lmol = None
        if lg["denticity"] >= 2:
            # ideal donor positions at the assigned polyhedron vertices (metal at origin)
            # -> the constrained metallacycle embed pins the chelate bite to match them.
            _dent = lg["denticity"]
            # CONFIG-FAITHFUL: arm index i -> the i-th CANONICAL (element-sorted)
            # donor, so the enumerator's arm->vertex assignment is realised the
            # same way for every instance of a ligand type (see _canonical_arm_order).
            _dons_d = _canonical_arm_order(lg, _dent)
            _vts = [v for v, arm in sorted(va, key=lambda x: x[1])]
            # element of the i-th canonical arm (for the ideal M-D target distance)
            _delems = [lg["mol"].GetAtomWithIdx(int(di)).GetSymbol() for di in _dons_d]
            try:
                _dtp = [ref[_vts[i]] / np.linalg.norm(ref[_vts[i]])
                        * MSB.md_distance(metal, _delems[i]) for i in range(_dent)]
            except Exception:
                _dtp = None
            ring_confs = _embed_metallacycle(lg["mol"], _dons_d, metal, donor_target_pos=_dtp)
        if ring_confs is not None:
            lsyms, coords_list = ring_confs
        else:
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
        best_Q, best_clash = None, 1e18
        for lP in coords_list:
            if lg["denticity"] == 1:
                v = va[0][0]
                Vunit = ref[v] / np.linalg.norm(ref[v])
                md = MSB.md_distance(metal, lsyms[dons[0]])
                if len(lsyms) == 1:
                    Q = (Vunit * md).reshape(1, 3)
                else:
                    lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, dons[0])
                    Q = (lPv - lPv[dons[0]]) @ _rot_align(lp, -Vunit).T + Vunit * md
            else:
                # chelate (bi-/tri-dentate): the d assigned vertices host the d donors;
                # arm index i (config order) seats on verts[i] in FIXED correspondence
                # (config-faithful), using the canonical element-sorted arm order so the
                # enumeration bijects onto the distinct stereoisomers.
                dent = lg["denticity"]
                dons_d = _canonical_arm_order(lg, dent)
                verts = [v for v, arm in sorted(va, key=lambda x: x[1])]
                targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                           * MSB.md_distance(metal, lsyms[dons_d[i]]) for i in range(dent)]
                # config-faithful seating is only meaningful for ASYMMETRIC chelates
                # (distinct donor elements); for symmetric chelates every arm seating
                # is the same isomer, so keep the best-fit (lower ring strain).
                _de = [lsyms[di] for di in dons_d]
                _asym = len(set(_de)) > 1
                Q = None
                if ring_confs is not None:
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=_asym)
                    # iter-32e (YILNUF oxalate-collapse class): the DG metallacycle
                    # embed sometimes produces a backbone with collapsed C-C / C-O bonds.
                    # _orient_chelate_to_vertices preserves that (rigid Kabsch+rescale),
                    # so the resulting Q poisons _build_is_clean and the whole complex
                    # falls back to legacy.  Reject such Q here → the bidentate rigid
                    # fallback below picks up; for tridentate we skip the config (the
                    # build-gate would reject anyway).  Env-gated default OFF.
                    if Q is not None and _has_collapsed_heavy_bonds(lsyms, Q):
                        Q = None
                if Q is None and dent == 2:         # embed/orient failed -> rigid fallback (bidentate only)
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1], targets[0], targets[1])
                if Q is None:                       # tridentate embed failure -> skip this config
                    continue
            cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
            if cl < best_clash:
                best_clash, best_Q = cl, Q
            if cl == 0:
                break
        if best_Q is None:                  # no conformer could be placed -> bail to legacy
            return None
        out_syms += lsyms
        for row in best_Q:
            placed.append(row)
        placed_syms += lsyms
        relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))   # bonds for internal relax
        # freeze the donor atoms at their vertices
        for d in dons:
            fixed.add(pos + d)
        pos += len(lsyms)
    P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
    if refine:
        # FUNDAMENTAL internal relaxation (division-of-labor doctrine): build the
        # ligands-only mol (NO metal) at the placed coords and UFF-relax it with the
        # DONORS FIXED + inter-fragment vdW ON.  All organic internals (bond lengths,
        # angles, funcgroup/aromatic planarity, H geometry) recover to MM-ideal and
        # inter-ligand clashes resolve, while the constructed coordination is preserved
        # (donors pinned -> M-D invariant; metal never enters the force field).  This
        # replaces the weak defect-count refine(), which left ETKDG-rough internals.
        # INTERNAL finishing.  THREE-WAY RACE design:
        #   Track 3 (FF-FREE, DEFAULT): pure construction -- GEOMETRIC correctors only
        #     (sp2-flatten projection + defect-count clash-relief).  No force field.
        #   Track 2 (ligand-FF, opt-in DELFIN_FFFREE_LIGANDFF=1): additionally a constrained
        #     UFF relax of the organic internals (metal + donors frozen) before the geometric
        #     correctors.  Track 1 (UFF) is the legacy path (DELFIN_FFFREE_BUILDER=0).
        # The complex mol (metal + ligand bonds) is the sp2-flatten template (and the optional
        # UFF mol); metal + donors stay frozen so the constructed coordination is preserved.
        try:
            cm = Chem.RWMol()
            cm.AddAtom(Chem.Atom(out_syms[0]))
            for frag, _off in relax_frags:
                base = cm.GetNumAtoms()
                for a in frag.GetAtoms():
                    cm.AddAtom(Chem.Atom(a.GetAtomicNum()))
                for b in frag.GetBonds():
                    cm.AddBond(b.GetBeginAtomIdx() + base, b.GetEndAtomIdx() + base, b.GetBondType())
            donor_globals = sorted(fixed - {0})
            for dg in donor_globals:
                cm.AddBond(0, dg, Chem.BondType.SINGLE)
            if cm.GetNumAtoms() == len(out_syms):
                conf = Chem.Conformer(cm.GetNumAtoms())
                for kk in range(cm.GetNumAtoms()):
                    x = P[kk]; conf.SetAtomPosition(kk, [float(x[0]), float(x[1]), float(x[2])])
                cm.AddConformer(conf, assignId=True)
                try:
                    Chem.SanitizeMol(cm, catchErrors=True)
                except Exception:
                    pass
                if os.environ.get("DELFIN_FFFREE_LIGANDFF", "0") == "1":   # Track 2: FF relax
                    if _constrained_uff_relax(cm, [0] + donor_globals):
                        Pn = np.array(cm.GetConformer().GetPositions(), float)
                        if Pn.shape == P.shape:
                            P = Pn
                            c = cm.GetConformer()
                            for kk in range(len(P)):
                                c.SetAtomPosition(kk, [float(P[kk][0]), float(P[kk][1]), float(P[kk][2])])
                # FF-FREE geometric planarity corrector (BOTH tracks): project sp2 atoms onto
                # their neighbour plane (aromatic/amide/funcgroup planarity).  Pure geometry.
                try:
                    from delfin.smiles_converter import _flatten_sp2_atoms_xyz
                    xyz_str = "\n".join(f"{s} {float(p[0]):.6f} {float(p[1]):.6f} {float(p[2]):.6f}"
                                        for s, p in zip(out_syms, P))
                    flat = _flatten_sp2_atoms_xyz(xyz_str, cm)
                    if flat:
                        newP = np.array([[float(x) for x in ln.split()[1:4]]
                                         for ln in flat.splitlines() if ln.strip()], float)
                        if newP.shape == P.shape:
                            P = newP
                except Exception:
                    pass
        except Exception:
            pass
        # FF-free geometric clash-relief (both tracks)
        try:
            from delfin.fffree.refine import refine as _refine
            P = _refine(out_syms, P, fixed)
        except Exception:
            pass
    donors = sorted(fixed - {0})              # global indices of the constructed donor atoms
    return out_syms, P, donors


def assemble_heteroleptic(metal: str, geometry: str, vertex_specs):
    """vertex_specs[i] = (ligand_smiles, donor_idx) for polyhedron vertex i.
    Heteroleptic monodentate assembly (different ligand per vertex) — the basis
    for building enumerated coordination isomers (e.g. cis/trans MA4B2)."""
    ref = MSB._ref_vectors(geometry)
    if len(vertex_specs) != len(ref):
        raise ValueError("vertex_specs count != vertices")
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    for i, (smi, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        lsyms, lP, lmol = _ligand_3d(smi)
        md = MSB.md_distance(metal, lsyms[di])
        vertex = Vunit * md
        if len(lsyms) == 1:                       # monatomic ligand (e.g. Cl-)
            out_syms += lsyms; blocks.append(vertex.reshape(1, 3)); continue
        lp = _donor_and_lp(lsyms, lP, lmol, di)
        R = _rot_align(lp, -Vunit)
        Q = (lP - lP[di]) @ R.T + vertex
        out_syms += lsyms; blocks.append(Q)
    return out_syms, np.vstack(blocks)


def generate_complex_conformers(metal, ligand_smiles, donor_idx, geometry,
                                 max_lig_conf=6):
    """Wire L3 into the complex: enumerate the LIGAND's conformers (conformer_enum),
    build the homoleptic complex from each, dedup. Demonstrates conformers-per-isomer
    (the L3 dimension of generate-gate-floor)."""
    from delfin.fffree import conformer_enum as CE
    _, confs = CE.enumerate_conformers(ligand_smiles)
    ref = MSB._ref_vectors(geometry); n = len(ref)
    out = []
    for e, m in confs[:max_lig_conf]:
        lsyms = [a.GetSymbol() for a in m.GetAtoms()]
        lP = m.GetConformer().GetPositions()
        md = MSB.md_distance(metal, lsyms[donor_idx])
        lp = _donor_and_lp(lsyms, lP, m, donor_idx)
        syms = [metal]; blocks = [np.zeros((1, 3))]
        for i in range(n):
            Vunit = ref[i] / np.linalg.norm(ref[i])
            R = _rot_align(lp, -Vunit)
            Q = (lP - lP[donor_idx]) @ R.T + Vunit * md
            syms += lsyms; blocks.append(Q)
        out.append((round(e, 2), syms, np.vstack(blocks)))
    return out
