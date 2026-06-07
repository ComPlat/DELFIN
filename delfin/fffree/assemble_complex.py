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

# Track 3 pure FF-free: skip MMFF in ligand 3D embed; rely on ETKDG initial
# coordinates + downstream refine.py (FF-free defect-loss gradient descent).
# DELFIN_FFFREE_PURE_TRACK3=1 → no FF anywhere in fffree pipeline.
# Default OFF (byte-identical when unset). Universal, deterministic, no template.
_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"


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
    if not _PURE_TRACK3:                 # Track 3 pure FF-free: skip MMFF
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
        #
        # G15c (2026-06-01, User pure-construction-mandate): auto-enable this path
        # under PT3 with RELAXED donor-donor tolerance.  The old CHELATE_BITE used
        # tight donor-donor (+/- 0.05 A), which forces the bipyridine's natural 2.7 A
        # bite up to the cis-octahedron 2.9 A vertex-distance and stretches M-D
        # asymmetrically (e.g. Pt-N3 2.07, Pt-N14 2.46 -- the 28-iter-gate-axis
        # signature).  Relaxing donor-donor to +/- 0.20 A under PT3 lets the
        # bipy stay at its natural bite while M-D lands on the COD-empirical ideal,
        # closing the M-X distortion gap at the SOURCE (no post-hoc polish brittle-
        # ness).  M-D bound stays tight (+/- 0.05) because M-D length is the
        # quantity we want preserved.  Feasibility: triangle smoothing checks
        # this is satisfiable; otherwise the existing unconstrained-embed fallback
        # still fires.
        cids = []
        _G15C_DD_RELAX = (
            os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
            or os.environ.get("DELFIN_FFFREE_DD_RELAX", "0") == "1"
        )
        if (donor_target_pos is not None and len(donor_target_pos) == len(donor_idxs)
                and (os.environ.get("DELFIN_FFFREE_CHELATE_BITE", "0") == "1"
                     or _G15C_DD_RELAX)):
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
                # M-D tolerance stays tight (M-D length is the conserved quantity).
                # Donor-donor tolerance is the new G15c relaxation knob.
                md_tol = 0.05
                dd_tol = 0.20 if _G15C_DD_RELAX else 0.05
                for i in range(nA):                       # reset metal->non-donor permissive
                    if i != mi and i not in dset:
                        bm[min(i, mi)][max(i, mi)] = 100.0
                        bm[max(i, mi)][min(i, mi)] = 1.2
                for a, da in enumerate(donor_idxs):       # pin M-D + donor-donor to ideal vertices
                    _setb(mi, int(da), float(np.linalg.norm(tp[a])), tol=md_tol)
                    for b in range(a + 1, len(donor_idxs)):
                        _setb(int(da), int(donor_idxs[b]),
                              float(np.linalg.norm(tp[a] - tp[b])), tol=dd_tol)
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


def _orient_chelate_to_vertices(lP, donor_idxs, targets):
    """Rotate a metal-centered chelate conformer (from _embed_metallacycle) so its
    donor directions best-fit the target vertex directions, then uniformly rescale
    to the ideal M-donor distance.  Tries EVERY donor<->target correspondence and
    keeps the lowest-residual Kabsch rotation, so a bi-/tri-dentate auto-seats on
    the matching cis-edge / mer / fac vertex arrangement without the enumerator
    needing to know it.  The ring geometry (backbone clears the metal) is preserved
    as a rigid body.  Works for any denticity."""
    dvecs = [np.asarray(lP[d], float) for d in donor_idxs]
    nrm = [float(np.linalg.norm(v)) for v in dvecs]
    if any(x < 1e-6 for x in nrm):
        return None
    u = np.array([v / x for v, x in zip(dvecs, nrm)])
    Vt = [np.asarray(T, float) / np.linalg.norm(T) for T in targets]
    tgt_md = [float(np.linalg.norm(T)) for T in targets]
    best = None
    for perm in itertools.permutations(range(len(targets))):
        Varr = np.array([Vt[p] for p in perm])
        R = _kabsch_rot(u, Varr)
        resid = float(np.sum((u @ R.T - Varr) ** 2))
        if best is None or resid < best[0]:
            best = (resid, R, perm)
    Q = lP @ best[1].T
    # Per-donor RADIAL placement at the exact ideal M-donor distance (NOT a uniform
    # scale, which preserved the ETKDG embed's M-D asymmetry -> over-contracted donors,
    # the FEKZON CCDC defect).  Keep each donor's Kabsch-rotated DIRECTION (so the embed's
    # natural bite angle is preserved) and set only its radius to md.  The constrained
    # relax (donors fixed here) then pulls the backbone into consistency.
    perm = best[2]
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
    if not _PURE_TRACK3:                 # Track 3 pure FF-free: skip MMFF
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
    if not _PURE_TRACK3:                 # Track 3 pure FF-free: skip MMFF
        try:
            AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=1)
        except Exception:
            pass
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return syms, [np.array(m.GetConformer(c).GetPositions(), float) for c in cids], m


def _hh_clash_include_active_assemble() -> bool:
    """``True`` iff ``DELFIN_FFFREE_HH_CLASH_INCLUDE`` is on (default OFF).

    Wired-in 2026-06-04 (healing-wiring cleanup).  Default OFF ->
    :func:`_clash_count` returns the byte-identical legacy 0.70× heavy/H
    count.  When ON, the count is augmented by cross-block H-H pairs
    below the Bondi vdW floor (0.85 × 2.40 Å), catching methyl-methyl
    eclipsing that the loose 0.70× factor passes through clean.
    """
    import os as _os_loc
    raw = _os_loc.environ.get("DELFIN_FFFREE_HH_CLASH_INCLUDE", "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def _clash_count(Q, existing, syms_Q, syms_ex):
    """# heavy/H pairs between block Q and existing atoms closer than 0.7*(vdW sum).

    H-H clash augmentation (2026-06-04, wiring cleanup)
    --------------------------------------------------
    When ``DELFIN_FFFREE_HH_CLASH_INCLUDE=1`` is set, the count is augmented
    by cross-block H-H pairs below the Bondi vdW floor (0.85 × 2 × 1.20 Å
    = 2.04 Å) -- catches methyl-methyl eclipsing between Q and the already-
    placed atoms.  Default OFF -> byte-identical with HEAD b195dba.
    """
    if len(existing) == 0:
        return 0
    from delfin.fffree.refine import _vdw
    c = 0
    for a in range(len(Q)):
        for b in range(len(existing)):
            d = float(np.linalg.norm(Q[a] - existing[b]))
            if d < 0.70 * (_vdw(syms_Q[a]) + _vdw(syms_ex[b])):
                c += 1
    # H-H augmentation (env-gated, default OFF byte-identical).
    if _hh_clash_include_active_assemble():
        try:
            from delfin.fffree.hh_clash_detector import (
                H_VDW_RADIUS,
                hh_clash_factor,
            )
            fac = hh_clash_factor()
            d_min = fac * (H_VDW_RADIUS + H_VDW_RADIUS)
            for a in range(len(Q)):
                if str(syms_Q[a]) != "H":
                    continue
                Qa = Q[a]
                for b in range(len(existing)):
                    if str(syms_ex[b]) != "H":
                        continue
                    d = float(np.linalg.norm(Qa - existing[b]))
                    # Skip exact coincidence + already-counted (d < 0.70 floor).
                    if d < 1e-9:
                        continue
                    if d < 0.70 * (_vdw(syms_Q[a]) + _vdw(syms_ex[b])):
                        continue  # already in c
                    if d < d_min:
                        c += 1
        except Exception:
            # Defence in depth: detector failure must never break assembly.
            pass
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
    # Phase B Task #63 (2026-06-03): multi-metal / cluster dispatcher.
    # When DELFIN_FFFREE_MULTI_METAL=1 (or DELFIN_FFFREE_PURE_TRACK3=1) AND the
    # decomposed config carries a multi-metal mol (>=2 metal atoms in the graph),
    # dispatch to the cluster orchestrator BEFORE the single-metal assembly.
    # Hard-rollback on contract violation: the orchestrator returns None and we
    # silently continue into the legacy single-metal path -> byte-identical when
    # the env flag is unset.
    try:
        _mm_mol = None
        if isinstance(config, dict):
            _mm_mol = config.get("__multi_metal_mol__")
        if _mm_mol is None and ligands:
            for _lg in ligands:
                if isinstance(_lg, dict):
                    _cand = _lg.get("__parent_mol__") or _lg.get("mol")
                    if _cand is not None:
                        try:
                            n_metal = sum(1 for a in _cand.GetAtoms()
                                          if _bd._is_metal(a.GetSymbol()))
                        except Exception:
                            n_metal = 0
                        if n_metal >= 2:
                            _mm_mol = _cand
                            break
        if _mm_mol is not None:
            from delfin.fffree import multi_metal_assemble as _MMA
            if _MMA.should_dispatch_multi_metal(_mm_mol):
                _mm_res = _MMA.assemble_multi_metal(_mm_mol, ligands=ligands,
                                                    metals=None, config=config)
                if _mm_res is not None:
                    return _mm_res
                # _mm_res is None -> silent fall-through to single-metal path.
    except Exception:
        # Any error in the multi-metal branch is non-fatal: we drop to the
        # legacy single-metal path. Production safety contract.
        pass

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
        ring_confs_unconstrained = None     # alt embed for the pick-better gate
        lmol = None
        if lg["denticity"] >= 2:
            # ideal donor positions at the assigned polyhedron vertices (metal at origin)
            # -> the constrained metallacycle embed pins the chelate bite to match them.
            _dent = lg["denticity"]
            _vts = [v for v, arm in sorted(va, key=lambda x: x[1])]
            _delems = lg.get("donor_elems") or [None] * _dent
            try:
                _dtp = [ref[_vts[i]] / np.linalg.norm(ref[_vts[i]])
                        * MSB.md_distance(metal, _delems[i]) for i in range(_dent)]
            except Exception:
                _dtp = None
            ring_confs = _embed_metallacycle(lg["mol"], dons[:_dent], metal, donor_target_pos=_dtp)
            # CHELATE_BITE pick-better gate (Task #55-56, 2026-06-04).
            # When DELFIN_FFFREE_CHELATE_BITE=1 + DELFIN_FFFREE_CHELATE_BITE_PICK_BETTER=1:
            # ALSO compute the unconstrained embed (donor_target_pos=None) so we can
            # later choose, per-chelate, whichever placement gives lower local CShM.
            # Filters out the catastrophic regressions (D-YOMHIR / D-DIMXUS / D-JAHZUQ)
            # while keeping the wins (ABUMET 11.97 -> 0.055 etc).  Env-gated default-OFF
            # byte-identical -- the alt embed is NEVER computed unless BOTH env flags
            # are set, so HEAD path is unchanged.
            try:
                from delfin.fffree.chelate_bite_pick_better import pick_better_active
                if pick_better_active():
                    ring_confs_unconstrained = _embed_metallacycle(
                        lg["mol"], dons[:_dent], metal, donor_target_pos=None,
                    )
            except Exception:
                ring_confs_unconstrained = None
        if ring_confs is not None:
            lsyms, coords_list = ring_confs
        else:
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
        # Phase G6 (2026-05-31): hapto-π piano-stool placement.
        # When the ligand is hapto (η3/η4/η5/η6/η7/η8 Cp/arene/allyl/diene/COT),
        # is_hapto=True from decompose.py. Place the ring centroid at the
        # polyhedron vertex, with the ring perpendicular to the M-centroid axis.
        # Universal across hapto modes via hapto_modes.m_centroid_distance().
        is_hapto_lig = lg.get("is_hapto", False)
        eta = lg.get("hapto_eta", 0)
        best_Q, best_clash = None, 1e18
        for lP in coords_list:
            if is_hapto_lig and eta >= 3:
                # Piano-stool placement: metal at vertex, ring centroid at distance
                # _md (from hapto_modes table) along the vertex unit vector.
                v = va[0][0]
                Vunit = ref[v] / np.linalg.norm(ref[v])
                # Get all hapto-donor positions (the ring carbons)
                ring_indices = lg.get("donor_local_idxs", [dons[0]])
                ring_pos = np.array([lP[i] for i in ring_indices])
                centroid = ring_pos.mean(axis=0)
                # M-centroid distance from hapto_modes (env-gated import; fallback to 1.7)
                try:
                    from delfin.fffree.hapto_modes import m_centroid_distance
                    _mc = m_centroid_distance(metal, eta)
                except Exception:
                    _mc = 1.7
                # Rotate ring so its normal aligns with Vunit, then place centroid at Vunit*_mc
                centered = ring_pos - centroid
                _, _, Vt = np.linalg.svd(centered, full_matrices=False)
                ring_normal = Vt[-1]
                # Ensure normal points outward from metal
                if np.dot(ring_normal, Vunit) < 0:
                    ring_normal = -ring_normal
                # Rotate all ligand atoms so ring-normal -> Vunit
                R = _rot_align(ring_normal, Vunit)
                lP_aligned = (lP - centroid) @ R.T
                # Translate so ring centroid sits at Vunit * _mc
                Q = lP_aligned + Vunit * _mc
            elif lg["denticity"] == 1:
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
                # the exact donor<->vertex seating is found geometrically in orient.
                dent = lg["denticity"]
                dons_d = dons[:dent]
                verts = [v for v, arm in sorted(va, key=lambda x: x[1])]
                targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                           * MSB.md_distance(metal, lsyms[dons_d[i]]) for i in range(dent)]
                Q = None
                if ring_confs is not None:
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets)
                    # iter-32e (YILNUF oxalate-collapse class): the DG metallacycle
                    # embed sometimes produces a backbone with collapsed C-C / C-O bonds.
                    # _orient_chelate_to_vertices preserves that (rigid Kabsch+rescale),
                    # so the resulting Q poisons _build_is_clean and the whole complex
                    # falls back to legacy.  Reject such Q here → the bidentate rigid
                    # fallback below picks up; for tridentate we skip the config (the
                    # build-gate would reject anyway).  Env-gated default OFF.
                    if Q is not None and _has_collapsed_heavy_bonds(lsyms, Q):
                        Q = None
                # CHELATE_BITE pick-better gate (Task #55-56, 2026-06-04).
                # When BOTH env flags set, orient the unconstrained embed too and
                # pick whichever placement gives lower local CShM (vs the assigned
                # polyhedron geometry).  Eliminates catastrophic regressions while
                # keeping wins.  Env-gated default-OFF byte-identical: the
                # ``ring_confs_unconstrained`` reference is None unless the alt
                # embed was computed above, so this block is a no-op.
                if (ring_confs_unconstrained is not None
                        and Q is not None
                        and len(lP) == len(ring_confs_unconstrained[1][0])):
                    try:
                        from delfin.fffree.chelate_bite_pick_better import (
                            pick_better as _pb,
                            md_invariant_preserved as _md_ok,
                        )
                        _alt_lsyms, _alt_coords_list = ring_confs_unconstrained
                        # Use the first unconstrained conformer (matches the current
                        # outer-loop iteration on ``lP`` from coords_list); both embeds
                        # share atom indexing because metallacycle keeps the donor
                        # local idxs invariant.
                        _alt_lP = _alt_coords_list[0]
                        _alt_Q = _orient_chelate_to_vertices(_alt_lP, dons_d, targets)
                        if _alt_Q is not None and _has_collapsed_heavy_bonds(_alt_lsyms, _alt_Q):
                            _alt_Q = None
                        if _alt_Q is not None:
                            _picked, _src, _c, _u = _pb(
                                Q, _alt_Q, dons_d, geometry,
                                metal_pos=np.zeros(3),
                            )
                            # Defence-in-depth: only swap when M-D invariant holds
                            # (orient already places donors on the targets; this
                            # guards against any future change).
                            if (_picked is not None
                                    and _md_ok(Q, _picked, dons_d, metal_pos=np.zeros(3))):
                                Q = _picked
                    except Exception:
                        # Any failure -> keep the constrained Q (legacy CHELATE_BITE
                        # behaviour).  Build never regresses on import errors.
                        pass
                if Q is None and dent == 2:         # embed/orient failed -> rigid fallback (bidentate only)
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1], targets[0], targets[1])
                if Q is None:                       # tridentate embed failure -> skip this config
                    continue
                # Iter-32e deeper fix (YILNUF oxalate-5-ring planar projector,
                # 2026-06-05).  When DELFIN_FFFREE_OXALATE_5RING_PLANAR=1, detect
                # M-O-C(=O)-C(=O)-O 5-ring oxalate chelates in this block and
                # project the 5 ring atoms onto a planar ideal (bite ~78°,
                # C-C ~1.55, C-O endo ~1.27, ring planarity <0.05 A).  Exocyclic
                # C=O is projected onto the same plane (carboxylate conjugation
                # preserved).  Universal, geometry-only, deterministic.  Env-gated
                # default-OFF byte-identical (helper returns input copy).
                try:
                    from delfin.fffree.oxalate_5ring_planar import (
                        flag_active as _ox5_active,
                        apply as _ox5_apply,
                    )
                    if _ox5_active() and dent == 2:
                        Q, _ox5_applied = _ox5_apply(
                            lsyms, Q, dons_d,
                            metal_pos=np.zeros(3), metal_sym=metal,
                        )
                except Exception:
                    # Any failure -> keep the existing Q (legacy path), build
                    # never regresses on the env-gated projector's import errors.
                    pass
            cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
            # CONSTRUCTION-FIX #3 (User 2026-06-03): build-time clash gate.
            # Penalise candidates that internally contain collapsed bonds (X-H
            # or heavy-heavy below the structqual floor) so the selection loop
            # naturally picks the conformer with the fewest intra-Q defects.
            # Env-gated default-OFF byte-identical (helper returns 0 when no
            # flag set).
            try:
                from delfin.fffree.build_time_clash_gate import (
                    _flag_active as _bcg_active, collapse_count as _bcg_count,
                )
                if _bcg_active():
                    cl = cl + _bcg_count(lsyms, Q)
            except Exception:
                pass
            # SURGICAL FIX 1 (User 2026-06-03 F24-interlig forensik):
            # inter-ligand non-bonded vdW penalty.  Adds a quadratic overlap
            # score between the candidate ``Q`` and already-placed atoms
            # ``placed`` so configurations with cis-CO axis collisions /
            # bulky-phosphine methyl clashes / π-π-too-close are
            # deprioritised.  Pure read-only score; coordinates unchanged.
            # Env-gated default-OFF byte-identical -- helper returns
            # (0.0, 0) when no flag set.
            try:
                from delfin.fffree.build_time_clash_gate import (
                    interlig_penalty_active as _ilp_active,
                    _interlig_penalty_for_pair as _ilp_pair,
                    _interlig_weight as _ilp_w,
                )
                if _ilp_active() and len(placed) > 1:
                    _ilp_score, _ilp_n = _ilp_pair(
                        lsyms, Q, placed_syms, np.array(placed),
                    )
                    cl = cl + _ilp_w() * _ilp_score
            except Exception:
                pass
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
                # CONSTRUCTION-FIX #2 (User 2026-06-03): M-D direction alignment.
                # Before sp2-flatten so the lone-pair-anti rotation re-orients
                # aromatic donors / sp3-anti donors INTO their correct M-D
                # geometry; the sp2-flatten then sees the already-correct
                # neighbour plane.  Pure axis rotation around M-D preserves
                # the donor distance EXACTLY (M-D invariant by construction).
                # Env-gated default-OFF byte-identical.
                try:
                    from delfin.fffree.mddir_alignment import (
                        align_donor_lonepairs, _flag_active as _mda_active,
                    )
                    if _mda_active():
                        P_pre = P.copy()
                        P_md, _nf, _ns = align_donor_lonepairs(
                            out_syms, P, metal_idx=0, donor_idxs=donor_globals,
                        )
                        if (P_md is not None and isinstance(P_md, np.ndarray)
                            and P_md.shape == P.shape
                            and np.all(np.isfinite(P_md))):
                            # M-D invariant defence-in-depth
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(P_pre[_dg] - P_pre[0]))
                                _d_new = float(np.linalg.norm(P_md[_dg] - P_md[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_md
                                conf = cm.GetConformer()
                                for kk in range(len(P)):
                                    conf.SetAtomPosition(kk, [float(P[kk][0]),
                                                              float(P[kk][1]),
                                                              float(P[kk][2])])
                except Exception:
                    pass
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
                # SURGICAL FIX 2 (User 2026-06-03 F24-interlig forensik):
                # terminal CO/NO/CN ligand torsion stagger.  When two
                # terminal sp-ligands occupy adjacent polyhedron vertices
                # the default placement collides their O/N atoms at d ~
                # 1.85 A; a rigid 60 deg spin of the SECOND ligand about
                # its own M-donor axis preserves r(M-donor) EXACTLY
                # while opening the terminal-terminal distance.  Donors
                # and the metal are NEVER moved.  Env-gated default-OFF
                # byte-identical (helper returns (P_copy, 0) when no
                # flag set).
                try:
                    from delfin.fffree.terminal_ligand_stagger import (
                        apply_stagger, stagger_active as _tls_active,
                    )
                    if _tls_active():
                        P_pre_t = P.copy()
                        P_t, _ntls = apply_stagger(
                            out_syms, P, metal_idx=0,
                            donor_idxs=donor_globals,
                        )
                        if (P_t is not None and isinstance(P_t, np.ndarray)
                            and P_t.shape == P.shape
                            and np.all(np.isfinite(P_t))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(P_pre_t[_dg] - P_pre_t[0]))
                                _d_new = float(np.linalg.norm(P_t[_dg] - P_t[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_t
                except Exception:
                    pass
                # CONSTRUCTION-FIX #1 (User 2026-06-03): amide-VSEPR template.
                # AFTER sp2-flatten so it cleans up the amide-N residue the
                # geometric flatten misses (sp2 angle window 100-135 deg leaves
                # some pyramidal amide-N untouched).  Pure projection of N
                # onto the 3-substituent plane.  Donors NEVER moved (M-D
                # invariant by construction).  Env-gated default-OFF.
                try:
                    from delfin.fffree.amide_vsepr_template import (
                        enforce_amide_planarity, _flag_active as _amd_active,
                    )
                    if _amd_active():
                        P_pre_a = P.copy()
                        P_a, _nfa = enforce_amide_planarity(
                            out_syms, P, metal_idx=0, donor_idxs=donor_globals,
                        )
                        if (P_a is not None and isinstance(P_a, np.ndarray)
                            and P_a.shape == P.shape
                            and np.all(np.isfinite(P_a))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(P_pre_a[_dg] - P_pre_a[0]))
                                _d_new = float(np.linalg.norm(P_a[_dg] - P_a[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_a
                except Exception:
                    pass
                # CONSTRUCTION-FIX #4 (User 2026-06-03): oxoanion-VSEPR
                # template (re-activation of iter-32a-3).  Voll-pool
                # f8c9905 verdict showed nitrate_pct_files_with_violation
                # +139 % and nitrate_pct_no3_broken +93 % vs pgcorr-v3.
                # This hook projects detected NO3-/ClO4-/SO4^2-/PO4^3-
                # oxygens onto their ideal D3h/Td template AFTER the
                # amide-VSEPR step (templates are independent --
                # amide-N never overlaps an oxoanion-X centre).  Donors
                # NEVER moved (M-D invariant by construction).
                # Env-gated default-OFF, byte-identical when unset.
                try:
                    from delfin.fffree.oxoanion_vsepr_template import (
                        enforce_oxoanion_vsepr,
                        _flag_active as _oxo_active,
                    )
                    if _oxo_active():
                        P_pre_o = P.copy()
                        P_o, _nfo = enforce_oxoanion_vsepr(
                            out_syms, P, metal_idx=0, donor_idxs=donor_globals,
                        )
                        if (P_o is not None and isinstance(P_o, np.ndarray)
                            and P_o.shape == P.shape
                            and np.all(np.isfinite(P_o))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(P_pre_o[_dg] - P_pre_o[0]))
                                _d_new = float(np.linalg.norm(P_o[_dg] - P_o[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_o
                except Exception:
                    pass
                # CONSTRUCTION-FIX (hapto-honest, 2026-06-03): rigid-block
                # correction of hapto units (η²-η⁸).  After all σ-side
                # construction fixes and BEFORE GRIP/refine, snap each
                # detected π-cluster centroid to the empirical M-centroid
                # distance for (metal, η) AND align the ring axis to the
                # M-centroid line.  Metal + σ donors NEVER moved (the
                # σ-only metric stays byte-invariant).  Defence-in-depth
                # M-D validator + rigid-block validator inside the helper.
                # Env-gated default-OFF byte-identical.  Pairs with the
                # standalone hapto-only-CShM metric (see
                # ``agent_workspace/quality_framework/scripts/
                # hapto_only_cshm.py``).
                try:
                    from delfin.fffree.hapto_honest_construction import (
                        apply_hapto_honest, honest_active as _hho_active,
                    )
                    if _hho_active():
                        P_pre_h = P.copy()
                        P_h, _nfh = apply_hapto_honest(
                            out_syms, P, metal_idx=0,
                            donor_idxs=donor_globals,
                        )
                        if (P_h is not None and isinstance(P_h, np.ndarray)
                            and P_h.shape == P.shape
                            and np.all(np.isfinite(P_h))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(P_pre_h[_dg] - P_pre_h[0]))
                                _d_new = float(np.linalg.norm(P_h[_dg] - P_h[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_h
                except Exception:
                    pass
                # CONSTRUCTION-FIX (sp³-H umbrella, 2026-06-04, GUPPY
                # ea26dcb): rebuild every sp³-C/N H umbrella that is
                # degenerate (H-X-H 90°/180° pattern observed in
                # 29-Ni_pincer-tBu-imid_.xyz where all 3 methyls had H
                # atoms on orthogonal cartesian axes).  L-BFGS angle-loss
                # cannot move 180° H pairs (saddle of cos(angle)), so this
                # correction MUST live in construction.  Metal + σ donors
                # NEVER moved; centre + heavy neighbours NEVER moved; only
                # H positions are rewritten on a Td template anchored on
                # the heavy-atom direction.  Defence-in-depth M-D check.
                # Env-gated default-OFF, byte-identical when unset.
                try:
                    from delfin.fffree.sp3_h_umbrella import (
                        enforce_all_sp3_umbrella,
                        umbrella_active as _sp3_umbrella_active,
                    )
                    if _sp3_umbrella_active():
                        P_pre_u = P.copy()
                        P_u, _rep_u = enforce_all_sp3_umbrella(
                            P, out_syms, mol=None,
                            preserve_bond_lengths=True,
                            only_degenerate=True,
                            enable_clash_rollback=True,
                        )
                        if (P_u is not None and isinstance(P_u, np.ndarray)
                            and P_u.shape == P.shape
                            and np.all(np.isfinite(P_u))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(
                                    P_pre_u[_dg] - P_pre_u[0]))
                                _d_new = float(np.linalg.norm(
                                    P_u[_dg] - P_u[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_u
                except Exception:
                    pass
                # Overnight aromatic-bond fix (2026-06-04, LUHMOT 1.13 Å
                # bug): per-bond pull of aromatic edges toward CCDC empirical
                # means (delfin.fffree.aromatic_bond_targets). Catches the
                # case where aromatic_ring_scale's mean-based / uniform pull
                # misses intra-ring variance AND topology_healing ignores
                # compressed bonds. Frozen-set: metal-coordinated + hapto-π
                # + donors. M-D invariance check before accept. Default-OFF
                # byte-identical (env DELFIN_FFFREE_AROMATIC_BONDS=1; auto
                # under PURE_TRACK3).
                try:
                    from delfin.fffree.aromatic_bond_enforcement import (
                        enforce_aromatic_bonds as _enforce_aromatic_bonds,
                        is_enabled as _aromatic_bonds_enabled,
                    )
                    if _aromatic_bonds_enabled():
                        P_pre_ab = P.copy()
                        _, P_ab = _enforce_aromatic_bonds(
                            out_syms, P, mol=None,
                            fixed=set(donor_globals) | {0},
                        )
                        if (P_ab is not None and isinstance(P_ab, np.ndarray)
                            and P_ab.shape == P.shape
                            and np.all(np.isfinite(P_ab))):
                            _md_ok = True
                            for _dg in donor_globals:
                                _d_old = float(np.linalg.norm(
                                    P_pre_ab[_dg] - P_pre_ab[0]))
                                _d_new = float(np.linalg.norm(
                                    P_ab[_dg] - P_ab[0]))
                                if abs(_d_old - _d_new) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = P_ab
                except Exception:
                    pass  # silent fallback to pre-aromatic P
                # Phase 4 (SPEC_GRIP §4.2): env-gated GRIP polish.  When
                # DELFIN_FFFREE_GRIP=1 we run the CCDC-grounded L-BFGS
                # refinement on the post-sp2-flatten geometry, BEFORE the
                # final ff-free clash-relief refine.  grip_polish is
                # internally safe (accept-if-better, M-D/topology/chiral
                # validators, NaN-guard), but we add a defence-in-depth
                # M-D check and a NaN gate here so any unexpected failure
                # cleanly rolls back to the pre-GRIP P.  Default-OFF =
                # byte-identical to HEAD (the entire block is skipped).
                if os.environ.get("DELFIN_FFFREE_GRIP", "0") == "1":
                    try:
                        from delfin.fffree.grip_polish import grip_polish
                        from delfin.fffree.grip_mogul_lookup import GripLibrary
                        _grip_lib = GripLibrary.get_default()  # may be None
                        _md_targets = [float(np.linalg.norm(P[d] - P[0])) for d in donor_globals]
                        Pg = grip_polish(
                            P, cm, metal=0, donors=donor_globals,
                            geom=str(geometry or ""), mogul_lib=_grip_lib,
                        )
                        if (
                            Pg is not None
                            and isinstance(Pg, np.ndarray)
                            and Pg.shape == P.shape
                            and np.all(np.isfinite(Pg))
                        ):
                            # Defence-in-depth M-D check (±0.05 Å, SPEC §11).
                            _md_ok = True
                            for _di, _dg in enumerate(donor_globals):
                                _d_now = float(np.linalg.norm(Pg[_dg] - Pg[0]))
                                if abs(_d_now - _md_targets[_di]) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = Pg
                    except Exception:
                        pass  # silent fallback to existing P
                # Phase 5 (User 2026-06-02, post-GRIP-hapto-VOLLPOOL forensik):
                # post-GRIP geometric correctors -- env-gated, default OFF =
                # byte-identical.  Three independent fixes that each verify the
                # M-D invariant before accepting.  See
                # delfin.fffree.post_grip_corrector for fix details (F20 sp2
                # planarity, haxis donor-H rotation, mdshort secondary-contact
                # push-back).  Activates when any of the per-fix env flags is set
                # or DELFIN_FFFREE_POST_GRIP_ALL=1.
                try:
                    from delfin.fffree.post_grip_corrector import (
                        post_grip_corrections, any_active as _pg_any_active,
                    )
                    if _pg_any_active():
                        Pc, _ = post_grip_corrections(
                            out_syms, P, metal_idx=0, donor_idxs=donor_globals,
                        )
                        if (
                            Pc is not None
                            and isinstance(Pc, np.ndarray)
                            and Pc.shape == P.shape
                            and np.all(np.isfinite(Pc))
                        ):
                            # Defence-in-depth M-D check (mirrors GRIP block).
                            _md_ok = True
                            _md_targets_pg = [
                                float(np.linalg.norm(P[d] - P[0])) for d in donor_globals
                            ]
                            for _di, _dg in enumerate(donor_globals):
                                _d_now = float(np.linalg.norm(Pc[_dg] - Pc[0]))
                                if abs(_d_now - _md_targets_pg[_di]) > 0.05:
                                    _md_ok = False
                                    break
                            if _md_ok:
                                P = Pc
                except Exception:
                    pass  # silent fallback to pre-corrector P
        except Exception:
            pass
        # FF-free geometric clash-relief (both tracks)
        try:
            from delfin.fffree.refine import refine as _refine
            P = _refine(out_syms, P, fixed)
        except Exception:
            pass
    donors = sorted(fixed - {0})              # global indices of the constructed donor atoms
    # CONSTRUCTION-SANITY GATE (2026-06-07, Bug-class 2+3):
    # Post-embedding sanity check rejects builds that contain Pauli-floor
    # violations, severely stretched/compressed bonds, drifted M-D distances
    # or a catastrophic coordination-polyhedron CShM.  When the check fires we
    # return ``None`` so the upstream caller (smiles_converter) can retry with
    # a different config / fall back to the legacy path.  Env-gated default
    # OFF byte-identical (``DELFIN_FFFREE_CONSTRUCTION_SANITY=1``).
    try:
        from delfin.fffree.construction_sanity import (
            sanity_active as _cs_active,
            assert_construction_sane as _cs_assert,
        )
        if _cs_active():
            # Derive bond list from the (placed) ligand mols + M-D bonds.
            # When the relax block built ``cm`` we re-use it; otherwise we
            # infer bonds from the per-ligand mols and the (metal, donor)
            # pairs in ``fixed``.
            _bonds_for_sanity = []
            try:
                # Per-ligand internal bonds (offset by lig_offset in the
                # global frame).
                for _frag_mol, _off in relax_frags:
                    for _b in _frag_mol.GetBonds():
                        _i = _b.GetBeginAtomIdx() + _off + 1
                        _j = _b.GetEndAtomIdx() + _off + 1
                        _bonds_for_sanity.append((_i, _j))
            except Exception:
                # relax_frags may not exist when an early return path is taken;
                # the sanity check still works without bond-stretch detection.
                _bonds_for_sanity = []
            # Metal-donor bonds
            for _dg in donors:
                _bonds_for_sanity.append((0, _dg))
            _ok, _viols = _cs_assert(
                P, out_syms, _bonds_for_sanity,
                metal_idx=0,
                donor_idxs=donors,
                geometry=str(geometry or ""),
            )
            if not _ok:
                # Reject the build -- caller falls back to the legacy path.
                return None
    except Exception:
        # Any error in the sanity gate is non-fatal: keep the build (the
        # gate is opt-in, byte-identical when disabled).
        pass
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
