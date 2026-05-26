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


def _embed_metallacycle(lmol, donor_idxs, metal_sym, k=6):
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
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, randomSeed=SEED,
                                               numThreads=1, useRandomCoords=False))
        if not cids:
            if AllChem.EmbedMolecule(mh, randomSeed=SEED, useRandomCoords=True) != 0:
                return None
            cids = [0]
        keep = [i for i in range(mh.GetNumAtoms()) if i != mi]      # drop placeholder metal
        lsyms = [mh.GetAtomWithIdx(i).GetSymbol() for i in keep]
        metal_pos = {cid: np.array(mh.GetConformer(cid).GetPositions()[mi], float)
                     for cid in cids}
        # MMFF-relax the ligand INTERNALS with the donors FIXED.  The metallacycle embed
        # gets the chelating ring SHAPE right but leaves ETKDG-rough internals (funcgroup
        # planarity, bond lengths, H geometry); MMFF cannot touch the placeholder metal
        # (no params), so remove it and relax the metal-free ligand while pinning the
        # donor atoms -> internals recover to MMFF-ideal WITHOUT losing the chelating
        # donor arrangement.  Donor indices (< mi) are unchanged by the metal removal,
        # and the ligand atom order then equals `keep` / lsyms.
        ligrw = Chem.RWMol(mh); ligrw.RemoveAtom(mi); ligmol = ligrw.GetMol()
        coords = []
        for cid in cids:
            try:
                props = AllChem.MMFFGetMoleculeProperties(ligmol)
                ff = AllChem.MMFFGetMoleculeForceField(ligmol, props, confId=cid) if props else None
                if ff is not None:
                    for di in donor_idxs:
                        ff.AddFixedPoint(int(di))
                    ff.Minimize(maxIts=200)
            except Exception:
                pass
            P = np.array(ligmol.GetConformer(cid).GetPositions(), float)
            coords.append(P - metal_pos[cid])                       # metal -> origin
        return lsyms, coords
    except Exception:
        return None


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
    best = None
    for perm in itertools.permutations(range(len(targets))):
        Varr = np.array([Vt[p] for p in perm])
        R = _kabsch_rot(u, Varr)
        resid = float(np.sum((u @ R.T - Varr) ** 2))
        if best is None or resid < best[0]:
            best = (resid, R)
    Q = lP @ best[1].T
    obs = sum(float(np.linalg.norm(Q[d])) for d in donor_idxs)
    tgt = sum(float(np.linalg.norm(T)) for T in targets)
    return Q * (tgt / obs if obs > 1e-6 else 1.0)


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
        for i in fixed_idx:
            ff.AddFixedPoint(i)
        ff.Initialize()
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
    for li, va in by_lig.items():
        lg = ligands[li]
        dons = lg["donor_local_idxs"]
        # Chelates: embed the metallacycle (ligand + dative-bonded metal) so the
        # ring forms with correct geometry (backbone clears the metal) instead of
        # rigid-fitting a non-chelating free-ligand conformer.  Fall back to the
        # free-ligand path if the metallacycle embed fails.
        ring_confs = None
        lmol = None
        if lg["denticity"] >= 2:
            ring_confs = _embed_metallacycle(lg["mol"], dons[:lg["denticity"]], metal)
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
                # the exact donor<->vertex seating is found geometrically in orient.
                dent = lg["denticity"]
                dons_d = dons[:dent]
                verts = [v for v, arm in sorted(va, key=lambda x: x[1])]
                targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                           * MSB.md_distance(metal, lsyms[dons_d[i]]) for i in range(dent)]
                Q = None
                if ring_confs is not None:
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets)
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
        # freeze the donor atoms at their vertices
        for d in dons:
            fixed.add(pos + d)
        pos += len(lsyms)
    P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
    if refine:
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
