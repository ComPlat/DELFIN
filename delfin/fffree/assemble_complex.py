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


def _planar_mer_cn5_enabled() -> bool:
    """Geometry-aware MERIDIONAL bite for a RIGID PLANAR tridentate on CN5 (TBP-5 /
    SPY-5).  The OC-6 meridional fix (DELFIN_FFFREE_PLANAR_MER) gives such a ligand the
    correct meridional vertex ASSIGNMENT, but on CN5 the rigid-Kabsch orient cannot
    OPEN the embed's ~110deg folded bite onto the meridional axis.  When enabled, the
    metallacycle embed is forced through the distance-geometry bite constraint (pinning
    donor-donor distances to the meridional vertices) so the embedded conformer ALREADY
    carries the meridional ~150-165deg bite.  Default OFF -> byte-identical (the forced
    bite is never applied; the SPY-5 trans-basal triples are never emitted)."""
    return os.environ.get("DELFIN_FFFREE_PLANAR_MER_CN5", "0") == "1"


def _embed_metallacycle(lmol, donor_idxs, metal_sym, k=6, donor_target_pos=None,
                        harden=False, force_bite=False):
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
        # CHELATE-BACKBONE hardening (DELFIN_FFFREE_CHELATE_BACKBONE, default OFF -> byte-id):
        # PHASE 0 of the polydentate project (K4_MACROCYCLE_DESIGN_2026_06_17.md §3.2).  For
        # a large/strained backbone (BIQCOV-class) the CHELATE_BITE HARD donor-donor pin
        # (tol 0.05) over-constrains the bounds matrix -> Triangle-Smoothing fails -> the
        # unconstrained fallback embed buckles the backbone INTO the coordination shell
        # (H-H / C-H collapses, the self-gate then drops the whole complex to legacy).  The
        # fix: keep the M-D distances HARD-pinned (so M-D stays invariant), but relax the
        # donor-donor distances to a WIDE SOFT window [d-tol, d+tol] (DELFIN_FFFREE_POLY_BB_TOL,
        # default 0.4 A) so the strained ring has folding freedom, AND embed MORE conformers
        # (DELFIN_FFFREE_POLY_K_CONFS, default 12) so the assembler can pick a collapse-free
        # backbone pose.  Activated by EITHER flag; CHELATE_BACKBONE => soft windows + more
        # confs, CHELATE_BITE-only => the historic hard pins (byte-id).  Per-donor M-D radial
        # rescale downstream (_orient_chelate_to_vertices) re-sets M-D exactly, so the soft
        # window never moves a donor off its ideal radius.
        cids = []
        # harden == True only for the newly-admitted LARGE backbones (set by the caller,
        # per-arm > historic cap, flag-gated).  Existing chelates (harden=False) keep the
        # exact historic embed (k=6, hard donor-donor pins ONLY under CHELATE_BITE).
        # force_bite (DELFIN_FFFREE_PLANAR_MER_CN5, rigid-planar CN5): activate the bite
        # constraint for THIS embed regardless of the global CHELATE_BITE flag, so the
        # rigid planar tridentate embeds with the meridional bite pinned to the assigned
        # CN5 vertices (else the rigid-Kabsch orient keeps the folded ~110deg).
        _chel_bite = os.environ.get("DELFIN_FFFREE_CHELATE_BITE", "0") == "1"
        # OC-6 mixed-donor vertex fix (DELFIN_FFFREE_OC6_VERTEX, default OFF -> byte-id):
        # ROOT CAUSE — the metallacycle DG embed bonds the placeholder metal ONLY to the
        # donor atoms, so a NON-DONOR heavy atom that shares a donor's parent (the canonical
        # case: a carboxylate's non-coordinating carbonyl O in an N,O-glycinate/amino-acid
        # chelate, but also any pendant heteroatom on a donor's first shell) has NO lower
        # bound to the metal and DG happily collapses it ONTO the metal (~0.7-1.1 A).  That
        # in-shell intruder makes _build_is_clean reject the whole complex (over-coordination
        # / overlap) -> it falls back to legacy UFF, whose mixed-donor OC-6 buckles the donors
        # together ("kein OC-6 Polyeder, donors bunched" — eye-flagged ADOHOT/ADOROD/ZUSBEU/
        # ASEBAC class).  The donor VERTEX assignment itself is already a perfect octahedron
        # (12 cis 90deg / 3 trans 180deg, verified); only the embed's non-donor collapse
        # poisons it.  Fix: a metal->non-donor-heavy EXCLUSION floor in the DG bounds matrix
        # (no non-donor heavy atom may enter the metal's coordination shell), applied via the
        # constrained-DG path -- and turned on for EVERY chelate embed (not just the
        # bite-pinned ones) so the collapse is prevented at construction.  Universal,
        # graph/geometry-only (donor set comes from the cleave), deterministic, never raises.
        _oc6_vertex = os.environ.get("DELFIN_FFFREE_OC6_VERTEX", "0") == "1"
        # exclusion floor: keep non-donor heavy atoms outside the donor shell.  Default 2.4 A
        # (~ the shortest realistic non-bonded M...heavy contact, above any M-D bond) so a
        # genuine bridging/agostic contact is not forced but the carbonyl-O collapse is barred.
        _excl = float(os.environ.get("DELFIN_FFFREE_OC6_NONDONOR_EXCL", "2.4"))
        if harden:
            k = max(int(k), int(os.environ.get("DELFIN_FFFREE_POLY_K_CONFS", "12")))
        _dd_tol = float(os.environ.get("DELFIN_FFFREE_POLY_BB_TOL", "0.4")) if harden else 0.05
        # The bite-pinned DG path needs donor targets; the OC6_VERTEX exclusion path does not
        # (it only adds metal->non-donor lower bounds, leaving donor distances to RDKit's own
        # covalent bounds, which the downstream per-donor radial rescale re-sets exactly).
        _have_targets = (donor_target_pos is not None
                         and len(donor_target_pos) == len(donor_idxs))
        _use_bm = (_have_targets and (_chel_bite or harden or force_bite)) or _oc6_vertex
        if _use_bm:
            try:
                from rdkit.Chem import rdDistGeom as _DG
                from rdkit import DistanceGeometry as _DGs
                nA = mh.GetNumAtoms()
                bm = _DG.GetMoleculeBoundsMatrix(mh)
                dset = {int(d) for d in donor_idxs}
                tp = ([np.asarray(p, float) for p in donor_target_pos]
                      if _have_targets else None)

                def _setb(i, j, dist, tol=0.05):
                    lo, hi = (i, j) if i < j else (j, i)
                    bm[lo][hi] = float(dist + tol)
                    bm[hi][lo] = float(max(dist - tol, 0.0))
                # metal->non-donor lower bound.  Without OC6_VERTEX this is the historic
                # permissive 1.2 A (byte-identical); with it the exclusion floor (no
                # non-donor heavy atom inside the coordination shell).  H atoms keep the
                # permissive floor (an X-H may legitimately point near the metal).
                _md_lo = _excl if _oc6_vertex else 1.2
                for i in range(nA):                       # reset metal->non-donor bounds
                    if i != mi and i not in dset:
                        is_h = mh.GetAtomWithIdx(i).GetAtomicNum() == 1
                        lo, hi = min(i, mi), max(i, mi)
                        this_lo = 1.2 if is_h else _md_lo
                        # raise the UPPER bound above the floor so lo<=hi stays consistent
                        bm[lo][hi] = max(float(bm[lo][hi]), this_lo + 0.1, 100.0
                                         if not _oc6_vertex else this_lo + 2.0)
                        bm[hi][lo] = float(this_lo)
                if tp is not None:
                    for a, da in enumerate(donor_idxs):   # M-D HARD; donor-donor (soft for backbone)
                        _setb(mi, int(da), float(np.linalg.norm(tp[a])))
                        for b in range(a + 1, len(donor_idxs)):
                            _setb(int(da), int(donor_idxs[b]),
                                  float(np.linalg.norm(tp[a] - tp[b])), tol=_dd_tol)
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


def _planar_polydentate_place_enabled() -> bool:
    """In-plane (coplanar-metal) PLACEMENT for a RIGID PLANAR polydentate.

    A rigid (aromatic/conjugated) planar polydentate — terpyridine-class tridentate,
    pincer — holds its >=3 donors coplanar in the ligand's own flat plane, so the
    coordinated metal physically MUST lie IN that donor plane (donors + M coplanar =
    a meridional in-plane arrangement; that is what a planar tridentate demands).

    The historic metallacycle embed (``_embed_metallacycle``) places M bonded to all
    donors but lets the rigid backbone FOLD so the metal lifts ~1.0 A OUT of the
    donors' own plane (a tripod buckle).  A post-hoc rigid rotation of a frozen
    non-coplanar donor set CANNOT fix this — rotating a rigid body never changes the
    metal's signed distance to the donor plane (``DELFIN_FFFREE_PI_COPLANAR_M``
    proved this).  So the fix is at PLACEMENT time (this flag): construct a
    metal-centered conformer in which the metal is solved INTO the rigid donor plane
    before orienting onto the assigned meridional vertices.

    Default OFF -> byte-identical (the coplanar conformer is never built; the
    historic folded embed is used)."""
    return os.environ.get("DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE", "0") == "1"


def _coplanar_metal_centered_conformer(mol, donor_idxs, metal_sym, mds, k=8):
    """Build a metal-centered conformer of a RIGID PLANAR polydentate in which the
    metal is COPLANAR with the donor set (the in-plane meridional pose a flat
    conjugated tridentate physically requires), returning the SAME contract as
    ``_embed_metallacycle``: ``(lsyms, [coords, ...])`` where ``lsyms`` + each
    ``coords`` exclude the placeholder metal, match ``AddHs(mol)`` atom order
    (donor indices preserved), and are RECENTERED so the metal sits at the ORIGIN
    (donor positions are then the M->donor vectors).  Returns ``None`` on failure.

    Method (universal, geometry/graph-only, no SMILES/refcode knowledge):
      1. Embed FREE-ligand conformers (the rigid backbone keeps the donors flat) and
         MMFF-optimise the organic internals.  The free embed gives the ligand's
         NATURAL flat donor geometry (~150-165 deg outer-outer span), NOT the
         metallacycle's folded ~110 deg.
      2. For each near-planar conformer, fit the donor/backbone mean-plane (normal
         ``n``) and solve, by 2-D Gauss-Newton IN that plane, for the metal position
         that best matches every ideal M-donor distance ``mds[i]`` simultaneously.
         The solution is IN the plane by construction => M is coplanar with the
         donors (out-of-plane = 0 exactly).
      3. Pick the conformer with the SMALLEST distance residual (the least-splayed
         flat pose — the donors that can be reached at ~equal M-D bonds), recenter so
         M is at the origin, and emit it.

    The downstream rigid Kabsch-orient + per-donor RADIAL rescale
    (``_orient_chelate_to_vertices``) then seats the donors on the assigned in-plane
    meridional vertices and sets each M-D to exactly its ideal — radial rescale moves
    each donor only ALONG its (in-plane) M->donor ray, so the metal stays coplanar
    and the meridional span is preserved.  Deterministic (fixed seed, single thread);
    never raises."""
    try:
        mh = Chem.AddHs(mol)
        cids = list(AllChem.EmbedMultipleConfs(
            mh, numConfs=max(int(k), 1), randomSeed=SEED,
            numThreads=1, useRandomCoords=False))
        if not cids:
            if AllChem.EmbedMolecule(mh, randomSeed=SEED, useRandomCoords=True) != 0:
                return None
            cids = [0]
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mh, numThreads=1)
        except Exception:
            pass
        dons = [int(x) for x in donor_idxs]
        if len(dons) < 3 or len(mds) != len(dons):
            return None
        # backbone atoms = the heavy atoms on every donor-donor shortest path (the
        # rigid conjugated framework that holds the donors coplanar)
        ring = set()
        for a in range(len(dons)):
            for b in range(a + 1, len(dons)):
                sp = Chem.GetShortestPath(mol, dons[a], dons[b])
                if not sp:
                    return None
                ring.update(int(i) for i in sp)
        ring = sorted(ring)
        if len(ring) < 4:
            return None

        def _lp_in(P, di, nrm):
            nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(di).GetNeighbors()
                    if nb.GetAtomicNum() > 1]
            if not nbrs:
                return None
            v = P[di] - np.mean([P[n] for n in nbrs], axis=0)
            w = v - float(np.dot(v, nrm)) * nrm
            nn = float(np.linalg.norm(w))
            return w / nn if nn > 1e-6 else None

        best = None                       # (residual, recentered coords)
        keep = list(range(mh.GetNumAtoms()))
        for cid in cids:
            P = np.array(mh.GetConformer(cid).GetPositions(), float)
            cen = P[ring].mean(axis=0)
            B = P[ring] - cen
            try:
                _, sv, Vt = np.linalg.svd(B)
            except Exception:
                continue
            if sv[0] < 1e-9 or (sv[2] / sv[0]) > 0.18:   # backbone not flat -> skip
                continue
            nrm = Vt[2]; e1 = Vt[0]; e2 = Vt[1]

            def _to2d(p):
                q = p - cen
                return np.array([float(np.dot(q, e1)), float(np.dot(q, e2))])

            D2 = [_to2d(P[di]) for di in dons]
            # initial M guess = mean of (donor + md * inward in-plane lone-pair)
            seeds = []
            for i, di in enumerate(dons):
                lp = _lp_in(P, di, nrm)
                if lp is None:
                    lp = -D2[i] / max(np.linalg.norm(D2[i]), 1e-6)  # fallback: toward 2-D centroid
                    lp3 = e1 * lp[0] + e2 * lp[1]
                    lp = np.array([float(np.dot(lp3, e1)), float(np.dot(lp3, e2))])
                    seeds.append(D2[i] + mds[i] * lp)
                    continue
                lp2 = np.array([float(np.dot(lp, e1)), float(np.dot(lp, e2))])
                seeds.append(D2[i] + mds[i] * lp2)
            u = np.mean(seeds, axis=0)
            for _ in range(200):           # 2-D Gauss-Newton: f_i = |u - D2_i| - md_i
                J = []; r = []
                for i in range(len(dons)):
                    diff = u - D2[i]; nn = max(float(np.linalg.norm(diff)), 1e-9)
                    r.append(nn - mds[i]); J.append(diff / nn)
                J = np.array(J); r = np.array(r)
                try:
                    step = np.linalg.lstsq(J, -r, rcond=None)[0]
                except Exception:
                    break
                u = u + step
                if float(np.linalg.norm(step)) < 1e-9:
                    break
            resid = float(np.linalg.norm(r))
            M = cen + u[0] * e1 + u[1] * e2          # metal IN the donor plane
            coords = P[keep] - M                     # recenter: M -> origin
            if not np.all(np.isfinite(coords)):
                continue
            if best is None or resid < best[0]:
                best = (resid, coords)
        if best is None:
            return None
        lsyms = [mh.GetAtomWithIdx(i).GetSymbol() for i in keep]
        return lsyms, [best[1]]
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


def _collapsed_heavy_bonds_strict(syms, P, factor=0.82):
    """True if any BONDED heavy-heavy non-metal pair sits below ``factor`` × the
    covalent-sum ideal — same logic as ``_has_collapsed_heavy_bonds`` but NOT env-
    gated (always active).  Used to reject the few DG-metallacycle conformers of a
    RIGID PLANAR tridentate that carry a collapsed donor-backbone bond, so a clean
    conformer is selected from the pool.  Universal, geometry-only, deterministic."""
    n = len(syms)
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            ideal = _bd._ideal_bond(syms[i], syms[j])
            if d > 1.30 * ideal:                 # only check pairs that ARE bonded
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
    # RIGID PLANAR tridentate (terpy / pincer, DELFIN_FFFREE_PLANAR_MER): the
    # enumerator seats arm index 1 on the meridian's CENTRAL vertex and arms 0/2 on
    # the outer (antipodal) vertices.  So order the arms [outer, CENTRAL, outer]
    # where the CENTRAL donor = the one lying on the backbone path between the other
    # two (graph-central).  This makes the central pyridyl-N seat on the central
    # vertex -> coplanar meridional placement, outer-outer ~158deg.  Only when the
    # ligand was tagged rigid_planar (flag ON), else byte-identical below.
    if lg.get("rigid_planar") and dent == 3:
        c = _rigid_planar_central_arm(lg["mol"], dons)
        if c is not None:
            outer = [dons[i] for i in range(3) if i != c]
            return [outer[0], dons[c], outer[1]]
    elems = lg.get("donor_elems") or [None] * len(dons)
    elems = list(elems)[:dent]
    try:
        order = sorted(range(dent),
                       key=lambda i: (str(elems[i]) if i < len(elems) else "",
                                      int(dons[i])))
        return [dons[i] for i in order]
    except Exception:
        return dons


def _rigid_planar_central_arm(mol, dons):
    """Index (0/1/2) into ``dons`` of the CENTRAL donor of a rigid planar tridentate
    = the donor that lies ON the backbone shortest path between the OTHER two donors
    (terpy's central pyridyl-N, a pincer's central donor).  Graph-only,
    deterministic; returns None if no single such donor (then the caller falls back
    to the element-sorted order)."""
    try:
        for c in range(3):
            others = [dons[i] for i in range(3) if i != c]
            sp = Chem.GetShortestPath(mol, int(others[0]), int(others[1]))
            if sp and int(dons[c]) in [int(x) for x in sp]:
                return c
    except Exception:
        return None
    return None


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


# --- donor-local VSEPR bend for under-coordinated bent-capable donors ----------
_CHALCOGENS = frozenset(("O", "S", "Se", "Te"))
_PNICTOGENS = frozenset(("N", "P", "As", "Sb"))
# VSEPR ideal M-D-X angles for a SINGLE-substituent donor that keeps its lone
# pairs (chalcogen 2-coord ether/thioether/selenoether/selenolate ~100 deg; a
# pyramidal pnictogen ~107 deg).  CCDC-sane: H2Se 91, R2Se ~96-98, R2Te ~95,
# H2O 104.5, R2O ~111, R3N/R3P ~107.
_DONOR_BEND_DEG = {"O": 109.0, "S": 100.0, "Se": 98.0, "Te": 95.0,
                   "N": 107.0, "P": 100.0, "As": 96.0, "Sb": 95.0}


def _donor_bend_angle(mol, atom):
    """If a SINGLE-ligand-substituent donor ``atom`` (so M + this one substituent
    => 2-coordinate) is a CHALCOGEN or a BENT (pyramidal) PNICTOGEN that retains
    lone pairs, return its VSEPR ideal M-D-X angle in degrees; else ``None`` (=
    keep the linear placement).  Graph/hybridisation-only, no coordinates.

    Genuinely-linear donors are NOT bent: an sp-hybridised nitrogen (nitrile
    N#C, azo/diazo, azide-terminal N), a terminal double-bonded oxo / carbonyl
    O (M=O, M-O#... ), or any donor whose single neighbour is reached by a
    triple bond / allene-type sp centre.  These keep the metal antiperiplanar
    to the substituent (180 deg)."""
    sym = atom.GetSymbol()
    deg = _DONOR_BEND_DEG.get(sym)
    if deg is None:
        return None
    nbrs = list(atom.GetNeighbors())
    if len(nbrs) != 1:
        return None                       # only the M + one-substituent (2-coord) case
    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbrs[0].GetIdx())
    bt = bond.GetBondTypeAsDouble() if bond is not None else 1.0
    # sp donor / multiply-bonded terminal donor => genuinely linear, do not bend.
    hyb = str(atom.GetHybridization())
    if hyb == "SP":
        return None
    if sym in _PNICTOGENS:
        # bend only a *pyramidal* (sp3-ish single-bonded) pnictogen; a
        # double/triple-bonded terminal N (imido/nitrido/diazo) or an aromatic
        # sp2 N stays linear-to-substituent (its lone pair is already the donor
        # axis and bending would distort the multiple bond).
        if bt >= 2.0 or atom.GetIsAromatic():
            return None
        if hyb not in ("SP3", "UNSPECIFIED", "S"):
            return None
    else:  # chalcogen
        # a terminal oxo/chalcogenide double bond (M=O, =S) is linear (the lone
        # pairs sit perpendicular; the donor axis is the pi bond) -> no bend.
        if bt >= 2.0:
            return None
    return float(deg)


def _donor_c_angle(mol, atom):
    """Ideal M-D-R angle (deg) for a SINGLE-heavy-substituent donor ``atom`` whose
    correct local geometry is TETRAHEDRAL or TRIGONAL but which is otherwise built
    LINEAR (180 deg) -- the carbon-donor (and general sp3/sp2 donor) analogue of
    ``_donor_bend_angle``.  Returns ``None`` to keep the linear placement.

    Root cause this addresses: an alkyl / Grignard-type carbanion donor ``M-CH2-R``
    (and ``M-CH3``) loses its donor-carbon hydrogens in the placement graph (the
    fragment is ``[H]C([H])([H])[C]`` with the donor carbon carrying ZERO H and a
    single heavy neighbour), so the donor carbon reaches ``_vsepr_reconstruct`` as a
    k==1 atom and falls through to the linear branch -> a near-linear M-C-C angle
    where VSEPR demands ~109.5 deg.  ``_donor_bend_angle`` only rescues 2-coordinate
    chalcogen / pnictogen donors; this fills the gap for CARBON and any other donor
    whose hybridisation says the metal must sit off the substituent axis.

    Hybridisation-only (graph-derived, no coordinates):
      * sp3 single-substituent donor  -> 109.47 deg (tetrahedral vacancy)
      * sp2 single-substituent donor  -> 120.0  deg (trigonal vacancy)
      * sp  donor                     -> None  (genuinely linear: M-C#O carbonyl,
                                                M-C#N isocyanide, allene/cumulene C,
                                                kept antiperiplanar at 180 deg)
    The genuinely-linear discriminator is HYBRIDISATION (sp), not bond order: a
    kekulised sigma-vinyl donor ``[H]C([H])=[C]`` is sp2 and trigonal (120 deg)
    even though its single substituent is reached by a double bond.  Only an sp3
    donor double/triple-bonded to its substituent is held linear (the double bond
    contradicts sp3 -> ambiguous, keep the conservative 180 deg).  Donors already
    handled by ``_donor_bend_angle`` (chalcogen / pnictogen) return ``None`` here
    so the two helpers never both fire on the same donor."""
    sym = atom.GetSymbol()
    if sym in _CHALCOGENS or sym in _PNICTOGENS:
        return None                       # owned by _donor_bend_angle
    nbrs = list(atom.GetNeighbors())
    if len(nbrs) != 1:
        return None                       # only the M + one-substituent (k==1) case
    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbrs[0].GetIdx())
    bt = bond.GetBondTypeAsDouble() if bond is not None else 1.0
    hyb = str(atom.GetHybridization())
    # The genuinely-linear case is SP hybridisation: a cumulene / vinylidene donor
    # carbon (M=C=CR2), an isocyanide carbon (M-C#N-R), a terminal carbyne.  These
    # keep the metal on the substituent axis (180 deg).
    if hyb == "SP":
        return None
    if hyb == "SP2":
        # trigonal donor (sigma-vinyl/aryl carbanion, sp2 carbene): 120 deg.  A
        # double bond to the *substituent* (kekulised sigma-vinyl [H]C([H])=[C]) is
        # fine -- the donor is still trigonal, only the metal-facing vacancy moves.
        return 120.0
    if hyb in ("SP3", "UNSPECIFIED", "S"):
        # tetrahedral donor (alkyl carbanion).  A double/triple bond to the single
        # substituent contradicts sp3 -> defer to the linear default (do not bend).
        if bt >= 2.0:
            return None
        return 109.47
    return None                           # hypervalent / unknown -> keep linear


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
        # Default: linear (metal antiperiplanar to the single substituent).  With
        # DELFIN_FFFREE_DONOR_BEND=1, a 2-coordinate BENT-CAPABLE donor (chalcogen
        # selenoether/thioether/ether, or a pyramidal pnictogen) gets its real
        # VSEPR M-D-X angle instead of the colinear 180 deg: place the metal
        # vacancy at angle theta from the substituent in an arbitrary (but
        # deterministic) lone-pair plane.  Substituent stays put; the caller's
        # _rot_align(lp, -Vunit) makes M-D-X == theta.  Genuinely-linear donors
        # (sp nitrile/azo N, terminal oxo, =S, M=N) return None -> stay 180 deg.
        if os.environ.get("DELFIN_FFFREE_DONOR_BEND", "0") == "1":
            bend = _donor_bend_angle(lmol, atom)
            if bend is not None:
                s = u[0]                       # donor->substituent unit vector
                # deterministic perpendicular to s (lone-pair plane in-plane axis)
                tmp = (np.array([1.0, 0.0, 0.0]) if abs(s[0]) < 0.9
                       else np.array([0.0, 1.0, 0.0]))
                p = tmp - s * float(np.dot(tmp, s))
                np_ = np.linalg.norm(p)
                if np_ > 1e-6:
                    p = p / np_
                    th = np.radians(bend)
                    # vacancy a with angle(a, s) == theta: cos(theta) along s,
                    # sin(theta) along the perpendicular p.
                    a = np.cos(th) * s + np.sin(th) * p
                    na = np.linalg.norm(a)
                    if na > 1e-6 and np.all(np.isfinite(a)):
                        return lP, a / na
        # DELFIN_FFFREE_DONOR_C_ANGLE=1: the CARBON-donor (and general sp3/sp2
        # single-heavy-substituent donor) analogue of DONOR_BEND.  An alkyl /
        # Grignard carbanion donor M-CH2-R (fragment [H]C([H])([H])[C], donor C
        # with 0 H + 1 heavy neighbour) reaches here as k==1 and would otherwise be
        # placed LINEAR (M-C-C 180 deg) -- VSEPR demands ~109.5 deg (sp3) / 120 deg
        # (sp2).  Offset the metal vacancy off the substituent axis by the
        # hybridisation-ideal angle, IDENTICAL technique to DONOR_BEND above: the
        # substituent subtree stays put on the donor vertex, the caller's
        # _rot_align(lp, -Vunit) makes M-C-R == theta.  Genuinely-linear sp donors
        # (M-C#O carbonyl, M-C#N isocyanide, =C= cumulene) return None -> 180 deg.
        if os.environ.get("DELFIN_FFFREE_DONOR_C_ANGLE", "0") == "1":
            cbend = _donor_c_angle(lmol, atom)
            if cbend is not None:
                s = u[0]                       # donor->substituent unit vector
                tmp = (np.array([1.0, 0.0, 0.0]) if abs(s[0]) < 0.9
                       else np.array([0.0, 1.0, 0.0]))
                p = tmp - s * float(np.dot(tmp, s))
                np_c = np.linalg.norm(p)
                if np_c > 1e-6:
                    p = p / np_c
                    th = np.radians(cbend)
                    a = np.cos(th) * s + np.sin(th) * p
                    na = np.linalg.norm(a)
                    if na > 1e-6 and np.all(np.isfinite(a)):
                        return lP, a / na
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


def _diatomic_orient_enabled() -> bool:
    """DELFIN_FFFREE_DIATOMIC_ORIENT — post-placement orientation guard for linear
    diatomic donors (M-C#O carbonyl, M-C#N cyanide, M-N=O nitrosyl).  Default OFF
    => the guard is never invoked => byte-identical."""
    return os.environ.get("DELFIN_FFFREE_DIATOMIC_ORIENT", "0") == "1"


def _diatomic_donor_partner(lg):
    """For a strictly diatomic (exactly TWO heavy atoms) monodentate ligand, return
    ``(donor_elem, partner_elem)`` where the DONOR is the atom the SMILES bonds to the
    metal (``lg['donor_local_idxs'][0]`` in the ligand graph) and the PARTNER is the
    other heavy atom — but ONLY when the two heavy atoms are DIFFERENT elements (so the
    correct orientation is unambiguous: C#O, C#N, N=O).  Returns ``None`` otherwise
    (not diatomic, polydentate, or homonuclear N#N/etc. where no flip is detectable).

    CONNECTIVITY-ONLY: the donor element comes from the molecular graph, never from a
    hardcoded "C is the donor" rule -> covers C-donor carbonyl/cyanide AND N-donor
    nitrosyl correctly.  Used by the post-placement diatomic-orientation guard."""
    try:
        mol = lg["mol"]
        dons = lg.get("donor_local_idxs", [])
        if int(lg.get("denticity", 0)) != 1 or len(dons) != 1:
            return None
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        if len(heavy) != 2:
            return None                       # not a diatomic
        di = int(dons[0])
        if di not in heavy:
            return None
        partner = heavy[0] if heavy[1] == di else heavy[1]
        de = mol.GetAtomWithIdx(di).GetSymbol()
        pe = mol.GetAtomWithIdx(partner).GetSymbol()
        if de == pe:
            return None                       # homonuclear -> orientation symmetric
        return de, pe
    except Exception:
        return None


def _orient_diatomic_block(Q, lsyms, donor_elem, partner_elem, metal_pos, vertex):
    """Re-orient a placed diatomic ligand block ``Q`` (atoms ordered by ``lsyms``) so
    the SMILES-bonded DONOR atom faces the metal: donor at the coordination vertex,
    partner pointing OUTWARD along the metal->vertex axis.

    The donor / partner atoms are located in the placed block BY ELEMENT (the diatomic
    is heteronuclear, see ``_diatomic_donor_partner``), so the guard is robust even when
    the placed-block atom ordering differs from the ligand-graph ordering (e.g. a shared
    conformer-cache entry built from another instance of the same ligand type — the
    actual root cause of the M-O-C isocarbonyl flip).

    No-op (returns ``Q`` unchanged) when the donor is ALREADY closer to the metal than
    the partner.  Rigid: the donor-partner bond length is preserved exactly.  Pure
    geometry, deterministic, never raises (returns ``Q`` on any failure)."""
    try:
        di = [i for i, s in enumerate(lsyms) if s == donor_elem]
        pi = [i for i, s in enumerate(lsyms) if s == partner_elem]
        if len(di) != 1 or len(pi) != 1:
            return Q                          # ambiguous (extra atoms) -> leave as built
        di, pi = di[0], pi[0]
        d_pos = np.asarray(Q[di], float)
        p_pos = np.asarray(Q[pi], float)
        m = np.asarray(metal_pos, float)
        if not (np.all(np.isfinite(d_pos)) and np.all(np.isfinite(p_pos))):
            return Q
        d_md = float(np.linalg.norm(d_pos - m))
        p_md = float(np.linalg.norm(p_pos - m))
        if d_md <= p_md:
            return Q                          # already donor-bound -> byte-identical
        # flipped (partner closer to metal): rebuild the rigid 2-atom unit with the
        # donor at the vertex and the partner outward along the metal->vertex axis.
        bond = float(np.linalg.norm(p_pos - d_pos))
        if bond < 1e-6:
            return Q
        vtx = np.asarray(vertex, float)
        axis = vtx - m
        na = float(np.linalg.norm(axis))
        if na < 1e-6:
            return Q
        out = axis / na                       # metal -> vertex = outward direction
        newQ = np.array(Q, float)
        newQ[di] = vtx                         # donor seats on the vertex
        newQ[pi] = vtx + out * bond            # partner points outward, bond preserved
        if not np.all(np.isfinite(newQ)):
            return Q
        return newQ
    except Exception:
        return Q


def assemble_monodentate(metal: str, ligand_smiles: str, donor_idx: int,
                         geometry: str, thetas=None) -> Tuple[List[str], np.ndarray]:
    """Assemble; thetas[i] = rotation of ligand i about its own M-D axis (free DOF
    used for clash relief — geometric, metal-FF-free)."""
    ref = MSB._ref_vectors(geometry)
    n = len(ref)
    lsyms, lP, lmol = _ligand_3d(ligand_smiles)
    donor_elem = lsyms[donor_idx]
    md = MSB.md_distance(metal, donor_elem,
                         atom=lmol.GetAtomWithIdx(donor_idx), mol=lmol)
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
        T1 = ref[verts[0]] / np.linalg.norm(ref[verts[0]]) * MSB.md_distance(
            metal, lsyms[d1], atom=lmol.GetAtomWithIdx(d1), mol=lmol)
        T2 = ref[verts[1]] / np.linalg.norm(ref[verts[1]]) * MSB.md_distance(
            metal, lsyms[d2], atom=lmol.GetAtomWithIdx(d2), mol=lmol)
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
    md1 = MSB.md_distance(metal, lsyms[d1], atom=lmol.GetAtomWithIdx(d1), mol=lmol)
    md2 = MSB.md_distance(metal, lsyms[d2], atom=lmol.GetAtomWithIdx(d2), mol=lmol)
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


# Process-local conformer memo (deterministic embed -> cacheable).  Keyed by the
# fragment's canonical SMILES + k, so the ENSEMBLE builder re-uses one ligand embed
# across all its variant builds instead of re-running ETKDG per variant (the embed of
# a large flexible η-substituent dominates the per-build cost).  Transparent: the same
# input always returns the same (deterministic) result, so behaviour is unchanged —
# this is a speedup, not a semantic change.  Bounded to avoid unbounded growth.
_CONF_CACHE = {}
_CONF_CACHE_MAX = 256


def _ligand_confs_from_mol(frag_mol, k=10):
    """UNIVERSAL multi-conformer generation for a ligand (deterministic): K diverse
    ETKDG conformers (fixed seed, single-thread) + MMFF.  Returns (syms, [coords],
    mol).  Used to pick the clash-minimal conformer per ligand at placement — a
    fundamental Layer-3 mechanism applied to every ligand, not a per-case patch.

    Result is memoised by canonical SMILES + k (deterministic embed): repeated calls
    for the same ligand (e.g. the ensemble builder's variant loop) skip the costly
    re-embed.  The cached coords/mol are NOT mutated by any caller."""
    key = None
    try:
        key = (Chem.MolToSmiles(frag_mol), int(k))
    except Exception:
        key = None
    if key is not None and key in _CONF_CACHE:
        return _CONF_CACHE[key]
    m = Chem.AddHs(frag_mol)
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                           numThreads=1))
    if not cids:
        if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
            if key is not None:
                _CONF_CACHE[key] = None
            return None
        cids = [0]
    try:
        AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=1)
    except Exception:
        pass
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    out = (syms, [np.array(m.GetConformer(c).GetPositions(), float) for c in cids], m)
    if key is not None and len(_CONF_CACHE) < _CONF_CACHE_MAX:
        _CONF_CACHE[key] = out
    return out


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


def _ligand_block_bonds(lmol, offset, donor_local):
    """True connectivity of one ligand block for the #308 torsion relaxer.

    Returns ``(offset, [(li, lj), ...], donor_local)`` where the local (i,j) bond
    pairs come directly from the ligand mol (atom order preserved through assembly,
    AddHs included), so the relaxer never has to GUESS bonds from distance on a
    crowded complex (where two ligands at a fortuitous bonding distance would be
    mis-read as covalently bonded).  Returns ``None`` on any failure (relaxer then
    falls back to geometric perception)."""
    try:
        lb = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in lmol.GetBonds()]
        return (int(offset), lb, int(donor_local))
    except Exception:
        return None


def _torsion_relax_frame(out_syms, P, fixed, block_specs):
    """Apply the env-gated #308 torsion-space clash relax to one assembled frame,
    threading the true per-ligand connectivity (``block_specs`` = list of
    ``(offset, lmol, donor_local)``).  No-op when the flag is unset; never raises."""
    try:
        from delfin.fffree import torsion_relax as _TR
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _TR.bonds_from_blocks(0, blocks)
        return np.asarray(_TR.relax_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P


def _joint_declash_frame(out_syms, P, fixed, block_specs):
    """Apply the env-gated JOINT global INTER-LIGAND heavy-heavy declash to one
    assembled frame (``DELFIN_FFFREE_JOINT_DECLASH``), threading the true
    per-ligand connectivity (``block_specs`` = list of ``(offset, lmol,
    donor_local)``).  Runs AFTER #308 torsion-relax and BEFORE the self-gate so a
    declashed class-B build passes ``_build_is_clean``.  No-op when the flag is
    unset; never raises."""
    try:
        from delfin.fffree import joint_declash as _JD
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _JD._TR.bonds_from_blocks(0, blocks)
        return np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P


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
    block_specs = []                  # (offset, lmol, donor_local) for the torsion relax
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag)
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
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
                # diatomic-orientation guard (env-gated, default OFF -> byte-identical):
                # keep the SMILES-bonded donor atom (not its partner) at the vertex.
                if _diatomic_orient_enabled():
                    _dp = _diatomic_donor_partner(
                        {"mol": frag, "donor_local_idxs": [di], "denticity": 1})
                    if _dp is not None:
                        Q = _orient_diatomic_block(
                            Q, lsyms, _dp[0], _dp[1], np.zeros(3), vertex)
            cl = _clash_count(Q, np.array(placed_P), lsyms, placed_syms)
            if cl < best_clash:
                best_clash, best_Q = cl, Q
            if cl == 0:
                break
        block_specs.append((pos, lmol, di))
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
        # #308 whole-complex torsion-space clash relax (env-gated, default-OFF
        # byte-id): when rigid M-D-axis selection is not enough and ligand-internal
        # rotation is needed, jointly optimise all rotatable single bonds of the
        # assembled complex (metal + donors = `fixed` frozen).  Torsion-only -> bond
        # lengths/angles preserved exactly; never-worse.  No-op when flag unset.
        P = _torsion_relax_frame(out_syms, P, fixed, block_specs)
        # JOINT inter-ligand declash (env-gated, default-OFF byte-id): whole-ligand
        # M-D-axis rotations + internal torsions jointly minimising the GLOBAL
        # inter-ligand heavy-heavy clash (the self-gate blocker for class-B), core
        # frozen.  Runs after #308; no-op when DELFIN_FFFREE_JOINT_DECLASH unset.
        P = _joint_declash_frame(out_syms, P, fixed, block_specs)
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
    block_specs = []          # (offset, lmol, donor_local) for the #308 torsion relax
    pos = 1
    metal_P = np.zeros((1, 3))
    metal_sym = [metal]
    for i, (frag, di) in enumerate(vertex_specs):
        Vunit = ref[i] / np.linalg.norm(ref[i])
        confs = _ligand_confs_from_mol(frag, k=max(per_lig_confs, 1))
        if confs is None:
            return None
        lsyms, coords_list, lmol = confs
        block_specs.append((pos, lmol, di))
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
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
                # diatomic-orientation guard (env-gated, default OFF -> byte-identical)
                if _diatomic_orient_enabled():
                    _dp = _diatomic_donor_partner(
                        {"mol": frag, "donor_local_idxs": [di], "denticity": 1})
                    if _dp is not None:
                        Q = _orient_diatomic_block(
                            Q, lsyms, _dp[0], _dp[1], np.zeros(3), vertex)
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

    # ---- inter-ligand-clash-aware selection (env-gated, default OFF = byte-id) ----
    # ROOT FIX (#306 inter-ligand): the per-vertex candidate scoring above clashes vs
    # the METAL ONLY, and this combo loop orders by metal-clash RANK SUM and NEVER
    # scores the assembled frame for inter-ligand overlap.  For homoleptic bulky
    # ligands every vertex then independently picks the same metal-optimal pose ->
    # systematic inter-ligand clash, identical in every frame (e.g. Zr(C6H5)6: min
    # inter-aryl H-H = 0.73 A in all frames).  The axial-spin candidates that COULD
    # interleave the ligands already exist in per_vertex_cands -- they are just never
    # selected for.  When DELFIN_FFFREE_INTERLIG_RANK=1: (a) lead frame 0 with a
    # GREEDY sequential placement that picks, at each vertex, the candidate minimising
    # _clash_count vs the already-placed ligands (mirrors assemble_heteroleptic_from_
    # mols), and (b) emit the remaining frames in ASCENDING total inter-ligand clash.
    # RIGID only: we re-select among existing conformer/spin candidates; the core
    # (vertex directions + M-D distances) stays fixed.  Deterministic, geometry-only.
    _interlig = os.environ.get("DELFIN_FFFREE_INTERLIG_RANK", "0") == "1"
    MAX_EVAL = 64                                    # bound the build work

    def _interlig_clash(cb):
        """Total inter-ligand _clash_count for combo cb: sum over each ligand block
        of its clash vs the union of PREVIOUSLY-placed ligand blocks (metal at origin
        excluded; same H-inclusive measure used everywhere).  Rigid-block geometry."""
        placed = []                       # accumulated already-placed ligand atoms
        placed_syms = []
        total = 0
        for vi, ci in enumerate(cb):
            Q = per_vertex_cands[vi][ci][0]
            lsyms = vertex_lsyms[vi]
            if placed:
                total += _clash_count(Q, np.array(placed), lsyms, placed_syms)
            for row in Q:
                placed.append(row)
            placed_syms += lsyms
        return total

    if _interlig:
        MAX_EVAL = 128                               # modest, capped bump when ON
        # (a) greedy inter-ligand-minimal lead combo: fixed deterministic vertex order,
        #     each vertex takes the candidate minimising clash vs already-placed blocks.
        greedy = []
        gplaced = []; gplaced_syms = []
        for vi in range(len(per_vertex_cands)):
            cands_vi = per_vertex_cands[vi]
            lsyms = vertex_lsyms[vi]
            best_ci, best_cl = 0, None
            for ci, (Q, _mc) in enumerate(cands_vi):
                if not gplaced:
                    cl = 0
                else:
                    cl = _clash_count(Q, np.array(gplaced), lsyms, gplaced_syms)
                # strictly-less keeps the FIRST (lowest-index = lowest metal-clash)
                # candidate on ties -> deterministic, and =all-best when no clash.
                if best_cl is None or cl < best_cl:
                    best_cl, best_ci = cl, ci
                if cl == 0:
                    break
            greedy.append(best_ci)
            Qsel = cands_vi[best_ci][0]
            for row in Qsel:
                gplaced.append(row)
            gplaced_syms += lsyms
        greedy = tuple(greedy)
        # rank the bounded combo window by total inter-ligand clash (ascending), then
        # by the original metal-rank ordering as a stable deterministic tiebreak; put
        # the greedy lead combo first (dedup so it is not double-evaluated).
        window = combos[:MAX_EVAL]
        if greedy not in window:
            window = [greedy] + window[:MAX_EVAL - 1]
        scored = [(_interlig_clash(cb), sum(cb), cb) for cb in window]
        # greedy combo forced to the front via a sentinel key; rest by (clash, rank).
        scored.sort(key=lambda t: (0 if t[2] == greedy else 1, t[0], t[1], t[2]))
        eval_order = [t[2] for t in scored]
    else:
        eval_order = combos[:MAX_EVAL]

    frames = []                                      # (syms, P) kept (deduped)
    for cb in eval_order:
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
            # #308 whole-complex torsion-space clash relax (env-gated, default-OFF
            # byte-id); torsion-only, never-worse, metal+donors (`fixed`) frozen.
            P = _torsion_relax_frame(out_syms, P, fixed, block_specs)
            # JOINT inter-ligand declash (env-gated, default-OFF byte-id): global
            # inter-ligand heavy-heavy minimisation, core frozen.  After #308.
            P = _joint_declash_frame(out_syms, P, fixed, block_specs)
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
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
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


def _eta_centroid_distance(metal, eta_n):
    """Crystallographic metal→ring-centroid distance (Å) for an η-face.  Pulls the
    open-source literature-averaged table from smiles_converter (a hard-coded table
    of published averages; reads no proprietary database file) and falls back to its
    geometric estimator.  Deterministic; always finite."""
    try:
        from delfin.smiles_converter import _target_mc_dist
        d = float(_target_mc_dist(metal, int(eta_n)))
        if np.isfinite(d) and d > 0.5:
            return d
    except Exception:
        pass
    # last-resort: covalent M-C sum minus a ring-radius slip (always finite)
    rmc = PLY.COV.get(metal, 1.5) + PLY.COV.get("C", 0.76)
    if eta_n >= 3:
        rr = 1.40 / (2.0 * np.sin(np.pi / max(eta_n, 3)))
        v = rmc * rmc - rr * rr
        return float(np.sqrt(v)) if v > 0.25 else 0.80 * rmc
    return 0.85 * rmc


def _rot_about_axis(axis, theta):
    """Rodrigues rotation matrix: rotate by ``theta`` (rad) about unit ``axis``."""
    axis = np.asarray(axis, float)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        return np.eye(3)
    axis = axis / n
    c, s = float(np.cos(theta)), float(np.sin(theta))
    x, y, z = axis
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) * c + s * K + (1.0 - c) * np.outer(axis, axis)


def _piano_leg_tilt_enabled() -> bool:
    """Piano-stool σ-leg TILT correction (default OFF -> byte-identical).

    A half-sandwich / piano-stool (ONE η-face + a small set of σ-legs) is NOT a
    regular tetrahedron: the η-ring is a fat 3-electron-pair FACE donor, and the legs
    splay DOWN, away from the ring, so the real (ring-centroid)-M-L angle is ~120-130
    deg (not the ~90-110 deg a generic CN4/CN5 polyhedron places the legs at).  When
    ON, the legs are rigidly tilted to that piano-stool angle; OFF, the block is
    skipped and the build is byte-identical to the historic polyhedron placement."""
    return os.environ.get("DELFIN_FFFREE_PIANO_LEG_TILT", "0") == "1"


def _piano_leg_target_deg() -> float:
    """Target (η-ring-centroid)-M-(σ-leg) tilt angle in degrees for the piano-stool
    leg-tilt correction.  Default 125 deg -- the canonical crystallographic value for
    half-sandwich piano-stools: (arene)Cr(CO)3 centroid-Cr-C ~125.7 deg, CpMn(CO)3
    centroid-Mn-C ~125 deg, CpFe(CO)2X ~121-126 deg.  Override via
    DELFIN_FFFREE_PIANO_LEG_DEG; the legs splay DOWN (away from the ring), so values
    >90 deg point the legs to the FAR hemisphere from the centroid axis."""
    try:
        v = float(os.environ.get("DELFIN_FFFREE_PIANO_LEG_DEG", "125.0"))
    except (TypeError, ValueError):
        return 125.0
    # clamp to a sane piano-stool window so a bad env value can never flip the tripod
    return min(160.0, max(95.0, v))


def _tilt_piano_legs(P, prim_V, sigma_blocks, target_deg):
    """Rigidly tilt each σ-leg block to the piano-stool (ring-centroid)-M-L angle.

    ``P`` is the assembled coordinate array with the metal at the origin (row 0).
    ``prim_V`` is the UNIT M->(primary η-ring centroid) direction.  ``sigma_blocks``
    is the list of (start, end) atom-row half-open ranges of the σ-co-ligand blocks
    (the legs).  ``target_deg`` is the desired centroid-M-L angle.

    For each leg, the representative donor (the block atom CLOSEST to the metal) gives
    the leg's current M->L direction.  The whole leg block is rotated RIGIDLY about an
    axis through the metal (the origin) that is PERPENDICULAR to the plane spanned by
    ``prim_V`` and the leg direction, by exactly the angle needed to bring the
    centroid-M-L angle to ``target_deg``.  Because the rotation passes through the
    metal it preserves every atom's metal distance (M-L INVARIANT) and the leg's
    internal geometry; it only swings the leg DOWN, away from the η-ring.  Each leg
    keeps its own azimuth about ``prim_V`` (the rotation axis is radial-tangential),
    so a 3-fold-symmetric tripod stays 3-fold symmetric.  The η-ring is never touched.

    Universal over leg count (1..n) and ring type; geometry-only, deterministic,
    never raises (returns ``P`` unchanged on any degeneracy).  Returns the modified
    array (in-place edits applied to ``P``)."""
    prim_V = np.asarray(prim_V, float)
    nV = np.linalg.norm(prim_V)
    if nV < 1e-9 or not sigma_blocks:
        return P
    prim_V = prim_V / nV
    tgt = np.radians(float(target_deg))
    M = P[0]
    for (s, e) in sigma_blocks:
        if e <= s:
            continue
        blk = P[s:e]
        # representative leg donor = block atom nearest the metal (the σ-donor).
        rel = blk - M
        dists = np.linalg.norm(rel, axis=1)
        di = int(np.argmin(dists))
        L = rel[di]
        nL = np.linalg.norm(L)
        if nL < 1e-9:
            continue
        Lu = L / nL
        cur = float(np.arccos(np.clip(np.dot(prim_V, Lu), -1.0, 1.0)))
        delta = tgt - cur
        if abs(delta) < 1e-6:
            continue
        # rotation axis = prim_V x leg-direction (perpendicular to the swing plane).
        # Sign of `delta` then moves the leg toward/away from the centroid axis; a
        # positive delta (target > current) opens the angle -> legs splay DOWN.
        ax = np.cross(prim_V, Lu)
        nax = np.linalg.norm(ax)
        if nax < 1e-8:
            # leg is (anti)parallel to the centroid axis -> tilt direction undefined;
            # use any axis perpendicular to prim_V (deterministic) so it still splays.
            ref = np.array([1.0, 0.0, 0.0]) if abs(prim_V[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ax = np.cross(prim_V, ref)
            nax = np.linalg.norm(ax)
            if nax < 1e-8:
                continue
        ax = ax / nax
        R = _rot_about_axis(ax, delta)
        # rotate the WHOLE leg block about the metal at the origin (no recentre).
        P[s:e] = blk @ R.T
    return P


def _place_eta_ring(metal, lsyms, lP, eta_idxs, Vunit, mc_dist,
                    spin=0.0, pucker=0.0):
    """Rigidly seat an η-face so its ring CENTROID sits on the polyhedron vertex
    direction ``Vunit`` at distance ``mc_dist`` from the metal (origin), with the
    ring plane PERPENDICULAR to the M→centroid axis.  The whole ligand (ring +
    substituents) is moved as one rigid body, so substituents are dragged with the
    ring and the η-carbons sit at the ring radius AROUND the centroid (never on the
    metal).  Returns the transformed coords, or None on degeneracy.

    Geometry-only, deterministic.  The in-plane spin is fixed canonically (first
    ring atom placed at azimuth 0) so two builds of the same input are identical.

    ``spin`` (rad): additional rigid rotation of the whole ligand ABOUT the
    M→centroid axis — generates the η-ring ROTAMER orientations (rotates substituent
    positions around the ring while the ring itself maps onto itself up to symmetry).
    ``pucker`` (Å): out-of-plane displacement amplitude of alternating ring atoms
    (Cremer-Pople-style envelope/twist mode) applied to NON-aromatic η-rings before
    seating, so diene / allyl / cyclohexadienyl rings emit puckered conformers.  Both
    default to 0.0 -> byte-identical to the canonical single rigid build."""
    lP = np.asarray(lP, float).copy()
    ring = [int(i) for i in eta_idxs]
    if len(ring) < 2:
        return None
    R = lP[ring]
    cen = R.mean(axis=0)                                   # ring centroid (ligand frame)
    X = R - cen
    # ring-plane normal via SVD (smallest singular direction); robust for any n>=2.
    try:
        _, _, Vt = np.linalg.svd(X)
    except Exception:
        return None
    if len(ring) == 2:
        # η2 (alkene): "plane normal" is ambiguous; use the C=C bond direction as the
        # in-plane axis and pick a normal perpendicular to it (deterministic).
        nrm = Vt[2] if Vt.shape[0] > 2 else np.cross(Vt[0], np.array([0.0, 0.0, 1.0]))
    else:
        nrm = Vt[2]
    nn = np.linalg.norm(nrm)
    if nn < 1e-8:
        return None
    nrm = nrm / nn
    # Cremer-Pople-style pucker: displace ring atoms along the ring NORMAL with an
    # alternating sign (B/T-like envelope), in the LIGAND frame, before seating.  Only
    # the ring atoms move; substituents stay rigid relative to the (unpuckered) ring
    # carbon they hang off — this is a deterministic conformer, not an MMFF relax.
    if abs(pucker) > 1e-9 and len(ring) >= 4:
        for k, ai in enumerate(ring):
            lP[ai] = lP[ai] + nrm * (float(pucker) * (1.0 if (k % 2 == 0) else -1.0))
        # recompute centroid after the symmetric in-out displacement (≈ unchanged)
        cen = lP[ring].mean(axis=0)
    Vunit = np.asarray(Vunit, float)
    Vunit = Vunit / np.linalg.norm(Vunit)
    target_centroid = Vunit * float(mc_dist)
    # Rotate so the ring-plane normal aligns with the M→centroid axis: the metal
    # then sits along the ring normal (face-on coordination), exactly as in a real
    # piano-stool / metallocene.  Sign chosen so the ring is on the FAR side of the
    # metal (centroid points away from origin along +Vunit).
    Rrot = _rot_align(nrm, Vunit)
    Q = (lP - cen) @ Rrot.T + target_centroid
    if abs(spin) > 1e-9:
        # rotate the seated ligand about the M→centroid axis through the centroid:
        # spins substituents to a distinct rotamer; the bare ring maps onto itself.
        Spin = _rot_about_axis(Vunit, float(spin))
        Q = (Q - target_centroid) @ Spin.T + target_centroid
    if not np.all(np.isfinite(Q)):
        return None
    return Q


def assemble_hapto(metal, geometry, d, variant=None):
    """Build a hapto complex on the FF-free path: η-faces are placed as RIGID rings
    on their centroid vertices; σ-donors are oriented onto the remaining vertices
    (homoleptic/heteroleptic monodentate).  ``d`` is a rigid-hapto decompose dict
    (has per-ligand 'is_eta'/'eta_local_idxs'/'eta_n').  v1: chelating σ-arms are
    placed by the existing rigid chelate block.  Returns (syms, P, donors) where
    ``donors`` are the global indices fffree treats as the coordination shell (one
    representative atom per η-face + every σ-donor), or None on failure.

    Universal, geometry-only, deterministic (fixed seeds, canonical in-plane spin).

    Returns the tuple (syms, P, donors, exempt_pairs) — exempt_pairs are the global
    index pairs of GENUINE multiple bonds (C≡O carbonyls, C≡N nitriles, C=C, …) whose
    short length is chemically correct, so the self-gate does not mistake them for a
    collapse.

    ``variant`` (optional dict, default None = byte-identical canonical build): a
    per-build perturbation spec used by ``assemble_hapto_ensemble`` to enumerate the
    rigid ENSEMBLE.  Recognised keys::

        variant["eta"][li] = {"spin": rad, "pucker": Å, "eta_idxs": [...], "eta_n": k}

    where ``spin`` rotates the seated η-ligand about the M→centroid axis (η-ring
    ROTAMER), ``pucker`` applies a Cremer-Pople envelope to non-aromatic rings, and
    ``eta_idxs``/``eta_n`` (optional) override the ring-atom set + hapticity for a
    deterministic ring-SLIP (η6→η4→η2 / η5→η3→η1) isomer.  All variation is geometric
    + deterministic; an absent/None variant reproduces the historic single build."""
    variant = variant or {}
    eta_var = variant.get("eta", {})
    # Hebel A (DELFIN_FFFREE_HAPTO_AXIS_ROT): rigid rotation (rad) of every NON-
    # primary-η atom block (the CO-tripod / co-ligands / a 2nd ring) about the
    # PRIMARY η-face M→centroid axis.  The metal sits at the origin, so any rotation
    # about an axis through the origin preserves every atom's distance to the metal
    # (M-centroid + all M-D distances INVARIANT); only the ring-vs-rest azimuthal
    # CLOCK changes — the η-ring internal geometry is untouched.  0.0 = byte-identical
    # to the historic build (the rotation is skipped entirely).
    _axis_rot = float(variant.get("axis_rot", 0.0) or 0.0)
    ref = MSB._ref_vectors(geometry)
    n_vert = len(ref)
    ligands = d["ligands"]
    # how many vertices each ligand occupies (η-face = 1, σ = denticity)
    occ = [(1 if lg.get("is_eta") else lg["denticity"]) for lg in ligands]
    if sum(occ) != n_vert:
        return None
    # Deterministic vertex assignment: η-faces first (lowest vertex indices), then
    # σ-donors in ligand order.  Simple + reproducible; isomer enumeration is a
    # later increment (v1 emits ONE faithful build per hapto complex).
    order = sorted(range(len(ligands)), key=lambda i: (0 if ligands[i].get("is_eta") else 1, i))
    out_syms = [metal]
    placed = [np.zeros(3)]
    placed_syms = [metal]
    donors = []
    relax_frags = []
    exempt_pairs = []        # global index pairs that are genuine double/triple bonds
    pos = 1
    vi = 0
    fixed = {0}
    # Hebel-A bookkeeping: per-η-face (atom-block-start, atom-block-end, centroid-unit-
    # direction); and the contiguous atom blocks of every NON-η ligand (CO/σ co-ligands).
    eta_blocks = []          # [(start, end, Vunit), ...] in emission (= placement) order
    nonprimary_eta_blocks = []   # (start, end) of η faces 2..n  -> rotated with the rest
    sigma_blocks = []        # (start, end) of σ ligands

    def _collect_exempt(frag_mol, lig_offset):
        """Local heavy-atom double/triple bonds -> global index pairs (the AddHs
        ligand block starts at lig_offset+1 in the assembled coords)."""
        for b in frag_mol.GetBonds():
            bt = b.GetBondTypeAsDouble()
            if bt >= 1.5:                       # aromatic(1.5)/double(2)/triple(3)
                a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                if (frag_mol.GetAtomWithIdx(a1).GetAtomicNum() > 1
                        and frag_mol.GetAtomWithIdx(a2).GetAtomicNum() > 1):
                    g1 = lig_offset + 1 + a1
                    g2 = lig_offset + 1 + a2
                    exempt_pairs.append((min(g1, g2), max(g1, g2)))
    for li in order:
        lg = ligands[li]
        lig_offset = pos - 1
        if lg.get("is_eta"):
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            Vunit = ref[vi] / np.linalg.norm(ref[vi])
            # per-build variation (ensemble): η-ring spin/pucker + optional ring-slip.
            _ev = eta_var.get(li, {})
            _spin = float(_ev.get("spin", 0.0))
            _pucker = float(_ev.get("pucker", 0.0))
            _eta_idxs = _ev.get("eta_idxs", lg["eta_local_idxs"])
            _eta_n = int(_ev.get("eta_n", lg["eta_n"]))
            mc = _eta_centroid_distance(metal, _eta_n)
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                Q = _place_eta_ring(metal, lsyms, lP, _eta_idxs, Vunit, mc,
                                    spin=_spin, pucker=_pucker)
                if Q is None:
                    continue
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            # representative donor = ring atom closest to the metal (for the self-gate).
            # Uses the ACTIVE ring set (_eta_idxs) so ring-slip isomers freeze + report
            # the η-carbons they actually coordinate, not the full pre-slip face.
            ring_globals = [lig_offset + 1 + int(j) for j in _eta_idxs]
            rep = min(ring_globals, key=lambda g: float(np.linalg.norm(placed[g])))
            donors.append(rep)
            fixed.update(lig_offset + 1 + int(j) for j in _eta_idxs)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            _blk = (lig_offset + 1, lig_offset + 1 + len(lsyms))
            eta_blocks.append((_blk[0], _blk[1], Vunit.copy()))
            if len(eta_blocks) > 1:               # secondary η face -> rotates with rest
                nonprimary_eta_blocks.append(_blk)
            pos += len(lsyms)
            vi += 1
        elif lg["denticity"] == 1:
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            di = lg["donor_local_idxs"][0]
            Vunit = ref[vi] / np.linalg.norm(ref[vi])
            md = MSB.md_distance(metal, lsyms[di],
                                 atom=lmol.GetAtomWithIdx(di), mol=lmol)
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                if len(lsyms) == 1:
                    Q = (Vunit * md).reshape(1, 3)
                else:
                    lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, di)
                    Q = (lPv - lPv[di]) @ _rot_align(lp, -Vunit).T + Vunit * md
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            donors.append(lig_offset + 1 + di)
            fixed.add(lig_offset + 1 + di)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            sigma_blocks.append((lig_offset + 1, lig_offset + 1 + len(lsyms)))
            pos += len(lsyms)
            vi += 1
        else:
            # chelating σ-arm: rigid-fit onto the next `denticity` vertices.
            dent = lg["denticity"]
            confs = _ligand_confs_from_mol(lg["mol"])
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
            dons_d = _canonical_arm_order(lg, dent)
            verts = list(range(vi, vi + dent))
            targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                       * MSB.md_distance(metal, lsyms[dons_d[i]],
                                         atom=lmol.GetAtomWithIdx(dons_d[i]), mol=lmol)
                       for i in range(dent)]
            best_Q, best_clash = None, 1e18
            for lP in coords_list:
                if dent == 2:
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1],
                                             targets[0], targets[1])
                else:
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=True)
                if Q is None:
                    continue
                cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
                if cl < best_clash:
                    best_clash, best_Q = cl, Q
                if cl == 0:
                    break
            if best_Q is None:
                return None
            out_syms += lsyms
            for row in best_Q:
                placed.append(row)
            placed_syms += lsyms
            for di in dons_d:
                donors.append(lig_offset + 1 + di)
                fixed.add(lig_offset + 1 + di)
            relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))
            _collect_exempt(lg["mol"], lig_offset)
            sigma_blocks.append((lig_offset + 1, lig_offset + 1 + len(lsyms)))
            pos += len(lsyms)
            vi += dent
    P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
    # --- Hebel A: rigid η-axis rotation of the non-primary-η atoms --------------
    # Rotate the CO-tripod / co-ligands / a secondary ring as ONE rigid body about
    # the PRIMARY η-face M→centroid axis (through the metal at the origin).  This is
    # the missing relative DOF for piano-stool / sandwich systems: the η-ring sits
    # right but the tripod/2nd-ring CLOCK is not yet at the crystal rotamer.  Because
    # the axis passes through the metal, every atom's metal distance is preserved
    # (M-centroid + M-D INVARIANT); only the azimuth changes.  Strictly additive: the
    # canonical build (axis_rot==0) is byte-identical (this block is skipped).
    if abs(_axis_rot) > 1e-9 and eta_blocks:
        prim_V = np.asarray(eta_blocks[0][2], float)
        nV = np.linalg.norm(prim_V)
        if nV > 1e-9:
            prim_V = prim_V / nV
            Rax = _rot_about_axis(prim_V, _axis_rot)
            rot_blocks = list(sigma_blocks) + list(nonprimary_eta_blocks)
            for (s, e) in rot_blocks:
                P[s:e] = P[s:e] @ Rax.T            # metal at origin -> no recentre
            if not np.all(np.isfinite(P)):
                return None
    # --- Piano-stool σ-leg TILT (DELFIN_FFFREE_PIANO_LEG_TILT, default OFF) -------
    # A half-sandwich (EXACTLY one η-face + σ-legs, no second ring) is not a regular
    # polyhedron: the η-ring is a fat FACE donor and the legs splay DOWN, so the real
    # (ring-centroid)-M-L angle is ~125 deg, not the ~110 deg (T-4) / ~90 deg (SP-4)
    # the generic CNn polyhedron places them at.  Rigidly tilt each σ-leg block about
    # the metal so that angle is correct; the η-ring + every M-L distance are invariant.
    # Skipped (byte-identical) when the flag is off, when there is not exactly one
    # η-face (full sandwich / metallocene untouched), or when there are no σ-legs.
    if (_piano_leg_tilt_enabled() and len(eta_blocks) == 1 and sigma_blocks):
        try:
            prim_V = np.asarray(eta_blocks[0][2], float)
            P = _tilt_piano_legs(P, prim_V, sigma_blocks, _piano_leg_target_deg())
        except Exception:
            pass
        if not np.all(np.isfinite(P)):
            return None
    # FF-free geometric clash-relief (η-ring + σ-donors all frozen so the rigid
    # ring + constructed coordination are preserved; periphery relaxes only).
    try:
        from delfin.fffree.refine import refine as _refine
        P = _refine(out_syms, P, fixed)
    except Exception:
        pass
    if not np.all(np.isfinite(P)):
        return None
    return out_syms, P, sorted(set(donors)), exempt_pairs


# ---------------------------------------------------------------------------
# Rigid-hapto ENSEMBLE (deterministic, RMSD-deduplicated) — used by the FF-free
# converter backend when DELFIN_FFFREE_RIGID_HAPTO=1.  Every member is a fully
# rigid build (no ring collapse); the variation is purely geometric + enumerated
# deterministically, so two runs of the same SMILES emit byte-identical frames.
# ---------------------------------------------------------------------------

# Symmetry-distinct in-plane spin counts per hapticity for a BARE ring (the ring
# carbons map onto themselves under the cyclic group C_n; a SUBSTITUTED ring breaks
# that symmetry, so spinning to n_fold orientations samples the substituent rotamers
# without emitting bare-ring duplicates — the RMSD dedup removes any that coincide).
_ETA_NFOLD = {2: 2, 3: 3, 4: 4, 5: 5, 6: 6}


def _ring_spins(eta_n):
    """Deterministic, symmetry-reduced list of η-ring spin angles (rad), starting
    at the canonical 0.  For an η-n face the carbons recur every 2π/n; we sample one
    full inter-carbon sector in n_fold steps (so a substituted ring visits each
    distinct substituent-vs-coligand clock position exactly once).  Duplicate frames
    (symmetric rings) are removed later by the RMSD dedup."""
    nf = _ETA_NFOLD.get(int(eta_n), max(2, int(eta_n)))
    sector = 2.0 * np.pi / float(nf)
    # sample the sector at nf sub-steps -> the substituent sweeps a full inter-carbon
    # gap; canonical 0 first so member[0] is the historic single build.
    return [sector * (k / float(nf)) for k in range(nf)]


def _axis_rot_angles(n_rot, eta_n):
    """Deterministic η-axis rotation angles (rad) for Hebel A: rotate the non-η
    block (CO-tripod / co-ligands / 2nd ring) about the primary η M→centroid axis.

    The canonical 0 is NOT included here (the base ensemble already carries it as
    axis_rot==0); these are the ADDED frames.  We sweep the full 0–2π in ``n_rot``
    equal steps and drop the 0 step, so the extra frames are a uniform azimuthal
    sweep of the ring-vs-rest clock.  Bare-ring / symmetric duplicates are removed
    downstream by the RMSD dedup; the sweep is intentionally NOT symmetry-reduced
    here because the relevant symmetry is the PRODUCT of the ring C_n and the tripod
    C_m (system-dependent), and over-reducing would skip the crystal rotamer."""
    n_rot = max(2, int(n_rot))
    step = 2.0 * np.pi / float(n_rot)
    return [step * k for k in range(1, n_rot)]


def _slip_modes(lg):
    """Deterministic ring-SLIP isomers for one η-ligand, gated by valence sanity.
    Returns a list of (eta_idxs, eta_n) the ring may bind through, lowest-disruption
    first and ALWAYS including the full face.  An η-n face can ring-slip to a
    CONTIGUOUS sub-arc of the ring (η6→η4→η2 for arenes/dienes, η5→η3→η1 for Cp/allyl)
    — the canonical odd/even hapticity ladder.  We only emit a slip when the sub-arc
    atoms are mutually contiguous through ring bonds (graph-checked), so no impossible
    mode is invented; substituents stay rigidly attached (whole ligand still rigid)."""
    full = list(lg["eta_local_idxs"])
    n = len(full)
    modes = [(full, n)]
    if n < 3:
        return modes
    mol = lg["mol"]
    # ring-bond adjacency among the η carbons (graph-only)
    fs = set(full)
    adj = {i: [] for i in full}
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a1 in fs and a2 in fs:
            adj[a1].append(a2)
            adj[a2].append(a1)
    # canonical ring traversal order (follow adjacency from the lowest index)
    seq = [full[0]]
    seen = {full[0]}
    while len(seq) < n:
        nxt = next((x for x in adj.get(seq[-1], []) if x not in seen), None)
        if nxt is None:
            return modes                      # not a simple ring path -> no slip
        seq.append(nxt)
        seen.add(nxt)
    # hapticity ladder: drop 2 carbons at a time (one from each end of the arc) so the
    # bound arc stays contiguous + centred; preserve odd/even parity (η6→η4→η2,
    # η5→η3→η1).  Stop at η1 (single σ-C, no longer a face — handled as edge case).
    k = n - 2
    while k >= 1:
        drop = (n - k) // 2
        arc = seq[drop: drop + k] if k > 1 else [seq[n // 2]]
        if len(arc) == k and k >= 1:
            modes.append((sorted(arc), k))
        k -= 2
    return modes


def _cp_pucker_amps(eta_n, aromatic):
    """Cremer-Pople-style out-of-plane pucker amplitudes (Å) for a NON-aromatic
    η-ring (diene / allyl / cyclohexadienyl).  Aromatic faces (Cp, arene) are planar
    -> no pucker (amp 0 only).  Returns [0.0] for planar faces (single conformer) and
    a small symmetric set for puckerable rings."""
    if aromatic or int(eta_n) < 4:
        return [0.0]
    return [0.0, 0.15, -0.15]


def _eta_is_aromatic(lg):
    """True if every η-ring atom is aromatic (planar face: Cp / arene)."""
    mol = lg["mol"]
    try:
        return all(mol.GetAtomWithIdx(int(j)).GetIsAromatic()
                   for j in lg["eta_local_idxs"])
    except Exception:
        return True


def _rmsd_aligned(A, B):
    """Heavy-rigid coordinate RMSD between two SAME-ORDERING coordinate sets after a
    proper-rotation Kabsch superposition (both already metal-centred at origin).
    Deterministic; used only to dedup ensemble members of ONE complex."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    if A.shape != B.shape or A.shape[0] < 3:
        return 1e9
    Ac = A - A.mean(0); Bc = B - B.mean(0)
    R = _kabsch_rot(Ac, Bc)
    D = Ac @ R.T - Bc
    return float(np.sqrt((D * D).sum() / len(A)))


def _dedup_builds(builds, rmsd_tol=0.25):
    """RMSD-deduplicate a list of (syms, P, donors, exempt) builds (same complex, so
    identical atom ordering).  Keeps the FIRST occurrence (emission order = canonical
    build first), dropping any later build within ``rmsd_tol`` Å of a kept one.
    Deterministic.  Distinct-by-atom-count builds (ring-slip changes nothing in the
    atom list, so counts always match) are compared directly."""
    kept = []
    for b in builds:
        P = b[1]
        dup = False
        for kb in kept:
            if kb[1].shape == P.shape and _rmsd_aligned(kb[1], P) < rmsd_tol:
                dup = True
                break
        if not dup:
            kept.append(b)
    return kept


def _hapto_base_variants(d):
    """Deterministic base variant grid for the rigid-hapto ensemble: η-ring rotamers
    (symmetry-reduced spin), valence-gated η/σ ring-slip isomers, Cremer-Pople pucker
    (non-aromatic faces).  Returns ([variant, ...], eta_lis) with the canonical
    all-base variant first, or (None, None) if the decompose dict carries no η-face.
    Graph-only, deterministic; factored out so Hebel A can reuse the same base grid."""
    ligands = d["ligands"]
    eta_lis = [i for i, lg in enumerate(ligands) if lg.get("is_eta")]
    if not eta_lis:
        return None, None
    per_eta_choices = {}
    for li in eta_lis:
        lg = ligands[li]
        aromatic = _eta_is_aromatic(lg)
        slips = _slip_modes(lg)                       # [(eta_idxs, eta_n), ...]
        choices = []
        for (eidx, en) in slips:
            spins = _ring_spins(en)
            for sp in spins:
                for pk in _cp_pucker_amps(en, aromatic):
                    choices.append({"spin": sp, "pucker": pk,
                                    "eta_idxs": eidx, "eta_n": en})
        per_eta_choices[li] = choices
    # Vary ONE η-ligand at a time off the canonical base; canonical all-base leads.
    base_choice = {li: per_eta_choices[li][0] for li in eta_lis}
    variants = [{"eta": dict(base_choice)}]           # member 0 = canonical
    for li in eta_lis:
        for ch in per_eta_choices[li][1:]:
            ev = dict(base_choice)
            ev[li] = ch
            variants.append({"eta": ev})
    return variants, eta_lis


def assemble_hapto_ensemble(metal, geometry, d, max_builds=30):
    """Deterministic RIGID-hapto ENSEMBLE: enumerate η-ring rotamers (symmetry-
    reduced spin), η/σ ring-slip isomers (valence-gated), and Cremer-Pople ring
    pucker (non-aromatic faces), assembling each as a fully rigid build (no collapse),
    then RMSD-deduplicate.  Returns a list of (syms, P, donors, exempt_pairs) builds
    (canonical build first), or None on total failure.

    Universal, graph-only, deterministic (fixed seeds; canonical ordering); the
    canonical build (member 0) is byte-identical to ``assemble_hapto(...)``."""
    variants, eta_lis = _hapto_base_variants(d)
    if variants is None:
        return None

    builds = []
    for var in variants:
        if len(builds) >= max_builds * 3:             # generous cap before dedup
            break
        try:
            b = assemble_hapto(metal, geometry, d, variant=var)
        except Exception:
            b = None
        if b is None:
            continue
        builds.append(b)
    if not builds:
        return None
    builds = _dedup_builds(builds)
    return builds[:max_builds] if builds else None


def assemble_hapto_axis_rotants(metal, geometry, d, n_axis=8, max_builds=60):
    """Hebel A (DELFIN_FFFREE_HAPTO_AXIS_ROT) — STRICTLY ADDITIVE η-axis rotamers.

    For every base variant (same grid as ``assemble_hapto_ensemble``) emit the rigid
    builds that additionally rotate the non-primary-η block (CO-tripod / co-ligands /
    a 2nd ring) about the PRIMARY η M→centroid axis by a uniform 0–2π sweep (the
    canonical 0 step is EXCLUDED — those are exactly the base ensemble builds, which
    the caller emits separately).  The rotation axis runs through the metal at the
    origin, so M-centroid + every M-D distance is INVARIANT and the η-ring internal
    geometry is untouched; only the ring-vs-rest azimuthal clock changes.

    Returns a list of (syms, P, donors, exempt_pairs) builds (NEW orientations only),
    or [] if no η-face / no rotatable block.  RMSD-deduplicated.  Deterministic."""
    variants, eta_lis = _hapto_base_variants(d)
    if variants is None:
        return []
    prim_n = int(d["ligands"][eta_lis[0]].get("eta_n", 6))
    angles = _axis_rot_angles(n_axis, prim_n)
    builds = []
    for var in variants:
        if len(builds) >= max_builds * 3:
            break
        for a in angles:
            v = {"eta": dict(var.get("eta", {})), "axis_rot": float(a)}
            try:
                b = assemble_hapto(metal, geometry, d, variant=v)
            except Exception:
                b = None
            if b is not None:
                builds.append(b)
            if len(builds) >= max_builds * 3:
                break
    if not builds:
        return []
    builds = _dedup_builds(builds)
    return builds[:max_builds]


def assemble_from_config(metal, geometry, config, ligands, refine=True,
                         n_frames=1, per_lig_confs=6, rmsd_dedup=0.5,
                         planar_bite=None, planar_coplanar=None):
    """Build a 3D complex from a chelate-isomer config (vertex -> (ligand_idx,
    arm_idx)) and the decomposed ligand list.  Chelating ligands are Kabsch-fit
    onto their two assigned vertices; monodentate ligands are oriented onto their
    vertex.  Multi-conformer selection per ligand + constrained refine.

    ``n_frames``: with the default ``1`` the single clash-minimal frame is returned
    (byte-identical to the historic behaviour) as ``(syms, P, donors)``.  With
    ``n_frames > 1`` an RMSD-deduped ENSEMBLE of full-complex frames is returned as
    a list ``[(syms, P, donors), ...]`` (canonical clash-minimal frame first) that
    varies ligand internal / chelate-ring conformation while keeping the
    constructed coordination (vertex directions + M-D distances) fixed -- the same
    best-of-ensemble lever proven for CN2 / rigid-hapto, generalised to the chelate
    σ sub-path (Task A.1).  Returns ``None`` on failure (in either mode).
    Deterministic (fixed SEED, single-thread); every emitted frame is internally
    refined; ETKDG / metallacycle conformer pool already samples chelate-ring
    pucker (Cremer-Pople) across conformers, so distinct puckers fall out naturally.
    """
    ref = MSB._ref_vectors(geometry)
    # group config by ligand instance: lig_idx -> [(vertex, arm), ...]
    by_lig = {}
    for v, (li, arm) in config.items():
        by_lig.setdefault(li, []).append((v, arm))
    out_syms = [metal]; placed = [np.zeros(3)]; placed_syms = [metal]
    fixed = {0}; pos = 1
    relax_frags = []          # (AddHs(lg.mol), ligands-only offset) for the internal relax
    # ENSEMBLE refactor (Task A.1): collect per-ligand candidate placements
    # (Q, clash_vs_metal) deduped by intra-ligand RMSD instead of only the single
    # clash-minimal pick, then enumerate full-complex combinations.  n_frames==1
    # keeps only the all-best combo -> byte-identical to the historic single path.
    ensemble = int(n_frames) > 1
    _k_confs = max(int(per_lig_confs), 1) if ensemble else 6
    per_lig_cands = []        # per ligand: [(Q, clash_vs_metal), ...] sorted best-first
    per_lig_syms = []         # per ligand: lsyms (for combo assembly)
    metal_P = np.zeros((1, 3)); metal_sym = [metal]
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
                        * MSB.md_distance(metal, _delems[i],
                                          atom=lg["mol"].GetAtomWithIdx(int(_dons_d[i])),
                                          mol=lg["mol"]) for i in range(_dent)]
            except Exception:
                _dtp = None
            # CHELATE-BACKBONE hardening is SCOPED to the newly-admitted LARGE backbones
            # (per-arm heavy > the historic cap of 8) so every existing native chelate
            # (per-arm <= 8: en, acac, bipy, terpy, ...) keeps the EXACT historic embed
            # -> byte-identical to OFF.  Only an over-cap arm (BIQCOV/ABEZAJ-class, only
            # reachable at all because the flag lifted the decompose cap) gets the soft
            # donor-donor windows + more conformers.  Flag-gated; OFF => harden=False
            # everywhere => byte-identical.
            _nheavy_lg = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
            _harden = (os.environ.get("DELFIN_FFFREE_CHELATE_BACKBONE", "0") == "1"
                       and _nheavy_lg / max(_dent, 1) > 8.0)
            # RIGID PLANAR tridentate on CN5 (TBP-5 / SPY-5): force the DG bite
            # constraint so the embed carries the meridional bite (the rigid-Kabsch
            # orient cannot open a folded ~110deg bite onto the TBP/SPY meridian).
            # Scoped to CN5 + rigid_planar + dent 3.  planar_bite overrides the env
            # gate (the never-worse caller builds BOTH ways); when None the env flag
            # DELFIN_FFFREE_PLANAR_MER_CN5 decides (else byte-identical: force_bite False).
            _eligible_rp_cn5 = bool(
                lg.get("rigid_planar") and _dent == 3
                and geometry in ("TBP-5 trigonal bipyramid", "SPY-5 square pyramid"))
            if planar_bite is None:
                _force_bite = bool(_eligible_rp_cn5 and _planar_mer_cn5_enabled())
            else:
                _force_bite = bool(_eligible_rp_cn5 and planar_bite)
            ring_confs = _embed_metallacycle(lg["mol"], _dons_d, metal,
                                             donor_target_pos=_dtp, harden=_harden,
                                             force_bite=_force_bite)
            # RIGID PLANAR polydentate (terpy / pincer, dent >= 3): the metallacycle
            # embed FOLDS the rigid backbone so the metal lifts ~1.0 A OUT of the
            # donors' own plane.  When enabled, build a metal-COPLANAR conformer
            # (the metal solved INTO the rigid donor plane) so the in-plane meridional
            # pose is realised at placement time.  planar_coplanar overrides the env
            # gate (the never-worse caller builds BOTH ways); when None the env flag
            # DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE decides (else byte-identical:
            # the coplanar conformer is never built).  Scoped to rigid_planar dent>=3;
            # falls back to the folded metallacycle embed if the coplanar solve fails.
            _eligible_rp_cop = bool(lg.get("rigid_planar") and _dent >= 3
                                    and _dtp is not None)
            _use_cop = (bool(_eligible_rp_cop and _planar_polydentate_place_enabled())
                        if planar_coplanar is None
                        else bool(_eligible_rp_cop and planar_coplanar))
            if _use_cop:
                _mds = [float(np.linalg.norm(p)) for p in _dtp]
                _cop = _coplanar_metal_centered_conformer(
                    lg["mol"], _dons_d, metal, _mds)
                if _cop is not None:
                    ring_confs = _cop
            elif (_dent >= 2 and _dtp is not None
                  and os.environ.get("DELFIN_FFFREE_PI_RIGID_PLACE", "0") == "1"):
                # PHASE 2 — FLEXIBLE multi-arm polydentate (NOT rigid_planar, e.g. a
                # tripodal with pyridyl arms on an sp3 backbone = 91 % of the metal-
                # out-of-plane defects).  Choose the backbone conformer that lets the
                # metal sit COPLANAR with EVERY aromatic arm's ring (per-ring solve;
                # _multiarm_coplanar.best_conformer).  Self-gating: returns None when
                # the donors are not in aromatic rings -> falls back to the metallacycle
                # embed.  Default-OFF / byte-identical when the flag is unset.
                try:
                    from delfin._multiarm_coplanar import embed_coplanar_multiarm as _ma_emb
                    _ma = _ma_emb(
                        lg["mol"], _dons_d, metal,
                        [float(np.linalg.norm(p)) for p in _dtp],
                        donor_target_pos=_dtp)
                    if _ma is not None:
                        ring_confs = _ma
                except Exception:
                    pass
        if ring_confs is not None:
            lsyms, coords_list = ring_confs
        else:
            confs = _ligand_confs_from_mol(lg["mol"], k=_k_confs)
            if confs is None:
                return None
            lsyms, coords_list, lmol = confs
        # Build every valid placement for this ligand across the conformer pool.
        # Single mode: pick the clash-minimal vs already-placed atoms (historic,
        # byte-identical).  Ensemble mode: keep distinct low-clash-vs-metal
        # candidates (placement-order-independent + diverse), deduped intra-ligand.
        best_Q, best_clash = None, 1e18
        cands = []                                  # (Q, clash_vs_metal) for ensemble
        seen_local = []                             # intra-ligand RMSD dedup
        for lP in coords_list:
            if lg["denticity"] == 1:
                v = va[0][0]
                Vunit = ref[v] / np.linalg.norm(ref[v])
                md = MSB.md_distance(metal, lsyms[dons[0]],
                                     atom=lg["mol"].GetAtomWithIdx(int(dons[0])),
                                     mol=lg["mol"])
                if len(lsyms) == 1:
                    Q = (Vunit * md).reshape(1, 3)
                else:
                    lPv, lp = _vsepr_reconstruct(lsyms, lP, lmol, dons[0])
                    Q = (lPv - lPv[dons[0]]) @ _rot_align(lp, -Vunit).T + Vunit * md
                    # Post-placement diatomic-orientation guard (env-gated, default OFF
                    # => byte-identical): for a linear diatomic donor (M-C#O carbonyl /
                    # M-C#N cyanide / M-N=O nitrosyl) make sure the SMILES-bonded DONOR
                    # atom — not its partner — sits at the vertex.  A shared conformer
                    # cache (keyed by canonical SMILES) can hand back a block whose atom
                    # ordering differs from this ligand's graph, flipping the unit to an
                    # M-O-C isocarbonyl; the guard reflects it donor-first.  Connectivity-
                    # derived donor element, no per-case logic.
                    if _diatomic_orient_enabled():
                        _dp = _diatomic_donor_partner(lg)
                        if _dp is not None:
                            Q = _orient_diatomic_block(
                                Q, lsyms, _dp[0], _dp[1], np.zeros(3), Vunit * md)
            else:
                # chelate (bi-/tri-dentate): the d assigned vertices host the d donors;
                # arm index i (config order) seats on verts[i] in FIXED correspondence
                # (config-faithful), using the canonical element-sorted arm order so the
                # enumeration bijects onto the distinct stereoisomers.
                dent = lg["denticity"]
                dons_d = _canonical_arm_order(lg, dent)
                verts = [v for v, arm in sorted(va, key=lambda x: x[1])]
                targets = [ref[verts[i]] / np.linalg.norm(ref[verts[i]])
                           * MSB.md_distance(metal, lsyms[dons_d[i]],
                                             atom=lg["mol"].GetAtomWithIdx(int(dons_d[i])),
                                             mol=lg["mol"]) for i in range(dent)]
                # config-faithful seating is only meaningful for ASYMMETRIC chelates
                # (distinct donor elements); for symmetric chelates every arm seating
                # is the same isomer, so keep the best-fit (lower ring strain).
                _de = [lsyms[di] for di in dons_d]
                _asym = len(set(_de)) > 1
                # RIGID PLANAR tridentate: the CENTRAL donor MUST seat on the central
                # meridian vertex (arm 1 -> targets[1]); the best-fit permutation would
                # scramble central/outer (a flat terpy is element-symmetric all-N), so
                # force config-faithful seating to realise the meridional arrangement.
                if lg.get("rigid_planar") and dent == 3:
                    _asym = True
                _rigid_planar = bool(lg.get("rigid_planar")) and dent == 3
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
                    # RIGID PLANAR tridentate (DELFIN_FFFREE_PLANAR_MER): the DG
                    # metallacycle embed yields chelating conformers (outer-outer ~150-
                    # 158deg) but a few carry a collapsed donor-backbone bond.  Reject
                    # such conformers UNCONDITIONALLY (not env-gated) so a CLEAN one is
                    # selected from the pool -> the meridional placement survives the
                    # self-gate instead of bailing the whole complex to legacy.
                    if Q is not None and _rigid_planar and _collapsed_heavy_bonds_strict(lsyms, Q):
                        Q = None
                if Q is None and dent == 2:         # embed/orient failed -> rigid fallback (bidentate only)
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1], targets[0], targets[1])
                if Q is None:                       # tridentate embed failure -> skip this config
                    continue
            if not np.all(np.isfinite(Q)):
                continue
            cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
            if cl < best_clash:
                best_clash, best_Q = cl, Q
            if ensemble:
                # dedup conformers of THIS ligand by intra-ligand RMSD (identity corr.)
                dup = False
                for Qs in seen_local:
                    if Qs.shape == Q.shape:
                        c0 = Q - Q.mean(axis=0); c1 = Qs - Qs.mean(axis=0)
                        if float(np.sqrt(((c0 - c1) ** 2).sum(axis=1).mean())) < rmsd_dedup:
                            dup = True
                            break
                if not dup:
                    seen_local.append(Q)
                    cands.append((Q, _clash_count(Q, metal_P, lsyms, metal_sym)))
            elif cl == 0:
                break
        if best_Q is None:                  # no conformer could be placed -> bail to legacy
            return None
        if ensemble:
            if not cands:
                return None
            cands.sort(key=lambda t: t[1])          # low clash-vs-metal first
            per_lig_cands.append(cands)
            per_lig_syms.append(lsyms)
        out_syms += lsyms
        placed_syms += lsyms
        relax_frags.append((Chem.AddHs(lg["mol"]), lig_offset))   # bonds for internal relax
        # in single mode the placed coords come from best_Q; the ensemble assembles
        # placed coords per combo below (placed[] is only the single-mode buffer)
        if not ensemble:
            for row in best_Q:
                placed.append(row)
        # freeze the donor atoms at their vertices
        for d in dons:
            fixed.add(pos + d)
        pos += len(lsyms)
    donors = sorted(fixed - {0})              # global indices of the constructed donor atoms
    if not ensemble:
        P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
        P = _finish_config_frame(out_syms, P, fixed, relax_frags, refine)
        return out_syms, P, donors

    # ENSEMBLE assembly (Task A.1): enumerate the Cartesian product of per-ligand
    # conformer choices, greedily ordered (sum of candidate ranks -> the best
    # combinations first = frame 0 == the single-path pick), refine each, RMSD-dedup
    # at the complex level, keep up to n_frames.  Capped product keeps it bounded.
    import itertools as _it
    rank_lists = [list(range(len(c))) for c in per_lig_cands]
    combos = list(_it.product(*rank_lists))
    combos.sort(key=lambda cb: (sum(cb), cb))      # deterministic; frame 0 = all-best
    MAX_EVAL = 64
    frames = []                                    # (syms, P) kept (deduped)
    for cb in combos[:MAX_EVAL]:
        blocks = [np.zeros((1, 3))]
        ok = True
        for vi, ci in enumerate(cb):
            Q = per_lig_cands[vi][ci][0]
            blocks.append(Q)
        Pc = np.vstack(blocks)
        if not np.all(np.isfinite(Pc)):
            continue
        Pc = _finish_config_frame(out_syms, Pc, fixed, relax_frags, refine)
        if not np.all(np.isfinite(Pc)):
            continue
        dup = False
        for _, Pk in frames:
            if Pk.shape == Pc.shape and _complex_rmsd(out_syms, Pc, Pk) < rmsd_dedup:
                dup = True
                break
        if dup:
            continue
        frames.append((list(out_syms), Pc))
        if len(frames) >= int(n_frames):
            break
    if not frames:
        return None
    return [(syms, P, donors) for syms, P in frames]


def _finish_config_frame(out_syms, P, fixed, relax_frags, refine=True):
    """Shared finishing tail for a single assembled chelate-config frame: the
    FUNDAMENTAL internal relaxation (division-of-labor doctrine) — build the
    ligands-only mol (NO metal) at the placed coords and UFF-relax it with the
    DONORS FIXED + inter-fragment vdW ON.  All organic internals (bond lengths,
    angles, funcgroup/aromatic planarity, H geometry) recover to MM-ideal and
    inter-ligand clashes resolve, while the constructed coordination is preserved
    (donors pinned -> M-D invariant; metal never enters the force field).  This
    replaces the weak defect-count refine(), which left ETKDG-rough internals.
    INTERNAL finishing.  THREE-WAY RACE design:
      Track 3 (FF-FREE, DEFAULT): pure construction -- GEOMETRIC correctors only
        (sp2-flatten projection + defect-count clash-relief).  No force field.
      Track 2 (ligand-FF, opt-in DELFIN_FFFREE_LIGANDFF=1): additionally a constrained
        UFF relax of the organic internals (metal + donors frozen) before the geometric
        correctors.  Track 1 (UFF) is the legacy path (DELFIN_FFFREE_BUILDER=0).
    The complex mol (metal + ligand bonds) is the sp2-flatten template (and the optional
    UFF mol); metal + donors stay frozen so the constructed coordination is preserved.
    Returns the (possibly relaxed) coordinate array; never raises."""
    if not refine:
        return P
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
    # #308 whole-complex torsion-space clash relax (env-gated, default-OFF byte-id):
    # joint multi-axis torsion of all rotatable single bonds, metal + donors (`fixed`)
    # frozen.  Chelate ring + M-D arms are ring bonds -> kept rigid by construction;
    # only the σ-arms/substituents rotate.  Torsion-only, never-worse.  No-op when unset.
    try:
        from delfin.fffree import torsion_relax as _TR
        bp = None
        try:
            # true connectivity: ligand-internal bonds (lig_offset shifts the AddHs
            # mol's local indices to global = lig_offset+1, atom 0 = metal) + M-D bonds
            # to every constructed donor (so the relaxer keeps chelate metallacycles
            # rigid via ring detection and never mis-perceives a crowded contact).
            pairs = []
            for frag, lig_off in relax_frags:
                base = lig_off + 1                 # global index of the frag's atom 0
                for b in frag.GetBonds():
                    pairs.append((base + b.GetBeginAtomIdx(), base + b.GetEndAtomIdx()))
            for dg in sorted(fixed - {0}):
                pairs.append((0, dg))
            bp = pairs or None
        except Exception:
            bp = None
        P = np.asarray(_TR.relax_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                       dtype=float)
        # JOINT inter-ligand declash (env-gated, default-OFF byte-id): global
        # inter-ligand heavy-heavy minimisation, core frozen.  After #308.  Reuses
        # the same true-connectivity bond list.  No-op when the flag is unset.
        from delfin.fffree import joint_declash as _JD
        P = np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                       dtype=float)
    except Exception:
        pass
    return P


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
        md = MSB.md_distance(metal, lsyms[di],
                             atom=lmol.GetAtomWithIdx(di), mol=lmol)
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
        md = MSB.md_distance(metal, lsyms[donor_idx],
                             atom=m.GetAtomWithIdx(donor_idx), mol=m)
        lp = _donor_and_lp(lsyms, lP, m, donor_idx)
        syms = [metal]; blocks = [np.zeros((1, 3))]
        for i in range(n):
            Vunit = ref[i] / np.linalg.norm(ref[i])
            R = _rot_align(lp, -Vunit)
            Q = (lP - lP[donor_idx]) @ R.T + Vunit * md
            syms += lsyms; blocks.append(Q)
        out.append((round(e, 2), syms, np.vstack(blocks)))
    return out
