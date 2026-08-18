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
import math
import os
import itertools
from typing import List, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from delfin.manta import metal_sphere_builder as MSB
import delfin.manta._bond_decollapse as _bd

SEED = 42


def _rot_align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix mapping unit vector a -> unit vector b (Rodrigues)."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(np.dot(a, b))
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        # ANTIPARALLEL (a -> -a).  `-np.eye(3)` maps a to b correctly but has det = -1: it is an
        # INVERSION, not a rotation.  Applied to a CHIRAL ligand it silently emits the ENANTIOMER --
        # a different molecule that no dedup or mirror filter downstream catches, and that shows up
        # in the eye as a lost CCDC isomer rather than as the placement bug it is.  The proper
        # replacement is a 180 deg rotation about ANY axis perpendicular to a: R = 2*n n^T - I with
        # n _|_ a.  It maps a -> -a exactly as well, but det = +1, so handedness is preserved.
        # Env-gated (DELFIN_FFFREE_PROPER_ANTIPARALLEL, default OFF -> byte-identical) so the
        # correction gets its own A/B; watch ccdc_isomer_lost, where the mirrored ligand surfaces.
        if os.environ.get("DELFIN_FFFREE_PROPER_ANTIPARALLEL", "0") != "1":
            return -np.eye(3)
        _perp = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(_perp, a))) > 0.9:
            _perp = np.array([0.0, 1.0, 0.0])
        _n = np.cross(a, _perp)
        _n = _n / np.linalg.norm(_n)
        return 2.0 * np.outer(_n, _n) - np.eye(3)
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


def _ring_bounds_enabled() -> bool:
    """Tighten the aromatic 5-/6-ring INTERIOR ANGLE in the metallacycle DG embed.

    Eye-find QEBLOC (2026-06-24): the κ4 cage cap's rigid 5-ring embeds with its
    interior angle OPEN (~115deg vs the ideal ~108deg) — the unconstrained ETKDG
    embed + kekulization-skip (KEKULIZE_SPLIT) leaves the ring-atom 1-3 bounds loose,
    so even rigid-body seating inherits the ~6deg ring splay.  This flag adds an
    interior-angle 1-3 bound (regular-polygon ideal: 108deg for 5-rings, 120deg for
    6-rings) so the ring embeds flat-and-correct.  Geometry-only, feasibility-safe
    (falls back to the unconstrained embed).  Default OFF -> byte-identical."""
    return os.environ.get("DELFIN_FFFREE_RING_BOUNDS", "0") == "1"


def _bite_free_bounds() -> bool:
    """THE one place DELFIN_FFFREE_BITE_FREE is read (default OFF -> byte-identical).

    THE BITE IS A RING CLOSURE, NOT AN INPUT.  Measured 2026-07-31 over 18 metals:
        cos(bite) = (r1^2 + r2^2 - d_DD^2) / (2 r1 r2)      RMS 0.24 deg, no fitted parameter.
    Two M-D radii and ONE donor-donor distance already determine the bite exactly.  So a
    bounds matrix that pins BOTH the radii AND the donor-donor separation is not stating one
    law twice -- it is OVER-DETERMINING a triangle, and when the second statement comes from
    the ideal polyhedron while the ligand's backbone says otherwise, the two contradict.
    Triangle smoothing spreads that contradiction over the whole matrix and the embedder
    returns the least-bad compromise: that is the splay described at the TRILATERATED BOUNDS
    comment below, and the source of chelate bites beyond 92 deg (real chelates: 98 % below
    90, none above 95).

    This flag keeps the half that is measured and real -- the M-D radii -- and DROPS the
    donor-donor entries, leaving them to RDKit's own covalent bounds, i.e. to the separations
    the ligand actually has.  The bite then comes OUT of the embed as the consequence of the
    ring closure instead of being asked for, and there is nothing left to correct afterwards.

    Distinct from the trilateration at line 273: that one KEEPS the donor-donor entries and
    moves the targets instead, which needs a pre-existing conformer to trilaterate from (see
    the GetConformer guard there).  This needs nothing.  It wins over the bite pin when both
    are on, because it is the more specific statement about the very same entries."""
    return os.environ.get("DELFIN_FFFREE_BITE_FREE", "0") == "1"


def _tighten_ring_bounds(bm, mh, tol=0.06):
    """Set the 1-3 distance bound of every aromatic 5-/6-ring to the regular-polygon
    interior angle, using the bounds matrix's OWN 1-2 midpoints (so RDKit's bond
    lengths are respected and only the ANGLE is enforced).  In place; never raises."""
    try:
        ri = mh.GetRingInfo()
    except Exception:
        return

    def _mid(i, j):
        lo, hi = (i, j) if i < j else (j, i)
        return 0.5 * (float(bm[lo][hi]) + float(bm[hi][lo]))
    for ring in ri.AtomRings():
        n = len(ring)
        if n not in (5, 6):
            continue
        ct = float(np.cos(np.radians(180.0 * (n - 2) / n)))   # 108deg (5) / 120deg (6)
        for a in range(n):
            i, j, kk = ring[a], ring[(a + 1) % n], ring[(a + 2) % n]
            d12 = _mid(i, j); d23 = _mid(j, kk)
            d13 = float(np.sqrt(max(d12 * d12 + d23 * d23 - 2 * d12 * d23 * ct, 0.0)))
            lo, hi = (i, kk) if i < kk else (kk, i)
            bm[lo][hi] = float(d13 + tol)
            bm[hi][lo] = float(max(d13 - tol, 0.0))


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
        _sops = (Chem.SanitizeFlags.SANITIZE_ALL
                 ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        # Kekulize-robust (DELFIN_FFFREE_KEKULIZE_SPLIT, default OFF -> byte-id): the
        # cleaved aromatic-N⁺ cap (triazolide / pyridinium scorpionate) is unkekulizable,
        # so the strict sanitize fails and the metallacycle embed returns None -> legacy.
        # Skip kekulization (aromaticity retained) so the DG embed can proceed.
        if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") == "1":
            _sops ^= Chem.SanitizeFlags.SANITIZE_KEKULIZE
        try:
            Chem.SanitizeMol(m, sanitizeOps=_sops)
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
        _ring_b = _ring_bounds_enabled()
        # BITE_FREE opens the bounds path on its own: it needs the M-D radii out of the
        # targets, so _have_targets, but NOT the bite pin -- dropping the donor-donor entries
        # is the whole point of it (see _bite_free_bounds).
        _bite_free = _bite_free_bounds()
        _use_bm = (_have_targets and (_chel_bite or harden or force_bite or _bite_free)) \
            or _oc6_vertex or _ring_b
        if _use_bm:
            try:
                from rdkit.Chem import rdDistGeom as _DG
                from rdkit import DistanceGeometry as _DGs
                nA = mh.GetNumAtoms()
                bm = _DG.GetMoleculeBoundsMatrix(mh)
                dset = {int(d) for d in donor_idxs}
                tp = ([np.asarray(p, float) for p in donor_target_pos]
                      if _have_targets else None)
                # RING_BOUNDS-only (rigid cage): do NOT pin donor-donor to the ideal
                # polyhedron (that re-introduces the cap splay the rigid seating fixes);
                # only enforce the ring interior angle below.  Keep the cap's natural
                # donor geometry -> coordination stays emergent.
                # BITE_FREE belongs in this list too, and for the opposite reason: it wants
                # the M-D pins that only tp carries, and it drops the donor-donor entries by
                # itself further down.  Dropping tp here would silently take its radii away
                # as well, so RING_BOUNDS+BITE_FREE would quietly become RING_BOUNDS alone.
                if _ring_b and not (_chel_bite or harden or force_bite or _bite_free):
                    tp = None

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
                    # TRILATERATED BOUNDS (DELFIN_FFREE_TRILAT_BOUNDS=1, default OFF).
                    #
                    # This block already IS "specify and solve": bounds matrix -> triangle
                    # smoothing -> embed.  The solver is not the problem.  Its INPUT is:
                    # tp are the IDEAL POLYHEDRON vertex positions, and the donor-donor
                    # entries derived from them ask for a separation the ligand does not
                    # have.  Triangle smoothing then propagates that impossible request
                    # through the whole matrix, and the embedder returns the least-bad
                    # compromise -- which is where the splayed rings and the 29 % of
                    # chelate bites beyond 92 deg come from (real chelates: 98 % below 90,
                    # none above 95).
                    #
                    # So keep the M-D radii, which are measured and real, and replace the
                    # donor-donor entries with the separations the LIGAND ACTUALLY HAS.
                    # Same statement as the cosine law for two donors, generalised: the
                    # coordination angles become a consequence of the bounds instead of an
                    # input to them, and the enumerated vertex assignment is untouched
                    # because trilateration starts from those very targets.
                    _tp_use = tp
                    if os.environ.get("DELFIN_FFREE_TRILAT_BOUNDS", "0") == "1":
                        try:
                            _lp_src = mh.GetConformer().GetPositions()
                        except Exception:
                            _lp_src = None
                        if _lp_src is not None:
                            _t2 = _trilaterate_donor_targets(
                                _lp_src, [int(d) for d in donor_idxs], list(tp))
                            if _t2 is not None:
                                _tp_use = _t2
                    for a, da in enumerate(donor_idxs):   # M-D HARD; donor-donor (soft for backbone)
                        _setb(mi, int(da), float(np.linalg.norm(_tp_use[a])))
                        if _bite_free:
                            # RING CLOSURE: the two radii above already fix the bite once the
                            # ligand's own donor-donor distance is known, so writing that
                            # distance from the ideal polyhedron as well over-determines the
                            # triangle.  Leave it to RDKit's covalent bounds and let the bite
                            # come out.  Deliberately AFTER the M-D line, not instead of it.
                            continue
                        for b in range(a + 1, len(donor_idxs)):
                            _setb(int(da), int(donor_idxs[b]),
                                  float(np.linalg.norm(_tp_use[a] - _tp_use[b])),
                                  tol=_dd_tol)
                if _DGs.DoTriangleSmoothing(bm):
                    ep = _DG.EmbedParameters()
                    ep.randomSeed = SEED
                    ep.useRandomCoords = True
                    ep.SetBoundsMat(bm)
                    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, params=ep))
            except Exception:
                cids = []
        if not cids:                                      # default / fallback: unconstrained embed
            try:
                cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, randomSeed=SEED,
                                                       numThreads=1, useRandomCoords=False))
            except Exception:
                if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") != "1":
                    raise                                 # byte-identical default
                # κ4 chain: an unkekulizable aromatic-N⁺ scorpionate cap makes the
                # embedder throw -> metallacycle returns None -> cage to legacy.
                # Neutralise the artefact ring-N⁺ charges (geometry-only) and retry so
                # the chelate-aware metallacycle embed succeeds (clean OC-6 seating).
                cids = []
                try:
                    rw2 = Chem.RWMol(m)
                    for a in rw2.GetAtoms():
                        if (a.GetIsAromatic() and a.GetSymbol() == "N"
                                and a.GetFormalCharge() > 0):
                            a.SetFormalCharge(0)
                    m2 = rw2.GetMol()
                    Chem.SanitizeMol(m2)
                    mh = Chem.AddHs(m2)
                    cids = list(AllChem.EmbedMultipleConfs(
                        mh, numConfs=k, randomSeed=SEED, numThreads=1,
                        useRandomCoords=False))
                except Exception:
                    cids = []
            # EMBED-FLOOR (DELFIN_FFFREE_EMBED_FLOOR, default-OFF -> byte-id): the
            # unconstrained fallback embed has NO metal->non-donor lower bound, so a
            # NON-DONOR heavy atom (a carboxylate carbonyl-O/C, a pendant heteroatom on
            # a donor's first shell) collapses ONTO the metal (~0.7-1.1 A) -> that
            # in-shell intruder makes _build_is_clean reject the whole complex -> it
            # falls to legacy UFF whose buckled donors are the generic-sigma set-gap
            # (#1 RMSD lever, ~2/3 of it).  Generate ADDITIONAL conformers with a
            # metal->non-donor EXCLUSION bounds matrix and POOL them with the
            # unconstrained ones (never replace).  ADDITIVE = never-worse by
            # construction: the caller's clash-selection keeps the best, and the
            # original unconstrained conformers stay in the pool (so a system that was
            # already fine, e.g. MEZFOM, keeps its good frame — this is why we pool,
            # not swap: raw OC6_VERTEX swapped and regressed MEZFOM 2.05->2.54).
            # License-clean (covalent bounds + donor set from the cleave), deterministic.
            if os.environ.get("DELFIN_FFFREE_EMBED_FLOOR", "0") == "1":
                try:
                    from rdkit.Chem import rdDistGeom as _DGf
                    from rdkit import DistanceGeometry as _DGfs
                    _bmf = _DGf.GetMoleculeBoundsMatrix(mh)
                    _dsetf = {int(d) for d in donor_idxs}
                    _exclf = float(os.environ.get("DELFIN_FFFREE_OC6_NONDONOR_EXCL", "2.4"))
                    for _ai in range(mh.GetNumAtoms()):
                        if _ai != mi and _ai not in _dsetf:
                            _ish = mh.GetAtomWithIdx(_ai).GetAtomicNum() == 1
                            _lof = 1.2 if _ish else _exclf
                            _a2, _b2 = min(_ai, mi), max(_ai, mi)
                            _bmf[_a2][_b2] = max(float(_bmf[_a2][_b2]), _lof + 2.0)
                            _bmf[_b2][_a2] = float(_lof)
                    if _DGfs.DoTriangleSmoothing(_bmf):
                        _epf = _DGf.EmbedParameters()
                        _epf.randomSeed = SEED
                        _epf.useRandomCoords = True
                        _epf.SetBoundsMat(_bmf)
                        _fcids = list(AllChem.EmbedMultipleConfs(mh, numConfs=k, params=_epf))
                        cids = list(cids) + [c for c in _fcids if c not in cids]
                except Exception:
                    pass
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


def _rigid_ligand_cavity_conformer(mol, donor_idxs, metal_sym, mds, k=12):
    """LIGAND-GEOMETRY-FIRST universal seating for a RIGID POLYDENTATE (macrocycle /
    cage / conjugated pincer): embed the FREE ligand (its own internal geometry is the
    hard INVARIANT — the conjugated/macrocyclic backbone holds the donors in their
    natural arrangement), then solve, for each conformer, the 3-D metal position that
    best matches EVERY ideal M-donor distance ``mds[i]`` SIMULTANEOUSLY (Gauss-Newton
    on the over-determined distance system).  Pick the conformer whose donor cavity
    most consistently admits the metal at the ideal radii (smallest distance residual
    = least-strained seating that respects the ligand), recenter so the metal sits at
    the ORIGIN, and emit it.

    Same return contract as ``_embed_metallacycle`` / ``_coplanar_metal_centered_
    conformer``: ``(lsyms, [coords, ...])`` excluding the placeholder metal, matching
    ``AddHs(mol)`` atom order (donor indices preserved), recentered on the metal.
    Returns ``None`` on failure.

    Unlike ``_coplanar_metal_centered_conformer`` this does NOT constrain the metal to
    a backbone PLANE and does NOT require the backbone to be flat, so it is UNIVERSAL:
    a 3-D cage solves to the cavity centre; a planar κ4 macrocycle solves to the
    in-plane donor centre AUTOMATICALLY (the equidistant point of coplanar donors lies
    in their plane).  Paired with the rigid-body orient (``_orient_chelate_to_vertices
    (rigid=True)``: Kabsch fit, NO per-donor radial rescale) the ligand's internal
    bonds AND angles are preserved EXACTLY and the coordination polyhedron comes out
    EMERGENT — the "ligand-first → polyhedron emergent" construction order.  Universal,
    geometry/graph-only, deterministic (fixed seed, single thread); never raises."""
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
        if len(dons) < 2 or len(mds) != len(dons):
            return None
        keep = list(range(mh.GetNumAtoms()))
        best = None                       # (residual, recentered coords)
        for cid in cids:
            P = np.array(mh.GetConformer(cid).GetPositions(), float)
            D = P[dons]
            # 3-D Gauss-Newton: minimise sum_i (|u - D_i| - md_i)^2, seed = donor centroid
            u = D.mean(axis=0)
            r = np.zeros(len(dons))
            for _ in range(200):
                J = []; r = []
                for i in range(len(dons)):
                    diff = u - D[i]; nn = max(float(np.linalg.norm(diff)), 1e-9)
                    r.append(nn - mds[i]); J.append(diff / nn)
                J = np.array(J); r = np.array(r)
                try:
                    step = np.linalg.lstsq(J, -r, rcond=None)[0]
                except Exception:
                    break
                u = u + step
                if float(np.linalg.norm(step)) < 1e-9:
                    break
            resid = float(np.sqrt(np.mean(np.asarray(r) ** 2)))
            coords = P[keep] - u              # recenter: metal -> origin
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


def _donor_follow_weights(syms, P, donor_idxs, span):
    """Per-atom (owner donor, weight) for letting a donor's NEIGHBOURHOOD follow its move.

    THE DEFECT THIS ADDRESSES.  The per-donor radial placement further down sets each donor
    to its exact ideal M-D radius and moves NOTHING else -- the donor slides by delta while
    its own bonded neighbour stays put, so the donor-backbone bond is stretched or compressed
    by the full delta (measured median ~0.18 A).  The code says "the constrained relax then
    pulls the backbone into consistency", and that relax sits behind LIGANDFF, which is not a
    champion flag: in the shipped build nothing ever pulls it back.  That is the largest
    single source of GATE_COLLAPSED_BOND (171 of 215 CHELATE_EMPTY systems die there).

    THE FIX IS A DECAY, NOT A REPAIR.  Instead of ending the displacement abruptly at the
    donor, let it fall off linearly over the BOND GRAPH: an atom k bonds away moves by
    (1 - k/span) of its donor's delta.  The donor-neighbour bond then distorts by delta/span
    instead of delta -- with the default span of 3, a third -- and the distortion keeps
    shrinking outward instead of piling onto one bond.  No force field, no iteration, no
    reference: one BFS over a covalent graph read off the geometry.

    EVERY DONOR KEEPS ITS EXACT RADIUS.  Atoms are assigned to their NEAREST donor and a tie
    moves nothing, so no field ever reaches another donor -- the M-D lengths the placement
    just fixed, and the bite between them, are untouched by construction.
    """
    n = len(syms)
    dons = [int(d) for d in donor_idxs]
    dset = set(dons)
    adj = [[] for _ in range(n)]
    for i in range(n):
        if _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _bd._is_metal(syms[j]):
                continue
            if float(np.linalg.norm(P[i] - P[j])) <= 1.30 * _bd._ideal_bond(syms[i], syms[j]):
                adj[i].append(j)
                adj[j].append(i)
    INF = 10 ** 6
    dist = [INF] * n
    owner = [-1] * n
    frontier = []
    for d in dons:
        dist[d] = 0
        owner[d] = d
        frontier.append(d)
    k = 0
    while frontier and k < span:
        k += 1
        nxt = []
        claim = {}
        for a in frontier:
            for b in adj[a]:
                if dist[b] <= k or b in dset:
                    continue                      # already closer, or a donor: never moved
                claim.setdefault(b, set()).add(owner[a])
        for b, owners in claim.items():
            dist[b] = k
            owner[b] = next(iter(owners)) if len(owners) == 1 else -1   # tie -> no move
            nxt.append(b)
        frontier = nxt
    out = {}
    for a in range(n):
        if a in dset or owner[a] < 0 or dist[a] >= span:
            continue
        out[a] = (owner[a], 1.0 - float(dist[a]) / float(span))
    return out


def _collapsed_heavy_bonds_strict(syms, P, factor=None):   # None -> _bd.COLLAPSE_FLOOR (war 0.82)
    """True if any BONDED heavy-heavy non-metal pair sits below ``factor`` × the
    covalent-sum ideal — same logic as ``_has_collapsed_heavy_bonds`` but NOT env-
    gated (always active).  Used to reject the few DG-metallacycle conformers of a
    RIGID PLANAR tridentate that carry a collapsed donor-backbone bond, so a clean
    conformer is selected from the pool.  Universal, geometry-only, deterministic."""
    if factor is None:
        factor = _bd.COLLAPSE_FLOOR   # EINE Quelle fuer den Kollaps-Boden, siehe _bond_decollapse
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


def _hydrogens_riding_on(syms, P, idx, cut=1.35):
    """Indices of the hydrogens covalently attached to atom ``idx``.

    Geometric, so it needs no molecule object and no vocabulary: an H whose distance to
    ``idx`` is inside a covalent X-H range belongs to it.  1.35 A clears every real X-H
    (C-H 1.09, N-H 1.03, O-H 0.99, B-H 1.19, Si-H 1.48 is the only common one above it and
    a silane donor is not a case this path seats) while staying far below any H...X contact.
    """
    out = []
    if syms is None:
        return out
    p = np.asarray(P[idx], float)
    for j in range(len(syms)):
        if syms[j] != "H" or j == idx:
            continue
        if float(np.linalg.norm(np.asarray(P[j], float) - p)) <= cut:
            out.append(j)
    return out


_H_CONTACT_FLOOR = 1.6      # inter-ligand H...H stays above this in crystals (eye: hhclash)


def _riders_that_may_move(syms, P, riders, delta, parent):
    """Which riding hydrogens may follow their parent WITHOUT tightening a contact.

    MEASURED CAUSE, KIQNUT 2026-08-01.  H_FOLLOW moved every rider along the M-D radial
    direction with no clash check at all.  On a crowded CN6 hexaamine that drove
    inter-ligand hydrogens into each other -- H13..H44 went 1.34 -> 1.21 A where crystals
    stay above 1.6 -- and the tightened contact cost the topology match: topo_correct
    true -> false, broken_frac 0.0 -> 1.0.  That ONE system was the cap_LOST that blocked
    the entire A/B, which was otherwise a win (12 affected, valid 5->5, mean_delta -0.836:
    the eye read the rest as BETTER).

    The heavy-atom graph was never the problem -- BOTH arms carry the identical N-C
    compression (1.32 / 1.34 / 1.39 A vs crystal 1.47), so the rescale is not what broke
    it.  The hydrogens were.

    RULE: a rider moves only if the move does not leave its closest contact both TIGHTER
    than before and below the crystal floor.  Never-worse by construction: the outcome is
    either today's champion behaviour (the rider simply stays) or a move that does not
    tighten anything past what crystals show.

    ⛔ MEASURED NOT TO FIX KIQNUT, and the reason is structural.  isoH_FOLLOW2 came back
    identical to isoH_FOLLOW down to the MD5 of the built file: this guard moved nothing.
    It runs inside _orient_chelate_to_vertices, which works on ONE ligand's coordinates in
    the metal frame, while the clash is INTER-ligand (H13 on C2 against H44 on N31).  At
    rescale time the neighbouring ligand is not in the array, so the guard cannot see the
    contact it was written to catch -- checkable before building, and I did not check it.

    The guard is kept: it is correct for intra-ligand contacts and costs nothing where it
    cannot fire.  But H_FOLLOW's real defect is one step earlier -- the DONOR is pushed onto
    its ideal radius with no knowledge of where its hydrogen will then point, and without
    H_FOLLOW that same move merely stretches the X-H bond instead.  Both are wrong in
    different ways and neither repairs the other.  A next attempt must run on the ASSEMBLED
    complex, or better, the seating must account for the neighbour before it places the donor.
    """
    if not len(riders):
        return riders
    A = np.asarray(P, float)
    keep = []
    for h in riders:
        others = [j for j in range(len(A)) if j != h and j != parent and j not in riders]
        if not others:
            keep.append(h); continue
        O = A[others]
        d_before = float(np.min(np.linalg.norm(O - A[h], axis=1)))
        d_after = float(np.min(np.linalg.norm(O - (A[h] + delta), axis=1)))
        if d_after >= d_before or d_after >= _H_CONTACT_FLOOR:
            keep.append(h)
    return keep


# ===== TRILATERATION AS A RESCUE RUNG, NOT AS THE PRIMARY PATH ===============================
#
# Measured 2026-08-02 (trilatAB2, 995 systems).  DELFIN_FFREE_TRILATERATE as a PRIMARY path:
#     12 losses -- 12 of 12 were topo_correct BEFORE
#     11 gains  --  0 of 11 had a valid frame BEFORE
# Not one borderline case in either direction: it repairs what was broken and damages what was
# whole.  PLANAR_MER measured identically.  A switch with that signature is not a better way to
# place ligands, it is a SECOND way -- and a second way belongs where the first one already
# failed.  converter_backend already HAS that ladder (_maybe_decollapse -> _seat_via_conformers
# -> legacy); this adds a rung to it rather than a parallel mechanism.
#
# Reached only after the self-gate has rejected the rigid build, the rescue is additive BY
# CONSTRUCTION: a clean frame can never be replaced by it, so never-worse holds structurally
# instead of having to be re-measured.
#
# The env read lives HERE and only here (converter_backend asks through trilat_rescue_enabled /
# trilaterate_rescue), so the flag has one home and cannot drift out of sync with its callers.
_DP_RESID_SLACK = 0.15          # A, how much donor-to-vertex residual beta may buy back
_DP_STEPS = 72                  # 5 deg scan; the objective is smooth, no optimiser needed


def _donor_plane_beta(P, syms, d):
    """beta at donor d: angle between M->D and the plane of d's own substituents.

    The metal is NOT part of the plane fit -- that is the whole point (pyramid_root.py):
    a three-coordinate donor has one Walsh angle and it does not say WHICH of the three
    partners is displaced.  Fitting the plane WITHOUT the metal makes the question
    answerable: if the substituents stay flat and the metal is off, the seating is at fault.
    Returns None when d has fewer than two substituents (no plane exists to be out of).
    """
    _d = np.asarray(P[d], float)
    nb = [i for i in range(len(syms)) if i != d and syms[i] != "H"
          and float(np.linalg.norm(np.asarray(P[i], float) - _d)) < 1.95]
    if len(nb) < 2:
        return None
    A = np.array([np.asarray(P[i], float) - _d for i in nb], float)
    try:                                        # plane normal = smallest singular vector
        nrm = np.linalg.svd(A)[2][-1]
    except Exception:
        return None
    v = np.asarray(P[0], float) * 0.0 - _d      # M sits at the origin in this frame
    n_ = float(np.linalg.norm(v))
    if n_ < 1e-6:
        return None
    return abs(math.degrees(math.asin(max(-1.0, min(1.0, float(np.dot(nrm, v / n_)))))))


_BETA_BAND = 4.8        # deg -- the crystals' own upper beta, measured over clean CCDC
                        # structures: 3.4 monodentate, 3.7 bidentate, 4.8 tetradentate.
                        # NOT a tuning knob: below it a donor is as flat as real chemistry
                        # gets, so there is nothing to win by preferring a flatter conformer.


_COMBO_MATERIALISE_MAX = 200_000    # above this the full product is never built at all


def _ranked_combos(rank_lists, k):
    """The first ``k`` index-combinations in (sum, lexicographic) order -- WITHOUT ever
    materialising the Cartesian product.

    WHY THIS EXISTS.  Both ensemble paths did

        combos = list(itertools.product(*rank_lists))
        combos.sort(key=lambda cb: (sum(cb), cb))
        ... combos[:MAX_EVAL]

    i.e. they built and sorted the WHOLE product to use 64 of it.  Eight ligands with ten
    conformers each is 10^8 tuples -- per system, times every parallel worker.  That is the
    measured cause of four consecutive OOM kills of the sigma-ensemble path (journal:
    "Failed with result 'oom-kill'"), which is why the biggest conformer lever in the tree
    has never once produced a verdict.  Note the shape of the mistake: the memory blows up
    in the SELECTION, not in the chemistry -- a single build peaks near 0.25 G.

    EXACTLY ORDER-EQUIVALENT to the sort it replaces.  Best-first over the index lattice:
    pop the smallest (sum, tuple), push its one-step increments.  Every combination is
    reachable by incrementing coordinates from all-zeros, and the heap key IS the sort key,
    so the k-th element out is the k-th element of the sorted product.  Memory O(k * n).

    Below _COMBO_MATERIALISE_MAX the caller keeps the historic path verbatim, so nothing
    changes for the small cases that always worked."""
    import heapq as _hq
    n = len(rank_lists)
    if n == 0:
        return []
    start = tuple(0 for _ in range(n))
    heap = [(0, start)]
    seen = {start}
    out = []
    while heap and len(out) < k:
        s, cb = _hq.heappop(heap)
        out.append(cb)
        for i in range(n):
            if cb[i] + 1 < len(rank_lists[i]):
                nxt = cb[:i] + (cb[i] + 1,) + cb[i + 1:]
                if nxt not in seen:
                    seen.add(nxt)
                    _hq.heappush(heap, (s + 1, nxt))
    return out


def _beta_score(syms, Q, donor_idxs):
    """Sum of SQUARED out-of-plane angles over the donors that HAVE a plane (degrees^2).

    The metal sits at the origin in Q, which is exactly what _donor_plane_beta assumes.
    Donors with fewer than two heavy substituents have no plane to be out of and simply do
    not contribute -- a carboxylate O is not scored here, and must not be: its in-plane
    statement is a TORSION, a different quantity.

    ONLY THE EXCESS OVER THE CRYSTAL BAND IS SCORED.  Real complexes are not flat either:
    measured against clean crystals, beta sits at 3.4 deg for monodentate donors, 3.7 for
    bidentate and 4.8 for tetradentate.  A donor already inside that band is RIGHT, and
    preferring an even flatter conformer over it buys nothing while disturbing a pick that
    the historic clash order had made for a reason.

    Measured 2026-08-02, and this is why the band is here rather than a raw sum: scoring the
    raw beta (betasel) moved pyramidal_sp2 from 19.05 % to 12.70 % and the hard-finding rate
    from 50.8 % to 42.5 %, but cost 4 capabilities against 3 gained.  Pairing it with the
    collapse criterion (betacsel) made it WORSE, not better -- 6 lost against 2 -- so the four
    losses are not collapse-related; the lever is simply too eager.  _BETA_BAND is not a fitted
    knob: it is the crystals' own upper figure."""
    s = 0.0
    for d in donor_idxs:
        try:
            b = _donor_plane_beta(Q, syms, int(d))
        except Exception:
            b = None
        if b is not None:
            e = float(b) - _BETA_BAND
            if e > 0.0:
                s += e * e
    return s


def _donor_plane_relax(Q, syms, donor_idxs, tgt, tmu):
    """Rotate the whole ligand about the donor-centroid axis to put the metal into the
    donor planes.  Rigid -> bonds, internal angles and the BITE are untouched by
    construction.  Returns the improved coordinates or None if nothing beat the input."""
    def _beta_sum(X):
        s = 0.0
        for d in donor_idxs:
            b = _donor_plane_beta(X, syms, d)
            if b is not None:
                s += b * b
        return s

    def _resid(X):
        return float(np.sqrt(np.mean(np.sum(
            (np.array([np.asarray(X[d], float) for d in donor_idxs]) - tgt) ** 2, axis=1))))

    # THE FREE AXIS IS THE ONE THROUGH THE DONORS THEMSELVES, not metal -> centroid.
    #
    # Measured 2026-08-02 (donorplane, rc=3): with the metal->centroid axis the lever had
    # ZERO reach on 187 systems -- the loop's own probe refused the A/B before it could
    # report "no effect" and let me mistake a wiring fault for a verdict.  The reason is
    # geometry, not code: rotating about metal->centroid swings every donor on a CONE, away
    # from the vertex it was just fitted to, so the residual guard rightly killed every
    # candidate.  That rotation is free only for a MONODENTATE -- and the monodentate path
    # already sits at beta 0.9 deg, better than the crystals.
    #
    # Rotating about the line THROUGH the donors leaves the donors themselves on that line
    # and therefore on their targets: for a bidentate the donor-donor axis IS the bite, so
    # the bite is untouched by construction and both donors stay put to first order.  What
    # swings is the backbone -- and with it the donor planes, which is exactly the quantity
    # we want to move.  Nothing that is already placed pays for it.
    # ONLY WHERE THE ROTATION IS PROVABLY FREE: one or two donors.
    #
    # With two donors the axis is the line THROUGH both, so both stay exactly where the
    # radial reset just put them and the donor-donor distance -- the BITE -- is invariant.
    # With one donor it is the M-D bond itself: same argument, trivially.
    #
    # From THREE donors on the claim fails: the principal direction is a best-fit line that
    # passes through none of them, so every donor swings off the position it was just given.
    # Measured (donorplane3, 187 systems): reach 56, valid 32->35, 4 systems gained, but
    # OVAVEO lost and the TOPOLOGY floor broke -- exactly the polydentate case where "free"
    # was never true.  The lever keeps the part it can prove and drops the part it cannot;
    # kappa3+ needs the donors placed GLOBALLY, not one ligand rotated after the fact.
    _D = np.array([np.asarray(Q[d], float) for d in donor_idxs])
    if len(_D) > 2:
        return None
    if len(_D) == 2:
        ax = _D[1] - _D[0]                                # the line through both donors
    else:
        ax = _D[0] - np.zeros(3)                          # monodentate: the M-D bond itself
    n_ = float(np.linalg.norm(ax))
    if n_ < 1e-6:
        return None
    ax = ax / n_
    tmu = _D.mean(0)                            # rotate about the donor centroid, on-axis
    base_b, base_r = _beta_sum(Q), _resid(Q)
    best = None
    K = np.array([[0.0, -ax[2], ax[1]], [ax[2], 0.0, -ax[0]], [-ax[1], ax[0], 0.0]])
    for i in range(1, _DP_STEPS):
        th = 2.0 * math.pi * i / _DP_STEPS
        R = np.eye(3) + math.sin(th) * K + (1.0 - math.cos(th)) * (K @ K)
        X = (np.asarray(Q, float) - tmu) @ R.T + tmu
        b = _beta_sum(X)
        if b >= base_b or _resid(X) > base_r + _DP_RESID_SLACK:
            continue
        if best is None or b < best[0]:
            best = (b, X)
    return None if best is None else best[1]


_TRILAT_RESCUE = False


def trilat_rescue_enabled() -> bool:
    """THE one place DELFIN_FFREE_TRILAT_RESCUE is read (default OFF -> byte-identical)."""
    return os.environ.get("DELFIN_FFREE_TRILAT_RESCUE", "0") == "1"


class trilaterate_rescue:
    """Re-run one assembly with trilaterated donor targets instead of ideal vertices."""

    def __enter__(self):
        global _TRILAT_RESCUE
        self._prev = _TRILAT_RESCUE
        _TRILAT_RESCUE = True
        return self

    def __exit__(self, *_exc):
        global _TRILAT_RESCUE
        _TRILAT_RESCUE = self._prev
        return False


def _trilat_targets_on() -> bool:
    """Primary-path flag (legacy A/B, default OFF) OR an active rescue re-build."""
    return _TRILAT_RESCUE or os.environ.get("DELFIN_FFREE_TRILATERATE", "0") == "1"


def _orient_chelate_to_vertices(lP, donor_idxs, targets, asym=True, rigid=False, lsyms=None):
    """Rotate a metal-centered chelate conformer (from _embed_metallacycle) so its
    donors seat onto the target vertex directions, then per-donor rescale to the
    ideal M-donor distance.  The ring geometry (backbone clears the metal) is
    preserved as a rigid body.  Works for any denticity.

    LIGAND-GEOMETRY-FIRST (``rigid=True``, DELFIN_FFFREE_LIGAND_RIGID): for a RIGID
    polydentate (clathrochelate cage cap / macrocycle / conjugated pincer) the
    ligand's OWN internal geometry is the hard constraint and the coordination
    polyhedron is EMERGENT.  The historic per-donor radial rescale moves each donor
    INDEPENDENTLY onto its exact ideal radius, which — on a rigid backbone — splays
    the ring open (eye-find QEBLOC: ring N-N-N 110deg->129deg, B-N 1.54->1.52,
    N-N compressed) because the donors are forced apart while the backbone cannot
    follow.  Instead apply a SINGLE UNIFORM scale (metal at origin) so EVERY M-X
    distance scales by the same factor: all bond ANGLES (ligand-internal AND the
    coordination bite) are preserved EXACTLY, donors land at ~the target radius, and
    the coordination polyhedron comes out slightly twisted but PHYSICALLY REAL
    (the cage cavity dictates it).  This reverses the construction order from
    "coordination ideal -> ligand distorts" to "ligand geometry first -> coordination
    emergent".  Default ``rigid=False`` -> byte-identical (per-donor rescale).

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

    # TRILATERATED TARGETS (DELFIN_FFREE_TRILATERATE=1, default OFF -> byte-identical).
    #
    # The two branches below force a choice: per-donor radial gives the right M-D lengths
    # and splays the ring, rigid keeps the ring and lets M-D come out emergent.  The choice
    # only exists because the TARGETS are ideal polyhedron vertices, which the ligand
    # generally cannot reach.  Move the targets first -- onto points that carry both the
    # M-D radii and the ligand's own donor separations -- and the rigid fit lands on them
    # with a small residual, so both hold at once.  Bidentate already does exactly this via
    # _bite_aware_targets (the cosine law); this is the same statement for any denticity,
    # and it is the only place polydentates have ever had it: _bite_aware_targets is called
    # from _place_chelate_block alone, which runs only for dent == 2.
    if _trilat_targets_on():
        try:
            _tri = _trilaterate_donor_targets(
                lP, list(donor_idxs), [targets[perm[i]] for i in range(len(donor_idxs))])
        except Exception:
            _tri = None
        if _tri is not None:
            _don = np.array([np.asarray(lP[d], float) for d in donor_idxs])
            _tgt = np.array(_tri, float)
            _dmu, _tmu = _don.mean(0), _tgt.mean(0)
            try:
                _Rk = _kabsch_rot(_don - _dmu, _tgt - _tmu)
                _Qt = (np.asarray(lP, float) - _dmu) @ _Rk.T + _tmu
                if np.all(np.isfinite(_Qt)):
                    return _Qt
            except Exception:
                pass

    if rigid:
        # LIGAND-GEOMETRY-FIRST: rigid-body best-fit (rotation + translation, NO scale,
        # NO per-donor radial move) of the donor set onto the target POINTS.  Both the
        # ligand-internal geometry (bonds AND angles) AND the coordination bite are
        # preserved EXACTLY as embedded; the metal-donor distances come out EMERGENT
        # (the rigid cage cavity dictates them).  A uniform scale was rejected here:
        # it preserves angles but, because the metallacycle embed's M-D runs short
        # (~1.95 vs 2.13), the M-D correction factor inflates EVERY ligand bond
        # (B-N 1.54->1.66).  The pure rigid fit keeps the embedded bonds untouched.
        don = np.array([np.asarray(lP[d], float) for d in donor_idxs])
        tgt = np.array([np.asarray(targets[perm[i]], float)
                        for i in range(len(donor_idxs))])
        dmu = don.mean(0); tmu = tgt.mean(0)
        Rk = _kabsch_rot(don - dmu, tgt - tmu)
        Qr = (np.asarray(lP, float) - dmu) @ Rk.T + tmu
        # BETA IN THE SEATING (DELFIN_FFFREE_DONOR_PLANE, default OFF -> byte-identical).
        #
        # beta = angle between the M->D bond and the plane of the donor's OWN conjugated
        # environment (0 deg = metal in plane).  Measured against crystals:
        #     monodentate 0.9 deg (BETTER than the 3.4 crystals allow)
        #     bidentate   9.3      tetradentate 15.9   (crystals: 3.7 / 4.8)
        # The monodentate path is right because it aligns the metal along the lone pair
        # (_donor_and_lp + _rot_align).  THIS path never looks at the donor plane at all --
        # it fits donors rigidly onto ideal vertices, and beta comes out as a by-product.
        # The break is a PATH CHANGE, not physics.  pyramidal_sp2 is our largest INDEPENDENT
        # defect (15.7 % of our frames vs 1.0 % of crystals; sole hard finding on 37.9 % of
        # its firings), and trilatresc2 lost its verdict to exactly this axis (OVAVEO,
        # pyramid_frame_regressed + root_defects_increased).
        #
        # THE FREE DEGREE OF FREEDOM: a rigid rotation about the axis through the donor
        # centroid changes NOTHING that is already right -- every bond, every internal
        # angle and the coordination BITE are preserved exactly (it is a rigid motion of
        # the whole ligand, and the bite is an internal distance).  It only trades vertex
        # alignment, which the polyhedron never had exactly anyway.  So beta can be reduced
        # at zero cost to the two quantities we already get right.
        #
        # Loss-guarded: the rotation is kept ONLY if it reduces the summed beta AND the
        # donor-to-target residual does not grow beyond _DP_RESID_SLACK.  Otherwise the
        # rigid fit stands unchanged -- a candidate that helps neither is never taken.
        if os.environ.get("DELFIN_FFFREE_DONOR_PLANE", "0") == "1" and lsyms is not None:
            _q = _donor_plane_relax(Qr, lsyms, list(donor_idxs), tgt, tmu)
            if _q is not None:
                Qr = _q
        # CAGE-CRUSH GUARD (DELFIN_FFFREE_CAGE_MD_GUARD, default-OFF -> byte-id):
        # the rigid fit preserves the embed's ligand geometry EXACTLY, so it also
        # preserves a CRUSHED embed (ETKDG produced a cavity far too small for the
        # metal) -> the emergent M-D comes out physically IMPOSSIBLE (eye/CCDC:
        # clathrochelate Co-O 1.49 vs ideal ~1.95, a 24% collapse).  A crushed
        # coordination bond is UNRELAXABLE downstream (breaks topology / QM), whereas
        # the per-donor radial placement gives correct M-D at the cost of (soft,
        # UFF-relaxable) angle splay.  So when the emergent coordination sphere is
        # crushed below a physical fraction of ideal, the "ligand-geometry-first"
        # premise has FAILED for this embed -> DON'T trust it: fall through to the
        # per-donor radial placement (ideal M-D).  Healthy rigid cages (emergent M-D
        # ~0.9 of ideal, e.g. QEBLOC) stay on the rigid path untouched.  General:
        # keys only on the physical M-D fraction, never on any SMILES/refcode.
        if os.environ.get("DELFIN_FFFREE_CAGE_MD_GUARD", "0") == "1":
            try:
                _frac = float(os.environ.get("DELFIN_FFFREE_CAGE_MD_FRAC", "0.85"))
            except Exception:
                _frac = 0.85
            # PER-DONOR IN-PLACE reset (NOT mean-vs-mean, NOT fall-through-to-rescale):
            # the crush is typically UNEVEN — a subset of donors (the 3 O of an N3O3
            # clathrochelate, or the 2 short Cd-N of an over-coordinated diimine)
            # collapse below ideal while the rest sit near ideal.  For EACH donor whose
            # emergent M-D (|Qr[d]|, metal at origin) is crushed below _frac of ITS OWN
            # ideal M-D, reset it radially to ideal (keep direction); healthy donors +
            # the backbone stay on the rigid fit.  This fixes the EMITTED frame DIRECTLY
            # — a fall-through to the full per-donor rescale changes conformer SELECTION
            # and the clash metric then prefers the still-crushed conformer (YECSUW:
            # emitted stayed 1.785).  The downstream constrained relax (donors frozen)
            # pulls the backbone to follow.  General, keys only on the physical M-D
            # ratio (never SMILES/refcode); byte-id OFF; healthy cages (ratio ~0.9,
            # QEBLOC/FEKZON) have no donor below threshold -> untouched.
            _nfix = 0
            _hf = (lsyms is not None
                   and os.environ.get("DELFIN_FFREE_H_FOLLOW", "0") == "1")
            for _i, _d in enumerate(donor_idxs):
                _r = float(np.linalg.norm(Qr[_d]))
                _idl = float(tgt_md[perm[_i]])
                if _r > 1e-6 and _idl > 1e-6 and _r < _frac * _idl:
                    _nv = Qr[_d] / _r * _idl
                    if _hf:                      # hydrogens ride with their parent, as above
                        _rd = _hydrogens_riding_on(lsyms, Qr, _d)
                        _dl = _nv - Qr[_d]
                        for _h in _riders_that_may_move(lsyms, Qr, _rd, _dl, _d):
                            Qr[_h] = Qr[_h] + _dl
                    Qr[_d] = _nv
                    _nfix += 1
            if os.environ.get("DELFIN_CAGE_DEBUG", "0") == "1" and _nfix:
                os.write(2, ("[CAGE_GUARD] dent=%d reset %d/%d donors to ideal M-D\n"
                             % (len(donor_idxs), _nfix, len(donor_idxs))).encode())
            return Qr
        else:
            return Qr
    # Per-donor RADIAL placement at the exact ideal M-donor distance (NOT a uniform
    # scale, which preserved the ETKDG embed's M-D asymmetry -> over-contracted donors,
    # the FEKZON CCDC defect).  Keep each donor's Kabsch-rotated DIRECTION (so the embed's
    # natural bite angle is preserved) and set only its radius to md.  The constrained
    # relax (donors fixed here) then pulls the backbone into consistency.
    # HYDROGENS RIDE WITH THEIR PARENT (DELFIN_FFREE_H_FOLLOW=1, default OFF -> byte-id).
    #
    # The loop below moves the DONOR ATOM and nothing else.  A donor that carries hydrogens
    # -- an amine N-H, a hydroxyl O-H, an agostic C-H -- therefore has its heavy atom
    # displaced radially while its H stay where the embed put them, which corrupts both the
    # X-H length and its direction by exactly the displacement.
    #
    # MEASURED, reference-free (weddell/tools/h_geometry.py over 296696 X-H bonds of the
    # champion archive): the MEDIAN X-H length is textbook-correct everywhere -- aromatic
    # C-H 1.080, methyl 1.109, N-H 1.034, O-H 0.990 -- so the placement RULE is right and
    # must not be touched.  What is wrong is a TAIL: C|2|2|1 has p10 = 0.948 A, C|1|3|1
    # spans 0.886 to 1.237, and the direction deviation reaches p90 = 58 deg where its own
    # median is 6.  A rule that produced correct medians does not produce that tail; being
    # left behind by a later move does.
    #
    # (The CCDC comparison cannot referee this and was nearly a trap: crystal C-H sit at
    # p10/p50/p90 = 0.930/0.949/0.960 with a 0.45 deg direction width, i.e. a RIDING MODEL,
    # not a measurement.  "Fixing" our lengths towards it would have broken correct H.)
    #
    # This is additive in the strict sense: it invents no geometry, it preserves the X-H
    # geometry the embed already had.  Nothing about heavy-atom placement changes.
    _hfollow = (lsyms is not None
                and os.environ.get("DELFIN_FFREE_H_FOLLOW", "0") == "1")
    # LET THE NEIGHBOURHOOD FOLLOW (DELFIN_FFFREE_DONOR_FOLLOW, default OFF -> byte-identical).
    # The weights are read off the geometry BEFORE any donor moves, so the graph is the one
    # the embed produced; see _donor_follow_weights for why a decay and not a repair.
    _dfollow = None
    if lsyms is not None and os.environ.get("DELFIN_FFFREE_DONOR_FOLLOW", "0") == "1":
        try:
            _span = max(2, int(os.environ.get("DELFIN_FFFREE_DONOR_FOLLOW_SPAN", "3")))
            _dfollow = _donor_follow_weights(lsyms, Q, list(donor_idxs), _span)
        except Exception:
            _dfollow = None
    _ddelta = {}
    for i, di in enumerate(donor_idxs):
        r = float(np.linalg.norm(Q[di]))
        if r > 1e-6:
            _new = Q[di] / r * tgt_md[perm[i]]
            if _hfollow:
                _riders = _hydrogens_riding_on(lsyms, Q, di)   # BEFORE the move
                _delta = _new - Q[di]
                for _h in _riders_that_may_move(lsyms, Q, _riders, _delta, di):
                    Q[_h] = Q[_h] + _delta
            _ddelta[int(di)] = _new - Q[di]
            Q[di] = _new
    if _dfollow:
        # ONLY WHERE THE PLACEMENT ACTUALLY BROKE SOMETHING (default; the global form is
        # DELFIN_FFFREE_DONOR_FOLLOW_ALWAYS=1).
        #
        # Measured on 187 systems, applied to EVERY seating: reach 92, capability +13 and
        # valid 56 -> 66 -- the largest capability gain of the day -- but cap_LOST 3 and
        # sixteen red terms (pyramid_frame_regressed 23, smiles_ccdc_regressed 13,
        # isomers_lost 7).  The root is right and the scope was wrong: it also moved the
        # backbones of frames that were already clean, and those had everything to lose.
        #
        # A frame that already carries a collapsed bond has nothing to lose, so restricting
        # the decay to exactly those frames cannot cost a capability by construction -- the
        # same argument that every lever which landed here rests on.  Where the placement was
        # clean, this is byte-identical.
        # KEEP IT ONLY WHERE THE GRADED BOND AXIS STRICTLY IMPROVES.
        #
        # A first attempt gated on "the frame already collapsed" and was wrong TWICE, both
        # measured on 187 systems:
        #   * Reach fell from 92 to 3.  The displacement is about 0.18 A on a C-N bond of
        #     1.47 A, i.e. 0.88 x ideal -- ABOVE the 0.82 collapse threshold, so the gate
        #     almost never fired.  The decay PREVENTS a collapse; it does not repair one.
        #   * On those 3 it still did damage (isomers_lost 2, ccdc_arrangement_lost 2), so
        #     "a collapsed frame has nothing to lose" is false at FRAME level -- that argument
        #     holds for a system with broken_frac 1.0, not for one frame among many, which can
        #     still be the only carrier of an isomer.
        #
        # The graded quantity is the one the eye reads (org_bond), so that is what decides:
        # keep the decayed frame only if the WORST relative bond deviation strictly drops.
        # Where it does not, the frame is left exactly as it was -- never-worse on the axis
        # this lever exists to improve, by construction rather than by hope.
        _always = os.environ.get("DELFIN_FFFREE_DONOR_FOLLOW_ALWAYS", "0") == "1"

        def _worst_bond_dev(_s, _P):
            _P = np.asarray(_P, float)
            _w = 0.0
            for _i in range(len(_s)):
                if _s[_i] == "H" or _bd._is_metal(_s[_i]):
                    continue
                for _j in range(_i + 1, len(_s)):
                    if _s[_j] == "H" or _bd._is_metal(_s[_j]):
                        continue
                    _id = _bd._ideal_bond(_s[_i], _s[_j])
                    _dd = float(np.linalg.norm(_P[_i] - _P[_j]))
                    if _id <= 0 or _dd > 1.30 * _id:
                        continue
                    _dv = abs(_dd - _id) / _id
                    if _dv > _w:
                        _w = _dv
            return _w

        _Qf = np.array(Q, float)
        for _a, (_own, _w) in _dfollow.items():
            _d = _ddelta.get(int(_own))
            if _d is not None:
                _Qf[_a] = _Qf[_a] + _w * _d
        try:
            if _always or _worst_bond_dev(lsyms, _Qf) < _worst_bond_dev(lsyms, Q) - 1e-9:
                Q = _Qf
        except Exception:
            pass
    # BETA IN THE SETTING -- ON THE PATH THAT ACTUALLY RUNS.
    #
    # The first two attempts hooked this into the `if rigid:` branch above and had ZERO
    # reach on 187 systems, twice (donorplane, donorplane2, both rc=3).  Cause, found only
    # after the second refusal: line ~3654 sets
    #     _rigid_seat = dent >= 3 and (LIGAND_RIGID or RIGID_LIGAND_SEAT)
    # and BOTH of those flags are dark.  The rigid branch never executes in the champion, so
    # the lever was not in the wrong path -- it was in NO path.  Third case of that class in
    # one day (POLY6 sat in the parked functional), hence the rule: before building into a
    # branch, grep the ENCLOSING CONDITION, not just the function.
    #
    # This is the live default path: every donor has just been reset radially to its ideal
    # M-D length, so r(M-D) is exactly right and beta is whatever the embed happened to give.
    # A rotation about the line through the donors leaves the donors on that line -- for a
    # bidentate the axis IS the donor-donor separation, i.e. the bite, so the bite and the
    # just-corrected M-D lengths both survive untouched.  Only the backbone swings, and with
    # it the donor planes.  Loss-guarded inside _donor_plane_relax.
    if os.environ.get("DELFIN_FFFREE_DONOR_PLANE", "0") == "1" and lsyms is not None:
        try:
            _tg = np.array([np.asarray(targets[perm[i]], float)
                            for i in range(len(donor_idxs))], float)
            _q = _donor_plane_relax(Q, lsyms, list(donor_idxs), _tg, _tg.mean(0))
            if _q is not None:
                Q = _q
        except Exception:
            pass
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


def _straighten_sp_chain(lsyms, lP, lmol, di, flag="DELFIN_FFREE_SP_LINEAR"):
    """An SP centre in the donor's substituent chain is LINEAR.  Make it so.

    (DELFIN_FFREE_SP_LINEAR=1, default OFF -> byte-identical.)

    MEASURED, and the measurement is what located it.  smiles_sp-not-linear fires on 437
    findings of the champion archive, 148 of them on frames where it is the ONLY hard finding,
    and on 1500 clean CCDC crystals it fires ZERO times.  Its distance to the nearest metal is
    a razor-thin band -- p10/p50/p90 = 2.89 / 3.10 / 3.36 A, 0.2 % beyond 5 A -- i.e. always
    exactly ONE BOND beyond the coordination sphere.  Every case is a cumulated pseudohalide,
    M-N=C=S / M-N=C=Se (UJAZUD02, JEJSID, ADITIT, QUHWAT, QINJAB), built at 105-116 deg where
    the SMILES itself says 180.

    WHY THE EXISTING MECHANISM MISSES IT -- a scope, not a blindness.  _vsepr_reconstruct
    reads the hybridisation OF THE DONOR and sets the M-D-substituent angle; for these ligands
    the donor is the N (correctly placed linear, k==1) and the sp atom is the CARBON one bond
    further out, which nothing in that function ever touches.  67-90 % of each affected
    system's frames carry it, so the geometry IS reachable -- the build simply does not
    insist on it.

    No table, no fit, no crystal reference: a centre carrying a triple bond, or two double
    bonds, is linear by definition, and the input graph already states the bond orders.  The
    far subtree is rotated rigidly about the sp atom, so every bond LENGTH and every angle
    inside that subtree is preserved exactly -- only the one angle that was wrong changes.
    """
    if os.environ.get(flag, "0") != "1":
        return lP
    try:
        P = np.array(lP, float).copy()
        for a in lmol.GetAtomWithIdx(int(di)).GetNeighbors():
            ai = int(a.GetIdx())
            nb = [n.GetIdx() for n in a.GetNeighbors()]
            if len(nb) != 2 or a.IsInRing():
                continue
            orders = sorted(float(b.GetBondTypeAsDouble()) for b in a.GetBonds())
            # sp iff a triple bond, or two doubles (a cumulene) -- read off the graph.
            if not (orders[-1] >= 2.9 or (len(orders) == 2 and orders[0] >= 1.9
                                          and orders[1] >= 1.9)):
                continue
            far = int(nb[0]) if int(nb[1]) == int(di) else int(nb[1])
            v1 = P[int(di)] - P[ai]; v2 = P[far] - P[ai]
            n1 = float(np.linalg.norm(v1)); n2 = float(np.linalg.norm(v2))
            if n1 < 1e-6 or n2 < 1e-6:
                continue
            c = float(np.dot(v1 / n1, v2 / n2))
            ang = math.degrees(math.acos(max(-1.0, min(1.0, c))))
            if ang > 170.0:
                continue                       # already linear
            axis = np.cross(v2, v1)
            na = float(np.linalg.norm(axis))
            if na < 1e-6:
                continue
            # rotate the FAR subtree (never the donor side) onto the straight continuation
            grp = _subtree(lmol, far, ai)
            R = _axis_rot(axis / na, math.radians(180.0 - ang))
            for g in grp:
                P[g] = (P[g] - P[ai]) @ R.T + P[ai]
        return P
    except Exception:
        return lP


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
    lP = _straighten_sp_chain(lsyms, lP, lmol, di)
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


def _trilaterate_donor_targets(lP, donor_idxs, targets, iters=400, damp=0.5):
    """Donor targets that honour the ligand's OWN donor-donor distances AND the M-D radii.

    THE GENERALISATION.  ``_bite_aware_targets`` solves exactly this for two donors, with
    the cosine law: given the two M-D radii and the donor separation the ligand actually
    has, the angle between them FOLLOWS.  For k donors the same statement is a
    trilateration --

        |x_i| = r_i                    (the measured / referenced M-D length)
        |x_i - x_j| = D_ij             (the separation the LIGAND already has)

    -- and the coordination angles are again a CONSEQUENCE, never an input.  No ideal
    polyhedron appears anywhere in it, no metal is named, and k = 2 reduces to the cosine
    law, so bidentate and macrocycle are one rule.

    WHY IT MATTERS HERE.  This seating documents a trade-off it treats as unavoidable: the
    per-donor radial placement gives correct M-D lengths at the price of splayed ring
    angles, while the rigid-body fit keeps the ligand exact but lets M-D come out emergent
    and sometimes impossible (a clathrochelate Co-O at 1.49 against ~1.95, a 24 % collapse).
    Trilateration is the resolution: it satisfies BOTH as far as geometry allows, and where
    they are genuinely incompatible the residual states by how much instead of silently
    picking a side.  Measured motivation: over 995 built systems, 29 % of chelate bites come
    out beyond 92 deg, a region real chelates essentially do not occupy (49 measured bins,
    98 % below 90, none above 95).

    Solved by alternating projection -- project every donor back onto its own radius, then
    correct each pair toward the ligand's separation, damped, repeat.  Deterministic, no
    random start, no minimiser.  It BEGINS at the enumerated vertex targets, so donor i
    stays in the region of vertex i and the isomer assignment the enumeration made is
    preserved; this only moves the targets to where the ligand can actually reach them.
    """
    k = len(donor_idxs)
    if k < 2:
        return None
    P = np.array([np.asarray(targets[i], float) for i in range(k)])
    r = np.array([float(np.linalg.norm(P[i])) for i in range(k)])
    if np.any(r < 1.0e-6):
        return None
    L = np.array([np.asarray(lP[d], float) for d in donor_idxs])
    D = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            D[i, j] = D[j, i] = float(np.linalg.norm(L[i] - L[j]))
    if not np.all(np.isfinite(D)):
        return None
    for _ in range(int(iters)):
        # 1) back onto each donor's own sphere around the metal
        for i in range(k):
            n = float(np.linalg.norm(P[i]))
            if n > 1.0e-9:
                P[i] *= r[i] / n
        # 2) pull/push every pair toward the separation the ligand actually has
        for i in range(k):
            for j in range(i + 1, k):
                v = P[i] - P[j]
                d = float(np.linalg.norm(v))
                if d < 1.0e-9 or D[i, j] < 1.0e-9:
                    continue
                corr = damp * 0.5 * (d - D[i, j]) * (v / d)
                P[i] -= corr
                P[j] += corr
    for i in range(k):                      # radii are the hard side: end on them
        n = float(np.linalg.norm(P[i]))
        if n > 1.0e-9:
            P[i] *= r[i] / n
    if not np.all(np.isfinite(P)):
        return None
    return [P[i] for i in range(k)]


def _lone_pair_dir(P, mol, idx):
    """Direction of the donor's lone pair, from its own bonds only.

    Minus the sum of the unit vectors to ALL its neighbours -- hydrogens included, because
    for an sp3 amine the two N-H bonds are what fix where the lone pair points; dropping
    them would leave only the C-N axis and give the wrong direction entirely.

    One expression covers every case, which is the test this project applies to any rule:
    an aromatic pyridine N (two ring neighbours) gets the in-plane outward bisector, an sp3
    amine (three neighbours) gets the fourth tetrahedral direction, a linear nitrile N (one
    neighbour) gets the bond axis reversed.  No hybridisation label is read, no functional
    group is named, and it reads the same for a donor nobody has classified yet.
    """
    v = np.zeros(3)
    try:
        for nb in mol.GetAtomWithIdx(int(idx)).GetNeighbors():
            d = P[int(nb.GetIdx())] - P[int(idx)]
            n = float(np.linalg.norm(d))
            if n > 1.0e-9:
                v += d / n
    except Exception:
        return None
    n = float(np.linalg.norm(v))
    if n < 1.0e-6:
        return None                      # symmetric surroundings: no defined direction
    return -v / n


def _lp_orient_seated_bidentate(Q, mol, d1, d2, lsyms=None):
    """Set the ONE rotational freedom a SEATED bidentate still has, from the lone pairs.

    WHERE THE TILT ACTUALLY COMES FROM.  The live path seats a chelate with
    _orient_chelate_to_vertices, a Kabsch fit of the donor directions onto the target vertex
    directions.  With TWO donors that fit is UNDERDETERMINED: rotating the ligand about the
    donor-donor axis leaves both donors exactly where they are, so it changes nothing the fit
    can see, and whatever the SVD happens to return decides it.  That leftover rotation is
    precisely the one that decides whether the metal lies in a conjugated donor's pi plane --
    i.e. beta, the largest single realism gap this project has (18-19 % of frames against
    0.98 % of crystals).  It was being set by a numerical tie-break.

    The law that fixes it already existed and was already measured -- over 995 built systems
    the tilt tracks co-ligand contact (median 12.1 deg where vdW shells overlap, 0.9 deg where
    they do not) -- but it sat in _place_chelate_block, which assemble_from_config reaches ONLY
    as a fallback after the embed has already failed.  Three independent flags in that function
    (LP_ORIENT, BITE_LAW, BITE_MEASURED) all measured ZERO reach on a 187-system sweep, while
    flags elsewhere in the same file measured 123 and 7 -- the function is live code on a path
    the champion does not take.

    SAFE BY CONSTRUCTION for the seating: both donors lie ON the rotation axis, so they do not
    move at all.  The M-D distances, the bite angle and the arm-to-vertex correspondence are
    all preserved EXACTLY; only the backbone turns.  Returns None when no lone-pair direction
    is defined, or when the turn would introduce a collapsed bond the seated frame did not
    have -- the caller then keeps what it had.
    """
    Q = np.asarray(Q, float)
    axis = Q[d2] - Q[d1]
    if float(np.linalg.norm(axis)) < 1.0e-9:
        return None
    mid = 0.5 * (Q[d1] + Q[d2])
    th = _lp_aligned_angle(Q, mol, d1, d2, axis, mid)
    if th is None:
        return None
    Qr = (Q - mid) @ _axis_rot(axis, th).T + mid
    if not np.all(np.isfinite(Qr)):
        return None
    # NO STERIC ESCAPE HERE, AND THAT IS THE POINT.  A first version of this rejected the
    # derived orientation whenever it brought a backbone atom near the metal, i.e. it fell
    # back to the old sweep exactly where the sweep and the law disagree -- and measured
    # affected=0 on 187 systems, changing nothing at all.  The docstring of _lp_aligned_angle
    # records the same experiment and the same outcome one layer down: "it is why the first
    # version changed almost nothing".  Twice is enough.  The M-D sigma bond goes THROUGH the
    # lone pair, so the plane is a BONDING constraint and a clash is not a licence to break
    # it; a real molecule relieves the clash by turning substituents, opening the bite or
    # twisting the polyhedron.  If a collision remains it now STAYS VISIBLE for the self-gate
    # and the clash axis to report, instead of being hidden inside a tilt.
    return Qr


def _lp_aligned_angle(Q, mol, d1, d2, axis, mid):
    """The rotation about the donor-donor axis that points BOTH lone pairs at the metal.

    THE DEFECT THIS REPLACES.  The orientation used to be picked by a 36-step sweep that
    maximised the minimum metal-distance of the backbone -- collision avoidance, i.e. a
    STERIC heuristic where the constraint is a BONDING one.  Three consequences, and all
    three look like the flapping the manifold is full of: the criterion is simply the wrong
    one, so a planar chelate's plane tilts out; 36 steps quantise to 10 deg; and scoring by
    the MINIMUM over backbone atoms lets a single atom decide, so a small change in the
    ligand flips the chosen step and the orientation jumps between similar systems.

    A conjugated donor's lone pair lies IN its pi plane, so the metal lies in that plane
    too -- there is no such thing as a tilted aromatic chelate.  That is a hard geometric
    condition, not a preference, and it has a closed form.  Both donors sit ON the rotation
    axis, so they do not move and their target directions t_i are constant; each lone pair
    rotates by Rodrigues, and the total alignment is

        sum_i lp_i(theta) . t_i  =  A cos(theta) + B sin(theta) + C

    maximised at ``theta = atan2(B, A)``.  One atan2 instead of 36 samples: exact, no
    quantisation, and no discontinuity to jump across.

    Returns ``None`` when no lone-pair direction is defined at either donor, so the caller
    keeps the old sweep rather than inventing an orientation.
    """
    k = axis / float(np.linalg.norm(axis))
    A = B = C = 0.0
    used = 0
    for d in (d1, d2):
        v = _lone_pair_dir(Q, mol, d)
        if v is None:
            continue
        t = mid * 0.0 - Q[d]             # metal sits at the origin
        nt = float(np.linalg.norm(t))
        if nt < 1.0e-9:
            continue
        t = t / nt
        kv, kt = float(np.dot(k, v)), float(np.dot(k, t))
        A += float(np.dot(v, t)) - kv * kt
        B += float(np.dot(np.cross(k, v), t))
        C += kv * kt
        used += 1
    if used == 0 or (abs(A) < 1.0e-12 and abs(B) < 1.0e-12):
        return None
    return math.atan2(B, A)


def _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=None):
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

    # LONE-PAIR ORIENTATION (DELFIN_FFREE_LP_ORIENT, default OFF -> byte-identical).
    #
    # THE PLANE WINS, WITHOUT A STERIC ESCAPE.  A first version kept a floor that fell back
    # to the sweep when the derived orientation pushed a backbone atom near the metal.  That
    # is the wrong hierarchy and it is why the first version changed almost nothing: the
    # M-D sigma bond goes THROUGH the lone pair, so a conjugated donor holds the metal in
    # its pi plane, and a real molecule relieves a clash by turning substituents, opening
    # the bite or twisting the polyhedron -- never by tipping the metal out of the plane,
    # which would destroy the overlap that binds it.  Measured over 995 built systems: the
    # tilt tracks the co-ligand contact (Pearson -0.20; median tilt 12.1 deg where the vdW
    # shells already overlap, 0.9 deg where they do not), i.e. the seat is spending the
    # ONE degree of freedom it has on steric relief and paying for it with the bond.
    #
    # So the plane is set, unconditionally.  If a collision remains it now STAYS VISIBLE
    # for the clash axis to report, instead of being hidden inside a tilt -- and it has to
    # be resolved where reality resolves it (cis-edge choice, conformer), not here.
    if mol is not None and os.environ.get("DELFIN_FFREE_LP_ORIENT", "0") == "1":
        _th = _lp_aligned_angle(Q, mol, d1, d2, axis, mid)
        if _th is not None:
            return (Q - mid) @ _axis_rot(axis, _th).T + mid

    best, bestQ = -1e9, Q
    for k in range(36):
        Qr = (Q - mid) @ _axis_rot(axis, 2 * np.pi * k / 36).T + mid
        score = min(float(np.linalg.norm(Qr[i])) for i in body_idx)   # closest backbone atom to metal
        if score > best:
            best, bestQ = score, Qr
    return bestQ


_BITE_MAX_RING = 7      # above this a "shared ring" is a macrocycle, not a bite


def bite_close(u, v, theta_deg):
    """Rotate two unit directions symmetrically in their own plane until the angle between
    them is ``theta_deg``, keeping their bisector fixed.

    WHY.  An ideal octahedron puts two cis vertices 90 deg apart, and a five-membered
    chelate ring measured over real crystals bites at 70.6 deg (p10 69.4, p90 72.6 -- a
    3 deg wide band, one of the best defined quantities we have).  So every five-ring
    chelate is seated 19.4 deg too open, and chelates are most of our systems.  This is the
    quantified form of the earthquake root, "rigid polydentates forced onto ideal vertices",
    and of poly_cshm_vs_ccdc being the worst axis in all seven UFF-constraint A/Bs.

    Keeping the bisector fixed matters: the chelate must still occupy the SAME edge of the
    polyhedron that the isomer enumeration assigned to it.  Only the opening changes, so
    the isomer identity is untouched and this stays a pure seating correction.
    """
    u = np.asarray(u, float); v = np.asarray(v, float)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-9 or nv < 1e-9:
        return u, v
    u = u / nu; v = v / nv
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    cur = math.degrees(math.acos(c))
    if abs(cur - theta_deg) < 1e-6:
        return u, v
    bis = u + v
    nb = np.linalg.norm(bis)
    if nb < 1e-9:                       # exactly opposed: no bisector, nothing to close on
        return u, v
    bis = bis / nb
    perp = u - float(np.dot(u, bis)) * bis
    np_ = np.linalg.norm(perp)
    if np_ < 1e-9:
        return u, v
    perp = perp / np_
    half = math.radians(theta_deg) / 2.0
    return (math.cos(half) * bis + math.sin(half) * perp,
            math.cos(half) * bis - math.sin(half) * perp)


def _chelate_ring_size(lmol, d1, d2):
    """Size of the chelate ring the two donors would close WITH the metal.

    The ligand mol has no metal, so the donors are not yet in a common ring: the chelate
    ring is the shortest donor-to-donor path through the ligand, plus the metal itself.
    Path of k bonds -> ring of k+2 atoms (k+1 ligand atoms and the metal).
    """
    try:
        from rdkit import Chem
        path = Chem.GetShortestPath(lmol, int(d1), int(d2))
        return (len(path) + 1) if path else 0
    except Exception:
        return 0


def _measured_bite(metal, cn, ring_size, e1, e2):
    """Measured p50 bite angle (deg) for this chelate, or None.

    Only for SMALL rings: shared rings of 11-14 come out 68-82 deg WIDE, because two
    donors of a macrocycle can sit cis or trans while still sharing a ring.  Beyond
    _BITE_MAX_RING there is no single bite to aim at, and forcing one would be the
    bimodal-p50 mistake this table was designed to avoid.
    """
    if not ring_size or ring_size > _BITE_MAX_RING:
        return None
    path = os.environ.get("DELFIN_FFREE_DMD_BANDS", "")
    if not path:
        return None
    try:
        from delfin.manta.polyhedra import _load_band_tsv
        tbl = _load_band_tsv(path)
    except Exception:
        return None
    ea, eb = sorted((e1, e2))
    for k in (f"{metal},{cn}|{ea},{eb}|{ring_size}",
              f"{metal},{cn}|*|{ring_size}",
              f"{metal}|*|{ring_size}",
              f"*,{cn}|*|{ring_size}"):
        row = tbl.get(k)
        if row is not None:
            return float(row[1])                     # p50
    return None


def assemble_multichelate(metal: str, geometry: str, chelate_specs):
    """chelate_specs: list of (smiles, [d1,d2], [v1,v2]).  Place several chelating
    ligands on their assigned cis vertex-pairs (e.g. tris-chelate [M(en)3])."""
    ref = MSB._ref_vectors(geometry)
    out_syms = [metal]; blocks = [np.zeros((1, 3))]
    for smi, dons, verts in chelate_specs:
        lsyms, lP, lmol = _ligand_3d(smi)
        d1, d2 = dons
        u1 = ref[verts[0]] / np.linalg.norm(ref[verts[0]])
        u2 = ref[verts[1]] / np.linalg.norm(ref[verts[1]])
        # MEASURED BITE (DELFIN_FFREE_BITE_MEASURED, default OFF -> byte-identical).
        # The ideal polyhedron separates two cis vertices by 90 deg; a five-membered
        # chelate measured over real crystals bites at 70.6 deg with a 3 deg wide band.
        # Closing the pair onto the measured angle -- bisector fixed, so the ligand keeps
        # the very edge the isomer enumeration assigned it -- is a pure SETTING: no weight,
        # no gate, no minimiser.
        _r1 = MSB.md_distance(metal, lsyms[d1],
                              atom=lmol.GetAtomWithIdx(d1), mol=lmol)
        _r2 = MSB.md_distance(metal, lsyms[d2],
                              atom=lmol.GetAtomWithIdx(d2), mol=lmol)
        # DERIVED BITE (DELFIN_FFREE_BITE_LAW, default OFF -> byte-identical).
        #
        # A chelate closes a ring THROUGH the metal, so the bite is not an independent
        # quantity at all -- it is ring closure:
        #
        #     cos(bite) = (r1^2 + r2^2 - d_DD^2) / (2 r1 r2)
        #
        # with d_DD the donor-donor distance the LIGAND ALREADY HAS, read straight off its
        # own FF-free 3D construction.  Verified against 131769 crystals without building
        # anything: inverting this on the measured D-M-D table returns chemically exact
        # donor-donor distances (1.402 A for an eta2 C=C, 2.614 A for an en/bipy 5-ring)
        # and reproduces the measured bite to 0.24 deg across 18 metals.  Nine of eleven
        # ligand types sit at or below the precision of the M-D input itself.
        #
        # WHY THIS AND NOT THE TABLE.  Not precision -- the error in r enters either way.
        # COUPLING.  mdAB measured what happens when they are set independently: the M-D
        # length was corrected while the vertex angle stayed at the ideal 90 deg, which
        # forces d_DD = |u1*r1 - u2*r2| to a value the ligand DOES NOT HAVE, so the ligand
        # absorbs the difference (valid 155 -> 146).  Deriving the bite from d_DD and r
        # makes that contradiction impossible by construction: the ligand keeps its own
        # donor separation, whatever r turns out to be.
        #
        # It also extrapolates where the table has holes (bins below n=200 -- precisely the
        # exotic systems carrying our defects), separates cis from trans by itself in a
        # macrocycle (they simply ARE two different d_DD), and is a formula rather than a
        # dataset.
        if os.environ.get("DELFIN_FFREE_BITE_LAW", "0") == "1":
            _dd = float(np.linalg.norm(lP[d1] - lP[d2]))
            if _dd > 1.0e-6 and _r1 > 1.0e-6 and _r2 > 1.0e-6:
                _c = (_r1 * _r1 + _r2 * _r2 - _dd * _dd) / (2.0 * _r1 * _r2)
                # Outside [-1, 1] the triangle does not close: the ligand's own donor
                # separation cannot be reached at these M-D lengths.  Forcing a clamped
                # angle there would silently invent a geometry, so leave the vertex
                # direction alone and let the existing path handle it.
                if -1.0 <= _c <= 1.0:
                    u1, u2 = bite_close(u1, u2, math.degrees(math.acos(_c)))
        elif os.environ.get("DELFIN_FFREE_BITE_MEASURED", "0") == "1":
            # MEASURED BITE.  The ideal polyhedron separates two cis vertices by 90 deg; a
            # five-membered chelate measured over real crystals bites at 70.6 deg with a
            # 3 deg wide band.  Closing the pair onto the measured angle -- bisector fixed,
            # so the ligand keeps the very edge the isomer enumeration assigned it -- is a
            # pure SETTING: no weight, no gate, no minimiser.
            _ring = _chelate_ring_size(lmol, d1, d2)
            _bite = _measured_bite(metal, len(ref), _ring, lsyms[d1], lsyms[d2])
            if _bite is not None:
                u1, u2 = bite_close(u1, u2, _bite)
        T1 = u1 * _r1
        T2 = u2 * _r2
        Q = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=lmol)
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
    bestQ = _place_chelate_block(metal, lsyms, lP, d1, d2, T1, T2, mol=lmol)
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


def _relax_confs_ffree(m, cids):
    """Minimise each conformer with U_total instead of MMFF.

    Routed through ``variational_refine`` rather than a second copy of the minimiser:
    ONE minimiser, one place.  The XYZ round-trip is cheap against L-BFGS.
    ``enable_global_pg=False`` because Tier D is the point group of the WHOLE molecule
    and this is only a cut-out fragment; Tier B/C (Morgan equivalence + graph orbits)
    stay on and are exactly what keeps a symmetric ligand's conformer symmetric.
    A conformer the refiner declines is left untouched -- never-worse per conformer.
    """
    try:
        from delfin.manta._variational_refiner import variational_refine
    except Exception:
        return
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    for cid in cids:
        try:
            conf = m.GetConformer(cid)
            pos = conf.GetPositions()
            xyz = f"{len(syms)}\nconf\n" + "\n".join(
                f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}"
                for s, p in zip(syms, pos))
            new_xyz, rep = variational_refine(xyz, m, class_label="no_metal",
                                              enable_global_pg=False)
            if rep.get("fallback_used", True):
                continue
            n_set = 0
            for ln in new_xyz.splitlines():
                parts = ln.split()
                if len(parts) == 4 and n_set < len(syms):
                    conf.SetAtomPosition(n_set, (float(parts[1]), float(parts[2]),
                                                 float(parts[3])))
                    n_set += 1
        except Exception:
            continue


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
    try:
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                               numThreads=1))
    except Exception:
        # Default (byte-identical): re-raise the original unguarded behaviour.
        if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") != "1":
            raise
        # Kekulize-robust retry (same root as the decompose split): a cleaved
        # aromatic-N⁺ ligand (triazolide / pyridinium / scorpionate cap) carries
        # ARTEFACT positive charges on ring N from the metal-dative SMILES encoding,
        # which leave the ring unkekulizable -> EmbedMultipleConfs raises and the
        # WHOLE complex falls to legacy.  Neutralise those aromatic-N⁺ formal charges
        # (geometry-only: connectivity + donor indices unchanged) so a full sanitize
        # kekulizes and the fragment embeds.  Verified: scorpionate cap embeds.
        cids = []
        try:
            fm = Chem.RWMol(frag_mol)
            for a in fm.GetAtoms():
                if (a.GetIsAromatic() and a.GetSymbol() == "N"
                        and a.GetFormalCharge() > 0):
                    a.SetFormalCharge(0)
            m = fm.GetMol()
            Chem.SanitizeMol(m)
            m = Chem.AddHs(m)
            cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                                   numThreads=1))
        except Exception:
            cids = []
    if not cids:
        try:
            _emb_ok = AllChem.EmbedMolecule(m, randomSeed=SEED) == 0
        except Exception:
            if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") != "1":
                raise                  # byte-identical: original unguarded behaviour
            _emb_ok = False            # kekulize-robust: clean None instead of crash
        if not _emb_ok:
            if key is not None:
                _CONF_CACHE[key] = None
            return None
        cids = [0]
    # FF-FREE CONFORMER RELAX (DELFIN_FFREE_CONF_RELAX, default OFF -> byte-identical).
    # MMFF is the last real force field on the conformer axis, and it has NO metal
    # parameters -- it relaxes the conformers of a COORDINATED ligand with a model that
    # does not know the metal exists, the same defect class as the
    # "UFFTYPER: Unrecognized atom type: Pd+2" this build prints.  The fragment is cut
    # metal-free here, so the functional's existing "no_metal" preset fits exactly:
    # k_topology = 0, k_A = 0, and bond / signature-angle / torsion / clash / symmetry
    # carry the geometry.  ETKDG above is NOT touched: it is distance geometry, not a
    # force field, and it stays the generator.
    if os.environ.get("DELFIN_FFREE_CONF_RELAX", "0") == "1":
        _relax_confs_ffree(m, cids)
    else:
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
    from delfin.manta.refine import _vdw
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
        from delfin.manta import torsion_relax as _TR
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


def _joint_declash_frame(out_syms, P, fixed, block_specs, geom=None):
    """Apply the env-gated JOINT global INTER-LIGAND heavy-heavy declash to one
    assembled frame (``DELFIN_FFFREE_JOINT_DECLASH``), threading the true
    per-ligand connectivity (``block_specs`` = list of ``(offset, lmol,
    donor_local)``).  Runs AFTER #308 torsion-relax and BEFORE the self-gate so a
    declashed class-B build passes ``_build_is_clean``.  No-op when the flag is
    unset; never raises."""
    try:
        from delfin.manta import joint_declash as _JD
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _JD._TR.bonds_from_blocks(0, blocks)
        return np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, geom=geom, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P


def _refine_guarded(out_syms, P, fixed):
    """``refine()`` mit der Zusicherung der Setzung davor und dahinter.

    ⚠️ WARUM DIESE FUNKTION EXISTIERT -- ein Fehler von mir, am 18.08. gemessen und
    hier festgehalten, damit ihn niemand wiederholt.  Ich hatte den Schutz zuerst
    INLINE an EINE Aufrufstelle geschrieben (``assemble_heteroleptic_from_mols``) und
    danach auf den 19 gemessenen OC-6-Faellen geprueft: 19 von 19 byte-identisch.  Das
    sah aus wie "die Zusicherung haelt".  Die Positivkontrolle hat es widerlegt: mit
    Toleranz 0,0001 Angstroem -- wo JEDE Relaxation anschlagen muss -- blieb der Bau
    ebenfalls identisch.  Der Block lief also nie.  Der FF-freie Chelatbauer geht durch
    ``assemble_from_config`` -> ``_finish_config_frame``, eine ANDERE Funktion mit einer
    EIGENEN refine-Aufrufstelle.
    Ich hatte die Zeile auf Erreichbarkeit geprueft und die FUNKTION nicht -- dieselbe
    Bauart wie ``ISOLATED_SEAT``, das ich am selben Tag bei anderen dokumentiert habe.
    ⇒ Der Schutz gehoert an ALLE vier refine-Aufrufstellen, also in EINE Funktion.

    Vorgabe AUS -> byte-identisch: ohne den Schalter ist dies exakt der alte
    ``try: P = refine(...) except: pass``-Block.
    """
    _assert_on = os.environ.get("DELFIN_FFFREE_ASSERT_ENFORCE", "0") == "1"
    _assertion, _P_before = None, None
    if _assert_on:
        try:
            from delfin.manta import _frame_assertions as _FA
            # Die Menge des BAUERS, nicht meine Rekonstruktion davon: refine() bekommt
            # `fixed` als Zusage, also ist genau das der Vertrag, den es halten muss.
            _assertion = _FA.derive((list(out_syms), P), frozen=fixed)
            _P_before = P.copy()
        except Exception:
            _assertion = _P_before = None
    try:
        from delfin.manta.refine import refine as _refine
        P = _refine(out_syms, P, fixed)
    except Exception:
        pass
    if _assertion is not None and _P_before is not None:
        try:
            from delfin.manta import _frame_assertions as _FA
            _v = _FA.violations(_assertion, (list(out_syms), P))
            # ⚠️ EINE SPUR, WEIL EIN BYTE-VERGLEICH HIER NICHT ENTSCHEIDET.
            # Der erste Rauchtest zeigte "identisch" -- und das hat drei mit blossem
            # Auge ununterscheidbare Ursachen: (a) der Block laeuft nicht, (b)
            # derive() liefert None, (c) refine() bewegt nichts, dann ist die
            # Ruecknahme ein No-op.  Genau diese Verwechslung hat am 14.08. den
            # Feuerzensus Befunde erfinden lassen.  Die Spur trennt sie:
            #   derived=1 sagt (b) ab, moved=... sagt (c) ab, broke=1 ist der Treffer.
            # DELFIN_ASSERT_TRACE=<pfad>, sonst still und kostenlos.
            _tp = os.environ.get("DELFIN_ASSERT_TRACE", "")
            if _tp and _tp != "0":
                try:
                    _mv = float(np.max(np.linalg.norm(P - _P_before, axis=1)))
                except Exception:
                    _mv = -1.0
                try:
                    with open(_tp, "a") as _fh:
                        _fh.write("[ASSERT] n=%d derived=1 moved=%.4f broke=%d %s\n"
                                  % (len(out_syms), _mv,
                                     1 if (_v and _v.get("any_broken")) else 0,
                                     "" if not _v else
                                     "md=%d planar=%d frozen=%d trans=%d" % (
                                         _v.get("md_broken", 0), _v.get("planar_broken", 0),
                                         _v.get("frozen_moved", 0),
                                         _v.get("trans_lost_metals", 0))))
                except Exception:
                    pass
            if _v is not None and _v.get("any_broken"):
                P = _P_before              # Ruecknahme: die Behauptung wiegt schwerer
        except Exception:
            pass
    return P


def _sphere_flex_frame(out_syms, P, fixed, block_specs):
    """Apply the env-gated soft coordination-sphere clash relax to one assembled
    frame (``DELFIN_FFFREE_SPHERE_FLEX``).  Donors are soft-restrained (not frozen)
    so the sphere can BREATHE a few hundredths of an Angstrom to open the residual
    inter-ligand heavy-heavy contacts that a frozen-donor refine (and pure M-D-axis
    rotation) cannot.  Threads the true connectivity; runs AFTER joint-declash and
    BEFORE the self-gate.  No-op when the flag is unset; never raises."""
    try:
        from delfin.manta import sphere_flex as _SF
        from delfin.manta import torsion_relax as _TR
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _TR.bonds_from_blocks(0, blocks)
        return np.asarray(_SF.flex_if_enabled(out_syms, P, fixed, bond_pairs=bp),
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
        # ===== DIE SETZUNG UEBERGIBT DER RELAXATION IHRE ZUSICHERUNG ==============
        # Gemessen 2026-08-18 auf 30921 Systemen: netto 988 Systeme fliessen vom
        # Oktaeder ins trigonale Prisma (McNemar X2 = 860,8), der Bauer erzeugt 2,99
        # mal zu viele Prismen, waehrend jede andere Form zwischen 0,86 und 1,29
        # bleibt.  Es sind aber KEINE Prismen -- die CShM-Masse liegt unimodal bei 8
        # bis 12 statt bei 16,7, also ein HALBER Bailar-Twist.  Und es ist keine
        # Auswahl: poly_match ist in 1061 von 1061 Faellen false, obwohl das Auge den
        # besten Frame ueber den GANZEN Manifold liest -- im ganzen Manifold gibt es
        # kein Oktaeder.  Das Signal ist die Verzahnung, monoton von 2,49 % bei null
        # Chelatringen auf 16,12 % bei fuenf; Metall und d-Zahl sind flach.
        #
        # ⇒ Die Setzung stellt das Polyeder richtig, und der Chelatzug dreht es danach
        # heraus.  Genau dafuer ist das Zusicherungsprotokoll gebaut: die Konstruktion
        # sagt, WAS sie behauptet, und die Relaxation darf es nicht brechen.
        #
        # Die vorhandene Gegenmassnahme im Quelltext (smiles_converter.py:38074, die
        # UFF-Winkelziele auf die gegenueberliegenden Donorpaare) haengt ueber :27451
        # an apply_uff und liegt hinter dem FF-freien Return -- sie lief auf diesem
        # Pfad NIE.  Dies hier ist ihr FF-freies Gegenstueck, und es korrigiert nicht,
        # es VERBIETET: bricht die Relaxation die Zusicherung, gilt der Frame VOR der
        # Relaxation.  Never-worse per Konstruktion, kein Zielwert, keine Schwelle,
        # die sich auf einen Pool feintunen liesse.
        #
        # ⚠ WARUM DAS KEIN REPARATEUR IST.  Der Modulzensus vom 18.08. hat gemessen,
        # dass 20 von 21 Reparateuren ohnehin nichts tun und der eine verbleibende
        # (die unbedingte Nachrelaxation) TRAGEND ist -- ohne sie wird jede Achse
        # schlechter (uffoffE).  Die Antwort ist also nicht "keine Optimierung",
        # sondern "keine BLINDE Optimierung".  Dieser Block nimmt der Relaxation
        # nichts weg; er gibt ihr nur, was sie bisher nicht wusste.
        #
        # DELFIN_FFFREE_ASSERT_ENFORCE (Vorgabe 0 -> byte-identisch).
        P = _refine_guarded(out_syms, P, fixed)
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
        P = _joint_declash_frame(out_syms, P, fixed, block_specs, geom=geometry)
        # SOFT coordination-sphere flex (env-gated, default-OFF byte-id): let the
        # donors breathe a hard-bounded amount to open the residual MILD inter-ligand
        # heavy-heavy clashes that rotation/frozen-donor refine cannot (clash forensik
        # 2026-06-29).  Restored toward the ideal vertex+M-D; never-worse on clash.
        P = _sphere_flex_frame(out_syms, P, fixed, block_specs)
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
    _nprod = 1
    for _rl in rank_lists:
        _nprod *= max(1, len(_rl))
        if _nprod > _COMBO_MATERIALISE_MAX:
            break
    if _nprod > _COMBO_MATERIALISE_MAX:      # never materialise a 10^8 product to use 64
        combos = _ranked_combos(rank_lists, 256)
    else:
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
            P = _refine_guarded(out_syms, P, fixed)
            # #308 whole-complex torsion-space clash relax (env-gated, default-OFF
            # byte-id); torsion-only, never-worse, metal+donors (`fixed`) frozen.
            P = _torsion_relax_frame(out_syms, P, fixed, block_specs)
            # JOINT inter-ligand declash (env-gated, default-OFF byte-id): global
            # inter-ligand heavy-heavy minimisation, core frozen.  After #308.
            P = _joint_declash_frame(out_syms, P, fixed, block_specs, geom=geometry)
            # SOFT coordination-sphere flex (env-gated, default-OFF byte-id): donors
            # breathe a bounded amount to open residual mild inter-ligand clashes.
            P = _sphere_flex_frame(out_syms, P, fixed, block_specs)
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
    from delfin.manta import polyhedra as PLY
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

    # LENGTH-GATED HAPTO EXEMPTION (DELFIN_FFFREE_HAPTO_EXEMPT_LENGTH, default OFF -> byte-identical).
    # The hapto path collects every heavy-heavy bond of order >= 1.5 into a LIST, and
    # converter_backend._is_exempt reads the list form as an UNCONDITIONAL pass ("if key in _ex:
    # return True") -- so once a bond is multiple/aromatic, ANY length survives the collapse
    # self-gate, however crushed.  Measured: AJUWUY ships a C-C at 0.495 A through this door.
    # The dict form is already implemented on the reading side (_ex_len + DELFIN_FFFREE_MULTIBOND_TOL
    # 0.15) and gates on `d >= ideal - tol`, which still passes every genuine short multiple bond
    # (C=O 1.13, C#N 1.16, aromatic C~C 1.39) while catching the collapses.  Class: 908 systems,
    # 86.6% of them no_valid -- it tracks the W 69% / Re 53% / Rh 48% / Mo 44% failure rates.
    # Pairs with no known ideal map to 0.0 = unconditional, i.e. exactly today's behaviour.
    _hapto_len_gate = os.environ.get("DELFIN_FFFREE_HAPTO_EXEMPT_LENGTH", "0") == "1"
    if _hapto_len_gate:
        exempt_pairs = {}

    def _collect_exempt(frag_mol, lig_offset):
        """Local heavy-atom double/triple bonds -> global index pairs (the AddHs
        ligand block starts at lig_offset+1 in the assembled coords)."""
        for b in frag_mol.GetBonds():
            bt = b.GetBondTypeAsDouble()
            if bt >= 1.5:                       # aromatic(1.5)/double(2)/triple(3)
                a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                at1 = frag_mol.GetAtomWithIdx(a1)
                at2 = frag_mol.GetAtomWithIdx(a2)
                if at1.GetAtomicNum() > 1 and at2.GetAtomicNum() > 1:
                    g1 = lig_offset + 1 + a1
                    g2 = lig_offset + 1 + a2
                    key = (min(g1, g2), max(g1, g2))
                    if _hapto_len_gate:
                        from delfin.manta.converter_backend import multibond_ideal as _mbi
                        exempt_pairs[key] = _mbi(at1.GetSymbol(), at2.GetSymbol(), bt) or 0.0
                    else:
                        exempt_pairs.append(key)
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
                                             targets[0], targets[1], mol=lg["mol"])
                else:
                    _rigid_seat = (os.environ.get("DELFIN_FFFREE_LIGAND_RIGID", "0") == "1"
                                   and dent >= 3)
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=True,
                                                    rigid=_rigid_seat, lsyms=lsyms)
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
    P = _refine_guarded(out_syms, P, fixed)
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


# ===== GLOBAL DONOR SEATING: the distance nobody in the construction is looking at ===========
#
# WHAT IS MISSING.  Every seating quantity this builder computes is computed for ONE ligand at
# a time: _orient_chelate_to_vertices takes ONE ligand's coordinates, fits ITS donors onto ITS
# assigned vertices and resets ITS M-D radii.  The bite is an intra-ligand distance, r(M-D) is
# a radius, beta is a local plane.  Not one term in the whole seating has a distance between
# TWO DIFFERENT ligands in it.  The conformer PICK sees the neighbour (_clash_count vs
# `placed`), but it is greedy and sequential -- ligand 1 is chosen while ligands 2 and 3 do not
# exist yet -- and once a conformer is picked, nothing can move it relative to its neighbours.
#
# MEASURED, three times independently on 2026-08-02:
#   KEJCUZ dissected: C4...O28 and O15...C39 both sit at 1.85 A = 57 % of the vdW sum, and BOTH
#     pairs are INTER-LIGAND (donor of one chelate against backbone carbon of another).  The
#     per-ligand bite is blind to them by construction.
#   donorplane3 (187 systems): rotating a seated ligand about its donor BEST-FIT line breaks
#     the topology from THREE donors on -- the line passes through none of them.
#   donorplane4 (187 systems): even at TWO donors the backbone swing cost 12 systems
#     (ccdc_arrangement_lost / ccdc_backbone_lost) when the objective was beta.
#
# The lesson of the last two is NOT "never move a seated ligand".  It is "do not move it for a
# reason unrelated to what is broken": beta fired on EVERY system, so it paid the backbone cost
# everywhere and collected a benefit only sometimes.  This fires ONLY where an inter-ligand
# contact is actually below the floor -- a frame with no such contact is returned untouched and
# unchanged, which is most frames.
#
# WHAT IS SOLVED FOR, and why the two things we already get right cannot be traded away.
# The unknowns are the positions of ALL donors of ALL ligands at once.  The two constraints the
# seating already satisfies are enforced by the PARAMETERISATION, not as penalty terms that an
# optimiser could sell:
#     r(M-D)      is invariant under any rotation about the METAL (the metal is the origin in
#                 this frame), so every M-D band stays exactly where md_distance put it;
#     the BITE    -- and with it every intra-ligand distance -- is invariant under ANY rigid
#                 motion of the whole ligand body.  d_DD does not change, so by the ring-closure
#                 law cos(bite) = (r1^2 + r2^2 - d_DD^2)/(2 r1 r2) the bite ANGLE cannot change
#                 either.  Nothing is re-solved; it is carried.
# What is left free is precisely the quantity that was never in any objective: how the ligands
# sit relative to EACH OTHER.  Both invariants are re-CHECKED numerically per move (_gd_move_ok)
# rather than asserted -- one instrument that nobody holds against a second says nothing.
#
# WHY THIS IS NOT joint_declash.  delfin/manta/joint_declash.py already minimises an
# inter-ligand heavy-heavy objective, but it runs with the metal AND ALL DONORS FROZEN, and
# torsion_relax.identify_dofs drops any rotation whose moving half contains a frozen atom other
# than the pivot.  For a CHELATE the M-D1 whole-body spin moves D2, which is frozen -> the DOF
# is dropped, and the D1-D2 line is not a bond so it is not a torsion axis either.  A rigid
# chelate (acac, oxalate, bipy, phen) therefore has NO degree of freedom there at all -- exactly
# the KEJCUZ case.  This pass exists because the donors have to be allowed to move, together.
_GD_CLASH_F = 0.75      # the same clash factor the self-gate and joint_declash use
_GD_H_W = 0.05          # H contacts are a tie-breaker only; the gate rejects on heavy-heavy
_GD_RESID = 0.05        # A, furthest a donor may travel from the vertex it was seated on
# 2026-08-02: 0.25 A was MEASURED and lost.  n=53 affected, valid 25->27, cap_gained=4 but
# cap_LOST=2 (KEGMEP, KIQNUT) and good_regr=1 -> topology_floor=False.  0.25 A at r=2.1 A is
# ~7 deg, which is enough for a donor to leave the polyhedron vertex it was enumerated onto,
# and the topology floor is exactly the term that notices.  0.05 A is ~1.4 deg: the bounded
# stage survives only as a nudge, and what is left is essentially the axis on which the
# donors PROVABLY do not move.  Set DELFIN_FFFREE_GD_RESID to sweep it; 0 disables the
# bounded stage outright (free axis only).
_GD_FREE_STEPS = 24     # 15 deg grid on the zero-cost axis; the objective is smooth
_GD_BOUND_STEPS = 6     # steps each way inside the residual cap
_GD_PASSES = 4          # coordinate-descent sweeps over the ligands
_GD_EPS = 1e-6          # numerical slack on "did not move" -- a rigid rotation about an axis
                        # THROUGH a donor leaves it put only to float precision, so resid=0
                        # must still admit that, or the free axis rejects itself.


def _global_donor_seat_enabled() -> bool:
    """THE one place DELFIN_FFFREE_GLOBAL_DONORS is read (default OFF -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_GLOBAL_DONORS", "0") == "1"


def _gd_resid() -> float:
    """THE one place DELFIN_FFFREE_GD_RESID is read.  Only ever consulted from inside
    _global_donor_seat, i.e. only when DELFIN_FFFREE_GLOBAL_DONORS is on -- with the flag off
    this is dead code and the frame is byte-identical whatever the variable says."""
    try:
        return max(0.0, float(os.environ["DELFIN_FFFREE_GD_RESID"]))
    except Exception:
        return _GD_RESID


def _gd_loss(X, mh, ml, fl):
    """(loss, worst inter-ligand heavy contact).  loss = sum of squared shortfalls below the
    vdW floor over inter-ligand pairs, heavy-heavy dominant.  Same shape as the self-gate's
    own measure, so a move that lowers this is a move toward passing the gate."""
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    over = np.clip(fl - D, 0.0, None)
    L = float((over[mh] ** 2).sum() + _GD_H_W * (over[ml] ** 2).sum())
    return L, (float(D[mh].min()) if mh.any() else float("inf"))


def _gd_move_ok(T, X0, blocks, resid):
    """Re-measure the two invariants instead of trusting the parameterisation.

    r(M-D): the metal is the origin, so |x_d| must be unchanged to numerical precision.
    The BITE: every donor-donor distance INSIDE a ligand must be unchanged likewise.
    Plus the one quantity this pass is allowed to spend -- how far a donor has drifted from
    the vertex the enumeration seated it on -- capped at ``resid`` against the ORIGINAL
    frame (not the previous step), so repeated sweeps cannot accumulate a walk-away.
    ``resid`` is passed in rather than read here: this runs once per candidate angle per
    donor, and an environment lookup in that loop would be the most expensive line in it."""
    cap = resid + _GD_EPS
    for _st, _ln, dn in blocks:
        for a in range(len(dn)):
            da = dn[a]
            if abs(float(np.linalg.norm(T[da])) - float(np.linalg.norm(X0[da]))) > 1e-6:
                return False                                   # M-D band broken
            if float(np.linalg.norm(T[da] - X0[da])) > cap:
                return False                                   # drifted off its vertex
            for b in range(a + 1, len(dn)):
                db = dn[b]
                if abs(float(np.linalg.norm(T[da] - T[db]))
                       - float(np.linalg.norm(X0[da] - X0[db]))) > 1e-6:
                    return False                               # bite broken
    return True


def _global_donor_seat(syms, P, blocks):
    """Move every ligand body -- donors included -- against every OTHER ligand, jointly.

    ``blocks``: one ``(start, n_atoms, [donor global indices])`` per placed ligand; the metal
    is index 0 and never moves.  Returns the improved frame, or None when there is nothing
    below the floor / nothing improved (the caller then keeps the seated frame verbatim).

    Two kinds of rigid motion, tried in this order because the first one is FREE:

      1) rotation about the ligand's OWN donor axis -- the M-D bond for a monodentate, the
         line THROUGH both donors for a bidentate.  The donors lie ON the axis, so they do
         not move at all: r(M-D), the bite AND the vertex alignment are all untouched, and
         only the backbone swings.  Cost to everything already right: exactly zero.
         NOT offered from three donors on: there the "axis" would be a best-fit line through
         none of them, which is the donorplane3 failure verbatim.

      2) rotation of the whole ligand about the METAL, capped so no donor leaves its vertex
         by more than _GD_RESID.  This is the step where the donors genuinely re-seat, and it
         is the only one available to a rigid tri-/tetradentate -- and it is also the step
         that lost the 0.25 A measurement, because a donor that leaves its vertex is exactly
         what the topology floor is watching for.  At resid 0 it is not offered at all and
         only (1) remains.

    Coordinate descent, fixed ligand order, fixed angular grid, no RNG, accept-only-if-better,
    with a never-worse floor on the WORST inter-ligand heavy contact so the sum objective
    cannot buy three mild reliefs by crushing one comfortable pair.
    """
    try:
        X0 = np.array(P, float)
    except Exception:
        return None
    n = len(syms)
    if len(blocks) < 2 or X0.shape != (n, 3) or not np.all(np.isfinite(X0)):
        return None
    try:
        from delfin.manta.refine import _vdw
    except Exception:
        return None
    lig = np.full(n, -1, dtype=int)
    for bi, (st, ln, _dn) in enumerate(blocks):
        if st < 1 or st + ln > n:
            return None                                  # bookkeeping mismatch -> do nothing
        lig[st:st + ln] = bi
    vdw = np.array([_vdw(s) for s in syms], float)
    fl = _GD_CLASH_F * (vdw[:, None] + vdw[None, :])
    isH = np.array([s == "H" for s in syms])
    inter = ((lig[:, None] != lig[None, :]) & (lig[:, None] >= 0) & (lig[None, :] >= 0)
             & np.triu(np.ones((n, n), bool), 1))
    mh = inter & ~isH[:, None] & ~isH[None, :]
    ml = inter & (isH[:, None] | isH[None, :])
    if not mh.any() and not ml.any():
        return None
    L0, h0 = _gd_loss(X0, mh, ml, fl)
    if L0 <= 1e-12:
        return None            # nothing inter-ligand under the floor -> the frame is returned
    hfloor = h0 - 1e-6         # ... unchanged, which is what makes this cheap on clean frames
    resid = _gd_resid()        # read ONCE, outside every loop
    Xc = X0.copy()
    bestL = L0
    for _p in range(_GD_PASSES):
        improved = False
        for st, ln, dn in blocks:
            if ln < 1 or not dn:
                continue
            sl = slice(st, st + ln)
            cands = []                         # (origin, axis, angles) in try-order
            # 1) THE FREE AXIS (see above): donors stay exactly put.  Needs a BODY to swing
            #    -- for a monatomic ligand (Cl-, the atom IS the donor) it moves nothing, so
            #    such a ligand only gets the bounded stage below.
            if len(dn) == 1 and ln >= 2:
                cands.append((np.zeros(3), np.asarray(Xc[dn[0]], float),
                              [2.0 * math.pi * k / _GD_FREE_STEPS
                               for k in range(1, _GD_FREE_STEPS)]))
            elif len(dn) == 2:
                cands.append((np.asarray(Xc[dn[0]], float),
                              np.asarray(Xc[dn[1]], float) - np.asarray(Xc[dn[0]], float),
                              [2.0 * math.pi * k / _GD_FREE_STEPS
                               for k in range(1, _GD_FREE_STEPS)]))
            # 2) THE BOUNDED AXES through the metal: the donors move, together, by at most
            #    the chord _GD_RESID allows.  The chord of a given rotation grows with the
            #    radius, so it is the LONGEST M-D in this ligand that decides the angle for
            #    the whole body -- taking the shortest would let the outer donors overrun the
            #    cap (_gd_move_ok would then throw those candidates away, silently).
            rmax = max(float(np.linalg.norm(Xc[d])) for d in dn)
            if rmax > 1e-6 and resid > _GD_EPS:
                tmax = 2.0 * math.asin(min(1.0, resid / (2.0 * rmax)))
                angs = [s * m for s in
                        [tmax * k / _GD_BOUND_STEPS for k in range(1, _GD_BOUND_STEPS + 1)]
                        for m in (1.0, -1.0)]
                for e in (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]),
                          np.array([0.0, 0.0, 1.0])):
                    cands.append((np.zeros(3), e, angs))
            for org, ax, angs in cands:
                na = float(np.linalg.norm(ax))
                if na < 1e-6:
                    continue
                a_ = np.asarray(ax, float) / na
                locL, locX = bestL, None
                for th in angs:
                    R = _axis_rot(a_, th)
                    T = Xc.copy()
                    T[sl] = (Xc[sl] - org) @ R.T + org
                    if not np.all(np.isfinite(T)):
                        continue
                    Lt, ht = _gd_loss(T, mh, ml, fl)
                    if Lt < locL - 1e-9 and ht >= hfloor and _gd_move_ok(T, X0, blocks, resid):
                        locL, locX = Lt, T
                if locX is not None:
                    Xc, bestL, improved = locX, locL, True
        if not improved:
            break
    if bestL >= L0 - 1e-9 or not np.all(np.isfinite(Xc)):
        return None
    return Xc


def assemble_from_config(metal, geometry, config, ligands, refine=True,
                         n_frames=1, per_lig_confs=6, rmsd_dedup=0.5,
                         planar_bite=None, planar_coplanar=None, prefer_beta=False,
                         lp_orient=False):
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
    lig_blocks = []           # (start, n_atoms, [donor global idxs]) per placed ligand --
                              # which atoms form ONE rigid body, for the global donor seating
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
                    from delfin.manta._multiarm_coplanar import embed_coplanar_multiarm as _ma_emb
                    _ma = _ma_emb(
                        lg["mol"], _dons_d, metal,
                        [float(np.linalg.norm(p)) for p in _dtp],
                        donor_target_pos=_dtp)
                    if _ma is not None:
                        # ADD the coplanar-M conformer to the ensemble (preserve
                        # completeness — do NOT replace the existing conformers;
                        # V16 ADD-not-REPLACE lesson).  _ma is only returned when
                        # genuinely coplanar (accept-if-better gate in the module),
                        # so no regression: ranking/clean-gate pick the best frame.
                        if (ring_confs is not None and ring_confs[0] == _ma[0]):
                            ring_confs = (ring_confs[0],
                                          list(ring_confs[1]) + list(_ma[1]))
                        else:
                            ring_confs = _ma
                except Exception:
                    pass
            # UNIVERSAL LIGAND-INVARIANT cavity seating (DELFIN_FFFREE_RIGID_LIGAND_SEAT,
            # default OFF -> byte-identical).  For ANY rigid polydentate (dent>=3:
            # κ3 pincer / κ4 macrocycle / κ6 cage) replace the metallacycle/coplanar
            # conformer with the FREE-ligand conformer whose donor cavity best admits the
            # metal at the ideal M-D radii (3-D distance solve; NO backbone-flatness
            # requirement -> works where _coplanar_metal_centered_conformer bails).  Paired
            # with the rigid-body orient below (_rigid_seat forced True for this ligand) the
            # ligand internal geometry is preserved EXACTLY and the polyhedron is EMERGENT.
            # Overrides _use_cop / PI_RIGID_PLACE when set (the more general construction).
            # Self-gating: None on solve failure -> keeps the metallacycle embed (never-worse).
            _rls = (_dent >= 3 and _dtp is not None
                    and os.environ.get("DELFIN_FFFREE_RIGID_LIGAND_SEAT", "0") == "1")
            if _rls:
                _mds = [float(np.linalg.norm(p)) for p in _dtp]
                _cav = _rigid_ligand_cavity_conformer(lg["mol"], _dons_d, metal, _mds)
                if _cav is not None:
                    ring_confs = _cav
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
        #
        # COLLAPSE-AWARE SELECTION (DELFIN_FFFREE_COLLAPSE_AWARE_SELECT, default OFF ->
        # byte-identical).  The pick below ranks conformers by CLASH alone, and is blind to
        # the one criterion the self-gate will later reject the whole complex for: a bonded
        # heavy-heavy pair below 0.82 x the covalent sum.  Measured 2026-08-02: 171 of the
        # 215 CHELATE_EMPTY systems die on exactly that, and the mode of the surviving
        # frames' tightest bond sits in the first bin ABOVE the floor -- the gate cuts
        # through the middle of the distribution, it does not trim a tail.
        #
        # _collapsed_heavy_bonds_strict is the gate's OWN predicate, same 0.82, same
        # bonded-pair rule, already unconditional -- it is merely scoped to rigid-planar
        # tridentates further down.  This widens the KNOWLEDGE of it to every chelate
        # without widening the REJECTION: a collapsed conformer is only DEPRIORITISED, never
        # discarded.  If every conformer collapses, the same one wins as today, so the set
        # of buildable systems cannot shrink -- never-worse by construction, not by measurement.
        _csel = os.environ.get("DELFIN_FFFREE_COLLAPSE_AWARE_SELECT", "0") == "1"
        # BETA-AWARE SELECTION (DELFIN_FFFREE_BETA_AWARE_SELECT, default OFF -> byte-id).
        #
        # pyramidal_sp2 -- the metal sitting OUT of its donor's own substituent plane -- is the
        # largest single defect class we carry: 18-19 % of frames against 0.98 % in 509 clean
        # crystals, a 19-fold over-representation.  And it is a PATH difference, not physics:
        # the monodentate path aligns the metal onto the donor's lone-pair axis by construction
        # (_vsepr_reconstruct + _rot_align) and reaches 0.9 deg, BETTER than the crystals; the
        # chelate path never looks at the donor plane at all and reaches 9.3 / 15.9 deg.
        #
        # WHY THIS IS A SELECTION AND NOT A CONSTRAINT.  beta is SECOND ORDER in distance:
        # with d^2 = r^2 + a^2 + 2*r*a*cos(beta)*cos(half) the first derivative vanishes at
        # beta=0, so d changes only as beta^2.  A 0.10 A bound still admits ~37 deg, a 0.05 A
        # bound ~26 deg, and our worst 15.9 deg costs all of 0.018 A.  Holding beta to 5 deg
        # would need a 0.002 A tolerance, which triangle smoothing will not survive.  NO set of
        # distances to atoms IN the donor plane can pin beta to first order -- the distance
        # from a point to coplanar anchors is stationary against perpendicular displacement.
        # (The in-plane azimuth theta IS first order, 0.014 A/deg; that one belongs in bounds.)
        #
        # WHY IT REACHES THE FRAME.  The per-donor radial rescale below scales a donor along
        # M->D, a line that lies IN the plane, so plane(D',X1,X2) is the same plane and still
        # contains M: beta is preserved.  The Kabsch fit is rigid, likewise.  So whatever beta
        # the chosen conformer has is the beta that ships -- the lever belongs here and nowhere
        # after.
        #
        # beta decides ONLY among conformers of EQUAL clash.  It cannot displace the criterion
        # that decides today, so it cannot introduce a clash regression at all.
        # ``prefer_beta`` lets a CALLER ask for the beta-optimal pick without touching the
        # environment.  THE POINT IS THAT IT IS ADDITIVE: the caller builds the frame twice
        # and appends the second as a SIBLING, so the primary is untouched by construction.
        #
        # That shape is not a preference, it is what the record says works.  Every flag that
        # LANDED in this project adds something -- D8_SQ_ADD "a PURELY ADDITIVE sibling ...
        # the PRIMARY frame is untouched", CN6_OH_ADD "adds the OC-6 isomer as a PURELY
        # ADDITIVE sibling", STEREOCENTER_ENUM "additively builds every buildable fold",
        # CN4_BOTH "native-additive".  Everything measured on 2026-08-02 that CHOSE instead
        # of ADDING died: beta-as-replacement cost 3-4 capabilities, collapse-as-replacement
        # died on all three pools, and five reach levers that replaced legacy's frame lost
        # 78/52/52/15/3.  The one lever still alive that night -- the ring pucker -- is the
        # only one that appends and leaves the primary byte-identical.
        _bsel = bool(prefer_beta) or os.environ.get(
            "DELFIN_FFFREE_BETA_AWARE_SELECT", "0") == "1"
        best_Q, best_clash = None, 1e18
        best_coll = True                            # a collapsed pick loses to a clean one
        best_beta = float("inf")                    # ... and among equals, the flatter donor
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
                # κ4 chain: when the chelate's metallacycle embed failed (ring_confs
                # is None, e.g. an unkekulizable aromatic-N⁺ scorpionate cap) the
                # tridentate was NEVER seated -> skipped -> whole complex to legacy.
                # Seat the FREE-ligand conformer too (Kabsch-fit the rigid donor
                # triangle onto the assigned vertices).  Gated default-OFF byte-id.
                _orient_free = (ring_confs is None
                                and os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") == "1")
                if ring_confs is not None or _orient_free:
                    # rigid-body orient (no per-donor rescale) when LIGAND_RIGID OR the
                    # universal cavity seating (RIGID_LIGAND_SEAT) is active for a dent>=3
                    # rigid polydentate -> ligand internal geometry preserved, polyhedron
                    # emergent.  Default OFF on both -> byte-identical (per-donor rescale).
                    _rigid_seat = (dent >= 3 and (
                        os.environ.get("DELFIN_FFFREE_LIGAND_RIGID", "0") == "1"
                        or os.environ.get("DELFIN_FFFREE_RIGID_LIGAND_SEAT", "0") == "1"))
                    Q = _orient_chelate_to_vertices(lP, dons_d, targets, asym=_asym,
                                                    rigid=_rigid_seat, lsyms=lsyms)
                    # THE LEFTOVER ROTATION, SET BY THE LONE PAIRS INSTEAD OF BY THE SVD.
                    # Seating two donors onto two vertices leaves exactly one rotational
                    # freedom -- about the donor-donor axis -- and the Kabsch fit cannot see
                    # it, because both donors lie ON that axis and do not move.  See
                    # _lp_orient_seated_bidentate for why that one angle IS beta, and for the
                    # 995-system measurement behind the law.  Donors, M-D distances, bite and
                    # arm-to-vertex correspondence are all untouched; only the backbone turns.
                    #
                    # AS A SEATING THIS IS MEASURED NEGATIVE, AS A SIBLING IT IS OPEN.
                    # lpseat2 (187 systems, the rotation applied unconditionally): reach 49,
                    # but cap_LOST 2, valid 28->27, mean +0.181 and ten red terms.  The turned
                    # backbones do collide, the self-gate drops those colorings, and the whole
                    # system hands over to legacy.  So the source's own claim at
                    # _lp_aligned_angle -- "the plane wins, without a steric escape" -- is
                    # refuted AT THE LIVE SEATING: a bonding argument does not survive contact
                    # with a co-ligand that is also real.
                    # A guarded version is no answer either: with the collapse guard the same
                    # run measured affected=0, i.e. exactly the "changed almost nothing" the
                    # docstring already recorded once.  Guarded it does nothing, unguarded it
                    # costs capability.
                    # Hence ``lp_orient``: the CALLER asks for the lp-oriented frame and
                    # appends it as a SIBLING, leaving the primary untouched.  Then the plane
                    # can win in a frame of its own without any config ever being dropped,
                    # which is the one shape that has ever landed here.
                    if (Q is not None and dent == 2
                            and (lp_orient
                                 or os.environ.get("DELFIN_FFFREE_LP_SEAT", "0") == "1")):
                        _Ql = _lp_orient_seated_bidentate(Q, lg.get("mol"),
                                                          dons_d[0], dons_d[1], lsyms=lsyms)
                        if _Ql is not None:
                            Q = _Ql
                    # Collapse-rigid-fallback (DELFIN_FFFREE_CHELATE_RIGID_FALLBACK,
                    # default OFF -> byte-id).  The per-donor RADIAL rescale moves each
                    # donor independently onto its exact ideal radius, which on many
                    # chelates distorts the donor-backbone bonds into a collapse near
                    # the metal (the #1 build-coverage gap: 15% selfgate-collapse, 77%
                    # within 3A of the metal).  When the rescaled seating collapsed,
                    # retry the RIGID-body Kabsch fit (no rescale): ligand geometry is
                    # preserved exactly, the M-D distances come out emergent (cavity-
                    # dictated, refined downstream), and the backbone does NOT collapse
                    # -> the config is recovered instead of skipped to legacy.
                    if (Q is not None and not _rigid_seat
                            and os.environ.get("DELFIN_FFFREE_CHELATE_RIGID_FALLBACK", "0") == "1"
                            and _has_collapsed_heavy_bonds(lsyms, Q)):
                        _Qr = _orient_chelate_to_vertices(lP, dons_d, targets,
                                                          asym=_asym, rigid=True, lsyms=lsyms)
                        if (_Qr is not None and np.all(np.isfinite(_Qr))
                                and not _has_collapsed_heavy_bonds(lsyms, _Qr)):
                            Q = _Qr
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
                    Q = _place_chelate_block(metal, lsyms, lP, dons_d[0], dons_d[1],
                                             targets[0], targets[1], mol=lg["mol"])
                if Q is None:                       # tridentate embed failure -> skip this config
                    continue
            if not np.all(np.isfinite(Q)):
                continue
            cl = _clash_count(Q, np.array(placed), lsyms, placed_syms)
            # (collapsed, clash) beats (clash) alone: a clean conformer outranks a colliding
            # one, and among equals the historic clash order is untouched.  With the flag OFF
            # _coll is False for every candidate, so the tuple compare degenerates to the
            # historic `cl < best_clash` exactly -- byte-identical.
            _coll = bool(_csel and _collapsed_heavy_bonds_strict(lsyms, Q))
            _bta = _beta_score(lsyms, Q, dons) if _bsel else 0.0
            if (_coll, cl, _bta) < (best_coll, best_clash, best_beta):
                best_coll, best_clash, best_beta, best_Q = _coll, cl, _bta, Q
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
            elif cl == 0 and not _coll and not _bsel:
                # the early exit must not fire on a conformer that is clash-free but carries
                # a bond the gate will kill the whole complex for -- that is the exact trade
                # this flag exists to stop.  Flag OFF -> _coll is False -> historic break.
                # It must also not fire under BETA_AWARE_SELECT at all: taking the FIRST
                # clash-free conformer is precisely what makes the beta tie-break unable to
                # see the others.  At most per_lig_confs candidates, so the cost is bounded.
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
        # ... and record which atoms belong to this one ligand.  Pure bookkeeping: nothing
        # reads it unless DELFIN_FFFREE_GLOBAL_DONORS is on, so the frame is unaffected.
        lig_blocks.append((pos, len(lsyms), sorted(pos + int(d) for d in dons)))
        pos += len(lsyms)
    donors = sorted(fixed - {0})              # global indices of the constructed donor atoms
    if not ensemble:
        P = np.vstack([np.zeros((1, 3))] + [np.array(placed[1:], float)])
        # GLOBAL DONOR SEATING (default OFF -> byte-identical).  Every ligand up to here was
        # seated ALONE; this is the first and only point where all of them exist at once, so
        # it is the first point where an inter-ligand distance can even be written down.  It
        # runs BEFORE _finish_config_frame on purpose: the relax there pins the donors, so
        # whatever this leaves them at is what the coordination is finished around.
        if _global_donor_seat_enabled():
            _g = _global_donor_seat(out_syms, P, lig_blocks)
            if _g is not None:
                P = _g
        P = _finish_config_frame(out_syms, P, fixed, relax_frags, refine, geom=geometry)
        return out_syms, P, donors

    # ENSEMBLE assembly (Task A.1): enumerate the Cartesian product of per-ligand
    # conformer choices, greedily ordered (sum of candidate ranks -> the best
    # combinations first = frame 0 == the single-path pick), refine each, RMSD-dedup
    # at the complex level, keep up to n_frames.  Capped product keeps it bounded.
    import itertools as _it
    rank_lists = [list(range(len(c))) for c in per_lig_cands]
    _nprod = 1
    for _rl in rank_lists:
        _nprod *= max(1, len(_rl))
        if _nprod > _COMBO_MATERIALISE_MAX:
            break
    if _nprod > _COMBO_MATERIALISE_MAX:      # never materialise a 10^8 product to use 64
        combos = _ranked_combos(rank_lists, 256)
    else:
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
        # same global seating for the ensemble combos; the per-ligand spans are identical
        # across combos (same conformer atom counts), so lig_blocks applies unchanged.
        if _global_donor_seat_enabled():
            _g = _global_donor_seat(out_syms, Pc, lig_blocks)
            if _g is not None:
                Pc = _g
        Pc = _finish_config_frame(out_syms, Pc, fixed, relax_frags, refine, geom=geometry)
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


def _finish_config_frame(out_syms, P, fixed, relax_frags, refine=True, geom=None):
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
    P = _refine_guarded(out_syms, P, fixed)
    # #308 whole-complex torsion-space clash relax (env-gated, default-OFF byte-id):
    # joint multi-axis torsion of all rotatable single bonds, metal + donors (`fixed`)
    # frozen.  Chelate ring + M-D arms are ring bonds -> kept rigid by construction;
    # only the σ-arms/substituents rotate.  Torsion-only, never-worse.  No-op when unset.
    try:
        from delfin.manta import torsion_relax as _TR
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
        from delfin.manta import joint_declash as _JD
        P = np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, geom=geom, bond_pairs=bp),
                       dtype=float)
        # SOFT coordination-sphere radial flex (env-gated, default-OFF byte-id):
        # let crowded monodentate ligands translate radially outward a bounded
        # amount to open residual mild inter-ligand clashes (real-crystal 0.85*vdw
        # target).  This is the assemble_from_config path (the main metal-complex
        # builder); never-worse on clash.  No-op when DELFIN_FFFREE_SPHERE_FLEX unset.
        from delfin.manta import sphere_flex as _SF
        P = np.asarray(_SF.flex_if_enabled(out_syms, P, fixed, bond_pairs=bp),
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
    from delfin.manta import conformer_enum as CE
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


# ---------------------------------------------------------------------------
# Self-test — does the DERIVED bite actually preserve the ligand's own d_DD?
# ---------------------------------------------------------------------------

def _self_test_bite_law() -> None:
    """The claim is not "the angle is nicer" -- it is that d_DD is PRESERVED.

    mdAB measured what happens when it is not: correcting the M-D length while the vertex
    angle stays at the ideal 90 deg forces d_DD = |u1*r1 - u2*r2| to a value the ligand
    does not have, and the ligand absorbs the difference (valid 155 -> 146).  So the test
    that matters compares the donor-donor distance the ligand HAS against the one the
    seating ASKS FOR, with the law off and on.  A lever that leaves that gap unchanged
    would be dead code that looks alive.
    """
    cases = [
        ("NCCN",              [0, 3], "ethylenediamine, sp3 backbone"),
        ("c1ccc(-c2ccccn2)nc1", None, "2,2'-bipyridine, aromatic backbone"),
        ("CC(=O)CC(C)=O",     None,   "acetylacetone, 6-ring"),
    ]
    # WHERE THE MISMATCH ACTUALLY SHOWS.  _place_chelate_block fits the ligand RIGIDLY, so
    # d_DD is preserved either way -- the first version of this test measured that and
    # learned nothing.  A rigid body cannot satisfy two vertex directions AND two M-D
    # lengths when the ligand's own bite disagrees with the vertex spacing, so the
    # compromise lands in the M-D DISTANCES: the donors end up nearer or further than the
    # length that was requested.  That is also the sharper reading of mdAB's failure --
    # setting r is pointless if the rigid fit then moves the donor off it again.
    print("derived bite -- does the seat DELIVER the M-D length it asked for?")
    print(f"{'ligand':>34} | {'d_DD':>6} | {'M-D err off':>11} | {'M-D err on':>10} | "
          f"{'angle off':>9} {'angle on':>8} | OK")
    metal_sym = "Ni"
    for smi, dons, label in cases:
        try:
            lsyms, lP, lmol = _ligand_3d(smi)
        except Exception as exc:
            print(f"{label:>34} | SKIP ({type(exc).__name__})")
            continue
        if dons is None:                      # first two N/O heavy donors in the SMILES
            dons = [a.GetIdx() for a in lmol.GetAtoms()
                    if a.GetSymbol() in ("N", "O")][:2]
        if len(dons) != 2:
            print(f"{label:>34} | SKIP (no donor pair)")
            continue
        d1, d2 = dons
        dd_lig = float(np.linalg.norm(lP[d1] - lP[d2]))
        _want1 = MSB.md_distance(metal_sym, lsyms[d1],
                                 atom=lmol.GetAtomWithIdx(d1), mol=lmol)
        _want2 = MSB.md_distance(metal_sym, lsyms[d2],
                                 atom=lmol.GetAtomWithIdx(d2), mol=lmol)
        errs, angs = [], []
        for flag in ("0", "1"):
            prev = os.environ.get("DELFIN_FFREE_BITE_LAW")
            os.environ["DELFIN_FFREE_BITE_LAW"] = flag
            try:
                # OC-6 vertices 0/1 are TRANS ([1,0,0] / [-1,0,0]); a chelate needs a CIS
                # edge, so 0/2 ([1,0,0] / [0,1,0], 90 deg apart).
                syms, P = assemble_multichelate(metal_sym, "OC-6 octahedron",
                                                [(smi, [d1, d2], [0, 2])])
                # metal sits at row 0; the ligand block follows in ligand order
                v1, v2 = P[1 + d1] - P[0], P[1 + d2] - P[0]
                g1, g2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
                errs.append(max(abs(g1 - _want1), abs(g2 - _want2)))
                angs.append(math.degrees(math.acos(max(-1.0, min(1.0,
                            float(np.dot(v1, v2)) / (g1 * g2))))))
            except Exception as exc:
                errs.append(float("nan"))
                angs.append(float("nan"))
                print(f"    ({label}: {type(exc).__name__}: {exc})")
            finally:
                if prev is None:
                    os.environ.pop("DELFIN_FFREE_BITE_LAW", None)
                else:
                    os.environ["DELFIN_FFREE_BITE_LAW"] = prev
        # NEVER-WORSE is the criterion, not "always changes something".  Three outcomes are
        # all correct: the law repairs a real gap (en), it reproduces what an already-exact
        # rigid fit does (bipy -- a confirmation, not a no-op), or it DECLINES because the
        # triangle does not close (acac, whose freely embedded conformer is the open form,
        # not the chelating one -- d_DD must come from the chelating basin).
        _c = ((_want1 ** 2 + _want2 ** 2 - dd_lig ** 2) / (2.0 * _want1 * _want2)
              if _want1 > 0 and _want2 > 0 else 2.0)
        note = ("declined (cos=%.2f)" % _c if not -1.0 <= _c <= 1.0
                else ("repairs" if errs[0] - errs[1] > 1.0e-6 else "reproduces"))
        ok = errs[1] <= errs[0] + 1.0e-9
        print(f"{label:>34} | {dd_lig:6.3f} | {errs[0]:11.4f} | {errs[1]:10.4f} | "
              f"{angs[0]:9.2f} {angs[1]:8.2f} | {str(ok):>5} {note}")


def _run_self_tests() -> None:
    _self_test_bite_law()
    _self_test_lp_orient()
    _self_test_trilateration()


def _self_test_lp_orient() -> None:
    """How far do the donor lone pairs point PAST the metal?

    That angle IS the flapping.  A conjugated donor's lone pair lies in its pi plane, so a
    correctly seated planar chelate has the metal IN that plane and the angle is zero;
    there is no such thing as a tilted aromatic chelate.  The old orientation was chosen by
    a 36-step sweep maximising the minimum backbone-to-metal distance -- collision
    avoidance where the constraint is bonding -- so this measures what that costs.
    """
    cases = [
        ("c1ccc(-c2ccccn2)nc1",        "2,2'-bipyridine (planar, aromatic)"),
        ("c1cnc2c(c1)ccc1cccnc12",     "1,10-phenanthroline (fused, rigid)"),
        ("NCCN",                       "ethylenediamine (sp3)"),
        ("CC(=O)[O-]",                 "acetate (carboxylate O donors)"),
    ]
    print("\nlone-pair orientation -- degrees the lone pairs miss the metal by")
    print("  'converge' = angle between the two lone pairs measured on the FREE ligand.")
    print("  A chelating conformer aims both at one point, so the two directions must")
    print("  CONVERGE (obtuse, near 180-bite).  If they diverge, no rotation of the whole")
    print("  ligand can aim them at a metal -- the conformer itself is the wrong one, and")
    print("  the seat can only flap it into place.")
    print(f"{'ligand':>36} | {'conv':>5} | {'--':>6} {'lp':>6} {'bite':>6} "
          f"{'both':>6} | {'best':>7} | OK")
    for smi, label in cases:
        try:
            lsyms, lP, lmol = _ligand_3d(smi)
        except Exception as exc:
            print(f"{label:>36} | SKIP ({type(exc).__name__})")
            continue
        dons = [a.GetIdx() for a in lmol.GetAtoms() if a.GetSymbol() in ("N", "O")][:2]
        if len(dons) != 2:
            print(f"{label:>36} | SKIP (no donor pair)")
            continue
        d1, d2 = dons
        _l1, _l2 = _lone_pair_dir(lP, lmol, d1), _lone_pair_dir(lP, lmol, d2)
        conv = (math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(_l1, _l2))))))
                if (_l1 is not None and _l2 is not None) else float("nan"))
        # Four configurations, because the sweep was only ONE candidate root.  A rigid
        # ligand like phenanthroline cannot have a wrong conformer, so if its lone pairs
        # still miss, the cause is upstream: the donors are pinned at an angle and a length
        # the ligand does not aim at.  That is precisely what the derived bite removes, so
        # it has to be in the comparison.
        miss = []
        for flag, bite in (("0", "0"), ("1", "0"), ("0", "1"), ("1", "1")):
            prev = os.environ.get("DELFIN_FFREE_LP_ORIENT")
            prevb = os.environ.get("DELFIN_FFREE_BITE_LAW")
            os.environ["DELFIN_FFREE_LP_ORIENT"] = flag
            os.environ["DELFIN_FFREE_BITE_LAW"] = bite
            try:
                syms, P = assemble_multichelate("Ni", "OC-6 octahedron",
                                                [(smi, [d1, d2], [0, 2])])
                worst = 0.0
                for d in (d1, d2):
                    lp = _lone_pair_dir(P[1:], lmol, d)
                    if lp is None:
                        continue
                    t = P[0] - P[1 + d]
                    nt = float(np.linalg.norm(t))
                    if nt < 1.0e-9:
                        continue
                    c = max(-1.0, min(1.0, float(np.dot(lp, t / nt))))
                    worst = max(worst, math.degrees(math.acos(c)))
                miss.append(worst)
            except Exception as exc:
                miss.append(float("nan"))
                print(f"    ({label}: {type(exc).__name__}: {exc})")
            finally:
                for _k, _v in (("DELFIN_FFREE_LP_ORIENT", prev),
                               ("DELFIN_FFREE_BITE_LAW", prevb)):
                    if _v is None:
                        os.environ.pop(_k, None)
                    else:
                        os.environ[_k] = _v
        best = min(m for m in miss if m == m)
        ok = best <= miss[0] + 1.0e-6             # never-worse is the criterion
        print(f"{label:>36} | {conv:4.0f}d | {miss[0]:6.1f}d {miss[1]:6.1f}d "
              f"{miss[2]:6.1f}d {miss[3]:6.1f}d | {miss[0]-best:+6.1f}d | {ok}")




def _self_test_trilateration() -> None:
    """Does trilateration really hold BOTH sides, where the two old branches each drop one?

    The seating documents the trade-off as unavoidable: per-donor radial keeps the M-D
    lengths and splays the ring, rigid keeps the ring and lets M-D drift.  So the test
    reports exactly those two errors for all three placements.  A lever that merely moves
    the error from one column to the other would be a rename, not a fix.
    """
    from delfin.manta import polyhedra as _PH
    cases = [
        ("tridentate, 2.1 A radii", "OC-6 octahedron", [0, 2, 4], 2.10),
        ("tridentate, 2.4 A radii", "OC-6 octahedron", [0, 2, 4], 2.40),
        ("tetradentate",            "OC-6 octahedron", [0, 2, 4, 1], 2.10),
    ]
    print("\ntrilateration -- both constraints at once?  (errors in Angstrom)")
    print(f"{'case':>24} | {'placement':>14} | {'M-D err':>8} | {'D-D err':>8}")
    for label, geom, verts, rad in cases:
        ref = _PH.ref_vectors(geom)
        # a synthetic rigid donor set: an equilateral-ish arm spacing the ideal vertices
        # cannot satisfy, which is the situation a real polydentate is in
        L = np.array([[0.0, 0.0, 0.0], [2.55, 0.0, 0.0], [1.27, 2.30, 0.0],
                      [1.27, 0.77, 2.10]])[:len(verts)]
        T = [ref[v] / np.linalg.norm(ref[v]) * rad for v in verts]
        D = {(i, j): float(np.linalg.norm(L[i] - L[j]))
             for i in range(len(verts)) for j in range(i + 1, len(verts))}

        def _err(X):
            md = max(abs(float(np.linalg.norm(X[i])) - rad) for i in range(len(verts)))
            dd = max(abs(float(np.linalg.norm(X[i] - X[j])) - D[(i, j)])
                     for (i, j) in D)
            return md, dd

        # (a) ideal vertices, per-donor radial: radii exact, ligand distances ignored
        a_md, a_dd = _err(np.array(T))
        # (b) rigid body onto the ideal vertices: ligand exact, radii drift
        dmu, tmu = L.mean(0), np.array(T).mean(0)
        Rk = _kabsch_rot(L - dmu, np.array(T) - tmu)
        Xb = (L - dmu) @ Rk.T + tmu
        b_md, b_dd = _err(Xb)
        # (c) trilateration, then the same rigid fit onto the moved targets
        tri = _trilaterate_donor_targets(L, list(range(len(verts))), T)
        if tri is None:
            print(f"{label:>24} | trilateration returned None")
            continue
        tgt = np.array(tri)
        Rk2 = _kabsch_rot(L - dmu, tgt - tgt.mean(0))
        Xc = (L - dmu) @ Rk2.T + tgt.mean(0)
        c_md, c_dd = _err(Xc)
        for nm, (md, dd) in (("radial(ideal)", (a_md, a_dd)),
                             ("rigid(ideal)", (b_md, b_dd)),
                             ("TRILATERATED", (c_md, c_dd))):
            print(f"{label:>24} | {nm:>14} | {md:8.3f} | {dd:8.3f}")


if __name__ == "__main__":       # pragma: no cover
    _run_self_tests()
