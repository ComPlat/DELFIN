"""Ligand embedding for the FF-free assembly: rigid alignment, Kabsch, ring-bound tightening, metallacycle embedding, planar polydentate placement and rigid cavity conformers.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.manta.assemble_donor_plane import (
    _trilaterate_donor_targets,
)


SEED = 42


# Bin width of the polyhedron-fidelity selection (18.08.2026).  See the long
# justification at the selection site: unbinned, the floating-point number would decide
# every comparison and the landed beta term would never get its turn again.
_POLY_FIDELITY_BIN = 0.05


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
        # together ("no OC-6 polyhedron, donors bunched" — eye-flagged ADOHOT/ADOROD/ZUSBEU/
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
