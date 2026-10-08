"""Orientation of chelates onto polyhedron vertices, donor bend angles, sp-chain straightening, VSEPR reconstruction and diatomic donor orientation of the FF-free assembly.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import itertools
import math
import numpy as np
import os

import delfin.manta._bond_decollapse as _bd
from delfin.manta.assemble_donor_plane import (
    _donor_follow_weights,
    _donor_plane_relax,
    _ffree_flag,
    _hydrogens_riding_on,
    _riders_that_may_move,
    _trilat_targets_on,
    _trilaterate_donor_targets,
)
from delfin.manta.assemble_ligand_embed import (
    _kabsch_rot,
    _rot_align,
)


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
                   and _ffree_flag("H_FOLLOW"))
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
                and _ffree_flag("H_FOLLOW"))
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
