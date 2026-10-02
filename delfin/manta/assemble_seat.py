"""Seating passes of the FF-free constructor: global donor seat, OC-6 twist seat and ligand-DOF reseat.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import numpy as np
import os

from delfin.manta.assemble_ligand_confs import (
    _clash_count,
)
from delfin.manta.assemble_ligand_embed import (
    _kabsch_rot,
    _rot_align,
)
from delfin.manta.assemble_orient import (
    _axis_rot,
)


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


# ===== THE HALF BAILAR TWIST, TURNED BACK FF-FREE ===========================
# Measured 18.08.2026 on 30921 systems: net +988 systems flow from the octahedron
# into the trigonal prism (McNemar X2 = 860.8 on 1 df).  TPR-6 is built 2.99 times
# as often as it really occurs, while EVERY other shape lies between 0.86 and 1.29
# -- the largest single defect of the polyhedron axis, 3.7 times the mass of the
# second-largest pair.
#
# THREE MEASUREMENTS SAY WHAT IT IS NOT:
#   * It is NOT a ligand-field effect.  Metal flat, d-count flat; the signal is
#     solely the INTERLOCKING by chelate rings, monotone from 2.49 % at zero rings
#     to 16.12 % at five (6.5-fold).  So geometry, not chemistry.
#   * They are NOT real prisms.  CShM(OC-6) lies at a median of 11.03 instead of
#     16.7, as an ideal TPR would have -- a HALF Bailar twist, stopped
#     halfway.
#   * It is NOT a selection question.  poly_match is false in 1061 of 1061 cases,
#     although the eye reads poly_build as the minimum over ALL realistic frames.
#     In the whole manifold there is no octahedron -- so none is built.
#
# ⛔ WHY THE EXISTING REPAIR IS NOT ENOUGH.  It exists twice, both in
# smiles_converter.py: DELFIN_FFFREE_CN6_OH_ADD (:27452) begins with `apply_uff and`,
# and DELFIN_FFFREE_CN6_OH_ANGLES (:38153) sits in
# _build_coordination_constraints_from_xyz, i.e. in the UFF constraint machinery.  Both
# give UFF octahedral angle targets (90/180 degrees) so that UFF relaxes the twist out.
# On the FF-free path no UFF runs ⇒ reach zero.  They cannot be
# wired in; they must arise anew on the construction side, and that is this block.
#
# ⚠ WHAT IS DELIBERATELY CARRIED OVER HERE, AND WHY.  The three safeguards of the UFF
# version are not trimmings, they are the reason it is isomer-safe:
#   1) The three trans pairs come from the FRAME ITSELF (greedy: each donor with its
#      most opposite one), NOT from an enumerator permutation.  The
#      PERM variant was measured on 2026-07-14 and COLLAPSED isomers
#      (VOYWUD lost all-trans + trans-OH) -- it forces the wrong trans set on some
#      arrangements.  What the frame already has stays: fac stays fac,
#      cis stays cis.  It is a twist correction, not an arrangement change.
#   2) Only if ALL three pairs are clearly trans (min > 120 degrees).  A valid
#      TPR/OC frame sits at 140-180; an ambiguous one does not -> skipped, so that
#      nothing collapses.
#   3) Only at CN 6 and _PREFERRED_CN6_GEOMETRY.get(metal, 'OH') == 'OH'.
#
# ⚠ AND WHAT IS DIFFERENT -- that is the reason why this can land.  By the cost law
# measured on 18.08. at three points, ordering/selection costs about 0,
# isometry +0.98 pp, a rigid rotation with NEW conformation +6.57 pp and a
# re-embedding +11.9 pp.  What ADDS, without inventing new geometry, lands.
# That is why here the WHOLE ligand arm is rotated RIGIDLY about the metal and not the
# single donor atom shifted: a rotation about the metal leaves r(M-D) exact and
# the bite exact, it invents no conformation.  Shifting single donors
# would tear bonds -- exactly the re-embedding that was measured to be the most expensive.
#
# A chelate with a 78-degree bite CANNOT give a perfect octahedron, and it is not
# meant to: the Kabsch fit puts its donors as close to the ideal directions as
# its bite allows, and the bite wins.  The goal is not CShM 0, the goal
# is "no more half twist".
_OC6_TRANS_MIN = 120.0   # degrees; below it the pair is not unambiguously trans -> abort


_OC6_AXIS_MIN = 0.20     # smallest singular value of the three axes: below it they are


                         # almost coplanar and the orthonormalisation would be guessed
_OC6_RIGID_TOL = 1e-6    # Angstrom; re-measurement of r(M-D) and bite AFTER the rotation


# ⚠ ONE EXIT, ONE LINE -- and `call` right at the front.  The function has ten ways
# to come back with None, and every single one means something different: "not my case",
# "ambiguous, hands off", "rotated, but it achieved nothing".  Without the denominator
# they would all be the same zero in the report -- the mistake that was made here on
# 10.08. and once more on 14.08.  Only recorded if the corrector is called at all,
# i.e. only behind oc6_twist=True.  Read by _self_test_oc6_twist and by the
# sibling self-test in converter_backend.
_OC6_SEAT_CENSUS = dict.fromkeys(
    ("call", "not_oc6", "metal_pref", "shape", "book", "not_cn6", "zero_md",
     "ambiguous", "coplanar_axes", "rot_bad", "rigid_broken", "cshm_flat", "ok"), 0)


# (CShM before, CShM after, smallest trans angle) per call -- the raw numbers from
# which to read WHETHER the seating is twisted at all.  A counter alone could not
# say that: "not improved" means either "already correct" or "too
# bad to rescue", and those are opposite findings.  Only behind oc6_twist.
_OC6_SEAT_CSHM = []


_OC6_CSHM_KEEP = 4096   # cap; the counters above stay complete, only the


                        # raw-value list stops growing at some point


def _oc6_twist_seat_enabled() -> bool:
    """THE one read site of DELFIN_FFFREE_OC6_TWIST_SEAT (default 0 -> byte-identical).

    It does NOT decide about the primary frame.  ``assemble_from_config`` executes the
    correction exclusively on the keyword ``oc6_twist=True``, which
    is False by default -- so the built primary frame is byte-identical, no
    matter what is in the environment.  This switch only says whether the caller
    additionally builds a SIBLING FRAME."""
    return os.environ.get("DELFIN_FFFREE_OC6_TWIST_SEAT", "0") == "1"


def _oc6_trans_pairs(u, donors):
    """The three trans pairs from the frame itself: each donor with its most
    opposite one, greedy in fixed index order (deterministic, no RNG).

    ``u``: dict donor index -> unit vector from the metal.  Returns
    ``(pairs, smallest_trans_angle_deg)`` or ``(None, 0.0)``.  Word for word the
    pairing of the UFF version in smiles_converter.py:38160 -- not out of convenience,
    but because EXACTLY this pairing is the isomer-safe one (see head note)."""
    rem = list(donors)
    pairs = []
    min_trans = 180.0
    while len(rem) >= 2:
        a = rem[0]
        b = min(rem[1:], key=lambda x: float(np.dot(u[a], u[x])))
        c = max(-1.0, min(1.0, float(np.dot(u[a], u[b]))))
        min_trans = min(min_trans, math.degrees(math.acos(c)))
        pairs.append((a, b))
        rem.remove(a)
        rem.remove(b)
    if len(pairs) != 3:
        return None, 0.0
    return pairs, min_trans


def _oc6_ideal_axes(u, pairs):
    """The three measured trans axes, pulled onto the NEAREST orthonormal triad
    (polar decomposition, ``A = U S Vt`` -> ``U Vt``).

    WHY POLAR DECOMPOSITION AND NOT GRAM-SCHMIDT: Gram-Schmidt is order-dependent
    -- the first axis would stay untouched, the third would carry the whole error.  The
    polar decomposition minimises the sum of squares over all three simultaneously and
    is thereby independent of which pair was found first.  That matters,
    because the greedy pairing above has an index order, but the chemistry does not.

    The handedness is NOT corrected.  What is sought are three mutually perpendicular
    unit vectors; {±e1, ±e2, ±e3} is the same octahedron whether the triad is right-
    or left-handed.  A det correction would be no protection here, but an
    additional, unnecessary rotation.  Returns ``(E, smallest_singular_value)``."""
    A = []
    for a, b in pairs:
        ax = u[a] - u[b]
        na = float(np.linalg.norm(ax))
        if na < 1e-9:
            return None, 0.0
        A.append(ax / na)
    A = np.asarray(A, float)
    try:
        U, S, Vt = np.linalg.svd(A)
    except Exception:
        return None, 0.0
    if not np.all(np.isfinite(U)) or not np.all(np.isfinite(Vt)):
        return None, 0.0
    return U @ Vt, float(S[-1])


def _oc6_twist_seat(syms, P, blocks, metal, geometry):
    """Rotate the half twist out -- rigidly, ligand by ligand, about the metal.

    ``blocks``: per seated ligand one ``(start, n_atoms, [global donor indices])``,
    the same bookkeeping ``_global_donor_seat`` uses; atom 0 is the metal.
    Returns: the corrected frame, or ``None`` if one of the safeguards
    triggers OR the correction does not measurably reduce the twist.  The caller
    then keeps the seated frame verbatim.

    Procedure:
      1) CN 6, OC-6 requested, metal prefers OH -- otherwise nothing.
      2) trans pairs from the frame, all three clearly trans (> 120 degrees).
      3) orthonormalise the three axes -> the octahedron NEAREST to the frame.
         Not the lab-axes octahedron: the nearest is the one requiring the least
         movement, and movement is exactly what is paid for under the cost
         law.
      4) per ligand ONE rigid rotation about the metal that puts its donors onto
         their target directions in the Kabsch sense.  A monodentate ligand gets
         the minimal rotation (Rodrigues), from two donors on the Kabsch fit -- which
         is NOT degenerate for two points, because the metal at the origin is held
         along and the covariance thus has rank 2, whose null direction is unique.
      5) RE-MEASURE instead of trust: r(M-D) and every intra-ligand donor-donor
         distance must be unchanged to 1e-6, and CShM(OC-6) must have STRICTLY
         fallen.  Both are measurements on the result, not thresholds tuned to a
         pool -- a rotation that does not reduce the twist
         is rejected, instead of being dressed up."""
    def _no(reason):
        _OC6_SEAT_CENSUS[reason] += 1
        return None

    _OC6_SEAT_CENSUS["call"] += 1
    if not str(geometry).startswith("OC-6"):
        return _no("not_oc6")              # a REQUESTED TPR-6 stays a TPR-6
    # Metal preference.  Deferred import as in _finish_config_frame; if the
    # module fails, the default 'OH' applies -- the same the table itself gives.
    try:
        from delfin.smiles_converter import _PREFERRED_CN6_GEOMETRY as _PCN6
        if _PCN6.get(str(metal), 'OH') != 'OH':
            return _no("metal_pref")
    except Exception:
        pass
    try:
        X0 = np.asarray(P, float)
    except Exception:
        return _no("shape")
    n = len(syms)
    if X0.shape != (n, 3) or not np.all(np.isfinite(X0)) or not blocks:
        return _no("shape")
    donors = []
    for st, ln, dn in blocks:
        if st < 1 or st + ln > n:
            return _no("book")             # bookkeeping does not fit -> do nothing
        donors += [int(x) for x in dn]
    if len(donors) != 6 or len(set(donors)) != 6:
        return _no("not_cn6")              # CN 6, and every donor exactly once
    M = X0[0].copy()
    u = {}
    r = {}
    for d in sorted(donors):
        v = X0[d] - M
        nv = float(np.linalg.norm(v))
        if nv < 1e-6:
            return _no("zero_md")
        u[d] = v / nv
        r[d] = nv
    pairs, min_trans = _oc6_trans_pairs(u, sorted(donors))
    if pairs is None or min_trans <= _OC6_TRANS_MIN:
        return _no("ambiguous")            # ambiguous -> skip, nothing collapses
    E, smin = _oc6_ideal_axes(u, pairs)
    if E is None or smin < _OC6_AXIS_MIN:
        return _no("coplanar_axes")        # almost coplanar axes -> the triad would be guessed
    tgt = {}
    for i, (a, b) in enumerate(pairs):
        e = np.asarray(E[i], float)
        ne = float(np.linalg.norm(e))
        if ne < 1e-9:
            return _no("coplanar_axes")
        e = e / ne
        if float(np.dot(u[a], e)) < 0.0:
            e = -e                         # the axis points towards a, not away from a
        tgt[a] = e * r[a]
        tgt[b] = -e * r[b]
    Xc = X0.copy()
    for st, ln, dn in blocks:
        dn = [int(x) for x in dn]
        if ln < 1 or not dn:
            continue
        obs = np.asarray([X0[d] - M for d in dn], float)
        tar = np.asarray([tgt[d] for d in dn], float)
        if len(dn) == 1:
            R = _rot_align(obs[0], tar[0])
        else:
            R = _kabsch_rot(obs, tar)
        if R is None or not np.all(np.isfinite(R)):
            return _no("rot_bad")
        Xc[st:st + ln] = (X0[st:st + ln] - M) @ R.T + M
    if not np.all(np.isfinite(Xc)):
        return _no("rot_bad")
    # 5a) RE-MEASURE the two invariants.  A rotation about the metal holds them
    #     mathematically; they are measured anyway, because a degenerate Kabsch matrix
    #     could silently smuggle in a mirroring exactly here.
    for st, ln, dn in blocks:
        dn = [int(x) for x in dn]
        for i, da in enumerate(dn):
            if abs(float(np.linalg.norm(Xc[da] - M)) - r[da]) > _OC6_RIGID_TOL:
                return _no("rigid_broken")             # r(M-D) broken
            for db in dn[i + 1:]:
                if abs(float(np.linalg.norm(Xc[da] - Xc[db]))
                       - float(np.linalg.norm(X0[da] - X0[db]))) > _OC6_RIGID_TOL:
                    return _no("rigid_broken")         # bite broken
    # 5b) and the one number that matters.  If it does not fall, the block has nothing
    #     to offer and returns the frame unchanged (the caller keeps it).
    try:
        from delfin.manta import polyhedra as _PH
        before = _PH.cshm([X0[d] - M for d in sorted(donors)], "OC-6 octahedron")
        after = _PH.cshm([Xc[d] - M for d in sorted(donors)], "OC-6 octahedron")
        if len(_OC6_SEAT_CSHM) < _OC6_CSHM_KEEP:      # capped: a 30-hour run
            _OC6_SEAT_CSHM.append(                    # must not drag a list along
                (float(before), float(after), float(min_trans)))
    except Exception:
        return _no("cshm_flat")
    if not (after < before - 1e-9):
        return _no("cshm_flat")
    _OC6_SEAT_CENSUS["ok"] += 1
    return Xc


# ===== THE LIGAND'S LAST DEGREE OF FREEDOM ===================================
# THE FINDING THAT FORCES THIS BLOCK (18.08., 492 clean remaining systems).
# `org_bond` -- the largest defect mass of the organic geometry -- does NOT hang on
# the bond class (within a molecule the localisation is massive, but identical between
# hits and failures, max |delta| 0.14).  Two quantities genuinely enrich:
# `worst_n` (number of simultaneously bent organic bonds) 3.83x, and the
# COORDINATION NUMBER 1.62x (CN>=5 51.3 % against CN<=4 31.7 %).  And CN acts THROUGH
# `worst_n`: mean `worst_n` rises 3.46 -> 7.07 (CN 2..6) at practically
# constant molecule size, with a CORRECTLY built polyhedron 3.6 against 7.1.
# ⇒ At equal ligand size a CN-6 centre bends twice as many
# organic bonds as a CN-4 centre.  The polyhedron is enforced, the ligand pays.
#
# THE BUILDER'S DOF BALANCE (re-measured, see _self_test_ligand_dof):
#   * A ligand block is seated RIGIDLY -- 6 rigid-body DOF.
#   * Monodentate: `_rot_align(lp, -Vunit)` fixes 5 of them (3 translation via
#     `+ Vunit*md`, 2 direction via the Rodrigues rotation).  The SIXTH -- the
#     azimuth about the M-D axis -- is chemically completely undetermined and is
#     nonetheless nailed down, namely to the random value that the MINIMAL rotation of
#     `_rot_align` happens to deliver.  It is NEVER sampled on the configuration path.
#   * Bidentate: the Kabsch fit in `_orient_chelate_to_vertices` puts both donors
#     onto their vertices; exactly ONE rotation remains -- about the donor-donor axis.
#     The fit cannot see it (both donors lie ON the axis).  Today it is
#     used only by `_lp_orient_seated_bidentate` to optimise a BOND angle
#     (default OFF, measured negative as a seating) -- for the PACKING nobody
#     has ever used it.
#   * Tridentate and higher: three non-collinear donors fix the rigid body
#     COMPLETELY.  ZERO free rotation.  The register is right: 0 DOF.
#
# WHAT THIS AXIS PRESERVES -- and that is the reason it is built.  Both
# cases rotate about a line that CONTAINS EVERY donor of this block.  The donors
# are thereby POINTWISE fixed: M-D distance, M-D DIRECTION, bite, vertex angles and CShM
# do not change by one bit.  It is an isometry of the ligand block on an
# already built conformation -- no new embedding, no new conformation.
# By the cost law measured at four points (ordering/selection ~0, isometry
# +0.98 pp, rigid rotation with NEW conformation +6.57 pp, re-embedding +11.9 pp)
# that is the cheapest class that can change anything about the packing at all.
#
# ⚠ DISTINCTION FROM DELFIN_FFFREE_LIGAND_SWING (joint_declash.py:305).  The swing is
# a DIFFERENT movement at a DIFFERENT place: it rotates the whole ligand about the
# axis M -> donor centroid, AFTER assembly, in the declash pass, and in doing so THE
# DONORS MOVE -- that is why it needs a 3-degree cap and a CShM budget, and
# that is why it can damage the polyhedron at all (at 8 degrees it did).  This
# axis here does not move the donors, therefore needs no angle cap and may sample the
# full circle; and it acts at the SELECTION, where the candidate is still being
# decided, not afterwards on the finished frame.  Extending the swing would mean
# rebuilding it into a movement it is not -- and in a file that does not belong to this
# territory.  Hence a second axis, not an extension.
_LIGAND_DOF_AXIS_TOL = 1.0e-6      # collinearity of the donors / axis proximity (Angstrom)


def _ligand_dof_seat_enabled() -> bool:
    """THE one place where DELFIN_FFFREE_LIGAND_DOF_SEAT is read
    (default 0 -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_LIGAND_DOF_SEAT", "0") == "1"


def _ligand_dof_seat_steps() -> int:
    """Number of sampled angles on the full circle (default 6, like the already
    existing CN2 axis rotation DELFIN_FFFREE_CN2_SPINS).  Only read when the
    switch above is on."""
    try:
        n = int(os.environ.get("DELFIN_FFFREE_LIGAND_DOF_SEAT_N", "6"))
    except Exception:
        n = 6
    return max(2, min(n, 36))


def _free_rigid_axis(Q, donor_locals, metal_pos):
    """The ONLY line about which this seated ligand block may still be rotated
    rigidly without changing the coordination sphere by even one bit -- or ``None``.

    Condition: the axis must contain EVERY donor of the block, then all donors are
    pointwise fixed and thus M-D distance, M-D direction, bite and CShM exactly preserved.
      * 1 donor  -> the line M--D (among all lines through the donor the only one that
                    additionally preserves EVERY M-X distance of the ligand, because the
                    metal then lies on the axis itself -- the ligand cannot swing into
                    the metal).
      * >=2 donors -> the line through the donors, but ONLY if they are collinear.
      * otherwise -> ``None``.  Three non-collinear donors fix the rigid body
                    completely; rotating anything here would mean touching the polyhedron.
    Returns ``(origin, unit_axis)``; ``origin`` is a donor point, so that this donor
    stays BIT-EXACT under the rotation.
    """
    try:
        d = [np.asarray(Q[int(i)], float) for i in donor_locals]
    except Exception:
        return None
    if not d or any(not np.all(np.isfinite(x)) for x in d):
        return None
    if len(d) == 1:
        a = d[0] - np.asarray(metal_pos, float)
        n = float(np.linalg.norm(a))
        if n < _LIGAND_DOF_AXIS_TOL:
            return None
        return (d[0], a / n)
    a = d[-1] - d[0]
    n = float(np.linalg.norm(a))
    if n < _LIGAND_DOF_AXIS_TOL:
        return None
    a = a / n
    for p in d[1:-1]:                       # Kollinearitaet ALLER Donoren
        w = p - d[0]
        if float(np.linalg.norm(w - float(np.dot(w, a)) * a)) > _LIGAND_DOF_AXIS_TOL:
            return None
    return (d[0], a)


def _free_dof_reseat(Q, lsyms, donor_locals, existing, existing_syms, base_clash,
                     metal_pos=None):
    """Samples the one free degree of freedom of this ligand block and returns
    ``(Q_rotated, clash)`` -- or ``None`` if nothing is STRICTLY better.

    NEVER-WORSE BY CONSTRUCTION, at three places:
      1) Switch off -> immediately ``None``, the caller sees nothing.  Byte-identical.
      2) The same quantity is measured that the selection reads anyway as its first
         key (`_clash_count` against the already seated atoms) -- no
         new criterion that could win against the old one.
      3) Only a STRICTLY smaller clash is taken; on a tie the
         historic pose stays.  And only a block that ALREADY collides
         (`base_clash > 0`) is touched at all -- a collision-free seating stays
         untouched, there is nothing to gain in the packing there.
    Deterministic: fixed angle list, smallest index wins on a tie.
    """
    if not _ligand_dof_seat_enabled():
        return None
    if base_clash <= 0 or len(Q) < 2:
        return None
    ax = _free_rigid_axis(Q, donor_locals,
                          np.zeros(3) if metal_pos is None else metal_pos)
    if ax is None:
        return None
    o, u = ax
    W = np.asarray(Q, float) - o
    perp = W - np.outer(W @ u, u)
    if float(np.max(np.linalg.norm(perp, axis=1))) < _LIGAND_DOF_AXIS_TOL:
        return None                     # everything lies ON the axis -> rotation = identity
    n = _ligand_dof_seat_steps()
    best = None
    for k in range(1, n):
        try:
            R = _axis_rot(u, 2.0 * np.pi * float(k) / float(n))
        except Exception:
            continue
        Qk = W @ R.T + o
        if not np.all(np.isfinite(Qk)):
            continue
        ck = _clash_count(Qk, existing, lsyms, existing_syms)
        if ck < base_clash and (best is None or ck < best[1]):
            best = (Qk, ck)
    if best is None:
        return None
    # THE ASSERTION IS RE-MEASURED, NOT CLAIMED.  Mathematically a rotation about an
    # axis through all donors holds every donor point fixed -- it is measured
    # anyway, because exactly here a degenerate axis could silently do something
    # else, and the price would be the coordination sphere.  If the assertion breaks,
    # the rotation is REJECTED (the caller keeps its historic pose).
    _mp = np.zeros(3) if metal_pos is None else np.asarray(metal_pos, float)
    for _d in donor_locals:
        _d = int(_d)
        if float(np.linalg.norm(best[0][_d] - Q[_d])) > 1.0e-9:
            return None                                   # donor has moved
        if abs(float(np.linalg.norm(best[0][_d] - _mp))
               - float(np.linalg.norm(np.asarray(Q[_d], float) - _mp))) > 1.0e-9:
            return None                                   # r(M-D) broken
    _tp = os.environ.get("DELFIN_LIGAND_DOF_TRACE", "")
    if _tp and _tp != "0":
        # ⚠️ TRACE, BECAUSE "byte-identical" HAS SEVERAL CAUSES HERE: the block does
        # not run, there is no free axis, nothing collided, or the rotation found
        # no better pose.  Without this line one could not say which applies.
        try:
            with open(_tp, "a") as _fh:
                _fh.write("[LIGDOF] ndon=%d nat=%d clash %d -> %d\n"
                          % (len(list(donor_locals)), len(Q), int(base_clash),
                             int(best[1])))
        except Exception:
            pass
    return best
