"""delfin.manta.joint_declash — JOINT global INTER-LIGAND heavy-heavy declash.

The deepest FF-free recall lever for the "class-B" bulk: large / multi-ligand
complexes whose COORDINATION CORE is already IDEAL (donor-angle-RMSD ~0-8°, sane
M-D) but which fall to the distorted legacy-UFF fallback purely because the
ligand BODIES clash — the self-gate ``_build_is_clean`` rejects on genuine
INTER-LIGAND heavy-heavy overlaps (~1.9-2.2 A) even though the first coordination
shell is perfect.

Why a dedicated pass on top of #308 ``torsion_relax``
----------------------------------------------------
#308 minimises a SINGLE global clash SUM that is H-inclusive and counts BOTH
intra- and inter-ligand contacts equally.  On a crowded class-B complex the H-H
and intra-ligand terms DOMINATE that sum, so the coordinate descent spends its
moves relieving cheap H-H / intra clashes and the never-worse-on-MIN floor is set
by an H-H pair — leaving the load-bearing HEAVY-HEAVY inter-ligand overlap (the
ONLY thing the self-gate actually rejects on) under-optimised.  This module makes
the objective the thing the gate measures: a GLOBAL INTER-LIGAND HEAVY-HEAVY
clash sum (H terms kept only as a light secondary tie-breaker), so every move is
spent opening the contacts that block the gate.

Degrees of freedom (identical kinematics to #308 -> identical safety property)
------------------------------------------------------------------------------
For every ligand, jointly:
  * its WHOLE-BODY RIGID ROTATION about the M-donor axis (anchor = metal, pivot =
    donor; the donor lies ON the axis so the M-D distance is exactly preserved),
  * its INTERNAL rotatable-bond torsions (rigid distal sub-tree about each
    non-ring single-bond axis).
Both are pure rigid rotations of a sub-tree about a bond axis: they change ONLY
the dihedral / orientation, never a bond length or a bond angle.  The metal +
ALL donor atoms are FROZEN, so the coordination polyhedron is invariant by
construction (proven in the test-suite: bond-length + bond-angle RMSD before/
after == 0).

Method
------
Internal-coordinate coordinate descent over the joint DOF set (whole-ligand M-D
rotations FIRST — they move the most mass and relieve the inter-ligand overlap
most directly — then internal torsions), each DOF swept on a fixed deterministic
angular grid to its objective-minimising angle.  Accept-if-better (strictly
decreasing), never-worse-on-MIN floor on the inter-ligand heavy minimum, hard
M-D invariant guard (<= md_tol, default 0.05 A; any violation rolls the whole
relaxation back to the input).  Deterministic: fixed DOF order, fixed grid, no
RNG.  Never returns non-finite; on any exception the input frame is returned.

Integration
-----------
Env-gated ``DELFIN_FFFREE_JOINT_DECLASH`` (default ``"0"`` => byte-identical: the
pass is never invoked AND the coupled decompose gate-lift stays at 8).  Runs as a
post-placement pass on each assembled frame, AFTER #308 (if also on) and BEFORE
the self-gate, so a declashed class-B build now PASSES ``_build_is_clean``.

License: open-source / pure geometry only (Bondi-style vdW radii reused from the
FF-free refiner).  No CSD/CCDC data.
"""
from __future__ import annotations

import math
import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin.manta._bond_decollapse as _bd
from delfin.manta import torsion_relax as _TR
from delfin.manta.refine import _vdw
from delfin.manta import polyhedra as _PLY

# --- Metalloid-donor awareness (env-gated, default OFF) ----------------------
# _bond_decollapse._METALS mis-classifies the heavy metalloid sigma-donors
# Sb/Sn/Pb/Ge/Bi as METALS, so this inter-ligand declash treats them as
# coordination CENTRES: the len(metals)==1 guard no-ops and the SbPh3/AsPh3
# ligand BODY is mis-partitioned, so the bulky metalloid ligand is never spun to
# relieve the inter-ligand clash that the short metalloid M-D distance
# (DELFIN_FFFREE_METALLOID_MD_LEN) exposes.  When
# DELFIN_FFFREE_DECLASH_METALLOID_LIGAND=1 a metalloid DONOR counts as a LIGAND
# atom (not a centre), recovering the whole-ligand M-D-axis spin DOF.  The declash
# spins about the M-D axis, so the M-D distance stays EXACTLY invariant -> this can
# never re-detach.  Default OFF -> _is_center == _bd._is_metal (byte-identical).
try:
    from delfin.manta.decompose import _METALLOID_DONORS as _MLD
except Exception:  # pragma: no cover - defensive
    _MLD = frozenset({"Sb", "As", "Bi", "Te", "Se", "Ge", "Sn", "Pb"})


def _is_center(sym: str) -> bool:
    """Coordination-CENTRE test for the declash partition.  Same as ``_bd._is_metal``
    except that, with DELFIN_FFFREE_DECLASH_METALLOID_LIGAND=1, a heavy metalloid
    sigma-donor is a LIGAND atom, not a centre."""
    if sym in _MLD and os.environ.get("DELFIN_FFFREE_DECLASH_METALLOID_LIGAND", "0") == "1":
        return False
    return _bd._is_metal(sym)


# Geometric clash factor — identical to the self-gate / #308 / the spec's f.
_CLASH_F = 0.75

# H-H and X-H contacts kept as a LIGHT secondary tie-breaker: the self-gate
# rejects on HEAVY-HEAVY only, so heavy-heavy must dominate the objective.
#
# ⚠ 25.08.2026 -- THE JUSTIFICATION IN THE LINE ABOVE IS REFUTED.  It says the
# self-gate rejects only on heavy-heavy, so heavy-heavy must dominate.
# Measured over ALL 7631 clash pairs from `archive_aromrad6k_off`:
#     heavy-heavy pairs                            1103  (min 1.733 A, median 2.231 A)
#       of which below GROSS_OVERLAP (0.60 x Sigma_cov)   0   (0.00 %)
#       of which below the refine floor (0.78 x Sigma)    0   (0.00 %)
#       of which below the collapse floor (0.82 x Sigma)  0   (0.00 %)
#     H-involving pairs                            6528  (self-gate skips H)
#   ⇒ reachable by the self-gate: 0 of 7631.
# So the self-gate does NOT reject on heavy-heavy -- it does not reject at
# all.  Its thresholds sit on the BOND scale (C-C 0.92 A), the detector's on
# the VDW scale (C-C 2.38 A).  With that, the premise for 1/20 no longer holds.
#
# WHAT THIS COSTS, measured: 85.5 % of the clash pairs carry at least one H
# (C-H 3751, H-H 1723); of the 3341 frames with a clash, 2507 = 75.0 % have
# EXCLUSIVELY H-involving pairs.  On three quarters of the affected frames this
# declasher sees the whole defect at one twentieth of the weight, and its
# rollback floor `hdmin` (heavy only) notices nothing of it at all.  Pool-wide
# that is around 9.98 % of ALL build frames.
# The crystal says 0.00 % to this on 509 clean structures -- for the H pairs
# too; an H...H below 1.68 A does not exist in real crystals.  The 1/20
# weighting is a builder assumption, not chemistry.
#
# ⛔ DEFAULT UNCHANGED 0.05 -> byte-identical.  Live only with
# `DELFIN_FFFREE_DECLASH_H_FULL=1`.  Reason for the switch: JOINT_DECLASH is ON
# in the champion and has three call sites on the FF-free path -- a silent
# change would be in every running build at once and separable in no A/B.
_H_WEIGHT = 0.05


def _jd_h_voll() -> bool:
    """Do H contacts count in full -- and does the rollback floor count them too?

    Read at CALL time, not at import: a switch whose effect depends on the import
    order has already ended up as a dark switch twice in this project (most
    recently PLANAR_KEEP, two runs with reach 0/24)."""
    return os.environ.get("DELFIN_FFFREE_DECLASH_H_FULL", "0") == "1"


def _jd_h_gewicht() -> float:
    return 1.0 if _jd_h_voll() else _H_WEIGHT

# Default coordinate-descent controls (env-overridable; bounded + deterministic).
_DEF_GRID = 24          # angular grid steps per DOF
_DEF_PASSES = 8         # max coordinate-descent passes over all DOFs
_DEF_MD_TOL = 0.05      # A, hard M-D invariant guard
_DEF_MAX_DOFS = 96      # cap DOFs optimised (bulkiest first; whole-ligand spins kept)


def _enabled() -> bool:
    return os.environ.get("DELFIN_FFFREE_JOINT_DECLASH", "0") == "1"


# ---------------------------------------------------------------------------
# Per-atom ligand membership (so we can score INTER-ligand contacts only)
# ---------------------------------------------------------------------------


def _ligand_of_atom(n: int, syms: Sequence[str],
                    bond_pairs: Optional[Sequence[Tuple[int, int]]],
                    P: np.ndarray) -> np.ndarray:
    """Per-atom ligand id (0..k-1) for the assembled complex; the metal gets -1.

    The ligand graph is the bond graph with the metal-donor coordination bonds
    REMOVED — each connected component on the non-metal atoms is one ligand.
    Using the threaded ``bond_pairs`` (true connectivity) is strongly preferred so
    two interpenetrating ligands at a fortuitous bonding distance are not fused
    into one component (which would hide their inter-ligand clash from the
    objective).  Falls back to geometric perception when ``bond_pairs`` is None.
    """
    adj, _bonds = _TR._adjacency(syms, P, bond_pairs)
    metals = {i for i in range(n) if _is_center(syms[i])}
    lig = np.full(n, -1, dtype=int)
    cur = 0
    for start in range(n):
        if start in metals or lig[start] != -1:
            continue
        # BFS over non-metal atoms only (coordination bonds to the metal are not
        # traversed -> ligands stay separate components).
        stack = [start]
        lig[start] = cur
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b in metals or lig[b] != -1:
                    continue
                lig[b] = cur
                stack.append(b)
        cur += 1
    return lig


# ---------------------------------------------------------------------------
# Inter-ligand clash objective (heavy-heavy dominant, H light tie-breaker)
# ---------------------------------------------------------------------------


def _inter_mask(syms: Sequence[str], lig: np.ndarray,
                excl: Sequence[Set[int]]) -> Tuple[np.ndarray, np.ndarray]:
    """Two upper-triangular (n,n) boolean masks of counted pairs:

      * ``heavy``: i<j, DIFFERENT ligands, neither H, neither metal, not 1-2/1-3
        (the self-gate-relevant inter-ligand HEAVY-HEAVY contacts), and
      * ``light``: i<j, DIFFERENT ligands, at least one H, not 1-2/1-3 (kept only
        as a small secondary tie-breaker so the descent does not introduce gross
        H clashes while opening heavy ones).

    Intra-ligand contacts are NOT counted (this is the *inter-ligand* declash;
    intra-ligand strain is #308's / the refiner's job and is geometry-fixed by the
    rigid sub-tree kinematics anyway).  Pairs sharing the metal (i.e. a ligand
    atom vs the metal) are excluded — the metal is frozen and on the axis.
    """
    n = len(syms)
    heavy = np.zeros((n, n), dtype=bool)
    light = np.zeros((n, n), dtype=bool)
    is_h = np.array([s == "H" for s in syms])
    is_m = np.array([_is_center(s) for s in syms])
    iu, ju = np.triu_indices(n, k=1)
    for k in range(len(iu)):
        i, j = int(iu[k]), int(ju[k])
        if is_m[i] or is_m[j]:
            continue
        if lig[i] == lig[j]:                      # same ligand (or both -1) -> intra
            continue
        if lig[i] < 0 or lig[j] < 0:              # safety: unassigned -> skip
            continue
        if j in excl[i]:                          # 1-2/1-3 (cannot happen inter, but safe)
            continue
        if is_h[i] or is_h[j]:
            light[i, j] = True
        else:
            heavy[i, j] = True
    return heavy, light


def _objective(P: np.ndarray, heavy: np.ndarray, light: np.ndarray,
               rsum: np.ndarray) -> Tuple[float, float]:
    """Return ``(L, heavy_dmin)``.

    ``L = sum_{heavy pairs} over^2 + _H_WEIGHT * sum_{light pairs} over^2`` where
    ``over = max(0, f*(vdw_i+vdw_j) - d)``.  ``heavy_dmin`` is the minimum
    inter-ligand HEAVY-HEAVY distance (the quantity the self-gate gates on).
    """
    diff = P[:, None, :] - P[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    over_h = np.where(heavy, rsum - dist, 0.0)
    over_h = np.where(over_h > 0.0, over_h, 0.0)
    over_l = np.where(light, rsum - dist, 0.0)
    over_l = np.where(over_l > 0.0, over_l, 0.0)
    loss = float((over_h * over_h).sum()) + _jd_h_gewicht() * float((over_l * over_l).sum())
    # ⚠ THE ROLLBACK FLOOR MUST FOLLOW, otherwise the weighting has no effect.
    # `hdmin` is the minimum over the HEAVY pairs; a step that crushes an H
    # contact never drops below this floor and is accepted.
    # A higher H weight in the objective whose guard does not know H would be
    # half a mechanism -- exactly the construction that has already been exposed
    # four times today.  With the switch, the floor counts ALL non-bonded pairs.
    hd = dist[heavy] if not _jd_h_voll() else dist[heavy | light]
    hdmin = float(hd.min()) if hd.size else float("inf")
    return loss, hdmin


# ---------------------------------------------------------------------------
# Core: joint inter-ligand coordinate-descent declash
# ---------------------------------------------------------------------------


def declash(syms: Sequence[str], P, frozen: Iterable[int],
            grid: int = _DEF_GRID, passes: int = _DEF_PASSES,
            md_tol: float = _DEF_MD_TOL, max_dofs: int = _DEF_MAX_DOFS,
            bond_pairs: Optional[Sequence[Tuple[int, int]]] = None,
            geom: Optional[str] = None) -> np.ndarray:
    """Joint global INTER-LIGAND heavy-heavy declash of an assembled complex.

    ``frozen``: indices that must not move (metal + all donor atoms).
    ``bond_pairs`` (optional, strongly preferred): the TRUE connectivity threaded
    from the builder (see :func:`torsion_relax._adjacency`) — both for robust DOF
    identification and for correct ligand-membership partition.

    DOFs: each ligand's whole-body rotation about its M-D axis PLUS its internal
    rotatable-bond torsions (reusing #308's :func:`identify_dofs`, which already
    treats the M-D bond as a rotatable axis with the metal on the anchor side).
    The whole-ligand M-D spins are floated FIRST (they move the most mass and
    relieve the inter-ligand overlap most directly), then internal torsions.

    Returns the relaxed coordinate array.  Accept-if-better, never-worse-on
    heavy-MIN, hard M-D invariant; deterministic; never raises (input on failure).
    """
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return np.array(P, dtype=float)
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return P0
    frozen_set = set(int(x) for x in frozen)

    dofs = _TR.identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                             bond_pairs=bond_pairs)
    if not dofs:
        return P0

    # Reorder: whole-ligand M-D spins (anchor == metal) first, then internal
    # torsions; both already sorted bulkiest-first within #308's identify_dofs,
    # so this is a stable partition preserving determinism.
    metals = {i for i in range(n) if _is_center(syms[i])}
    md_spins = [d for d in dofs if d["anchor"] in metals]
    internal = [d for d in dofs if d["anchor"] not in metals]
    ordered = md_spins + internal

    adj, _bonds = _TR._adjacency(syms, P0, bond_pairs)
    excl = _TR._excl_1_2_3(adj, n)
    lig = _ligand_of_atom(n, syms, bond_pairs, P0)

    # ===== LIGAND SWING (DELFIN_FFFREE_LIGAND_SWING, default OFF -> byte-identical) =====
    #
    # THE OVER-DETERMINATION.  This pass freezes the metal AND ALL DONORS as ATOMS; the
    # polyhedron is thereby "invariant by construction".  That is exactly why a RIGID
    # chelate (acac, bipy, phen) has ZERO degrees of freedom here: its two donors fix
    # the ligand position completely, there are no internal torsions, and identify_dofs
    # rejects every rotation whose moving half contains a frozen donor.
    # Whatever clashes at seating thus clashes forever.
    #
    # MEASURED, 185 869 frames of the 1000-pool: intclash_pair 18.01 % against 0.00 % in
    # the crystal -- the largest single item of the whole output, and in the crystal it
    # does not exist.  Real crystals resolve it by BENDING the polyhedron a few
    # degrees.  We cannot, because we hold the ATOMS fixed instead of the QUANTITY
    # that is chemically truly invariant.
    #
    # THE MISSING MOVE TYPE.  A rigid-body rotation of a WHOLE ligand about the
    # METAL CENTRE preserves EVERY M-D distance exactly (rotation about M leaves all
    # radii unchanged) and changes exclusively the M-D DIRECTIONS.  With that, the
    # chemically hard quantity -- the bond length -- stays exactly held, while the
    # soft quantity -- the vertex angle -- may yield within a band.
    # For a monodentate this is the already existing M-D axis rotation; for a
    # CHELATE it is new: it swings the chelate plane about the axis M -> donor centroid.
    #
    # The swing is tightly capped (LIGAND_SWING_DEG).  The default was 8 degrees and the
    # comment here itself said that this was "not a measured calibration".
    #
    # ⚠ 2026-08-09: IT IS NOW MEASURED, SO IT STANDS HERE.  The curve from 07.08.
    # (8 -> 3 -> 1 degrees, same pool, same eye):
    #     8 degrees   poly_cshm_regressed 2 · poly_lost 1 · poly_type_lost 1
    #     3 degrees   all three ZERO, capability +3/-0, 21:8, mean -0.517
    #     1 degree    too tight, the degree of freedom no longer contributes
    # So the 8 did not merely stand there unsupported, it was REFUTED -- and whoever
    # turns the swing on without DELFIN_FFFREE_LIGAND_SWING_DEG=3 gets the worse number.
    # Exactly the class of error that MONO_REACH_18 was: a threshold beside the distribution.
    # Default therefore 3; byte-identical as long as the switch is off.
    #
    # The two CShM bounds deliberately stay per env and without a default: they are
    # derived from CCDC crystals and must not go into this repo (license).
    if _TR._env_int("DELFIN_FFFREE_LIGAND_SWING", 0, 0, 1):
        _swing_deg = _TR._env_int("DELFIN_FFFREE_LIGAND_SWING_DEG", 3, 1, 30)
        # ===== THE POLYHEDRON BOUND, MEASURED INSTEAD OF GUESSED (2026-08-07) =====
        # The angle cap above is a GUESSED number, and the first run refuted it:
        # at 8 degrees poly_cshm_regressed 2, poly_lost 1, poly_type_lost 1 tripped; at 3
        # degrees all three were gone (capability +3/-0, 21:8, mean -0.517).  Tighter was
        # better on EVERY axis -- the cap was the problem, not the move type.
        #
        # The right quantity is not an angle but THE ONE THE GATE MEASURES: the
        # deviation of the donor set from the ideal polyhedron (CShM), because
        # `poly_cshm_regressed` is the term that tripped.  And the builder can compute
        # it itself -- polyhedra.cshm is pure geometry, no reference data.
        #
        # THE BOUND: the swing may not distort the polyhedron any FURTHER than it already
        # is.  Default 0.0 = strictly never-worse, without any guessed number.  A larger
        # budget can be set per env and is then a MEASUREMENT QUESTION, not a feeling.
        #
        # ⚠ On the calibration, and why it is NOT here: 553 crystals themselves do not
        # sit on the ideal (p50 0.37 · p90 3.65 · p99 8.0 CShM against their own ideal
        # polyhedron).  A budget of this order of magnitude would therefore be
        # chemically covered -- but the number is CCDC-derived and does not belong in
        # this repo.  It lives in the private workspace; here the default stays at 0.
        _swing_cshm = float(os.environ.get("DELFIN_FFFREE_LIGAND_SWING_CSHM", "0") or 0.0)
        _swings = []
        for _lid in sorted({int(x) for x in lig if int(x) >= 0}):
            _atoms = [i for i in range(n) if int(lig[i]) == _lid]
            _don = [i for i in _atoms if i in frozen_set]
            if len(_don) < 2:
                continue          # monodentate: the M-D axis rotation already covers it
            _m = next((i for i in metals), None)
            if _m is None:
                continue
            _c = P0[_don].mean(axis=0) - P0[_m]
            if float(np.linalg.norm(_c)) < 1e-6:
                continue
            _swings.append({"anchor": int(_m), "pivot": -1, "axis_vec": _c,
                            "rotating": _atoms, "max_deg": int(_swing_deg),
                            "cshm_budget": _swing_cshm,
                            "score": 10_000 + len(_atoms)})
        # Swings FIRST: they move the most mass and relieve inter-ligand overlap
        # most directly -- the same reasoning that puts the M-D spins up front.
        ordered = _swings + ordered
    heavy, light = _inter_mask(syms, lig, excl)
    if not heavy.any() and not light.any():
        return P0                                     # nothing inter-ligand to declash

    base_md = _TR._md_pairs(syms, P0)
    vdw = np.array([_vdw(s) for s in syms], dtype=float)
    rsum = _CLASH_F * (vdw[:, None] + vdw[None, :])

    Pcur = P0.copy()
    best_loss, base_hdmin = _objective(Pcur, heavy, light, rsum)
    if best_loss <= 1e-12:
        return Pcur
    # never-worse-on-heavy-MIN floor: a move is rejected if it would push the worst
    # inter-ligand HEAVY contact below the input frame's worst (a SUM objective
    # could otherwise relieve several mild clashes by crushing one comfortable
    # heavy pair -> lower L but a WORSE self-gate-relevant minimum).
    hdmin_floor = base_hdmin - 1e-6

    angles = [2.0 * math.pi * k / float(grid) for k in range(grid)]
    for _ in range(max(1, passes)):
        improved = False
        for dof in ordered:
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            rot = dof["rotating"]
            origin = Pcur[anchor]
            # A ligand swing brings its axis along as a VECTOR (M -> donor centroid);
            # it cannot be expressed as an atom pair, because the centroid is not an
            # atom.  All remaining DOFs stay atom-pair-defined, unchanged.
            _av = dof.get("axis_vec")
            axis = np.asarray(_av, dtype=float) if _av is not None else (Pcur[pivot] - Pcur[anchor])
            if float(np.linalg.norm(axis)) < 1e-9:
                continue
            # Capped swing: sample only the narrow angle window, not the full circle.
            _md = dof.get("max_deg")
            if _md:
                _step = max(1, int(_md) // 4)
                _dofs_angles = [math.radians(d) for d in
                                range(-int(_md), int(_md) + 1, _step) if d != 0]
            else:
                _dofs_angles = angles
            local_best_loss = best_loss
            local_best_P = None
            for ang in _dofs_angles:
                if abs(ang) < 1e-12:
                    continue
                trial = _TR._rotate_subtree(Pcur, origin, axis, ang, rot)
                if not np.all(np.isfinite(trial)):
                    continue
                if not _TR._md_ok(base_md, trial, md_tol):
                    continue
                # POLYHEDRON BOUND for the ligand swing: the same quantity the gate
                # measures.  The M-D DISTANCE is exactly preserved by the rotation about M,
                # but the DIRECTIONS change -- and exactly that tripped poly_cshm_regressed
                # at 8 degrees.  So check here instead of finding out afterwards.
                _cb = dof.get("cshm_budget")
                if _cb is not None and geom:
                    try:
                        _don_idx = [i for i in range(n) if i in frozen_set
                                    and not _is_center(syms[i])]
                        if _don_idx:
                            _m0 = next((i for i in metals), 0)
                            _c_before = _PLY.cshm([Pcur[i] - Pcur[_m0] for i in _don_idx],
                                                  geom)
                            # ===== NO SWING ON AN UNSAFE STARTING SHAPE (2026-08-08) =====
                            # MEASURED, ligswing1k on the 1000-pool: of the 9 damaged
                            # systems, 5 (56 %) had already built the WRONG polyhedron,
                            # in the undamaged comparison group only 7 of 43 (16 %) --
                            # a 3.5-fold enrichment.  (On the smaller 180-pool the signal
                            # was still 43 % against 32 % and thus not robust; only the
                            # larger sample separates.)
                            #
                            # Physically this is inescapable: whoever has already built the
                            # wrong shape optimises the swing WITHIN a shape that is not right.
                            # Every move then leads away from the crystal just as likely as
                            # towards it -- the clash objective knows nothing about that.
                            #
                            # The builder can check this without a crystal: if the donor set
                            # already sits FURTHER from its own ideal polyhedron than a real
                            # crystal at the p90 (CShM 3.65 over 553 crystals), the
                            # starting shape is no trustworthy basis.  The bound comes
                            # per env, because the number is CCDC-derived; default 0 = off.
                            _c_floor = float(os.environ.get(
                                "DELFIN_FFFREE_LIGAND_SWING_MAX_START_CSHM", "0") or 0.0)
                            if _c_floor > 0.0 and _c_before > _c_floor:
                                continue          # unsafe starting shape -> do not swing at all
                            _c_after = _PLY.cshm([trial[i] - trial[_m0] for i in _don_idx],
                                                 geom)
                            if _c_after > _c_before + float(_cb) + 1e-9:
                                continue          # would distort the polyhedron further
                    except Exception:
                        pass                      # not assessable -> the old safeguards apply
                tl, thd = _objective(trial, heavy, light, rsum)
                if thd < hdmin_floor:
                    continue                          # would worsen worst heavy contact
                if tl < local_best_loss - 1e-12:
                    local_best_loss = tl
                    local_best_P = trial
            if local_best_P is not None:
                Pcur = local_best_P
                best_loss = local_best_loss
                improved = True
            if best_loss <= 1e-12:
                break
        if not improved or best_loss <= 1e-12:
            break

    # final hard guard: roll the whole declash back to the input on any M-D
    # violation / non-finite / worse-than-input loss or heavy-min (never-worse).
    if not np.all(np.isfinite(Pcur)) or not _TR._md_ok(base_md, Pcur, md_tol):
        return P0
    fin_loss, fin_hdmin = _objective(Pcur, heavy, light, rsum)
    if fin_loss > best_loss + 1e-9 or fin_hdmin < hdmin_floor:
        return P0
    return Pcur


def declash_if_enabled(syms: Sequence[str], P, frozen: Iterable[int],
                       geom: Optional[str] = None,
                       bond_pairs: Optional[Sequence[Tuple[int, int]]] = None):
    """Wire-in entry: when ``DELFIN_FFFREE_JOINT_DECLASH=1`` run :func:`declash`,
    else return ``P`` unchanged (byte-identical default-OFF).  ``bond_pairs``
    (optional) is the true connectivity from the builder.  Never raises."""
    if not _enabled():
        return P
    try:
        grid = _TR._env_int("DELFIN_FFFREE_JOINT_GRID", _DEF_GRID, 4, 72)
        passes = _TR._env_int("DELFIN_FFFREE_JOINT_PASSES", _DEF_PASSES, 1, 32)
        md_tol = _TR._env_float("DELFIN_FFFREE_JOINT_MD_TOL", _DEF_MD_TOL, 0.0, 2.0)
        max_dofs = _TR._env_int("DELFIN_FFFREE_JOINT_MAX_DOFS", _DEF_MAX_DOFS, 1, 256)
        return declash(syms, P, frozen, grid=grid, passes=passes, geom=geom,
                       md_tol=md_tol, max_dofs=max_dofs, bond_pairs=bond_pairs)
    except Exception:
        return P


# ===== THE ADDITIVE VERSION OF THE M-D ROTATION ==============================
# (DELFIN_FFFREE_MD_SPIN_SIBLINGS, default 0)
#
# WHY THIS FUNCTION EXISTS.  `declash` above CAN already do the rotation about M-D --
# `md_spins` (:259) collects exactly the degrees of freedom whose anchor is the metal,
# and puts them up front.  It is ON in the champion (DELFIN_FFFREE_JOINT_DECLASH).  But
# it is a REPAIRER:
#   * it fires only when an inter-ligand clash is already present, and
#   * it REPLACES the pose instead of appending a second one.
# A monodentate whose azimuth is set merely arbitrarily, but clash-free, thus
# produces NO additional frame -- and the eye's backbone bin stays empty.
#
# MEASURED 20.08. on rows_HIST1KV2_268f120a (969 systems): of 384 systems that miss
# ccdc_backbone, 128 (33.3 %) fail EXCLUSIVELY on metal-containing torsions.  Of the
# 60 among them that carry a coordination number at all, 44 (73 %) are at CN 4/5/6.
# So it is not the MOVE that is missing, it is the OUTPUT FORM.
#
# ⚠ WHY ADDITIVE AND NOT "REPAIR BETTER".  The register is unambiguous on this
# point: what ADDS lands, what CHOOSES dies (03.08.).  And the cost law (four points,
# 19.08.) calls the rotation about an axis THROUGH the metal an ISOMETRY: it leaves
# every distance to M exactly unchanged -- numerically checked, largest |M-D| change
# 0.000000 A.  Price +0.98 pp, the cheapest appending class.
#
# MODEL that is copied here instead of reinvented: `_cn2_spins`
# (assemble_complex.py:2839) -- fixed azimuth steps, RMSD-deduplicated, primary frame
# untouched.  With that, `cap_lost` is impossible by construction: nothing is taken
# away, only placed alongside.
#
# ⚠ BYTE-IDENTITY IS NOT A CLAIM HERE BUT STRUCTURE: this function has ZERO call
# sites in the whole tree.  It cannot change anything as long as nobody calls it.
# The switch below is for the day a call site is added -- that belongs at the place
# where `_cn2_spins` also hands off its siblings, NOT here.
def md_spin_siblings(syms, P, frozen, bond_pairs=None, n_steps=6, max_dofs=4):
    """Sibling poses by rotating whole ligands about their M-D axis.

    Returns a LIST of additional poses (without the input pose).  Empty list if
    the switch is off, no M-D axis exists, or every rotation is degenerate.

    NO clash precondition -- that is the whole difference from `declash`.  The
    azimuth of a monodentate is underdetermined even when nothing clashes; exactly
    these cases are missing from the manifold today.

    The selection of which pose is fit is NOT made by this function but by the
    caller's self-gate -- as with every sibling.  Whoever filters here already
    builds a chooser again.
    """
    if os.environ.get("DELFIN_FFFREE_MD_SPIN_SIBLINGS", "0") != "1":
        return []
    try:
        P0 = np.array(P, dtype=float)
    except Exception:
        return []
    n = len(syms)
    if n < 3 or P0.shape != (n, 3) or not np.all(np.isfinite(P0)):
        return []
    try:
        frozen_set = set(int(x) for x in frozen)
        dofs = _TR.identify_dofs(syms, P0, frozen_set, max_dofs=max_dofs,
                                 bond_pairs=bond_pairs)
        if not dofs:
            return []
        metals = {i for i in range(n) if _is_center(syms[i])}
        # ONLY the M-D axes.  Internal torsions are a different axis with a different
        # price (rigid rotation with a new conformation, +6.57 pp) and do not belong
        # in the same output -- otherwise a verdict cannot be attributed afterwards.
        spins = [d for d in dofs if d.get("anchor") in metals]
        if not spins:
            return []
        step = 360.0 / max(2, int(n_steps))
        out, seen = [], [P0]
        for d in spins:
            anchor, pivot = int(d["anchor"]), int(d["pivot"])
            rot = d.get("rot") or d.get("rotating")
            if rot is None:
                continue
            rot = [int(x) for x in rot]
            origin = P0[anchor]
            axis = P0[pivot] - P0[anchor]
            if float(np.linalg.norm(axis)) < 1e-9:
                continue
            for k in range(1, int(n_steps)):
                ang = step * k
                try:
                    trial = _TR._rotate_subtree(P0, origin, axis, ang, rot)
                except Exception:
                    continue
                if trial is None or not np.all(np.isfinite(trial)):
                    continue
                # RMSD dedup against ALL previous ones, not only the input pose:
                # for a C2-symmetric ligand the half rotation coincides with the
                # starting position, and a duplicate with a label is not a frame.
                if any(float(np.sqrt(np.mean(np.sum((trial - q) ** 2, axis=1)))) < 0.25
                       for q in seen):
                    continue
                seen.append(trial)
                out.append(trial)
        return out
    except Exception:
        return []
