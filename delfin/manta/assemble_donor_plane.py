"""Donor-plane relaxation of the FF-free assembly: collapsed-bond checks, donor follow weights, arm order, riding hydrogens, beta scoring and the trilateration rescue switch.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import numpy as np
import os
from rdkit import Chem

import delfin.manta._bond_decollapse as _bd


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
        factor = _bd.COLLAPSE_FLOOR   # ONE source for the collapse floor, see _bond_decollapse
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


def _ffree_flag(name: str) -> bool:
    """ONE SWITCH, TWO SPELLINGS -- and only one of them can ever land.

    ===== THE LETTER ON WHICH 50.8 PERCENT HANG (26.08.2026) ======================

    `cli_manta.py:293` sets the champion like this:

        os.environ["DELFIN_FFFREE_" + f] = "1"        # THREE F

    But there is a whole class of switches that is read as `DELFIN_FFREE_` -- TWO F.
    No switch of this class can ever get into the champion, no matter how good its
    verdict turns out.  Affected, with MEASURED reach:

        TRILATERATE     50.8 %      MD_MEASURED   19.9 %      H_FOLLOW  6.4 %
        CONF_RELAX       4.0 %      TRILAT_RESCUE  3.7 %

    🔑 And `TRILAT_RESCUE` is not just any of them: it is the LAST RUNG of the
       rescue ladder below the chelate self-gate, which rejects 2822 of 3597 chelate
       configs (78.5 %).  Line 1 of `_trilat_rescue` is a return on exactly this
       switch -- and it is read at ONE place and set at NONE.

    WHY NOT SIMPLY RENAME.  Measured: 31 of 1286 axis files set the OLD name.
    A rename would make these 31 archives unreproducible -- they would have been
    measured with a name the code no longer reads.  That is exactly the class of
    silent failure this campaign is built against.

    ⇒ So read BOTH.  The tree already does this elsewhere: `_mirror_enum`
      reads `DELFIN_FFFREE_MIRROR_ENUM` or `DELFIN_MIRROR_ENUM` in one line.

    ⛔ THIS IS NOT A LANDING, BUT ITS PRECONDITION.  Before, the mechanism was
      unlandable; now it is landABLE.  Reach is not benefit -- it still needs
      its own verdict.
    ⛔ Default stays OFF under BOTH names -> byte-identical.
    """
    return (os.environ.get("DELFIN_FFFREE_" + name, "0") == "1"
            or os.environ.get("DELFIN_FFREE_" + name, "0") == "1")


def trilat_rescue_enabled() -> bool:
    """THE one place TRILAT_RESCUE is read (default OFF -> byte-identical).

    Reads both spellings, see `_ffree_flag`.  Without that, this rung of the
    rescue ladder can never get into the champion, because the setter writes three F."""
    return _ffree_flag("TRILAT_RESCUE")


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
    return _TRILAT_RESCUE or _ffree_flag("TRILATERATE")


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
