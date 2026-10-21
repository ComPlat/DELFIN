"""Fold fingerprints of the FF-free constructor: ring folds from blocks, the amplitude axis, metallacycle arms and complex RMSD for the dedup.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import numpy as np
import os

from delfin.manta.assemble_donor_plane import (
    _canonical_arm_order,
)
from delfin.manta.assemble_ligand_embed import (
    _kabsch_rot,
)


# ===== THE MISSING EYE OF THE DEDUP: THE FOLD FINGERPRINT =========================
#
# WHAT WAS MEASURED (26.08.2026, harness/faltung_dedup_schwelle.py, raw data under
# results/FALTUNG_DEDUP_2026_08_26/).  Every dedup in the build computes RMSD over
# ALL heavy atoms.  But a ring fold moves only the RING ATOMS -- their
# median share of the heavy atoms is 0.1316.  Measured counterfactually, on real
# archive frames folded with the build's own generator (`_ring_pucker._set_pucker`
# via `_pucker_candidates`, substituents ride along rigidly):
#
#   archive gkfam6kb_on (5681 files, 500 drawn, 199 evaluable, 459 rings)
#       4539 of 5754 fold candidates below 0.50 A total RMSD        =  78.88 %
#       3311 of 5754                    below 0.30 A               =  57.54 %
#       2672 of 5754                    below 0.25 A               =  46.44 %
#   archive hplacegate6k_on (4313 files, 136 evaluable, 327 rings)
#       3117 of 4206 fold candidates below 0.50 A                  =  74.11 %
#   chair against COUNTER-chair, the most classic flip of all:
#         94 of  294 even rings         below 0.50 A               =  31.97 %
#
# ⛔ DO NOT LOWER THE THRESHOLD.  That would be a knob, gameable, and it would hit
#    EVERY other axis as well -- rotamers, seating variants, combinations.  Instead,
#    the PREDICATE itself is tightened:
#
#        rmsd < thr        ->        rmsd < thr  AND  same fold
#
# ⚠ WHY THIS CANNOT EXPLODE HERE, although a CP dedup on 26.08. in
#   `_ring_pucker` failed on EXACTLY THAT.  There it had been measured: without CP
#   n=5 converges at 3/3/3 states, with CP it explodes to 9/13/14 -- all
#   thirteen at theta=90 and Q=0.300, distinguished ONLY by phi.  That is the
#   pseudorotation of an unsubstituted ring; phi hangs on the ATOM NUMBERING,
#   not on the chemistry.  TFD folds in the topological symmetry, CP cannot
#   do that.  CP ALONE is no substitute for TFD.
#     The difference here is the CONJUNCTION.  There CP stood alone and DECIDED;
#   here it stands BEHIND the RMSD and can only RESCUE pairs that already lie below
#   the RMSD threshold -- i.e. are almost congruent atom by atom.  Two
#   pseudorotamers 60 degrees apart are precisely NOT that: half the ring set
#   switches side, the RMSD separates them anyway.  The explosion mode lies by
#   construction outside the reach of this term, and the number of additionally
#   kept frames is capped from above by `n_frames` / `max_builds`.
#
# RESOLUTION NOT GUESSED, BUT TAKEN FROM THE GENERATOR.  `_pucker_candidates`
# samples the CP equator with K = max(8, 2n) phases; the finest spacing the
# generator itself draws is therefore 360/K degrees.  The tolerance below is HALF the
# sampling spacing, 180/K: below it two states lie within the granularity of the
# generator, above it they are two different candidates for the generator itself.
#     n=4  22.5   n=5  18.0   n=6  15.0   n=7  12.86   n=8  11.25 degrees
# No free parameter, no turned knob.
#
# ⚠ COST, MEASURED instead of assumed (harness/faltung_fp_kosten.py):
#       _cp_theta_phi   n=5 75.2 us · n=6 85.9 us · n=7 87.9 us   per ring and frame
#       _cp_abstand      6.9 us                                   per ring and PAIR
#       _tfd          1435.8 us                                   per PAIR (35 atoms)
#   Three consequences from that:
#     1. TFD is ruled out at this site.  It is by construction a PAIR
#        quantity and cannot be cached per frame; in `_dedup_builds` with
#        up to 180 candidates against up to 60 kept, that is 10 800 pairs
#        = 15.5 s per complex.  That is not an instrument verdict, that is a price.
#     2. The fingerprint is computed and memorised PER FRAME, never per pair.  Naively
#        per pair it would be 21 600 CP evaluations instead of 180 at the same site --
#        exactly the quadratic trap.
#     3. It is moreover computed LAZILY: only when an RMSD proximity
#        occurs at all.  Where nothing is deduplicated, it costs zero.
#   Upper bound thereby: 16 rings x 86 us = 1.4 ms per frame, x 180 frames = 0.25 s per
#   complex in the worst case; measured on average were 459/199 = 2.31 foldable
#   rings per system, so the cap practically never bites.
#
# ⚠ DEFERRED IMPORT, and the reason is NOT circularity.  Measured (both
#   load orders, each in a fresh interpreter): `_ring_pucker` pulls in on
#   loading only `delfin`, `delfin.manta`, `delfin.manta._ring_pucker` -- NO
#   `assemble_complex`.  A module-level import would therefore be permitted.  It stays
#   in the function nonetheless, because otherwise it would be a LOAD-TIME coupling: a
#   momentarily broken `_ring_pucker` would then take down the whole builder, even with
#   the fingerprint switched off.  With default OFF the import never runs.
#
# Switch: DELFIN_FFFREE_DEDUP_FOLD_FP (default 0 -> predicate byte-identical).

_FOLD_FP_QMIN = 0.075        # A -- half the smallest amplitude the build produces


                             # (`_cp_pucker_amps` gives +/-0.15 A for eta faces).
                             # Below it the ring is flat and theta/phi are noise.
_FOLD_FP_RINGMAX = 12        # largest ring size considered


_FOLD_FP_MAXRING = 16        # cost cap: rings per frame (deterministically sorted)


def _fold_fp_enabled():
    return os.environ.get("DELFIN_FFFREE_DEDUP_FOLD_FP", "0") == "1"


def _fold_rings_from_blocks(blocks, syms):
    """Global ring index lists of the potentially FOLDABLE rings.

    ``blocks`` is a sequence of ``(global_offset, mol)``, where
    ``global_offset + local_index`` is the index in the frame -- the convention the
    builder itself writes (``_collect_exempt``: "the AddHs ligand block starts
    at lig_offset+1 in the assembled coords").

    The filter is deliberately COARSE: not fully aromatic, size 4..12.  The
    actual gate is the amplitude floor at run time -- a rigid, flat
    ring has Q ~ 0 in BOTH frames and thus counts as equally folded anyway.
    The filter only saves compute time, it decides nothing.

    ⚠ RETURNS ``None`` AS SOON AS THE OFFSET MODEL DOES NOT WORK OUT.  No verdict
      is better than a wrong one: with ``None`` the predicate falls back to the old,
      pure RMSD.  Every ring atom is checked against its element symbol
      in the frame.  ⚠ What that does NOT catch: ``_ligand_confs_from_mol`` memoises
      its conformer pools by canonical SMILES; two constitutionally identical
      ligands with different internal atom order get the same pool,
      and then local j can be a DIFFERENT atom of the same name.  That is a
      pre-existing trait of the builder (lines 3705/3722/3744 mix the same two
      index spaces); here the consequence would at most be one additionally kept
      near-duplicate, never a lost frame."""
    rings = []
    try:
        for off, mol in blocks:
            if mol is None:
                continue
            off = int(off)
            for r in mol.GetRingInfo().AtomRings():
                n = len(r)
                if n < 4 or n > _FOLD_FP_RINGMAX:
                    continue
                if all(mol.GetAtomWithIdx(int(j)).GetIsAromatic() for j in r):
                    continue                     # flat, rigid face: no axis
                g = [off + int(j) for j in r]
                if min(g) < 0 or max(g) >= len(syms):
                    return None
                for j, gj in zip(r, g):
                    if syms[gj] != mol.GetAtomWithIdx(int(j)).GetSymbol():
                        return None
                rings.append(tuple(g))
    except Exception:
        return None
    if not rings:
        return None
    rings.sort()                                 # deterministic, independent of SSSR
    return rings[:_FOLD_FP_MAXRING]


def _fold_fp(P, rings):
    """The Cremer-Pople state per ring, computed with the instrument of the
    FOLD GENERATOR itself (``_ring_pucker._cp_theta_phi``) instead of re-implemented.
    Per entry ``(n, Q, theta, phi)``.  ``None`` = no verdict possible."""
    try:
        from delfin.manta._ring_pucker import _cp_theta_phi     # deferred, see above
    except Exception:
        return None
    out = []
    try:
        for r in rings:
            Q, th, ph = _cp_theta_phi(P, list(r))
            out.append((len(r), float(Q), float(th), float(ph)))
    except Exception:
        return None
    return out


# ===== THE THIRD AXIS OF THE SPHERE: THE AMPLITUDE (26.08.2026) ==================
#
# The fingerprint above read (Q, theta, phi) and compared ONLY (theta, phi) of it.
# Q stood in the tuple and was held exclusively against the floor -- the
# DIRECTION of the fold decided everything, its DEPTH nothing.  Two gaps follow
# from that, and both are measured, not suspected
# (`harness/faltung_fp_rettung_metallacyclus.py`, `archive_gkfam6kb_on`, seed 11,
#  500 systems -- the same draw as all numbers of this block):
#
#   (1) FLAT AGAINST FOLDED.  If ONE ring lies below the floor and the other
#       above, the great-circle distance judges by the theta/phi of the flat one --
#       and that is the phase of a vanishing displacement, i.e. noise.
#       Measured 1191 such pairs, of which 146 (12.26 %) judged as SAME.
#       That is the `5M:planar` state, which the eye counts as its own basin:
#       a flat and a folded ring passed as the same fold.
#   (2) FLAT AGAINST DEEP.  Two rings with the same fold direction, but
#       Q = 0.15 against Q = 0.60, have great-circle distance ZERO -- a hinted
#       and a pronounced boat counted as the same state.
#
# ⚖ THE TOLERANCE IS DERIVED, NOT CHOSEN -- by THE law by which the
#   two already existing numbers of this instrument are built:
#
#       TOLERANCE = HALF THE GENERATOR'S SAMPLING SPACING ON THIS AXIS.
#
#   * angle: `_pucker_candidates` draws K = max(8, 2n) phases around the equator,
#     spacing 360/K, half spacing 180/max(8,2n)  -- exactly the line below.
#   * floor:  the smallest amplitude the build sets at all is 0.15 A
#     (`_cp_pucker_amps` for eta faces); the two neighbouring states there
#     are 0 and 0.15, half spacing 0.075  -- exactly `_FOLD_FP_QMIN`.
#   * amplitude: `_pucker_candidates` emits `q_scale` from {0.0, 1.0} (the zero under
#     DELFIN_FFFREE_PUCKER_PLANAR) and multiplies by `_ring_pucker._amp(n)`.
#     So the ladder is {0, _amp(n)}, its spacing _amp(n), the half spacing
#
#         Q_TOL(n) = _amp(n) / 2      n=4 0.175 · n=5 0.200 · n=6 0.315
#                                     n=7 0.360 · n=8 0.400 · n=12 0.49
#
#   ⇒ The floor is thereby no longer a second parameter, but the same expression
#     for the eta ladder: 0.15/2.  Two constants, ONE law, NO free knob.
#
# 🔬 THE SECOND INSTRUMENT, and it says almost the same.  A derived number without
#   a second measurement would be a claim -- so the cross-check:  `_ring_pucker`
#   has long deduplicated its OWN fold candidates on BOTH axes,
#
#       `_cp_abstand(_cp,_k) < _cp_tol  AND  |Q_cp - Q_k| < _cp_qtol`   (:1144)
#
#   with `_cp_tol` = 15 degrees and `_cp_qtol` = 0.15 A -- the latter NOT from the
#   grid, but crystallographically justified (coordinate esd 0.002..0.01 A
#   as lower, the disorder limit 0.3..0.5 A as upper bound).
#   ⇒ Two independent routes, ONE predicate form.  The fingerprint here was
#     until today the WEAKER half of it: it demanded only the angle.
#   ⇒ On the angle the two even meet exactly: 180/max(8,2*6) = 15.0 degrees.
#     On the amplitude the number derived here is by a factor 1.17 (n=4) to
#     2.67 (n=8) COARSER than the crystallographic one.  The direction of this
#     deviation is the safe one (see below); its size is the rescue this
#     draft deliberately leaves on the table.
#
# ⚠ WHY THE COARSE LADDER AND NOT THE FINE ONE.  `_pucker_space_grid` divides
#   the same span into `n_amp` steps (default 2) and would give _amp(n)/4.  The
#   LARGER tolerance is taken, and for the same reason the coarse angle tolerance
#   above stays in place: `_fold_same` returns True on "within the
#   tolerance", and True means DUPLICATE -- the old behaviour.
#   A too-coarse tolerance can only leave rescue on the table, never wrongly keep a
#   frame.  The measured rescue numbers are thereby LOWER BOUNDS.
#
# ⚠ (1) NEEDS NO TOLERANCE AT ALL, and therefore gets none.  The floor itself
#   already says what a ring below it IS: flat, without a fold axis.  A
#   flat and a folded ring are, by exactly this definition, not
#   the same fold -- that is no new number, but the missing half
#   of a case distinction that until now had only its symmetric branch
#   ("BOTH flat -> same").
#
# ⚠ THE CHANGE IS MONOTONE, and that is its most important property.  Every
#   pair that was previously called DIFFERENT stays different; only pairs that
#   were called SAME can flip.  So no rescue is lost, and the
#   byte proof at default OFF stays the same (`_fold_same` is unreachable
#   as long as `fold_rings is None` -- all three call sites bail out before).
#   ⚠ Not only argued, but RECOUNTED: in the re-measurement the rescue
#     rises in EVERY bucket of both axes or stays equal, in none does it fall
#     -- organic n=4..7 plus chair-against-counter-chair, metallacycle n=4..12
#     plus chair-against-counter-chair, plus the three overall tables.
#
# 🎯 THAT IT IS NOT NOISE IS SHOWN BY THE FOUR-RING -- and analytically, not
#   statistically.  An N-ring has N-3 fold degrees of freedom; for the FOUR-RING that
#   is EXACTLY ONE, and that one is the amplitude.  In `_cp_theta_phi` the term
#   `q2s = -sqrt(2/n) * sum z_j sin(pi j)` drops identically to zero for n=4 (sin(pi j) = 0
#   for integer j), so phi is constantly 0 or 180, and q2/q3 stand in a
#   fixed ratio, so theta too takes only two values.
#   ⇒ For a four-ring, (theta, phi) NEVER carried information about the fold, only
#     its SIGN.  The old fingerprint was blind there by construction.
#   ⇒ And exactly there the jump is largest -- measured, not expected:
#         organic    n=4   32 -> 78 of 94 eaten            34.04 % -> 82.98 %
#         MC         n=4   67 -> 205 of 356                18.82 % -> 57.58 %
#     while the six-ring, where the direction really carries two degrees of freedom,
#     barely gains (organic 3083 -> 3111 of 3136, 98.31 % -> 99.20 %).
#   A noise term would have scattered evenly across all ring sizes.  This one
#   hits the class that was known beforehand to HAVE to be blind.
#
# ⚠ IN WHICH DIRECTION THIS DRAFT ERRS, if it errs: it SPLITS TOO MUCH.
#   Two rings just either side of the floor (0.074 against 0.076) are practically
#   both flat and are still called different -- a spurious variant, one frame
#   too many.  That is the more expensive, but the right error direction: a
#   surplus conformer costs compute time, an eaten state is
#   irretrievable, and the north star is completeness.  The error in
#   the opposite direction persists anyway (coarse tolerance, see above).


def _fold_same(fa, fb):
    """True = the SAME fold (or no verdict possible -> old behaviour).

    Comparison uses the GREAT-CIRCLE DISTANCE on the CP sphere
    (``_ring_pucker._cp_abstand``), not |dtheta|+|dphi|.  At the pole (theta 0
    or 180 -- chair and counter-chair) phi is meaningless; the naive metric
    holds two identical chairs with phi=136 and phi=339 to be 200 degrees
    apart.  The great-circle distance solves that geometrically, without a special rule.

    The AMPLITUDE enters as well -- derivation in the block above.  Order of the
    gates: first the floor (is the sphere responsible at all?), then the
    amplitude (how DEEP), then the direction (where to)."""
    if not fa or not fb or len(fa) != len(fb):
        return True
    try:
        from delfin.manta._ring_pucker import _cp_abstand, _amp   # deferred, see above
    except Exception:
        return True
    for a, b in zip(fa, fb):
        n = int(a[0])
        if int(b[0]) != n:
            return True                          # ring lists do not match: no verdict
        if max(a[1], b[1]) < _FOLD_FP_QMIN:
            continue                             # both flat -> no fold axis
        if min(a[1], b[1]) < _FOLD_FP_QMIN:
            return False                         # ONE flat, ONE folded: the floor
                                                 # itself calls that two states
        if abs(a[1] - b[1]) > 0.5 * _amp(n):
            return False                         # depth: half the generator's
                                                 # amplitude step ({0, _amp(n)})
        if _cp_abstand(a[1:], b[1:]) > 180.0 / max(8, 2 * n):
            return False
    return True


# ===== THE FINGERPRINT'S GAP: THE RING AROUND THE METAL (26.08.2026) ==============
#
# The fingerprint above takes its rings from ``mol.GetRingInfo()`` of the LIGAND MOL.
# There is no metal there, so there is no CHELATE RING there either.  The fold
# of a metallacycle -- the "step" of a salen ring, the flip of an
# ethylenediamine five-ring -- until now ran unhindered into the same complex RMSD and
# was eaten there as a duplicate.
#
# ⚠ THIS IS NOT AN EDGE CASE, AND THE NUMBER IS MEASURED, NOT ESTIMATED
#   (`harness/faltung_fp_rettung_metallacyclus.py`, `archive_gkfam6kb_on`,
#    seed 11, the same 500 drawn systems as for the organic measure):
#       reach        314 of 500 systems carry a metallacycle (62.8 %)
#       rings        816 found, of which 706 after the aromatics rule below
#                    (n=5: 336 · n=6: 307 · n=4: 33 · n>=7: 30)
#       candidates   8520  -- MORE than the 5754 organic ones of the predecessor measure
#       eaten        5550 of 8520 lie below 0.50 A complex heavy-atom RMSD
#       RESCUED      3630 of 5550 (65.41 %) = 42.61 % of all 8520
#   For comparison the organic measure of the same day: 4257 of 4539 (93.79 %).
#   So the metallacycle is rescued LESS OFTEN, but there are more candidates.
#   ⚠ Without the aromatics rule it would be 816 rings / 9620 candidates / 4445 of
#     6621 rescued.  The 110 dropped rings are the bipyridine class
#     (fully aromatic apart from the metal); they are NOT in the
#     headline number here, because the generator does not fold them in the first
#     place -- see below.
#   ⚠ THESE FOUR NUMBERS ARE THE STATE BEFORE THE AMPLITUDE.  The Q term in the block
#     above `_fold_same` raised them on the same day; the same tools,
#     the same archive, the same seed 11, the same 500 systems -- only the predicate
#     is sharper (`results/FALTUNG_FP_Q_2026_08_26/`):
#         metallacycle   3630 -> 3913 of 5550   65.41 % -> 70.50 %
#         the same without the aromatics rule
#                        4445 -> 4780 of 6621   67.13 % -> 72.19 %
#         organic        4257 -> 4433 of 4539   93.79 % -> 97.66 %
#         chair against counter-chair on the metallacycle
#                          60 ->   79 of   80   75.00 % -> 98.75 %
#     The denominators stand still because they measure the RMSD and not the verdict --
#     exactly by that the change is recognisable as a pure predicate tightening.
#
# ⚠ WHY THE RING LIST AND NOT THE MATHEMATICS WAS THE PROBLEM.  Cremer-Pople
#   needs only a cyclic order; the metal is a ring atom like any
#   other, and ``_cp_theta_phi`` is purely geometric.  That is checked, not
#   assumed (read-back test over 9620 seatings: the formula delivers a
#   well-defined (Q, theta, phi) for every chelate ring).
#   ⚠ What the read-back test ALSO shows and what one needs to know: the SET
#     amplitude is not reached -- |Q_read - Q_set| median 0.155 A
#     against a set ``_amp(5)`` = 0.40.  That is NOT a metal effect in the
#     CP calculation, but the frozen coordination sphere: in a
#     five-ring M-D-X-Y-D THREE of the five ring atoms are in ``frozen``, only
#     X and Y move.  For the fingerprint that does not matter -- it reads
#     the ACHIEVED geometry, not the desired one.
#
# ⚠ THE TOLERANCE STAYS 180/max(8,2n), AND THAT AFTER A MEASUREMENT THAT SPEAKS
#   AGAINST IT.  Its derivation is half the sampling spacing of the equator that
#   ``_pucker_candidates`` draws (K = max(8,2n) phases).  That presumes that the
#   generator also REACHES the requested phase spacing.  On the metallacycle it
#   does not -- measured the CP distance of NEIGHBOURING equator candidates,
#   read back:
#       n=5  median 11.5 degrees against tolerance 18.0  ->  2678 of 3969 below
#       n=6  median 18.6 degrees against tolerance 15.0  ->   879 of 3377 below
#   The softness sits precisely in the M-D bonds, and those are frozen; the
#   generator samples more tightly than it believes.
#   ⇒ The tolerance is too COARSE for the metallacycle.  It is nevertheless NOT
#     tightened, and the reason is the DIRECTION of the error: ``_fold_same`` returns
#     True on "within the tolerance", and True means DUPLICATE,
#     i.e. exactly the old behaviour.  A too-coarse tolerance can only leave RESCUE
#     on the table, never wrongly keep a frame.  A tightened number,
#     by contrast, would be a fitted knob without derivation -- and the 4445 above are
#     measured with the coarse tolerance, hence a LOWER BOUND.
#
# ✅ WHAT THIS TERM COULD NOT DO -- CLOSED, NOT LEFT STANDING.  Until today this
#   stood here: ``_fold_same`` does not compare Q at all, only (theta, phi); if
#   ONE ring lies below the amplitude floor and the other above, the
#   great-circle distance judges by a meaningless phi.  Measured were 1191 such
#   pairs, of which 146 (12.26 %) judged as SAME -- the planar chelate-ring state
#   (`5M:planar`), which the eye counts as its own basin.
#   The amplitude block above ``_fold_same`` closes that: RE-MEASURED with
#   the same tool, archive, seed and draw, now **0 of 1191** (0.00 %) stand
#   on SAME.  It hit, as predicted here, the organic axis as well -- which is why
#   it is not a silent change, but one with newly collected numbers on BOTH
#   axes (see the block above ``_fold_same`` and the numbers above).
#
# ⛔ AND NOW THE UNCOMFORTABLE REACH QUESTION, CHECKED INSTEAD OF ASSUMED.
#   The occasion for this term was: "DELFIN_FFFREE_PUCKER_MC produces
#   metallacycle folds that subsequently run into the same complex RMSD."
#   THAT IS NOT TRUE.  Looked up in the emitter instead of believed:
#     `converter_backend._append_ffree_ring_puckers` calls `_ring_pucker.generate`
#     (PUCKER_MC sits there) and appends every surviving fold with
#     `results.append(...)` DIRECTLY -- :537.  These frames see none of the three
#     dedups; the only dedup over `results` is an EXACT
#     string comparison against the primary frame (:2066/:2740/:3065).  And the
#     call sites (:2098/:2160/:2800/:3103) stand AFTER the build, i.e. after the
#     dedup.
#   ⇒ This term does NOT protect the folds of the fold emitter.  It protects
#     the metallacycle folds that are there ALREADY BEFORE the dedup: those from
#     the conformer pool per ligand (`_ligand_confs_from_mol` -> the combination
#     product in `assemble_from_config`) and the eta variants in `_dedup_builds`.
#     That is a real axis -- the 8520 candidates above lie on it --, but
#     it is NOT the axis that triggered the task.
#   ⚠ And on the small byte battery the effect is ZERO, honestly counted
#     (`harness/faltung_fp_mc_vergleich.sh`, census of both runs):
#         MC=0   rings per build  12,  one ring list without verdict (None)
#         MC=1   rings per build  26,  no ring list without verdict any more
#                14 chelate rings in 7 ring lists, all n=5
#         `_fold_same` in BOTH runs: 15 SAME / 2 DIFFERENT
#     So the switch is not dark -- it doubles the ring set and gives
#     a build a verdict in the first place --, but NO frame changes and
#     not a single `_fold_same` verdict flips.  On those four chelate systems
#     the conformer product varies the ligand periphery, not the
#     chelate ring fold -- which fits the P3 measurement above: with nailed-down
#     donors the achieved fold spread on the metallacycle stays small.
#   ⇒ WHOEVER WANTS TO SWITCH THIS ON must measure it on a 5000-system A/B,
#     not on a battery -- and the number they may expect is the one from
#     the archive measurement, not the one from the byte run.
#
# Switch: DELFIN_FFFREE_DEDUP_FOLD_FP_MC (default 0).  Its OWN switch, although
# the parent switch is OFF anyway -- only that way do the two findings
# (organic / metallacycle) stay separately measurable.  ⚠ The two percentages that
# stood here -- 73.98 and 65.41 -- are superseded since the Q term; the valid
# numbers with their denominators stand above in ONE block, so that they cannot
# drift apart at three places.

_FOLD_FP_MCMAX = 8           # cost cap: metallacycles per frame (deterministically


                             # sorted).  Measured 816/314 = 2.6 per system.


def _fold_fp_mc_enabled():
    return (_fold_fp_enabled()
            and os.environ.get("DELFIN_FFFREE_DEDUP_FOLD_FP_MC", "0") == "1")


def _fold_mc_arms(lg):
    """The donors with which THIS ligand spans a chelate ring -- or
    ``None`` if it spans none.

    ⚠ eta ligands are EXCLUDED, and that is not a reservation but
      geometry: an eta face is not a sequence of sigma donors, but ONE
      pi bond.  ``assemble_hapto`` also freezes the whole ring set there
      (``fixed.update(...eta_local_idxs)``) and reports only ONE representative as
      donor (:3902).  A "ring" M-C1-C2-C3 would be the cut-open Cp ring, not a
      metallacycle -- and its middle atom is itself metal-bound, which the
      minimality test in ``_fold_mc_rings`` catches anyway.
    ⚠ ``_canonical_arm_order`` instead of ``donor_local_idxs``: that is the list the
      builder ACTUALLY seats on vertices (:3954 in the hapto branch, :5340 in the
      configuration branch).  ``donor_local_idxs`` can be longer than the denticity
      -- that would produce rings nobody has coordinated."""
    try:
        if lg.get("is_eta"):
            return None
        dent = int(lg.get("denticity") or 0)
        if dent < 2:
            return None                          # one donor spans no ring
        arms = [int(x) for x in _canonical_arm_order(lg, dent)]
        return arms if len(arms) >= 2 else None
    except Exception:
        return None


def _lig_path(mol, a, b, maxlen):
    """Shortest path ``a`` -> ``b`` IN THE LIGAND GRAPH, as a list of local indices.

    The ligand mol carries no metal, so the path cannot shortcut via the metal
    -- exactly the property that defines a chelate ring.  Neighbours are
    visited in ascending order -> deterministic for equally long paths.  H is skipped
    as an intermediate atom (it is terminal and can never lie on a
    shortest path between two heavy atoms)."""
    a, b = int(a), int(b)
    if a == b:
        return None
    vor = {a: None}
    dq = [(a, 1)]
    head = 0
    while head < len(dq):
        cur, tiefe = dq[head]
        head += 1
        if tiefe >= maxlen:
            continue
        for nb in sorted(int(x.GetIdx()) for x in mol.GetAtomWithIdx(cur).GetNeighbors()):
            if nb in vor:
                continue
            if nb != b and mol.GetAtomWithIdx(nb).GetSymbol() == "H":
                continue
            vor[nb] = cur
            if nb == b:
                weg = [nb]
                while vor[weg[-1]] is not None:
                    weg.append(vor[weg[-1]])
                weg.reverse()
                return weg
            dq.append((nb, tiefe + 1))
    return None


def _fold_mc_rings(blocks, syms, metal_idx=0):
    """Global ring index lists of the CHELATE RINGS, in cyclic order.

    ``blocks`` = sequence of ``(global_offset, mol, donor_locals)``.  The ring is
    ``[metal_idx] + shortest ligand path(d1 -> d2)`` -- that is cyclic order
    by construction (M-d1, the path bonds, d2-M), i.e. exactly what
    ``_cp_theta_phi`` demands.

    ⚠ MINIMALITY INSTEAD OF ALL PAIRS.  A tridentate has three
      donor pairs, but only two chelate rings: the third pair runs via the
      middle donor and is the CONCATENATION of the two.  A path that touches a
      third donor of the same ligand is therefore rejected.

    ⚠ THE METAL INDEX IS CHECKED, NOT BELIEVED.  All callers write the
      metal to index 0 (``out_syms = [metal]`` :3262/:5042, ``off = 1`` in
      ``_hapto_fold_rings``), but an offset model is an assumption -- so
      ``_elements.is_metal`` stands in front here, and every ring atom is checked,
      as in ``_fold_rings_from_blocks``, against its element symbol in the frame.
      If anything does not fit: ``None``, and the predicate falls back to the pure
      RMSD.

    ⚠ AROMATIC IS A DIFFERENT CRITERION HERE THAN ON THE ORGANIC RING.  There
      ONE aromatic ring atom tips the whole ring (flat, rigid face).
      A salen chelate ring has aromatic phenolate carbons and folds
      nonetheless -- at the M-D bonds.  Therefore only what is ENTIRELY
      aromatic apart from the metal is excluded -- bipyridine, terpyridine,
      metallabenzene.  That is literally the rule that
      ``_ring_pucker._is_puckerable`` applies under DELFIN_FFFREE_PUCKER_MC,
      not a second opinion on it, and it is THE RIGHT WAY ROUND HERE: the
      generator never folds such rings, so there is nothing to rescue there either.
      A fingerprint that carried them would be stricter than the generator and
      would keep frames that differ in nothing.
      MEASURED what the rule costs (`harness/faltung_fp_mc_selbsttest.py` and
      the 500-system measurement): 110 of 816 chelate rings drop out, all of the
      bipyridine type; salicylaldiminate (M-O-C(ar)-C(ar)-C=N) stays in, because N
      and O are not aromatic -- exactly the case the `_ring_pucker`
      comment calls "step or umbrella fold"."""
    try:
        from delfin.manta import _elements as _EL
    except Exception:
        return None
    try:
        mi = int(metal_idx)
        if mi < 0 or mi >= len(syms) or not _EL.is_metal(syms[mi]):
            return None                          # no metal there: no verdict
        ringe = []
        gesehen = set()
        for off, mol, dons in blocks:
            if mol is None or not dons or len(dons) < 2:
                continue
            off = int(off)
            dset = {int(d) for d in dons}
            arme = sorted(dset)
            for ai in range(len(arme)):
                for bi in range(ai + 1, len(arme)):
                    weg = _lig_path(mol, arme[ai], arme[bi], _FOLD_FP_RINGMAX)
                    if not weg or len(weg) < 3:
                        continue                 # no path / ring size < 4
                    if len(weg) + 1 > _FOLD_FP_RINGMAX:
                        continue
                    if any(j in dset for j in weg[1:-1]):
                        continue                 # third donor: not minimal
                    if all(mol.GetAtomWithIdx(int(j)).GetIsAromatic() for j in weg):
                        continue                 # metallabenzene: planar-rigid
                    g = [mi] + [off + int(j) for j in weg]
                    if min(g) < 0 or max(g) >= len(syms):
                        return None
                    for j, gj in zip(weg, g[1:]):
                        if syms[gj] != mol.GetAtomWithIdx(int(j)).GetSymbol():
                            return None
                    schl = frozenset(g)
                    if schl in gesehen:
                        continue
                    gesehen.add(schl)
                    ringe.append(tuple(g))
    except Exception:
        return None
    if not ringe:
        return None
    ringe.sort()                                 # deterministic
    return ringe[:_FOLD_FP_MCMAX]


def _fold_rings_with_mc(blocks, syms, metal_idx=0):
    """The fingerprint's ring set: organic rings as before, plus the
    chelate rings when ``DELFIN_FFFREE_DEDUP_FOLD_FP_MC`` is on.

    ``blocks`` = ``(offset, mol, donor_locals)``; ``donor_locals=None`` switches
    the metallacycle off for this block.  ⛔ MC switch OFF -> the return
    value is literally ``_fold_rings_from_blocks(...)``, i.e. byte-identical
    to the state in which the organic measure was collected."""
    org = _fold_rings_from_blocks([(o, m) for (o, m, _d) in blocks], syms)
    if not _fold_fp_mc_enabled():
        return org
    mc = _fold_mc_rings(blocks, syms, metal_idx)
    if not mc:
        return org
    if org is None:
        return mc
    # ⚠ APPENDED, NOT SORTED IN.  `_fold_same` runs with `zip` over two
    #   fingerprints that come from THE SAME ring list; any stable
    #   order will do, but it must be the same between the frames.
    return org + mc


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
