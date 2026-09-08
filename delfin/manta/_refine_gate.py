"""_refine_gate.py — THE ROLLBACK GATE.  A FILTER, not a corrector.

WHY (24.08.2026, from the user's question).  A refiner that moves atoms wins on many systems
and loses on a few.  The landing gate knows 33 forbidding terms, all on the LOSS SIDE, and
`capability_gained` is computed, printed -- and feeds into not a single blocking decision.
A net-positive refiner can therefore never land.

THE ANSWER IS NOT TO LOOSEN THE GATE.  A mean cannot see a dead class, and the eye itself
changes -- never-worse is the only criterion that survives a change of eye.  The answer is to
TAKE THE LOSS AWAY from the refiner:

    after the chain, NO frame is worse than before.

Then "move atoms so that they sit better" is exactly what happens -- without collateral
damage, and never-worse holds BY CONSTRUCTION instead of on average.

⛔ THREE RULES THAT DETERMINE THE DESIGN (user, 24.08.):

 1. **NO CORRECTOR OF A CORRECTOR.**  This module NEVER writes coordinates.  It chooses
    between two frames that already exist -- result or original.  A corrector that fixes up
    a corrector is a sign that one of the two is wrong; then that one should be made right,
    not a third one placed beside it.

 2. **NO BLOAT.**  ONE public function, ONE module, and TWO lines per refinement chain --
    not per refiner.  There are more than twenty `_apply_*` sites; wrapping each one
    individually would be exactly the bloat that is ruled out here.
    And NO radii, thresholds or graphs are redefined: everything comes from
    `_h_placement`, where it is calibrated and justified against the detectors.

 3. **SPEED.**  The point at which this becomes cheap: a refiner leaves most frames
    BYTE-IDENTICAL.  A string comparison up front costs nothing and excludes them.  Only
    what has really changed is evaluated -- typically a handful of frames per system
    instead of all of them.  Without this pre-filter the gate would be more expensive than
    the refiners it guards.

THE CRITERION -- reference-free, local, TWO INTEGERS, both "smaller is better":

    n_clash      atom pairs that fall below their floor (H...H 1.50 A; H...heavy
                 0.85 x vdW sum; heavy...heavy 0.70 x vdW sum), only for truly
                 non-bonded pairs from graph distance `_NONBONDED_MIN_HOPS` on
    n_bond_out   covalent bonds outside [0.85 .. 1.15] x sum of covalent radii

A frame counts as NOT WORSE if NEITHER of the two numbers rises.  Number against number,
i.e. dimensionally clean -- no comparison of a severity with a frame count, of the kind that
has led this repo astray three times already.

⚠ THE CRITERION IS NOT THE EYE, and that is deliberate.  Whoever rebuilds the eye inside the
builder hands Goodhart the keys: the builder then optimizes the measurement instead of the
geometry.  Only physical invariants stand here, which the builder can check from itself.

⚠ NEW LABELS PASS THROUGH UNTOUCHED.  Enumerators (mirror, fold, atropisomer) APPEND;
they cannot lose anything by construction, and an appended sibling MAY be worse than the
original -- that is completeness, not a regression.  The matching runs over the LABEL, so
that this module keeps the two mechanism classes apart without knowledge of the individual
passes.

⚠ NAMING `_rg_`: the watchdog `check_exists_first` reported `_enabled` and `_score` as a
collision (`joint_declash`, `sphere_flex`, `smiles_converter`).  The same lesson as with
`_hp_` and `_atrop_` -- an unambiguous name makes the copy visible instead of disguising it.

Default OFF (``DELFIN_FFFREE_REFINE_GATE``) -> byte-identical.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

# NO radii, floors or graphs of our own: everything from `_h_placement`, where it is calibrated.
from delfin.manta._h_placement import (
    _HH_FLOOR,
    _H_HEAVY_FRAC,
    _NONBONDED_MIN_HOPS,
    _XH_FIRE_HI,
    _XH_FIRE_LO,
    _hp_cov,
    _hp_graph,
    _hp_metal,
    _hp_read,
    _hp_vdw,
    _hp_within,
)

_LOG = logging.getLogger(__name__)

FLAG = "DELFIN_FFFREE_REFINE_GATE"

# Heavy-against-heavy is the ONLY floor that `_h_placement` does not need (it only moves
# H).  0.70 x vdW sum is deliberately LOOSE: a real bond lies far below it and is excluded
# via the graph distance anyway; only gross interpenetrations are meant.
_RG_HEAVY_HEAVY_FRAC = 0.70


def _rg_on() -> bool:
    """THE one read site.  A switch that is read twice drifts."""
    return os.environ.get(FLAG, "0") == "1"


# ── THE THIRD METRIC: sp2 centres that were built pyramidal ──────────────────────
#
# WHY (25.08.2026).  The defect ranking over 79 495 frames (aromrad6k_off, 1618
# systems) says that the two LARGEST classes are the same chemistry:
#     smiles_hyb-angle   22.18 %   crystal 0.00 %   mass 3911
#     pyramidal_sp2      20.02 %   crystal 0.98 %   mass 3028
# Together 6939 -- more than the next three classes combined, and the crystal
# sits at zero: pure build error, no detector noise.
#
# `_rg_score` was BLIND to both: it knew only clashes and bond lengths.
# So the rollback gate could protect every repair EXCEPT the one on the
# largest class.  A flattener that straightens 300 centres and bends 30 was
# waved through by the gate and struck down by the landing gate.
#
# ⚠ THE MEASURE IS THE EYE'S, NOT MY OWN.  At first I wanted to take
# `_frame_assertions._oop` -- the distance of the centre from the plane of its
# three neighbours, in ANGSTROM.  But the eye measures the WALSH ANGLE in DEGREES
# (`find_pyramidalization.py:515`):
#       walsh_deg = atan2(d_oop, r_mean)
# i.e. exactly `_oop`, NORMALIZED to the mean centre-neighbour bond length.  An
# Angstrom measure with a degree threshold would again have been a label instead
# of a measurement -- the mistake that kept the atropisomer enumerator crippled
# until 23.08.  That is why the eye's formula stands here, character for character.
#
# ⚠ WHAT IS DELIBERATELY COARSER HERE THAN THE EYE, and why that suffices: the eye
# chooses its ceiling ENVIRONMENT-DEPENDENTLY (`_PYR_ENV`, p99.9 from 307k crystals).
# This gate, by contrast, COMPARES the same frame before and after the chain -- the
# same atoms, the same environment.  For a "has it become worse" a CONSISTENT
# measure suffices; the absolute calibration only decides which centres count at all.
# The floor 12.0 degrees is nevertheless the eye's (`find_graph_geometry._PYR_BUILD_MIN`),
# so that the numbers remain comparable.
#
# ⚠ mu-BRIDGES ARE NOT PYRAMIDAL BUT BRIDGING -- the eye's R3 rule
# (`full_verdict.py:514-525`, established on CCDC).  A light atom bonded to TWO or more
# metals is skipped; otherwise every mu2-oxo corner would count as a defect.
#
# ⚠ SEPARATE SWITCH, and that is not caution but obligation: `picopgate6k` and
# `hplacegate6k` have been running SINCE NOON TODAY with the two-tuple, and the builder
# re-reads its file for EVERY system.  Without a separate switch this line would have
# altered two running A/Bs mid-run.  Default OFF -> two-tuple, byte-identical.
FLAG_PYR = "DELFIN_FFFREE_REFINE_GATE_PYR"
_RG_PYR_FLOOR = 12.0                      # degrees -- find_graph_geometry._PYR_BUILD_MIN
_RG_PYR_ELEMS = frozenset(("C", "N", "O", "S", "P"))


def _rg_pyr_on() -> bool:
    """THE one read site for the third metric."""
    return os.environ.get(FLAG_PYR, "0") == "1"


def _rg_walsh(P, c: int, nb: Sequence[int]) -> Optional[float]:
    """Walsh angle in degrees, formula from `find_pyramidalization.py:501-515`.
    Planar sp2 -> ~0; the more pyramidal, the larger."""
    a, b, d = int(nb[0]), int(nb[1]), int(nb[2])
    nrm = np.cross(P[b] - P[a], P[d] - P[a])
    ln = float(np.linalg.norm(nrm))
    if ln < 1e-9:
        return None
    d_oop = abs(float(np.dot(P[c] - P[a], nrm / ln)))
    r_mean = (float(np.linalg.norm(P[c] - P[a]))
              + float(np.linalg.norm(P[c] - P[b]))
              + float(np.linalg.norm(P[c] - P[d]))) / 3.0
    if r_mean < 1e-6:
        return None
    return math.degrees(math.atan2(d_oop, r_mean))


def _rg_pyr_count(syms: Sequence[str], P, adj) -> int:
    """How many THREE-BONDED light-atom centres are pyramidalized beyond the floor?

    Three-bonded, because an sp2 centre has exactly three sigma neighbours -- H counts,
    as with the eye ("including every sp2 C-H", `full_verdict.py:507`).  Without bond
    orders, "three neighbours on C/N/O/S/P" is the best available sp2 approximation;
    a true sp3 centre has four and thus drops out by itself."""
    n_pyr = 0
    for i in range(len(syms)):
        if syms[i] not in _RG_PYR_ELEMS:
            continue
        nb = list(adj[i])
        if len(nb) != 3:
            continue
        if sum(1 for j in nb if _hp_metal(syms[j])) >= 2:
            continue                      # mu-bridge, not a mis-planarization (R3)
        w = _rg_walsh(P, i, nb)
        if w is not None and w > _RG_PYR_FLOOR:
            n_pyr += 1
    return n_pyr


def _rg_score(xyz: str) -> Tuple[int, ...]:
    """(n_clash, n_bond_out[, n_pyr]) -- all "smaller is better".

    The third position is added only with `DELFIN_FFFREE_REFINE_GATE_PYR=1`; without it
    the return value is the old two-tuple.  The comparison in `keep_better` is a tuple
    comparison and carries both lengths -- as long as BOTH sides come from the same
    run, and they do by construction (same process, same environment).
    `(-1, -1)` or `(-1, -1, -1)` means unreadable."""
    _pyr = _rg_pyr_on()
    try:
        syms, P, _lines = _hp_read(xyz)
    except Exception:
        return (-1, -1, -1) if _pyr else (-1, -1)
    n = len(syms)
    if n < 2:
        return (0, 0, 0) if _pyr else (0, 0)
    adj = _hp_graph(syms, P)
    n_bond_out = 0
    for i in range(n):
        for j in adj[i]:
            if j <= i:
                continue
            tgt = _hp_cov(syms[i]) + _hp_cov(syms[j])
            if tgt <= 1e-6:
                continue
            r = float(np.linalg.norm(P[i] - P[j])) / tgt
            if r < _XH_FIRE_LO or r > _XH_FIRE_HI:
                n_bond_out += 1
    # Clashes: only truly non-bonded pairs, otherwise every normal 1,3-contact
    # reports a violation.
    nah: List[Set[int]] = [_hp_within(adj, i, _NONBONDED_MIN_HOPS - 1) for i in range(n)]
    n_clash = 0
    for i in range(n):
        si = syms[i]
        if _hp_metal(si):
            continue                          # the M-D sphere is not a clash question
        for j in range(i + 1, n):
            if j in nah[i]:
                continue
            sj = syms[j]
            if _hp_metal(sj):
                continue
            if si == "H" and sj == "H":
                thr = _HH_FLOOR
            elif si == "H" or sj == "H":
                thr = _H_HEAVY_FRAC * (_hp_vdw("H") + _hp_vdw(sj if si == "H" else si))
            else:
                thr = _RG_HEAVY_HEAVY_FRAC * (_hp_vdw(si) + _hp_vdw(sj))
            if float(np.linalg.norm(P[i] - P[j])) < thr:
                n_clash += 1
    if _pyr:
        return (n_clash, n_bond_out, _rg_pyr_count(syms, P, adj))
    return (n_clash, n_bond_out)


def keep_better(before: Sequence, after: Sequence):
    """After the refinement chain, NO frame is worse than before.

    `before` is the list BEFORE the chain, `after` the one after it; both (xyz, label).
    Returned is a list of the same length as `after`:
      * label existed before AND the text has changed AND the score is
        worse -> the ORIGINAL,
      * otherwise the result of the chain, unchanged.

    Byte-identical frames are excluded via a string comparison BEFORE anything is
    computed -- the reason why this gate costs almost nothing in the build.
    """
    if not _rg_on() or not before or not after:
        return after

    def _rg_lbl(e):
        try:
            return e[1] if len(e) > 1 else ""
        except Exception:
            return ""

    def _rg_xyz(e):
        try:
            return e[0]
        except Exception:
            return None

    vorher: Dict[str, str] = {}
    for e in before:
        l, x = _rg_lbl(e), _rg_xyz(e)
        if l and x is not None and l not in vorher:
            vorher[l] = x

    out = list(after)
    ZAEHLER["keep_better_gelaufen"] += 1
    n_rueck = n_geprueft = 0
    for k, e in enumerate(out):
        l, x = _rg_lbl(e), _rg_xyz(e)
        alt = vorher.get(l)
        if alt is None or x is None or alt == x:
            continue                          # new (enumerator) or unchanged -> free
        n_geprueft += 1
        s_alt, s_neu = _rg_score(alt), _rg_score(x)
        # ⚠ BOTH LINES WERE NAILED TO THE TWO-TUPLE (noticed 25.08., before the
        # third metric ever ran).  `== (-1, -1)` is FALSE for `(-1, -1, -1)`,
        # so an unreadable frame would no longer have been skipped; and the
        # comparison read only index 0 and 1, so `n_pyr` would have been computed
        # and then THROWN AWAY -- a counter without an effect line, exactly the
        # design that has ended up as a "dark switch" here several times already.
        # Now length-independent.
        if -1 in s_alt or -1 in s_neu:
            continue                          # unreadable -> do not judge, pass through
        if len(s_alt) != len(s_neu):
            continue                          # can only happen on a switch change in the
            # MIDDLE of a run; then no comparison is possible and passing through is right.
        if any(b > a for a, b in zip(s_alt, s_neu)):
            try:
                out[k] = ((alt,) + tuple(e[1:])) if isinstance(e, tuple) else ([alt] + list(e[1:]))
            except Exception:
                continue
            n_rueck += 1
    ZAEHLER["keep_better_zurueck"] += n_rueck
    ZAEHLER["keep_better_geprueft"] += n_geprueft
    if n_rueck:
        # WARNING, NOT INFO (01.09.2026).  Measured: `refine-gate` appears in
        # **0 of 2031** run logs -- the builder's INFO level does not reach the
        # run log; `mirror_enum` is visible there only because it uses `warning`.
        # This gate was therefore UNOBSERVABLE since it was built: its effect
        # was assumed, never seen.
        # ⚠ And the consequence was worse than missing curiosity: with `addroot3` I
        #   could NOT tell from the silence whether the gate found nothing or did
        #   not run at all.  Exactly the distinction that `loop.py:_fire_out` raises
        #   to a rule with "measured and hit nothing is a DIFFERENT statement than
        #   not measured".
        # ⇒ An event is reported; from now on silence means "nothing done".
        _LOG.warning("refine-gate: %d von %d veraenderten Frames zurueckgenommen "
                     "(Kollisionen oder Bindungslaengen wurden schlechter)",
                     n_rueck, n_geprueft)
    return out


FLAG_ADD = "DELFIN_FFFREE_ADD_NEVER_REPLACE"

# ── COUNTERS INSTEAD OF LOG (01.09.2026) ─────────────────────────────────────────
# MEASURED: `refine-gate` appears in 0 of 2031 run logs.  I thereupon raised the
# messages from INFO to WARNING -- and re-measured on `addroot4`: NO builder
# warning at all reaches the run log, not even `mirror_enum`, which certainly
# fired.  The builder runs as a subprocess whose log stream is discarded; the
# log is the wrong channel, no matter at which level.
#
# THE RIGHT CHANNEL is the JSON line the builder writes per system anyway
# (loop.py:406).  Next to it stands `_fire_out()` with exactly this justification:
# "ALWAYS pass it along, even empty: 'measured and hit nothing' is a DIFFERENT
#  statement than 'not measured'".
#
# ⚠ `gelaufen` (ran) is counted INDEPENDENTLY of `getroffen` (hit).  Without that a
#   zero is not readable -- and that is exactly where the question "did the gate find
#   nothing, or did it not run at all?" failed on `addroot3` and `addroot4`.
ZAEHLER = {
    "keep_better_gelaufen": 0,      # calls in which the gate was ON
    "keep_better_zurueck": 0,       # frames actually rolled back
    # How many frames were actually COMPARED. Measured 2026-09-08 on picopret6k: the
    # gate ran 905 times and rolled back nothing, and from the outside "compared 500,
    # none worse" was indistinguishable from "compared none, the gate never saw the
    # axis at all". It was the second: the axis it was meant to guard runs in the
    # caller, after this window has already closed. A gate that cannot say how much
    # it looked at cannot be told apart from one that is wired past its subject.
    "keep_better_geprueft": 0,      # frames compared (label known, text changed)
    "keep_all_gelaufen": 0,         # calls in which the gate was ON
    "keep_all_wieder": 0,           # frames actually restored
    # ── Dual-parse union (smiles_converter.py:~35092), 01.09.2026 ────────────────
    # For metal SMILES with `canonical != input` the build runs TWICE; the two
    # result sets are united via an XYZ signature.  This signature LEAVES OUT
    # HYDROGENS and sorts the heavy-atom lines -- two frames that differ only in
    # H positions (stereocentre against mirror image) get the same key.  The
    # assignment keeps the LATER one.  This is the only stage found with a
    # reversed preference rule and an H-blind key.
    # ⚠ HERE ONLY COUNTING HAPPENS, NOTHING IS CHANGED.  The thesis is proven only
    #   once `dual_sig_kollision` is non-zero on ABUSAU/JEJROI AND disappears
    #   with the mirror.
    "dual_parse_gelaufen": 0,       # second build executed at all
    "dual_sig_kollision": 0,        # frames that OVERWRITE an existing key
    # ── Pairwise gate for appended frames (converter_backend, 01.09.2026) ──────
    # `gelaufen` (ran) separate from `verworfen` (rejected), for the same reason as
    # above: a zero with gelaufen>0 means "checked, nothing to reject", a zero with
    # gelaufen==0 means "the gate did not run".  Without the separation it is unreadable.
    "pairgate_gelaufen": 0,         # calls of `_append_reembed` with gate ON
    "pairgate_verworfen": 0,        # appended frames that brought a NEW too-close pair
}


def _ka_on() -> bool:
    """THE one read site -- like `_rg_on`.  A switch, read twice, drifts."""
    return os.environ.get(FLAG_ADD, "0") == "1"


def keep_all(before: Sequence, after: Sequence):
    """ADD, NEVER REPLACE: no frame from `before` may be MISSING at the end.

    Sibling of `keep_better`.  That one enforces *"no frame is worse than
    before"*, this one *"no frame is GONE"* -- two halves of the same
    contract, and the second was not enforced until today.

    OCCASION (01.09.2026, measured on `mirrleg6k`).  `harness/frame_keys_additiv.py`
    checked the eight blockers of the best landing candidate at CONTENT LEVEL:

        LIYGAC QOYTEE XUFHEM XUFHIQ FEDDEA VAPNEI   nur_basis 0   strictly additive
        ABUSAU  58 -> 59 frames                     nur_basis 1   CONTENT GONE
        JEJROI  89 -> 90 frames                     nur_basis 1   CONTENT GONE

    On ABUSAU and JEJROI ONE frame each disappears, although `_mirror_enum.
    expand_results` is additive by construction (`:370  return list(results) +
    added`; the only other exit returns `results` unchanged).
    So the loss arises DOWNSTREAM, in a stage that SELECTS.
    The labels say what it hits:

        ABUSAU  gone: ...Br-N2-D-conf4_stereo-u    new: ...Br-N2-L-conf3_mirror (2x)
        JEJROI  gone: all-cis-L-conf3-2_stereo-uu  new: Isomer 2_mirror         (2x)

    In both cases a STEREOCENTRE frame falls while mirrors are added --
    and mirroring inverts Lambda/Delta, so the new frames are close relatives
    of the fallen ones.  Two additive passes displace each other through a
    selection stage.  The same design cost `trans208` a CCDC isomer (note at
    the mirror call site, smiles_converter.py:32621).

    WHAT THIS GATE DOES NOT DO.  It does NOT say which stage takes the frame -- it
    restores it.  The root remains open and is tracked as a separate task; a
    contract that is enforced only at the end of the chain is a SEAM and
    not a cure.

    AND IT CAN BRING BACK A FRAME THAT A SELECTION DELIBERATELY REJECTED.
    That is intended: the contract says an ENUMERATOR may not cost a frame.
    Whoever wants a selection must make it BEFORE the enumerator, not after.
    Whoever sees it differently leaves the gate off -- it is default OFF.

    ⚠️ CORRECTION 01.09.2026, later in the day: `_gfnff_ensemble_rank_filter` stood
       here as the suspect, "on in the champion".  THAT WAS WRONG.  `_CHAMPION_FLAGS`
       (cli_manta.py:42-281) has 20 entries, `GFNFF_RANK` is not one of them; my
       hit came from comment text next to the definition.  The culprit is
       still UNKNOWN -- and `harness/kettenzaehler.py` has meanwhile ruled out the
       whole OUTER chain (14 stages): the frame count falls nowhere there.
       So what this gate does does not depend on that conjecture.

    MATCHING OVER THE LABEL MULTISET, not over the text: stages in between
    reformat coordinates; a string comparison then reported "gone" where it was
    only printed differently, and the gate appended DUPLICATES.  Labels are stable
    (passes append suffixes, they do not rewrite).  MULTISET, because duplicate
    labels are the normal case -- ABUSAU carries 20 of them in the base arm; a
    SET instead of a counter would never miss the second frame.

    Byte-identical as long as `DELFIN_FFFREE_ADD_NEVER_REPLACE != 1`.
    """
    if not _ka_on() or not before or not after:
        return after
    ZAEHLER["keep_all_gelaufen"] += 1
    try:
        def _lbl(e):
            try:
                return e[1] if len(e) > 1 else ""
            except Exception:
                return ""

        n_vor: Dict[str, int] = {}
        for e in before:
            n_vor[_lbl(e)] = n_vor.get(_lbl(e), 0) + 1
        n_nach: Dict[str, int] = {}
        for e in after:
            n_nach[_lbl(e)] = n_nach.get(_lbl(e), 0) + 1

        fehlt = {l: n - n_nach.get(l, 0) for l, n in n_vor.items()
                 if n - n_nach.get(l, 0) > 0}
        if not fehlt:
            return after

        # WHICH occurrence to restore?  The label says HOW MANY are missing,
        # not WHICH.  If a label appears twice in `before` and once in
        # `after`, and one simply takes the first, one appends a DUPLICATE of
        # the frame that is present -- and the one really missing stays lost.
        # (Exactly this is what the first version failed on in the self-test.)
        # So: contents that `after` already carries are struck off once; what
        # remains afterwards is restored preferentially.  The text comparison is
        # fit for use HERE because it only SELECTS -- whether anything is missing
        # at all has already been decided by the label multiset.
        vorhanden: Dict[str, int] = {}
        for e in after:
            try:
                _k = str(e[0])
            except Exception:
                continue
            vorhanden[_k] = vorhanden.get(_k, 0) + 1

        rest = list(after)
        offen = dict(fehlt)
        nachrang = []
        for e in before:
            l = _lbl(e)
            if offen.get(l, 0) <= 0:
                continue
            try:
                _k = str(e[0])
            except Exception:
                _k = None
            if _k is not None and vorhanden.get(_k, 0) > 0:
                vorhanden[_k] -= 1          # this content is already there
                nachrang.append(e)          # keep only as a fallback
                continue
            rest.append(e)
            offen[l] -= 1
        # Fallback: if slots remain open (all candidates were identical in content),
        # the count must still be honoured -- otherwise the gate would be silent
        # depending on the data.
        for e in nachrang:
            l = _lbl(e)
            if offen.get(l, 0) > 0:
                rest.append(e)
                offen[l] -= 1
        # NO SILENT RESTORATION, and WARNING instead of INFO -- for the same
        # measured reason as with `keep_better` above: INFO from the builder does not
        # reach the run log (0 of 2031).  With `addroot3` exactly that robbed me of
        # the answer whether this gate found nothing or did not run.
        ZAEHLER["keep_all_wieder"] += sum(fehlt.values())
        _LOG.warning("ADD-never-replace: %d Frame(s) wiederhergestellt, die die "
                     "Kette verloren hatte (%d Etikett(en): %s)",
                     sum(fehlt.values()), len(fehlt), ", ".join(sorted(fehlt)[:4]))
        return rest
    except Exception as exc:                  # pragma: no cover - fail-safe
        _LOG.warning("ADD-never-replace nicht angewandt (%s) -- es wird NICHTS "
                     "wiederhergestellt", type(exc).__name__)
        return after


def _rg_selbsttest() -> int:  # pragma: no cover - tool, not a production path
    """This module's claims against numbers, not against confidence.

        python -m delfin.manta._refine_gate

    What is checked is what CAN go wrong, not what is convenient:
      1. Default OFF  -> `_rg_score` returns a TWO-tuple.  Without that, the
         running A/Bs `picopgate6k` and `hplacegate6k` would have been altered mid-run.
      2. Switch ON    -> THREE-tuple, and the third number DISTINGUISHES flat from
         pyramidal.  A counter that says the same on both forms is a
         table constant and not a detector.
      3. The Walsh angle matches the EYE's numbers: planar formaldehyde ~0 degrees,
         and the formula is `atan2(d_oop, r_mean)`, not `d_oop` alone.
      4. `keep_better` ROLLS BACK a worsening of the third number -- otherwise
         it would be computed and thrown away.
    """
    import numpy as _np

    def _xyz(rows):
        return "%d\ntest\n" % len(rows) + "\n".join(
            f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, x, y, z in rows) + "\n"

    # formaldehyde-like: C with three neighbours, exactly planar (z = 0 for all)
    flach = _xyz([("C", 0.0, 0.0, 0.0), ("O", 0.0, 1.21, 0.0),
                  ("H", 0.94, -0.54, 0.0), ("H", -0.94, -0.54, 0.0)])
    # the same C, but pulled 0.35 A out of the plane -> clearly pyramidal
    pyr = _xyz([("C", 0.0, 0.0, 0.35), ("O", 0.0, 1.21, 0.0),
                ("H", 0.94, -0.54, 0.0), ("H", -0.94, -0.54, 0.0)])

    fehler = 0

    def _urteil(name, ok, zusatz=""):
        nonlocal fehler
        print(f"  {'OK  ' if ok else 'FEHL'} {name}{('  ' + zusatz) if zusatz else ''}")
        if not ok:
            fehler += 1

    _alt = os.environ.get(FLAG_PYR)
    _alt_haupt = os.environ.get(FLAG)
    try:
        os.environ[FLAG_PYR] = "0"
        s_aus = _rg_score(flach)
        _urteil("Vorgabe AUS liefert ein Zwei-Tupel", len(s_aus) == 2, f"-> {s_aus}")

        os.environ[FLAG_PYR] = "1"
        s_flach = _rg_score(flach)
        s_pyr = _rg_score(pyr)
        _urteil("Schalter AN liefert ein Drei-Tupel", len(s_flach) == 3, f"-> {s_flach}")
        _urteil("flaches sp2 zaehlt NICHT als pyramidal",
                len(s_flach) == 3 and s_flach[2] == 0, f"n_pyr={s_flach[2] if len(s_flach)>2 else '?'}")
        _urteil("gezogenes sp2 zaehlt SEHR WOHL",
                len(s_pyr) == 3 and s_pyr[2] >= 1, f"n_pyr={s_pyr[2] if len(s_pyr)>2 else '?'}")
        _urteil("die beiden Formen sind UNTERSCHEIDBAR (kein konstanter Zaehler)",
                len(s_flach) == 3 and len(s_pyr) == 3 and s_flach[2] != s_pyr[2])

        # Walsh angle recomputed against the eye's formula
        P = _np.array([[0.0, 0.0, 0.35], [0.0, 1.21, 0.0],
                       [0.94, -0.54, 0.0], [-0.94, -0.54, 0.0]])
        w = _rg_walsh(P, 0, [1, 2, 3])
        _urteil("Walsh-Winkel liegt ueber dem Boden", w is not None and w > _RG_PYR_FLOOR,
                f"{w:.1f} Grad, Boden {_RG_PYR_FLOOR}")
        Pf = _np.array([[0.0, 0.0, 0.0], [0.0, 1.21, 0.0],
                        [0.94, -0.54, 0.0], [-0.94, -0.54, 0.0]])
        wf = _rg_walsh(Pf, 0, [1, 2, 3])
        _urteil("planar ergibt ~0 Grad", wf is not None and wf < 1.0, f"{wf:.2f} Grad")

        # keep_better must roll back the worsening of the THIRD number.
        # ⚠ TWO SWITCHES, and the first test draft set only ONE -- the check
        # failed and looked like a code error.  `keep_better` bails out at line
        # 258 FIRST on `_rg_on()`, i.e. on the MAIN SWITCH; `FLAG_PYR` alone
        # does nothing.  This is not test cosmetics: a run that sets only
        # DELFIN_FFFREE_REFINE_GATE_PYR measures NOTHING and would report it as
        # "no effect".  The queue line must carry BOTH.
        os.environ[FLAG] = "1"
        vor = [(flach, "f0")]
        nach = [(pyr, "f0")]
        zurueck = keep_better(vor, nach)
        _urteil("keep_better nimmt die Pyramidalisierung ZURUECK",
                bool(zurueck) and zurueck[0][0] == flach)

        # and when OFF it must NOT do exactly that (the proof that it was due to the
        # third number and not to clashes or bond lengths)
        os.environ[FLAG_PYR] = "0"
        zurueck_aus = keep_better(vor, nach)
        _urteil("bei Vorgabe AUS bleibt dieselbe Aenderung STEHEN",
                bool(zurueck_aus) and zurueck_aus[0][0] == pyr,
                "sonst kaeme die Ruecknahme woanders her")
    finally:
        for _f, _v in ((FLAG_PYR, _alt), (FLAG, _alt_haupt)):
            if _v is None:
                os.environ.pop(_f, None)
            else:
                os.environ[_f] = _v

    # ── ADD, NEVER REPLACE ──────────────────────────────────────────────────────
    # ABUSAU is rebuilt: duplicate labels are the normal case (20 in the
    # base arm), ONE frame disappears, two are added.  A gate that can only handle
    # the simple case falls over exactly here.
    _alt_add = os.environ.get(FLAG_ADD)
    try:
        vor = [("A", "L-conf3"), ("B", "L-conf3"), ("C", "D-conf4_stereo-u")]
        nach = [("A", "L-conf3"), ("B", "L-conf3"),
                ("M1", "L-conf3_mirror"), ("M2", "L-conf3_mirror")]

        os.environ[FLAG_ADD] = "0"
        _urteil("bei Vorgabe AUS ist es ein NO-OP (byte-identisch)",
                keep_all(vor, nach) is nach)

        os.environ[FLAG_ADD] = "1"
        her = keep_all(vor, nach)
        _urteil("das verschwundene Etikett ist wieder da",
                any(l == "D-conf4_stereo-u" for _x, l in her))
        _urteil("und zwar mit dem ORIGINALINHALT, nicht mit einem Ersatz",
                any(x == "C" for x, _l in her))
        _urteil("die Spiegel bleiben unangetastet",
                sum(1 for _x, l in her if l == "L-conf3_mirror") == 2)
        _urteil("nichts wird doppelt angehaengt (Multiset, nicht Menge)",
                len(her) == 5, "4 vorhandene + genau 1 wiederhergestellter")

        # THE COUNTER-CHECK, without which the number would be worth nothing: if
        # NOTHING is missing, the gate must not touch anything either -- otherwise
        # the manifold grows on every call.
        _urteil("fehlt nichts, bleibt die Liste unveraendert",
                keep_all(vor, list(vor)) == list(vor))

        # And the case on which a SET instead of a counter would fail:
        # a label appears twice before and only once after.
        vor2 = [("A", "dup"), ("B", "dup")]
        nach2 = [("A", "dup")]
        her2 = keep_all(vor2, nach2)
        _urteil("ein von ZWEI gleichen Etiketten verlorener Frame wird bemerkt",
                len(her2) == 2 and any(x == "B" for x, _l in her2),
                "eine Menge statt eines Zaehlers saehe hier nichts")
    finally:
        if _alt_add is None:
            os.environ.pop(FLAG_ADD, None)
        else:
            os.environ[FLAG_ADD] = _alt_add

    print(f"\n  {'ALLE PROBEN BESTANDEN' if fehler == 0 else str(fehler) + ' PROBE(N) GESCHEITERT'}")
    return 1 if fehler else 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys
    _sys.exit(_rg_selbsttest())
