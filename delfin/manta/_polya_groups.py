"""Polyhedron symmetry-group registry — single source-of-truth for
prescribed-isomer enumeration (Welle-5n-Pre, 2026-05-18).

This module is the **enumeration-layer** facade over the existing
:mod:`delfin.manta._burnside_groups` symmetry-builder.  It adds the
graph-/CN-driven *selection logic* and the *positional-descriptor*
helper needed by :mod:`delfin.manta._prescribed_isomer_enumerator` to label
each orbit representative with a chemistry-faithful tag (``fac`` /
``mer`` / ``cis`` / ``trans`` / ``ccc`` / ``axial`` …) — all derived
from the polyhedron geometry vectors, never from SMILES patterns.

Doctrine (per ``feedback_universal_fundamental_doctrine`` 2026-05-18):
the entire dispatch is keyed on (coordination number, graph features),
so **no element allowlist, no SMILES regex, no refcode prefix** appears.
Every coordination centre — d-block, f-block, main-group, sigma, hapto
— goes through the same selection.

Public API
----------
``polyhedra_for_cn(cn)``
    Return the canonical list of polyhedron codes (as registered in
    :mod:`delfin.manta._burnside_groups._GEO`) for a given coordination
    number.  Order = chemistry-preferred first (e.g. ``["OH", "TPR"]``
    for CN6).
``trans_positions(geom)``
    Return the canonical list of ``(i, j)`` index pairs at 180° in the
    polyhedron — used for chelate-cis enforcement during enumeration
    and for cis/trans positional-descriptor computation.
``positional_descriptor(geom, types, perm, chelate_pairs)``
    Return a short human-readable positional tag for a single orbit
    representative (e.g. ``"fac"``, ``"mer"``, ``"cis"``, ``"trans"``,
    ``"axial"``, ``"all-cis"``, ``"ccc"``, or ``""`` for trivial /
    single-orbit cases).  Pure geometry/graph signal.
``polyhedron_geometry(geom)``
    Tuple of vertex coordinates (mirror of ``_burnside_groups._GEO``).
"""

from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple

from delfin.manta._burnside_groups import (
    _GEO as _BURNSIDE_GEO,
    burnside_canonical_key,
    get_groups,
    group_size,
)


# ---------------------------------------------------------------------------
# CN → polyhedra preference list.
#
# Chemistry-preferred first; alternatives second.  The selection matches
# the SMILES-converter's ``_enumerate_topological_isomers`` dispatch so
# that the *prescribed*-layer agrees on what is reachable per CN.  CN10-
# 12 (lanthanide/actinide) polyhedra are returned unconditionally — the
# prescribed layer never gates by element block.  Caller can subset via
# the dispatcher (e.g. f-block-only for the realisation layer).
# ---------------------------------------------------------------------------
_CN_TO_POLYHEDRA: Dict[int, Tuple[str, ...]] = {
    2:  ("LIN",),
    3:  ("TP", "TS"),
    4:  ("TH", "SQ", "SS"),
    5:  ("TBP", "SP"),
    6:  ("OH", "TPR"),
    7:  ("PBP", "COH"),
    8:  ("SAP", "DD"),
    9:  ("TTP",),
    10: ("BCSAP", "PAP"),
    11: ("CPAP",),
    12: ("ICOS", "CUBO", "HBP"),
}


# ---------------------------------------------------------------------------
# Trans-position tables — mirror of ``_TOPO_TRANS_POSITIONS`` in
# ``smiles_converter.py``.  Duplicated here intentionally so that the
# enumeration-layer is import-cycle-free of the heavy converter; both
# tables MUST stay in lockstep (audited by tests).
# ---------------------------------------------------------------------------
_TRANS_POSITIONS: Dict[str, Tuple[Tuple[int, int], ...]] = {
    'LIN':   ((0, 1),),
    'TP':    (),
    'TS':    ((0, 1),),
    'OH':    ((0, 1), (2, 3), (4, 5)),
    'SQ':    ((0, 1), (2, 3)),
    'TH':    (),
    'TBP':   ((0, 1),),
    'SP':    ((1, 3), (2, 4)),
    'PBP':   ((0, 1),),
    'SAP':   ((0, 6), (1, 7), (2, 4), (3, 5)),
    'DD':    ((0, 2), (1, 3), (4, 6), (5, 7)),
    'SS':    ((0, 1),),
    'TPR':   (),
    'COH':   ((0, 1), (2, 3), (4, 5)),
    'TTP':   (),
    'BCSAP': ((0, 2), (1, 3), (4, 6), (5, 7), (8, 9)),
    'PAP':   ((0, 7), (1, 8), (2, 9), (3, 5), (4, 6)),
    'CPAP':  (),
    'ICOS':  ((0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)),
    'CUBO':  ((0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11)),
    'HBP':   ((0, 3), (1, 4), (2, 5), (6, 9), (7, 10), (8, 11)),
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def polyhedra_for_cn(cn: int) -> Tuple[str, ...]:
    """Return canonical polyhedra for a coordination number.

    Parameters
    ----------
    cn : int
        Observed coordination number (2..12 covered; outside that range
        returns an empty tuple).

    Returns
    -------
    Tuple of polyhedron codes (e.g. ``("OH", "TPR")`` for CN=6).
    """
    return _CN_TO_POLYHEDRA.get(int(cn), ())


def trans_positions(geom: str) -> Tuple[Tuple[int, int], ...]:
    """Return the (a, b) pairs of polyhedron vertex indices at 180°.

    Empty tuple for polyhedra without strict trans pairs (TH, TP, TPR,
    TTP, CPAP).  Used both for chelate-cis enforcement (chelate may not
    occupy a trans pair) and for cis/trans descriptor generation.
    """
    return _TRANS_POSITIONS.get(geom, ())


def polyhedron_geometry(geom: str) -> Tuple[Tuple[float, float, float], ...]:
    """Return the polyhedron geometry vectors (same as ``_burnside_groups._GEO``)."""
    return _BURNSIDE_GEO.get(geom, ())


# ---------------------------------------------------------------------------
# Positional descriptor — pure geometry signal.
#
# For each polyhedron we recognise a small set of *named* positional
# isomers determined by the donor-class multiset and the trans-pair
# graph.  This is the place where ``fac`` / ``mer`` / ``cis`` / ``trans``
# / ``axial`` etc. labels originate.  The logic is graph/geometry-only:
#
#  * For a multiset with k distinct classes, we walk every donor class
#    and ask: "does this class form an all-cis triangle (= fac) or a
#    meridional row (= mer)?"  Both predicates are CN-agnostic.
#  * Two-class A2B2 / A3B3 / A4B2 etc.: cis / trans / fac / mer fall out
#    of the same trans-pair / triangle analysis.
#  * A single-element multiset / a singleton orbit collapse to "" (no
#    further qualifier needed beyond the polyhedron + multiset key).
# ---------------------------------------------------------------------------
def _positions_of(types: Sequence[str], cls: str) -> List[int]:
    """Indices of donor positions whose donor-class equals ``cls``."""
    return [i for i, t in enumerate(types) if t == cls]


def _all_trans(positions: Sequence[int], trans_pairs: Sequence[Tuple[int, int]]) -> bool:
    """True iff every pair in ``positions`` is a trans-pair of the geom."""
    pos = list(positions)
    if len(pos) % 2 != 0:
        return False
    if not pos:
        return False
    trans_set = {frozenset(p) for p in trans_pairs}
    matched = set()
    for i in range(len(pos)):
        if i in matched:
            continue
        partnered = False
        for j in range(i + 1, len(pos)):
            if j in matched:
                continue
            if frozenset((pos[i], pos[j])) in trans_set:
                matched.add(i)
                matched.add(j)
                partnered = True
                break
        if not partnered:
            return False
    return len(matched) == len(pos)


def _any_trans(positions: Sequence[int], trans_pairs: Sequence[Tuple[int, int]]) -> bool:
    """True iff ANY pair within ``positions`` is at 180°."""
    if len(positions) < 2:
        return False
    s = set(positions)
    for a, b in trans_pairs:
        if a in s and b in s:
            return True
    return False


def _is_fac(positions: Sequence[int], trans_pairs: Sequence[Tuple[int, int]]) -> bool:
    """A class of 3 positions is "fac" iff none of the 3 pairs are trans."""
    if len(positions) != 3:
        return False
    if not _any_trans(positions, trans_pairs):
        return True
    return False


def _is_mer(positions: Sequence[int], trans_pairs: Sequence[Tuple[int, int]]) -> bool:
    """A class of 3 positions is "mer" iff EXACTLY ONE of the 3 pairs is trans.

    On an octahedron, 3-class A3 mer = one trans + two cis (≡ T-shape on
    the polyhedron).  Other CNs without strict trans pairs cannot
    register mer.
    """
    if len(positions) != 3:
        return False
    s = set(positions)
    n_trans = 0
    for a, b in trans_pairs:
        if a in s and b in s:
            n_trans += 1
    return n_trans == 1


def _chelate_position_signature(
    perm: Sequence[int],
    chelate_pairs: Iterable[FrozenSet[int]],
    trans_pairs: Sequence[Tuple[int, int]],
) -> Tuple[int, int]:
    """Return (n_chelate_pairs, n_chelate_at_trans_positions).

    Helper for chelate-aware positional descriptors.  Even when the
    multiset is homoleptic (single donor class), chelate-pair placement
    distinguishes physically distinct isomers.  For tris-bidentate on
    OH, all chelates sit on cis pairs → distinguishing Δ/Λ propellers.
    """
    n_total = 0
    n_trans = 0
    trans_set = {frozenset(p) for p in trans_pairs}
    perm_list = list(perm)
    for cp in chelate_pairs:
        donors = sorted(cp)
        if len(donors) != 2:
            continue
        try:
            pa = perm_list.index(donors[0])
            pb = perm_list.index(donors[1])
        except ValueError:
            continue
        n_total += 1
        if frozenset((pa, pb)) in trans_set:
            n_trans += 1
    return n_total, n_trans


def positional_descriptor(
    geom: str,
    types: Sequence[str],
    perm: Sequence[int],
    chelate_pairs: Iterable[FrozenSet[int]] = (),
) -> str:
    """Return a positional tag for one orbit representative.

    Inputs
    ------
    geom : str
        Polyhedron code (``"OH"`` …).
    types : sequence of donor-class strings, length = CN
        ``types[i]`` is the donor class at polyhedron vertex i under
        this orbit representative.
    perm : permutation
        Required by the signature for parallel use with the chirality
        classifier (``perm[pos] = donor_list_index``); the
        positional-descriptor itself reads only the resulting ``types``.
    chelate_pairs : iterable of frozenset(donor_list_index pairs)
        Currently unused in positional logic (chelates affect chirality
        but the positional tag is derived from ``types`` alone).

    Tag grammar
    -----------
    Empty string ``""`` means "no further qualifier needed beyond the
    polyhedron + multiset key".  Non-empty tags follow this hierarchy:

    *  ``"only-isomer"`` — homoleptic (single donor class).
    *  ``"trans"``      — A2B(CN-2) with the A-pair across a trans
                          line; or A2B2 / A2B2 etc. with all A-A and
                          B-B at 180°.
    *  ``"cis"``        — A2... with the A-pair adjacent (default
                          opposite of ``"trans"`` for 2-of-N classes).
    *  ``"all-cis"``    — every multi-class pair is at 90° (e.g.
                          octahedral A2B2C2 with no trans pair).
    *  ``"fac"``        — a 3-class occupies an all-cis triangle.
    *  ``"mer"``        — a 3-class spans one trans pair (the meridian).
    *  ``"axial"``      — TBP / PBP class sits on both axial vertices.
    *  ``"equatorial"`` — TBP / PBP class sits on equatorial only.
    *  ``"ccc-trans"``  — A3B3 / A2B2C2 ``trans-trans-trans`` (all
                          three classes occupy trans pairs).
    *  empty            — orbit alone determines the isomer fully
                          (multiset key already canonical).
    """
    tps = _TRANS_POSITIONS.get(geom, ())
    type_list = list(types)
    counts: Dict[str, int] = {}
    for t in type_list:
        counts[t] = counts.get(t, 0) + 1

    # ---- Homoleptic chelate-aware path ---------------------------------
    chelate_list = [cp for cp in chelate_pairs if len(cp) == 2]
    if len(counts) <= 1:
        if len(chelate_list) >= 2:
            # Multi-chelate homoleptic — Λ/Δ helicity distinguishes orbits.
            # Tag by chelate-trans count: ``tris-chelate`` for all-cis-pair
            # tris-bidentate (OH tris-en classic), ``chelate-mixed`` otherwise.
            n_total, n_trans = _chelate_position_signature(
                perm, chelate_list, tps
            )
            if n_total >= 2 and n_trans == 0:
                return "chelate-propeller"
            return "chelate-mixed"
        return "only-isomer"

    # ---- Bis-chelate path (2 chelates) -- distinguish trans/cis arrangement
    if len(chelate_list) == 2 and tps:
        n_total, n_trans = _chelate_position_signature(perm, chelate_list, tps)
        # Both chelates are bidentate (length-2 chelate pairs).  At each
        # pair, the donors must be cis (we filter that earlier in the
        # enumerator).  The remaining freedom is how the chelates relate
        # to each other: parallel "side-by-side" vs "perpendicular".
        # Detect via: do the chelate position-sets share a common axis?
        cp_sets = [set(c) for c in chelate_list]
        # If chelates share a vertex (impossible if donors disjoint) → skip.
        # Otherwise: check if the two chelates lie on a common "plane"
        # of the polyhedron (= same Burnside-equivalence class of position
        # pairs).  For OH bis-chelate: two cis pairs in same plane = cis;
        # two cis pairs in orthogonal planes = "trans" (chelates opposite).
        # We just return "" here and let the multiset_key disambiguate;
        # the bis-chelate descriptor is best derived from chelate-position
        # tuples which are already on the entry.

    # ---- Octahedral-family (OH, COH, SAP, DD, BCSAP, ICOS, CUBO, HBP) ----
    # Recognise mer/fac for 3-occupancy classes; trans/cis for 2-occ.
    if tps:
        # Pre-scan class positions.
        by_cls = {c: _positions_of(type_list, c) for c in counts}

        # All three-classes that fit fac/mer?
        three_classes = [c for c, n in counts.items() if n == 3]
        for c in three_classes:
            if _is_mer(by_cls[c], tps):
                return "mer"
            if _is_fac(by_cls[c], tps):
                return "fac"

        # All two-classes occupying trans pairs → trans-iso descriptor
        two_classes = [c for c, n in counts.items() if n == 2]
        if two_classes:
            # If at least one 2-class is fully trans → label it "trans"
            for c in two_classes:
                if _all_trans(by_cls[c], tps):
                    # Multi-2-class: ccc-trans if EVERY 2-class is trans
                    if all(_all_trans(by_cls[t], tps) for t in two_classes):
                        if len(two_classes) >= 3:
                            return "ccc-trans"
                        return "trans"
                    return "trans"
            # No 2-class is fully trans → cis (default for A2B2-style)
            # All-cis catch (A2B2C2 OH all-cis = special name)
            n_total = sum(counts.values())
            if (
                n_total == 6
                and len(two_classes) == 3
                and all(len(by_cls[c]) == 2 and not _all_trans(by_cls[c], tps)
                        for c in two_classes)
            ):
                return "all-cis"
            return "cis"

        # 4-occ class on SQ-like: trans = both pairs collinear, cis = adjacent
        four_classes = [c for c, n in counts.items() if n == 4]
        for c in four_classes:
            if _all_trans(by_cls[c], tps):
                return "trans"
            if _any_trans(by_cls[c], tps):
                # mixed: at least one trans pair, but not all → mixed cis/trans
                return "cis"

    # ---- TBP / PBP — axial / equatorial split -------------------------
    if geom in ("TBP", "PBP"):
        axial = {0, 1}
        for c, n in counts.items():
            if n == 2:
                positions = set(_positions_of(type_list, c))
                if positions == axial:
                    return "axial"
                if positions.isdisjoint(axial):
                    return "equatorial"
            if n == len(type_list) - 2:
                # one class fills the equator, partner sits on axes
                positions = set(_positions_of(type_list, c))
                if positions.issuperset({i for i in range(len(type_list)) if i not in axial}):
                    return "equatorial-major"
        return ""

    # ---- TH / TP / TS / TPR / TTP / CPAP — no trans constraints -------
    # Tetrahedral / trigonal-planar: every distinct multiset = one isomer.
    return ""


__all__ = [
    "polyhedra_for_cn",
    "trans_positions",
    "polyhedron_geometry",
    "positional_descriptor",
    "burnside_canonical_key",
    "get_groups",
    "group_size",
]
