"""Orbit enumeration of topological coordination isomers (Burnside / Polya) in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

from typing import Dict, FrozenSet, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _all_polyhedra_codes,
    _class_conditional_flag,
    _delfin_env_int,
)
from delfin.manta.isomer_labels import (
    _TOPO_CANONICAL_FNS,
    _TOPO_GEOMETRY_VECTORS,
    _TOPO_TRANS_POSITIONS,
)
from delfin.manta.ml_tables import (
    _PREFERRED_CN4_GEOMETRY,
    _PREFERRED_CN5_GEOMETRY,
    _PREFERRED_CN6_GEOMETRY,
    _classify_cn5_geometry,
    _classify_cn5_geometry_from_labels,
    _cn5_d_electron_count,
)

logger = get_logger("delfin.smiles_converter")


def _enumerate_orbits_topo(
    geom_name: str,
    donor_labels: List[str],
    n_coord: int,
    chelate_pairs: List[FrozenSet],
    canonical_fn,
    trans_pos: List[Tuple[int, int]],
    chiral: bool,
    helicity_aware_pairs=None,
):
    """Burnside-orbit enumeration for one polyhedron (P5 — DELFIN_ORBIT_ENUM).

    Replaces the factorial ``itertools.permutations(range(n_coord))`` sweep
    with an enumeration over the *distinct chelate-extended label tuples*
    (``set(permutations(extended_labels))``), orbit-reduced by the proper
    rotation group of ``geom_name`` (CODE keys OH/SAP/DD/TTP/… from
    :mod:`delfin.manta._burnside_groups`, which cover CN<=12).  For a label
    multiset with multiplicities the number of distinct label tuples is
    ``n!/prod(mult!)`` — far smaller than ``n!`` — so high-CN cases
    (CN8 SAP/DD, CN9 TTP) enumerate completely without ever hitting the
    factorial cap that silently dropped isomers in the legacy loop.

    Returns a list of ``(orbit_key, canonical_form, perm_list)`` — ONE
    entry per proper-rotation orbit (the complete Cauchy-Frobenius count)
    — in deterministic order.  ``orbit_key`` is the lex-min orbit
    representative (used by the caller for cross-geometry dedup);
    ``perm_list`` is the lex-min donor-list-index assignment that realises
    the orbit representative's label arrangement.  Same-class donors are
    interchangeable, so the lex-min assignment is a canonical, reproducible
    choice.  The completeness basis is the proper rotation group (Λ/Δ
    enantiomer-distinct), matching ``fffree.polya_isomer_count.count_isomers``.

    Graph/group-theory only — no SMILES, refcode, or atom-index keying.
    Only invoked when ``DELFIN_ORBIT_ENUM=1``; the legacy factorial loop
    (with its shared ``perm_count`` and cap-return) is untouched when the
    flag is unset, preserving byte-identity (incl. the CN8 truncation).
    """
    import itertools as _it

    from delfin.manta._burnside_groups import get_groups as _get_groups

    proper, _full = _get_groups(geom_name)
    # Chelate-extended labels: append a per-chelate colour suffix to each
    # chelate donor so the orbit reduction is chelate-aware (mirrors
    # ``_build_extended_label`` in _prescribed_isomer_enumerator.py).
    extended: List[str] = list(donor_labels)
    for k, cp in enumerate(chelate_pairs):
        for d in cp:
            if 0 <= d < n_coord:
                extended[d] = f"{extended[d]}@c{k}"

    trans_set = frozenset(frozenset((a, b)) for a, b in trans_pos)

    # Map each distinct extended-label to the sorted list of donor-list
    # indices carrying it (used to recover a lex-min perm per arrangement).
    label_to_idxs: Dict[str, List[int]] = {}
    for di in range(n_coord):
        label_to_idxs.setdefault(extended[di], []).append(di)
    for lst in label_to_idxs.values():
        lst.sort()

    def _lexmin_perm(arrangement: Tuple[str, ...]) -> Optional[List[int]]:
        """Lex-smallest perm (positions->donor idx) realising ``arrangement``
        (extended[perm[pos]] == arrangement[pos]); deterministic."""
        used = [False] * n_coord
        perm: List[int] = []
        for pos in range(n_coord):
            want = arrangement[pos]
            chosen = -1
            for cand in label_to_idxs.get(want, ()):  # already sorted asc
                if not used[cand]:
                    chosen = cand
                    break
            if chosen < 0:
                return None
            used[chosen] = True
            perm.append(chosen)
        return perm

    # Return ONE (cf, perm) per Burnside orbit (no internal cf-dedup): the
    # caller applies the final dedup granularity (coarse ``cf`` when
    # DELFIN_BURNSIDE_FULL is off, fine ``(geom, bk, chir)`` when on),
    # exactly mirroring the legacy factorial loop's dedup semantics.
    results: List[Tuple[tuple, List[int]]] = []
    # Orbit-min over distinct extended-label arrangements.  Per-geometry
    # local scope: no shared counter, no cap-return across geometries.
    # Deterministic order: sort the distinct arrangements lexicographically.
    seen_orbit: set = set()
    for arrangement in sorted(set(_it.permutations(extended))):
        # Orbit-min key under the proper rotation group (chiral) — the
        # complete Burnside invariant of the label arrangement's orbit.
        best = arrangement
        for g in proper:
            if len(g) != n_coord:
                continue
            cand = tuple(arrangement[g[i]] for i in range(n_coord))
            if cand < best:
                best = cand
        if best in seen_orbit:
            continue
        seen_orbit.add(best)

        # Recover the lex-min concrete perm for the orbit representative
        # ``best`` so the chelate-trans validity filter and the canonical
        # form are evaluated on a real donor assignment.
        perm = _lexmin_perm(best)
        if perm is None:
            continue

        # Chelate-trans validity filter on the representative.
        valid = True
        for chelate in chelate_pairs:
            chelate_list = sorted(chelate)
            if len(chelate_list) != 2:
                continue
            try:
                pos_i = perm.index(chelate_list[0])
                pos_j = perm.index(chelate_list[1])
            except ValueError:
                continue
            if frozenset((pos_i, pos_j)) in trans_set:
                valid = False
                break
        if not valid:
            continue

        types = tuple(donor_labels[perm[pos]] for pos in range(n_coord))
        cf = canonical_fn(types)
        if chiral and helicity_aware_pairs is not None:
            try:
                cf = helicity_aware_pairs(cf, perm, chelate_pairs, geom_name)
            except Exception as _ce_exc:
                logger.debug(
                    "Chirality enumerator no-op (orbit, %s, perm=%s): %s",
                    geom_name, perm, _ce_exc,
                )
        # ===== THEOREM-D for ASYMMETRIC bidentates (wired 2026-08-06) =====
        # _theorem_d_asymmetric_bidentate.py (327 lines) had been sitting in the tree since
        # 18.05., written, documented and with its own env switch -- and was NEVER IMPORTED.
        # Zero import sites, no dynamic imports; one of eight such modules.
        #
        # WHAT IT SOLVES.  The existing universal classifier sorts the chelate vectors
        # by DONOR LIST INDEX.  For an asymmetric (A,B) chelate its sign is thereby
        # meaningless: an (A,B) and a (B,A) chelate yield opposite signs for the same
        # stereochemistry, the sum cancels, and the helicity collapses to ''.
        # Consequence, documented in the module on X10-YIVROM: Fe(III) with three
        # asymmetric (O,S) chelates has, by Polya, FOUR stereoisomers
        # (fac-Delta, fac-Lambda, mer-Delta, mer-Lambda) -- THREE get built.
        #
        # Theorem-D establishes a chemically meaningful orientation and attaches its OWN tag
        # ('chir_td'), not 'chir'.  That way both serve as independent split keys when a
        # permutation is at once Lambda-by-legacy and Delta-by-Theorem-D.
        #
        # BUILD FORM: purely ADDITIVE.  The wrapper returns cf unchanged when an argument is
        # missing, the chelate set does not pass the asymmetry gate, or the classifier yields
        # ''.  So it can only split a previously MERGED permutation, never remove an
        # existing one -- the same build form as the four flags that have each landed.
        # Not coupled to `chiral`: the legacy path collapses precisely on these cases,
        # the module's own gate (is_asymmetric_bidentate_set) is the right filter.
        # Default OFF -> byte-identical.
        if _delfin_env_int("DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE", 0):
            try:
                from delfin.manta._theorem_d_asymmetric_bidentate import (
                    theorem_d_aware_pairs as _td_pairs)
                cf = _td_pairs(cf, perm, chelate_pairs, donor_labels, geom_name)
            except Exception as _td_exc:
                logger.debug(
                    "Theorem-D no-op (orbit, %s, perm=%s): %s",
                    geom_name, perm, _td_exc,
                )
        # (orbit_key, cf, perm): orbit_key = the proper-rotation orbit
        # representative ``best``.  One entry per rotation orbit → the
        # complete Cauchy-Frobenius count (e.g. CN9 TTP N5O4 → 24).
        results.append((best, cf, list(perm)))

    # ===== REALISABILITY: what is chemically unbuildable is not enumerated in the first place =====
    # (wired 2026-08-06.  _realisability.py, 513 lines, in the tree since Welle-3, NEVER imported.)
    #
    # THE GOAL, in the user's words: all chemically realisable frames -- i.e. those that
    # can be built WITHOUT an anomaly.  Isomers and conformers that are buildable only with
    # chemically unrealistic defects we do not need.
    #
    # Polya/Burnside generates ALL orbit-distinct colourings under the point group.  Many
    # of them are combinatorially valid and chemically IMPOSSIBLE: a five-membered chelate ring
    # cannot span a 180-degree trans pair; two large sigma donors (P, As, Sb, wide-cone NHC)
    # do not fit onto adjacent vertices without getting below 2 x r_vdW; an annulated
    # aromatic donor cannot span two vertices whose angle deviates far from the planar
    # aryl ideal; d-electron count and donor sigma pairs can be incompatible.
    #
    # The builder tries them today anyway and delivers defective frames.  Those then count
    # doubly harmful: they drag every quality term (worst_sev reads the WORST
    # frame over the manifold) AND they inflate the denominator of the isomer coverage.  Leaving
    # them out here thus lowers the defect count and at the same time makes the completeness number honest.
    #
    # THE GATE FOR THIS IS ALREADY BUILT (loop.py:1541, user 2026-07-19): quality-weighted
    # completeness -- "all REALISTIC frames AT quality, ANCHORED by the CCDC isomer, NOT every
    # Polya isomer".  The generic isomers_lost floor is therefore soft; all CCDC-anchored
    # floors stay HARD.  Whoever removes an orbit here that the crystal actually is
    # fails at ccdc_isomer_lost -- exactly the right lock.
    #
    # Default OFF (DELFIN_REALISABILITY=0) -> bit-exact HEAD.
    if _delfin_env_int("DELFIN_REALISABILITY", 0) and results:
        try:
            from delfin.manta import _realisability as _realis
            _verts = _TOPO_GEOMETRY_VECTORS.get(geom_name)
            if _verts:
                _rep = _realis.filter_isomer_labels(
                    [(cf, perm) for _ok, cf, perm in results],
                    geom_name, _verts, donor_labels, n_coord, chelate_pairs)
                _kept_perms = {tuple(p) for _c, p in getattr(_rep, "kept", [])}
                _before = len(results)
                _filtered = [t for t in results if tuple(t[2]) in _kept_perms]
                if _filtered:                      # never wipe out the whole manifold
                    results = _filtered
                if len(results) != _before:
                    logger.debug("Realisability: %d von %d Orbits behalten (%s)",
                                 len(results), _before, geom_name)
        except Exception as _rl_exc:
            logger.debug("Realisability no-op (%s): %s", geom_name, _rl_exc)
    return results


def _enumerate_topological_isomers(
    donor_labels: List[str],
    n_coord: int,
    chelate_pairs: List[FrozenSet],
    metal_symbol: str = '',
    metal_formal_charge: int = 0,
    mol=None,
    metal_idx: Optional[int] = None,
    donor_indices: Optional[List[int]] = None,
) -> List[Tuple[tuple, List[int]]]:
    """Return all unique (canonical_form, permutation) pairs.

    ``perm[position] = donor_list_index`` — permutations where chelate-
    constrained donor pairs are never placed in trans positions.

    Args:
        donor_labels: Element symbol per donor atom (index = position in
            donor_indices list passed by the caller).
        n_coord: Coordination number (2–8).
        chelate_pairs: frozensets of donor-list indices that must stay cis.
        metal_symbol: Metal element symbol (used for CN=4/5 geometry
            preference).
        metal_formal_charge: Formal charge on the metal centre.  Optional,
            consumed only when ``DELFIN_CN5_GEOM_AWARE=1`` to compute the
            d-electron count for CN=5 polyhedron selection.  Defaults to 0
            (matches legacy uncharged-metal assumption).
        mol: Optional RDKit ``Mol`` for the rich CN=5 classifier.  When
            provided together with ``metal_idx`` + ``donor_indices`` and
            ``DELFIN_CN5_GEOM_AWARE=1``, the graph-based
            :func:`_classify_cn5_geometry` is consulted; otherwise the
            label-only facade is used.
        metal_idx: 0-based atom index of the metal centre in ``mol``.
        donor_indices: 0-based atom indices of the donor atoms in ``mol``
            in the same order as ``donor_labels``.

    Returns:
        List of (canonical_form_tuple, perm_list) for every unique isomer.
    """
    import itertools
    import os

    # Permutation-cap to prevent exponential explosion on high-CN + multi-chelate
    # cases. CN=8 gives 80640 perms, CN=9 gives 362880 — without cap, a single
    # call can run for minutes and consume 17-22 cores via numpy/MKL via RDKit
    # (observed: 2026-04-29 INSIGHTS_LOG ~18:51 UTC).  Default 50000 covers all
    # reasonable CN<=7 cases and most CN=8 cases; CN=9+ returns reduced set.
    max_perms_cap = int(os.environ.get('DELFIN_MAX_PERMS_CAP', '50000'))

    if n_coord == 2:
        geometries = ['LIN']
    elif n_coord == 3:
        geometries = _all_polyhedra_codes(3, metal_symbol, ['TP', 'TS'])
    elif n_coord == 4:
        # Metal-aware geometry ordering: preferred geometry first
        pref = _PREFERRED_CN4_GEOMETRY.get(metal_symbol, 'SQ')
        other = 'TH' if pref == 'SQ' else 'SQ'
        geometries = _all_polyhedra_codes(4, metal_symbol, [pref, other, 'SS'])
    elif n_coord == 5:
        # CN=5 polyhedron preference: chemistry-aware (DELFIN_CN5_GEOM_AWARE=1)
        # or legacy element-only map (default OFF -> bit-exact HEAD).
        if _delfin_env_int("DELFIN_CN5_GEOM_AWARE", 0):
            if (
                mol is not None
                and metal_idx is not None
                and donor_indices is not None
                and len(donor_indices) == 5
            ):
                pref5 = _classify_cn5_geometry(
                    mol, int(metal_idx), list(donor_indices)
                )
            else:
                pref5 = _classify_cn5_geometry_from_labels(
                    metal_symbol=metal_symbol,
                    formal_charge=int(metal_formal_charge or 0),
                    donor_labels=donor_labels,
                    chelate_pairs=chelate_pairs,
                )
        else:
            pref5 = _PREFERRED_CN5_GEOMETRY.get(metal_symbol, 'TBP')
        other5 = 'SP' if pref5 == 'TBP' else 'TBP'
        geometries = _all_polyhedra_codes(5, metal_symbol, [pref5, other5])
    elif n_coord == 6:
        pref6 = _PREFERRED_CN6_GEOMETRY.get(metal_symbol, 'OH')
        other6 = 'TPR' if pref6 == 'OH' else 'OH'
        _geoms6 = [pref6, other6]
        # REALISM (DELFIN_FFFREE_CN6_TPR_SUPPRESS, default-OFF -> byte-identical):
        # trigonal-prismatic CN6 is realistic ONLY for d0-d2 early-TM (dithiolene /
        # tris-S).  For d3-d10 the octahedron is overwhelmingly ligand-field-favoured
        # -> a TPR "isomer" is UNREALISTIC junk that twists rigid meridional ligands
        # out of the metal's pi-plane (the "skewed pi-systems" bloat) and never
        # matches the crystal.  Drop TPR for d3-d10 so the manifold holds only the
        # realistic (octahedral) isomers.  d0-d2 / unknown keep BOTH (never lose a
        # genuine prism -> the CCDC-isomer HARD floor stays safe).
        if _delfin_env_int("DELFIN_FFFREE_CN6_TPR_SUPPRESS", 0):
            _d6 = _cn5_d_electron_count(metal_symbol, int(metal_formal_charge or 0))
            if _d6 is not None and 3 <= _d6 <= 10:
                _geoms6 = [g for g in _geoms6 if g != 'TPR'] or ['OH']
        geometries = _all_polyhedra_codes(6, metal_symbol, _geoms6)
    elif n_coord == 7:
        geometries = _all_polyhedra_codes(7, metal_symbol, ['PBP', 'COH'])
    elif n_coord == 8:
        geometries = _all_polyhedra_codes(8, metal_symbol, ['SAP', 'DD'])
    elif n_coord == 9:
        geometries = ['TTP']
    elif n_coord == 10:
        # Iter-2D — DELFIN_CN_HIGH_ENABLE=1 enables CN10 lanthanide/actinide.
        # Default OFF (=0) → bit-exact HEAD (returns []).
        if _delfin_env_int("DELFIN_CN_HIGH_ENABLE", 0):
            geometries = ['BCSAP', 'PAP']
        else:
            return []
    elif n_coord == 11:
        if _delfin_env_int("DELFIN_CN_HIGH_ENABLE", 0):
            geometries = ['CPAP']
        else:
            return []
    elif n_coord == 12:
        if _delfin_env_int("DELFIN_CN_HIGH_ENABLE", 0):
            geometries = ['ICOS', 'CUBO', 'HBP']
        else:
            return []
    else:
        return []

    results: List[Tuple[tuple, List[int]]] = []
    seen_canonical: set = set()
    perm_count = 0

    # Iter-2 Λ/Δ helicity gate — hoisted out of the inner loop.  Active only
    # when DELFIN_CHIRAL_ENUM=1 and ≥ 2 chelate pairs (chirality undefined
    # otherwise).  Bit-exact when env-flag is off.
    _chir_enabled = (
        _delfin_env_int('DELFIN_CHIRAL_ENUM', 0)
        and len(chelate_pairs) >= 2
        # When the geometric ENANTIOMER_MIRROR post-pass is active it is the SOLE enantiomer source
        # (exact mirror -> dedup-able + coincidence-based elimination).  CHIRAL_ENUM here would build a
        # SECOND, INDEPENDENTLY-embedded Λ/Δ set (not exact mirrors) -> unmergeable duplicates (ATOZUG).
        # Disable it so the two mechanisms never collide.
        and not _delfin_env_int('DELFIN_FFFREE_ENANTIOMER_MIRROR', 0)
    )
    _helicity_aware_pairs = None
    if _chir_enabled:
        try:
            from delfin.manta._chirality_enumerator import helicity_aware_pairs as _helicity_aware_pairs  # type: ignore
        except Exception as _ce_exc:
            logger.debug("Chirality enumerator import failed: %s", _ce_exc)
            _chir_enabled = False

    # Iter-3 Pólya-Burnside completeness gate.  The achiral-by-construction
    # ``_canonical_<poly>`` heuristics collapse several orbit-distinct
    # stereoisomer permutations into a single bucket on PBP / SAP / DD /
    # COH / TPR / SS (see ``results/polya_audit.csv``).  When
    # ``DELFIN_BURNSIDE_FULL=1`` is set, an extra Burnside orbit-min key is
    # appended to the dedup tuple: this strictly *adds* enumeration buckets
    # without ever collapsing legitimate ones.
    #
    # Welle-5l T1.3 (formerly Welle-5j Agent C) — class-conditional wrap,
    # per-class default-OFF.  Default OFF → bit-exact HEAD; awaits the
    # ITER-adaptive_timeout_audit timeout-budget before flipping any
    # class default-ON.  Operator overrides:
    #   DELFIN_BURNSIDE_FULL=1               → enable globally
    #   DELFIN_BURNSIDE_FULL_CLASSES=sigma,hapto  → per-class opt-in
    # Welle-5i Agent G measured on 339 SMILES (gold7 + 332 class-balanced):
    # global default-flip FAILS (+28 frames, 7 catastrophic), but per-class:
    # sigma +203 named (77 imp / 9 reg), hapto +40 (4/0), multi-sigma −11
    # catastrophic, multi-hapto 0.  Re-arm class defaults via the operator
    # CLASSES override once timeout budget allows full enumeration.
    _burnside_enabled = _class_conditional_flag(
        "DELFIN_BURNSIDE_FULL", mol, default=0,
    )
    _burnside_canonical_key = None
    if _burnside_enabled:
        try:
            from delfin.manta._burnside_groups import burnside_canonical_key as _burnside_canonical_key  # type: ignore
        except Exception as _be_exc:
            logger.debug("Burnside enumerator import failed: %s", _be_exc)
            _burnside_enabled = False

    # ------------------------------------------------------------------
    # P5 — DELFIN_ORBIT_ENUM (default 0): Burnside-orbit enumeration that
    # replaces the factorial permutation sweep + shared ``perm_count`` +
    # whole-function cap-return.  The legacy loop initialises ``perm_count``
    # ONCE before the geometry loop and ``return results`` on cap-hit, so a
    # high-CN geometry that exhausts the cap (e.g. CN8 SAP at 40320) starves
    # every later geometry (DD) → isomers silently dropped (verified: CN9
    # N5O4 yielded 4 instead of the Burnside count 24).  The orbit path
    # enumerates the *distinct* chelate-extended label tuples per geometry
    # (n!/prod(mult!) << n!), orbit-reduces by the proper rotation group,
    # and resets the work per geometry — complete and cap-free for CN<=12.
    # Flag unset → fall through to the unchanged factorial loop below →
    # byte-identical (incl. the CN8 truncation).
    if _delfin_env_int("DELFIN_ORBIT_ENUM", 0):
        for geom_name in geometries:
            canonical_fn = _TOPO_CANONICAL_FNS[geom_name]
            trans_pos = _TOPO_TRANS_POSITIONS[geom_name]
            orbit_isomers = _enumerate_orbits_topo(
                geom_name=geom_name,
                donor_labels=donor_labels,
                n_coord=n_coord,
                chelate_pairs=chelate_pairs,
                canonical_fn=canonical_fn,
                trans_pos=trans_pos,
                chiral=_chir_enabled,
                helicity_aware_pairs=(
                    _helicity_aware_pairs if _chir_enabled else None
                ),
            )
            for orbit_key, cf, perm in orbit_isomers:
                # The orbit path is inherently Burnside-complete (one entry
                # per proper-rotation orbit), so it does NOT route through
                # the coarse ``cf``-dedup or the additive DELFIN_BURNSIDE_FULL
                # key — both of which exist to recover orbits the factorial
                # loop would otherwise collapse.  Cross-geometry dedup uses
                # ``(geom_name, orbit_key)`` so equal cf across distinct
                # polyhedra never collide and never collapse.
                dedup_key = (geom_name, orbit_key)
                if dedup_key not in seen_canonical:
                    seen_canonical.add(dedup_key)
                    results.append((cf, list(perm)))
        return results

    for geom_name in geometries:
        canonical_fn = _TOPO_CANONICAL_FNS[geom_name]
        trans_pos = _TOPO_TRANS_POSITIONS[geom_name]

        for perm in itertools.permutations(range(n_coord)):
            perm_count += 1
            if perm_count > max_perms_cap:
                # Early exit: return whatever isomers we found so far rather
                # than blocking subprocess for minutes on exhaustive enumeration.
                logger.warning(
                    "_enumerate_topological_isomers cap %d reached (CN=%d, "
                    "n_chelates=%d, geom=%s); returning %d isomers",
                    max_perms_cap, n_coord, len(chelate_pairs),
                    geom_name, len(results),
                )
                return results

            # Validate chelate constraints: paired donors must not sit trans
            valid = True
            for chelate in chelate_pairs:
                chelate_list = sorted(chelate)
                # Find the geometry positions of these two donors.
                # chelate_list contains donor-list indices (0..n_coord-1),
                # which are always present in perm (a permutation of the
                # same range), so .index() will always succeed.
                pos_i = perm.index(chelate_list[0])
                pos_j = perm.index(chelate_list[1])
                for ta, tb in trans_pos:
                    if (pos_i == ta and pos_j == tb) or (pos_i == tb and pos_j == ta):
                        valid = False
                        break
                if not valid:
                    break
            if not valid:
                continue

            # Build canonical form from donor element symbols at each position
            types = tuple(donor_labels[perm[pos]] for pos in range(n_coord))
            cf = canonical_fn(types)
            # Iter-2 Λ/Δ helical-isomer enumeration (env-gated, additive):
            # When DELFIN_CHIRAL_ENUM=1 and ≥ 2 chelate pairs are present,
            # tag the canonical form with helicity ('L'/'D') so Λ and Δ
            # enantiomeric permutations survive as separate buckets instead
            # of collapsing via the achiral-by-construction _canonical_*.
            # When env-flag is off (default) cf is unchanged → bit-exact HEAD.
            if _chir_enabled and _helicity_aware_pairs is not None:
                try:
                    cf = _helicity_aware_pairs(
                        cf, perm, chelate_pairs, geom_name,
                    )
                except Exception as _ce_exc:
                    logger.debug(
                        "Chirality enumerator no-op (%s, perm=%s): %s",
                        geom_name, perm, _ce_exc,
                    )
            # Iter-3 Pólya-Burnside dedup key (env-gated; bit-exact when off).
            # Verified against 105-row `results/polya_audit.csv`: the
            # Burnside orbit-min `types`-tuple is a *complete* invariant
            # of the donor-multiset orbit under the polyhedron's point
            # group.  We therefore use it ALONE as dedup key — using
            # ``(cf, bk)`` would over-count whenever ``_canonical_<poly>``
            # is finer than the true point-group orbit, e.g. for SAP D4d
            # where the trans-pair-sort heuristic is finer than Pólya.
            # Geometry tag is prefixed so PBP-bucket keys never collide
            # with COH-bucket keys for the same multiset.
            if _burnside_enabled and _burnside_canonical_key is not None:
                try:
                    bk = _burnside_canonical_key(
                        geom_name, types, chiral=_chir_enabled,
                    )
                    # When chirality enumeration is on, augment with the
                    # helicity tag (already baked into cf via Iter-2 helper)
                    # so Λ/Δ partners stay separate even though their
                    # achiral-Burnside bk collides.
                    chir_tag = ''
                    if _chir_enabled:
                        for _itm in cf:
                            if (isinstance(_itm, tuple) and len(_itm) == 2
                                    and _itm[0] == 'chir'):
                                chir_tag = _itm[1]
                                break
                    dedup_key = (geom_name, bk, chir_tag)
                except Exception as _be_exc:
                    logger.debug(
                        "Burnside enumerator no-op (%s, perm=%s): %s",
                        geom_name, perm, _be_exc,
                    )
                    dedup_key = cf
            else:
                dedup_key = cf
            if dedup_key not in seen_canonical:
                seen_canonical.add(dedup_key)
                results.append((cf, list(perm)))

    return results
