"""Welle-5l T6.2 — Theorem-D asymmetric-bidentate Δ/Λ enumeration.

Universal helicity classifier for **asymmetric** bidentate chelates.  Solves the
collapse problem of the achiral universal scalar-triple-product
(``_chirality_enumerator._classify_helicity``) on OH/SP/TBP-symmetric
tris-bidentate(A,B) complexes such as the user-flagged X10-YIVROM:

    Fe(III)·[(O,S)·(O,S)·(O,S)] — 3 asymmetric (O,S) chelates,
    Pólya expectation 4 distinct stereoisomers
    (fac-Δ, fac-Λ, mer-Δ, mer-Λ); pre-Iter-T6.2 emitted only 3.

The pre-existing universal classifier sorts chelate-pair vectors by
**donor-list index** which makes its sign meaningless when the chelate is
asymmetric: an (A,B) chelate and a (B,A) chelate would contribute opposite
signs even though they describe the same stereochemistry.  Theorem-D
reinstates a *chemically* meaningful orientation by sorting each chelate's
primary→secondary endpoint by the (atomic-number, donor-class-key)
**donor label** so that for an (O,S) chelate the vector ALWAYS runs
``O → S`` regardless of donor-list ordering.

Universal-fundamental doctrine (see ``feedback_universal_fundamental_doctrine``):

  - Trigger logic is **graph-feature**-driven (≥ 2 chelates with distinct
    donor labels) — no SMILES patterns, no element allowlists, no refcodes.
  - Works for ANY asymmetric (A,B) bidentate set: (N,O), (N,S), (O,S),
    (P,O), (P,N), (Cp-η5, X), (C-aryl, N), etc.
  - Δ/Λ is the right-hand rule on the chelate-twist direction, computed
    via the SO(3)-invariant mixed product
    ``Σ_{i<j} (c_i × c_j) · (p_j − p_i)`` where ``c_k = pos(B_k) − pos(A_k)``
    and ``p_k = pos(A_k)``.  A → B is the label-oriented primary→secondary
    direction; the second factor breaks the centroid-symmetry zero that
    sinks the legacy classifier on OH-tris-(A,B).
  - Bit-exact off:  caller must check
    ``DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE`` before invoking.

Stand-alone helper file — no smiles_converter imports, only the polyhedron
geometry vectors which are duplicated from ``_chirality_enumerator`` (kept
in lockstep with smiles_converter's ``_TOPO_GEOMETRY_VECTORS``).
"""
from __future__ import annotations

import os
from typing import FrozenSet, Iterable, List, Optional, Sequence, Tuple


def _delfin_env_int(name: str, default: int) -> int:
    """Local copy of ``smiles_converter._delfin_env_int`` (avoid import cycle)."""
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


# ---------------------------------------------------------------------------
# Trigger gate
# ---------------------------------------------------------------------------


def is_asymmetric_bidentate_set(
    chelate_pairs: Iterable[FrozenSet],
    donor_labels: Sequence[str],
) -> bool:
    """Return True when the chelate set warrants Theorem-D enumeration.

    Trigger conditions (ALL must hold):

      1. ≥ 2 bidentate chelates are present.
      2. ≥ 2 of those chelates have **distinct** donor labels
         (i.e. ``donor_labels[a] != donor_labels[b]``).

    This is graph-feature-driven only — no element allowlist, no SMILES
    pattern, no refcode lookup.  Works for any (A,B) chelate where the
    two donors differ in element OR donor-class key (Morgan environment).

    Symmetric chelates (en, oxalate, acac) → both donors carry the same
    label and contribute ``0`` to the asymmetric count.  A complex with
    only one asymmetric chelate cannot have helical handed-ness, so the
    gate also refuses single-asymmetric cases.
    """
    n_asym = 0
    n_total = 0
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        n_total += 1
        a, b = pair[0], pair[1]
        if a < 0 or b < 0 or a >= len(donor_labels) or b >= len(donor_labels):
            continue
        if donor_labels[a] != donor_labels[b]:
            n_asym += 1
    return n_total >= 2 and n_asym >= 2


# ---------------------------------------------------------------------------
# Polyhedron geometry vectors — duplicated from ``_chirality_enumerator``
# (which mirrors ``smiles_converter._TOPO_GEOMETRY_VECTORS``).  Kept here
# to remove the import-cycle risk.  Lockstep is enforced by tests.
# ---------------------------------------------------------------------------

_OH_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (2.0, 0.0, 0.0), (-2.0, 0.0, 0.0),
    (0.0, 2.0, 0.0), (0.0, -2.0, 0.0),
    (0.0, 0.0, 2.0), (0.0, 0.0, -2.0),
)
_TPR_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (2.0, 0.0, 1.0), (-1.0, 1.732, 1.0), (-1.0, -1.732, 1.0),
    (2.0, 0.0, -1.0), (-1.0, 1.732, -1.0), (-1.0, -1.732, -1.0),
)
_TBP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (0.0, 0.0, 2.0), (0.0, 0.0, -2.0),
    (2.0, 0.0, 0.0), (-1.0, 1.732, 0.0), (-1.0, -1.732, 0.0),
)
_SP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (0.0, 0.0, 2.2),
    (2.0, 0.0, 0.4), (0.0, 2.0, 0.4), (-2.0, 0.0, 0.4), (0.0, -2.0, 0.4),
)
_SAP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (1.414, 0.0, 0.8), (0.0, 1.414, 0.8), (-1.414, 0.0, 0.8), (0.0, -1.414, 0.8),
    (1.0, 1.0, -0.8), (-1.0, 1.0, -0.8), (-1.0, -1.0, -0.8), (1.0, -1.0, -0.8),
)
_DD_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (1.414, 0.0, 1.0), (0.0, 1.414, -1.0), (-1.414, 0.0, 1.0), (0.0, -1.414, -1.0),
    (0.9, 0.9, 0.0), (-0.9, 0.9, 0.0), (-0.9, -0.9, 0.0), (0.9, -0.9, 0.0),
)
_COH_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (2.0, 0.0, 0.0), (-2.0, 0.0, 0.0),
    (0.0, 2.0, 0.0), (0.0, -2.0, 0.0),
    (0.0, 0.0, 2.0), (0.0, 0.0, -2.0),
    (1.155, 1.155, 1.155),
)
_TTP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (1.633, 0.0, 1.155), (-0.816, 1.414, 1.155), (-0.816, -1.414, 1.155),
    (1.633, 0.0, -1.155), (-0.816, 1.414, -1.155), (-0.816, -1.414, -1.155),
    (1.0, 1.732, 0.0), (-2.0, 0.0, 0.0), (1.0, -1.732, 0.0),
)
_PBP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (0.0, 0.0, 2.0), (0.0, 0.0, -2.0),
    (2.0, 0.0, 0.0), (0.618, 1.902, 0.0),
    (-1.618, 1.176, 0.0), (-1.618, -1.176, 0.0), (0.618, -1.902, 0.0),
)
_TH_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (1.155, 1.155, 1.155),  (-1.155, -1.155, 1.155),
    (-1.155, 1.155, -1.155), (1.155, -1.155, -1.155),
)
_SS_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (0.0, 0.0, 2.0), (0.0, 0.0, -2.0),
    (2.0, 0.0, 0.4), (-2.0, 0.0, 0.4),
)


_GEOM_VECTORS = {
    'OH':  _OH_VECTORS,
    'TPR': _TPR_VECTORS,
    'TBP': _TBP_VECTORS,
    'SP':  _SP_VECTORS,
    'SAP': _SAP_VECTORS,
    'DD':  _DD_VECTORS,
    'COH': _COH_VECTORS,
    'TTP': _TTP_VECTORS,
    'PBP': _PBP_VECTORS,
    'TH':  _TH_VECTORS,
    'SS':  _SS_VECTORS,
}

_EPS = 1e-3


def _cross(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


# ---------------------------------------------------------------------------
# Theorem-D classifier
# ---------------------------------------------------------------------------


def classify_helicity_asymmetric(
    perm: Sequence[int],
    chelate_pairs: Iterable[FrozenSet],
    donor_labels: Sequence[str],
    geometry: str,
) -> str:
    """Return 'L' (Λ), 'D' (Δ), or '' (achiral / undefined).

    Theorem-D scalar (label-oriented mixed-product):

        For each asymmetric chelate ``k = (a_k, b_k)`` with
        ``donor_labels[a_k] < donor_labels[b_k]``:

            p_k  =  polyhedron vertex at perm-position of a_k   (PRIMARY)
            s_k  =  polyhedron vertex at perm-position of b_k   (SECONDARY)
            c_k  =  s_k − p_k                                   (chelate vector)

        Symmetric chelates (same donor label on both ends) are skipped — they
        do not contribute helicity to Theorem-D because their primary→secondary
        direction is undefined.

        Helicity scalar:

            H  =  Σ_{i<j}  (c_i × c_j) · (p_j − p_i)

        ``sign(H)`` is invariant under the polyhedron's point-group rotations
        and flips exactly once under any reflection — the signature of a
        helical pseudo-scalar.  The ``p_j − p_i`` factor (instead of
        ``m_i + m_j`` from the legacy classifier) breaks the centroid-symmetry
        zero that collapses the legacy helicity to 0 for OH-tris-(A,B).

    Returns:
        'L' if H > +EPS, 'D' if H < −EPS, '' otherwise (< 2 asymmetric
        chelates, unknown geometry, or numerical collapse to zero).
    """
    vectors = _GEOM_VECTORS.get(geometry)
    if vectors is None:
        return ''

    perm_list = list(perm)
    # Build label-oriented pair data: (c_vec, primary_pos) for each ASYM chelate.
    pair_data: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        a, b = pair[0], pair[1]
        if a < 0 or b < 0 or a >= len(donor_labels) or b >= len(donor_labels):
            continue
        la, lb = donor_labels[a], donor_labels[b]
        if la == lb:
            continue  # symmetric chelate skipped — Theorem-D scope
        # Primary = donor with smaller label (lexicographic).  This is the
        # universal-fundamental orientation rule: same primary→secondary
        # convention is applied for every chelate so the resulting
        # chelate-vector ``c_k`` is in a chemically meaningful direction.
        pri, sec = (a, b) if la < lb else (b, a)
        try:
            pos_p = perm_list.index(pri)
            pos_s = perm_list.index(sec)
        except ValueError:
            continue
        if pos_p >= len(vectors) or pos_s >= len(vectors):
            continue
        vp = vectors[pos_p]
        vs = vectors[pos_s]
        c = _sub(vs, vp)
        pair_data.append((c, vp))

    if len(pair_data) < 2:
        return ''

    # Canonical sort by primary-position vector (round-stable) so swapping the
    # chelate-iteration order does not flip the sign — required for permutation
    # symmetry invariance across equivalent perms.
    pair_data.sort(key=lambda cp: (round(cp[1][0], 3),
                                    round(cp[1][1], 3),
                                    round(cp[1][2], 3)))

    total = 0.0
    n = len(pair_data)
    for i in range(n):
        for j in range(i + 1, n):
            c_i, p_i = pair_data[i]
            c_j, p_j = pair_data[j]
            cross_cc = _cross(c_i, c_j)
            d_pp = _sub(p_j, p_i)
            total += _dot(cross_cc, d_pp)

    if total > _EPS:
        return 'L'
    if total < -_EPS:
        return 'D'
    return ''


def theorem_d_aware_pairs(
    canonical_pairs: tuple,
    perm: Optional[Sequence[int]],
    chelate_pairs: Optional[Iterable[FrozenSet]],
    donor_labels: Optional[Sequence[str]],
    geometry: str,
) -> tuple:
    """Append a Theorem-D helicity tag (``('chir_td', 'L'|'D')``) to a
    canonical-form tuple when the chelate set satisfies the asymmetric-
    bidentate gate AND the classifier returns a non-empty handed-ness.

    Returns ``canonical_pairs`` unchanged when:

      - any required argument is None,
      - the chelate set fails ``is_asymmetric_bidentate_set`` gate,
      - the classifier returns '' (numerical achirality / unknown geometry).

    The tag key is ``'chir_td'`` (not ``'chir'``) so it is **distinct** from
    the legacy Iter-2 helicity tag.  This lets the dedup pipeline use BOTH
    tags as independent split-keys without collisions when a permutation is
    simultaneously Λ-under-legacy AND Δ-under-Theorem-D (orthogonal axes).

    The caller is responsible for env-flag gating (the helper is pure-function
    no I/O).  Recommended gate:
    ``DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE``.
    """
    if perm is None or chelate_pairs is None or donor_labels is None:
        return canonical_pairs
    cp_list = list(chelate_pairs)
    if not is_asymmetric_bidentate_set(cp_list, donor_labels):
        return canonical_pairs
    h = classify_helicity_asymmetric(perm, cp_list, donor_labels, geometry)
    if not h:
        return canonical_pairs
    return canonical_pairs + (('chir_td', h),)


__all__ = [
    'is_asymmetric_bidentate_set',
    'classify_helicity_asymmetric',
    'theorem_d_aware_pairs',
]
