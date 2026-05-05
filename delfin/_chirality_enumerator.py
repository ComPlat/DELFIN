"""Iter-2 Λ/Δ Helical-Isomer Enumeration Helper.

Computes the helicity (Λ / Δ / —) of a coordination-polyhedron permutation
given its chelate-pair graph.  The doctrine matches Iter-12 / Iter-14:

    - Universal: works for any polyhedron with ≥ 2 chelate pairs (OH, TPR,
      TBP, SP, SAP, DD, COH).  Driven by the polyhedron geometry vectors
      and the chelate-pair list — no element / SMILES / refcode shortcuts.
    - Geometry-driven: helicity is the sign of a scalar triple-product /
      pseudo-vector built from chelate-spanning vectors. No hardcoded
      angle thresholds beyond a small ``EPS`` floor for achiral cases.
    - Env-flag-gated: callers must check ``DELFIN_CHIRAL_ENUM`` themselves;
      this helper only computes the classification given the inputs.
    - Side-effect-free: returns 'L', 'D', or '' — no I/O, no rng, no globals.

Stufe-3-Innovation pattern (BECHUQ Λ/Δ gap, 2026-05-04):

    Multi-chelate Polyhedra have first-class chirality (Λ vs Δ helical
    enantiomers) that the achiral-by-construction `_canonical_*` (double
    `sorted()`) destroys.  This helper restores the missing degree of
    freedom by computing a deterministic helicity tag from the donor
    permutation and the chelate-pair list, which can then be appended
    to the canonical-form tuple to keep both hands as separate buckets.
"""
from __future__ import annotations

from typing import FrozenSet, Iterable, List, Optional, Sequence, Tuple

# Polyhedron geometry vectors (must match _TOPO_GEOMETRY_VECTORS in
# smiles_converter.py).  Duplicated here intentionally to keep the helper
# import-cycle-free; the tuple of tuples is a tiny constant.
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
# Tricapped trigonal prism (D3h, CN=9): 6 prism vertices (top z=+1.155,
# bottom z=-1.155, eclipsed) + 3 equatorial caps at z=0 between prism edges.
# Identical to _TOPO_GEOMETRY_VECTORS['TTP'] in smiles_converter.py - kept in
# lockstep so chirality classifier and isomer enumerator agree on positions.
_TTP_VECTORS: Tuple[Tuple[float, float, float], ...] = (
    (1.633, 0.0, 1.155), (-0.816, 1.414, 1.155), (-0.816, -1.414, 1.155),
    (1.633, 0.0, -1.155), (-0.816, 1.414, -1.155), (-0.816, -1.414, -1.155),
    (1.0, 1.732, 0.0), (-2.0, 0.0, 0.0), (1.0, -1.732, 0.0),
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


def _norm2(a: Sequence[float]) -> float:
    return _dot(a, a)


# ---------------------------------------------------------------------------
# Iter-2C polyhedron-specific Λ/Δ classifiers (DD, SAP, TPR)
#
# The universal scalar-triple-product (used for OH/TBP/SP/COH below) collapses
# to zero for several D2d/D4d/D3h-symmetric chelate-pair configurations, so
# DD CN8, SAP CN8, TPR CN6 receive dedicated classifiers driven by the
# *signed projection of the chelate-pair cross-product onto the principal
# improper-rotation axis* of the polyhedron (z by construction in our
# canonical vector-set).  No SMILES-, refcode-, or element-specific logic.
# ---------------------------------------------------------------------------


def _signed_angle_around_z(va: Sequence[float], vb: Sequence[float]) -> float:
    """Signed planar angle from va->vb measured CCW around the +z axis.

    Returns a number in (-pi, +pi].  Used as a tier-internal handedness
    discriminator for SAP and TPR where the principal axis is z.
    """
    import math
    ang_a = math.atan2(va[1], va[0])
    ang_b = math.atan2(vb[1], vb[0])
    d = ang_b - ang_a
    while d > math.pi:
        d -= 2 * math.pi
    while d < -math.pi:
        d += 2 * math.pi
    return d


def _helicity_dd(perm: Sequence[int],
                 chelate_pairs: Iterable[FrozenSet]) -> str:
    """DD (D2d, CN=8) - interpenetrating-trapezoid helicity.

    Vertex layout (`_DD_VECTORS`):
      Trap-A:   0,2 at (+/-x,0,+1); 1,3 at (0,+/-y,-1)   - alternating tier
      Trap-B:   4-7 at (+/-0.9,+/-0.9, 0)                - middle ring

    S4 axis along z.  Helicity scalar combines (1) the z-component of the
    chelate-spanning cross-product and (2) a tier-tilt term distinguishing
    Trap-A vs Trap-B (vertices 0..3 vs 4..7).
    """
    vectors = _GEOM_VECTORS['DD']
    perm_list = list(perm)
    total = 0.0
    pair_count = 0
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        try:
            pos_a = perm_list.index(pair[0])
            pos_b = perm_list.index(pair[1])
        except ValueError:
            continue
        if pos_a >= 8 or pos_b >= 8:
            continue
        va = vectors[pos_a]
        vb = vectors[pos_b]
        # z-projected angular-momentum term
        cross_z = va[0] * vb[1] - va[1] * vb[0]
        # tier-tilt term: distinguishes Trap-A (pos<4) vs Trap-B (pos>=4)
        tier_a = 1 if pos_a < 4 else -1
        tier_b = 1 if pos_b < 4 else -1
        tier_tilt = 0.5 * (tier_a - tier_b) * (va[2] + vb[2])
        total += cross_z + tier_tilt
        pair_count += 1
    if pair_count < 2:
        return ''
    if total > _EPS:
        return 'L'
    if total < -_EPS:
        return 'D'
    return ''


def _helicity_sap(perm: Sequence[int],
                  chelate_pairs: Iterable[FrozenSet]) -> str:
    """SAP (D4d, CN=8) - square-antiprismatic helicity.

    Vertex layout (`_SAP_VECTORS`):
      Top square z=+0.8: 0..3 at angles {0, 90, 180, 270} deg
      Bot square z=-0.8: 4..7 at angles {45, 135, 225, 315} deg

    S8 axis along z. Helicity scalar = z-cross of cross-tier chelate pairs
    + signed-angle of within-tier chelate pairs.
    """
    vectors = _GEOM_VECTORS['SAP']
    perm_list = list(perm)
    total = 0.0
    pair_count = 0
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        try:
            pos_a = perm_list.index(pair[0])
            pos_b = perm_list.index(pair[1])
        except ValueError:
            continue
        if pos_a >= 8 or pos_b >= 8:
            continue
        va = vectors[pos_a]
        vb = vectors[pos_b]
        cross_z = va[0] * vb[1] - va[1] * vb[0]
        tier_a = 0 if pos_a < 4 else 1
        tier_b = 0 if pos_b < 4 else 1
        if tier_a != tier_b:
            total += cross_z
        else:
            total += 0.5 * _signed_angle_around_z(va, vb)
        pair_count += 1
    if pair_count < 2:
        return ''
    if total > _EPS:
        return 'L'
    if total < -_EPS:
        return 'D'
    return ''


def _helicity_tpr(perm: Sequence[int],
                  chelate_pairs: Iterable[FrozenSet]) -> str:
    """TPR (D3h, CN=6) - eclipsed trigonal-prism helicity.

    Vertex layout (`_TPR_VECTORS`):
      Top triangle z=+1: 0,1,2 at angles {0, 120, 240} deg
      Bot triangle z=-1: 3,4,5 at angles {0, 120, 240} deg (eclipsed)
    """
    vectors = _GEOM_VECTORS['TPR']
    perm_list = list(perm)
    total = 0.0
    pair_count = 0
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        try:
            pos_a = perm_list.index(pair[0])
            pos_b = perm_list.index(pair[1])
        except ValueError:
            continue
        if pos_a >= 6 or pos_b >= 6:
            continue
        va = vectors[pos_a]
        vb = vectors[pos_b]
        cross_z = va[0] * vb[1] - va[1] * vb[0]
        tier_a = 0 if pos_a < 3 else 1
        tier_b = 0 if pos_b < 3 else 1
        if tier_a != tier_b:
            total += cross_z
        else:
            total += 0.7 * _signed_angle_around_z(va, vb)
        pair_count += 1
    if pair_count < 2:
        return ''
    if total > _EPS:
        return 'L'
    if total < -_EPS:
        return 'D'
    return ''


_GEOM_HELICITY_FN = {
    'DD':  _helicity_dd,
    'SAP': _helicity_sap,
    'TPR': _helicity_tpr,
}


def _classify_helicity(
    perm: Sequence[int],
    chelate_pairs: Iterable[FrozenSet],
    geometry: str,
) -> str:
    """Return 'L' (Λ), 'D' (Δ), or '' (achiral / undefined).

    Universal scalar-triple-product method for ≥ 2 chelate pairs:

      For each chelate pair (a, b) in canonical (sorted) donor-list-index
      order, the chelate-spanning vector at the polyhedron is

          c_k = pos(perm.index(b)) - pos(perm.index(a))

      where ``pos(i)`` is the i-th polyhedron geometry vector.  The chelate
      *centroid-vector* is

          m_k = pos(perm.index(b)) + pos(perm.index(a))   (×0.5 implicit)

      Helicity-pseudovector:

          P = Σ_k  c_k × m_k_unit

      where ``m_k_unit`` is the centroid-vector normalised.  ``P`` points
      along the helical axis.  The handed-ness comes from the sign of the
      mixed product

          sign(  (c_1 × c_2) · (m_1 + m_2)  )

      which is invariant under the relabelling symmetries of the polyhedron
      (mirroring through any plane flips the sign exactly once).  For ≥ 3
      chelate pairs we extend by summing over consecutive pair-pairs.

    Returns:
        'L' for Λ (right-handed, propeller +), 'D' for Δ (left-handed, −),
        '' if helicity undefined: < 2 chelate pairs, geometry without
        registered vectors, all triples ≤ EPS (achiral by construction).
    """
    # Iter-2C: dispatch DD / SAP / TPR to their polyhedron-specific
    # classifiers driven by the principal-improper-axis projection.  The
    # universal scalar-triple-product fall-through covers OH, TBP, SP, COH
    # (and any future polyhedron registered only in _GEOM_VECTORS without
    # a dedicated helicity_fn).
    fn = _GEOM_HELICITY_FN.get(geometry)
    if fn is not None:
        return fn(perm, chelate_pairs)

    vectors = _GEOM_VECTORS.get(geometry)
    if vectors is None:
        return ''

    # Build chelate-pair geometric data in canonical donor-list order.
    pair_data: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
    perm_list = list(perm)
    for cp in chelate_pairs:
        pair = sorted(cp)
        if len(pair) != 2:
            continue
        try:
            pos_a = perm_list.index(pair[0])
            pos_b = perm_list.index(pair[1])
        except ValueError:
            continue
        if pos_a >= len(vectors) or pos_b >= len(vectors):
            continue
        va = vectors[pos_a]
        vb = vectors[pos_b]
        c = _sub(vb, va)               # chelate-spanning vector
        m = (va[0] + vb[0], va[1] + vb[1], va[2] + vb[2])  # centroid (×2)
        pair_data.append((c, m))

    if len(pair_data) < 2:
        return ''

    # Sort pair_data canonically by centroid magnitude+direction so the
    # output is permutation-symmetry stable: swapping the pair-iteration
    # order must not flip helicity (would create non-determinism between
    # equivalent perms that just relabel chelate-pair indices).
    pair_data.sort(key=lambda cm: (round(_norm2(cm[1]), 3),
                                    round(cm[1][0], 3),
                                    round(cm[1][1], 3),
                                    round(cm[1][2], 3)))

    # Mixed product: sum over all i<j pair-pairs of
    #   ((c_i × c_j) · (m_i + m_j))
    # This is the SO(3)-invariant helicity scalar.  For 2 chelate pairs
    # this collapses to the single triple product; for 3 pairs (tris-bidentate
    # Δ/Λ classic) it is the standard tris-chelate handedness sign.
    total = 0.0
    n = len(pair_data)
    for i in range(n):
        for j in range(i + 1, n):
            c_i, m_i = pair_data[i]
            c_j, m_j = pair_data[j]
            cross_cc = _cross(c_i, c_j)
            m_sum = (m_i[0] + m_j[0], m_i[1] + m_j[1], m_i[2] + m_j[2])
            total += _dot(cross_cc, m_sum)

    if total > _EPS:
        return 'L'
    if total < -_EPS:
        return 'D'
    return ''


def helicity_aware_pairs(canonical_pairs: tuple,
                         perm: Optional[Sequence[int]],
                         chelate_pairs: Optional[Iterable[FrozenSet]],
                         geometry: str) -> tuple:
    """Append helicity tag to a canonical-form pair-tuple if defined.

    Returns ``canonical_pairs`` unchanged if helicity is undefined ('') —
    so achiral isomers keep their bit-exact pre-Iter-2 canonical form.
    """
    if perm is None or chelate_pairs is None:
        return canonical_pairs
    h = _classify_helicity(perm, chelate_pairs, geometry)
    if not h:
        return canonical_pairs
    return canonical_pairs + (('chir', h),)


__all__ = ['_classify_helicity', 'helicity_aware_pairs']
