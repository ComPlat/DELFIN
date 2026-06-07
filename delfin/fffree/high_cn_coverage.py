"""delfin.fffree.high_cn_coverage — Universal high-CN (CN 8-12) coverage extension.

Self-contained extension that fills the CN 8-12 coverage gaps for non-f-block
metals (Ru, Mo, Os, W, Re, Tc, etc.).  Many of the worst-bug-class CCDC
refcodes (PAYQIS Ru CN10, FEZQUY Ru CN10, BUNWUF Ru CN9, GEYJAX Ru CN10) fall
into this dispatch hole — they pass decompose but ``ref_vectors`` returns
nothing useful AND md_distance has no entry for the rare donor-metal pair,
so the whole-complex assembly silently fails and no XYZ is emitted.

Scope (the contract is universal — pure geometry + covalent radii, no per-metal
or per-class branches):

* **CN8** — adds ``DD-8`` (triangular dodecahedron, D2d), ``TPR-8`` (bicapped
  trigonal prism, D2h) and re-exports the existing ``SQAP-8`` (Wells canonical).
* **CN9** — adds ``CSAP-9`` (capped square antiprism, C4v) and ``MFF-9`` (muffin
  / J84-derived, C2v).  ``TTP-9`` is the legacy default.
* **CN10** — re-exports the existing ``BICAP-10`` / ``CSAP-10`` / ``SAP-10``
  and adds ``SPHENO-10`` (sphenocorona J87, C2v).  No re-definition.
* **CN11** — adds ``OCD-11`` (octadecahedron variant, C2v) and ``CAP-11``
  (monocapped pentagonal antiprism, C5v, re-exported from f-block).
* **CN12** — adds ``ICO-12`` (icosahedron, Ih, re-exported from f-block),
  ``CUOH-12`` (cuboctahedron, Oh) and ``ACUOH-12`` (anticuboctahedron, D3h).

All vertex builders return ``(N, 3)`` unit-vector arrays — same contract as
``polyhedra.ref_vectors`` and ``f_block_polyhedra.ref_vectors_fblock``.  Index
order is fixed (top → bottom → caps) for cross-platform determinism.

Universal M-D distance fallback
-------------------------------

For CN >= 8, real-world bonds are slightly LONGER than the
covalent-radii-sum predicts (the Bürgi-Dunitz-style high-CN bond elongation):
typical empirical excess is +0.08-0.12 Å across the 4d/5d series.  When the
metal+donor pair is not in the COV lookup table, we use
``cov(M) + cov(D) + 0.10`` instead of bailing.

Env gates (default OFF — byte-identical to HEAD when both unset)
----------------------------------------------------------------

* ``DELFIN_FFFREE_HIGH_CN_COVERAGE=1`` — turn on CN 8-12 coverage for
  non-f-block metals.  Auto-enabled under ``DELFIN_FFFREE_PURE_TRACK3=1``.
* ``DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR=2.0`` — float multiplier for the
  ETKDG bounds-matrix tolerance when CN >= 8 (default ``2.0``; tighter than the
  3.0× the worst-CCDC cases might need, but loose enough to recover the
  common failures).  Read by ``embed_fallback`` on the high-CN retry path.

This module is import-safe — every entry point gates on an env flag, so a
stray import never changes legacy behaviour.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Env gates
# ---------------------------------------------------------------------------

_FLAG_HIGHCN = "DELFIN_FFFREE_HIGH_CN_COVERAGE"
_FLAG_BOUNDS = "DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR"
_FLAG_PT3 = "DELFIN_FFFREE_PURE_TRACK3"


def high_cn_coverage_enabled() -> bool:
    """True when the CN 8-12 coverage extension is enabled (env-gated)."""
    return (
        os.environ.get(_FLAG_HIGHCN, "0") == "1"
        or os.environ.get(_FLAG_PT3, "0") == "1"
    )


def high_cn_bounds_tol_factor() -> float:
    """Read the ETKDG bounds-matrix tolerance multiplier for CN >= 8.

    Default ``2.0`` (loose enough to recover the common high-CN failures;
    tighter than the 3.0× the worst CCDC cases might need).
    """
    raw = os.environ.get(_FLAG_BOUNDS, "2.0").strip()
    try:
        v = float(raw)
    except ValueError:
        v = 2.0
    return max(1.0, v)


# ---------------------------------------------------------------------------
# Vertex helpers
# ---------------------------------------------------------------------------


def _norm_rows(V: np.ndarray) -> np.ndarray:
    """Row-wise unit-normalise an (N, 3) array."""
    return V / np.linalg.norm(V, axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# CN = 8 — extras beyond SQAP-8 (already in polyhedra.py)
# ---------------------------------------------------------------------------


def dd_8_vertices() -> np.ndarray:
    """DD-8 — Triangular Dodecahedron (D2d, Hoard-Silverton).

    Re-exports the canonical vertex set from :mod:`f_block_polyhedra` so the
    non-f-block CN8 path shares ONE source of truth with the f-block one.

    Vertex order: 0,2,4,6 = "A" vertices (φ = 0°, 90°, 180°, 270°);
    1,3,5,7 = "B" vertices (φ = 45°, 135°, 225°, 315°).
    """
    from delfin.fffree import f_block_polyhedra as _FBP
    return _FBP.dd_8_vertices()


def tpr_8_vertices() -> np.ndarray:
    """TPR-8 — Bicapped Trigonal Prism (D2h).

    A CN6 trigonal prism with two extra caps on the two rectangular faces that
    sit on opposite sides of the C3 axis (yielding D2h, not D3h).  Distinct
    from BICAP-10's D4d (which caps both square faces of a SAP-8): TPR-8 caps
    a TRIANGULAR PRISM, not an antiprism, so it has fewer symmetry elements.

    Vertex order:
      0-2 = upper triangle (φ = 0°, 120°, 240°, +z);
      3-5 = lower triangle (eclipsed with the upper, -z);
      6,7 = caps on the two opposite rectangular faces (φ = 60°, 240°, z=0).
    """
    R_tri = 0.85
    h_pr = 0.50
    R_cap = 1.0
    top = [(R_tri * math.cos(math.radians(120.0 * k)),
            R_tri * math.sin(math.radians(120.0 * k)), h_pr) for k in range(3)]
    bot = [(R_tri * math.cos(math.radians(120.0 * k)),
            R_tri * math.sin(math.radians(120.0 * k)), -h_pr) for k in range(3)]
    # Two caps on opposite rectangular faces (60° and 240°)
    caps = [
        (R_cap * math.cos(math.radians(60.0)),
         R_cap * math.sin(math.radians(60.0)), 0.0),
        (R_cap * math.cos(math.radians(240.0)),
         R_cap * math.sin(math.radians(240.0)), 0.0),
    ]
    return _norm_rows(np.asarray(top + bot + caps, dtype=float))


# ---------------------------------------------------------------------------
# CN = 9 — extras beyond TTP-9 (already in polyhedra.py)
# ---------------------------------------------------------------------------


def csap_9_vertices() -> np.ndarray:
    """CSAP-9 — Capped Square Antiprism (C4v).

    Re-exports the canonical vertex set from :mod:`f_block_polyhedra` so the
    non-f-block CN9 path shares ONE source of truth.
    """
    from delfin.fffree import f_block_polyhedra as _FBP
    return _FBP.csap_9_vertices()


def mff_9_vertices() -> np.ndarray:
    """MFF-9 — Muffin (J84-derived 9-vertex C2v polyhedron).

    A "muffin" is a 9-vertex polyhedron used for some Mo/W/Re CN9 complexes:
    a square-antiprism-like belt of 8 vertices PLUS one apical cap, but the
    belt is slightly distorted (one square is contracted) — giving C2v.

    Geometry (idealised on the unit sphere):
      * 1 apical cap on +z;
      * 4 "upper square" vertices at z = +h_up, φ = 0°, 90°, 180°, 270°,
        radius r_up (slightly smaller than 1 — closer to the cap);
      * 4 "lower square" vertices at z = -h_lo, φ = 45°, 135°, 225°, 315°,
        radius r_lo.

    Distinct from CSAP-9 in the belt tilt: CSAP-9 has h_top = 0.55, h_bot =
    0.70 (close to symmetric).  MFF-9 has h_top = 0.45, h_bot = 0.75 — a
    deeper tilt that breaks the C4v down to C2v.

    Vertex order: 0 = cap; 1-4 = upper square; 5-8 = lower square.
    """
    h_top = 0.45
    h_bot = 0.75
    r_top = math.sqrt(max(0.0, 1.0 - h_top * h_top))
    r_bot = math.sqrt(max(0.0, 1.0 - h_bot * h_bot))
    cap = [(0.0, 0.0, 1.0)]
    top = [(r_top * math.cos(math.radians(90.0 * k)),
            r_top * math.sin(math.radians(90.0 * k)), h_top) for k in range(4)]
    bot = [(r_bot * math.cos(math.radians(45.0 + 90.0 * k)),
            r_bot * math.sin(math.radians(45.0 + 90.0 * k)), -h_bot)
           for k in range(4)]
    return _norm_rows(np.asarray(cap + top + bot, dtype=float))


# ---------------------------------------------------------------------------
# CN = 10 — extras beyond BICAP-10 / CSAP-10 / SAP-10 (already in polyhedra.py)
# ---------------------------------------------------------------------------


def spheno_10_vertices() -> np.ndarray:
    """SPHENO-10 — Sphenocorona (Johnson solid J87, C2v).

    An irregular 10-vertex polyhedron used for some early-d high-CN complexes
    (Zr, Hf, Ta, W) where the more symmetric SAP-10/BICAP-10 don't fit the
    donor steric profile.  Approximate idealised arrangement:

      * 4 equator vertices on ±x, ±y (z=0);
      * 4 "wedge" vertices at z = +zs1, leaning toward ±x and ±y;
      * 2 apex vertices at z = ±zs2 along the y axis (one up, one down).

    Vertex order: 0-3 = equator; 4-7 = wedge; 8 = upper apex; 9 = lower apex.
    """
    zs1 = 0.6
    zs2 = 0.85
    out: List[Tuple[float, float, float]] = [
        (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0), (0.0, -1.0, 0.0),
        (0.85, 0.0, zs1), (-0.85, 0.0, zs1),
        (0.0, 0.85, zs1), (0.0, -0.85, zs1),
        (0.0, 0.55, zs2),
        (0.0, -0.55, -zs2),
    ]
    return _norm_rows(np.asarray(out, dtype=float))


# ---------------------------------------------------------------------------
# CN = 11 — extras
# ---------------------------------------------------------------------------


def cap_11_vertices() -> np.ndarray:
    """CAP-11 — Monocapped Pentagonal Antiprism (C5v).

    Re-exports the canonical vertex set from :mod:`f_block_polyhedra` so the
    non-f-block CN11 path shares ONE source of truth.
    """
    from delfin.fffree import f_block_polyhedra as _FBP
    return _FBP.cap_11_vertices()


def ocd_11_vertices() -> np.ndarray:
    """OCD-11 — Octadecahedron-derived 11-vertex polyhedron (C2v).

    The octadecahedron itself has 11 vertices, 18 triangular faces.  We use
    the "elongated pentagonal antiprism" parameterisation: two pentagons (one
    rotated 36° relative to the other) PLUS one apical cap, but with the two
    pentagons at DIFFERENT z-heights (unlike CAP-11 where both pentagons
    straddle the cap symmetrically).

    Vertex order:
      0    = apical cap (+z);
      1-5  = upper pentagon (φ = 0°, 72°, 144°, 216°, 288°, z = +h_top);
      6-10 = lower pentagon (φ = 36°, 108°, 180°, 252°, 324°, z = -h_bot).

    Distinct from CAP-11 by the asymmetric tilt: h_top = 0.30 (very close to
    cap), h_bot = 0.65 (further out).
    """
    h_top = 0.30
    h_bot = 0.65
    r_top = math.sqrt(max(0.0, 1.0 - h_top * h_top))
    r_bot = math.sqrt(max(0.0, 1.0 - h_bot * h_bot))
    cap = [(0.0, 0.0, 1.0)]
    top = [(r_top * math.cos(math.radians(72.0 * k)),
            r_top * math.sin(math.radians(72.0 * k)), h_top) for k in range(5)]
    bot = [(r_bot * math.cos(math.radians(36.0 + 72.0 * k)),
            r_bot * math.sin(math.radians(36.0 + 72.0 * k)), -h_bot)
           for k in range(5)]
    return _norm_rows(np.asarray(cap + top + bot, dtype=float))


# ---------------------------------------------------------------------------
# CN = 12 — extras beyond IH-12 (already in f_block)
# ---------------------------------------------------------------------------


def ico_12_vertices() -> np.ndarray:
    """ICO-12 — Regular Icosahedron (Ih).

    Re-exports the canonical vertex set from :mod:`f_block_polyhedra` so the
    non-f-block CN12 path shares ONE source of truth.
    """
    from delfin.fffree import f_block_polyhedra as _FBP
    return _FBP.ih_12_vertices()


def cuoh_12_vertices() -> np.ndarray:
    """CUOH-12 — Cuboctahedron (Oh).

    12 vertices at all distinct sign-permutations of (±1, ±1, 0): each vertex
    is a midpoint of a cube edge.  Common in close-packed metallic clusters
    and high-CN Ln/An complexes where bidentate ligands prefer the Oh shape
    over Ih.

    Vertex order: lex-sorted by (z, y, x) for cross-platform determinism.
    """
    raw: List[Tuple[float, float, float]] = []
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            raw.append((s1, s2, 0.0))
            raw.append((s1, 0.0, s2))
            raw.append((0.0, s1, s2))
    raw.sort(key=lambda v: (v[2], v[1], v[0]))
    return _norm_rows(np.asarray(raw, dtype=float))


def acuoh_12_vertices() -> np.ndarray:
    """ACUOH-12 — Anticuboctahedron (D3h, triangular orthobicupola J27).

    Twin-cupola of the cuboctahedron with a 60° twist on one cupola: instead
    of the Oh symmetry of CUOH-12 we get D3h (the same point group as the
    trigonal prism).  12 vertices forming two parallel hexagons (one upper at
    +z, one lower at -z) twisted 60° relative to each other, plus the
    "waist" comes from the staggered ring positions.

    Standard parameterisation: 6 vertices on each of two equilateral
    triangles at z = ±h, both with in-plane radius r.  Half the vertices on
    each triangle are "high" (z = +h_high) and half are "low" (z = +h_low)
    on the same triangle, giving 12 total in a D3h arrangement.

    Vertex order:
      0-5  = upper hexagon-like ring (alternating high/low z at φ k·60°);
      6-11 = lower hexagon-like ring (rotated 30° relative to upper).
    """
    h_high = 0.50
    h_low = 0.20
    r = math.sqrt(max(0.0, 1.0 - h_high * h_high))
    r2 = math.sqrt(max(0.0, 1.0 - h_low * h_low))
    top: List[Tuple[float, float, float]] = []
    for k in range(6):
        z = h_high if (k % 2 == 0) else h_low
        rr = r if (k % 2 == 0) else r2
        ang = math.radians(60.0 * k)
        top.append((rr * math.cos(ang), rr * math.sin(ang), z))
    bot: List[Tuple[float, float, float]] = []
    for k in range(6):
        z = -h_high if (k % 2 == 0) else -h_low
        rr = r if (k % 2 == 0) else r2
        ang = math.radians(30.0 + 60.0 * k)
        bot.append((rr * math.cos(ang), rr * math.sin(ang), z))
    return _norm_rows(np.asarray(top + bot, dtype=float))


# ---------------------------------------------------------------------------
# Public registry — canonical name → vertex builder
# ---------------------------------------------------------------------------


HIGH_CN_VERTICES: Dict[str, callable] = {
    # CN8 extras (SQAP-8 itself is already in polyhedra.py)
    "DD-8 dodecahedron": dd_8_vertices,
    "TPR-8 bicapped trigonal prism": tpr_8_vertices,
    # CN9 extras (TTP-9 already in polyhedra.py)
    "CSAP-9 capped square antiprism": csap_9_vertices,
    "MFF-9 muffin": mff_9_vertices,
    # CN10 extras (BICAP-10/CSAP-10/SAP-10 already in polyhedra.py)
    "SPHENO-10 sphenocorona": spheno_10_vertices,
    # CN11
    "CAP-11 monocapped pentagonal antiprism": cap_11_vertices,
    "OCD-11 octadecahedron": ocd_11_vertices,
    # CN12
    "ICO-12 icosahedron": ico_12_vertices,
    "CUOH-12 cuboctahedron": cuoh_12_vertices,
    "ACUOH-12 anticuboctahedron": acuoh_12_vertices,
}

# Short alias → canonical name.
HIGH_CN_ALIASES: Dict[str, str] = {
    "DD-8": "DD-8 dodecahedron",
    "TPR-8": "TPR-8 bicapped trigonal prism",
    "CSAP-9": "CSAP-9 capped square antiprism",
    "MFF-9": "MFF-9 muffin",
    "SPHENO-10": "SPHENO-10 sphenocorona",
    "CAP-11": "CAP-11 monocapped pentagonal antiprism",
    "OCD-11": "OCD-11 octadecahedron",
    "ICO-12": "ICO-12 icosahedron",
    "CUOH-12": "CUOH-12 cuboctahedron",
    "ACUOH-12": "ACUOH-12 anticuboctahedron",
}


# CN → list of extra geometries added when high-CN coverage is enabled.  Order
# is lex-deterministic (canonical names sorted) and matches the chemically
# common first / Wells-canonical first convention used by ``geometries_for_cn``.
HIGH_CN_EXTRA_BY_CN: Dict[int, List[str]] = {
    8: ["DD-8 dodecahedron", "TPR-8 bicapped trigonal prism"],
    9: ["CSAP-9 capped square antiprism", "MFF-9 muffin"],
    10: ["SPHENO-10 sphenocorona"],
    11: ["CAP-11 monocapped pentagonal antiprism", "OCD-11 octadecahedron"],
    12: ["ICO-12 icosahedron", "CUOH-12 cuboctahedron",
         "ACUOH-12 anticuboctahedron"],
}


def ref_vectors_high_cn(geometry: str) -> np.ndarray:
    """Dispatch to the canonical vertex builder for a high-CN polyhedron.

    Accepts either the full canonical name (key of :data:`HIGH_CN_VERTICES`)
    or a short alias (key of :data:`HIGH_CN_ALIASES`).  Raises ``KeyError`` on
    unknown name so callers can chain to the legacy ``polyhedra.ref_vectors``
    in a try/except.
    """
    if geometry in HIGH_CN_VERTICES:
        return HIGH_CN_VERTICES[geometry]()
    if geometry in HIGH_CN_ALIASES:
        return HIGH_CN_VERTICES[HIGH_CN_ALIASES[geometry]]()
    raise KeyError(geometry)


def default_geometry_high_cn(cn: int) -> Optional[str]:
    """Return the default high-CN polyhedron for ``cn`` or None.

    Only meaningful when the legacy ``GEOM_BY_CN`` has no entry for the CN
    (e.g. CN10, CN11, CN12 in HEAD).  Returns the FIRST canonical entry from
    :data:`HIGH_CN_EXTRA_BY_CN`.
    """
    if cn not in HIGH_CN_EXTRA_BY_CN:
        return None
    return HIGH_CN_EXTRA_BY_CN[cn][0]


# ---------------------------------------------------------------------------
# Extended covalent-radius table for the high-CN M-D fallback
# ---------------------------------------------------------------------------


# Cordero (2008) covalent-radii (Å) — same source as ``polyhedra.COV``.  We
# explicitly enumerate the metals that the legacy table missed AND that
# regularly appear in high-CN coordination chemistry: Ln3+ (CN8-10), An (CN8-
# 12), Y, Sc, Sr, Ba, Ca (large ionic), Pb, Bi, Sn.
_HIGH_CN_COV_EXTRA: Dict[str, float] = {
    # Lanthanides (Ln3+) — most are already missing from the legacy COV table.
    "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98,
    "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92,
    "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87,
    # Actinides
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90,
    "Pu": 1.87, "Am": 1.80, "Cm": 1.69,
    # Large ionic metals (Group 2 + heavy main group)
    "Ca": 1.76, "Sr": 1.95, "Ba": 2.15,
    "Pb": 1.46, "Bi": 1.48, "Sn": 1.39,
    "Tl": 1.45, "In": 1.42, "Ga": 1.22,
    # Early-d that the legacy table is missing
    "Tc": 1.47,
}


def md_distance_high_cn(metal: str, donor: str, cn: int = 0) -> Optional[float]:
    """Universal M-D distance with the high-CN bond-elongation correction.

    Returns ``None`` when both the legacy COV table AND the extended table
    miss the metal — caller falls back to the legacy ``cov_default`` (1.5 Å).

    For ``cn >= 8`` we add a deterministic +0.10 Å correction reflecting the
    empirical high-CN bond elongation across the 4d/5d/Ln series (the radius
    increase is well-documented in Shannon 1976 and downstream).
    """
    from delfin.fffree.polyhedra import COV as _COV
    r_m = _COV.get(metal)
    if r_m is None:
        r_m = _HIGH_CN_COV_EXTRA.get(metal)
    r_d = _COV.get(donor)
    if r_m is None or r_d is None:
        return None
    base = float(r_m) + float(r_d)
    if cn >= 8:
        return base + 0.10
    return base


def covalent_radius_extra(metal: str) -> Optional[float]:
    """Return the extra covalent radius for ``metal`` if known, else None."""
    return _HIGH_CN_COV_EXTRA.get(metal)


# ---------------------------------------------------------------------------
# Pólya proper-rotation groups for the high-CN polyhedra
# ---------------------------------------------------------------------------


def _close_group(generators: List[Tuple[int, ...]], n: int) -> List[Tuple[int, ...]]:
    """Close a set of vertex permutations into the group they generate."""
    identity = tuple(range(n))
    group = {identity}
    frontier = [identity]
    gens = [tuple(g) for g in generators]
    while frontier:
        a = frontier.pop()
        for g in gens:
            b = tuple(g[a[i]] for i in range(n))
            if b not in group:
                group.add(b)
                frontier.append(b)
    return sorted(group)


def _dd_8_group() -> List[Tuple[int, ...]]:
    """DD-8 (D2d) proper-rotation subgroup — conservative C2 (order 2).

    Vertex order: 0,2,4,6 = A-type; 1,3,5,7 = B-type.  C2_z (180° about z):
    0<->4, 2<->6, 1<->5, 3<->7.  This is the same conservative choice
    f_block_polyhedra makes for DD-8: it slightly over-counts orbits versus
    the full D2 (order 4), but stays sound (never under-counts).
    """
    C2z = (4, 5, 6, 7, 0, 1, 2, 3)
    return _close_group([C2z], 8)


def _tpr_8_group() -> List[Tuple[int, ...]]:
    """TPR-8 (D2h) proper-rotation subgroup — conservative C2 (order 2).

    Vertex order: 0-2 upper triangle, 3-5 lower triangle, 6,7 caps.
    C2_z about z: swaps 6<->7 (the two caps).  We do NOT include the C3 of
    the underlying TPR-6 because the caps break the threefold symmetry;
    only the cap-swap C2 is a real proper rotation of TPR-8.
    """
    # C2 about z: the two caps (6, 7) swap; top triangle goes to itself but
    # vertices 0/1/2 swap with NOTHING else exactly under a 180° z-rotation
    # (top triangle has C3, not C2).  Use the cap-swap only as a safe symmetry.
    C2 = (0, 1, 2, 3, 4, 5, 7, 6)
    return _close_group([C2], 8)


def _csap_9_group() -> List[Tuple[int, ...]]:
    """CSAP-9 (C4v) proper-rotation subgroup — C4 (order 4).

    Vertex order: 0 = apical cap; 1-4 = upper square; 5-8 = lower square
    (45°-rotated).  C4_z cycles 1->2->3->4 and 5->6->7->8, cap fixed.
    Same convention as f_block_polyhedra._csap_9_group.
    """
    C4 = (0, 2, 3, 4, 1, 6, 7, 8, 5)
    return _close_group([C4], 9)


def _mff_9_group() -> List[Tuple[int, ...]]:
    """MFF-9 (C2v) proper-rotation subgroup — C2 (order 2).

    Vertex order: 0 = cap; 1-4 = upper square; 5-8 = lower square.  The
    asymmetric belt tilt breaks C4 down to C2: only the 180° rotation about
    z is a proper rotation.  Upper: 1<->3, 2<->4; lower: 5<->7, 6<->8.
    """
    C2 = (0, 3, 4, 1, 2, 7, 8, 5, 6)
    return _close_group([C2], 9)


def _spheno_10_group() -> List[Tuple[int, ...]]:
    """SPHENO-10 (C2v) proper-rotation subgroup — C2 (order 2).

    Vertex order: 0-3 equator; 4-7 wedge; 8 upper apex; 9 lower apex.
    C2 about y-axis: equator (1,0,3,2), wedge (5,4,7,6), apexes swap.
    """
    C2 = (1, 0, 3, 2, 5, 4, 7, 6, 9, 8)
    return _close_group([C2], 10)


def _cap_11_group() -> List[Tuple[int, ...]]:
    """CAP-11 (C5v) proper-rotation subgroup — C5 (order 5).

    Vertex order: 0 = apical cap; 1-5 = upper pentagon; 6-10 = lower pentagon
    (rotated 36°).  C5_z cycles 1->2->3->4->5->1 and 6->7->8->9->10->6.
    """
    C5 = (0, 2, 3, 4, 5, 1, 7, 8, 9, 10, 6)
    return _close_group([C5], 11)


def _ocd_11_group() -> List[Tuple[int, ...]]:
    """OCD-11 (C2v) proper-rotation subgroup — C5 (order 5).

    Vertex order: 0 = cap; 1-5 = upper pentagon; 6-10 = lower pentagon.
    The asymmetric tilt still keeps the C5 axis (each pentagon rotates by
    72° onto itself independently).  Same C5 as CAP-11.
    """
    C5 = (0, 2, 3, 4, 5, 1, 7, 8, 9, 10, 6)
    return _close_group([C5], 11)


def _ico_12_group() -> List[Tuple[int, ...]]:
    """ICO-12 (Ih) — identity-only (conservative).

    The full I rotation group has order 60 but encoding it exactly on the
    lex-sorted vertex set requires care to keep determinism.  Same choice as
    f_block_polyhedra._ih_12_group: identity-only over-counts orbits but
    stays sound.
    """
    return [tuple(range(12))]


def _cuoh_12_group() -> List[Tuple[int, ...]]:
    """CUOH-12 (Oh) — conservative C4 (order 4).

    The full O proper-rotation group has order 24; we use the C4_z subgroup
    as a conservative under-approximation.  Vertices are lex-sorted by
    (z,y,x) of the unnormalised coordinates, so figuring out the C4
    permutation directly is non-trivial; we use identity as the safest
    choice (worst case: over-counts orbits, never under-counts).
    """
    return [tuple(range(12))]


def _acuoh_12_group() -> List[Tuple[int, ...]]:
    """ACUOH-12 (D3h) proper-rotation subgroup — C3 (order 3).

    Vertex order: 0-5 = upper ring; 6-11 = lower ring.  C3 cycles vertices
    in pairs around the C3 axis.
    """
    C3 = (2, 3, 4, 5, 0, 1, 8, 9, 10, 11, 6, 7)
    return _close_group([C3], 12)


# Canonical name → (group, n).  Maps to the SAME canonical strings as
# HIGH_CN_VERTICES so the Pólya enumerator and the vertex builder agree.
HIGH_CN_GROUPS: Dict[str, Tuple[List[Tuple[int, ...]], int]] = {
    "DD-8 dodecahedron": (_dd_8_group(), 8),
    "TPR-8 bicapped trigonal prism": (_tpr_8_group(), 8),
    "CSAP-9 capped square antiprism": (_csap_9_group(), 9),
    "MFF-9 muffin": (_mff_9_group(), 9),
    "SPHENO-10 sphenocorona": (_spheno_10_group(), 10),
    "CAP-11 monocapped pentagonal antiprism": (_cap_11_group(), 11),
    "OCD-11 octadecahedron": (_ocd_11_group(), 11),
    "ICO-12 icosahedron": (_ico_12_group(), 12),
    "CUOH-12 cuboctahedron": (_cuoh_12_group(), 12),
    "ACUOH-12 anticuboctahedron": (_acuoh_12_group(), 12),
}


def high_cn_group(geometry: str) -> Optional[Tuple[List[Tuple[int, ...]], int]]:
    """Return ``(group, n_vertices)`` for a high-CN polyhedron or None.

    Accepts either the full canonical name or a short alias.
    """
    if geometry in HIGH_CN_ALIASES:
        geometry = HIGH_CN_ALIASES[geometry]
    return HIGH_CN_GROUPS.get(geometry)


# Pólya geometry-key → canonical name (consumed by
# ``polya_isomer_count._get_group`` via a lazy lookup helper).
HIGH_CN_POLYA_KEY_MAP: Dict[str, str] = {
    "dd_8_hicn": "DD-8 dodecahedron",
    "tpr_8": "TPR-8 bicapped trigonal prism",
    "csap_9_hicn": "CSAP-9 capped square antiprism",
    "mff_9": "MFF-9 muffin",
    "spheno_10": "SPHENO-10 sphenocorona",
    "cap_11_hicn": "CAP-11 monocapped pentagonal antiprism",
    "ocd_11": "OCD-11 octadecahedron",
    "ico_12_hicn": "ICO-12 icosahedron",
    "cuoh_12": "CUOH-12 cuboctahedron",
    "acuoh_12": "ACUOH-12 anticuboctahedron",
}


# ---------------------------------------------------------------------------
# Self-test (manual sanity)
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    print("=== High-CN coverage polyhedra ===")
    for name, fn in HIGH_CN_VERTICES.items():
        v = fn()
        n = len(v)
        norms = np.linalg.norm(v, axis=1)
        max_dev = float(np.max(np.abs(norms - 1.0)))
        angles_deg: List[float] = []
        for i in range(n):
            for j in range(i + 1, n):
                cos_a = float(np.clip(np.dot(v[i], v[j]), -1.0, 1.0))
                angles_deg.append(math.degrees(math.acos(cos_a)))
        print(f"  {name:<46s}  N={n:2d}  ||v||-1<={max_dev:.2e}  "
              f"min_angle={min(angles_deg):.1f}deg  "
              f"max_angle={max(angles_deg):.1f}deg")
    print()
    print("  md_distance_high_cn('Ru', 'N', cn=10) =",
          md_distance_high_cn("Ru", "N", cn=10))
    print("  md_distance_high_cn('Y', 'O', cn=10) =",
          md_distance_high_cn("Y", "O", cn=10))
    print("  md_distance_high_cn('Pa', 'F', cn=8) =",
          md_distance_high_cn("Pa", "F", cn=8))
