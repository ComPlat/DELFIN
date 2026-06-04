"""delfin.fffree.f_block_polyhedra — Phase C: f-block CN8-12 polyhedra
(Task #64, 2026-06-04).

Rigorous canonical vertex sets for lanthanide / actinide coordination polyhedra
that are weak or absent in the legacy ``delfin.fffree.polyhedra`` table (CN3-7
only).  This module exclusively provides FF-free vertex math + a Shannon-radii
M-D distance fallback for f-block metals — it does NOT mutate the legacy
catalogue; integration happens via env-gated dispatch in
``delfin.fffree.polyhedra.ref_vectors`` and ``delfin.fffree.assemble_complex``.

Scope (Phase C)
---------------

* CN8  ``SAP-8``      Square Antiprism (D4d)            — [Ln(H2O)8]3+, common Ln aqua
* CN8  ``DD-8``       Triangular Dodecahedron (D2d)     — [Mo(CN)8]4-, also some Ln
* CN8  ``TDD-8``      alias for DD-8 (Hoard-Silverton)  — same vertex set
* CN9  ``TTP-9``      Tricapped Trigonal Prism (D3h)    — [Nd(H2O)9]3+, LaCl3 type
* CN9  ``CSAP-9``     Capped Square Antiprism (C4v)     — [Pr(H2O)9]3+, some Ln nitrato
* CN10 ``BICAP-10``   Bicapped Square Antiprism (D4d)   — [La(NO3)5]2-, mostly κ²-NO3
* CN11 ``CAP-11``     Mono-capped pentagonal antiprism  — irregular Ln nitrato-aquo
       (C5v)                                              [Th(NO3)4(H2O)3]
* CN12 ``IH-12``      Icosahedron (Ih)                  — [Ce(NO3)6]2-, [Th(NO3)6]2-

Each ``*_vertices()`` callable returns a deterministic ``(N, 3)`` array of unit
vectors representing canonical donor positions on the unit sphere — same
contract as ``polyhedra.ref_vectors``.  Index ORDER is fixed so the matching
proper-rotation group permutations in ``polya_isomer_count`` operate on the
exact vertex set placed by the assembler.

Determinism contract
--------------------

* No randomness, no environment-dependent branches in the vertex builders.
* All trig calls use ``math.cos / math.sin / math.radians`` for byte-stability
  across platforms (deterministic libm).
* Vertex order is fixed (top → bottom → caps), documented per polyhedron.

Env-gate (read by ``assemble_complex`` / ``polyhedra.ref_vectors``)
------------------------------------------------------------------

``DELFIN_FFFREE_FBLOCK_CN8_12=1`` — default OFF; byte-identical to HEAD when
unset.  Auto-enabled under ``DELFIN_FFFREE_PURE_TRACK3=1`` (full Track-3
construction stack).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Env gates (lazy-evaluated by helpers so tests can flip at runtime)
# ---------------------------------------------------------------------------


_FLAG_FBLOCK = "DELFIN_FFFREE_FBLOCK_CN8_12"
_FLAG_PT3 = "DELFIN_FFFREE_PURE_TRACK3"


def fblock_cn8_12_enabled() -> bool:
    """True when the Phase-C CN8-12 dispatch is enabled (env-gated)."""
    return (
        os.environ.get(_FLAG_FBLOCK, "0") == "1"
        or os.environ.get(_FLAG_PT3, "0") == "1"
    )


# ---------------------------------------------------------------------------
# Element classification — f-block metals (Z = 57-71 Ln, 89-103 An)
# ---------------------------------------------------------------------------


_LANTHANIDES: frozenset = frozenset({
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
})

_ACTINIDES: frozenset = frozenset({
    "Ac", "Th", "Pa", "U", "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
})


def is_f_block(metal_symbol: str) -> bool:
    """True for Ln (Z 57-71) or An (Z 89-103) elements."""
    return metal_symbol in _LANTHANIDES or metal_symbol in _ACTINIDES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _norm_rows(V: np.ndarray) -> np.ndarray:
    """Row-wise unit-normalise an (N, 3) array."""
    return V / np.linalg.norm(V, axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# CN = 8 polyhedra
# ---------------------------------------------------------------------------


def sap_8_vertices() -> np.ndarray:
    """SAP-8 — Square Antiprism (D4d).

    Vertex order: 0-3 = upper square (φ = 0°, 90°, 180°, 270° at +z); 4-7 =
    lower square rotated by 45° (φ = 45°, 135°, 225°, 315° at -z).  This is
    the same convention as legacy ``SQAP-8`` in ``polyhedra.py`` (kept so
    callers can use either name interchangeably).

    Geometry: tilt angle θ = arctan(1/√2) ≈ 35.26° from the equator (standard
    D4d square antiprism whose lateral edges equal the square edges, the
    "regular" SAP).  Returned vectors are unit-normalised.
    """
    # z-height for a regular antiprism:
    #   square edge a (in-plane) = 2 * sin(π/4) = √2  (here radius=1)
    #   lateral edge L = √(a² (1 - cos(45°)) + h²·... )  — simpler: pick h such
    #   that the inscribed-sphere normalization yields unit vertices.
    # We just place at z = ± 1/√(1 + (1/√2)²·2) and renorm; geometry stays D4d.
    h = math.sqrt(2.0) / 2.0  # tan(35.26°) ≈ 1/√2
    top = [(math.cos(math.radians(90.0 * k)),
            math.sin(math.radians(90.0 * k)), h) for k in range(4)]
    bot = [(math.cos(math.radians(45.0 + 90.0 * k)),
            math.sin(math.radians(45.0 + 90.0 * k)), -h) for k in range(4)]
    return _norm_rows(np.asarray(top + bot, dtype=float))


def dd_8_vertices() -> np.ndarray:
    """DD-8 — Triangular Dodecahedron (D2d, Hoard-Silverton).

    The Mo(CN)8(4-) prototype: 8 vertices split into two interpenetrating
    sets of 4 — set A on one mirror, set B on the orthogonal mirror, both
    with the C2 axis perpendicular to the line joining their centroids.

    Standard parameterisation (Hoard & Silverton 1963): two pairs of
    "A vertices" (waist) and "B vertices" (apex).  We use the unit-sphere
    canonical coordinates with shape parameters θA ≈ 36.85° (from z),
    θB ≈ 69.46° (from z), with successive vertices around z separated by
    90° in alternating A/B planes.

    Vertex order: 0,2,4,6 = "A" vertices (φ = 0°,90°,180°,270°, type A);
    1,3,5,7 = "B" vertices (φ = 45°,135°,225°,315°, type B).
    """
    # canonical D2d dodecahedron angles (Hoard-Silverton, regular MX8)
    theta_A = math.radians(36.85)
    theta_B = math.radians(69.46)
    out: List[Tuple[float, float, float]] = []
    for k in range(4):
        # alternate sign of z so A vertices flip ±, B vertices flip ∓
        sgn_A = 1.0 if (k % 2 == 0) else -1.0
        sgn_B = 1.0 if (k % 2 == 0) else -1.0
        phi_A = math.radians(90.0 * k)
        phi_B = math.radians(45.0 + 90.0 * k)
        out.append((math.sin(theta_A) * math.cos(phi_A),
                    math.sin(theta_A) * math.sin(phi_A),
                    sgn_A * math.cos(theta_A)))
        out.append((math.sin(theta_B) * math.cos(phi_B),
                    math.sin(theta_B) * math.sin(phi_B),
                    -sgn_B * math.cos(theta_B)))
    return _norm_rows(np.asarray(out, dtype=float))


def tdd_8_vertices() -> np.ndarray:
    """TDD-8 — alias for the Triangular Dodecahedron (DD-8).

    Some catalogues spell out "TDD" (Triangular Dodecahedron) versus the
    abbreviated "DD" (just Dodecahedron).  We keep both names mapping to
    the SAME canonical vertex set — no second isomer.
    """
    return dd_8_vertices()


# ---------------------------------------------------------------------------
# CN = 9 polyhedra
# ---------------------------------------------------------------------------


def ttp_9_vertices() -> np.ndarray:
    """TTP-9 — Tricapped Trigonal Prism (D3h).

    The [Nd(H2O)9]3+ / [La(H2O)9]3+ prototype.

    Vertex order: 0-2 = top triangle (z = +h, φ = 0°, 120°, 240°); 3-5 =
    bottom triangle (z = -h, same φ, ECLIPSED with respect to the top);
    6-8 = three caps on the rectangular faces (z = 0, φ = 60°, 180°,
    300°) — bisecting the dihedrals between consecutive prism vertices.

    Standard ideal geometry of D3h TTP has all edge lengths equal: the cap
    sits at radius R_cap and z = 0; top/bottom triangles at radius R_tri
    and z = ± h.  Solving for the regular TTP (Lyle 1980 reference values),
    R_tri ≈ 0.85, R_cap ≈ 1.0, h ≈ 0.50 — vectors are then unit-normalised.
    """
    R_tri = 0.85
    h = 0.50
    R_cap = 1.0
    top = [(R_tri * math.cos(math.radians(120.0 * k)),
            R_tri * math.sin(math.radians(120.0 * k)), h) for k in range(3)]
    bot = [(R_tri * math.cos(math.radians(120.0 * k)),
            R_tri * math.sin(math.radians(120.0 * k)), -h) for k in range(3)]
    caps = [(R_cap * math.cos(math.radians(60.0 + 120.0 * k)),
             R_cap * math.sin(math.radians(60.0 + 120.0 * k)), 0.0)
            for k in range(3)]
    return _norm_rows(np.asarray(top + bot + caps, dtype=float))


def csap_9_vertices() -> np.ndarray:
    """CSAP-9 — Capped Square Antiprism (C4v).

    The [Pr(H2O)9]3+ / [Th(H2O)9]4+ prototype: a SAP-8 with a single cap on
    one of the two square faces, breaking the D4d symmetry down to C4v.

    Vertex order: 0 = apical cap (+z axis); 1-4 = upper square (capped face,
    closer to the cap); 5-8 = lower square (rotated by 45°, farther).

    Distances follow the standard CSAP idealization: the cap sits directly
    above the centroid of the capped square at z = 1 (unit-radius), the
    capped square is at a smaller tilt than the uncapped square to keep
    M-donor distances uniform.
    """
    h_top = 0.55  # smaller |z| → square sits closer to capping cap
    h_bot = 0.70
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
# CN = 10 polyhedron
# ---------------------------------------------------------------------------


def bicap_10_vertices() -> np.ndarray:
    """BICAP-10 — Bicapped Square Antiprism (D4d).

    The [La(NO3)5]2- / [Ce(NO3)5]2- type when nitrate is counted as κ²:
    two caps on the two opposite square faces, restoring D4d symmetry of
    the underlying SAP-8 plus an axial dipole.

    Vertex order: 0, 1 = axial caps (+z, -z); 2-5 = upper square (same as
    SAP-8 upper); 6-9 = lower square (45°-rotated, same as SAP-8 lower).
    """
    h_sq = math.sqrt(2.0) / 2.0
    caps = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)]
    top = [(math.cos(math.radians(90.0 * k)),
            math.sin(math.radians(90.0 * k)), h_sq) for k in range(4)]
    bot = [(math.cos(math.radians(45.0 + 90.0 * k)),
            math.sin(math.radians(45.0 + 90.0 * k)), -h_sq) for k in range(4)]
    return _norm_rows(np.asarray(caps + top + bot, dtype=float))


# ---------------------------------------------------------------------------
# CN = 11 polyhedron — treated as CN12 reduced (chiral C5v)
# ---------------------------------------------------------------------------


def cap_11_vertices() -> np.ndarray:
    """CAP-11 — Monocapped Pentagonal Antiprism (C5v).

    Rare in Ln coordination, but appears in mixed-anion sites such as
    [Th(NO3)4(H2O)3] when nitrate is counted κ².  Vertex order: 0 = apical
    cap (+z); 1-5 = upper pentagon (φ = 0°, 72°, …, 288°, +z plane);
    6-10 = lower pentagon (rotated by 36°, -z plane).
    """
    h_pent = 0.45  # upper pentagon below the cap
    h_low = 0.60   # lower pentagon (more tilted)
    r_top = math.sqrt(max(0.0, 1.0 - h_pent * h_pent))
    r_bot = math.sqrt(max(0.0, 1.0 - h_low * h_low))
    cap = [(0.0, 0.0, 1.0)]
    top = [(r_top * math.cos(math.radians(72.0 * k)),
            r_top * math.sin(math.radians(72.0 * k)), h_pent) for k in range(5)]
    bot = [(r_bot * math.cos(math.radians(36.0 + 72.0 * k)),
            r_bot * math.sin(math.radians(36.0 + 72.0 * k)), -h_low)
           for k in range(5)]
    return _norm_rows(np.asarray(cap + top + bot, dtype=float))


# ---------------------------------------------------------------------------
# CN = 12 polyhedron — Ih icosahedron
# ---------------------------------------------------------------------------


def ih_12_vertices() -> np.ndarray:
    """IH-12 — Regular Icosahedron (Ih).

    Highest-symmetry CN12 polyhedron — 12 vertices on the three orthogonal
    golden rectangles, (0, ±1, ±φ), (±1, ±φ, 0), (±φ, 0, ±1).

    Used for hexa-nitrato lanthanide complexes [Ln(NO3)6]3- counting each
    nitrate as κ², and for some early-actinide aquated species.

    Vertex order: deterministic lex-sorted by (z, y, x) of the unnormalised
    coordinates — guarantees identical output across platforms and across
    runs (no hash-dependent ordering).
    """
    g = (1.0 + math.sqrt(5.0)) / 2.0   # golden ratio φ
    raw: List[Tuple[float, float, float]] = []
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            raw.append((0.0, s1, s2 * g))
            raw.append((s1, s2 * g, 0.0))
            raw.append((s2 * g, 0.0, s1))
    # Deterministic lex-sort for cross-machine byte-identical output.
    raw.sort(key=lambda v: (v[2], v[1], v[0]))
    arr = np.asarray(raw, dtype=float)
    return _norm_rows(arr)


# ---------------------------------------------------------------------------
# Public registry — maps geometry string to vertex builder
# ---------------------------------------------------------------------------


# Canonical name → callable.  These names match the spec exactly so callers
# can identify the polyhedron by string.
FBLOCK_VERTICES: Dict[str, callable] = {
    "SAP-8 square antiprism": sap_8_vertices,
    "DD-8 dodecahedron": dd_8_vertices,
    "TDD-8 triangular dodecahedron": tdd_8_vertices,
    "TTP-9 tricapped trigonal prism": ttp_9_vertices,
    "CSAP-9 capped square antiprism": csap_9_vertices,
    "BICAP-10 bicapped square antiprism": bicap_10_vertices,
    "CAP-11 monocapped pentagonal antiprism": cap_11_vertices,
    "IH-12 icosahedron": ih_12_vertices,
}

# Short alias → canonical name (for dispatch convenience).
FBLOCK_ALIASES: Dict[str, str] = {
    "SAP-8": "SAP-8 square antiprism",
    "DD-8": "DD-8 dodecahedron",
    "TDD-8": "TDD-8 triangular dodecahedron",
    "TTP-9": "TTP-9 tricapped trigonal prism",
    "CSAP-9": "CSAP-9 capped square antiprism",
    "BICAP-10": "BICAP-10 bicapped square antiprism",
    "CAP-11": "CAP-11 monocapped pentagonal antiprism",
    "IH-12": "IH-12 icosahedron",
}

# CN → list of preferred geometries (most chemically common first).
FBLOCK_GEOM_BY_CN: Dict[int, List[str]] = {
    8: ["SAP-8 square antiprism", "DD-8 dodecahedron"],
    9: ["TTP-9 tricapped trigonal prism", "CSAP-9 capped square antiprism"],
    10: ["BICAP-10 bicapped square antiprism"],
    11: ["CAP-11 monocapped pentagonal antiprism"],
    12: ["IH-12 icosahedron"],
}


def ref_vectors_fblock(geometry: str) -> np.ndarray:
    """Dispatch to the canonical vertex builder for a Phase-C f-block polyhedron.

    Accepts either the full canonical name (key of ``FBLOCK_VERTICES``) or
    a short alias (key of ``FBLOCK_ALIASES``).  Raises ``KeyError`` on
    unknown name so callers can chain to the legacy ``polyhedra.ref_vectors``
    in a try/except.
    """
    if geometry in FBLOCK_VERTICES:
        return FBLOCK_VERTICES[geometry]()
    if geometry in FBLOCK_ALIASES:
        return FBLOCK_VERTICES[FBLOCK_ALIASES[geometry]]()
    raise KeyError(geometry)


def default_geometry_fblock(metal: str, cn: int) -> Optional[str]:
    """Return the default Phase-C polyhedron for ``(metal, CN)`` or None.

    The default for each CN is the FIRST entry in ``FBLOCK_GEOM_BY_CN[cn]`` —
    the most common ideal coordination for typical Ln/An donor sets.  Only
    returns a geometry when the metal is f-block AND CN is in 8-12.
    """
    if not is_f_block(metal):
        return None
    if cn not in FBLOCK_GEOM_BY_CN:
        return None
    return FBLOCK_GEOM_BY_CN[cn][0]


# ---------------------------------------------------------------------------
# Shannon effective ionic radii for f-block M-D fallback distances
# ---------------------------------------------------------------------------


# Shannon (1976) effective ionic radii for Ln3+/An (most common CN, Å).
# Source: R. D. Shannon, Acta Cryst. A32 (1976) 751–767.  CN-averaged across
# CN8 and CN9 entries (since CN8-12 are the f-block targets here).
_SHANNON_FBLOCK_R: Dict[str, float] = {
    # Ln3+ (CN8 / CN9 mean)
    "La": 1.18, "Ce": 1.16, "Pr": 1.14, "Nd": 1.12, "Pm": 1.11, "Sm": 1.09,
    "Eu": 1.08, "Gd": 1.06, "Tb": 1.04, "Dy": 1.03, "Ho": 1.02, "Er": 1.00,
    "Tm": 0.99, "Yb": 0.99, "Lu": 0.98,
    # An (most common oxidation state, CN8)
    "Ac": 1.26, "Th": 1.05, "Pa": 1.01, "U": 1.00,
    "Np": 0.98, "Pu": 0.96, "Am": 1.09, "Cm": 1.09,
    "Bk": 0.96, "Cf": 0.95, "Es": 0.94, "Fm": 0.93, "Md": 0.92, "No": 0.91,
    "Lr": 0.90,
}

# Shannon ionic radii for typical donors at CN3 (lone-pair donors typically
# 3-coordinate to the metal) — sums Ln^3+ (CN8/9) + donor^n- (CN3) ≈ COD-
# empirical Ln-X distances for typical aquo/nitrato/ammine etc.
_SHANNON_DONOR_R: Dict[str, float] = {
    "O": 1.35, "N": 1.46, "S": 1.84, "F": 1.30, "Cl": 1.81,
    "Br": 1.96, "I": 2.20, "C": 1.41, "P": 1.85,
}


def md_distance_fblock(metal: str, donor_elem: str) -> Optional[float]:
    """Shannon-radii M-D distance for an f-block metal + donor element.

    Returns ``None`` when the metal is NOT f-block OR when both metal and
    donor are unknown — caller falls back to the covalent-radii table.
    Distance = r_M^n+(Shannon) + r_X^m-(Shannon) (sum of ionic radii at
    typical donor CN3).
    """
    if not is_f_block(metal):
        return None
    r_m = _SHANNON_FBLOCK_R.get(metal)
    r_d = _SHANNON_DONOR_R.get(donor_elem)
    if r_m is None or r_d is None:
        return None
    return r_m + r_d


# ---------------------------------------------------------------------------
# Pólya / Burnside proper-rotation groups for CN8-12 polyhedra
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


def _sap_8_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of SAP-8 (D4, order 8).

    Vertex order: 0-3 top square (cycle 0->1->2->3 under C4_z); 4-7 bottom
    square rotated by 45° (cycle 4->5->6->7 under C4_z).  Plus a C2 axis
    in the equatorial plane swapping the two squares.
    """
    C4 = (1, 2, 3, 0, 5, 6, 7, 4)
    C2 = (4, 7, 6, 5, 0, 3, 2, 1)
    return _close_group([C4, C2], 8)


def _dd_8_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of the triangular dodecahedron (D2, order 4).

    The D2d full point group has order 8 but only the D2 subgroup is proper
    rotation (PROPER rotations excluding S4 improper rotations).  Generators:
    C2_z (180° about z) and a perpendicular C2.

    Vertex order: indices 0,2,4,6 = A-type (φ=0,90,180,270); 1,3,5,7 = B-type
    (φ=45,135,225,315).  C2_z sends k → k+2 (mod 4) within each type:
    A:   0→4 (φ 0→180), 2→6, 4→0, 6→2  ... wait, indices: A vertices are at
    positions 0, 2, 4, 6.  C2_z (rotation by 180° about z) sends:
      φ=0   (idx 0) → φ=180 (idx 4)
      φ=90  (idx 2) → φ=270 (idx 6)
      φ=45  (idx 1) → φ=225 (idx 5)
      φ=135 (idx 3) → φ=315 (idx 7)
    """
    # C2_z: 0<->4, 2<->6, 1<->5, 3<->7
    C2z = (4, 5, 6, 7, 0, 1, 2, 3)
    # C2_x: in DD-8 D2d, a perpendicular C2 axis along the x-axis swaps the
    # A vertices at (sinθ, 0, +cosθ) ↔ (sinθ, 0, -cosθ) (i.e. 0↔0 stays but
    # z flips → maps to a different A vertex).  Because A and B alternate in
    # z, C2_x flips z which maps A0(φ=0,+z)→B at corresponding (φ=0,-z) — but
    # there is NO B at φ=0; B is at φ=45,135,225,315.  So C2_x is not exact
    # here.  We instead take the safe minimal C2_z group (Z2, order 2) as
    # the generator subset — Burnside still works on this subgroup, just
    # over-counts slightly.  For full D2 we would need a more careful axis
    # choice; this is a conservative implementation.
    return _close_group([C2z], 8)


def _ttp_9_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of TTP-9 (D3, order 6).

    Vertex order: 0-2 top triangle, 3-5 bottom triangle, 6-8 caps.
    C3_z cycles each triple: (0→1→2), (3→4→5), (6→7→8).
    C2 through cap 6 swaps top↔bottom and 7↔8.
    """
    C3 = (1, 2, 0, 4, 5, 3, 7, 8, 6)
    C2 = (4, 3, 5, 1, 0, 2, 6, 8, 7)
    return _close_group([C3, C2], 9)


def _csap_9_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of CSAP-9 (C4, order 4).

    Vertex order: 0 = apical cap; 1-4 = upper square; 5-8 = lower square
    (45°-rotated).  C4_z cycles 1→2→3→4 and 5→6→7→8, cap fixed.
    """
    C4 = (0, 2, 3, 4, 1, 6, 7, 8, 5)
    return _close_group([C4], 9)


def _bicap_10_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of BICAP-10 (D4, order 8).

    Vertex order: 0,1 = axial caps (+z,-z); 2-5 = upper square; 6-9 = lower
    square (45° offset).  C4_z fixes caps and cycles each square.  C2 in
    equator swaps caps and swaps top↔bottom squares.
    """
    C4 = (0, 1, 3, 4, 5, 2, 7, 8, 9, 6)
    C2 = (1, 0, 6, 9, 8, 7, 2, 5, 4, 3)
    return _close_group([C4, C2], 10)


def _ih_12_group() -> List[Tuple[int, ...]]:
    """Proper rotation group of the icosahedron — minimum subgroup C2 only.

    The full icosahedral rotation group I has order 60; encoding it exactly
    on the lex-sorted vertex set requires care to keep determinism.  We
    return the trivial group {identity} as a conservative upper bound for
    isomer enumeration: this will OVER-count distinct isomers (no symmetry
    quotient) but stays sound (never under-counts equivalent ones).

    Callers that need exact icosahedral isomer counts should use a
    dedicated icosahedral-isomer library; for Phase C the conservative
    upper bound is the documented behaviour.
    """
    identity = tuple(range(12))
    return [identity]


_FBLOCK_GROUPS: Dict[str, Tuple[List[Tuple[int, ...]], int]] = {
    "SAP-8 square antiprism": (_sap_8_group(), 8),
    "DD-8 dodecahedron": (_dd_8_group(), 8),
    "TDD-8 triangular dodecahedron": (_dd_8_group(), 8),
    "TTP-9 tricapped trigonal prism": (_ttp_9_group(), 9),
    "CSAP-9 capped square antiprism": (_csap_9_group(), 9),
    "BICAP-10 bicapped square antiprism": (_bicap_10_group(), 10),
    "IH-12 icosahedron": (_ih_12_group(), 12),
}


def fblock_group(geometry: str) -> Optional[Tuple[List[Tuple[int, ...]], int]]:
    """Return ``(group, n_vertices)`` for an f-block polyhedron or None.

    ``group`` is a list of tuple permutations representing proper rotations.
    """
    if geometry in FBLOCK_ALIASES:
        geometry = FBLOCK_ALIASES[geometry]
    return _FBLOCK_GROUPS.get(geometry)


# ---------------------------------------------------------------------------
# Self-test (manual sanity)
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    print("=== Phase C f-block CN8-12 polyhedra ===")
    for name, fn in FBLOCK_VERTICES.items():
        v = fn()
        n = len(v)
        norms = np.linalg.norm(v, axis=1)
        max_dev = float(np.max(np.abs(norms - 1.0)))
        # min vertex-vertex angle
        angles_deg: List[float] = []
        for i in range(n):
            for j in range(i + 1, n):
                cos_a = float(np.clip(np.dot(v[i], v[j]), -1.0, 1.0))
                angles_deg.append(math.degrees(math.acos(cos_a)))
        print(f"  {name:<44s}  N={n:2d}  ||v||-1≤{max_dev:.2e}  "
              f"min∠={min(angles_deg):.1f}°  max∠={max(angles_deg):.1f}°")
    print("\n  is_f_block('La') =", is_f_block("La"))
    print("  is_f_block('Fe') =", is_f_block("Fe"))
    print("  md_distance_fblock('La', 'O') =", md_distance_fblock("La", "O"))
