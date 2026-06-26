"""Symmetry-aware ideal donor unit-vector tables per (CN, geometry, isomer).

This module exposes a clean API providing the target unit vectors of an ideal
coordination polyhedron for any (coordination number, geometry, isomer-label)
triple. It is consumed by the Baustein-5 PBD optimizer for Hungarian-
assignment-based donor projection.

Independent re-implementation of the geometry tables — does NOT mutate the
legacy ``_TOPO_GEOMETRY_VECTORS`` dict in ``smiles_converter.py`` so that
parallel paths can co-exist during the Hybrid-Refactor.

All vectors returned by :func:`get_ideal_donor_vectors` are unit vectors
(``||v|| = 1``).  Per-geometry tables are stored as length-1 ``np.ndarray``
rows of shape ``(cn, 3)``.

Geometry codes (Schoenflies point group in parens):

* CN=2  ``linear_2``                — D∞h
* CN=3  ``trigonal_planar``         — D3h
* CN=4  ``Td``                      — Td
* CN=4  ``sqp_4``                   — D4h
* CN=4  ``see_saw``                 — C2v
* CN=5  ``tbp``                     — D3h
* CN=5  ``sqp_5``                   — C4v
* CN=6  ``Oh``                      — Oh
* CN=6  ``trig_prism``              — D3h
* CN=7  ``pbp``                     — D5h
* CN=7  ``capped_oct``              — C3v
* CN=8  ``sq_antiprism``            — D4d
* CN=8  ``cube``                    — Oh
* CN=8  ``dodecahedron``            — D2d
* CN=8  ``bicapped_trig_antiprism`` — D3d   (lanthanide / actinide preferred)
* CN=9  ``tricapped_tp``            — D3h
* CN=9  ``capped_sap``              — C4v   (lanthanide aquo / nitrato)
* CN=10 ``bicapped_sap``            — D4d
* CN=10 ``pentag_antiprism``        — D5d   (alternative lanthanide CN10)
* CN=10 ``sphenocorona``            — C2v   (Johnson solid J87, irregular CN10)
* CN=11 ``mono_capped_pentag_aprism``— C5v  (irregular CN11, lanthanide nitrato-aquo)
* CN=12 ``icosahedron``             — Ih
* CN=12 ``cuboctahedron``           — Oh    (lanthanide hexanitrato CN12)

Lanthanide / high-CN handling
-----------------------------

The CN=8-12 catalog is selectively expanded for f-block elements (La-Lu, U-Pu)
which prefer high coordination numbers with ionic (non-directional) bonding.
The metal-aware classifier :func:`classify_geometry_from_cn_donors_with_metal`
routes lanthanides and high-CN actinides to the polyhedra that ionic-radius
chemistry favours (e.g. CN=9 capped square antiprism for [La(H2O)9]3+, CN=8
bicapped trigonal antiprism for La(NTA), CN=10 bicapped SAP for La(NO3)5 type).

The new vector tables, point groups, and alias rows are emitted unconditionally
(they are static reference data, no runtime cost).  The metal-aware *selection*
logic is gated by the environment variable ``DELFIN_HIGH_CN_POLYHEDRA`` (read
by :mod:`delfin._symmetry_detection`).  Default OFF — bit-exact identical to
the pre-patch CN>=8 dispatch.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _u(*xyz: float) -> np.ndarray:
    """Return unit-normalised np.ndarray from xyz scalars."""
    v = np.asarray(xyz, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("zero-length vector cannot be unit-normalised")
    return v / n


def _stack(vectors: List[np.ndarray]) -> np.ndarray:
    """Stack list of unit vectors into shape (n,3) ndarray."""
    return np.vstack(vectors)


# ---------------------------------------------------------------------------
# Ideal vector tables
# ---------------------------------------------------------------------------

_GOLDEN: float = (1.0 + np.sqrt(5.0)) / 2.0  # φ ≈ 1.618


def _build_ideal_vectors() -> Dict[str, np.ndarray]:
    """Construct & cache the ideal-vector tables."""
    d: Dict[str, np.ndarray] = {}

    # ------------- CN = 2 -------------
    d["linear_2"] = _stack([
        _u(1, 0, 0),
        _u(-1, 0, 0),
    ])
    # Bent / water-like CN=2 (C2v, ~104.5°).  Mainly relevant for non-d10
    # CN=2 centres (e.g. Cu(II)/Ti where the d-shell is not full and a bent
    # ligand-field stabilisation is preferred).  d10 metals (Au(I), Ag(I),
    # Cu(I), Hg(II), Tl(I)) almost always sit at 180° linear and are routed
    # to ``linear_2`` by the metal-aware classifier below.
    _BENT_HALF_DEG: float = 104.5 / 2.0   # = 52.25°
    _bent_x = np.sin(np.deg2rad(_BENT_HALF_DEG))
    _bent_y = np.cos(np.deg2rad(_BENT_HALF_DEG))
    d["bent_2"] = _stack([
        _u(_bent_x,  _bent_y, 0.0),
        _u(-_bent_x, _bent_y, 0.0),
    ])

    # ------------- CN = 3 -------------
    d["trigonal_planar"] = _stack([
        _u(1, 0, 0),
        _u(-0.5,  np.sqrt(3) / 2.0, 0),
        _u(-0.5, -np.sqrt(3) / 2.0, 0),
    ])

    # ------------- CN = 4 -------------
    d["Td"] = _stack([
        _u( 1,  1,  1),
        _u( 1, -1, -1),
        _u(-1,  1, -1),
        _u(-1, -1,  1),
    ])
    d["sqp_4"] = _stack([
        _u( 1, 0, 0), _u(-1, 0, 0),
        _u(0,  1, 0), _u(0, -1, 0),
    ])
    # See-saw / C2v (e.g. SF4-type) — 2 axial + 2 equatorial bent
    # axial on z, equatorial in xz plane, equatorial bent ~30° down off +x/-x
    d["see_saw"] = _stack([
        _u(0, 0,  1), _u(0, 0, -1),
        _u(np.cos(np.deg2rad(15.0)),  0, -np.sin(np.deg2rad(15.0))),
        _u(-np.cos(np.deg2rad(15.0)), 0, -np.sin(np.deg2rad(15.0))),
    ])

    # ------------- CN = 5 -------------
    # Trigonal bipyramidal (D3h)
    d["tbp"] = _stack([
        _u(0, 0,  1), _u(0, 0, -1),          # axial
        _u(1, 0, 0),                          # equatorial 1
        _u(-0.5,  np.sqrt(3) / 2.0, 0),       # equatorial 2
        _u(-0.5, -np.sqrt(3) / 2.0, 0),       # equatorial 3
    ])
    # Square pyramidal (C4v): apical + 4 base, base tilted ~15° below xy-plane
    tilt = np.deg2rad(15.0)
    bz = -np.sin(tilt)
    br = np.cos(tilt)
    d["sqp_5"] = _stack([
        _u(0, 0, 1),
        _u(br,  0, bz),
        _u(-br, 0, bz),
        _u(0,  br, bz),
        _u(0, -br, bz),
    ])

    # ------------- CN = 6 -------------
    d["Oh"] = _stack([
        _u(1, 0, 0), _u(-1, 0, 0),
        _u(0, 1, 0), _u(0, -1, 0),
        _u(0, 0, 1), _u(0, 0, -1),
    ])
    # Trigonal prismatic (D3h) — eclipsed top/bottom triangles at z=±zhalf.
    # zhalf chosen so vertices land on unit sphere given xy radius 1 → can't
    # both be 1; pick zhalf=1/√3 so that x²+y²+z² = (2/3)+(1/3)=1 with xy
    # radius √(2/3) ≈ 0.8165.
    zhalf = 1.0 / np.sqrt(3.0)
    rho = np.sqrt(2.0 / 3.0)
    d["trig_prism"] = _stack([
        _u(rho * 1.0,                rho * 0.0,                 zhalf),
        _u(rho * -0.5,               rho * (np.sqrt(3) / 2.0),  zhalf),
        _u(rho * -0.5,               rho * -(np.sqrt(3) / 2.0), zhalf),
        _u(rho * 1.0,                rho * 0.0,                -zhalf),
        _u(rho * -0.5,               rho * (np.sqrt(3) / 2.0), -zhalf),
        _u(rho * -0.5,               rho * -(np.sqrt(3) / 2.0),-zhalf),
    ])

    # ------------- CN = 7 -------------
    # Pentagonal bipyramidal (D5h)
    pbp_rows: List[np.ndarray] = [
        _u(0, 0,  1), _u(0, 0, -1),
    ]
    for k in range(5):
        ang = 2.0 * np.pi * k / 5.0
        pbp_rows.append(_u(np.cos(ang), np.sin(ang), 0.0))
    d["pbp"] = _stack(pbp_rows)
    # Capped octahedral (C3v) — Oh base + one face cap on (1,1,1)/√3
    d["capped_oct"] = _stack([
        _u(1, 0, 0), _u(-1, 0, 0),
        _u(0, 1, 0), _u(0, -1, 0),
        _u(0, 0, 1), _u(0, 0, -1),
        _u(1, 1, 1),  # cap
    ])

    # ------------- CN = 8 -------------
    # Square antiprism (D4d): top square rotated 45° relative to bottom.
    # On unit sphere with xy-radius cos θ and z = ±sin θ; pick θ=π/4 so both
    # equal √2/2 → unit length.
    h = np.sqrt(2.0) / 2.0          # = sin(π/4) and = cos(π/4)
    sap_top: List[np.ndarray] = []
    sap_bot: List[np.ndarray] = []
    for k in range(4):
        ang_top = (np.pi / 2.0) * k + np.pi / 4.0  # rotated 45°
        ang_bot = (np.pi / 2.0) * k
        sap_top.append(_u(h * np.cos(ang_top), h * np.sin(ang_top),  h))
        sap_bot.append(_u(h * np.cos(ang_bot), h * np.sin(ang_bot), -h))
    d["sq_antiprism"] = _stack(sap_top + sap_bot)
    # Cube (Oh): 8 vertices at (±1,±1,±1)/√3
    cube_rows = [
        _u( 1,  1,  1), _u( 1,  1, -1),
        _u( 1, -1,  1), _u( 1, -1, -1),
        _u(-1,  1,  1), _u(-1,  1, -1),
        _u(-1, -1,  1), _u(-1, -1, -1),
    ]
    d["cube"] = _stack(cube_rows)
    # Dodecahedron / triangulated D2d (CN=8): two interpenetrating tetragonal
    # bisphenoids; pragmatic table (approximate D2d) with all unit vectors.
    a = np.deg2rad(36.85)
    b = np.deg2rad(69.46)
    r_a = np.sin(a); z_a = np.cos(a)
    r_b = np.sin(b); z_b = np.cos(b)
    d["dodecahedron"] = _stack([
        _u(r_a, 0,   z_a), _u(-r_a, 0,   z_a),
        _u(0,   r_a, -z_a), _u(0,   -r_a, -z_a),
        _u(r_b, 0,  -z_b), _u(-r_b, 0,  -z_b),
        _u(0,   r_b,  z_b), _u(0,   -r_b,  z_b),
    ])

    # ------------- CN = 9 -------------
    # Tricapped trigonal prism (D3h): TP base (6) + 3 equatorial caps.
    # Re-use trig_prism geometry then add 3 caps at z=0 perpendicular to TP
    # rectangular faces (between adjacent vertical edges).
    ttp_rows: List[np.ndarray] = []
    # 6 prism vertices on slightly compressed prism (z=±0.5, xy-r normalised)
    z9 = 0.5
    r9 = np.sqrt(1.0 - z9 * z9)  # ensure unit length
    for sign in (+1, -1):
        for k in range(3):
            ang = 2.0 * np.pi * k / 3.0
            ttp_rows.append(
                _u(r9 * np.cos(ang), r9 * np.sin(ang), sign * z9)
            )
    # 3 equatorial caps offset by π/3 relative to prism vertices
    for k in range(3):
        ang = 2.0 * np.pi * k / 3.0 + np.pi / 3.0
        ttp_rows.append(_u(np.cos(ang), np.sin(ang), 0.0))
    d["tricapped_tp"] = _stack(ttp_rows)

    # ------------- CN = 8 (extra, lanthanide / actinide) -------------
    # Bicapped trigonal antiprism (D3d) — common for 8-coordinate lanthanides
    # with small donors (e.g. La in NTA chelates, U in fluorides).
    # Geometry: trigonal antiprism (top triangle z=+z8, bottom rotated 60°,
    # z=-z8) + 2 axial caps on ±z.  Pick z8 so prism vertices sit on unit
    # sphere; caps are unit by construction.
    z8 = 1.0 / np.sqrt(3.0)
    r8 = np.sqrt(1.0 - z8 * z8)
    bcta_rows: List[np.ndarray] = []
    for k in range(3):
        ang_top = 2.0 * np.pi * k / 3.0
        ang_bot = 2.0 * np.pi * k / 3.0 + np.pi / 3.0  # rotated 60°
        bcta_rows.append(_u(r8 * np.cos(ang_top), r8 * np.sin(ang_top),  z8))
        bcta_rows.append(_u(r8 * np.cos(ang_bot), r8 * np.sin(ang_bot), -z8))
    bcta_rows.append(_u(0, 0,  1))
    bcta_rows.append(_u(0, 0, -1))
    d["bicapped_trig_antiprism"] = _stack(bcta_rows)

    # ------------- CN = 9 (extra, lanthanide aquo / nitrato) -------------
    # Capped square antiprism (C4v) — apex cap on +z plus SAP body offset
    # downward so all nine vectors are unit and roughly evenly spaced.
    # Body z = -tilt9, xy-radius r_body; top/bottom squares both below the
    # apex cap and rotated 45° relative to each other.
    tilt9_top = np.deg2rad(20.0)  # top square 20° below apex
    tilt9_bot = np.deg2rad(70.0)  # bottom square 70° below apex
    rt = np.sin(tilt9_top); zt = np.cos(tilt9_top)
    rb = np.sin(tilt9_bot); zb = np.cos(tilt9_bot)
    csap_rows: List[np.ndarray] = [_u(0, 0, 1)]  # apex
    for k in range(4):
        ang_t = (np.pi / 2.0) * k + np.pi / 4.0  # rotated 45°
        ang_b = (np.pi / 2.0) * k
        csap_rows.append(_u(rt * np.cos(ang_t), rt * np.sin(ang_t),  zt))
        csap_rows.append(_u(rb * np.cos(ang_b), rb * np.sin(ang_b), -zb))
    d["capped_sap"] = _stack(csap_rows)

    # ------------- CN = 10 -------------
    # Bicapped square antiprism (D4d): SAP + 2 axial caps along ±z.
    # Scale SAP equator so caps + equator are all unit; place equator at z=±h10
    # with xy-radius r10.
    h10 = 0.5
    r10 = np.sqrt(1.0 - h10 * h10)
    bcsap_rows: List[np.ndarray] = []
    for k in range(4):
        ang_top = (np.pi / 2.0) * k + np.pi / 4.0
        ang_bot = (np.pi / 2.0) * k
        bcsap_rows.append(_u(r10 * np.cos(ang_top), r10 * np.sin(ang_top),  h10))
        bcsap_rows.append(_u(r10 * np.cos(ang_bot), r10 * np.sin(ang_bot), -h10))
    bcsap_rows.append(_u(0, 0,  1))
    bcsap_rows.append(_u(0, 0, -1))
    d["bicapped_sap"] = _stack(bcsap_rows)

    # Pentagonal antiprism (D5d): two parallel pentagons at z=±zp10, the top
    # one rotated by 36° (π/5) relative to the bottom one.
    zp10 = 1.0 / np.sqrt(3.0)
    rp10 = np.sqrt(1.0 - zp10 * zp10)
    pap_rows: List[np.ndarray] = []
    for k in range(5):
        ang_t = 2.0 * np.pi * k / 5.0 + np.pi / 5.0
        ang_b = 2.0 * np.pi * k / 5.0
        pap_rows.append(_u(rp10 * np.cos(ang_t), rp10 * np.sin(ang_t),  zp10))
        pap_rows.append(_u(rp10 * np.cos(ang_b), rp10 * np.sin(ang_b), -zp10))
    d["pentag_antiprism"] = _stack(pap_rows)

    # Sphenocorona (J87, C2v) — irregular Johnson solid with 10 vertices.
    # Approximate idealised arrangement: 4 vertices forming a square belt at
    # z=0, 4 wedge vertices above (z=+zs1) tilted toward ±x, 2 apex vertices
    # at z=+zs2 along ±y.  All normalised to unit vectors.  This is a
    # pragmatic target rather than the exact Johnson-solid geometry, sized so
    # all directions are reasonably distinct.
    zs1 = 0.6
    zs2 = 0.85
    sphen_rows: List[np.ndarray] = [
        _u( 1.0,  0.0, 0.0), _u(-1.0,  0.0, 0.0),  # equator x
        _u( 0.0,  1.0, 0.0), _u( 0.0, -1.0, 0.0),  # equator y
        _u( 0.85,  0.0,  zs1), _u(-0.85,  0.0,  zs1),
        _u( 0.0,   0.85, zs1), _u( 0.0,  -0.85, zs1),
        _u( 0.0,   0.55, zs2),
        _u( 0.0,  -0.55, -zs2),  # one lobe descends to retain volume
    ]
    d["sphenocorona"] = _stack(sphen_rows)

    # ------------- CN = 11 (irregular, lanthanide nitrato-aquo) -------------
    # Monocapped pentagonal antiprism (C5v): apex cap on +z, then 2 pentagons
    # below in antiprism arrangement.
    z11_a = np.deg2rad(40.0)
    z11_b = np.deg2rad(75.0)
    ra11 = np.sin(z11_a); za11 = np.cos(z11_a)
    rb11 = np.sin(z11_b); zb11 = np.cos(z11_b)
    cn11_rows: List[np.ndarray] = [_u(0, 0, 1)]
    for k in range(5):
        ang_t = 2.0 * np.pi * k / 5.0
        ang_b = 2.0 * np.pi * k / 5.0 + np.pi / 5.0
        cn11_rows.append(_u(ra11 * np.cos(ang_t), ra11 * np.sin(ang_t),  za11))
        cn11_rows.append(_u(rb11 * np.cos(ang_b), rb11 * np.sin(ang_b), -zb11))
    d["mono_capped_pentag_aprism"] = _stack(cn11_rows)

    # ------------- CN = 12 -------------
    # Icosahedron (Ih): 12 unit vertices using golden ratio φ.
    phi = _GOLDEN
    raw = [
        ( 0,  1,  phi), ( 0,  1, -phi), ( 0, -1,  phi), ( 0, -1, -phi),
        ( 1,  phi, 0), ( 1, -phi, 0), (-1,  phi, 0), (-1, -phi, 0),
        ( phi, 0,  1), ( phi, 0, -1), (-phi, 0,  1), (-phi, 0, -1),
    ]
    d["icosahedron"] = _stack([_u(*xyz) for xyz in raw])

    # Cuboctahedron (Oh): 12 vertices at all permutations of (±1, ±1, 0).
    # Common ideal for [Ln(NO3)6]3- type CN=12 with bidentate nitrates.
    cubocta_raw = [
        ( 1,  1, 0), ( 1, -1, 0), (-1,  1, 0), (-1, -1, 0),
        ( 1, 0,  1), ( 1, 0, -1), (-1, 0,  1), (-1, 0, -1),
        ( 0,  1,  1), ( 0,  1, -1), ( 0, -1,  1), ( 0, -1, -1),
    ]
    d["cuboctahedron"] = _stack([_u(*xyz) for xyz in cubocta_raw])

    return d


_IDEAL_VECTORS: Dict[str, np.ndarray] = _build_ideal_vectors()


# Convenient alias map: (cn, geometry) → key in _IDEAL_VECTORS
_CN_GEOM_KEY: Dict[Tuple[int, str], str] = {
    (2, "linear"):           "linear_2",
    (2, "linear_2"):          "linear_2",
    (2, "LIN"):               "linear_2",
    (2, "bent"):              "bent_2",
    (2, "bent_2"):            "bent_2",
    (2, "BENT"):              "bent_2",
    (3, "trigonal_planar"):  "trigonal_planar",
    (3, "TP"):               "trigonal_planar",
    (4, "Td"):               "Td",
    (4, "tetrahedral"):      "Td",
    (4, "sqp"):              "sqp_4",
    (4, "sqp_4"):            "sqp_4",
    (4, "square_planar"):    "sqp_4",
    (4, "see_saw"):          "see_saw",
    (4, "SS"):               "see_saw",
    (5, "tbp"):              "tbp",
    (5, "trig_bipyramidal"): "tbp",
    (5, "sqp"):              "sqp_5",
    (5, "sqp_5"):            "sqp_5",
    (5, "square_pyramidal"): "sqp_5",
    (6, "Oh"):               "Oh",
    (6, "octahedral"):       "Oh",
    (6, "trig_prism"):       "trig_prism",
    (6, "TPR"):              "trig_prism",
    (7, "pbp"):              "pbp",
    (7, "pentag_bipyramidal"): "pbp",
    (7, "capped_oct"):       "capped_oct",
    (7, "COH"):              "capped_oct",
    (8, "sq_antiprism"):     "sq_antiprism",
    (8, "SAP"):              "sq_antiprism",
    (8, "cube"):             "cube",
    (8, "dodecahedron"):     "dodecahedron",
    (8, "DD"):               "dodecahedron",
    (8, "bicapped_trig_antiprism"): "bicapped_trig_antiprism",
    (8, "BCTA"):             "bicapped_trig_antiprism",
    (9, "tricapped_tp"):     "tricapped_tp",
    (9, "TTP"):              "tricapped_tp",
    (9, "capped_sap"):       "capped_sap",
    (9, "CSAP"):             "capped_sap",
    (10, "bicapped_sap"):    "bicapped_sap",
    (10, "BCSAP"):           "bicapped_sap",
    (10, "pentag_antiprism"): "pentag_antiprism",
    (10, "PAP"):             "pentag_antiprism",
    (10, "sphenocorona"):    "sphenocorona",
    (10, "J87"):             "sphenocorona",
    (11, "mono_capped_pentag_aprism"): "mono_capped_pentag_aprism",
    (11, "MCPAP"):           "mono_capped_pentag_aprism",
    (12, "icosahedron"):     "icosahedron",
    (12, "ICOS"):            "icosahedron",
    (12, "cuboctahedron"):   "cuboctahedron",
    (12, "COCT"):            "cuboctahedron",
}


# Point-group classification
_POINT_GROUPS: Dict[Tuple[int, str], str] = {
    (2, "linear_2"):         "D_inf_h",
    (2, "bent_2"):           "C2v",
    (3, "trigonal_planar"):  "D3h",
    (4, "Td"):               "Td",
    (4, "sqp_4"):            "D4h",
    (4, "see_saw"):          "C2v",
    (5, "tbp"):              "D3h",
    (5, "sqp_5"):            "C4v",
    (6, "Oh"):               "Oh",
    (6, "trig_prism"):       "D3h",
    (7, "pbp"):              "D5h",
    (7, "capped_oct"):       "C3v",
    (8, "sq_antiprism"):     "D4d",
    (8, "cube"):             "Oh",
    (8, "dodecahedron"):     "D2d",
    (8, "bicapped_trig_antiprism"): "D3d",
    (9, "tricapped_tp"):     "D3h",
    (9, "capped_sap"):       "C4v",
    (10, "bicapped_sap"):    "D4d",
    (10, "pentag_antiprism"): "D5d",
    (10, "sphenocorona"):    "C2v",
    (11, "mono_capped_pentag_aprism"): "C5v",
    (12, "icosahedron"):     "Ih",
    (12, "cuboctahedron"):   "Oh",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _resolve_key(cn: int, geometry: str) -> str:
    """Map (cn, geometry-alias) → canonical key in _IDEAL_VECTORS."""
    key = _CN_GEOM_KEY.get((cn, geometry))
    if key is not None:
        return key
    # tolerant fallback — accept canonical key directly
    if geometry in _IDEAL_VECTORS:
        vecs = _IDEAL_VECTORS[geometry]
        if vecs.shape[0] == cn:
            return geometry
    raise KeyError(
        f"no ideal-vector table for (cn={cn}, geometry={geometry!r})"
    )


def get_ideal_donor_vectors(
    cn: int,
    geometry: str,
    isomer_label: Optional[str] = None,
) -> np.ndarray:
    """Return ideal donor unit vectors for the requested polyhedron.

    Args:
        cn: coordination number (2-12).
        geometry: geometry name, e.g., ``'Oh'``, ``'Td'``, ``'sqp'``, ``'tbp'``.
            Both modern (``'sq_antiprism'``) and legacy DELFIN codes
            (``'SAP'``) are accepted.
        isomer_label: optional textual label such as ``'cis'``, ``'trans'``,
            ``'fac'``, ``'mer'``.  Currently consumed only by
            :func:`_arrange_donors_by_label` (no-op without ``donor_types``).

    Returns:
        np.ndarray of shape ``(cn, 3)`` — unit vectors pointing from the
        metal centre to the ideal donor positions.
    """
    key = _resolve_key(cn, geometry)
    vecs = _IDEAL_VECTORS[key]
    # Defensive copy so callers cannot mutate the table.
    return vecs.copy()


def get_target_point_group(
    cn: int,
    geometry: str,
    isomer_label: Optional[str] = None,
) -> str:
    """Return Schoenflies symbol of the ideal point group."""
    key = _resolve_key(cn, geometry)
    pg = _POINT_GROUPS.get((cn, key))
    if pg is None:
        return "C1"
    # mixed-ligand sets typically break the parent symmetry
    if isomer_label:
        return f"{pg}*"  # marker indicating reduced effective symmetry
    return pg


# ---------------------------------------------------------------------------
# d10 linear CN=2 special-case (Au(I), Ag(I), Cu(I), Hg(II), Tl(I))
# ---------------------------------------------------------------------------

# d10 closed-shell metals that *always* prefer linear (180°) geometry at CN=2.
# These ions have completely filled d-shells and no LF stabilisation from a
# bent arrangement; the σ + π back-bonding pattern strongly favours D∞h.
#
# Group 11 (Cu, Ag, Au) at +1 → d10
# Group 12 (Zn, Cd, Hg) at +2 → d10  — but only Hg routinely shows CN=2
# Group 13 (Tl) at +1         → s2 (post-transition) — Tl(I) is linear too
#
# Cu and Tl require an explicit oxidation-state check because Cu(II)/Tl(III)
# are not d10 and may show bent or other CN=2 patterns (rare but seen).
# Ag/Au/Hg are accepted in their nominal oxidation states (Ag+1, Au+1, Hg+2)
# without an additional charge gate because:
#   * Ag(II)/Au(III) almost never form CN=2 complexes (Ag(II) is rare; Au(III)
#     prefers square-planar CN=4)
#   * Hg(I) forms Hg-Hg dimers (one Hg-Hg + one external donor) and the
#     external linear geometry pattern is preserved anyway.
_D10_LINEAR_CN2_METALS_ALWAYS: frozenset = frozenset({"Ag", "Au", "Hg"})
_D10_LINEAR_CN2_METALS_CHARGED: Dict[str, int] = {
    "Cu": +1,
    "Tl": +1,
}


def _is_d10_linear_cn2(
    metal_symbol: Optional[str],
    cn: int,
    metal_formal_charge: Optional[int] = None,
) -> bool:
    """Return True iff this is a d¹⁰ CN=2 centre that prefers linear (180°).

    The d¹⁰ closed-shell electron configuration removes any ligand-field
    driving force for bending; consequently Au(I), Ag(I), Cu(I), Hg(II) and
    Tl(I) almost universally exhibit linear (D∞h) CN=2 coordination
    ([Au(CN)2]⁻, [Ag(NH3)2]⁺, [Cu(MeCN)2]⁺, [Hg(CN)2], [Tl(OR)2]⁻).

    The classifier deliberately avoids SMILES patterns and ligand-class
    keywords: only the metal element symbol + formal charge of the metal
    atom + coordination number are consulted.

    Args:
        metal_symbol: element symbol of the metal atom (e.g. ``'Au'``).
        cn: coordination number (donor count) at this metal.
        metal_formal_charge: optional integer formal charge on the metal
            atom (e.g. +1 for Cu(I)).  Required for Cu/Tl; ignored for
            Ag/Au/Hg.  ``None`` → conservatively only Ag/Au/Hg are treated
            as d¹⁰ linear.

    Returns:
        ``True`` if the metal is d¹⁰ and CN=2 (linear target should be
        forced).  ``False`` otherwise.
    """
    if cn != 2:
        return False
    if not metal_symbol:
        return False
    if metal_symbol in _D10_LINEAR_CN2_METALS_ALWAYS:
        return True
    required_q = _D10_LINEAR_CN2_METALS_CHARGED.get(metal_symbol)
    if required_q is None:
        return False
    if metal_formal_charge is None:
        # Charge unknown → cannot prove d¹⁰; default to "not d¹⁰".
        return False
    return int(metal_formal_charge) == required_q


def classify_geometry_from_cn_donors_with_metal(
    cn: int,
    donor_types: List[str],
    metal_symbol: Optional[str] = None,
    metal_formal_charge: Optional[int] = None,
) -> Tuple[str, Optional[str]]:
    """Metal-aware variant of :func:`classify_geometry_from_cn_donors`.

    Adds three universal special-cases (no SMILES patterns):

    * **CN=2 d¹⁰ linear** — Au(I), Ag(I), Cu(I), Hg(II), Tl(I) are routed to
      ``"linear_2"`` (180°, D∞h).  Already the default for CN=2 — explicit
      here so that the env-gated bent-CN=2 path below cannot accidentally
      flip them.

    * **CN=2 non-d¹⁰ bent** — env-gated by ``DELFIN_LINEAR_CN2`` (default 0).
      When enabled and the metal is *not* in the d¹⁰ linear-CN=2 set, CN=2
      centres are classified as ``"bent_2"`` (~104.5°, C2v).  This lets
      Cu(II)/Ti/early-TM CN=2 species relax to a bent target instead of
      being forced to 180°.

      With the default (env=0) this routine is bit-exact identical to
      :func:`classify_geometry_from_cn_donors` for every (cn, donor_types)
      pair — d¹⁰ metals and non-d¹⁰ metals both keep returning
      ``("linear_2", None)`` at CN=2.

    * **CN=8-12 high-CN ionic metals** — env-gated by
      ``DELFIN_HIGH_CN_POLYHEDRA`` (default 0).  When enabled and the metal
      is in the lanthanide / high-CN-actinide / large-ionic set
      (:data:`_LANTHANIDES`, :data:`_HIGH_CN_ACTINIDES`,
      :data:`_HIGH_CN_IONIC_OTHER`), CN=8-12 are routed to the
      lanthanide-preferred polyhedra (BCTA / CSAP / PAP-or-BCSAP / MCPAP /
      cuboctahedron) via :func:`_high_cn_lanthanide_geometry`.

      With the default (env=0) the CN=8-12 dispatch is bit-exact identical
      to :func:`classify_geometry_from_cn_donors`.

    Higher CN (3-12) delegates to the existing metal-agnostic classifier
    when none of the special-cases above apply.
    """
    import os

    if cn == 2:
        # Default OFF → bit-exact identical to legacy classifier.
        if os.environ.get("DELFIN_LINEAR_CN2", "0") not in ("1", "true", "True"):
            return ("linear_2", None)
        if _is_d10_linear_cn2(metal_symbol, cn, metal_formal_charge):
            return ("linear_2", None)
        return ("bent_2", None)

    # High-CN lanthanide / actinide / large-ionic override (env-gated).
    if cn >= 8:
        if os.environ.get("DELFIN_HIGH_CN_POLYHEDRA", "0") in ("1", "true", "True"):
            if _is_high_cn_ionic_metal(metal_symbol):
                override = _high_cn_lanthanide_geometry(cn, donor_types)
                if override is not None:
                    return override

    return classify_geometry_from_cn_donors(cn, donor_types)


def classify_geometry_from_cn_donors(
    cn: int,
    donor_types: List[str],
) -> Tuple[str, Optional[str]]:
    """Heuristic geometry + isomer-label classification from CN + donor list.

    Uses common chemistry defaults — no metal-specific d-electron logic.

    Args:
        cn: coordination number.
        donor_types: list of donor-atom element symbols (length should == cn).

    Returns:
        ``(geometry, isomer_label_or_None)``.
    """
    if cn == 2:
        return ("linear_2", None)
    if cn == 3:
        return ("trigonal_planar", None)
    if cn == 4:
        # default tetrahedral; sqp picked only via explicit context
        return ("Td", _maybe_cis_trans(donor_types, geometry="Td"))
    if cn == 5:
        return ("tbp", None)
    if cn == 6:
        iso = _maybe_octahedral_isomer(donor_types)
        return ("Oh", iso)
    if cn == 7:
        return ("pbp", None)
    if cn == 8:
        return ("sq_antiprism", None)
    if cn == 9:
        return ("tricapped_tp", None)
    if cn == 10:
        return ("bicapped_sap", None)
    if cn == 11:
        return ("mono_capped_pentag_aprism", None)
    if cn == 12:
        return ("icosahedron", None)
    return (f"undefined_cn{cn}", None)


# ---------------------------------------------------------------------------
# Metal-group taxonomy for high-CN polyhedron preference (lanthanides etc.)
# ---------------------------------------------------------------------------


# Lanthanides La-Lu (atomic numbers 57-71).  Ionic, prefer high CN, no strong
# directional d-electron preference.
_LANTHANIDES: frozenset = frozenset({
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
})

# Actinides up to Cm — those routinely seen in CCDC at high CN, ionic-like.
_HIGH_CN_ACTINIDES: frozenset = frozenset({
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
})

# Heavy / large alkaline-earth + Y/Sc/Pb/Bi that tolerate high CN with
# predominantly ionic bonding.
_HIGH_CN_IONIC_OTHER: frozenset = frozenset({
    "Y", "Sc", "Sr", "Ba", "Ca", "Pb", "Bi",
})


def _is_high_cn_ionic_metal(metal_sym: Optional[str]) -> bool:
    """Return True if metal prefers ionic / high-CN polyhedra.

    Used by :func:`classify_geometry_from_cn_donors_with_metal` to route
    CN=8-12 selections to lanthanide-style polyhedra, gated by the env
    variable ``DELFIN_HIGH_CN_POLYHEDRA``.
    """
    if not metal_sym:
        return False
    sym = metal_sym.strip()
    return (
        sym in _LANTHANIDES
        or sym in _HIGH_CN_ACTINIDES
        or sym in _HIGH_CN_IONIC_OTHER
    )


def _high_cn_lanthanide_geometry(
    cn: int, donor_types: List[str]
) -> Optional[Tuple[str, Optional[str]]]:
    """Return the lanthanide-preferred ``(geometry, isomer)`` for CN=8-12.

    Returns ``None`` if no override applies (CN outside 8-12, or no entry
    defined).  Otherwise the caller substitutes the returned tuple for the
    metal-agnostic default.

    Preferences (universal — driven by CN + donor-element composition only,
    never by SMILES or ligand-class strings):

    * CN=8  → ``bicapped_trig_antiprism`` (D3d) — closer to experimental
      Ln-NTA / U-fluoride / Th-fluoride 8-coord geometries than the
      late-d default ``sq_antiprism``.
    * CN=9  → ``capped_sap`` (C4v) — the canonical [Ln(H2O)9]3+ aquo
      geometry; also fits many Ln(NO3)3(L)x complexes better than D3h-TTP.
    * CN=10 → ``pentag_antiprism`` (D5d) for predominantly small (C/F/H)
      donor sets; otherwise ``bicapped_sap`` (D4d) is retained.
    * CN=11 → ``mono_capped_pentag_aprism`` (C5v) — common geometry for
      [Ln(NO3)4(H2O)3]- type species.
    * CN=12 → ``cuboctahedron`` (Oh) — preferred for hexa-bidentate
      species such as [Ln(NO3)6]3- where chelate bite forces 12 short
      M-O contacts in a near-Oh shell.
    """
    if cn == 8:
        return ("bicapped_trig_antiprism", None)
    if cn == 9:
        return ("capped_sap", None)
    if cn == 10:
        if donor_types:
            small = sum(1 for d in donor_types if d in ("F", "C", "H"))
            if small >= max(1, len(donor_types) // 2):
                return ("pentag_antiprism", None)
        return ("bicapped_sap", None)
    if cn == 11:
        return ("mono_capped_pentag_aprism", None)
    if cn == 12:
        return ("cuboctahedron", None)
    return None


# ---------------------------------------------------------------------------
# Isomer arrangement helpers
# ---------------------------------------------------------------------------


def _maybe_cis_trans(donor_types: List[str], geometry: str) -> Optional[str]:
    """Heuristic cis/trans assignment for CN=4 mixed sets — None by default."""
    if not donor_types:
        return None
    uniq = sorted(set(donor_types))
    if len(uniq) == 2 and donor_types.count(uniq[0]) == 2:
        return "cis"
    return None


def _maybe_octahedral_isomer(donor_types: List[str]) -> Optional[str]:
    """Heuristic isomer labelling for CN=6 mixed ligand sets.

    Returns one of ``'cis'``, ``'trans'``, ``'fac'``, ``'mer'`` or ``None``.
    """
    if not donor_types or len(donor_types) != 6:
        return None
    counts = {x: donor_types.count(x) for x in set(donor_types)}
    # Pattern A2B4 → cis vs trans (default cis = more common in CCDC stats)
    if sorted(counts.values()) == [2, 4]:
        return "cis"
    # Pattern A3B3 → fac vs mer (default fac slightly more common for
    # tridentate ligands; mer for linear-tridentate scaffolds)
    if sorted(counts.values()) == [3, 3]:
        return "fac"
    return None


def _arrange_donors_by_label(
    base_vectors: np.ndarray,
    donor_types: List[str],
    isomer_label: Optional[str],
) -> np.ndarray:
    """Reorder ``base_vectors`` so donor types match an isomer label.

    Currently implemented for the most common Oh isomers (A2B4 cis/trans,
    A3B3 fac/mer).  All other label inputs return ``base_vectors`` unchanged.
    """
    if isomer_label is None or donor_types is None:
        return base_vectors
    n = base_vectors.shape[0]
    if len(donor_types) != n:
        return base_vectors
    if n != 6:
        return base_vectors  # only Oh handled here

    counts = {x: donor_types.count(x) for x in set(donor_types)}
    minor = min(counts, key=counts.get)

    # Octahedral axis layout (matches _IDEAL_VECTORS["Oh"]):
    # 0=+x, 1=-x, 2=+y, 3=-y, 4=+z, 5=-z
    # cis: positions (0, 2)   adjacent (90°)
    # trans: positions (0, 1) opposite (180°)
    # fac: positions (0, 2, 4) — three mutually-cis vertices
    # mer: positions (0, 1, 2) — two trans + one cis

    target_minor_positions: Optional[List[int]] = None
    if sorted(counts.values()) == [2, 4]:
        if isomer_label == "cis":
            target_minor_positions = [0, 2]
        elif isomer_label == "trans":
            target_minor_positions = [0, 1]
    elif sorted(counts.values()) == [3, 3]:
        if isomer_label == "fac":
            target_minor_positions = [0, 2, 4]
        elif isomer_label == "mer":
            target_minor_positions = [0, 1, 2]

    if target_minor_positions is None:
        return base_vectors

    # Build permutation: minor donor type → target_minor_positions,
    # the rest fill remaining positions in input order.
    perm: List[int] = [-1] * n
    minor_target_iter = iter(target_minor_positions)
    rest_positions = [p for p in range(n) if p not in target_minor_positions]
    rest_iter = iter(rest_positions)
    for di, dt in enumerate(donor_types):
        slot = next(minor_target_iter) if dt == minor else next(rest_iter)
        perm[di] = slot

    out = np.zeros_like(base_vectors)
    for donor_idx, slot in enumerate(perm):
        out[donor_idx] = base_vectors[slot]
    return out


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


def _self_test() -> None:
    """Verify every vector in every table is unit-length."""
    print("== _polyhedron_targets._self_test ==")
    for key, vecs in _IDEAL_VECTORS.items():
        norms = np.linalg.norm(vecs, axis=1)
        np.testing.assert_allclose(
            norms, np.ones_like(norms), atol=1e-9,
            err_msg=f"non-unit vectors in {key!r}",
        )
        print(f"  {key:>18s}  cn={vecs.shape[0]:2d}  ||v||=ok")

    # smoke: every (cn, geom) alias resolves
    for (cn, alias), canon in _CN_GEOM_KEY.items():
        v = get_ideal_donor_vectors(cn, alias)
        assert v.shape == (cn, 3), (cn, alias, v.shape)
        pg = get_target_point_group(cn, alias)
        assert isinstance(pg, str) and pg
    print("== all aliases resolve, all vectors unit-length ==")

    # demo classifier
    for cn, donors in [
        (2, ["N", "N"]),
        (3, ["O", "O", "O"]),
        (4, ["N", "N", "N", "N"]),
        (4, ["N", "N", "Cl", "Cl"]),
        (5, ["N"] * 5),
        (6, ["N"] * 6),
        (6, ["N", "N", "Cl", "Cl", "Cl", "Cl"]),
        (6, ["N", "N", "N", "Cl", "Cl", "Cl"]),
        (7, ["O"] * 7),
        (8, ["O"] * 8),
        (9, ["O"] * 9),
        (10, ["O"] * 10),
        (12, ["O"] * 12),
    ]:
        geom, iso = classify_geometry_from_cn_donors(cn, donors)
        pg = get_target_point_group(cn, geom, iso)
        print(f"  cn={cn:2d}  donors={donors}  →  geom={geom!r}  iso={iso}  pg={pg}")

    # demo arrange-by-label (Oh A2B4 trans)
    base = get_ideal_donor_vectors(6, "Oh")
    donors = ["N", "N", "Cl", "Cl", "Cl", "Cl"]
    arr = _arrange_donors_by_label(base, donors, "trans")
    # first two donors are N → should be on opposite axes (180°)
    cos_NN = float(np.dot(arr[0], arr[1]))
    print(f"  trans-Cl2N4  N…N cos = {cos_NN:+.3f}  (expect -1.0)")
    np.testing.assert_allclose(cos_NN, -1.0, atol=1e-9)


if __name__ == "__main__":
    _self_test()
