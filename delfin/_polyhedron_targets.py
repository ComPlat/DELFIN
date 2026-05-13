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

* CN=2 ``linear_2``           — D∞h
* CN=3 ``trigonal_planar``    — D3h
* CN=4 ``Td``                 — Td
* CN=4 ``sqp_4``              — D4h
* CN=4 ``see_saw``            — C2v
* CN=5 ``tbp``                — D3h
* CN=5 ``sqp_5``              — C4v
* CN=6 ``Oh``                 — Oh
* CN=6 ``trig_prism``         — D3h
* CN=7 ``pbp``                — D5h
* CN=7 ``capped_oct``         — C3v
* CN=8 ``sq_antiprism``       — D4d
* CN=8 ``cube``               — Oh
* CN=8 ``dodecahedron``       — D2d
* CN=9 ``tricapped_tp``       — D3h
* CN=10 ``bicapped_sap``      — D4d
* CN=12 ``icosahedron``       — Ih
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

    # ------------- CN = 12 -------------
    # Icosahedron (Ih): 12 unit vertices using golden ratio φ.
    phi = _GOLDEN
    raw = [
        ( 0,  1,  phi), ( 0,  1, -phi), ( 0, -1,  phi), ( 0, -1, -phi),
        ( 1,  phi, 0), ( 1, -phi, 0), (-1,  phi, 0), (-1, -phi, 0),
        ( phi, 0,  1), ( phi, 0, -1), (-phi, 0,  1), (-phi, 0, -1),
    ]
    d["icosahedron"] = _stack([_u(*xyz) for xyz in raw])

    return d


_IDEAL_VECTORS: Dict[str, np.ndarray] = _build_ideal_vectors()


# Convenient alias map: (cn, geometry) → key in _IDEAL_VECTORS
_CN_GEOM_KEY: Dict[Tuple[int, str], str] = {
    (2, "linear"):           "linear_2",
    (2, "linear_2"):         "linear_2",
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
    (9, "tricapped_tp"):     "tricapped_tp",
    (9, "TTP"):              "tricapped_tp",
    (10, "bicapped_sap"):    "bicapped_sap",
    (10, "BCSAP"):           "bicapped_sap",
    (12, "icosahedron"):     "icosahedron",
    (12, "ICOS"):            "icosahedron",
}


# Point-group classification
_POINT_GROUPS: Dict[Tuple[int, str], str] = {
    (2, "linear_2"):         "D_inf_h",
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
    (9, "tricapped_tp"):     "D3h",
    (10, "bicapped_sap"):    "D4d",
    (12, "icosahedron"):     "Ih",
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
    if cn == 12:
        return ("icosahedron", None)
    return (f"undefined_cn{cn}", None)


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
