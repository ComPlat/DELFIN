"""Reference coordination polyhedra (unit vectors) + covalent radii for the
metal-FF-free builder (delfin.fffree).

Self-contained: ideal CN4/5/6 polyhedron vertex sets (tetrahedron, square planar,
trigonal bipyramid, square pyramid, octahedron, trigonal prism) and a covalent-radii
table.  Metal-donor distances use covalent-radii sums.
"""
from __future__ import annotations
import math
import numpy as np


def _norm_rows(V: np.ndarray) -> np.ndarray:
    return V / np.linalg.norm(V, axis=1, keepdims=True)


def _ref_polyhedra():
    R = {}
    t = 1 / math.sqrt(3)
    # Phase G (2026-05-31): CN2 L-2 linear D∞h symmetry.
    # 2 donors 180° apart on metal axis (Cu(I), Ag(I), Au(I), Hg(II)).
    R[("CN2", "L-2 linear")] = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    # CN3 (iter-32c, User 2026-05-28 ADUMOD: Pd CN3 built as Td=109.5°/linear=180°
    # instead of correct SP-3 trigonal-planar 120° or d8 T-shape 90°/180°).
    R[("CN3", "SP-3 trigonal planar")] = np.array(
        [[1.0, 0.0, 0.0],
         [-0.5, math.sqrt(3) / 2, 0.0],
         [-0.5, -math.sqrt(3) / 2, 0.0]])
    R[("CN3", "T-3 T-shape")] = np.array(
        [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    R[("CN4", "T-4 tetrahedron")] = np.array(
        [[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]])
    R[("CN4", "SP-4 square planar")] = np.array(
        [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
    R[("CN5", "TBP-5 trigonal bipyramid")] = np.array(
        [[0, 0, 1], [0, 0, -1], [1, 0, 0],
         [-0.5, math.sqrt(3) / 2, 0], [-0.5, -math.sqrt(3) / 2, 0]])
    # basal vertices in 90deg-cyclic order (45,135,225,315) so the proper-rotation
    # C4 in polya_isomer_count._spy_group (1->2->3->4) is a real geometric symmetry of
    # this vertex set -> chelate cis-edge enumeration & isomer dedup are consistent with
    # placement (the index space here IS the one assemble_from_config places into).
    R[("CN5", "SPY-5 square pyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [1, 1, 0.2], [-1, 1, 0.2], [-1, -1, 0.2], [1, -1, 0.2]], float))
    R[("CN6", "OC-6 octahedron")] = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    R[("CN6", "TPR-6 trigonal prism")] = _norm_rows(np.array(
        [[1, 0, 0.7], [-0.5, math.sqrt(3) / 2, 0.7], [-0.5, -math.sqrt(3) / 2, 0.7],
         [1, 0, -0.7], [-0.5, math.sqrt(3) / 2, -0.7], [-0.5, -math.sqrt(3) / 2, -0.7]]))
    # --- High-CN polyhedra (CN7-9). Vertex INDEX ORDER is chosen to match the
    # proper-rotation generators in polya_isomer_count (_pentagonal_bipyramid_group /
    # _square_antiprism_group / _tricapped_trigonal_prism_group) so isomer dedup is a
    # real geometric symmetry of this vertex set (same contract as SPY-5/TPR-6 above).
    # CN7 PB: idx 0,1 = axial (+z,-z); idx 2-6 = equatorial regular pentagon (0,72,...,288 deg).
    _pent = [[math.cos(2 * math.pi * k / 5), math.sin(2 * math.pi * k / 5), 0.0] for k in range(5)]
    R[("CN7", "PB-7 pentagonal bipyramid")] = _norm_rows(np.array(
        [[0, 0, 1], [0, 0, -1]] + _pent, float))
    # CN8 square antiprism: idx 0-3 = top square (0,90,180,270 deg, +z); idx 4-7 = bottom
    # square (45,135,225,315 deg, -z) -- the 45deg stagger that defines the antiprism.
    _h8 = 0.62
    _top = [[math.cos(math.pi * k / 2), math.sin(math.pi * k / 2), _h8] for k in range(4)]
    _bot = [[math.cos(math.pi * (k + 0.5) / 2 + 0.0), math.sin(math.pi * (k + 0.5) / 2), -_h8] for k in range(4)]
    # bottom at 45,135,225,315: angle = 45 + 90k
    _bot = [[math.cos(math.radians(45 + 90 * k)), math.sin(math.radians(45 + 90 * k)), -_h8] for k in range(4)]
    R[("CN8", "SQAP-8 square antiprism")] = _norm_rows(np.array(_top + _bot, float))
    # CN9 tricapped trigonal prism: idx 0-2 = top triangle (0,120,240 deg, +z); idx 3-5 =
    # bottom triangle (eclipsed, -z); idx 6-8 = caps on the 3 rectangular faces (60,180,300 deg, z=0).
    _tri_t = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), 0.7] for k in range(3)]
    _tri_b = [[math.cos(math.radians(120 * k)), math.sin(math.radians(120 * k)), -0.7] for k in range(3)]
    _caps = [[1.3 * math.cos(math.radians(60 + 120 * k)), 1.3 * math.sin(math.radians(60 + 120 * k)), 0.0] for k in range(3)]
    R[("CN9", "TTP-9 tricapped trigonal prism")] = _norm_rows(np.array(_tri_t + _tri_b + _caps, float))
    return {k: _norm_rows(v) for k, v in R.items()}


REFS = _ref_polyhedra()


# --- CN10 polyhedra for NON-f-block metals (Mission A2, 2026-06-05). -------
# When ``DELFIN_FFFREE_CN10_POLYHEDRA=1`` (or PURE_TRACK3=1) AND CN==10 AND the
# metal is NOT in the f-block, the assembler picks among three Wells-canonical
# CN10 polyhedra: BICAP-10 (D4d, bicapped square antiprism), CSAP-10 (C4v,
# capped square antiprism with an extra cap = monocapped SAP-9 capped twice),
# and SAP-10 (D5d, pentagonal antiprism, the "true" 10-vertex antiprism).
# Default OFF, byte-identical to HEAD when unset (non-f-block CN10 still falls
# back to legacy).  BICAP-10 reuses the f_block_polyhedra builder (one source).
def _bicap_10_unit() -> np.ndarray:
    """BICAP-10 vertex unit-vectors — re-exported from f_block_polyhedra so
    the legacy ``REFS`` table and the non-f-block CN10 path share ONE source.
    Vertex order: 0,1 = axial caps (+z,-z); 2-5 = upper square; 6-9 = lower
    square (45°-rotated)."""
    from delfin.fffree import f_block_polyhedra as _FBP
    return _FBP.bicap_10_vertices()


def _csap_10_unit() -> np.ndarray:
    """CSAP-10 — Capped Square Antiprism (with an extra cap on a square face).

    A 10-vertex C4v variant: one axial cap (+z) plus the SAP-8 square antiprism
    plus one EXTRA equatorial-belt vertex.  In the Wells / IUPAC catalogue this
    is also called the "monocapped square antiprism with a belt cap".  Vertex
    order: 0 = apical axial cap (+z); 1-4 = upper square (closer to cap);
    5-8 = lower square (45°-rotated); 9 = equatorial belt cap (φ=22.5°, z=0).
    """
    h_top = 0.55
    h_bot = 0.70
    r_top = math.sqrt(max(0.0, 1.0 - h_top * h_top))
    r_bot = math.sqrt(max(0.0, 1.0 - h_bot * h_bot))
    cap_axial = [(0.0, 0.0, 1.0)]
    top = [(r_top * math.cos(math.radians(90.0 * k)),
            r_top * math.sin(math.radians(90.0 * k)), h_top) for k in range(4)]
    bot = [(r_bot * math.cos(math.radians(45.0 + 90.0 * k)),
            r_bot * math.sin(math.radians(45.0 + 90.0 * k)), -h_bot)
           for k in range(4)]
    belt = [(math.cos(math.radians(22.5)), math.sin(math.radians(22.5)), 0.0)]
    return _norm_rows(np.asarray(cap_axial + top + bot + belt, dtype=float))


def _sap_10_unit() -> np.ndarray:
    """SAP-10 — Pentagonal Antiprism (D5d).

    The "true" 10-vertex antiprism: two parallel regular pentagons, the lower
    rotated by 36° relative to the upper.  Vertex order: 0-4 = upper pentagon
    (φ=0,72,144,216,288 at +z); 5-9 = lower pentagon (φ=36,108,180,252,324 at
    -z).  For the regular antiprism whose lateral edges equal the pentagon
    edges, h = sin(π/10) (in units where the in-plane radius is 1).  Vectors
    are then unit-normalised so all 10 vertices sit on the unit sphere.
    """
    h = math.sin(math.pi / 10.0)
    top = [(math.cos(math.radians(72.0 * k)),
            math.sin(math.radians(72.0 * k)), h) for k in range(5)]
    bot = [(math.cos(math.radians(36.0 + 72.0 * k)),
            math.sin(math.radians(36.0 + 72.0 * k)), -h) for k in range(5)]
    return _norm_rows(np.asarray(top + bot, dtype=float))


# Canonical CN10 geometry name → vertex builder.  Pólya-complete enumeration
# picks among these three; default is BICAP-10 (most common in CN10 X-ray).
CN10_VERTICES = {
    "BICAP-10 bicapped square antiprism": _bicap_10_unit,
    "CSAP-10 capped square antiprism": _csap_10_unit,
    "SAP-10 pentagonal antiprism": _sap_10_unit,
}

CN10_ALIASES = {
    "BICAP-10": "BICAP-10 bicapped square antiprism",
    "CSAP-10": "CSAP-10 capped square antiprism",
    "SAP-10": "SAP-10 pentagonal antiprism",
}


def cn10_polyhedra_enabled() -> bool:
    """True when the non-f-block CN10 dispatch is enabled (env-gated)."""
    import os as _os
    return (_os.environ.get("DELFIN_FFFREE_CN10_POLYHEDRA", "0") == "1"
            or _os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1")


GEOM_BY_CN = {
    2: ["L-2 linear"],                              # Phase G: Cu(I)/Ag(I)/Au(I)/Hg(II)
    3: ["SP-3 trigonal planar", "T-3 T-shape"],
    4: ["T-4 tetrahedron", "SP-4 square planar"],
    5: ["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"],
    6: ["OC-6 octahedron", "TPR-6 trigonal prism"],
    7: ["PB-7 pentagonal bipyramid"],
    8: ["SQAP-8 square antiprism"],
    9: ["TTP-9 tricapped trigonal prism"],
}


def geometries_for_cn(cn: int, metal: str = "") -> list:
    """CN → list of valid geometry names, optionally metal-aware.

    Phase C extension (Task #64): when ``metal`` is an f-block element AND
    ``DELFIN_FFFREE_FBLOCK_CN8_12=1`` (or PURE_TRACK3=1), CN 8-12 expand to
    include the dedicated f-block polyhedra (SAP-8, DD-8, TTP-9, CSAP-9,
    BICAP-10, CAP-11, IH-12) from :mod:`delfin.fffree.f_block_polyhedra`.

    Mission A2 (2026-06-05): when ``cn == 10`` AND the metal is NOT f-block
    AND ``DELFIN_FFFREE_CN10_POLYHEDRA=1`` (or PURE_TRACK3=1), the list also
    grows to include BICAP-10 / CSAP-10 / SAP-10 (Wells canonical CN10
    polyhedra) — so non-f-block CN10 (e.g. Y, early-d, large-anion) can build
    a real polyhedron instead of falling through to legacy.

    Default OFF — byte-identical to ``GEOM_BY_CN[cn]`` when unset.
    """
    import os as _os
    base = list(GEOM_BY_CN.get(cn, []))
    fb_on = (_os.environ.get("DELFIN_FFFREE_FBLOCK_CN8_12", "0") == "1"
             or _os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1")
    cn10_on = cn10_polyhedra_enabled()
    # CN10 non-f-block extension (Mission A2): inject 3 polyhedra here when
    # metal is non-f-block (or unknown).  F-block CN10 keeps its own dispatch
    # below, gated by the independent FBLOCK_CN8_12 flag.
    if cn == 10 and cn10_on:
        try:
            from delfin.fffree import f_block_polyhedra as _FBP_check
            is_fb = bool(metal) and _FBP_check.is_f_block(metal)
        except ImportError:
            is_fb = False
        if not is_fb:
            for g in CN10_VERTICES.keys():
                if g not in base:
                    base.append(g)
            return base
    if not fb_on:
        return base
    if not metal:
        return base
    try:
        from delfin.fffree import f_block_polyhedra as _FBP
    except ImportError:
        return base
    if not _FBP.is_f_block(metal):
        return base
    extra = _FBP.FBLOCK_GEOM_BY_CN.get(cn, [])
    # Append extras NOT already present (lex-deterministic).
    for g in extra:
        if g not in base:
            base.append(g)
    return base

# covalent radii (metal subset + donors); M-D = r(M) + r(D)
COV = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39, "Se": 1.20, "As": 1.19,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "La": 2.07,
}


def ref_vectors(geometry: str) -> np.ndarray:
    for (cn, shape), v in REFS.items():
        if shape == geometry:
            return v
    # Mission A2 (2026-06-05): non-f-block CN10 polyhedra dispatch.  Read-only
    # name lookup — callers only reach these geometries when env-gated.
    if geometry in CN10_VERTICES:
        return CN10_VERTICES[geometry]()
    if geometry in CN10_ALIASES:
        return CN10_VERTICES[CN10_ALIASES[geometry]]()
    # Phase C f-block CN8-12 dispatch (Task #64, 2026-06-04).  Env-gated
    # (``DELFIN_FFFREE_FBLOCK_CN8_12=1`` or ``DELFIN_FFFREE_PURE_TRACK3=1``)
    # — when unset, the new geometries are NEVER reached because callers
    # only pick them via ``decompose._default_geometry`` which itself
    # honours the env-gate.  We still attempt the f-block dispatch on the
    # error path so explicitly-named SAP-8/DD-8/TTP-9/… calls work even
    # outside the gate (read-only; geometry strings come from the caller).
    try:
        from delfin.fffree import f_block_polyhedra as _FBP
        return _FBP.ref_vectors_fblock(geometry)
    except (KeyError, ImportError):
        pass
    raise KeyError(geometry)


# Phase C: Shannon-radii M-D fallback for f-block metals.  When the metal is
# Ln/An AND the (metal, donor) pair is in the Shannon table, prefer that over
# the covalent-radii sum (Ln/An ionic bonding is poorly described by COV).
_FBLOCK_FALLBACK = None


def _get_fblock_fallback():
    global _FBLOCK_FALLBACK
    if _FBLOCK_FALLBACK is None:
        try:
            from delfin.fffree import f_block_polyhedra as _FBP
            _FBLOCK_FALLBACK = _FBP
        except ImportError:
            _FBLOCK_FALLBACK = False
    return _FBLOCK_FALLBACK if _FBLOCK_FALLBACK is not False else None


def md_distance(metal: str, donor: str) -> float:
    # Phase C: f-block metals get Shannon-radii distances (env-gated; default
    # path is byte-identical because the gate also controls whether f-block
    # CN8-12 geometries are reached at all).
    import os as _os
    if (_os.environ.get("DELFIN_FFFREE_FBLOCK_CN8_12", "0") == "1"
            or _os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"):
        _FBP = _get_fblock_fallback()
        if _FBP is not None:
            r = _FBP.md_distance_fblock(metal, donor)
            if r is not None:
                return r
    return COV.get(metal, 1.5) + COV.get(donor, 0.75)


def _kabsch_resid(P: np.ndarray, Q: np.ndarray) -> float:
    """Min mean-squared deviation aligning P onto Q by a proper rotation
    (Kabsch, determinant-corrected to forbid reflection)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = P @ R.T - Q
    return float((diff * diff).sum() / len(P))


def cshm(observed_vecs, geometry: str) -> float:
    """Continuous shape measure (0 = ideal, larger = worse) of observed donor
    unit-vectors against the ideal ``geometry`` polyhedron — scale-, rotation-
    and permutation-invariant (Kabsch over all vertex permutations, S =
    100·min_resid/obs_var).  Self-contained, deterministic; returns 0.0 for
    degenerate / size-mismatched input.  Used by the self-gate to reject
    catastrophic-coordination-shape outlier builds (#39)."""
    import itertools
    try:
        Q = ref_vectors(geometry)
    except KeyError:
        return 0.0
    P = np.asarray(observed_vecs, dtype=float)
    n = len(P)
    if n < 2 or n != len(Q):
        return 0.0

    def _unit_rms(a):
        rms = math.sqrt(float((a * a).sum()) / len(a))
        return a / rms if rms > 1e-9 else a

    P = _unit_rms(P)
    Q = _unit_rms(Q)
    obs_var = float((P * P).sum() / n)
    if obs_var < 1e-9:
        return 0.0
    best = float("inf")
    for perm in itertools.permutations(range(n)):
        r = _kabsch_resid(P, Q[list(perm)])
        if r < best:
            best = r
    return 100.0 * best / obs_var
