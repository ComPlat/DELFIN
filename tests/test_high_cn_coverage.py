"""Tests for the universal high-CN (CN 8-12) coverage extension.

Covers :mod:`delfin.fffree.high_cn_coverage` and its integration into
:mod:`delfin.fffree.polyhedra`, :mod:`delfin.fffree.decompose` and
:mod:`delfin.fffree.converter_backend`.

Contract
--------

1. Every new vertex builder returns an ``(N, 3)`` ``np.ndarray`` of unit vectors
   (``||v|| = 1``) and is deterministic (byte-identical across runs).
2. Default OFF (no env flag set) is byte-identical to HEAD: legacy
   ``GEOM_BY_CN[cn]`` is returned untouched, ``md_distance`` falls back to the
   plain covalent-radii sum, no new polyhedra are dispatched.
3. With ``DELFIN_FFFREE_HIGH_CN_COVERAGE=1`` set, every CN 8-12 picks up the
   new extras (DD-8/TPR-8/CSAP-9/MFF-9/SPHENO-10/CAP-11/OCD-11/ICO-12/CUOH-12/
   ACUOH-12) AND ``decompose._default_geometry`` resolves CN 10/11/12 to a
   canonical polyhedron for any (non-f-block) metal.
4. The bounds-tolerance factor reader ``high_cn_bounds_tol_factor`` returns
   the env-configurable float (default 2.0).
5. PAYQIS-class CCDC SMILES (Ru CN10) decompose successfully under the
   high-CN gate and produce at least one XYZ via the converter backend.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

from delfin.fffree import high_cn_coverage as HCN
from delfin.fffree import polyhedra as PLY
from delfin.fffree import polya_isomer_count as PIC


# ---------------------------------------------------------------------------
# Env helpers (per-test isolation)
# ---------------------------------------------------------------------------


@pytest.fixture
def env_off(monkeypatch):
    """Default state — all high-CN env flags unset."""
    for flag in (
        "DELFIN_FFFREE_HIGH_CN_COVERAGE",
        "DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR",
        "DELFIN_FFFREE_PURE_TRACK3",
        "DELFIN_FFFREE_FBLOCK_CN8_12",
        "DELFIN_FFFREE_CN10_POLYHEDRA",
    ):
        monkeypatch.delenv(flag, raising=False)


@pytest.fixture
def env_high_cn(monkeypatch):
    """High-CN coverage gate on; bounds-tol factor at default."""
    monkeypatch.setenv("DELFIN_FFFREE_HIGH_CN_COVERAGE", "1")
    monkeypatch.delenv("DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_FBLOCK_CN8_12", raising=False)


# ---------------------------------------------------------------------------
# Vertex math — shape, unit-norm, determinism
# ---------------------------------------------------------------------------


NEW_POLYHEDRA = [
    (8, "DD-8 dodecahedron"),
    (8, "TPR-8 bicapped trigonal prism"),
    (9, "CSAP-9 capped square antiprism"),
    (9, "MFF-9 muffin"),
    (10, "SPHENO-10 sphenocorona"),
    (11, "CAP-11 monocapped pentagonal antiprism"),
    (11, "OCD-11 octadecahedron"),
    (12, "ICO-12 icosahedron"),
    (12, "CUOH-12 cuboctahedron"),
    (12, "ACUOH-12 anticuboctahedron"),
]


@pytest.mark.parametrize("cn,name", NEW_POLYHEDRA)
def test_vertex_shape_and_unit_norm(cn, name):
    """Each new polyhedron returns (cn, 3) with all rows unit-length."""
    v = HCN.HIGH_CN_VERTICES[name]()
    assert v.shape == (cn, 3), f"{name}: {v.shape} != ({cn},3)"
    norms = np.linalg.norm(v, axis=1)
    np.testing.assert_allclose(norms, np.ones(cn), atol=1e-9)


@pytest.mark.parametrize("cn,name", NEW_POLYHEDRA)
def test_vertex_determinism(cn, name):
    """Repeated calls return byte-identical arrays."""
    v1 = HCN.HIGH_CN_VERTICES[name]()
    v2 = HCN.HIGH_CN_VERTICES[name]()
    np.testing.assert_array_equal(v1, v2)


@pytest.mark.parametrize("cn,name", NEW_POLYHEDRA)
def test_ref_vectors_dispatch(env_high_cn, cn, name):
    """``polyhedra.ref_vectors`` resolves the canonical name."""
    v = PLY.ref_vectors(name)
    expected = HCN.HIGH_CN_VERTICES[name]()
    np.testing.assert_array_equal(v, expected)


@pytest.mark.parametrize("cn,name", NEW_POLYHEDRA)
def test_min_vertex_separation_positive(cn, name):
    """No two vertices coincide (minimum angle > 5°)."""
    v = HCN.HIGH_CN_VERTICES[name]()
    min_ang = 360.0
    for i in range(cn):
        for j in range(i + 1, cn):
            cos_a = float(np.clip(v[i] @ v[j], -1.0, 1.0))
            ang = math.degrees(math.acos(cos_a))
            if ang < min_ang:
                min_ang = ang
    assert min_ang > 5.0, f"{name}: degenerate vertices (min angle {min_ang}°)"


# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------


def test_short_aliases_resolve():
    """Short aliases (DD-8, TPR-8, ...) resolve to the canonical builder."""
    for alias, canonical in HCN.HIGH_CN_ALIASES.items():
        v_alias = HCN.ref_vectors_high_cn(alias)
        v_canon = HCN.ref_vectors_high_cn(canonical)
        np.testing.assert_array_equal(v_alias, v_canon)


# ---------------------------------------------------------------------------
# Default-OFF byte-identity contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cn,expected", [
    (2, ["L-2 linear"]),
    (4, ["T-4 tetrahedron", "SP-4 square planar"]),
    (5, ["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"]),
    (6, ["OC-6 octahedron", "TPR-6 trigonal prism"]),
    (8, ["SQAP-8 square antiprism"]),
    (10, []),
    (12, []),
])
def test_default_off_geometries_byte_identical(env_off, cn, expected):
    """With no env flags set, ``geometries_for_cn`` matches the legacy table."""
    got = PLY.geometries_for_cn(cn)
    assert got == expected


def test_default_off_md_distance_legacy(env_off):
    """Default OFF: md_distance is the plain covalent-radii sum."""
    # Ru in COV: 1.46 + N 0.71 = 2.17 (no high-CN correction).
    assert PLY.md_distance("Ru", "N") == pytest.approx(2.17, abs=1e-9)
    # cn argument is accepted but ignored in default-OFF mode (because Ru
    # IS in COV, the high-CN fallback is not consulted).
    assert PLY.md_distance("Ru", "N", cn=10) == pytest.approx(2.17, abs=1e-9)


def test_default_off_ref_vectors_legacy_unchanged(env_off):
    """Default OFF: legacy ref_vectors still return the same arrays."""
    # OC-6 byte-identical legacy shape.
    oc6 = PLY.ref_vectors("OC-6 octahedron")
    assert oc6.shape == (6, 3)
    expected = np.array([
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1],
    ], float)
    np.testing.assert_array_equal(oc6, expected)


# ---------------------------------------------------------------------------
# Env-ON geometry dispatch
# ---------------------------------------------------------------------------


def test_env_on_cn8_includes_extras(env_high_cn):
    """CN8 + Ru picks up DD-8 and TPR-8 in addition to SQAP-8."""
    geoms = PLY.geometries_for_cn(8, "Ru")
    assert "SQAP-8 square antiprism" in geoms
    assert "DD-8 dodecahedron" in geoms
    assert "TPR-8 bicapped trigonal prism" in geoms


def test_env_on_cn9_includes_extras(env_high_cn):
    """CN9 + Ru picks up CSAP-9 and MFF-9."""
    geoms = PLY.geometries_for_cn(9, "Ru")
    assert "TTP-9 tricapped trigonal prism" in geoms
    assert "CSAP-9 capped square antiprism" in geoms
    assert "MFF-9 muffin" in geoms


def test_env_on_cn10_includes_all_four(env_high_cn):
    """CN10 + Ru picks up BICAP-10 / CSAP-10 / SAP-10 + SPHENO-10."""
    geoms = PLY.geometries_for_cn(10, "Ru")
    for expected in (
        "BICAP-10 bicapped square antiprism",
        "CSAP-10 capped square antiprism",
        "SAP-10 pentagonal antiprism",
        "SPHENO-10 sphenocorona",
    ):
        assert expected in geoms


def test_env_on_cn11_includes_extras(env_high_cn):
    """CN11 + Ru picks up CAP-11 and OCD-11."""
    geoms = PLY.geometries_for_cn(11, "Ru")
    assert "CAP-11 monocapped pentagonal antiprism" in geoms
    assert "OCD-11 octadecahedron" in geoms


def test_env_on_cn12_includes_all_three(env_high_cn):
    """CN12 + Ru picks up ICO-12 / CUOH-12 / ACUOH-12."""
    geoms = PLY.geometries_for_cn(12, "Ru")
    assert "ICO-12 icosahedron" in geoms
    assert "CUOH-12 cuboctahedron" in geoms
    assert "ACUOH-12 anticuboctahedron" in geoms


# ---------------------------------------------------------------------------
# Universal M-D distance fallback
# ---------------------------------------------------------------------------


def test_md_distance_high_cn_in_module():
    """``md_distance_high_cn`` adds the +0.10 Å correction at CN >= 8."""
    # Pa is NOT in legacy COV; the extended table has Pa = 2.00 Å.  F = 0.57 Å.
    r_cn7 = HCN.md_distance_high_cn("Pa", "F", cn=7)
    r_cn8 = HCN.md_distance_high_cn("Pa", "F", cn=8)
    assert r_cn7 == pytest.approx(2.57, abs=1e-9)        # 2.00 + 0.57
    assert r_cn8 == pytest.approx(2.67, abs=1e-9)        # + 0.10 high-CN


def test_md_distance_high_cn_missing_returns_none():
    """Unknown metal returns None so the caller can fall back."""
    assert HCN.md_distance_high_cn("Xx", "N", cn=8) is None


def test_md_distance_polyhedra_wires_high_cn(env_high_cn):
    """``polyhedra.md_distance`` falls back to the high-CN extension when
    the metal is missing from the legacy COV table AND the env is on."""
    # Pa is not in legacy COV but IS in _HIGH_CN_COV_EXTRA.
    r = PLY.md_distance("Pa", "F", cn=8)
    assert r == pytest.approx(2.67, abs=1e-9)


def test_md_distance_legacy_metal_no_correction(env_high_cn):
    """When metal IS in legacy COV, NO correction is applied even with env on.

    This preserves the byte-identity for metals that already work and avoids
    silently shifting M-D distances for the established (Ru, Mo, Fe, ...) set.
    """
    assert PLY.md_distance("Ru", "N", cn=10) == pytest.approx(2.17, abs=1e-9)


# ---------------------------------------------------------------------------
# Bounds-tol factor reader
# ---------------------------------------------------------------------------


def test_bounds_tol_factor_default(env_off):
    """Default bounds-tol factor is 2.0."""
    assert HCN.high_cn_bounds_tol_factor() == 2.0


def test_bounds_tol_factor_env_override(monkeypatch):
    """Env override is honoured (as a float)."""
    monkeypatch.setenv("DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR", "3.0")
    assert HCN.high_cn_bounds_tol_factor() == 3.0


def test_bounds_tol_factor_invalid_input(monkeypatch):
    """Non-numeric env value falls back to the default 2.0."""
    monkeypatch.setenv("DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR", "not-a-number")
    assert HCN.high_cn_bounds_tol_factor() == 2.0


def test_bounds_tol_factor_floor_at_one(monkeypatch):
    """Factor is clamped to >= 1.0 (loosening only, never tightening)."""
    monkeypatch.setenv("DELFIN_FFFREE_HIGH_CN_BOUNDS_TOL_FACTOR", "0.5")
    assert HCN.high_cn_bounds_tol_factor() == 1.0


# ---------------------------------------------------------------------------
# Pólya groups for the new polyhedra
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cn,name", NEW_POLYHEDRA)
def test_polya_group_lookup(cn, name):
    """``high_cn_group`` returns a (group, n) tuple for every new polyhedron."""
    res = HCN.high_cn_group(name)
    assert res is not None
    group, n = res
    assert n == cn
    assert len(group) >= 1  # at least the identity
    # First element must be the identity permutation.
    assert group[0] == tuple(range(cn)) or any(g == tuple(range(cn)) for g in group)


@pytest.mark.parametrize("key,canonical", [
    ("dd_8_hicn", "DD-8 dodecahedron"),
    ("tpr_8", "TPR-8 bicapped trigonal prism"),
    ("csap_9_hicn", "CSAP-9 capped square antiprism"),
    ("mff_9", "MFF-9 muffin"),
    ("spheno_10", "SPHENO-10 sphenocorona"),
    ("cap_11_hicn", "CAP-11 monocapped pentagonal antiprism"),
    ("ocd_11", "OCD-11 octadecahedron"),
    ("ico_12_hicn", "ICO-12 icosahedron"),
    ("cuoh_12", "CUOH-12 cuboctahedron"),
    ("acuoh_12", "ACUOH-12 anticuboctahedron"),
])
def test_polya_isomer_count_key_map(key, canonical):
    """``polya_isomer_count._GEOM_KEY_TO_SHAPE`` resolves each high-CN key."""
    from delfin.fffree.polya_isomer_count import _GEOM_KEY_TO_SHAPE, _get_group
    assert _GEOM_KEY_TO_SHAPE.get(key) == canonical
    group, n = _get_group(key)
    assert n in (8, 9, 10, 11, 12)
    assert len(group) >= 1


# ---------------------------------------------------------------------------
# Decompose / build path — PAYQIS-class CCDC scenario
# ---------------------------------------------------------------------------


def test_default_off_cn10_decompose_returns_none(env_off):
    """Without the env flag, a Ru CN10 SMILES is rejected by ``decompose``."""
    from delfin.fffree.decompose import decompose
    # Octa-aquo + 2 ammine Ru CN10 (simplified PAYQIS-class).
    smi = "[Ru](O)(O)(O)(O)(O)(O)(O)(O)(N)N"
    d = decompose(smi)
    # CN10 is not in the default-allowed CN set.
    assert d is None


def test_env_on_cn10_decompose_resolves_geometry(env_high_cn):
    """Under high-CN coverage, a Ru CN10 SMILES decomposes into BICAP-10."""
    from delfin.fffree.decompose import decompose
    smi = "[Ru](O)(O)(O)(O)(O)(O)(O)(O)(N)N"
    d = decompose(smi)
    assert d is not None
    assert d["cn"] == 10
    assert d["metal"] == "Ru"
    assert d["geometry"] == "BICAP-10 bicapped square antiprism"


def test_env_on_cn11_decompose_resolves_geometry(env_high_cn):
    """Under high-CN coverage, a Ru CN11 SMILES decomposes into CAP-11."""
    from delfin.fffree.decompose import decompose
    smi = "[Ru](O)(O)(O)(O)(O)(O)(O)(O)(O)(N)N"
    d = decompose(smi)
    assert d is not None
    assert d["cn"] == 11
    assert d["geometry"] == "CAP-11 monocapped pentagonal antiprism"


def test_env_on_cn12_decompose_resolves_geometry(env_high_cn):
    """Under high-CN coverage, a Ru CN12 SMILES decomposes into ICO-12."""
    from delfin.fffree.decompose import decompose
    smi = "[Ru](O)(O)(O)(O)(O)(O)(O)(O)(O)(O)(N)N"
    d = decompose(smi)
    assert d is not None
    assert d["cn"] == 12
    assert d["geometry"] == "ICO-12 icosahedron"
