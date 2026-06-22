"""Tests for the high-CN lanthanide / actinide polyhedron extension.

Adds vector tables, point-group labels and a metal-aware classifier for the
CN=8-12 polyhedra that lanthanides (La-Lu), high-CN actinides (Th, U, ...)
and other large ionic metals (Y, Sc, Sr, Ba, Ca, Pb, Bi) tend to adopt.

New polyhedra
-------------

* CN=8  ``bicapped_trig_antiprism`` (D3d)
* CN=9  ``capped_sap``              (C4v)
* CN=10 ``pentag_antiprism``        (D5d)
* CN=10 ``sphenocorona``            (C2v, Johnson J87, irregular)
* CN=11 ``mono_capped_pentag_aprism`` (C5v)
* CN=12 ``cuboctahedron``           (Oh)

Acceptance contract
-------------------

1. Every new vector table has shape ``(cn, 3)`` and consists of unit vectors.
2. Aliases (full name + short code) all resolve back to the same key.
3. Point-group lookup returns the correct Schoenflies symbol.
4. CN=8 La complex with default env: returns the legacy ``sq_antiprism``;
   with ``DELFIN_HIGH_CN_POLYHEDRA=1``: returns ``bicapped_trig_antiprism``.
5. CN=9 La complex with env=1: returns ``capped_sap``.
6. CN=10 La complex with predominantly small donors (C/F/H) and env=1:
   returns ``pentag_antiprism``; mixed/O-dominant donors stay on
   ``bicapped_sap``.
7. CN=12 La complex with env=1: returns ``cuboctahedron`` (not icosahedron).
8. Non-lanthanide metals (e.g. Ru, Mo, Fe) at CN=8 with env=1: NOT routed
   to the lanthanide overrides — they keep the legacy default.
9. Symmetry check — bicapped trigonal antiprism unit-vector sum is ≈ 0.
10. Default OFF — every (cn, donor_types) combination is bit-exact identical
    to :func:`classify_geometry_from_cn_donors` when the env flag is unset.
11. Malformed input — ``cn`` outside 8-12 or ``metal_sym=None`` falls back
    safely to the metal-agnostic classifier without raising.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from delfin._polyhedron_targets import (
    _CN_GEOM_KEY,
    _HIGH_CN_ACTINIDES,
    _HIGH_CN_IONIC_OTHER,
    _LANTHANIDES,
    _high_cn_lanthanide_geometry,
    _is_high_cn_ionic_metal,
    classify_geometry_from_cn_donors,
    classify_geometry_from_cn_donors_with_metal,
    get_ideal_donor_vectors,
    get_target_point_group,
)


# ---------------------------------------------------------------------------
# Env isolation
# ---------------------------------------------------------------------------


@pytest.fixture
def env_off(monkeypatch):
    """Default — DELFIN_HIGH_CN_POLYHEDRA unset.  Bit-exact legacy behaviour."""
    monkeypatch.delenv("DELFIN_HIGH_CN_POLYHEDRA", raising=False)
    monkeypatch.delenv("DELFIN_LINEAR_CN2", raising=False)
    return None


@pytest.fixture
def env_on(monkeypatch):
    """DELFIN_HIGH_CN_POLYHEDRA=1 — lanthanide overrides enabled."""
    monkeypatch.setenv("DELFIN_HIGH_CN_POLYHEDRA", "1")
    monkeypatch.delenv("DELFIN_LINEAR_CN2", raising=False)
    return None


# ---------------------------------------------------------------------------
# Vector tables — shape, unit-norm, alias resolution, point group
# ---------------------------------------------------------------------------


NEW_TABLES = [
    (8, "bicapped_trig_antiprism", "BCTA", "D3d"),
    (9, "capped_sap", "CSAP", "C4v"),
    (10, "pentag_antiprism", "PAP", "D5d"),
    (10, "sphenocorona", "J87", "C2v"),
    (11, "mono_capped_pentag_aprism", "MCPAP", "C5v"),
    (12, "cuboctahedron", "COCT", "Oh"),
]


@pytest.mark.parametrize("cn,key,short,pg", NEW_TABLES)
def test_new_polyhedron_table_shape_and_unit_norm(cn, key, short, pg):
    """Every new polyhedron returns ``(cn, 3)`` and only unit vectors."""
    vecs = get_ideal_donor_vectors(cn, key)
    assert vecs.shape == (cn, 3), f"{key}: expected ({cn},3), got {vecs.shape}"
    norms = np.linalg.norm(vecs, axis=1)
    np.testing.assert_allclose(norms, np.ones(cn), atol=1e-9)


@pytest.mark.parametrize("cn,key,short,pg", NEW_TABLES)
def test_short_alias_resolves(cn, key, short, pg):
    """Short alias (e.g. ``CSAP``) resolves to the same vectors."""
    full = get_ideal_donor_vectors(cn, key)
    by_short = get_ideal_donor_vectors(cn, short)
    np.testing.assert_allclose(full, by_short)


@pytest.mark.parametrize("cn,key,short,pg", NEW_TABLES)
def test_point_group_lookup(cn, key, short, pg):
    """Point-group label matches the expected Schoenflies symbol."""
    assert get_target_point_group(cn, key) == pg
    assert get_target_point_group(cn, short) == pg


def test_bicapped_trig_antiprism_centroid_is_origin():
    """8 vectors of bicapped trigonal antiprism sum to ~0 (centroid at metal)."""
    vecs = get_ideal_donor_vectors(8, "bicapped_trig_antiprism")
    centroid = vecs.sum(axis=0)
    np.testing.assert_allclose(centroid, np.zeros(3), atol=1e-9)


def test_tricapped_tp_centroid_is_origin():
    """Reference: 9 vectors of D3h TTP also sum to ~0."""
    vecs = get_ideal_donor_vectors(9, "tricapped_tp")
    centroid = vecs.sum(axis=0)
    np.testing.assert_allclose(centroid, np.zeros(3), atol=1e-9)


def test_pentag_antiprism_centroid_is_origin():
    """D5d pentagonal antiprism: 10 unit vectors must sum to ~0."""
    vecs = get_ideal_donor_vectors(10, "pentag_antiprism")
    centroid = vecs.sum(axis=0)
    np.testing.assert_allclose(centroid, np.zeros(3), atol=1e-9)


def test_cuboctahedron_centroid_is_origin():
    """Oh cuboctahedron: 12 unit vertices sum to 0."""
    vecs = get_ideal_donor_vectors(12, "cuboctahedron")
    centroid = vecs.sum(axis=0)
    np.testing.assert_allclose(centroid, np.zeros(3), atol=1e-9)


def test_cuboctahedron_vs_icosahedron_distinct():
    """Cuboctahedron must differ from icosahedron — different ideal CN=12."""
    cubo = get_ideal_donor_vectors(12, "cuboctahedron")
    ico = get_ideal_donor_vectors(12, "icosahedron")
    # arrays are different (cannot be element-wise close)
    assert not np.allclose(cubo, ico)


# ---------------------------------------------------------------------------
# Metal-group taxonomy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("sym", sorted(_LANTHANIDES))
def test_lanthanides_are_high_cn_ionic(sym):
    """All 15 lanthanides La..Lu count as high-CN ionic metals."""
    assert _is_high_cn_ionic_metal(sym), f"{sym} should be ionic"


@pytest.mark.parametrize("sym", sorted(_HIGH_CN_ACTINIDES))
def test_high_cn_actinides_are_classified(sym):
    """Th, U, Pa, ... routed to high-CN preference set."""
    assert _is_high_cn_ionic_metal(sym), f"{sym} should be ionic"


@pytest.mark.parametrize("sym", ["Y", "Sc", "Sr", "Ba", "Ca", "Pb", "Bi"])
def test_large_ionic_metals_classified(sym):
    """Y/Sc plus large alkaline-earth + Pb/Bi count as high-CN ionic."""
    assert _is_high_cn_ionic_metal(sym)


@pytest.mark.parametrize("sym", ["Fe", "Cu", "Ru", "Os", "Mo", "Ni", "Co", "Mn", "Ti"])
def test_d_block_metals_are_not_ionic(sym):
    """Late-d and early-d transition metals are NOT routed to ionic overrides."""
    assert not _is_high_cn_ionic_metal(sym), f"{sym} should not be ionic"


def test_empty_or_none_metal_is_not_ionic():
    """``None`` / empty string is safe — not ionic, no override."""
    assert _is_high_cn_ionic_metal(None) is False
    assert _is_high_cn_ionic_metal("") is False
    assert _is_high_cn_ionic_metal("   ") is False


# ---------------------------------------------------------------------------
# Override helper directly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cn,expected", [
    (8, "bicapped_trig_antiprism"),
    (9, "capped_sap"),
    (11, "mono_capped_pentag_aprism"),
    (12, "cuboctahedron"),
])
def test_high_cn_override_per_cn(cn, expected):
    """Override returns the lanthanide-preferred geometry for CN=8/9/11/12."""
    geom, iso = _high_cn_lanthanide_geometry(cn, ["O"] * cn)
    assert geom == expected
    assert iso is None


def test_high_cn_override_cn10_o_dominant_stays_bcsap():
    """CN=10 with all-O donors: lanthanide nitrato stays on bicapped SAP."""
    geom, _ = _high_cn_lanthanide_geometry(10, ["O"] * 10)
    assert geom == "bicapped_sap"


def test_high_cn_override_cn10_small_donors_picks_pap():
    """CN=10 with predominantly small (F/C/H) donors → pentagonal antiprism."""
    geom, _ = _high_cn_lanthanide_geometry(10, ["F"] * 8 + ["O", "O"])
    assert geom == "pentag_antiprism"


def test_high_cn_override_outside_range_returns_none():
    """CN<8 returns ``None`` — caller falls back to metal-agnostic classifier."""
    assert _high_cn_lanthanide_geometry(6, ["N"] * 6) is None
    assert _high_cn_lanthanide_geometry(7, ["O"] * 7) is None


# ---------------------------------------------------------------------------
# Metal-aware classifier — env-flag gating
# ---------------------------------------------------------------------------


def test_default_off_cn8_la_is_legacy_sap(env_off):
    """Env off: La CN=8 returns the metal-agnostic legacy default."""
    legacy = classify_geometry_from_cn_donors(8, ["O"] * 8)
    head = classify_geometry_from_cn_donors_with_metal(
        8, ["O"] * 8, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == legacy == ("sq_antiprism", None)


def test_env_on_cn8_la_routes_to_bcta(env_on):
    """Env on + La: CN=8 routed to bicapped trigonal antiprism."""
    head = classify_geometry_from_cn_donors_with_metal(
        8, ["O"] * 8, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("bicapped_trig_antiprism", None)


def test_env_on_cn9_la_aquo_routes_to_csap(env_on):
    """Env on + La aquo (CN=9, all O): capped square antiprism."""
    head = classify_geometry_from_cn_donors_with_metal(
        9, ["O"] * 9, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("capped_sap", None)


def test_env_on_cn10_la_nitrato_stays_bcsap(env_on):
    """Env on + La with 10× O donors (nitrato): bicapped SAP (not PAP)."""
    head = classify_geometry_from_cn_donors_with_metal(
        10, ["O"] * 10, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("bicapped_sap", None)


def test_env_on_cn10_la_fluoride_routes_to_pap(env_on):
    """Env on + La with F-dominant CN=10: pentagonal antiprism."""
    head = classify_geometry_from_cn_donors_with_metal(
        10, ["F"] * 10, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("pentag_antiprism", None)


def test_env_on_cn12_la_routes_to_cuboctahedron(env_on):
    """Env on + La CN=12 (hexa-bidentate nitrato): cuboctahedron, not Ih."""
    head = classify_geometry_from_cn_donors_with_metal(
        12, ["O"] * 12, metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("cuboctahedron", None)


def test_env_on_d_block_metal_keeps_legacy(env_on):
    """Env on + Ru (d-block, NOT ionic) at CN=8: legacy ``sq_antiprism``."""
    head = classify_geometry_from_cn_donors_with_metal(
        8, ["O"] * 8, metal_symbol="Ru", metal_formal_charge=2,
    )
    assert head == ("sq_antiprism", None)


def test_env_on_none_metal_keeps_legacy(env_on):
    """Env on but metal=None: no override → legacy classifier wins."""
    head = classify_geometry_from_cn_donors_with_metal(
        9, ["O"] * 9, metal_symbol=None, metal_formal_charge=None,
    )
    assert head == ("tricapped_tp", None)


def test_env_off_is_bit_exact_legacy(env_off):
    """Env off → every CN and metal pair equals the legacy classifier."""
    test_cases = [
        # (cn, donor_types, metal_sym)
        (4, ["N", "N", "N", "N"], "Fe"),
        (6, ["O"] * 6, "Ru"),
        (8, ["O"] * 8, "La"),
        (8, ["F"] * 8, "U"),
        (9, ["O"] * 9, "La"),
        (10, ["O"] * 10, "Ce"),
        (11, ["O"] * 11, "Nd"),
        (12, ["O"] * 12, "La"),
        (12, ["O"] * 12, "Fe"),
    ]
    for cn, donors, metal in test_cases:
        legacy = classify_geometry_from_cn_donors(cn, donors)
        head = classify_geometry_from_cn_donors_with_metal(
            cn, donors, metal_symbol=metal, metal_formal_charge=3,
        )
        assert head == legacy, (
            f"env-off mismatch cn={cn} metal={metal}: "
            f"legacy={legacy} head={head}"
        )


def test_env_on_actinide_uranium_cn8_routes_to_bcta(env_on):
    """U at CN=8: bicapped trigonal antiprism (UF8-type)."""
    head = classify_geometry_from_cn_donors_with_metal(
        8, ["F"] * 8, metal_symbol="U", metal_formal_charge=4,
    )
    assert head == ("bicapped_trig_antiprism", None)


# ---------------------------------------------------------------------------
# Malformed / edge-case input
# ---------------------------------------------------------------------------


def test_malformed_empty_donor_list_safe(env_on):
    """Empty donor list with high CN: still returns a valid (geom, iso) tuple."""
    head = classify_geometry_from_cn_donors_with_metal(
        8, [], metal_symbol="La", metal_formal_charge=3,
    )
    assert head == ("bicapped_trig_antiprism", None)


def test_malformed_unknown_cn_safe(env_on):
    """CN=13 with La + env on: falls through to legacy ``undefined_cn13``."""
    head = classify_geometry_from_cn_donors_with_metal(
        13, ["O"] * 13, metal_symbol="La", metal_formal_charge=3,
    )
    assert head[0].startswith("undefined")


def test_malformed_unknown_metal_safe(env_on):
    """Unknown metal symbol: no override, legacy classifier used."""
    head = classify_geometry_from_cn_donors_with_metal(
        8, ["O"] * 8, metal_symbol="Xx", metal_formal_charge=0,
    )
    assert head == ("sq_antiprism", None)


def test_legacy_classifier_cn11_returns_mcpap():
    """The metal-agnostic classifier also returns MCPAP for CN=11 (no env)."""
    geom, iso = classify_geometry_from_cn_donors(11, ["O"] * 11)
    assert geom == "mono_capped_pentag_aprism"
    assert iso is None


def test_alias_table_includes_new_entries():
    """Sanity: alias map has every new entry."""
    for key in [
        "bicapped_trig_antiprism",
        "capped_sap",
        "pentag_antiprism",
        "sphenocorona",
        "mono_capped_pentag_aprism",
        "cuboctahedron",
    ]:
        # Find at least one (cn, alias) pair whose value is `key`.
        assert any(v == key for v in _CN_GEOM_KEY.values()), (
            f"alias map missing {key}"
        )
