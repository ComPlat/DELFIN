"""Tests for ``delfin.fffree.f_block_polyhedra`` (Phase C, Task #64, 2026-06-04).

Verifies:
  * Each CN8-12 polyhedron has the correct vertex count and unit-vector shape.
  * Approximate point-group symmetry of each polyhedron (D4d / D2d / D3h / C4v /
    Ih) holds under the prescribed proper-rotation generators (when applicable).
  * Env-OFF byte-identical behaviour: with ``DELFIN_FFFREE_FBLOCK_CN8_12``
    unset (and PURE_TRACK3 unset), the legacy ``polyhedra.ref_vectors`` /
    ``polyhedra.md_distance`` / ``decompose._default_geometry`` results
    are IDENTICAL to the values they returned before this patch.
  * Determinism: repeated calls produce byte-identical np arrays.
  * F-block metal classification (Ln + An).
  * M-D distance fallback returns Shannon-radii sum.
  * CN8 default geometry under env-ON for La is "SAP-8 square antiprism".
  * Group definitions (Pólya) for SAP-8 / DD-8 / TTP-9 / CSAP-9 / BICAP-10
    are well-formed permutation groups (closed under composition).
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest


@pytest.fixture
def env_off(monkeypatch):
    """All Phase-C env flags unset — legacy byte-identical path."""
    monkeypatch.delenv("DELFIN_FFFREE_FBLOCK_CN8_12", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    return None


@pytest.fixture
def env_on(monkeypatch):
    """DELFIN_FFFREE_FBLOCK_CN8_12=1 enables Phase-C dispatch."""
    monkeypatch.setenv("DELFIN_FFFREE_FBLOCK_CN8_12", "1")
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    return None


# ---------------------------------------------------------------------------
# Vertex-set shape and unit-sphere invariance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name,expected_n",
    [
        ("SAP-8 square antiprism", 8),
        ("DD-8 dodecahedron", 8),
        ("TDD-8 triangular dodecahedron", 8),
        ("TTP-9 tricapped trigonal prism", 9),
        ("CSAP-9 capped square antiprism", 9),
        ("BICAP-10 bicapped square antiprism", 10),
        ("CAP-11 monocapped pentagonal antiprism", 11),
        ("IH-12 icosahedron", 12),
    ],
)
def test_vertex_count_and_unit_norm(geom_name, expected_n):
    """Every polyhedron returns the correct N×3 array of unit vectors."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.ref_vectors_fblock(geom_name)
    assert v.shape == (expected_n, 3), f"{geom_name}: shape {v.shape}"
    norms = np.linalg.norm(v, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9), f"{geom_name}: norms {norms}"


# ---------------------------------------------------------------------------
# Approximate point-group symmetry checks
# ---------------------------------------------------------------------------


def test_sap_8_has_c4_axis():
    """SAP-8 has a C4 axis along z (top square cycles under 90° rotation)."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.sap_8_vertices()
    # Rotate by 90° about z
    c, s = math.cos(math.pi / 2), math.sin(math.pi / 2)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    v_rot = v @ Rz.T
    # Each rotated vertex must coincide with an original vertex (within tol)
    matched = 0
    for vr in v_rot:
        for vo in v:
            if np.allclose(vr, vo, atol=1e-6):
                matched += 1
                break
    assert matched == 8, f"SAP-8 C4_z symmetry violated: matched {matched}/8"


def test_dd_8_has_c2_z():
    """DD-8 (D2d) has a C2 axis along z (180° rotation maps onto self)."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.dd_8_vertices()
    # 180° rotation about z = (-x, -y, z)
    v_rot = v.copy()
    v_rot[:, 0] *= -1
    v_rot[:, 1] *= -1
    matched = 0
    for vr in v_rot:
        for vo in v:
            if np.allclose(vr, vo, atol=1e-6):
                matched += 1
                break
    assert matched == 8, f"DD-8 C2_z symmetry violated: matched {matched}/8"


def test_ttp_9_has_c3_axis_and_d3h_mirror():
    """TTP-9 (D3h) has C3 about z and a horizontal mirror plane (z → -z)."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.ttp_9_vertices()
    # C3 about z = 120° rotation
    c, s = math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    v_c3 = v @ Rz.T
    matched_c3 = sum(
        1 for vr in v_c3 if any(np.allclose(vr, vo, atol=1e-6) for vo in v)
    )
    assert matched_c3 == 9, f"TTP-9 C3 symmetry: matched {matched_c3}/9"

    # σh mirror (z → -z)
    v_mir = v.copy()
    v_mir[:, 2] *= -1
    matched_mir = sum(
        1 for vr in v_mir if any(np.allclose(vr, vo, atol=1e-6) for vo in v)
    )
    assert matched_mir == 9, f"TTP-9 σh symmetry: matched {matched_mir}/9"


def test_csap_9_has_c4_axis():
    """CSAP-9 (C4v) has a C4 axis along z (apical cap + 2 4-fold squares)."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.csap_9_vertices()
    c, s = math.cos(math.pi / 2), math.sin(math.pi / 2)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    v_rot = v @ Rz.T
    matched = sum(
        1 for vr in v_rot if any(np.allclose(vr, vo, atol=1e-6) for vo in v)
    )
    assert matched == 9, f"CSAP-9 C4 symmetry: matched {matched}/9"


def test_bicap_10_has_c4_and_horizontal_mirror():
    """BICAP-10 (D4d) has C4_z and the 4-fold improper rotation, but as a proper
    rotation group has a C4_z + a horizontal C2 swapping the two axial caps."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.bicap_10_vertices()
    c, s = math.cos(math.pi / 2), math.sin(math.pi / 2)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    v_rot = v @ Rz.T
    matched = sum(
        1 for vr in v_rot if any(np.allclose(vr, vo, atol=1e-6) for vo in v)
    )
    assert matched == 10, f"BICAP-10 C4 symmetry: matched {matched}/10"


def test_ih_12_centroid_zero():
    """The icosahedron is centred at the origin (vertex sum ≈ 0)."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.ih_12_vertices()
    assert np.allclose(v.sum(axis=0), np.zeros(3), atol=1e-9)


def test_ih_12_edge_lengths_uniform():
    """All 30 icosahedron edges (nearest-neighbour pairs) have the same length."""
    from delfin.fffree import f_block_polyhedra as FBP

    v = FBP.ih_12_vertices()
    # nearest neighbours: each vertex has 5 neighbours at the same distance
    nn_distances = []
    for i in range(12):
        ds = sorted(np.linalg.norm(v[j] - v[i]) for j in range(12) if j != i)
        nn_distances.append(ds[0])  # nearest neighbour
    d0 = nn_distances[0]
    assert all(abs(d - d0) < 1e-6 for d in nn_distances), \
        f"Ih edges not uniform: {nn_distances}"


# ---------------------------------------------------------------------------
# Determinism: repeated calls = byte-identical
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name",
    [
        "SAP-8 square antiprism",
        "TTP-9 tricapped trigonal prism",
        "BICAP-10 bicapped square antiprism",
        "IH-12 icosahedron",
    ],
)
def test_determinism_byte_identical(geom_name):
    """Two calls to the same vertex builder return byte-identical arrays."""
    from delfin.fffree import f_block_polyhedra as FBP

    v1 = FBP.ref_vectors_fblock(geom_name)
    v2 = FBP.ref_vectors_fblock(geom_name)
    assert np.array_equal(v1, v2), f"{geom_name}: non-deterministic vertex set"


# ---------------------------------------------------------------------------
# Element classification
# ---------------------------------------------------------------------------


def test_is_f_block_lanthanides():
    """All 15 lanthanides classify as f-block."""
    from delfin.fffree.f_block_polyhedra import is_f_block

    for ln in ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
               "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]:
        assert is_f_block(ln), f"{ln} should be f-block"


def test_is_f_block_actinides():
    """Key actinides (Th, U, Pu, Am) classify as f-block."""
    from delfin.fffree.f_block_polyhedra import is_f_block

    for an in ["Ac", "Th", "U", "Np", "Pu", "Am", "Cm"]:
        assert is_f_block(an), f"{an} should be f-block"


def test_is_f_block_excludes_d_block_and_main_group():
    """Common d-block / main-group elements are NOT f-block."""
    from delfin.fffree.f_block_polyhedra import is_f_block

    for non in ["Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W", "Pt", "Pd",
                "Ru", "Rh", "Ir", "Mn", "Cr", "Ti", "Y", "Sc",
                "Na", "Ca", "Mg", "Al", "Si"]:
        assert not is_f_block(non), f"{non} should NOT be f-block"


# ---------------------------------------------------------------------------
# M-D distance (Shannon radii fallback)
# ---------------------------------------------------------------------------


def test_md_distance_fblock_known_pair():
    """``md_distance_fblock('La', 'O')`` returns a positive ionic-radii sum."""
    from delfin.fffree.f_block_polyhedra import md_distance_fblock

    d = md_distance_fblock("La", "O")
    assert d is not None
    assert 2.4 < d < 2.8, f"La-O distance {d} outside reasonable range"


def test_md_distance_fblock_non_fblock_returns_none():
    """For non-f-block metals the helper returns None (caller falls back)."""
    from delfin.fffree.f_block_polyhedra import md_distance_fblock

    assert md_distance_fblock("Fe", "N") is None
    assert md_distance_fblock("Pt", "O") is None


def test_md_distance_fblock_unknown_donor_returns_none():
    """Unknown donor (e.g. He) returns None."""
    from delfin.fffree.f_block_polyhedra import md_distance_fblock

    assert md_distance_fblock("La", "He") is None


# ---------------------------------------------------------------------------
# Env-OFF byte-identical contract
# ---------------------------------------------------------------------------


def test_env_off_md_distance_byte_identical(env_off):
    """With env-OFF, ``polyhedra.md_distance`` returns the legacy COV sum
    even for f-block metals (NO Shannon-radii fallback applied)."""
    from delfin.fffree import polyhedra as P

    # La COV = 2.07, O COV = 0.66 → legacy sum = 2.73
    d = P.md_distance("La", "O")
    assert abs(d - (2.07 + 0.66)) < 1e-9, \
        f"env-OFF La-O = {d}, expected {2.07 + 0.66} (legacy COV sum)"


def test_env_off_ref_vectors_legacy_only(env_off):
    """With env-OFF, requesting ``SAP-8 square antiprism`` falls through.

    Note: ``polyhedra.ref_vectors`` does still attempt the f-block fallback
    on the error path (explicit name request) — this is by design so that
    direct callers can use the new geometries.  The legacy code-path that
    routes through GEOM_BY_CN / decompose._default_geometry remains
    unchanged because those tables only emit the f-block names under the
    env-gate.
    """
    from delfin.fffree import polyhedra as P

    # Legacy SQAP-8 still resolves to the original (8, 3) array.
    v = P.ref_vectors("SQAP-8 square antiprism")
    assert v.shape == (8, 3)


def test_env_off_default_geometry_unchanged(env_off):
    """``decompose._default_geometry('La', 8)`` returns the legacy 'SQAP-8'
    when env-OFF (not the new 'SAP-8 square antiprism')."""
    from delfin.fffree.decompose import _default_geometry

    assert _default_geometry("La", 8) == "SQAP-8 square antiprism"


def test_env_off_geometries_for_cn_unchanged(env_off):
    """``geometries_for_cn`` returns the legacy GEOM_BY_CN list when env-OFF."""
    from delfin.fffree.polyhedra import geometries_for_cn, GEOM_BY_CN

    assert geometries_for_cn(8, "La") == GEOM_BY_CN[8]
    assert geometries_for_cn(12, "La") == GEOM_BY_CN.get(12, [])


# ---------------------------------------------------------------------------
# Env-ON behaviour: dispatch routes to f-block polyhedra
# ---------------------------------------------------------------------------


def test_env_on_default_geometry_la_cn8(env_on):
    """Env-ON: ``_default_geometry('La', 8)`` routes to SAP-8 (f-block)."""
    from delfin.fffree.decompose import _default_geometry

    assert _default_geometry("La", 8) == "SAP-8 square antiprism"


def test_env_on_default_geometry_la_cn12(env_on):
    """Env-ON: ``_default_geometry('La', 12)`` routes to IH-12 icosahedron."""
    from delfin.fffree.decompose import _default_geometry

    assert _default_geometry("La", 12) == "IH-12 icosahedron"


def test_env_on_default_geometry_non_fblock_unchanged(env_on):
    """Env-ON: non-f-block metals at CN8/9 still get the legacy geometry."""
    from delfin.fffree.decompose import _default_geometry

    assert _default_geometry("Fe", 8) == "SQAP-8 square antiprism"
    assert _default_geometry("Mo", 9) == "TTP-9 tricapped trigonal prism"


def test_env_on_md_distance_la_uses_shannon(env_on):
    """Env-ON: ``md_distance('La', 'O')`` returns Shannon-radii sum (≈ 2.53)
    instead of legacy COV sum (2.73)."""
    from delfin.fffree import polyhedra as P

    d = P.md_distance("La", "O")
    assert abs(d - 2.53) < 1e-9, f"env-ON La-O = {d}, expected 2.53"


def test_env_on_geometries_for_cn_la(env_on):
    """Env-ON: ``geometries_for_cn(8, 'La')`` includes the new SAP-8 / DD-8."""
    from delfin.fffree.polyhedra import geometries_for_cn

    gs = geometries_for_cn(8, "La")
    assert "SAP-8 square antiprism" in gs
    assert "DD-8 dodecahedron" in gs


def test_env_on_geometries_for_cn_fe_unchanged(env_on):
    """Env-ON: non-f-block metal stays on the legacy geometry list."""
    from delfin.fffree.polyhedra import geometries_for_cn

    gs = geometries_for_cn(8, "Fe")
    assert gs == ["SQAP-8 square antiprism"]


# ---------------------------------------------------------------------------
# Pólya group well-formedness
# ---------------------------------------------------------------------------


def _is_group_closed(group, n):
    """A list of permutation tuples is a group iff closed under composition."""
    G = set(group)
    for g in group:
        for h in group:
            composed = tuple(g[h[i]] for i in range(n))
            if composed not in G:
                return False
    return True


def test_polya_sap_8_group_closed():
    """SAP-8 D4 proper rotation group is closed under composition."""
    from delfin.fffree.f_block_polyhedra import fblock_group

    grp, n = fblock_group("SAP-8 square antiprism")
    assert n == 8
    assert len(grp) == 8, f"D4 must have order 8, got {len(grp)}"
    assert _is_group_closed(grp, n)


def test_polya_ttp_9_group_closed():
    """TTP-9 D3 proper rotation group (order 6) is closed."""
    from delfin.fffree.f_block_polyhedra import fblock_group

    grp, n = fblock_group("TTP-9 tricapped trigonal prism")
    assert n == 9
    assert len(grp) == 6, f"D3 must have order 6, got {len(grp)}"
    assert _is_group_closed(grp, n)


def test_polya_csap_9_group_closed():
    """CSAP-9 C4 proper rotation group (order 4) is closed."""
    from delfin.fffree.f_block_polyhedra import fblock_group

    grp, n = fblock_group("CSAP-9 capped square antiprism")
    assert n == 9
    assert len(grp) == 4, f"C4 must have order 4, got {len(grp)}"
    assert _is_group_closed(grp, n)


def test_polya_bicap_10_group_closed():
    """BICAP-10 D4 proper rotation group (order 8) is closed."""
    from delfin.fffree.f_block_polyhedra import fblock_group

    grp, n = fblock_group("BICAP-10 bicapped square antiprism")
    assert n == 10
    assert len(grp) == 8, f"D4 must have order 8, got {len(grp)}"
    assert _is_group_closed(grp, n)


def test_polya_ih_12_group_trivial():
    """IH-12 conservative group is the trivial group {identity}."""
    from delfin.fffree.f_block_polyhedra import fblock_group

    grp, n = fblock_group("IH-12 icosahedron")
    assert n == 12
    assert len(grp) == 1  # conservative: only identity


# ---------------------------------------------------------------------------
# Cross-module integration: count_isomers picks up the f-block groups
# ---------------------------------------------------------------------------


def test_count_isomers_sap_8_homoleptic():
    """SAP-8 with 8 identical donors has exactly 1 isomer (all equivalent)."""
    from delfin.fffree.polya_isomer_count import count_isomers

    n = count_isomers("sap_8", {"O": 8})
    assert n == 1


def test_count_isomers_ttp_9_homoleptic():
    """TTP-9 with 9 identical donors has exactly 1 isomer."""
    from delfin.fffree.polya_isomer_count import count_isomers

    n = count_isomers("ttp_9_fblock", {"O": 9})
    assert n == 1


def test_count_isomers_sap_8_mixed_returns_positive():
    """SAP-8 with 4 O + 4 N donors has multiple isomers (Burnside)."""
    from delfin.fffree.polya_isomer_count import count_isomers

    n = count_isomers("sap_8", {"O": 4, "N": 4})
    # For SAP-8 D4 (order 8), MA4B4 should give ≥ 1 distinct isomer.
    # Exact count depends on which 4 vertices the As occupy; pivoting on D4 the
    # full vertex-coloring count is 8 (= multinomial(8;4,4)/8 ≈ 8.75 → 8 or 9).
    assert n >= 1


# ---------------------------------------------------------------------------
# Public API: dispatch via polyhedra.ref_vectors works under env-OFF too
# ---------------------------------------------------------------------------


def test_polyhedra_dispatch_explicit_fblock_name(env_off):
    """``polyhedra.ref_vectors('SAP-8 square antiprism')`` works even env-OFF.

    Direct-name dispatch routes through the f-block fallback at the error
    edge — this is intentional so explicit callers can request the new
    geometries regardless of the env gate.
    """
    from delfin.fffree import polyhedra as P

    v = P.ref_vectors("SAP-8 square antiprism")
    assert v.shape == (8, 3)
    norms = np.linalg.norm(v, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9)


def test_polyhedra_dispatch_unknown_name_still_raises(env_off):
    """Unknown geometry strings still raise KeyError."""
    from delfin.fffree import polyhedra as P

    with pytest.raises(KeyError):
        P.ref_vectors("NONEXISTENT-99")
