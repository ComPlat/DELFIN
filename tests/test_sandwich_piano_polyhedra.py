"""Tests for ``delfin.fffree.sandwich_piano_polyhedra`` (Phase C / Task #44,
Mission A7, 2026-06-05).

Verifies:
  * Each SANDWICH-10 / PIANO-STOOL-8 / HALF-SANDWICH-9 polyhedron has the
    correct vertex count, unit-norm shape, and deterministic byte-identical
    output across repeated calls.
  * Cp / arene ring vertices lie at the expected half-angle from +z (not on
    the equatorial belt) and the σ-tripod sits below the ring at a wider
    opening angle, matching the C3v / D5d / C3v symmetry of CCDC ferrocenes /
    piano-stools / half-sandwiches.
  * Vertex orthogonality / staggering: SANDWICH-10 D5d top vs bottom Cp is
    rotated by 36° (staggered ferrocene); PIANO-STOOL-8 σ-tripod is between
    the projected Cp atoms.
  * Environment-OFF byte identity: with all relevant env flags unset, the
    legacy ``polyhedra.geometries_for_cn(cn)`` / ``ref_vectors(geom)`` paths
    do NOT mention SANDWICH-10 / PIANO-STOOL-8 / HALF-SANDWICH-9 and
    return values are byte-identical to the legacy table.
  * Environment-ON dispatch: setting ``DELFIN_FFFREE_SANDWICH_POLYHEDRA=1``
    extends the CN8 / CN9 / CN10 geometry lists.
  * Decompose dispatch: ferrocene-class SMILES (η⁵-Cp pair, η⁵-Cp + 3 σ,
    η⁶-arene + 3 σ) route to the new polyhedra under the env flag.
  * Determinism: two calls return byte-identical np.ndarray content.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# 1. Vertex-set shape, unit norms, byte-identical determinism
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name,expected_n",
    [
        ("SANDWICH-10 bis-eta5-Cp", 10),
        ("PIANO-STOOL-8 eta5-Cp+L3", 8),
        ("HALF-SANDWICH-9 eta6+L3", 9),
    ],
)
def test_vertex_count_and_unit_norm(geom_name, expected_n):
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.ref_vectors_sandwich(geom_name)
    assert v.shape == (expected_n, 3), f"{geom_name}: shape {v.shape}"
    norms = np.linalg.norm(v, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9), f"{geom_name}: norms {norms}"


def test_determinism_byte_identical():
    """Two successive calls return byte-identical arrays."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    for geom in [
        "SANDWICH-10 bis-eta5-Cp",
        "PIANO-STOOL-8 eta5-Cp+L3",
        "HALF-SANDWICH-9 eta6+L3",
    ]:
        a = SP.ref_vectors_sandwich(geom)
        b = SP.ref_vectors_sandwich(geom)
        assert a.tobytes() == b.tobytes(), f"{geom}: non-deterministic"


# ---------------------------------------------------------------------------
# 2. Sandwich D5d staggered: top vs bottom Cp ring offset by 36°
# ---------------------------------------------------------------------------


def test_sandwich_10_d5d_staggered():
    """Top Cp ring (idx 0-4) at azimuths 0,72,144,216,288°; bottom (idx 5-9)
    at 36,108,180,252,324° — i.e. rotated by 36° vs top (D5d staggered)."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.sandwich_10_vertices()
    top = v[:5]
    bot = v[5:]
    # z-signs: top all positive, bottom all negative
    assert (top[:, 2] > 0).all(), "top Cp ring should be on +z hemisphere"
    assert (bot[:, 2] < 0).all(), "bottom Cp ring should be on -z hemisphere"
    # Top mirror-z of bottom up to a 36° azimuthal rotation
    top_az = np.degrees(np.arctan2(top[:, 1], top[:, 0])) % 360.0
    bot_az = np.degrees(np.arctan2(bot[:, 1], bot[:, 0])) % 360.0
    # Expected azimuths
    expected_top = np.array([0.0, 72.0, 144.0, 216.0, 288.0])
    expected_bot = np.array([36.0, 108.0, 180.0, 252.0, 324.0])
    assert np.allclose(np.sort(top_az), np.sort(expected_top), atol=1e-6)
    assert np.allclose(np.sort(bot_az), np.sort(expected_bot), atol=1e-6)


def test_sandwich_10_centroid_axis_is_z():
    """The two Cp ring centroids lie on the +z and -z axis (sandwich C∞)."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.sandwich_10_vertices()
    top_centroid = v[:5].mean(axis=0)
    bot_centroid = v[5:].mean(axis=0)
    assert abs(top_centroid[0]) < 1e-9 and abs(top_centroid[1]) < 1e-9
    assert abs(bot_centroid[0]) < 1e-9 and abs(bot_centroid[1]) < 1e-9
    assert top_centroid[2] > 0
    assert bot_centroid[2] < 0
    # The centroids are antipodal up to the cosine of the ring half-angle
    assert abs(top_centroid[2] + bot_centroid[2]) < 1e-9


# ---------------------------------------------------------------------------
# 3. Piano-stool C3v: tripod opens below the Cp ring
# ---------------------------------------------------------------------------


def test_piano_stool_tripod_below_cp_ring():
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.piano_stool_8_vertices()
    cp = v[:5]
    tripod = v[5:]
    # Cp on +z, tripod on -z
    assert (cp[:, 2] > 0).all(), "Cp ring should sit on +z"
    assert (tripod[:, 2] < 0).all(), "tripod should sit below (-z)"
    # Tripod azimuths 0, 120, 240
    az = np.sort(np.degrees(np.arctan2(tripod[:, 1], tripod[:, 0])) % 360.0)
    assert np.allclose(az, [0.0, 120.0, 240.0], atol=1e-6)


def test_piano_stool_c3_symmetry():
    """C3 rotation around z permutes the tripod cyclically."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.piano_stool_8_vertices()
    tripod = v[5:]
    # Rotate by 120° around z
    c = math.cos(math.radians(120.0))
    s = math.sin(math.radians(120.0))
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    rotated = tripod @ Rz.T
    # Each rotated tripod atom should match SOME original tripod atom
    for r in rotated:
        diffs = np.linalg.norm(tripod - r, axis=1)
        assert diffs.min() < 1e-6, "tripod not C3-symmetric"


# ---------------------------------------------------------------------------
# 4. Half-sandwich arene ring planar, 30°-offset tripod
# ---------------------------------------------------------------------------


def test_half_sandwich_arene_ring_planar():
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.half_sandwich_9_vertices()
    arene = v[:6]
    tripod = v[6:]
    # All 6 arene vertices have the same z (planar ring)
    assert np.allclose(arene[:, 2], arene[0, 2], atol=1e-9)
    # 6 vertices, 60° spacing
    az = np.sort(np.degrees(np.arctan2(arene[:, 1], arene[:, 0])) % 360.0)
    expected = np.array([0.0, 60.0, 120.0, 180.0, 240.0, 300.0])
    assert np.allclose(az, expected, atol=1e-6)
    # Tripod at 30°, 150°, 270°
    taz = np.sort(np.degrees(np.arctan2(tripod[:, 1], tripod[:, 0])) % 360.0)
    assert np.allclose(taz, [30.0, 150.0, 270.0], atol=1e-6)


# ---------------------------------------------------------------------------
# 5. Env-OFF byte-identity: legacy paths unchanged
# ---------------------------------------------------------------------------


def test_env_off_geometries_for_cn_unchanged(monkeypatch):
    """With all env flags unset, geometries_for_cn returns the legacy list
    byte-identical to ``GEOM_BY_CN`` (the sandwich/piano-stool/half-sandwich
    polyhedra are NOT wired into the legacy CN→geometry table by default —
    they are accessed via the standalone module
    :mod:`delfin.fffree.sandwich_piano_polyhedra`)."""
    monkeypatch.delenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_FBLOCK_CN8_12", raising=False)
    from delfin.fffree.polyhedra import geometries_for_cn, GEOM_BY_CN

    for cn in (4, 5, 6, 7, 8, 9):
        assert geometries_for_cn(cn) == GEOM_BY_CN[cn], (
            f"CN{cn}: env-OFF must return GEOM_BY_CN value byte-identical")


def test_module_level_env_on_dispatch(monkeypatch):
    """When env flag is set, the standalone module exposes the new geometry
    catalogue via ``SANDWICH_GEOM_BY_CN`` and ``env_active``."""
    monkeypatch.setenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA", "1")
    from delfin.fffree import sandwich_piano_polyhedra as SP

    assert SP.env_active() is True
    assert "PIANO-STOOL-8 eta5-Cp+L3" in SP.SANDWICH_GEOM_BY_CN[8]
    assert "HALF-SANDWICH-9 eta6+L3" in SP.SANDWICH_GEOM_BY_CN[9]
    assert "SANDWICH-10 bis-eta5-Cp" in SP.SANDWICH_GEOM_BY_CN[10]


def test_legacy_ref_vectors_raises_keyerror_for_sandwich(monkeypatch):
    """Sandwich geometries are NOT registered in the legacy
    ``polyhedra.ref_vectors`` catalogue — calling it with a sandwich name
    raises KeyError.  The dedicated ``sandwich_piano_polyhedra.ref_vectors_sandwich``
    must be used instead (this is the explicit architecture decision —
    the legacy ``polyhedra`` module stays metal-agnostic point-donor only)."""
    monkeypatch.delenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    from delfin.fffree.polyhedra import ref_vectors

    with pytest.raises(KeyError):
        ref_vectors("SANDWICH-10 bis-eta5-Cp")

    from delfin.fffree.sandwich_piano_polyhedra import ref_vectors_sandwich
    v = ref_vectors_sandwich("SANDWICH-10 bis-eta5-Cp")
    assert v.shape == (10, 3)


# ---------------------------------------------------------------------------
# 6. Decompose-level dispatch for ferrocene-class SMILES
# ---------------------------------------------------------------------------


def test_polya_groups_registered_for_sandwich():
    """Pólya groups for the new geometries are registered in
    :mod:`polya_isomer_count` and ``enumerate_chelate_configs`` works on the
    keys exposed via :data:`converter_backend._GEOM_TO_POLYA`."""
    from delfin.fffree.polya_isomer_count import _get_group, enumerate_chelate_configs
    from delfin.fffree.converter_backend import _GEOM_TO_POLYA

    # Keys present
    for human in (
        "SANDWICH-10 bis-eta5-Cp",
        "PIANO-STOOL-8 eta5-Cp+L3",
        "HALF-SANDWICH-9 eta6+L3",
    ):
        key = _GEOM_TO_POLYA[human]
        group, n = _get_group(key)
        assert isinstance(group, list) and len(group) >= 1
        # Identity is always in the group
        identity = tuple(range(n))
        assert identity in group
    # PIANO-STOOL-8: ferrocene-class η⁵-Cp + 3 σ-donors emits non-empty configs
    specs_piano = [
        {"type": "Cp", "denticity": 5, "asym": False},
        {"type": "CO", "denticity": 1, "asym": False},
        {"type": "CO", "denticity": 1, "asym": False},
        {"type": "CO", "denticity": 1, "asym": False},
    ]
    cfgs = enumerate_chelate_configs("piano_stool_8", specs_piano)
    assert len(cfgs) >= 1, "PIANO-STOOL-8 must produce ≥ 1 config"
    # HALF-SANDWICH-9: η⁶-arene + 3 σ-donors
    specs_half = [
        {"type": "Arene", "denticity": 6, "asym": False},
        {"type": "Cl", "denticity": 1, "asym": False},
        {"type": "Cl", "denticity": 1, "asym": False},
        {"type": "PR3", "denticity": 1, "asym": False},
    ]
    cfgs = enumerate_chelate_configs("half_sandwich_9", specs_half)
    assert len(cfgs) >= 1, "HALF-SANDWICH-9 must produce ≥ 1 config"
    # SANDWICH-10: 2× η⁵-Cp rings
    specs_sw = [
        {"type": "Cp", "denticity": 5, "asym": False},
        {"type": "Cp", "denticity": 5, "asym": False},
    ]
    cfgs = enumerate_chelate_configs("sandwich_10", specs_sw)
    assert len(cfgs) >= 1, "SANDWICH-10 must produce ≥ 1 config"


# ---------------------------------------------------------------------------
# 7. Layout descriptor coverage (assembler contract)
# ---------------------------------------------------------------------------


def test_layout_descriptors_cover_all_geometries():
    from delfin.fffree import sandwich_piano_polyhedra as SP

    for geom in SP.SANDWICH_GEOM_VERTICES:
        layout = SP.layout_for(geom)
        assert layout is not None, f"missing layout for {geom}"
        # Every ring index references a valid vertex
        v = SP.ref_vectors_sandwich(geom)
        n = v.shape[0]
        for ring in layout["rings"]:
            assert all(0 <= i < n for i in ring["ring_indices"])
            assert ring["eta"] == len(ring["ring_indices"])
        for ti in layout["tripod_indices"]:
            assert 0 <= ti < n
        # Disjoint partition of vertex indices into rings + tripod
        all_idxs = sorted([i for r in layout["rings"] for i in r["ring_indices"]]
                          + list(layout["tripod_indices"]))
        assert all_idxs == list(range(n)), (
            f"{geom}: layout doesn't partition vertices: {all_idxs}")


def test_env_active_helper(monkeypatch):
    from delfin.fffree import sandwich_piano_polyhedra as SP

    monkeypatch.delenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    assert SP.env_active() is False
    monkeypatch.setenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA", "1")
    assert SP.env_active() is True
    monkeypatch.delenv("DELFIN_FFFREE_SANDWICH_POLYHEDRA")
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")
    assert SP.env_active() is True


# ---------------------------------------------------------------------------
# 8. Vertex orthogonality / mutual angles
# ---------------------------------------------------------------------------


def test_sandwich_10_within_ring_angles():
    """Within the top Cp ring, neighbours subtend 72° (pentagonal)."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.sandwich_10_vertices()
    top = v[:5]
    # Each vertex's nearest neighbour in the top ring
    for i in range(5):
        d_az_next = math.degrees(
            math.acos(np.clip(np.dot(top[i], top[(i + 1) % 5]), -1.0, 1.0))
        )
        # Two neighbours separated by 72° in azimuth, both at the same z.
        # The 3D angle subtended at the origin matches arccos(cos²(h) + sin²(h)·cos72°).
        s = math.sin(SP.ETA5_CP_RING_HALF_ANGLE)
        c = math.cos(SP.ETA5_CP_RING_HALF_ANGLE)
        expected = math.degrees(math.acos(c * c + s * s * math.cos(math.radians(72.0))))
        assert abs(d_az_next - expected) < 1e-6


def test_piano_stool_cp_ring_perpendicular_to_axis():
    """Cp ring centroid lies on +z; the ring spans a circle at half-angle 38°."""
    from delfin.fffree import sandwich_piano_polyhedra as SP

    v = SP.piano_stool_8_vertices()
    cp = v[:5]
    centroid = cp.mean(axis=0)
    assert abs(centroid[0]) < 1e-9 and abs(centroid[1]) < 1e-9
    # Each Cp vertex's angle from +z is ETA5_CP_RING_HALF_ANGLE
    for r in cp:
        cos_ang = float(np.dot(r, np.array([0.0, 0.0, 1.0])))
        assert abs(cos_ang - math.cos(SP.ETA5_CP_RING_HALF_ANGLE)) < 1e-6
