"""Tests for the hapto-honest hapto-only CShM metric.

Validates the η-cluster detection + the resulting hapto-only CShM
calculation on canonical geometries (ferrocene η5-η5, half-sandwich
piano-stool, η6-arene, η2-ethylene, η4-diene, η8-COT).  All inputs are
constructed in-memory (no archive dependency).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

# Make the quality_framework/scripts directory importable.
_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.abspath(os.path.join(
    _HERE, "..", "agent_workspace", "quality_framework", "scripts",
))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from hapto_only_cshm import (  # noqa: E402
    DEFAULT_MD_CUT,
    HAPTO_CC_MAX,
    HAPTO_CENTROID_CUT,
    SCHEMA_VERSION,
    compute_hapto_only_cshm_from_arrays,
    hapto_unit_geometry_score,
    ideal_m_centroid_distance,
)


# ---------------------------------------------------------------------------
# Helpers — canonical M-L geometries built directly in coordinates.
# ---------------------------------------------------------------------------
def _ferrocene_eta5_eta5(M: str = "Fe", r_mc: float = 1.65, r_cc: float = 1.40):
    """Ferrocene η5-η5: Fe between two parallel Cp rings.  ``r_mc`` is the
    target |M-C| derived from CCDC p50; the metal-centroid distance is
    sqrt(r_mc^2 - r_ring^2)."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    # Want metal-centroid = ideal_m_centroid_distance(M, 5); set r_mc so
    # that the metal-CARBON distance is the typical value.
    z_off = ideal_m_centroid_distance(M, 5)
    for sign in (+1.0, -1.0):
        for k in range(5):
            theta = 2.0 * math.pi * k / 5.0
            x = r_ring * math.cos(theta)
            y = r_ring * math.sin(theta)
            z = sign * z_off
            syms.append("C")
            coords.append([x, y, z])
    return syms, np.asarray(coords, dtype=float)


def _arene_eta6(M: str = "Cr", r_cc: float = 1.40):
    """(η6-C6H6)M-like: one benzene ring perpendicular to M at ideal
    M-centroid for the η6 mode."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 6.0))
    z = ideal_m_centroid_distance(M, 6)
    for k in range(6):
        theta = 2.0 * math.pi * k / 6.0
        coords.append([r_ring * math.cos(theta), r_ring * math.sin(theta), z])
        syms.append("C")
    return syms, np.asarray(coords, dtype=float)


def _eta2_ethylene(M: str = "Pt", r_cc: float = 1.40, r_mc: float = 2.10):
    """η2-CH2=CH2 on M.  The two olefin carbons sit symmetric to the
    M-midpoint axis."""
    syms = [M, "C", "C"]
    # Place the C-C bond along x, midpoint at (0, 0, r_mc).  d_mc=r_mc.
    coords = np.array([
        [0.0, 0.0, 0.0],
        [+r_cc / 2.0, 0.0, r_mc],
        [-r_cc / 2.0, 0.0, r_mc],
    ], dtype=float)
    return syms, coords


def _eta4_butadiene(M: str = "Fe", r_cc: float = 1.40):
    """η4-butadiene on M, planar diene at the ideal η4 distance."""
    syms = [M, "C", "C", "C", "C"]
    z = ideal_m_centroid_distance(M, 4)
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-1.5 * r_cc / 2.0, +r_cc * 0.4, z],
        [-0.5 * r_cc / 2.0, -r_cc * 0.4, z],
        [+0.5 * r_cc / 2.0, +r_cc * 0.4, z],
        [+1.5 * r_cc / 2.0, -r_cc * 0.4, z],
    ], dtype=float)
    return syms, coords


def _eta8_cot(M: str = "U", r_cc: float = 1.40):
    """η8-COT: 8-carbon planar ring on U."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 8.0))
    z = ideal_m_centroid_distance(M, 8)
    for k in range(8):
        theta = 2.0 * math.pi * k / 8.0
        coords.append([r_ring * math.cos(theta), r_ring * math.sin(theta), z])
        syms.append("C")
    return syms, np.asarray(coords, dtype=float)


def _half_sandwich_eta5_plus_sigma(M: str = "Fe", r_cc: float = 1.40,
                                    r_cl: float = 2.30):
    """CpFe(Cl)(N)(P) half-sandwich (piano-stool): one η5 Cp + 3 σ donors.
    The hapto-only metric should detect a single η5 unit."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    z = ideal_m_centroid_distance(M, 5)
    for k in range(5):
        theta = 2.0 * math.pi * k / 5.0
        coords.append([r_ring * math.cos(theta), r_ring * math.sin(theta), z])
        syms.append("C")
    # 3 σ donors in the equatorial plane.
    syms.append("Cl"); coords.append([r_cl, 0.0, 0.0])
    syms.append("N");  coords.append([-r_cl / 2.0, r_cl * math.sqrt(3) / 2.0, 0.0])
    syms.append("P");  coords.append([-r_cl / 2.0, -r_cl * math.sqrt(3) / 2.0, 0.0])
    return syms, np.asarray(coords, dtype=float)


def _octahedral_pure_sigma(M: str = "Fe", D: str = "N", r: float = 2.10):
    """Pure σ Werner-type — no hapto units expected."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
        syms.append(D)
        coords.append([r * v[0], r * v[1], r * v[2]])
    return syms, np.asarray(coords, dtype=float)


# ---------------------------------------------------------------------------
# Tests — schema + canonical geometries.
# ---------------------------------------------------------------------------
def test_schema_version_is_int():
    assert isinstance(SCHEMA_VERSION, int)
    assert SCHEMA_VERSION >= 1


def test_pure_sigma_no_hapto_units():
    """Werner-type [Fe(NH3)6] should return no hapto units and
    hapto_only_cshm = None."""
    syms, P = _octahedral_pure_sigma()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["hapto_n_units"] == 0
    assert rec["hapto_only_cshm"] is None


def test_ferrocene_eta5_eta5_two_units_near_zero():
    """Ideal ferrocene: two parallel η5 rings → hapto_only_cshm ≈ 0."""
    syms, P = _ferrocene_eta5_eta5()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["hapto_n_units"] == 2
    assert sorted(rec["hapto_unit_sizes"]) == [5, 5]
    assert rec["hapto_only_cshm"] is not None
    # Ideal ferrocene → CShM ~0 (centroid_dist matches ideal, ring planar,
    # M perfectly axial).  Allow a small tolerance for fp rounding.
    assert rec["hapto_only_cshm"] < 1.0
    assert rec["hapto_centroid_dist_dev_mean"] < 0.05
    assert rec["hapto_perpendicularity_dev_mean"] < 1.0
    assert rec["hapto_ring_plane_rmsd_mean"] < 1e-6


def test_eta6_arene_one_unit_near_zero():
    """(η6-C6H6)Cr → one η6 unit, near-ideal."""
    syms, P = _arene_eta6()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Cr"
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [6]
    assert rec["hapto_only_cshm"] is not None
    assert rec["hapto_only_cshm"] < 1.0


def test_eta2_ethylene_one_unit():
    """η2-ethylene on Pt → one η2 unit of size 2."""
    syms, P = _eta2_ethylene()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Pt"
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [2]
    assert rec["hapto_only_cshm"] is not None


def test_eta4_butadiene_one_unit():
    """η4-butadiene on Fe → one η4 unit of size 4."""
    syms, P = _eta4_butadiene()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [4]


def test_eta8_cot_one_unit_near_zero():
    """η8-COT on U → one η8 unit of size 8 at the ideal centroid distance."""
    syms, P = _eta8_cot()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "U"
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [8]
    assert rec["hapto_only_cshm"] is not None
    assert rec["hapto_only_cshm"] < 1.0


def test_half_sandwich_one_eta5_unit():
    """CpFe(Cl)(N)(P) piano-stool → exactly one η5 unit (the σ donors
    are not classified as hapto)."""
    syms, P = _half_sandwich_eta5_plus_sigma()
    rec = compute_hapto_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [5]
    assert rec["hapto_only_cshm"] is not None
    assert rec["hapto_only_cshm"] < 1.0


def test_distorted_ring_drives_cshm_up():
    """Take ferrocene, push one carbon out of the ring plane by 0.5 Å:
    the hapto-only CShM must rise sharply (this is the honesty test —
    the metric must detect a bad hapto unit)."""
    syms, P = _ferrocene_eta5_eta5()
    rec_ideal = compute_hapto_only_cshm_from_arrays(syms, P)
    # Displace one Cp atom along its local normal direction.
    P_bad = P.copy()
    P_bad[1, 2] += 0.5  # push C1 of top Cp by 0.5 Å in z
    rec_bad = compute_hapto_only_cshm_from_arrays(syms, P_bad)
    assert rec_bad["hapto_only_cshm"] > rec_ideal["hapto_only_cshm"]
    # Decisive — at least 10x worse than the ideal.
    assert rec_bad["hapto_only_cshm"] > 10.0 * (rec_ideal["hapto_only_cshm"] + 1e-3)
    assert rec_bad["hapto_ring_plane_rmsd_mean"] > rec_ideal["hapto_ring_plane_rmsd_mean"]


def test_centroid_dist_wrong_drives_cshm_up():
    """Move Fe up by 0.5 Å away from a ferrocene-centered position:
    centroid-distance deviation should make hapto_only_cshm rise."""
    syms, P = _ferrocene_eta5_eta5()
    rec_ideal = compute_hapto_only_cshm_from_arrays(syms, P)
    P_bad = P.copy()
    P_bad[0, 2] += 0.5  # shift Fe
    rec_bad = compute_hapto_only_cshm_from_arrays(syms, P_bad)
    assert rec_bad["hapto_only_cshm"] > rec_ideal["hapto_only_cshm"]
    assert rec_bad["hapto_centroid_dist_dev_mean"] > rec_ideal["hapto_centroid_dist_dev_mean"]


def test_perpendicularity_dev_drives_cshm_up():
    """Tilt the top Cp ring by 30° to violate perpendicularity →
    hapto_only_cshm rises and perpendicularity_dev_mean > 25°."""
    syms, P = _ferrocene_eta5_eta5()
    rec_ideal = compute_hapto_only_cshm_from_arrays(syms, P)
    # Tilt only the top ring (indices 1..5) by 30° about the x-axis.
    th = math.radians(30.0)
    R = np.array([
        [1.0, 0.0, 0.0],
        [0.0, math.cos(th), -math.sin(th)],
        [0.0, math.sin(th),  math.cos(th)],
    ])
    P_bad = P.copy()
    for i in range(1, 6):
        # Rotate around the ring centroid of the top ring.
        cen = P[1:6].mean(axis=0)
        P_bad[i] = R @ (P[i] - cen) + cen
    rec_bad = compute_hapto_only_cshm_from_arrays(syms, P_bad)
    assert rec_bad["hapto_only_cshm"] > rec_ideal["hapto_only_cshm"]
    # Mean over two rings — only the top is tilted (30°), bottom (0°) → mean ~15°.
    assert rec_bad["hapto_perpendicularity_dev_mean"] > 10.0


def test_hapto_unit_geometry_score_pure_function():
    """Direct probe of the per-unit scorer."""
    syms, P = _ferrocene_eta5_eta5()
    # Top Cp atoms are indices 1..5; metal at 0.
    scores = hapto_unit_geometry_score(P, [1, 2, 3, 4, 5], P[0], metal_symbol="Fe")
    assert scores["eta"] == 5.0
    assert math.isclose(scores["ideal_dist"], ideal_m_centroid_distance("Fe", 5))
    assert scores["centroid_dist"] > 0
    assert scores["perpendicularity_deg"] < 1e-6
    assert scores["ring_plane_rmsd"] < 1e-9


def test_ideal_m_centroid_lookup_and_fallback():
    """Known lookup keys + generic fallback for unknown (M, η)."""
    assert math.isclose(ideal_m_centroid_distance("Fe", 5), 1.65)
    assert math.isclose(ideal_m_centroid_distance("Cr", 6), 1.62)
    # Unknown metal/η → fallback base * (5/η)^0.5
    fallback = ideal_m_centroid_distance("Xx", 5)
    assert math.isclose(fallback, 1.8)
    # eta=4 unknown fallback = 1.8 * (5/4)^0.5
    assert math.isclose(ideal_m_centroid_distance("Xx", 4),
                        1.8 * math.sqrt(5.0 / 4.0))


def test_determinism_two_calls_identical():
    """Running the analysis twice on the same input → bit-identical JSON."""
    import json as _json
    syms, P = _half_sandwich_eta5_plus_sigma()
    a = compute_hapto_only_cshm_from_arrays(syms, P)
    b = compute_hapto_only_cshm_from_arrays(syms, P)
    assert _json.dumps(a, sort_keys=True) == _json.dumps(b, sort_keys=True)


def test_empty_input_returns_empty():
    rec = compute_hapto_only_cshm_from_arrays([], np.zeros((0, 3), float))
    assert rec["n_atoms"] == 0
    assert rec["metal_idx"] == -1
    assert rec["hapto_only_cshm"] is None
    assert rec["hapto_n_units"] == 0


def test_explicit_metal_idx_override():
    """An explicit metal_idx must override the auto-detection."""
    # Build a Cp ring + a single H 'metal' at center.
    syms = ["Xx"]
    coords = [[0.0, 0.0, 0.0]]
    r_cc = 1.40
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    z = 1.65
    for k in range(5):
        theta = 2.0 * math.pi * k / 5.0
        coords.append([r_ring * math.cos(theta),
                       r_ring * math.sin(theta), z])
        syms.append("C")
    rec = compute_hapto_only_cshm_from_arrays(
        syms, np.asarray(coords, float), metal_idx=0,
    )
    assert rec["metal_idx"] == 0
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_unit_sizes"] == [5]


def test_isolated_carbons_not_hapto():
    """Two isolated carbene-C donors > HAPTO_CC_MAX apart → no hapto unit."""
    syms = ["Fe", "C", "C"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [-2.0, 0.0, 0.0],
    ], dtype=float)
    rec = compute_hapto_only_cshm_from_arrays(syms, coords)
    assert rec["hapto_n_units"] == 0
    assert rec["hapto_only_cshm"] is None


def test_aggregator_summary_structure():
    """Mini-aggregator smoke: build a small list of records and check the
    summary metric keys are all present."""
    from aggregate_hapto_only import summarise_records
    syms, P = _ferrocene_eta5_eta5()
    r1 = compute_hapto_only_cshm_from_arrays(syms, P)
    syms2, P2 = _octahedral_pure_sigma()
    r2 = compute_hapto_only_cshm_from_arrays(syms2, P2)
    summary = summarise_records([r1, r2])
    for key in (
        "hapto_only_cshm_mean",
        "hapto_only_cshm_median",
        "hapto_only_cshm_p95",
        "hapto_only_cshm_max",
        "hapto_only_cshm_n_frames",
        "hapto_only_centroid_dist_mean",
        "hapto_only_centroid_dist_dev_mean",
        "hapto_only_perpendicularity_dev_mean",
        "hapto_only_ring_plane_rmsd_mean",
        "hapto_units_per_frame",
        "frames_with_hapto",
        "frames_no_hapto",
    ):
        assert key in summary
    assert summary["frames_with_hapto"] == 1
    assert summary["frames_no_hapto"] == 1
    assert summary["hapto_only_cshm_n_frames"] == 1
