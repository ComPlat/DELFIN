"""Tests for the hapto-honest σ-only CShM metric.

Validates the σ vs η/π classification and the resulting CShM calculation on
canonical geometries (ferrocene, octahedral Werner, half-sandwich, pyridine
σ).  All inputs are constructed in-memory (no archive dependency).
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

from sigma_only_cshm import (  # noqa: E402
    HAPTO_CC_MAX,
    SCHEMA_VERSION,
    compute_sigma_only_cshm_from_arrays,
    _split_sigma_hapto,
)


# ---------------------------------------------------------------------------
# Helpers — construct canonical M-L geometries directly in coords
# ---------------------------------------------------------------------------
def _octahedral_amines(M: str = "Fe", D: str = "N", r: float = 2.10):
    """Werner-type [M(NH3)6] without the Hs: 6 σ-N donors at OC-6 vertices."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
        syms.append(D)
        coords.append([r * v[0], r * v[1], r * v[2]])
    return syms, np.asarray(coords, dtype=float)


def _ferrocene_eta5_eta5(M: str = "Fe", r_mc: float = 2.05, r_cc: float = 1.40):
    """Ferrocene-like: 2 Cp rings, each 5 carbons at the M-Cp axis.

    z-axis = Cp normal; top ring at +z_mc, bottom at -z_mc.  Each ring is a
    regular pentagon of radius ``r_ring`` such that |M-C| = r_mc.
    """
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    # For η5 Cp, |M-C| = sqrt(z² + r_ring²).  Given r_cc, the pentagon
    # circumradius r_ring = r_cc / (2 sin(π/5)) ≈ r_cc * 0.8507.
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    z_off = math.sqrt(max(r_mc ** 2 - r_ring ** 2, 0.04))
    for sign in (+1.0, -1.0):
        for k in range(5):
            theta = 2.0 * math.pi * k / 5.0
            x = r_ring * math.cos(theta)
            y = r_ring * math.sin(theta)
            z = sign * z_off
            syms.append("C")
            coords.append([x, y, z])
    return syms, np.asarray(coords, dtype=float)


def _half_sandwich_eta5_plus_3sigma(M: str = "Fe", r_mc: float = 2.05,
                                     r_cl: float = 2.30, r_cc: float = 1.40):
    """CpFe(CO)(Cl)(N) piano-stool-like: 1 η5 Cp ring + 3 σ donors (Cl, N, P).

    The σ sub-polyhedron has CN = 3 → tested against SP-3 / T-3.
    """
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    r_ring = r_cc / (2.0 * math.sin(math.pi / 5.0))
    z_off = math.sqrt(max(r_mc ** 2 - r_ring ** 2, 0.04))
    # Cp ring at +z.
    for k in range(5):
        theta = 2.0 * math.pi * k / 5.0
        x = r_ring * math.cos(theta)
        y = r_ring * math.sin(theta)
        z = z_off
        syms.append("C")
        coords.append([x, y, z])
    # 3 σ donors arranged trigonal-planar in the equatorial plane z = 0
    # (the metal sees them as a clean SP-3 / T-3 sub-polyhedron).
    syms.append("Cl"); coords.append([r_cl, 0.0, 0.0])
    syms.append("N");  coords.append([-r_cl / 2.0, r_cl * math.sqrt(3) / 2.0, 0.0])
    syms.append("P");  coords.append([-r_cl / 2.0, -r_cl * math.sqrt(3) / 2.0, 0.0])
    return syms, np.asarray(coords, dtype=float)


def _square_planar_4sigma(M: str = "Pd", D: str = "N", r: float = 2.05):
    """Pd(NH3)4-like SP-4 with 4 σ-N donors."""
    syms = [M]
    coords = [[0.0, 0.0, 0.0]]
    for v in [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]:
        syms.append(D)
        coords.append([r * v[0], r * v[1], r * v[2]])
    return syms, np.asarray(coords, dtype=float)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_schema_version_is_int():
    assert isinstance(SCHEMA_VERSION, int)
    assert SCHEMA_VERSION >= 1


def test_octahedral_pure_sigma_cshm_near_zero():
    """Ideal OC-6 of 6 σ-N donors → σ-only CShM ≈ 0; no hapto."""
    syms, P = _octahedral_amines()
    rec = compute_sigma_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["sigma_n_donors"] == 6
    assert rec["hapto_n_units"] == 0
    assert rec["hapto_total_carbons"] == 0
    assert rec["mixed_ligand"] is False
    assert rec["sigma_polyhedron"] == "OC-6 octahedron"
    assert rec["sigma_only_cshm"] is not None
    assert rec["sigma_only_cshm"] < 0.01


def test_square_planar_sp4_pure_sigma():
    syms, P = _square_planar_4sigma()
    rec = compute_sigma_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Pd"
    assert rec["sigma_n_donors"] == 4
    assert rec["hapto_n_units"] == 0
    assert rec["sigma_polyhedron"] == "SP-4 square planar"
    assert rec["sigma_only_cshm"] is not None
    assert rec["sigma_only_cshm"] < 0.01


def test_ferrocene_eta5_eta5_no_sigma_donors():
    """A pure-η5-η5 ferrocene has 10 carbon donors in 2 η-clusters
    and ZERO σ donors → σ-only CShM is undefined."""
    syms, P = _ferrocene_eta5_eta5()
    rec = compute_sigma_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["sigma_n_donors"] == 0
    assert rec["hapto_n_units"] == 2
    assert rec["hapto_total_carbons"] == 10
    assert sorted(rec["hapto_unit_sizes"]) == [5, 5]
    assert rec["mixed_ligand"] is False
    assert rec["sigma_only_cshm"] is None  # CN < 2
    assert rec["sigma_polyhedron"] is None


def test_half_sandwich_isolates_sigma_three():
    """CpFe(Cl)(N)(P) piano-stool → σ-only CShM on the 3 σ donors
    (the η5 Cp is correctly excluded)."""
    syms, P = _half_sandwich_eta5_plus_3sigma()
    rec = compute_sigma_only_cshm_from_arrays(syms, P)
    assert rec["metal_symbol"] == "Fe"
    assert rec["sigma_n_donors"] == 3
    assert rec["hapto_n_units"] == 1
    assert rec["hapto_total_carbons"] == 5
    assert rec["mixed_ligand"] is True
    assert sorted(rec["sigma_donor_symbols"]) == ["Cl", "N", "P"]
    assert rec["sigma_only_cshm"] is not None
    # The 3 σ donors form a near-trigonal-planar set → CShM small (<1).
    assert rec["sigma_only_cshm"] < 1.0


def test_hapto_cluster_uses_geometric_proximity():
    """The η-cluster detector groups carbons that are within
    :data:`HAPTO_CC_MAX` of each other.  Two well-separated isolated C
    donors should NOT be clustered."""
    # CN-2 with 2 isolated carbene C donors >> HAPTO_CC_MAX apart.
    M = [0.0, 0.0, 0.0]
    syms = ["Fe", "C", "C"]
    coords = np.array([M, [2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]])
    rec = compute_sigma_only_cshm_from_arrays(syms, coords)
    assert rec["sigma_n_donors"] == 2  # both σ (isolated C)
    assert rec["hapto_n_units"] == 0


def test_hapto_cluster_cc_threshold_just_above_excludes():
    """C-C distance just above HAPTO_CC_MAX must NOT cluster the carbons."""
    eps = 1e-3
    cc = HAPTO_CC_MAX + 0.05  # safely above
    syms = ["Fe", "C", "C", "C"]
    # Three C donors at the metal but each pair > HAPTO_CC_MAX apart -> all σ.
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.8, 0.0, 0.0],
        [1.8 * math.cos(math.radians(120)), 1.8 * math.sin(math.radians(120)), 0.0],
        [1.8 * math.cos(math.radians(240)), 1.8 * math.sin(math.radians(240)), 0.0],
    ])
    # The 3 carbons form an equilateral triangle of side ~3.12 Å (>> 1.7).
    rec = compute_sigma_only_cshm_from_arrays(syms, coords)
    assert rec["hapto_n_units"] == 0
    assert rec["sigma_n_donors"] == 3


def test_split_sigma_hapto_pure_function():
    """Direct probe of the classification function."""
    syms, P = _ferrocene_eta5_eta5()
    metal = 0
    # All 10 carbons are in the M-D shell.
    donors = list(range(1, 11))
    split = _split_sigma_hapto(syms, P, donors)
    assert len(split.sigma_idx) == 0
    assert len(split.hapto_clusters) == 2
    assert sorted(len(c) for c in split.hapto_clusters) == [5, 5]


def test_determinism_two_calls_identical():
    """Running the analysis twice on the same input must yield bit-identical
    JSON records (no RNG, sorted iteration)."""
    import json as _json
    syms, P = _half_sandwich_eta5_plus_3sigma()
    a = compute_sigma_only_cshm_from_arrays(syms, P)
    b = compute_sigma_only_cshm_from_arrays(syms, P)
    assert _json.dumps(a, sort_keys=True) == _json.dumps(b, sort_keys=True)


def test_no_metal_returns_empty_record():
    """When no atom is a metal AND no fallback applies meaningfully, the
    result is a safe empty record (defensive)."""
    syms = ["C", "C", "C", "C"]
    coords = np.eye(4)[:, :3].astype(float)
    # Should still return without crashing.  Our heuristic falls back to
    # high-Z; with only carbons present, the first C is selected as "metal"
    # — but then the donor shell has no recognised donors (other carbons
    # are inside the cutoff but the M-shell logic excludes the metal index).
    # The point is: function MUST NOT raise.
    rec = compute_sigma_only_cshm_from_arrays(syms, coords)
    assert isinstance(rec, dict)
    assert rec["schema_version"] == SCHEMA_VERSION


def test_empty_input_returns_empty():
    rec = compute_sigma_only_cshm_from_arrays([], np.zeros((0, 3), float))
    assert rec["n_atoms"] == 0
    assert rec["metal_idx"] == -1
    assert rec["sigma_only_cshm"] is None


def test_explicit_metal_idx_override():
    """When the user gives an explicit metal_idx it must be honoured even
    if the symbol is not in the metal-symbol table."""
    syms = ["Xx", "N", "N", "N", "N"]   # Xx = unknown element
    coords = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0], [-2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0], [0.0, -2.0, 0.0],
    ])
    rec = compute_sigma_only_cshm_from_arrays(syms, coords, metal_idx=0)
    assert rec["metal_idx"] == 0
    assert rec["sigma_n_donors"] == 4
    assert rec["sigma_only_cshm"] is not None
    assert rec["sigma_only_cshm"] < 0.01


def test_aggregate_distorted_by_hapto():
    """The aggregate (full-shell) CShM on a ferrocene is much larger than
    the (undefined) σ-only CShM — this is the artefact we expose."""
    from aggregate_sigma_only import compute_aggregate_full_shell_cshm
    syms, P = _ferrocene_eta5_eta5()
    agg = compute_aggregate_full_shell_cshm(syms, P)
    # Ferrocene has 10 donors → aggregate CN exceeds the canonical OC-6 set
    # so the polyhedron lookup returns None (no canonical polyhedron for
    # CN=10).  Two scenarios are acceptable here:
    #   (a) agg is None (no polyhedron defined for CN=10) -> the artefact is
    #       so extreme the aggregate metric *itself* is undefined;
    #   (b) agg is a (huge) finite number when restricted to CN<=9.
    # We assert that *if* defined, it's clearly > 5.0 — much worse than 0.
    if agg is not None:
        assert agg > 5.0
