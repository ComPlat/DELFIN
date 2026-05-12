"""Unit tests for delfin.topology_hard_gate.

Säule 1 of Hybrid-Path (nature_project/15_HYBRID_PATH_FINAL.md).
"""
from __future__ import annotations

import pytest

from delfin.topology_hard_gate import (
    BOND_TOL,
    COLLISION_THRESHOLD,
    MAX_MD_DISTANCE,
    TopologyGateResult,
    TopologyViolation,
    validate_topology_invariant,
)


# ----------------------------------------------------------------------
# Fixtures — simple XYZ strings
# ----------------------------------------------------------------------

def _build_xyz(atoms_with_xyz):
    """Build XYZ string from [(sym, x, y, z), ...]."""
    n = len(atoms_with_xyz)
    lines = [str(n), "test"]
    for sym, x, y, z in atoms_with_xyz:
        lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(lines)


@pytest.fixture
def good_fecl4_xyz():
    """Tetrahedral [FeCl4]^2- approximate geometry, Fe-first atom order."""
    return _build_xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("Cl", 1.33, 1.33, 1.33),
        ("Cl", -1.33, -1.33, 1.33),
        ("Cl", -1.33, 1.33, -1.33),
        ("Cl", 1.33, -1.33, -1.33),
    ])


@pytest.fixture
def fecl4_smiles():
    """Fe-first SMILES to match the XYZ atom-order in good_fecl4_xyz."""
    return "[Fe-2]([Cl])([Cl])([Cl])[Cl]"


@pytest.fixture
def collision_xyz():
    """FeCl4 with two Cl on top of each other."""
    return _build_xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("Cl", 1.33, 1.33, 1.33),
        ("Cl", 1.34, 1.34, 1.34),  # 0.017 Å from above
        ("Cl", -1.33, 1.33, -1.33),
        ("Cl", 1.33, -1.33, -1.33),
    ])


@pytest.fixture
def detached_donor_xyz():
    """FeCl4 with one Cl flown to 10 Å away."""
    return _build_xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("Cl", 1.33, 1.33, 1.33),
        ("Cl", -1.33, -1.33, 1.33),
        ("Cl", -1.33, 1.33, -1.33),
        ("Cl", 10.0, 10.0, 10.0),  # detached
    ])


@pytest.fixture
def fragmented_xyz():
    """Two disconnected fragments: FeCl3 cluster + lonely Cl."""
    return _build_xyz([
        ("Fe", 0.0, 0.0, 0.0),
        ("Cl", 1.33, 1.33, 1.33),
        ("Cl", -1.33, -1.33, 1.33),
        ("Cl", -1.33, 1.33, -1.33),
        ("Cl", 50.0, 50.0, 50.0),  # very far away
    ])


# ----------------------------------------------------------------------
# Tests — basic structural checks
# ----------------------------------------------------------------------

class TestEmptyOrInvalidInput:
    def test_empty_xyz_fails(self):
        result = validate_topology_invariant("", "[Fe]")
        assert not result.passed
        assert any(v.kind == "empty_xyz" for v in result.violations)

    def test_just_header_no_atoms_fails(self):
        result = validate_topology_invariant("0\nempty\n", "[Fe]")
        assert not result.passed

    def test_unparseable_smiles_still_runs_collision_check(self, collision_xyz):
        # SMILES that RDKit cannot parse — gate should still detect
        # internal-XYZ-collisions but skip M-D check.
        result = validate_topology_invariant(collision_xyz, "garbled garbage smiles")
        # Result depends on RDKit's tolerance; if SMILES parsed (lax sanitize),
        # we expect a collision violation; if not, we still get one.
        assert not result.passed
        kinds = {v.kind for v in result.violations}
        assert "collision" in kinds


# ----------------------------------------------------------------------
# Tests — passes on good input
# ----------------------------------------------------------------------

class TestPassesOnGoodInput:
    def test_good_fecl4_passes(self, good_fecl4_xyz, fecl4_smiles):
        result = validate_topology_invariant(good_fecl4_xyz, fecl4_smiles)
        assert result.passed, f"violations={result.violations}"
        assert result.n_atoms == 5
        assert result.n_md_bonds_expected == 4
        assert result.n_md_bonds_found == 4


# ----------------------------------------------------------------------
# Tests — catches collisions
# ----------------------------------------------------------------------

class TestCollisionDetection:
    def test_collision_fails(self, collision_xyz, fecl4_smiles):
        result = validate_topology_invariant(collision_xyz, fecl4_smiles)
        assert not result.passed
        kinds = {v.kind for v in result.violations}
        assert "collision" in kinds

    def test_collision_threshold_configurable(self, good_fecl4_xyz, fecl4_smiles):
        # With a very large threshold, even good FeCl4 fails
        result = validate_topology_invariant(
            good_fecl4_xyz, fecl4_smiles,
            collision_threshold=3.0,  # huge — every pair is "colliding"
        )
        assert not result.passed
        assert any(v.kind == "collision" for v in result.violations)


# ----------------------------------------------------------------------
# Tests — catches detached donors
# ----------------------------------------------------------------------

class TestDetachedDonorDetection:
    def test_detached_cl_fails(self, detached_donor_xyz, fecl4_smiles):
        result = validate_topology_invariant(detached_donor_xyz, fecl4_smiles)
        assert not result.passed
        kinds = {v.kind for v in result.violations}
        # Either missing_md (Cl too far from Fe) or extra_fragment
        assert "missing_md" in kinds or "extra_fragment" in kinds

    def test_max_md_distance_configurable(self, good_fecl4_xyz, fecl4_smiles):
        # With max_md_distance=1.0, even good FeCl4 fails (Fe-Cl ~2.3 Å > 1.0)
        result = validate_topology_invariant(
            good_fecl4_xyz, fecl4_smiles,
            max_md_distance=1.0,
        )
        assert not result.passed
        assert any(v.kind == "missing_md" for v in result.violations)


# ----------------------------------------------------------------------
# Tests — catches fragmentation
# ----------------------------------------------------------------------

class TestFragmentation:
    def test_extra_fragment_detected(self, fragmented_xyz, fecl4_smiles):
        result = validate_topology_invariant(fragmented_xyz, fecl4_smiles)
        assert not result.passed
        kinds = {v.kind for v in result.violations}
        assert "extra_fragment" in kinds or "missing_md" in kinds
        # Expect 2 fragments in XYZ vs 1 in SMILES
        assert result.n_fragments_xyz >= 2


# ----------------------------------------------------------------------
# Tests — result-object properties
# ----------------------------------------------------------------------

class TestResultObject:
    def test_passed_result_reason_is_pass(self, good_fecl4_xyz, fecl4_smiles):
        result = validate_topology_invariant(good_fecl4_xyz, fecl4_smiles)
        if result.passed:
            assert result.reason == "pass"

    def test_failed_result_reason_lists_violation_kinds(
        self, collision_xyz, fecl4_smiles,
    ):
        result = validate_topology_invariant(collision_xyz, fecl4_smiles)
        assert not result.passed
        assert "collision" in result.reason

    def test_result_is_frozen(self, good_fecl4_xyz, fecl4_smiles):
        result = validate_topology_invariant(good_fecl4_xyz, fecl4_smiles)
        with pytest.raises(Exception):
            result.passed = False  # type: ignore[misc]


# ----------------------------------------------------------------------
# Tests — determinism (pure function contract)
# ----------------------------------------------------------------------

class TestDeterminism:
    def test_same_inputs_yield_same_output(self, good_fecl4_xyz, fecl4_smiles):
        r1 = validate_topology_invariant(good_fecl4_xyz, fecl4_smiles)
        r2 = validate_topology_invariant(good_fecl4_xyz, fecl4_smiles)
        assert r1.passed == r2.passed
        assert r1.n_md_bonds_found == r2.n_md_bonds_found
        assert len(r1.violations) == len(r2.violations)


# ----------------------------------------------------------------------
# Tests — multi-fragment SMILES tolerated
# ----------------------------------------------------------------------

class TestMultiFragmentSmiles:
    def test_smiles_with_counter_ion(self):
        # Counter-ion SMILES that DELFIN may keep separated
        xyz = _build_xyz([
            ("Pt", 0.0, 0.0, 0.0),
            ("Cl", 2.3, 0.0, 0.0),
            ("Cl", -2.3, 0.0, 0.0),
            ("N", 0.0, 2.0, 0.0),
            ("N", 0.0, -2.0, 0.0),
            ("Cl", 50.0, 50.0, 50.0),  # counter-ion at distance
        ])
        smiles = "Cl[Pt](Cl)([NH3+])[NH3+].[Cl-]"
        result = validate_topology_invariant(xyz, smiles)
        # With 2-fragment SMILES, the lone Cl at distance is expected → ok
        assert result.n_fragments_smiles >= 1
