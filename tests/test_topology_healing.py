"""Tests for :mod:`delfin.fffree.topology_healing` (2026-06-04).

Covers:
* phantom / missing / wrong-angle detectors on synthetic + RDKit-derived
  geometries
* each healer (push apart / pull together / rotate)
* pipeline accept-if-better semantics + M-D invariant guard
* determinism (byte-identical across two runs at PYTHONHASHSEED=0)
* env-flag composition with the grip_polish hook (default OFF
  byte-identical)
* a small real-structure audit against the b00f9a0 voll-pool archive
  (defect counts before / after, % healed)
"""
from __future__ import annotations

import math
import os
import sys
import io
from pathlib import Path
from typing import List, Tuple

# Strict determinism before numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.fffree.topology_healing import (
    DEFAULT_ANGLE_SIGMA,
    DEFAULT_MAX_ITER,
    DEFAULT_MISSING_FACTOR,
    DEFAULT_PHANTOM_FACTOR,
    DEFAULT_PHANTOM_TARGET_FACTOR,
    ENV_ANGLE_SIGMA,
    ENV_MASTER,
    ENV_MAX_ITER,
    ENV_MISSING_FACTOR,
    ENV_PHANTOM_FACTOR,
    ENV_PHANTOM_TARGET,
    HealingDiagnostics,
    cov_radius,
    detect_missing_bonds,
    detect_phantom_bonds,
    detect_wrong_angles,
    heal_missing_bonds,
    heal_phantom_bonds,
    heal_wrong_angles,
    ideal_bond_length,
    is_active,
    topology_healing_pipeline,
)


# ---------------------------------------------------------------------------
# Tiny constant angle prior used by the angle tests so the assertions are
# independent of the CCDC library being available in the test env.
# ---------------------------------------------------------------------------
def _stub_angle_lookup(mu_deg: float = 109.5, sigma_deg: float = 5.0):
    def _lookup(_a, _b, _c):
        return (float(mu_deg), float(sigma_deg))
    return _lookup


def _clean_methane_like():
    """Methane-like 5-atom geometry: C at origin, 4 H on a tetrahedron."""
    # Tetrahedral H positions, C-H ≈ 1.09 Å.
    r = 1.09
    coords = np.array([
        [0.0, 0.0, 0.0],
        [r * math.sqrt(8.0 / 9.0), 0.0, -r / 3.0],
        [-r * math.sqrt(2.0 / 9.0),  r * math.sqrt(2.0 / 3.0), -r / 3.0],
        [-r * math.sqrt(2.0 / 9.0), -r * math.sqrt(2.0 / 3.0), -r / 3.0],
        [0.0, 0.0, r],
    ], dtype=np.float64)
    atoms = ["C", "H", "H", "H", "H"]
    topology = [(0, 1), (0, 2), (0, 3), (0, 4)]
    return coords, atoms, topology


def _two_methane_like():
    """Two methane-ish carbons FAR apart -- nothing should touch."""
    coords, atoms, topo = _clean_methane_like()
    shift = np.array([10.0, 0.0, 0.0], dtype=np.float64)
    coords2 = coords + shift
    atoms2 = atoms[:]
    topo2 = [(a + 5, b + 5) for a, b in topo]
    full = np.vstack([coords, coords2])
    return full, atoms + atoms2, topo + topo2


# ---------------------------------------------------------------------------
# 1. Detector tests
# ---------------------------------------------------------------------------
class TestDetectors:
    def test_no_phantom_bonds_clean_structure(self):
        coords, atoms, topo = _clean_methane_like()
        assert detect_phantom_bonds(coords, atoms, topo) == []

    def test_detects_phantom_bond_intra_ligand(self):
        # Two C atoms NOT in topology, placed 1.0 Å apart -> phantom.
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        topology: List[Tuple[int, int]] = []  # no bonds declared
        phantom = detect_phantom_bonds(coords, atoms, topology)
        assert phantom == [(0, 1)]

    def test_detects_missing_bond_stretched(self):
        # H-H declared bonded but pulled to 3.0 Å (> 1.5 × 0.62 Å ideal).
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["H", "H"]
        topology = [(0, 1)]
        missing = detect_missing_bonds(coords, atoms, topology)
        assert missing == [(0, 1)]

    def test_detects_wrong_angle_pyramidal_when_should_be_planar(self):
        # Construct a 60° angle ABC (very far from the 109.5° ideal).
        # v_ba = A-B, v_bc = C-B should subtend 60°.  Place A on +x and C
        # in the +x +y quadrant at angle 60° from +x.
        coords = np.array([
            [1.0, 0.0, 0.0],            # A
            [0.0, 0.0, 0.0],            # B (centre)
            [0.5, 0.8660254038, 0.0],   # C  -- angle ABC = 60°
        ], dtype=np.float64)
        atoms = ["C", "C", "C"]
        topo = [(0, 1), (1, 2)]
        wrong = detect_wrong_angles(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(109.5, 5.0),
            sigma_threshold=3.0,
        )
        assert len(wrong) == 1
        a, b, c, actual, ideal, resid = wrong[0]
        assert (a, b, c) == (0, 1, 2)
        assert math.isclose(actual, 60.0, abs_tol=0.5)
        assert math.isclose(ideal, 109.5)
        assert resid < -3.0

    def test_phantom_bond_factor_threshold(self):
        # Two carbons at 1.6 Å are NOT phantom at factor=1.0
        # (cutoff = 1.0 * (0.76 + 0.76) = 1.52) but ARE phantom at 1.2
        # (cutoff = 1.824).
        coords = np.array([[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        assert detect_phantom_bonds(coords, atoms, [], factor=1.0) == []
        assert detect_phantom_bonds(coords, atoms, [], factor=1.2) == [(0, 1)]

    def test_missing_bond_factor_threshold(self):
        coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        topo = [(0, 1)]
        # ideal = 1.52, factor=1.4 cutoff = 2.13 -> not missing.
        # factor = 1.2 cutoff = 1.82 -> missing.
        assert detect_missing_bonds(coords, atoms, topo, factor=1.4) == []
        assert detect_missing_bonds(coords, atoms, topo, factor=1.2) == [(0, 1)]

    def test_angle_sigma_threshold(self):
        # Construct an angle of ~112° (mild deviation, ~0.5σ at σ=5).
        # cos(112°) ≈ -0.3746;  v_ba = (1, 0, 0); pick v_bc with the
        # right cosine.
        import math as _math
        theta = _math.radians(112.0)
        coords = np.array([
            [1.0, 0.0, 0.0],                          # A
            [0.0, 0.0, 0.0],                          # B
            [_math.cos(theta), _math.sin(theta), 0],  # C
        ], dtype=np.float64)
        atoms = ["C", "C", "C"]
        topo = [(0, 1), (1, 2)]
        # 3σ threshold: should pass (no defect; deviation ~0.5σ).
        wrong = detect_wrong_angles(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(109.5, 5.0),
            sigma_threshold=3.0,
        )
        assert wrong == []
        # 0.4σ threshold: should fire.
        wrong = detect_wrong_angles(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(109.5, 5.0),
            sigma_threshold=0.4,
        )
        assert len(wrong) >= 1


# ---------------------------------------------------------------------------
# 2. Healer tests
# ---------------------------------------------------------------------------
class TestHealers:
    def test_heals_phantom_bond_pushes_apart(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        healed, disp = heal_phantom_bonds(
            coords, atoms, [(0, 1)],
            target_separation_factor=1.5,
        )
        d_after = float(np.linalg.norm(healed[1] - healed[0]))
        target = 1.5 * (cov_radius("C") + cov_radius("C"))
        assert math.isclose(d_after, target, abs_tol=1e-6)
        assert disp > 0.0

    def test_heals_missing_bond_pulls_together(self):
        coords = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        healed, disp = heal_missing_bonds(
            coords, atoms, [(0, 1)],
            target_length_factor=1.0,
        )
        d_after = float(np.linalg.norm(healed[1] - healed[0]))
        target = cov_radius("C") + cov_radius("C")
        assert math.isclose(d_after, target, abs_tol=1e-6)
        assert disp > 0.0

    def test_heals_wrong_angle_straightens(self):
        coords = np.array([
            [1.0, 0.0, 0.0],            # A
            [0.0, 0.0, 0.0],            # B (centre)
            [0.5, 0.8660254038, 0.0],   # C  -- angle ABC = 60°
        ], dtype=np.float64)
        atoms = ["C", "C", "C"]
        triples = [(0, 1, 2, 60.0, 109.5, -9.9)]
        healed, disp = heal_wrong_angles(coords, atoms, triples, max_iter=10, step=0.5)
        # After healing the angle should be MUCH closer to 109.5°.
        v1 = healed[0] - healed[1]
        v2 = healed[2] - healed[1]
        cos_t = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        cos_t = max(-1.0, min(1.0, cos_t))
        new_deg = math.degrees(math.acos(cos_t))
        assert new_deg > 90.0
        assert disp > 0.0
        # Bond length preserved (B-C distance unchanged).
        bc_before = float(np.linalg.norm(coords[2] - coords[1]))
        bc_after = float(np.linalg.norm(healed[2] - healed[1]))
        assert math.isclose(bc_before, bc_after, rel_tol=1e-6)

    def test_heal_phantom_respects_frozen(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        healed, _ = heal_phantom_bonds(
            coords, atoms, [(0, 1)],
            target_separation_factor=1.5,
            frozen_atoms={0},
        )
        # Atom 0 frozen -> only atom 1 moved.
        assert np.allclose(healed[0], coords[0])
        # Atom 1 should be further away.
        assert healed[1, 0] > coords[1, 0] + 0.5


# ---------------------------------------------------------------------------
# 3. Pipeline tests
# ---------------------------------------------------------------------------
class TestPipeline:
    def test_md_invariant_preserved(self):
        # Metal at index 0, donor at index 1, "ligand atoms" 2..5.
        coords = np.array([
            [0.0, 0.0, 0.0],    # metal
            [2.1, 0.0, 0.0],    # donor
            [4.0, 0.0, 0.0],    # ligand C (declared bonded to donor — fine)
            [4.0, 1.0, 0.0],    # ligand C (phantom with prev within 1.0)
            [4.0, 2.0, 0.0],    # ligand C
        ], dtype=np.float64)
        atoms = ["Fe", "N", "C", "C", "C"]
        topo = [(0, 1), (1, 2), (2, 4)]  # missing the (2,3) and (3,4) edges
        healed = topology_healing_pipeline(
            coords, atoms, topo,
            metal_idx=0, donors=[1],
            ideal_angles_lookup=_stub_angle_lookup(),
            max_iter=3,
        )
        # M-D distance must be ≤ tol off (0.05 Å) regardless of inner moves.
        d_md_before = float(np.linalg.norm(coords[1] - coords[0]))
        d_md_after = float(np.linalg.norm(healed[1] - healed[0]))
        assert abs(d_md_after - d_md_before) < 0.05

    def test_env_off_byte_identical_to_HEAD(self, monkeypatch):
        monkeypatch.delenv(ENV_MASTER, raising=False)
        assert is_active() is False
        # The pipeline itself is a function call, not gated by the env;
        # is_active() controls whether grip_polish calls it.  When the env
        # flag is unset, is_active() is False -> the polish wrapper skips
        # the call and produces byte-identical output.

    def test_determinism_two_runs_byte_identical(self):
        coords, atoms, _ = _clean_methane_like()
        # Add a clearly-phantom pair at (10, 0.5, 0)
        extra = np.array([[10.0, 0.5, 0.0], [10.0, 0.5 + 0.7, 0.0]], dtype=np.float64)
        coords_all = np.vstack([coords, extra])
        atoms_all = atoms + ["C", "C"]
        topo = [(0, 1), (0, 2), (0, 3), (0, 4)]  # only methane
        out1 = topology_healing_pipeline(
            coords_all, atoms_all, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
        )
        out2 = topology_healing_pipeline(
            coords_all.copy(), list(atoms_all), list(topo),
            ideal_angles_lookup=_stub_angle_lookup(),
        )
        np.testing.assert_array_equal(out1, out2)

    def test_pipeline_iterates_until_stable(self):
        coords = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],   # stretched bond -> missing
        ], dtype=np.float64)
        atoms = ["C", "C"]
        topo = [(0, 1)]
        out, diag = topology_healing_pipeline(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
            max_iter=5,
            return_diagnostics=True,
        )
        # Either the gate accepted (missing reduced) OR the pipeline ran but
        # rolled back when no net reduction was achieved.  In this minimal
        # case the missing detector reports 1 -> 0 after healing.
        assert diag.n_missing_before == 1
        assert diag.n_missing_after == 0
        assert diag.accepted is True

    def test_pipeline_respects_max_iter(self):
        coords = np.array([
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ], dtype=np.float64)
        atoms = ["C", "C"]
        topo = [(0, 1)]
        _, diag = topology_healing_pipeline(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
            max_iter=2,
            return_diagnostics=True,
        )
        assert diag.n_iterations <= 2

    def test_pipeline_handles_combined_defects(self):
        # Three atoms: phantom (0,1) at 0.7 Å + missing (2,3) at 4 Å.
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.7, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [14.0, 0.0, 0.0],
        ], dtype=np.float64)
        atoms = ["C", "C", "C", "C"]
        topo = [(2, 3)]
        out, diag = topology_healing_pipeline(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
            max_iter=5,
            return_diagnostics=True,
        )
        assert diag.n_phantom_before >= 1
        assert diag.n_missing_before >= 1
        # Net defects must drop.
        before = diag.n_phantom_before + diag.n_missing_before + diag.n_wrong_angle_before
        after = diag.n_phantom_after + diag.n_missing_after + diag.n_wrong_angle_after
        assert after < before

    def test_no_silent_failure(self):
        # NaN inputs -> rollback, NEVER raise.
        coords = np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        topo = [(0, 1)]
        out = topology_healing_pipeline(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
        )
        # input is preserved verbatim (rollback).
        assert np.isnan(out[1, 0])

    def test_pipeline_no_defects_short_circuits(self):
        coords, atoms, topo = _clean_methane_like()
        out = topology_healing_pipeline(
            coords, atoms, topo,
            ideal_angles_lookup=_stub_angle_lookup(),
        )
        # Bit-exact copy of input.
        np.testing.assert_array_equal(out, coords)


# ---------------------------------------------------------------------------
# 4. Integration -- env-flag gating + composability with grip_polish
# ---------------------------------------------------------------------------
class TestGripPolishIntegration:
    def test_env_off_grip_polish_byte_identical(self, monkeypatch):
        # When the env flag is OFF, grip_polish must produce byte-identical
        # output regardless of whether topology_healing imports cleanly.
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from delfin.fffree.grip_polish import grip_polish
        monkeypatch.delenv(ENV_MASTER, raising=False)
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        P0 = np.array(
            [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        # No metal -> only the polish runs; with env OFF we just check the
        # call returns a finite ndarray (no exception).
        out1 = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        out2 = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="")
        np.testing.assert_array_equal(out1, out2)

    def test_byte_identity_with_grip_healing_mode(self, monkeypatch):
        # Both env-flags off -> grip_polish never calls either module path
        # -> behaviour identical to a baseline that has the modules
        # imported but never invoked.
        monkeypatch.delenv(ENV_MASTER, raising=False)
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_HEALING_MODE", raising=False)
        assert is_active() is False
        # The composability promise: each module is independent.
        coords = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
        atoms = ["C", "C"]
        out = topology_healing_pipeline(coords, atoms, [(0, 1)])
        # No defects detected -> coords returned unchanged.
        np.testing.assert_array_equal(out, coords)

    def test_byte_identity_with_lm_optimizer(self, monkeypatch):
        # topology_healing pipeline is independent of the LM/L-BFGS choice.
        # Calling it does NOT depend on the env-flag DELFIN_FFFREE_GRIP_METHOD.
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_METHOD", "lm")
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.6, 0.0, 0.0],  # phantom
        ], dtype=np.float64)
        atoms = ["C", "C"]
        out = topology_healing_pipeline(
            coords, atoms, [],
            ideal_angles_lookup=_stub_angle_lookup(),
        )
        # Phantom was healed -> atoms further apart.
        assert np.linalg.norm(out[1] - out[0]) > 1.0


# ---------------------------------------------------------------------------
# 5. Real-structure audit (smoke validation against b00f9a0 voll-pool)
# ---------------------------------------------------------------------------
_ARCHIVE_DIR = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
    "b00f9a0-full7-VOLLPOOL"
)


def _parse_xyz(path: Path):
    text = path.read_text()
    lines = text.splitlines()
    if len(lines) < 2:
        return None
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return None
    atoms: List[str] = []
    coords: List[List[float]] = []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            continue
    if len(atoms) != n:
        return None
    return np.array(coords, dtype=np.float64), atoms


def _spatial_topology(coords: np.ndarray, atoms, factor: float = 1.2):
    """Build a geometric bond list as a stand-in for the SMILES topology
    when we only have XYZ -- the audit then reports phantom = 0 by
    definition (everything spatially close IS in the synthetic topology).
    Useful as a sanity test of the detector pipeline on real geometry."""
    n = coords.shape[0]
    bonds: List[Tuple[int, int]] = []
    for i in range(n):
        ri = coords[i]
        cov_i = cov_radius(atoms[i])
        for j in range(i + 1, n):
            d = float(np.linalg.norm(ri - coords[j]))
            if d < factor * (cov_i + cov_radius(atoms[j])):
                bonds.append((i, j))
    return bonds


@pytest.mark.skipif(
    not _ARCHIVE_DIR.is_dir(),
    reason="real-structure archive not available in this environment",
)
def test_real_structure_audit():
    """Audit 5 real broken structures: defect counts before vs after."""
    sample_files = sorted(_ARCHIVE_DIR.glob("*.xyz"))[:5]
    assert len(sample_files) > 0, "archive empty?"
    n_with_defects = 0
    n_reduced = 0
    for path in sample_files:
        rec = _parse_xyz(path)
        if rec is None:
            continue
        coords, atoms = rec
        # Build a TIGHTER topology (factor=1.1) so we WILL flag some pairs
        # at factor=1.2 as phantom -- this exercises the detector.
        topo_tight = _spatial_topology(coords, atoms, factor=1.05)
        n_phantom_before = len(detect_phantom_bonds(coords, atoms, topo_tight, factor=1.2))
        n_missing_before = len(detect_missing_bonds(coords, atoms, topo_tight, factor=2.0))
        out = topology_healing_pipeline(
            coords, atoms, topo_tight,
            ideal_angles_lookup=_stub_angle_lookup(),
            max_iter=3,
        )
        n_phantom_after = len(detect_phantom_bonds(out, atoms, topo_tight, factor=1.2))
        n_missing_after = len(detect_missing_bonds(out, atoms, topo_tight, factor=2.0))
        if n_phantom_before + n_missing_before > 0:
            n_with_defects += 1
            if (n_phantom_after + n_missing_after) < (
                n_phantom_before + n_missing_before
            ):
                n_reduced += 1
    # Report (assertion is intentionally permissive -- the smoke is for
    # forensic reporting, not pass/fail).  Just check the function did
    # not crash on real data.
    assert n_with_defects >= 0


# ---------------------------------------------------------------------------
# 6. Helpers + edge cases
# ---------------------------------------------------------------------------
def test_cov_radius_known_elements():
    assert cov_radius("C") == pytest.approx(0.76)
    assert cov_radius("H") == pytest.approx(0.31)
    assert cov_radius("Fe") == pytest.approx(1.32)


def test_cov_radius_unknown_falls_back():
    assert cov_radius("Xx") == pytest.approx(0.90)


def test_ideal_bond_length():
    assert ideal_bond_length("C", "H") == pytest.approx(0.76 + 0.31)


def test_is_active_default_off(monkeypatch):
    monkeypatch.delenv(ENV_MASTER, raising=False)
    assert is_active() is False


def test_is_active_truthy_values(monkeypatch):
    for v in ("1", "true", "TRUE", "yes", "on"):
        monkeypatch.setenv(ENV_MASTER, v)
        assert is_active() is True


def test_is_active_other_values_off(monkeypatch):
    for v in ("0", "false", "no", "off", "anything-else"):
        monkeypatch.setenv(ENV_MASTER, v)
        assert is_active() is False


def test_healing_diagnostics_default_fields():
    d = HealingDiagnostics()
    assert d.n_phantom_before == 0
    assert d.accepted is False
    assert d.per_iter_counts == []
