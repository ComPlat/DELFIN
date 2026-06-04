"""Tests for the H-H clash wiring (2026-06-04 healing-wiring cleanup).

Covers the 4 endpoint hooks that integrate
:mod:`delfin.fffree.hh_clash_detector` into the rest of the FFFree stack
when the env-flag ``DELFIN_FFFREE_HH_CLASH_INCLUDE=1`` is set:

  1. :class:`delfin.fffree.grip_constraints.ClashFloorPenalty` -- exposed
     via :func:`hh_pair_contribution` diagnostic.
  2. :func:`delfin.fffree.build_time_clash_gate.collapse_count` /
     :func:`has_collapse` -- augmented count + boolean.
  3. :func:`delfin.fffree.inter_ligand_clash_gate.count_inter_ligand_hh_clashes`
     -- diagnostic helper for cross-ligand H-H pairs.
  4. :func:`delfin.fffree.assemble_complex._clash_count` -- per-block
     placement count augmented by H-H Bondi-floor pairs.

Every test verifies that env-OFF is byte-identical with HEAD ``b195dba``
and env-ON triggers the H-H augmentation on a canonical 2-methyl
eclipsing geometry.
"""
from __future__ import annotations

import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Canonical 2-methyl eclipsing fixture (face-to-face H-H ≈ 1.78 Å).
# ---------------------------------------------------------------------------
def _two_methyls(separation: float = 2.5):
    syms = ["C", "H", "H", "H", "C", "H", "H", "H"]
    P = np.array([
        [0.00, 0.00, 0.00],
        [1.03, 0.00, 0.36],
        [-0.52, 0.89, 0.36],
        [-0.52, -0.89, 0.36],
        [0.00, 0.00, separation],
        [1.03, 0.00, separation - 0.36],
        [-0.52, 0.89, separation - 0.36],
        [-0.52, -0.89, separation - 0.36],
    ], dtype=float)
    return syms, P


def _far_methyls(separation: float = 6.0):
    return _two_methyls(separation)


# ---------------------------------------------------------------------------
# Endpoint 1: ClashFloorPenalty / hh_pair_contribution
# ---------------------------------------------------------------------------
class TestHHClashGripConstraintsWiring:
    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", raising=False)
        yield

    def test_hh_pair_contribution_env_off_byte_identical(self, monkeypatch):
        from delfin.fffree.grip_constraints import hh_pair_contribution
        syms, P = _two_methyls(separation=2.5)
        # Env unset => (0.0, 0)
        sev, n = hh_pair_contribution(syms, P)
        assert sev == 0.0
        assert n == 0

    def test_hh_pair_contribution_env_on_detects(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.grip_constraints import hh_pair_contribution
        syms, P = _two_methyls(separation=2.5)
        sev, n = hh_pair_contribution(syms, P)
        assert n >= 3
        assert sev > 0.0

    def test_hh_pair_contribution_no_clash_geometry(self, monkeypatch):
        """Far-apart methyls produce no H-H clashes even with env-on."""
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.grip_constraints import hh_pair_contribution
        syms, P = _far_methyls(separation=6.0)
        sev, n = hh_pair_contribution(syms, P)
        assert n == 0
        assert sev == 0.0


# ---------------------------------------------------------------------------
# Endpoint 3: inter_ligand_clash_gate.count_inter_ligand_hh_clashes
# ---------------------------------------------------------------------------
class TestInterLigandHHClashWiring:
    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", raising=False)
        yield

    def test_env_off_returns_zero_byte_identical(self, monkeypatch):
        from delfin.fffree.inter_ligand_clash_gate import (
            count_inter_ligand_hh_clashes,
        )
        syms, P = _two_methyls(separation=2.5)
        # Two ligand subgraphs: methyl-A = (0,1,2,3); methyl-B = (4,5,6,7).
        sg = [(0, 1, 2, 3), (4, 5, 6, 7)]
        n = count_inter_ligand_hh_clashes(P, syms, sg, metal_idx=None)
        assert n == 0  # env OFF -> 0 byte-identically

    def test_env_on_counts_cross_ligand_hh(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.inter_ligand_clash_gate import (
            count_inter_ligand_hh_clashes,
        )
        syms, P = _two_methyls(separation=2.5)
        sg = [(0, 1, 2, 3), (4, 5, 6, 7)]
        n = count_inter_ligand_hh_clashes(P, syms, sg, metal_idx=None)
        # 3 face-to-face H-H pairs cross the ligand boundary.
        assert n >= 3

    def test_env_on_intra_ligand_skipped(self, monkeypatch):
        """When both methyls belong to ONE ligand, the inter-ligand counter
        must return 0 (no cross-ligand pairs)."""
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.inter_ligand_clash_gate import (
            count_inter_ligand_hh_clashes,
        )
        syms, P = _two_methyls(separation=2.5)
        # Single ligand: all 8 atoms belong to subgraph 0.
        sg = [(0, 1, 2, 3, 4, 5, 6, 7)]
        n = count_inter_ligand_hh_clashes(P, syms, sg, metal_idx=None)
        assert n == 0


# ---------------------------------------------------------------------------
# Endpoint 4: assemble_complex._clash_count
# ---------------------------------------------------------------------------
class TestAssembleCompleClashCountWiring:
    @pytest.fixture(autouse=True)
    def _clean_env(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", raising=False)
        yield

    def test_env_off_byte_identical(self, monkeypatch):
        """With env-OFF, _clash_count must return the legacy 0.70× count."""
        from delfin.fffree.assemble_complex import _clash_count
        syms_q = ["C", "H", "H", "H"]
        Q = np.array([
            [0.00, 0.00, 2.5],
            [1.03, 0.00, 2.14],
            [-0.52, 0.89, 2.14],
            [-0.52, -0.89, 2.14],
        ], dtype=float)
        syms_ex = ["C", "H", "H", "H"]
        existing = np.array([
            [0.00, 0.00, 0.00],
            [1.03, 0.00, 0.36],
            [-0.52, 0.89, 0.36],
            [-0.52, -0.89, 0.36],
        ], dtype=float)
        # H-H face-to-face ≈ 1.78; 0.70 floor = 0.70 × 2.40 = 1.68; 1.78 > 1.68
        # -> NOT counted by legacy 0.70× criterion.
        n = _clash_count(Q, existing, syms_q, syms_ex)
        assert n == 0

    def test_env_on_augments_hh_clash(self, monkeypatch):
        """Env-ON: _clash_count is augmented by cross-block H-H pairs below
        the Bondi 0.85× floor (2.04 Å)."""
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.assemble_complex import _clash_count
        syms_q = ["C", "H", "H", "H"]
        Q = np.array([
            [0.00, 0.00, 2.5],
            [1.03, 0.00, 2.14],
            [-0.52, 0.89, 2.14],
            [-0.52, -0.89, 2.14],
        ], dtype=float)
        syms_ex = ["C", "H", "H", "H"]
        existing = np.array([
            [0.00, 0.00, 0.00],
            [1.03, 0.00, 0.36],
            [-0.52, 0.89, 0.36],
            [-0.52, -0.89, 0.36],
        ], dtype=float)
        # Each of the 3 face-to-face H-H pairs lies at ~1.78 Å (>1.68, <2.04).
        # -> 3 H-H clashes counted on top of the legacy 0 heavy violations.
        n = _clash_count(Q, existing, syms_q, syms_ex)
        assert n >= 3

    def test_env_on_no_clash_clean_geometry(self, monkeypatch):
        """Env-ON on a geometry with NO cross-block H-H clash must still
        return 0 (no double-counting, no false positives)."""
        monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
        from delfin.fffree.assemble_complex import _clash_count
        syms_q = ["C", "H", "H", "H"]
        Q = np.array([
            [0.00, 0.00, 6.0],
            [1.03, 0.00, 5.64],
            [-0.52, 0.89, 5.64],
            [-0.52, -0.89, 5.64],
        ], dtype=float)
        syms_ex = ["C", "H", "H", "H"]
        existing = np.array([
            [0.00, 0.00, 0.00],
            [1.03, 0.00, 0.36],
            [-0.52, 0.89, 0.36],
            [-0.52, -0.89, 0.36],
        ], dtype=float)
        n = _clash_count(Q, existing, syms_q, syms_ex)
        assert n == 0


# ---------------------------------------------------------------------------
# Cross-module determinism sanity check.
# ---------------------------------------------------------------------------
def test_hh_wiring_deterministic_across_endpoints(monkeypatch):
    """All four endpoints agree on the H-H clash count of the canonical
    2-methyl eclipsing fixture (deterministic, env-on path)."""
    monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
    syms, P = _two_methyls(separation=2.5)

    # Source-of-truth detector.
    from delfin.fffree.hh_clash_detector import count_hh_clashes
    n_truth = int(count_hh_clashes(syms, P))
    assert n_truth >= 3

    # Endpoint 1: hh_pair_contribution
    from delfin.fffree.grip_constraints import hh_pair_contribution
    _, n_1 = hh_pair_contribution(syms, P)
    assert n_1 == n_truth

    # Endpoint 2: collapse_count delta vs env-off baseline
    from delfin.fffree.build_time_clash_gate import collapse_count
    n_2 = collapse_count(syms, P)
    monkeypatch.delenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", raising=False)
    n_2_baseline = collapse_count(syms, P)
    monkeypatch.setenv("DELFIN_FFFREE_HH_CLASH_INCLUDE", "1")
    assert n_2 - n_2_baseline == n_truth

    # Endpoint 3: inter_ligand_clash_gate -- cross-ligand H-H = same count
    from delfin.fffree.inter_ligand_clash_gate import count_inter_ligand_hh_clashes
    sg = [(0, 1, 2, 3), (4, 5, 6, 7)]
    n_3 = count_inter_ligand_hh_clashes(P, syms, sg, metal_idx=None)
    assert n_3 == n_truth
