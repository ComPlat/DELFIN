"""Tests for the GRACE Pareto-frontier + adaptive ensemble layer
(``delfin.fffree.grace_pareto``).

Covers:

  * env-flag default-OFF byte-identity (no Pareto reselect happens);
  * Pareto-dominance correctness on synthetic axis vectors;
  * frontier output is lex-sorted + deterministic across two runs;
  * dominator-count ranking of the dominated tail;
  * ``pareto_select`` fills exactly to ``max_keep`` when frontier is
    smaller than ``max_keep``;
  * adaptive plateau detection (no growth → ``should_stop`` True);
  * adaptive stop respects ``max_ensemble`` hard cap;
  * aromatic-planarity score is 0 for a perfectly flat benzene and
    non-zero for a puckered ring;
  * the legacy code path (``_rmsd_dedup`` with env OFF) is byte-identical
    to itself across two invocations.
"""
from __future__ import annotations

import os

# Strict determinism BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np
import pytest

from delfin.fffree import grace_pareto as GP


# ---------------------------------------------------------------------------
# Synthetic candidate stub (mirrors GraceCandidate's fields used by
# candidate_axes / pareto_select).
# ---------------------------------------------------------------------------
@dataclass
class _Stub:
    severity: float = 0.0
    clash_count: int = 0
    cshm: float = 0.0
    P: np.ndarray = field(default_factory=lambda: np.zeros((1, 3)))
    score: float = 0.0
    isomer_id: int = 0
    ring_id: int = 0
    rotamer_id: int = 0
    label: str = ""


# ---------------------------------------------------------------------------
# 1) Env-gate default OFF
# ---------------------------------------------------------------------------
def test_pareto_active_default_off(monkeypatch):
    monkeypatch.delenv(GP.ENV_PARETO_ENABLE, raising=False)
    assert GP.pareto_active() is False


def test_adaptive_active_default_off(monkeypatch):
    monkeypatch.delenv(GP.ENV_ADAPTIVE_ENABLE, raising=False)
    assert GP.adaptive_active() is False


def test_pareto_active_on(monkeypatch):
    monkeypatch.setenv(GP.ENV_PARETO_ENABLE, "1")
    assert GP.pareto_active() is True


def test_resolve_defaults(monkeypatch):
    monkeypatch.delenv(GP.ENV_MAX_ENSEMBLE, raising=False)
    monkeypatch.delenv(GP.ENV_ADAPTIVE_PATIENCE, raising=False)
    assert GP.resolve_max_ensemble() == GP.DEFAULT_MAX_ENSEMBLE
    assert GP.resolve_adaptive_patience() == GP.DEFAULT_ADAPTIVE_PATIENCE


def test_resolve_overrides(monkeypatch):
    monkeypatch.setenv(GP.ENV_MAX_ENSEMBLE, "17")
    monkeypatch.setenv(GP.ENV_ADAPTIVE_PATIENCE, "5")
    assert GP.resolve_max_ensemble() == 17
    assert GP.resolve_adaptive_patience() == 5


# ---------------------------------------------------------------------------
# 2) Dominance
# ---------------------------------------------------------------------------
def test_dominates_strict():
    # a dominates b when strictly better on one axis and not worse on others.
    a = (0.0, 0.0, 0.0, 0.0)
    b = (1.0, 0.0, 0.0, 0.0)
    assert GP.dominates(a, b)
    assert not GP.dominates(b, a)


def test_dominates_no_axis_strictly_better():
    a = (0.0, 0.0, 0.0, 0.0)
    b = (0.0, 0.0, 0.0, 0.0)
    # Identical vectors: neither dominates.
    assert not GP.dominates(a, b)
    assert not GP.dominates(b, a)


def test_dominates_mixed_better_worse():
    # a is better on first axis but worse on second: neither dominates.
    a = (0.0, 1.0, 0.0, 0.0)
    b = (1.0, 0.0, 0.0, 0.0)
    assert not GP.dominates(a, b)
    assert not GP.dominates(b, a)


def test_dominates_tolerance():
    # Differences below tol → treated as equal.
    a = (1.0, 1.0, 1.0, 1.0)
    b = (1.0 + 1e-16, 1.0, 1.0, 1.0)
    # a is not strictly better than b within tol.
    assert not GP.dominates(a, b)
    assert not GP.dominates(b, a)


# ---------------------------------------------------------------------------
# 3) Pareto frontier on synthetic stubs
# ---------------------------------------------------------------------------
def test_pareto_frontier_two_dominant_one_dominated():
    # Three candidates: two on the frontier, one strictly dominated.
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    b = _Stub(severity=2.0, clash_count=0, cshm=1.0, score=2.0, label="b")
    c = _Stub(severity=3.0, clash_count=1, cshm=3.0, score=3.0, label="c")
    dom, ddom = GP.compute_pareto_frontier([a, b, c], mol=None)
    labels_dom = sorted([d.label for d in dom])
    labels_ddom = [t[0].label for t in ddom]
    assert labels_dom == ["a", "b"]
    assert labels_ddom == ["c"]


def test_pareto_frontier_all_dominant():
    # Diagonal frontier — every point is non-dominated.
    cands = [
        _Stub(severity=float(i), clash_count=0, cshm=float(4 - i),
              score=float(i), label=f"x{i}")
        for i in range(4)
    ]
    dom, ddom = GP.compute_pareto_frontier(cands, mol=None)
    assert len(dom) == 4
    assert ddom == []


def test_pareto_frontier_lex_sorted_deterministic():
    # Output order must be deterministic across two runs.
    cands = [
        _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a"),
        _Stub(severity=2.0, clash_count=0, cshm=1.0, score=2.0, label="b"),
        _Stub(severity=1.0, clash_count=1, cshm=2.0, score=1.5, label="c"),
    ]
    dom1, _ = GP.compute_pareto_frontier(cands, mol=None)
    dom2, _ = GP.compute_pareto_frontier(list(reversed(cands)), mol=None)
    assert [c.label for c in dom1] == [c.label for c in dom2]


def test_pareto_select_fills_to_max_keep():
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    b = _Stub(severity=2.0, clash_count=0, cshm=1.0, score=2.0, label="b")
    c = _Stub(severity=3.0, clash_count=1, cshm=3.0, score=3.0, label="c")
    d = _Stub(severity=4.0, clash_count=2, cshm=4.0, score=4.0, label="d")
    sel = GP.pareto_select([a, b, c, d], max_keep=3, mol=None)
    assert len(sel) == 3
    # Frontier members must be retained.
    labels = sorted([s.label for s in sel])
    assert "a" in labels and "b" in labels


def test_pareto_select_empty():
    assert GP.pareto_select([], max_keep=3) == []


def test_pareto_select_zero_max_keep():
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    assert GP.pareto_select([a], max_keep=0) == []


def test_pareto_select_frontier_only_caps():
    # All four candidates are frontier; max_keep=2 must truncate.
    cands = [
        _Stub(severity=float(i), clash_count=0, cshm=float(4 - i),
              score=float(i), label=f"x{i}")
        for i in range(4)
    ]
    sel = GP.pareto_select(cands, max_keep=2, mol=None)
    assert len(sel) == 2


# ---------------------------------------------------------------------------
# 4) Adaptive plateau
# ---------------------------------------------------------------------------
def test_adaptive_plateau_no_growth_stops():
    # Three identical candidates after a single frontier point → plateau.
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    bad = _Stub(severity=2.0, clash_count=1, cshm=3.0, score=2.0, label="bad")
    state = GP.AdaptiveState(patience=3, max_ensemble=50)
    state.update(a)
    for _ in range(3):
        state.update(bad)
    assert state.should_stop() is True


def test_adaptive_max_ensemble_hard_cap():
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    state = GP.AdaptiveState(patience=999, max_ensemble=2)
    state.update(a)
    state.update(a)
    assert state.should_stop() is True


def test_adaptive_growth_resets_plateau():
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    b = _Stub(severity=0.5, clash_count=0, cshm=2.0, score=0.5, label="b")
    bad = _Stub(severity=10.0, clash_count=5, cshm=10.0, score=20.0, label="bad")
    state = GP.AdaptiveState(patience=3, max_ensemble=50)
    state.update(a)
    state.update(bad)  # plateau 1
    state.update(b)    # new frontier point → reset plateau
    assert state.plateau_run == 0
    assert state.should_stop() is False


def test_adaptive_plateau_check_helper():
    a = _Stub(severity=1.0, clash_count=0, cshm=2.0, score=1.0, label="a")
    bad = _Stub(severity=10.0, clash_count=5, cshm=10.0, score=20.0, label="bad")
    assert GP.adaptive_plateau_check([a, bad, bad, bad], patience=3) is True
    assert GP.adaptive_plateau_check([a, bad], patience=3) is False


# ---------------------------------------------------------------------------
# 5) Aromatic planarity score
# ---------------------------------------------------------------------------
def test_arom_plan_flat_benzene_zero():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=0)
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    conf = mol.GetConformer()
    P = np.asarray(conf.GetPositions(), dtype=np.float64)
    score = GP.compute_arom_plan_score(P, mol)
    # UFF-optimised benzene should be near-planar.
    assert score < 0.10


def test_arom_plan_no_aromatic_ring_zero():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=0)
    P = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    assert GP.compute_arom_plan_score(P, mol) == 0.0


def test_arom_plan_puckered_nonzero():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=0)
    conf = mol.GetConformer()
    P = np.asarray(conf.GetPositions(), dtype=np.float64)
    # Force a deterministic z-puck on the first aromatic atom.
    P[0, 2] += 0.5
    score = GP.compute_arom_plan_score(P, mol)
    assert score > 0.1


# ---------------------------------------------------------------------------
# 6) Determinism — two-run byte-identity
# ---------------------------------------------------------------------------
def test_pareto_two_run_byte_identical():
    cands = [
        _Stub(severity=float(i % 3), clash_count=i % 2,
              cshm=float((i * 7) % 5), score=float(i),
              label=f"x{i:03d}")
        for i in range(10)
    ]
    s1 = GP.pareto_select(cands, max_keep=5, mol=None)
    s2 = GP.pareto_select(list(reversed(cands)), max_keep=5, mol=None)
    assert [c.label for c in s1] == [c.label for c in s2]


# ---------------------------------------------------------------------------
# 7) Default-OFF: _rmsd_dedup byte-identity vs legacy truncation
# ---------------------------------------------------------------------------
def test_rmsd_dedup_default_off_byte_identity(monkeypatch):
    """When Pareto env-gate is OFF, _maybe_pareto_reselect must be a
    pure truncation pass (same as the legacy behaviour)."""
    monkeypatch.delenv(GP.ENV_PARETO_ENABLE, raising=False)
    cands = [
        _Stub(severity=float(i), clash_count=0, cshm=float(i),
              score=float(i), label=f"x{i:03d}")
        for i in range(8)
    ]
    from delfin.fffree.grace_ensemble import _maybe_pareto_reselect
    out = _maybe_pareto_reselect(cands, max_keep=3, mol=None)
    assert [c.label for c in out] == ["x000", "x001", "x002"]
