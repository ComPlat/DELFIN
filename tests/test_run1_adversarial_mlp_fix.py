"""Reviewer adversarial tests for the R1 fix (RUN1, nacht-s25).

Drive the real public path (MlpScreenWorkflow.run -> registry ->
contract-faithful fake mlp_rank adapter, registered AFTER lazy
discovery) and attack the fixed reporting logic with the shapes the
builder's own control does not cover: None/missing index, string and
bool energies, an all-crashed ranking.

A1: index=None finite entry — best_index must never be None.
A2: string energy (a decode artifact) — must be counted as non-finite
    and warned, not crash.
A3: all entries non-finite — no best_*, warning present; the result
    must not claim a clean success silently (warnings present).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from delfin.tools import _registry as reg
from delfin.tools._types import StepResult, StepStatus
from delfin.workflows.contrib.w21_mlp_screen import MlpScreenWorkflow


class _ShapeFake:
    name = "mlp_rank"
    description = "reviewer shape fake"

    def __init__(self, ranking):
        self.ranking = ranking

    def validate_params(self, **kwargs):
        return None

    def execute(self, work_dir, **kwargs):
        return StepResult(
            step_name="mlp_rank", status=StepStatus.SUCCESS,
            data={"ranking": self.ranking,
                  "n_structures": len(self.ranking)},
        )


def _run(tmp_path, monkeypatch, ranking):
    reg._discover_builtin_adapters()
    # monkeypatch.setitem is undone after each test — without it the
    # fake REPLACES the real mlp_rank adapter in reg._REGISTRY for the
    # rest of the CI process (CI #409 shard 3 / #410 shard 1: the
    # run_cli probe then succeeded on a missing file, and
    # test_tab_tools hit AttributeError on _ShapeFake.contract).
    monkeypatch.setitem(reg._REGISTRY, "mlp_rank", _ShapeFake(ranking))
    return MlpScreenWorkflow().run(
        config={}, ensemble="ens.xyz", work_dir=str(tmp_path / "wd"))


def test_a1_none_index_finite_energy_never_best(tmp_path, monkeypatch):
    res = _run(tmp_path, monkeypatch, [
        {"index": None, "energy_eV": -99.0},   # nonsense index
        {"index": 1, "energy_eV": -10.0},
    ])
    assert res.get("best_index") == 1, res
    assert res.get("best_energy_eV") == -10.0, res


def test_a2_string_energy_counts_nonfinite_and_warns(tmp_path, monkeypatch):
    res = _run(tmp_path, monkeypatch, [
        {"index": 0, "energy_eV": "-12.0"},    # string, not a number
        {"index": 1, "energy_eV": -10.0},
    ])
    assert res.get("best_index") == 1, res
    assert res.get("best_energy_eV") == -10.0, res
    assert res.get("warnings"), res


def test_a3_all_nonfinite_no_best_but_warned(tmp_path, monkeypatch):
    res = _run(tmp_path, monkeypatch, [
        {"index": 0, "energy_eV": float("inf")},
        {"index": 1, "energy_eV": float("nan")},
    ])
    assert "best_index" not in res, res
    assert "best_energy_eV" not in res, res
    assert res.get("warnings"), res
