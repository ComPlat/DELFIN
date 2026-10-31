"""RUN1 mlp_screen correctness controls (findings R1, R3 by nacht-s25).

Real call path: MlpScreenWorkflow.run() -> delfin.tools._registry
get("mlp_rank") -> adapter.execute(wd, ...).  The fake adapter is
contract-faithful: StepResult with data["ranking"] as the real mlp_rank
adapter produces it (structures whose MLP evaluation crashed come back
with float("inf") — delfin/mlp_tools/calculators.py:286-290).  The
workflow must not report an inf or index-invalid best silently.

Registry isolation: the fake adapter is registered into the GLOBAL
adapter registry, so every test here saves and restores the previous
"mlp_rank" entry — otherwise a fake leaks into later tests in the same
process (QS finding F-a: an order-dependent cross-file pollution made
test_zz_probe_run_cli_report green-when-lonely, red-in-set, because a
leaked fake never raised on a missing file).
"""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.tools import _registry as reg  # noqa: E402
from delfin.tools._types import StepResult, StepStatus  # noqa: E402
from delfin.workflows.contrib.w21_mlp_screen import (  # noqa: E402
    MlpScreenWorkflow,
)


class _FakeRankAdapter:
    """Contract-faithful mlp_rank fake returning a fixed ranking."""

    name = "mlp_rank"
    description = "test fake"

    def __init__(self, ranking):
        self._ranking = ranking

    def validate_params(self, backend="ani2x", charge=0, mult=1,
                        ensemble=None):
        assert ensemble, "validate_params must receive the ensemble path"

    def execute(self, work_dir, *, ensemble=None, backend="ani2x",
                charge=0, mult=1, device="cpu"):
        out = StepResult(
            step_name="mlp_rank", status=StepStatus.SUCCESS,
            data={"ranking": self._ranking,
                  "n_structures": len(self._ranking)},
        )
        (work_dir / "ranking.json").write_text(json.dumps(self._ranking))
        return out


@pytest.fixture
def _clean_registry():
    """Save the real mlp_rank entry; restore it after the test."""
    reg._discover_builtin_adapters()
    saved = reg._REGISTRY.get("mlp_rank")
    yield
    if saved is not None:
        reg.register(saved)


INF_RANKING = [
    {"index": 0, "energy_eV": float("inf")},
    {"index": 1, "energy_eV": -12.34},
    {"index": 2, "energy_eV": -10.0},
]

NONE_INDEX_RANKING = [
    {"index": None, "energy_eV": -99.0},
    {"index": 2, "energy_eV": -1.5},
]


def _run_with(ranking, tmp_path):
    reg.register(_FakeRankAdapter(ranking))
    return MlpScreenWorkflow().run(
        config={}, ensemble="ens.xyz",
        work_dir=str(tmp_path / "wd"))


def test_all_finite_ranking_does_not_warn(tmp_path, _clean_registry):
    """Negative control: a clean ranking must not gain a warning."""
    res = _run_with([
        {"index": 0, "energy_eV": -5.0},
        {"index": 1, "energy_eV": -7.5},
    ], tmp_path)
    assert res.get("warnings") is None
    assert res["best_index"] == 1 and res["best_energy_eV"] == -7.5


def test_silent_inf_is_reported_and_best_is_finite(tmp_path,
                                                   _clean_registry):
    res = _run_with(INF_RANKING, tmp_path)
    assert res.get("warnings"), (
        "no warning on an inf ranking energy: " + json.dumps(res))
    best = res.get("best_energy_eV")
    assert best is not None and best == best and best != float("inf"), (
        f"best_energy_eV not finite: {best!r}")
    assert res.get("best_index") == 1, (
        f"best_index {res.get('best_index')!r} is not the lowest "
        "finite-energy structure")


def test_none_index_entry_never_best(tmp_path, _clean_registry):
    res = _run_with(NONE_INDEX_RANKING, tmp_path)
    assert res["best_index"] == 2, res
    assert res["best_energy_eV"] == -1.5


def test_none_index_entry_is_warned(tmp_path, _clean_registry):
    res = _run_with(NONE_INDEX_RANKING, tmp_path)
    assert res.get("warnings"), res
