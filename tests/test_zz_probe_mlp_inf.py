"""R1 probe: mlp_screen reports a silent inf best_energy_eV.

Real call path: MlpScreenWorkflow.run() -> delfin.tools._registry.get("mlp_rank")
-> adapter.execute(wd, ...).  The fake adapter is contract-faithful:
StepResult with data["ranking"] as the real mlp_rank adapter produces it
(structures whose MLP evaluation crashed come back with float("inf") —
delfin/mlp_tools/calculators.py:286-290).  The workflow must not report
an inf best_energy_eV silently: it must warn, and best_index must point
at a structure that actually produced an energy.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.tools._registry import register  # noqa: E402
from delfin.tools._types import StepResult, StepStatus  # noqa: E402
from delfin.workflows.contrib.w21_mlp_screen import (  # noqa: E402
    MlpScreenWorkflow,
)


class _FakeRankAdapter:
    """Contract-faithful mlp_rank fake: one crashed structure -> inf."""

    name = "mlp_rank"
    description = "probe fake"

    def validate_params(self, backend="ani2x", charge=0, mult=1, ensemble=None):
        assert ensemble, "validate_params must receive the ensemble path"

    def execute(self, work_dir, *, ensemble=None, backend="ani2x",
                charge=0, mult=1, device="cpu"):
        ranking = [
            {"index": 0, "energy_eV": float("inf")},
            {"index": 1, "energy_eV": -12.34},
            {"index": 2, "energy_eV": -10.0},
        ]
        out = StepResult(
            step_name="mlp_rank", status=StepStatus.SUCCESS,
            data={"ranking": ranking, "n_structures": 3},
        )
        (work_dir / "ranking.json").write_text(json.dumps(ranking))
        return out


def test_all_finite_ranking_does_not_warn(tmp_path):
    """Negative control: a clean ranking must not gain a warning."""
    from delfin.tools import _registry as reg

    reg._discover_builtin_adapters()

    class _CleanAdapter(_FakeRankAdapter):
        def execute(self, work_dir, **kwargs):
            out = StepResult(
                step_name="mlp_rank", status=StepStatus.SUCCESS,
                data={"ranking": [{"index": 0, "energy_eV": -5.0},
                                  {"index": 1, "energy_eV": -7.5}],
                      "n_structures": 2},
            )
            return out

    reg.register(_CleanAdapter())
    res = MlpScreenWorkflow().run(config={}, ensemble="ens.xyz",
                                  work_dir=str(tmp_path / "wd"))
    assert res.get("warnings") is None
    assert res["best_index"] == 1 and res["best_energy_eV"] == -7.5

def test_silent_inf_is_reported_and_best_is_finite(tmp_path):
    # Register AFTER the workflow triggers lazy discovery, so the real
    # mlp_rank adapter (delfin/tools/adapters/mlp.py:324) does not
    # overwrite the fake.
    from delfin.tools import _registry as reg
    reg._discover_builtin_adapters()
    reg.register(_FakeRankAdapter())
    wf = MlpScreenWorkflow()
    res = wf.run(config={}, ensemble="ens.xyz",
                 work_dir=str(tmp_path / "wd"))

    assert res.get("warnings"), (
        f"no warning on an inf ranking energy: "
        f"{json.dumps({k: res.get(k) for k in res})}")
    best = res.get("best_energy_eV")
    assert best is not None and best == best and best != float("inf"), (
        f"best_energy_eV not finite: {best!r}")
    assert res.get("best_index") == 1, (
        f"best_index {res.get('best_index')!r} is not the lowest "
        "finite-energy structure")
