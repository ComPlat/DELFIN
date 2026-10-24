"""R3 probe (reviewer nacht-s25): best_* must only come from int-index entries.

Real call path: MlpScreenWorkflow.run() -> registry fake mlp_rank adapter
-> run() picks best over the "finite" filter.  The filter checks only
math.isfinite(energy_eV); an entry with a non-int "index" (e.g. None)
and a finite energy is then picked by min() and reported as best_index.
Expected: an index-invalid entry is treated like a non-finite one
(counted in the warning channel) and best_* comes only from entries
with an int index.  Run through the gate against UNFIXED 395696da: red.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.tools._types import StepResult, StepStatus  # noqa: E402
from delfin.workflows.contrib.w21_mlp_screen import (  # noqa: E402
    MlpScreenWorkflow,
)


class _FakeRankAdapter:
    name = "mlp_rank"
    description = "probe fake"

    def validate_params(self, backend="ani2x", charge=0, mult=1,
                        ensemble=None):
        assert ensemble, "validate_params must receive the ensemble path"

    def __init__(self, ranking):
        self._ranking = ranking

    def execute(self, work_dir, *, ensemble=None, backend="ani2x",
                charge=0, mult=1, device="cpu"):
        out = StepResult(
            step_name="mlp_rank", status=StepStatus.SUCCESS,
            data={"ranking": self._ranking,
                  "n_structures": len(self._ranking)},
        )
        (work_dir / "ranking.json").write_text(json.dumps(self._ranking))
        return out


def _run_with(ranking, tmp_path):
    from delfin.tools import _registry as reg

    reg._discover_builtin_adapters()
    reg.register(_FakeRankAdapter(ranking))
    return MlpScreenWorkflow().run(
        config={}, ensemble="ens.xyz",
        work_dir=str(tmp_path / "wd"))


def test_none_index_entry_never_best(tmp_path):
    res = _run_with([
        {"index": None, "energy_eV": -99.0},
        {"index": 2, "energy_eV": -1.5},
    ], tmp_path)
    assert res["best_index"] == 2, res
    assert res["best_energy_eV"] == -1.5


def test_none_index_entry_is_warned(tmp_path):
    res = _run_with([
        {"index": None, "energy_eV": -99.0},
        {"index": 2, "energy_eV": -1.5},
    ], tmp_path)
    assert res.get("warnings"), res
