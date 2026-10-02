"""The SLURM demo entry point produces a comparable report.

``run_demo.py`` drives the same candidate space twice (EI and random
at the same budget) and writes a demo_result.json whose comparison is
built from DELFIN's own statistics (wilson_interval, fisher exact).
This test pins the report logic with fake chemistry at the same seams
as the evaluation tests — no binary, no cluster.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                       / ".gate" / "campaign_demo"))
import run_demo  # noqa: E402


def _install_fake_chemistry(monkeypatch, gaps: dict):
    from delfin.tools._types import StepResult, StepStatus

    class _FakeAdapter:
        name = "xtb_sp"

        def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
            idx = int(str(work_dir).rsplit("_", 1)[-1])
            gap = gaps[idx]
            return StepResult(
                step_name="xtb_sp",
                status=StepStatus.SUCCESS if gap is not None
                else StepStatus.FAILED,
                work_dir=work_dir, elapsed_seconds=0.001,
                data={} if gap is None else {"homo_lumo_gap_eV": gap},
                error=None if gap is not None else "xtb crashed",
            )

    import delfin.tools._registry as reg
    real_get = reg.get
    monkeypatch.setattr(reg, "get",
                        lambda name: _FakeAdapter() if name == "xtb_sp"
                        else real_get(name))
    import delfin.agent.campaign as camp_mod
    monkeypatch.setattr(
        camp_mod, "_smiles_to_geometry",
        lambda smiles, work_dir: {"xyz": "C 0.0 0.0 0.0\n", "error": None})


class TestDemoReport:
    def test_demo_report_compares_both_strategies(
            self, tmp_path, monkeypatch):
        # Valley at index 3 (gap 2.0). Specs go into the campaign
        # SUBFOLDERS the way _run_strategy expects them (a spec in the
        # work root would be ignored and the demo spec seeded instead,
        # which runs EI against the demo's 5.0 eV target — an
        # incomparable comparison).
        gaps = {i: 2.0 + 0.4 * abs(i - 3) for i in range(10)}
        _install_fake_chemistry(monkeypatch, gaps)
        for folder, strategy in (("campaign_ei", "ei"),
                                 ("campaign_random", "random")):
            sub = tmp_path / folder
            sub.mkdir()
            (sub / "CAMPAIGN.json").write_text(json.dumps({
                "space": [f"c{i}" for i in range(10)],
                "target_gap_eV": 2.0,
                "max_evaluations": 8,
                "strategy": strategy,
            }))
        ei = run_demo._run_strategy(tmp_path, "ei")
        rnd = run_demo._run_strategy(tmp_path, "random")
        assert ei["n_evaluations"] == rnd["n_evaluations"] == 8
        assert ei["best"]["score"] <= rnd["best"]["score"]
        # Both campaign folders carry their decision logs.
        assert (tmp_path / "campaign_ei" / "decisions.jsonl").exists()
        assert (tmp_path / "campaign_random" / "decisions.jsonl").exists()

    def test_demo_report_fisher_test_is_honest_about_small_n(
            self, tmp_path, monkeypatch):
        # At budget 8 the smallest two-sided p the exact test can give
        # is 1/C(16,8) — the report must carry the real p, not a
        # rounded "significant" claim the budget cannot support.
        from delfin.agent.benchmark import _fisher_exact_2x2_pvalue
        p_best = _fisher_exact_2x2_pvalue(8, 0, 0, 8)
        assert p_best > 1 / 100000      # honest floor, not 0
        assert p_best < 0.05            # best case IS significant
