"""The campaign is a registered DELFIN workflow, drivable headless.

The workflow is DELFIN's public surface for the campaign: registered
under ``campaign`` like every contrib workflow (imag, co2, esd, …),
running the loop with a named search space and target, writing its
decision log and summary into the campaign folder, and returning the
summary so a caller (CLI, dashboard, agent turn) can report it.  The
SLURM demo path is pinned too: the campaign folder contains a
``CONTROL.txt``-companion (``CAMPAIGN.json``) describing space, target
and budget, and the summary lands in ``CAMPAIGN_summary.json``.

The workflow runs WITHOUT a language model: its decisions come from
the acquisition loop; the agent is only woken on failure, stagnation
or budget exhaustion via the scheduler contract tested earlier.
"""

import json

import pytest

from delfin.workflows.registry import get as get_workflow


# ---------------------------------------------------------------------------
# registration and contract
# ---------------------------------------------------------------------------


class TestRegistration:
    def test_campaign_workflow_is_registered(self):
        wf = get_workflow("campaign")
        assert wf is not None
        assert wf.name == "campaign"
        assert wf.description

    def test_run_cli_rejects_an_unknown_search_space_file(self, tmp_path):
        wf = get_workflow("campaign")
        rc = wf.run_cli([str(tmp_path / "missing" / "CAMPAIGN.json")])
        assert rc != 0


# ---------------------------------------------------------------------------
# run with fake chemistry (same seams as the evaluation tests)
# ---------------------------------------------------------------------------


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


class TestRun:
    def test_run_writes_log_and_summary_into_the_campaign_folder(
            self, tmp_path, monkeypatch):
        gaps = {i: 2.0 + 0.3 * (i - 4) ** 2 / 4 for i in range(6)}
        _install_fake_chemistry(monkeypatch, gaps)
        spec = {
            "space": [f"cand{i}" for i in range(6)],
            "target_gap_eV": 2.0,
            "max_evaluations": 4,
        }
        (tmp_path / "CAMPAIGN.json").write_text(json.dumps(spec))
        wf = get_workflow("campaign")
        result = wf.run(config={}, campaign_dir=str(tmp_path))
        assert result["n_evaluations"] == 4
        assert result["stopped_by"] in ("budget", "stagnation")
        assert (tmp_path / "decisions.jsonl").exists()
        summary = json.loads(
            (tmp_path / "CAMPAIGN_summary.json").read_text())
        assert summary["target_gap_eV"] == 2.0
        assert summary["n_evaluations"] == 4
        assert summary["best"]["value"] is not None

    def test_summary_ranks_the_best_candidate_first(
            self, tmp_path, monkeypatch):
        # A valley at index 3: the campaign's best must be that one.
        gaps = {i: 2.0 + 0.5 * abs(i - 3) for i in range(8)}
        _install_fake_chemistry(monkeypatch, gaps)
        spec = {
            "space": [f"cand{i}" for i in range(8)],
            "target_gap_eV": 2.0,
            "max_evaluations": 6,
        }
        (tmp_path / "CAMPAIGN.json").write_text(json.dumps(spec))
        wf = get_workflow("campaign")
        result = wf.run(config={}, campaign_dir=str(tmp_path))
        best = result["best"]
        assert best["score"] <= min(
            abs(v - 2.0) for v in
            [2.0 + 0.5 * abs(i - 3) for i in range(8)])
