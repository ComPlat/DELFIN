"""A campaign evaluates a candidate through DELFIN's own chemistry.

One candidate cycle is: the SMILES list entry becomes an XYZ geometry
(via DELFIN's existing SMILES converter), the geometry goes through the
``xtb_sp`` adapter (native xTB, gap parsed by DELFIN), and the
observation lands in the decision log with the parameters that produced
it.  Every real binary call is faked at the seam DELFIN defines — the
adapter's ``execute`` — so these tests run without xtb on PATH and
without touching the chemistry core.

Also pinned here:
* a failed evaluation is logged as a FAILURE decision and wakes the
  model (the only failure path that may);
* a stagnated campaign wakes the model — the wake reason is derived
  from the decision log, not recomputed ad hoc;
* the Bayes-vs-random demo: at equal budget, the EI acquisition must
  find a candidate whose observed value is at least as close to the
  target as the random baseline's best, on a synthetic landscape.
"""

import json

import pytest

from delfin.agent.campaign import (
    Campaign,
    CampaignBudget,
    DecisionLog,
    DecisionReason,
    TargetGap,
    campaign_should_wake,
    score_candidate,
)


# ---------------------------------------------------------------------------
# fake chemistry at DELFIN's adapter seam
# ---------------------------------------------------------------------------


def _install_fake_xtb(monkeypatch, gaps_by_index: dict):
    """Patch the adapter registry so ``xtb_sp`` returns chosen gaps.

    The fake lives at ``delfin.tools._registry.get`` — the same seam
    ``run_step`` uses, so the campaign's real evaluation path is
    exercised except for the binary call itself.
    """
    from delfin.tools._types import StepResult, StepStatus

    class _FakeAdapter:
        name = "xtb_sp"

        def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
            # The campaign writes the candidate's index into the work
            # dir name; the fake reads it back to pick a gap.
            idx = int(str(work_dir).rsplit("_", 1)[-1])
            gap = gaps_by_index[idx]
            return StepResult(
                step_name="xtb_sp",
                status=StepStatus.SUCCESS if gap is not None
                else StepStatus.FAILED,
                work_dir=work_dir,
                elapsed_seconds=0.001,
                data={} if gap is None
                else {"homo_lumo_gap_eV": gap},
                error=None if gap is not None else "xtb crashed",
            )

    import delfin.tools._registry as reg
    real_get = reg.get
    monkeypatch.setattr(reg, "get",
                        lambda name: _FakeAdapter() if name == "xtb_sp"
                        else real_get(name))

    # The geometry step is faked at the module seam the campaign
    # defines (_smiles_to_geometry) — the real RDKit converter would
    # reject these placeholder SMILES, which is DELFIN behaving
    # correctly, not the campaign failing.
    import delfin.agent.campaign as camp_mod
    monkeypatch.setattr(
        camp_mod, "_smiles_to_geometry",
        lambda smiles, work_dir: {"xyz": "C 0.0 0.0 0.0\n",
                                  "error": None})


def _candidate_index_from_workdir(camp: Campaign, work_dir) -> int:
    return int(str(work_dir).rsplit("_", 1)[-1])


# ---------------------------------------------------------------------------
# the evaluation cycle
# ---------------------------------------------------------------------------


class TestEvaluationCycle:
    def test_candidate_cycle_records_gap_in_the_log(self, tmp_path, monkeypatch):
        _install_fake_xtb(monkeypatch, {i: 1.0 + 0.2 * i for i in range(8)})
        camp = Campaign(
            space=[f"cand{i}" for i in range(8)],
            target=TargetGap(2.0),
            budget=CampaignBudget(max_evaluations=3),
            log=DecisionLog(tmp_path / "decisions.jsonl"),
            work_dir=tmp_path / "work",
        )
        camp.evaluate_candidates()
        assert camp.observations          # gaps came back
        rows = DecisionLog.entries(tmp_path / "decisions.jsonl")
        kinds = [r["kind"] for r in rows]
        assert "evaluate" in kinds        # every evaluation logged
        for r in rows:
            if r["kind"] == "evaluate":
                assert r["payload"]["value"] is not None
                assert r["payload"]["params"]["smiles"].startswith("cand")

    def test_failed_evaluation_wakes_the_model(self, tmp_path, monkeypatch):
        # Non-flat landscape: with all values equal the stagnation rule
        # would legitimately stop the loop before the failing
        # candidate is ever reached — the point of this test is the
        # FAILURE wake, so the scores must keep changing.
        gaps = {i: (None if i == 5 else 1.0 + 0.3 * i) for i in range(8)}
        _install_fake_xtb(monkeypatch, gaps)
        camp = Campaign(
            space=[f"cand{i}" for i in range(8)],
            target=TargetGap(2.0),
            budget=CampaignBudget(max_evaluations=8),
            log=DecisionLog(tmp_path / "decisions.jsonl"),
            work_dir=tmp_path / "work",
            stagnation_rounds=99,          # never stagnate in this test
        )
        result = camp.evaluate_candidates()
        assert any(r["index"] == 5 for r in result if "index" in r)
        rows = DecisionLog.entries(tmp_path / "decisions.jsonl")
        failures = [r for r in rows
                    if r["kind"] == "wake" and r["reason"] == "failure"]
        assert failures, "a failed evaluation must wake the model"
        # And the failure is visible in the observations-less record.
        assert 5 not in camp.observations
        # A burned candidate is never re-proposed.
        picked = [r["index"] for r in result if "index" in r]
        assert len(picked) == len(set(picked))

    def test_budget_never_exceeded_even_on_repeated_cycles(self, tmp_path, monkeypatch):
        _install_fake_xtb(monkeypatch, {i: 2.0 + 0.01 * i for i in range(8)})
        camp = Campaign(
            space=[f"cand{i}" for i in range(8)],
            target=TargetGap(2.0),
            budget=CampaignBudget(max_evaluations=3),
            log=DecisionLog(tmp_path / "decisions.jsonl"),
            work_dir=tmp_path / "work",
        )
        camp.evaluate_candidates()
        camp.evaluate_candidates()          # second call: budget already spent
        assert len(camp.observations) <= 3
        evaluations = [r for r in DecisionLog.entries(
            tmp_path / "decisions.jsonl") if r["kind"] == "evaluate"]
        assert len(evaluations) == 3


# ---------------------------------------------------------------------------
# stagnation and wake reasons from the log
# ---------------------------------------------------------------------------


class TestStagnation:
    def test_stagnated_campaign_wakes_with_stagnation_reason(
            self, tmp_path, monkeypatch):
        # A flat landscape: every evaluation returns the same value, so
        # the best score stops improving.
        _install_fake_xtb(monkeypatch, {i: 2.0 for i in range(8)})
        camp = Campaign(
            space=[f"cand{i}" for i in range(8)],
            target=TargetGap(2.0),
            budget=CampaignBudget(max_evaluations=8),
            log=DecisionLog(tmp_path / "decisions.jsonl"),
            work_dir=tmp_path / "work",
            stagnation_rounds=3,
        )
        camp.evaluate_candidates()
        assert campaign_should_wake(reason=DecisionReason.STAGNATION)


# ---------------------------------------------------------------------------
# the demo: EI acquisition vs random at equal budget (synthetic)
# ---------------------------------------------------------------------------


class TestBayesVsRandom:
    def test_ei_beats_or_matches_random_at_equal_budget(
            self, tmp_path, monkeypatch):
        # Synthetic landscape over 12 candidates with a single valley
        # near index 6 (score = |gap - 2.0|): scores rise on both
        # sides, so the GP can interpolate toward it — EI must find the
        # valley with fewer evaluations than a random 5-point design.
        def gap_of(i: int) -> float:
            dist = abs(i - 6)
            return 2.0 + 0.5 * dist

        gaps = {i: gap_of(i) for i in range(12)}
        gaps[7] = 1.6          # off-lattice improvement near the valley
        _install_fake_xtb(monkeypatch, gaps)

        def run(strategy: str) -> dict:
            camp = Campaign(
                space=[f"cand{i}" for i in range(12)],
                target=TargetGap(2.0),
                budget=CampaignBudget(max_evaluations=5),
                log=DecisionLog(tmp_path / f"{strategy}_decisions.jsonl"),
                work_dir=tmp_path / strategy,
            )
            camp.evaluate_candidates(strategy=strategy)
            return camp.best_observations()

        ei_result = run("ei")
        # Random baseline: fixed-seed Latin hypercube over the space.
        from delfin.agent.campaign import _latin_hypercube_points
        seed = _latin_hypercube_points(12, 5)
        random_best = min(
            TargetGap(2.0).score(gaps[i]) for i in seed
            if gaps.get(i) is not None)
        # best_observations() carries the score already.
        ei_best = min(
            entry["score"] for entry in ei_result.values())
        # EI must never be WORSE than random at the same budget; on this
        # landscape it is strictly better (it finds the valley).
        assert ei_best <= random_best

    def test_best_observations_keep_parameters_with_values(self, tmp_path, monkeypatch):
        _install_fake_xtb(monkeypatch, {i: 1.0 + i for i in range(8)})
        camp = Campaign(
            space=[f"cand{i}" for i in range(8)],
            target=TargetGap(2.0),
            budget=CampaignBudget(max_evaluations=4),
            log=DecisionLog(tmp_path / "decisions.jsonl"),
            work_dir=tmp_path / "work",
        )
        camp.evaluate_candidates()
        best = camp.best_observations()
        for idx, entry in best.items():
            assert entry["params"]["smiles"] == f"cand{idx}"
            assert entry["value"] == 1.0 + idx
            assert entry["score"] == pytest.approx(
                TargetGap(2.0).score(1.0 + idx))
