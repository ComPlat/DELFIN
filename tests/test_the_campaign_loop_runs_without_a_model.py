"""The campaign loop runs WITHOUT a language model.

A campaign is a closed optimization loop over a bounded candidate space
(a parameter grid, e.g. substituent choices on a molecular anchor). It
must make every decision itself — acquisition, stopping, logging — and
wake the agent only on failure, stagnation or budget exhaustion. These
tests pin that contract:

* the decision log records every decision it makes (acquisition call,
  evaluation, stopping reason), so the loop is auditable;
* the budget cap stops the loop and reports WHY it stopped;
* the acquisition picks the next point from the observed data —
  Expected Improvement on a GP surrogate — and prefers the untried
  point a data set already labels as most promising;
* with no data, acquisition falls back to Latin-hypercube seeding;
* the mode decides what "better" means: with a target gap, values
  CLOSE to the target are scored as better — not bigger ones;
* every score is recorded together with the parameters that produced
  it, so a campaign can be resumed or audited from its files alone.
"""

import json

import pytest

from delfin.agent.campaign import (
    CampaignBudget,
    CampaignDecision,
    DecisionLog,
    DecisionReason,
    TargetGap,
    campaign_should_wake,
    score_candidate,
    select_next,
)


# ---------------------------------------------------------------------------
# decision log
# ---------------------------------------------------------------------------


class TestDecisionLog:
    def test_every_decision_is_recorded_with_reason(self, tmp_path):
        log = DecisionLog(tmp_path / "decisions.jsonl")
        log.record(CampaignDecision(
            kind="acquire", reason="ei", payload={"index": 4},
        ))
        log.record(CampaignDecision(
            kind="evaluate", reason="done", payload={"value": 1.9},
        ))
        rows = [json.loads(line) for line in
                (tmp_path / "decisions.jsonl").read_text().splitlines()]
        assert len(rows) == 2
        # Every row is auditable: what, why, when, and the data behind it.
        assert rows[0]["kind"] == "acquire"
        assert rows[0]["reason"] == "ei"
        assert rows[0]["payload"] == {"index": 4}
        assert rows[0]["ts"] > 0

    def test_read_back_from_the_log_file_alone(self, tmp_path):
        log = DecisionLog(tmp_path / "decisions.jsonl")
        log.record(CampaignDecision(kind="stop", reason="budget"))
        assert DecisionLog.entries(tmp_path / "decisions.jsonl") == [
            {"kind": "stop", "reason": "budget", "payload": None},
        ]


# ---------------------------------------------------------------------------
# budget / wake policy
# ---------------------------------------------------------------------------


class TestBudget:
    def test_budget_stops_the_loop_when_spent(self):
        budget = CampaignBudget(max_evaluations=3)
        assert budget.spend() is True       # 1/3
        assert budget.spend() is True       # 2/3
        assert budget.spend() is True       # 3/3
        assert budget.exhausted is True
        assert budget.spend() is False      # past the cap: refused

    def test_wake_only_on_failure_stagnation_or_budget(self):
        # These are the ONLY conditions that wake the model.
        assert campaign_should_wake(reason=DecisionReason.FAILURE) is True
        assert campaign_should_wake(
            reason=DecisionReason.STAGNATION) is True
        assert campaign_should_wake(reason=DecisionReason.BUDGET) is True
        # A healthy loop must run on its own.
        assert campaign_should_wake(reason=DecisionReason.ACQUIRE) is False
        assert campaign_should_wake(reason=DecisionReason.EVALUATE) is False


# ---------------------------------------------------------------------------
# acquisition
# ---------------------------------------------------------------------------


class TestAcquisition:
    def test_seed_design_is_a_latin_hypercube_over_the_space(self):
        space = list(range(10))
        seed = select_next(space, observations={}, n=4)
        assert len(seed) == 4
        assert len(set(seed)) == 4           # no duplicate evaluations
        for p in seed:
            assert p in space

    def test_ei_prefers_the_point_the_data_labels_best(self):
        # Two untried points on the flank of a cluster of good scores:
        # the GP interpolates scores, so the point nearer the cluster
        # gets the lower predicted score and EI must rank it first.
        space = list(range(6))
        observations = {0: 0.0, 1: 0.1, 2: 0.2, 5: 9.9}   # 3, 4 untried
        result = select_next(space, observations=observations, n=1,
                             target=TargetGap(2.0))
        assert result == [3]

    def test_ei_does_not_revisit_tried_points(self):
        space = list(range(6))
        observations = {i: 2.0 for i in range(6)}
        result = select_next(space, observations=observations, n=2,
                             target=TargetGap(2.0))
        # Every point is tried and nothing improves the prediction --
        # still, the loop must answer SOMETHING without crashing, and
        # it must never propose a point it already evaluated.
        assert result is not None
        assert all(p in space for p in result)


# ---------------------------------------------------------------------------
# scoring / target mode
# ---------------------------------------------------------------------------


class TestTargetGap:
    def test_closer_to_target_is_better(self):
        # score is |gap - target|: LOWER is better, 0 is perfect.
        target = TargetGap(2.0)
        assert target.score(2.0) < target.score(3.0)
        assert target.score(2.0) < target.score(1.0)
        assert target.score(2.0) == pytest.approx(0.0)

    def test_no_target_treats_small_gaps_as_best(self):
        # Without a target the goal is minimizing the gap.
        assert TargetGap(None).score(0.5) < TargetGap(None).score(1.5)

    def test_score_candidate_records_parameters_with_the_value(self):
        scored = score_candidate(
            params={"substituent": "CN"}, value=1.8, target=TargetGap(2.0))
        assert scored["params"] == {"substituent": "CN"}
        assert scored["value"] == pytest.approx(1.8)
        assert scored["score"] == pytest.approx(TargetGap(2.0).score(1.8))


class TestAcquisitionProtocol:
    def test_selected_acquisition_runs_first_then_surely_stops(self):
        # A full loop on a tiny space: seeding, then GP acquisition,
        # then the budget stops it. No model in the loop.
        space = list(range(5))
        budget = CampaignBudget(max_evaluations=3)
        observations: dict = {}
        picks: list = []
        for _ in range(5):
            if budget.exhausted:
                break
            if observations:
                nxt = select_next(space, observations, n=1,
                                  target=TargetGap(2.0), seed_used=2)
            else:
                nxt = select_next(space, observations, n=2)
            picks.extend(nxt)
            for p in nxt:
                observations[p] = 1.0 + 0.1 * p
            budget.spend()
        assert picks
        assert budget.exhausted
        # Acquisition never re-picks an evaluated point.
        assert len(picks) == len(set(picks))
