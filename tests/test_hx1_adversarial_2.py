"""Adversarial tests for HX1 P3 run_round (reviewer nacht-s10)."""
import json

import pytest

from delfin.agent.harness_loop import (
    Candidate,
    RoundConfig,
    RoundReport,
    run_round,
)


def _cand(cid="c1", target="delfin/tools/guidance/tool_guidance.md"):
    return Candidate(
        kind="invariant",
        target=target,
        old_text="old line",
        new_text="new line",
        observation_task_id="task-7",
        observation="agent used the wrong workflow for the solvent request",
    )


def _config(**kw):
    base = dict(
        family_id="generalist",
        candidate_edits={"c1": _cand()},
        paired_tasks=["t1", "t2"],
        reference=[60.0, 70.0],
        probe_quality={"c1": [90.0, 100.0]},  # +30 mean, above floor-8 band
        spread=0.0,
        estimated_cost_usd=0.0,
        budget_usd=float("inf"),
    )
    base.update(kw)
    return RoundConfig(**base)


def test_empty_candidate_set_refused():
    # A round with no candidates is not a round; it must not return a
    # report that implies a measurement happened.
    with pytest.raises(ValueError):
        run_round(_config(candidate_edits={}), dry_run=True)


def test_probe_list_length_mismatch_is_loud():
    with pytest.raises(ValueError):
        run_round(_config(probe_quality={"c1": [90.0]}), dry_run=True)


def test_budget_equal_to_estimate_fits_above_refused():
    # Same "the estimate joins the day's spend before the instance
    # starts" semantics as bench_nightly: == fits, > refuses.
    report = run_round(_config(estimated_cost_usd=5.0, budget_usd=5.0),
                       dry_run=True)
    assert report["mode"] == "dry-run"
    with pytest.raises(RoundReport.BudgetRefused):
        run_round(_config(estimated_cost_usd=5.01, budget_usd=5.0),
                  dry_run=True)


def test_no_contrast_recorded_not_hidden():
    report = run_round(_config(
        probe_quality={"c1": [61.0, 71.0]}), dry_run=True)
    assert len(report["no_contrast"]) == 1
    assert report["composed"] == []
    assert report["no_contrast"][0]["mean_gain"] == pytest.approx(1.0)


def test_report_is_json_serialisable():
    report = run_round(_config(), dry_run=True)
    text = json.dumps(report)  # must not raise
    assert json.loads(text)["family_id"] == "generalist"


def test_every_decision_carries_observation_and_band():
    report = run_round(_config(), dry_run=True)
    for decision in report["composed"] + report["no_contrast"]:
        assert decision["observation_task_id"] == "task-7"
        assert decision["observation"]
        assert decision["band"] == report["noise_band"]


def test_lesson_evidence_kind_is_a_public_kind():
    # The lesson's provenance must use one of knowledge.py's public
    # evidence kinds (knowledge.py:96) — a lesson citing an invented
    # kind is unverifiable.
    from delfin.agent.knowledge import PUBLIC_EVIDENCE_KINDS
    report = run_round(_config(), dry_run=True)
    assert report["lessons"], "expected a composed lesson"
    for lesson in report["lessons"]:
        assert lesson["provenance"]["evidence_kind"] in PUBLIC_EVIDENCE_KINDS


def test_unprobed_candidate_is_loud_keyerror():
    with pytest.raises(KeyError):
        run_round(_config(candidate_edits={
            "c1": _cand(), "c2": _cand("c2")}), dry_run=True)
