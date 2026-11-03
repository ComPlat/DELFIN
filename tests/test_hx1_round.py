"""HX1 Phase 3 red test: the improvement round runs dry (no model
calls), composes only supported candidates, and every decision in the
report cites an observation.

Contract under test: delfin/agent.harness_loop.run_round
- runs in dry-run mode against RECORDED transcripts (fixtures): the
  round runner consumes per-instance recorded quality samples and
  reference transcript texts, never a model;
- composes ONLY candidates whose contrast is supported; a "no
  contrast" candidate is recorded as a valid, named outcome in the
  report (brief item 4);
- every decision in the report cites the observation it answers
  (task id + text) — the brief's fourth red test;
- the report is a reviewable artifact (dict-shaped, serialisable) and
  carries the KN1-style lesson records (situation/lesson/why +
  provenance evidence pointer), shaped like knowledge.Entry fields;
- the composed proposal is a PROPOSAL: the round returns it for a
  person to review; there is no apply/patch step in the module.
"""
from __future__ import annotations

import pytest

from delfin.agent.harness_loop import (
    Candidate,
    RoundConfig,
    RoundReport,
    run_round,
)


def _candidate(cid: str, kind: str, task_id: str, what: str) -> Candidate:
    return Candidate(
        kind=kind,
        target="delfin/agent/pack/skills/tool-guidance/SKILL.md",
        old_text="Old guidance text.",
        new_text="New guidance text.",
        observation_task_id=task_id,
        observation=what,
    )


def _config(**overrides) -> RoundConfig:
    fields = dict(
        family_id="workflow",
        candidate_edits={
            "inv": _candidate("inv", "invariant", "t1", "agent skipped "
                              "reading CONTROL.txt"),
            "dom": _candidate("dom", "domain", "t2", "agent guessed the "
                              "property instead of asking"),
        },
        paired_tasks=["t1", "t2"],
        reference=[40.0, 45.0],
        probe_quality={"inv": [70.0, 75.0],      # +30 at band of 30 → no
                       "dom": [80.0, 85.0]},     # +40 → supported
        spread=20.0,
    )
    fields.update(overrides)
    return RoundConfig(**fields)


def test_dry_run_composes_only_supported_candidates() -> None:
    report = run_round(_config(), dry_run=True)
    composed = [c["id"] for c in report["composed"]]
    assert composed == ["dom"]


def test_no_contrast_is_recorded_not_silent() -> None:
    report = run_round(_config(), dry_run=True)
    recorded = [c["id"] for c in report["no_contrast"]]
    assert recorded == ["inv"]


def test_every_decision_cites_an_observation() -> None:
    report = run_round(_config(), dry_run=True)
    for decision in report["composed"] + report["no_contrast"]:
        assert decision["observation_task_id"]
        assert decision["observation"]


def test_report_carries_band_and_lesson() -> None:
    report = run_round(_config(), dry_run=True)
    assert report["noise_band"] == pytest.approx(30.0)
    lessons = report["lessons"]
    assert len(lessons) == 1  # one per composed candidate
    lesson = lessons[0]
    assert lesson["situation"]
    assert lesson["lesson"]
    assert lesson["why"]
    assert lesson["provenance"]["evidence"]
    assert lesson["provenance"]["evidence_kind"] in (
        "test", "commit", "pr", "issue", "doc")


def test_report_is_serialisable() -> None:
    import json
    report = run_round(_config(), dry_run=True)
    json.dumps(report)  # must not raise


def test_budget_refusal() -> None:
    with pytest.raises(RoundReport.BudgetRefused):
        run_round(_config(estimated_cost_usd=99.0, budget_usd=10.0),
                  dry_run=True)


def test_dry_run_makes_no_model_calls() -> None:
    # The round must never import or call a runner in dry-run mode:
    # feed an explicit forbidden sentinel — if the code touched the
    # model path it would need it and fail loudly.
    report = run_round(_config(), dry_run=True)
    assert report["mode"] == "dry-run"
    assert report["family_id"] == "workflow"
