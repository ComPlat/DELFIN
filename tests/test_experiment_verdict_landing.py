"""Package G, Phase 5: the noise-gated verdict and the landing state machine.

A verdict "better" or "regressed" counts only when the effect clears the
noise gate (package D's ``compare_runs.significant``, read as the SINGLE
instrument -- not a threshold-only verdict string), every blocker is
classified (real regression vs measurement artefact, with reason) BEFORE
the verdict counts, and the experiment lands only through an explicit
human approval record.

The null run is part of this: a null comparison (same state twice) that
comes back significant by pure chance is exactly the spurious signal the
noise gate must not report as a real regression.  Landing is a state
machine -- verdict -> human review -> approval -> land -- whose ONLY path to
"approved" is a human approval record; an agent message never approves, a
double approval record is refused, and landing without approval is
impossible.
"""

import pytest

from delfin.agent.benchmark import compare_runs

from delfin.agent.experiment import (
    BlockerClassification,
    Experiment,
    ExperimentError,
    HumanApproval,
    InstrumentStamp,
    assert_same_stamp,
    can_land,
    instrument_stamp,
    land,
    pre_register,
    record_human_approval,
    record_measurement,
    status_of,
    submit_for_human_review,
    verdict_with_noise_gate,
)


def _row(task_id: str, n_pass: int, n: int) -> dict:
    """One result row in the shape compare_runs consumes."""
    flags = [True] * n_pass + [False] * (n - n_pass)
    return {
        "task_id": task_id,
        "success": n_pass * 2 >= n,
        "success_rate": n_pass / n,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 100 if n_pass else 10,
    }


def _exp(tmp_path, monkeypatch, *, status_ok=True) -> Experiment:
    import os
    import time as _t
    f = tmp_path / "code.py"
    f.write_text("switch default off", encoding="utf-8")
    monkeypatch.setenv("DELFIN_EXP_P5_ENV", "x")
    exp = Experiment(id="m5", hypothesis="h", switch="DELFIN_MODE",
                     why_chain=["one judge per comparison"])
    if status_ok:
        pre_register(exp, expectation="effect", reading="pass counts",
                     pool_size="large")
    return exp


def _stamp(tmp_path, monkeypatch, *, judge="j1", text="switch default off"):
    f1 = tmp_path / "code.py"
    f2 = tmp_path / "judge.py"
    f1.write_text(text, encoding="utf-8")
    f2.write_text("judge", encoding="utf-8")
    monkeypatch.setenv("DELFIN_EXP_P5_ENV", "x")
    return instrument_stamp(files=[str(f1), str(f2)],
                            env_keys=["DELFIN_EXP_P5_ENV"], judge=judge)


def _regression(exp, save=True):
    """A significant regression: 1/10 baseline vs 9/10 candidate (success
    FALLS) -- Fisher separates it; candidate < baseline."""
    base = [_row("c", 9, 10)]
    cand = [_row("c", 1, 10)]
    return base, cand


def _improvement(exp, save=True):
    base = [_row("c", 1, 10)]
    cand = [_row("c", 9, 10)]
    return base, cand


def _noise():
    base = [_row("c", 5, 10)]
    cand = [_row("c", 5, 10)]
    return base, cand


# ── Verdict refuses a cross-stamp or non-pre-registered experiment ───────


def test_verdict_refuses_different_stamps(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s1 = _stamp(tmp_path, monkeypatch, judge="v1")
    monkeypatch.setenv("DELFIN_EXP_P5_ENV", "x")
    s2 = _stamp(tmp_path, monkeypatch, judge="v2")
    base, cand = _noise()
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                stamp_baseline=s1, stamp_candidate=s2)


def test_verdict_refuses_unregistered_experiment(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch, status_ok=False)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _noise()
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                stamp_baseline=s, stamp_candidate=s)


# ── The noise gate is one instrument: significant, not a threshold string ─


def test_effect_inside_noise_is_not_better(tmp_path, monkeypatch):
    # 5/10 vs 5/10 is not significant; the verdict must NOT say "improved".
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _noise()
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    assert out["significant"] is False
    assert out["final"] == "noise"


def test_significant_improvement_is_improved(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    assert out["significant"] is True
    assert out["effect"] == "better"
    assert out["final"] == "improved"


# ── A regression needs a blocker classification before it counts ─────────


def test_regression_without_classification_refuses(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                stamp_baseline=s, stamp_candidate=s)


def test_regression_with_real_classification_counts(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    out = verdict_with_noise_gate(
        exp, baseline_rows=base, candidate_rows=cand,
        stamp_baseline=s, stamp_candidate=s,
        classifications=[BlockerClassification(
            blocker="quality dropped on every task",
            kind="real_regression", reason="reproduced without the switch on")])
    assert out["effect"] == "regression"
    assert out["final"] == "regressed"


def test_regression_classified_as_artefact_is_noise(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    out = verdict_with_noise_gate(
        exp, baseline_rows=base, candidate_rows=cand,
        stamp_baseline=s, stamp_candidate=s,
        classifications=[BlockerClassification(
            blocker="node went offline",
            kind="measurement_artefact", reason="same input re-ran clean")])
    assert out["final"] == "noise"


def test_bad_classification_kind_refuses(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(
            exp, baseline_rows=base, candidate_rows=cand,
            stamp_baseline=s, stamp_candidate=s,
            classifications=[BlockerClassification(
                blocker="b", kind="who_knows", reason="r")])


def test_classification_without_reason_refuses(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(
            exp, baseline_rows=base, candidate_rows=cand,
            stamp_baseline=s, stamp_candidate=s,
            classifications=[BlockerClassification(
                blocker="b", kind="real_regression", reason="")])


# ── A lucky-significant NULL run is not a real regression ────────────────


def test_null_run_significant_does_not_count_as_regression(tmp_path, monkeypatch):
    # A null comparison (the SAME state measured twice) that comes back
    # significant by chance (~5% at alpha 0.05) means the thresholds are
    # contaminated.  A main "regression" on top of a lucky-significant null
    # must NOT be reported as a real regression.
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    # The two null arms measure the SAME state, but noise produced 10/10 vs
    # 0/10 on the shared case -- a significant null by pure chance.
    null_a = [_row("c_null", 10, 10)]
    null_b = [_row("c_null", 0, 10)]
    with pytest.raises(ExperimentError):
        # no classification clearing the null -> refuse to count the
        # regression as real (would be reporting a lucky-significant null)
        verdict_with_noise_gate(
            exp, baseline_rows=base, candidate_rows=cand,
            stamp_baseline=s, stamp_candidate=s,
            null_rows=(null_a, null_b))


def test_null_run_significant_cleared_by_artefact_then_regression_counts(tmp_path, monkeypatch):
    # The lucky-significant null is EXPLAINED as a measurement artefact (the
    # null pool was contaminated), so the noise gate is trustworthy again and
    # the MAIN regression may be judged -- but it still needs its own
    # real_regression classification.
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _regression(exp)
    null_a = [_row("c_null", 10, 10)]
    null_b = [_row("c_null", 0, 10)]
    out = verdict_with_noise_gate(
        exp, baseline_rows=base, candidate_rows=cand,
        stamp_baseline=s, stamp_candidate=s,
        null_rows=(null_a, null_b),
        classifications=[
            BlockerClassification(
                blocker="null arm shows a spurious gap",
                kind="measurement_artefact",
                reason="the null pool was re-measured under a node that "
                       "re-ran a full pool; the lucky-significant null is a "
                       "chance draw"),
            BlockerClassification(
                blocker="quality dropped on every main task",
                kind="real_regression",
                reason="reproduced without the switch on"),
        ])
    assert out["null_significant"] is True
    assert out["final"] == "regressed"


# ── Landing: only a human approval record leads to "approved" / "landed" ─


def test_landing_refuses_without_human_approval(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    assert can_land(exp) is False
    with pytest.raises(ExperimentError):
        land(exp)
    assert status_of(exp) == "human_review"


def test_landing_requires_human_approval_record(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    record_human_approval(exp, HumanApproval(by="operator"))
    assert status_of(exp) == "approved"
    assert can_land(exp) is True
    land(exp)
    assert status_of(exp) == "landed"


def test_agent_message_is_not_approval(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    # An approval that is not a HumanApproval record is a messagelike "OK".
    with pytest.raises(ExperimentError):
        record_human_approval(exp, "OK from another agent")
    with pytest.raises(ExperimentError):
        record_human_approval(exp, {"by": "agent", "kind": "message"})
    assert status_of(exp) == "human_review"


def test_double_approval_record_is_refused(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    record_human_approval(exp, HumanApproval(by="operator"))
    assert status_of(exp) == "approved"
    # a second approval record is refused -- can't double-record
    with pytest.raises(ExperimentError):
        record_human_approval(exp, HumanApproval(by="operator"))
    land(exp)
    assert status_of(exp) == "landed"
    # landing twice is refused -- cannot double-land
    with pytest.raises(ExperimentError):
        land(exp)


# ---------------------------------------------------------------------------
# Reviewer findings A5.3 + A5.5 (P5 adversarial, from nacht-s6's
# tests/test_experiment_adversarial_verdict_landing.py @ 0e72df4e).  Both
# must be RED on the current code and GREEN after the fix.
# ---------------------------------------------------------------------------

def test_a53_empty_baseline_arm_is_a_refusal_not_a_verdict(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    _, cand = _improvement(exp)
    # an EMPTY baseline arm must refuse, not produce a verdict on old_n=0
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=[], candidate_rows=cand,
                                stamp_baseline=s, stamp_candidate=s)


def test_a55_significant_improvement_on_a_contaminated_null_is_refused(
        tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = _improvement(exp)   # significant improvement
    # null run of the SAME state that comes back significant by chance
    # (pooled 4/16 vs 16/16, p ~ 0.0006) contaminates the thresholds; no
    # artefact classification -> the improvement must be refused too.
    null_a = [_row("c", 4, 16)]
    null_b = [_row("c", 16, 16)]
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                stamp_baseline=s, stamp_candidate=s,
                                null_rows=(null_a, null_b))
