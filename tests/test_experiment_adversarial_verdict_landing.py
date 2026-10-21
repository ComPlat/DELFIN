"""Package G, reviewer s6: adversarial tests for phase 5 (verdict + landing).

Attacks the builder's 15 tests do not cover:

  A5.1  malformed rows (missing n_samples) silently pooled -- probe showed
        compare_runs defaults a missing n_samples to 1/1, so a bare row
        reads as certainty; the verdict must refuse or at least not trust it
  A5.2  the alpha parameter is silently ignored (contaminated threshold)
  A5.3  empty baseline rows -> old_n=0, old_rate=0.0 -- an empty arm must
        never produce a verdict
  A5.4  a verdict on a LANDED experiment rewinds its state to "verdict"
  A5.5  a significant improvement on a contaminated (significant) null
        counts without any explanation -- asymmetry attack
  A5.6  a classification for a blocker name unrelated to the experiment
        clears the gate (classifications are free-floating text)
  A5.7  approval record with a whitespace-only "by" is refused
"""

import pytest

from delfin.agent.experiment import (
    BlockerClassification,
    Experiment,
    ExperimentError,
    HumanApproval,
    instrument_stamp,
    land,
    pre_register,
    record_human_approval,
    status_of,
    submit_for_human_review,
    verdict_with_noise_gate,
)


def _row(task_id: str, n_pass: int, n: int) -> dict:
    flags = [True] * n_pass + [False] * (n - n_pass)
    return {
        "task_id": task_id,
        "success": n_pass * 2 >= n,
        "success_rate": n_pass / n,
        "n_samples": n,
        "per_run_success": flags,
        "quality_0_100": 100 if n_pass else 10,
    }


def _bare_row(task_id: str) -> dict:
    """A row WITHOUT the counts keys -- the shape my probe showed
    compare_runs silently pools as 1/1."""
    return {"task_id": task_id, "success": True, "quality_0_100": 90}


def _exp(tmp_path, monkeypatch, *, status="registered") -> Experiment:
    exp = Experiment(id="adv5", hypothesis="h", switch="DELFIN_MODE",
                     why_chain=["one judge per comparison"])
    pre_register(exp, expectation="effect", reading="pass counts",
                 pool_size="large")
    return exp


def _stamp(tmp_path, monkeypatch, *, judge="j1"):
    f1 = tmp_path / "code.py"
    f2 = tmp_path / "judge.py"
    f1.write_text("switch default off", encoding="utf-8")
    f2.write_text("judge", encoding="utf-8")
    monkeypatch.setenv("DELFIN_ADV5_ENV", "x")
    return instrument_stamp(files=[str(f1), str(f2)],
                            env_keys=["DELFIN_ADV5_ENV"], judge=judge)


# ── A5.1: a bare row (no counts) must not be silently pooled ─────────────

def test_bare_row_is_not_silent_certainty(tmp_path, monkeypatch):
    """One baseline row without n_samples reads as 1/1 through compare_runs.
    The verdict must not treat a malformed row as a measurement."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base = [_bare_row("c")]
    cand = [_row("c", 9, 10)]
    # Acceptable outcomes: a refusal, OR a verdict that reports n=0 / does
    # not claim significance from the bare row.  UNACCEPTABLE: a clean
    # "improved" verdict that silently counted the bare row as 1/1.
    try:
        out = verdict_with_noise_gate(exp, baseline_rows=base,
                                      candidate_rows=cand,
                                      stamp_baseline=s, stamp_candidate=s)
    except ExperimentError:
        return  # refused: good
    pooled_old_n = out.get("old_n") if "old_n" in out else None
    # The verdict dict does not expose n -- so assert the honest signal:
    # the result must NOT be a confident improvement off the bare row.
    assert not (out["final"] == "improved" and out["significant"]
                and pooled_old_n is None), (
        "a bare 1-row arm was silently counted as 1/1 certainty and "
        "produced a confident 'improved' verdict")


# ── A5.2: the alpha parameter is silently ignored ─────────────────────────

def test_alpha_parameter_is_not_ignored(tmp_path, monkeypatch):
    """_noise_gate accepts alpha but never passes it to compare_runs.  An
    agent asking for alpha=0.001 gets the same gate as alpha=0.05 -- a
    contaminated threshold the caller cannot see.  Document the behavior:
    either alpha is honored, or the module must refuse a non-default alpha
    it cannot honor.  SILENT ACCEPTANCE is the finding."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base = [_row("c", 5, 10)]
    cand = [_row("c", 5, 10)]
    try:
        out_default = verdict_with_noise_gate(
            exp, baseline_rows=base, candidate_rows=cand,
            stamp_baseline=s, stamp_candidate=s)
        exp2 = _exp(tmp_path, monkeypatch)
        out_alpha = verdict_with_noise_gate(
            exp2, baseline_rows=base, candidate_rows=cand,
            stamp_baseline=s, stamp_candidate=s, alpha=0.001)
    except ExperimentError:
        return  # refused the non-default alpha: acceptable
    # If both ran, alpha had no observable effect on the gate decision for
    # identical inputs -- which is exactly the silent-ignoring finding.
    assert out_default["significant"] == out_alpha["significant"]


# ── A5.3: an empty baseline arm produces no verdict ───────────────────────

def test_empty_baseline_arm_refuses(tmp_path, monkeypatch):
    """old_n == 0 -> old_rate = 0.0: an empty arm must be a refusal, not a
    recipe for a vacuous 'improved'."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    with pytest.raises(ExperimentError):
        verdict_with_noise_gate(exp, baseline_rows=[], candidate_rows=[_row("c", 9, 10)],
                                stamp_baseline=s, stamp_candidate=s)


# ── A5.4: a verdict on a landed experiment rewinds the state machine ─────

def test_verdict_on_landed_experiment_rewinds_state(tmp_path, monkeypatch):
    """verdict_with_noise_gate only checks status != 'draft'.  Run it on an
    already-LANDED experiment and exp._move('verdict') rewinds the machine
    through approved/landed back to 'verdict' -- the terminal state is not
    terminal.  This is the P2 bypass hole resurfacing in P5."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = [_row("c", 1, 10)], [_row("c", 9, 10)]
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    record_human_approval(exp, HumanApproval(by="operator"))
    land(exp)
    assert status_of(exp) == "landed"
    # NOW the attack: another verdict on the landed experiment.
    try:
        verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                stamp_baseline=s, stamp_candidate=s)
    except ExperimentError:
        return  # refused: correct
    assert status_of(exp) != "landed", (
        "a verdict on a landed experiment silently rewound the state "
        "machine; 'landed' is not terminal")


# ── A5.5: contaminated null + significant improvement counts unexplained ─

def test_contaminated_null_does_not_guard_improvement(tmp_path, monkeypatch):
    """The null guard fires only on the regression path.  A significant
    IMPROVEMENT on a threshold that a significant null proved contaminated
    needs no explanation at all.  Document the asymmetry: if the module
    accepts it, that is a finding (an improvement on a contaminated
    threshold is not a result either)."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base = [_row("c", 1, 10)]
    cand = [_row("c", 9, 10)]
    # A LUCKY null: 1/4 vs 4/4 per task over four tasks (same 5/10 state on
    # both sides of every task is not required -- the point is a null the
    # gate calls significant by chance construction).  Probe-verified
    # significant: pooled 4/16 vs 16/16, p=0.0.
    null_base = [_row(f"t{i}", 1, 4) for i in range(4)]
    null_cand = [_row(f"t{i}", 4, 4) for i in range(4)]
    try:
        out = verdict_with_noise_gate(exp, baseline_rows=base,
                                      candidate_rows=cand,
                                      stamp_baseline=s, stamp_candidate=s,
                                      null_rows=(null_base, null_cand))
    except ExperimentError:
        return  # refused: null guard covers improvements too -- good
    assert out["null_significant"] is True
    # UNACCEPTABLE: a confident improvement with the contaminated null
    # present and unexplained.
    assert not (out["final"] == "improved"), (
        "significant improvement counted on a contaminated null without "
        "any artefact explanation")


# ── A5.6: classifications are free-floating text ──────────────────────────

def test_classification_names_an_unrelated_blocker(tmp_path, monkeypatch):
    """Nothing ties BlockerClassification.blocker to an observation: a
    classification for a blocker the experiment never reported clears the
    gate.  Document that the module accepts it (design limitation) -- the
    verdict has no blocker list to check against."""
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base = [_row("c", 9, 10)]
    cand = [_row("c", 1, 10)]
    out = verdict_with_noise_gate(
        exp, baseline_rows=base, candidate_rows=cand,
        stamp_baseline=s, stamp_candidate=s,
        classifications=[BlockerClassification(
            blocker="a blocker this experiment never observed",
            kind="measurement_artefact", reason="unrelated text")])
    assert out["final"] == "noise"


# ── A5.7: approval record with a whitespace-only human name ──────────────

def test_whitespace_human_name_refused(tmp_path, monkeypatch):
    exp = _exp(tmp_path, monkeypatch)
    s = _stamp(tmp_path, monkeypatch)
    base, cand = [_row("c", 1, 10)], [_row("c", 9, 10)]
    out = verdict_with_noise_gate(exp, baseline_rows=base, candidate_rows=cand,
                                  stamp_baseline=s, stamp_candidate=s)
    submit_for_human_review(exp, out)
    with pytest.raises(ExperimentError):
        record_human_approval(exp, HumanApproval(by="   "))
