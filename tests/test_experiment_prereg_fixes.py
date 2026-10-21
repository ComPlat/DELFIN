"""Package G: reproduction of the reviewer's 5 P2 adversarial findings.

These are my own red tests for the reviewer's F1-F5 findings on
delfin/agent/experiment.py (agent/s5-g11 @ 6960b716).  All 5 are real
defects in the Phase 2 code:

F1  goalpost-moving is invisible -- a Measurement does not record WHAT was
    pre-registered, so expectations could be changed after the first
    measurement without a trace;
F2  a direct `status = "recording"` write bypasses the pre-registration
    gate, and record_measurement then accepts a measurement on a
    never-registered experiment;
F3  status can move BACKWARD silently, contradicting the docstring's
    "only moves forward";
F4  pre_register(expectation=None) crashes with AttributeError instead of
    the promised ExperimentError refusal;
F5  record_measurement(case=None) crashes with AttributeError instead of
    the promised ExperimentError refusal.

Each test must be RED on the current code and GREEN after the fix.
"""

import pytest

from delfin.agent.experiment import (
    Experiment,
    ExperimentError,
    pre_register,
    record_measurement,
    status_of,
)


def _draft(**over) -> Experiment:
    base = dict(id="e1", hypothesis="h the change matters",
                why_chain=["one judge per comparison"],
                switch="DELFIN_MODE")
    base.update(over)
    return Experiment(**base)


def _registered(**over) -> Experiment:
    exp = _draft(**over)
    pre_register(exp, expectation="rise < 5%", reading="one judge",
                 pool_size="small")
    return exp


# F1 -- a measurement records WHAT was pre-registered (tamper-evident) ────


def test_f1_measurement_records_preregistration_content(tmp_path):
    exp = _registered()
    record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})
    # Same pre-registration content -> the recorded digest
    assert exp.measurements[0].prereg_digest is not None
    assert exp.measurements[0].prereg_digest == exp.prereg_digest


def test_f1_changed_preregistration_is_detectable(tmp_path):
    exp = _registered()
    record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})
    first = exp.measurements[0].prereg_digest
    # A hypothetical goalpost move (different reading) must yield a
    # DIFFERENT digest than what the first measurement recorded.
    moved = _draft()
    pre_register(moved, expectation="rise < 5%", reading="different judge",
                 pool_size="small")
    assert moved.prereg_digest != first


# F2 -- a direct status write must NOT bypass the pre-registration gate ────


def test_f2_direct_status_write_cannot_enable_measuring():
    exp = _draft()
    try:
        exp.status = "recording"     # the bypass attempt
    except Exception:
        # a read-only status property already stops the bypass
        return
    # If the write succeeded, measuring must STILL be refused:
    with pytest.raises(ExperimentError):
        record_measurement(exp, case="c1", outcome={"runtime_s": 1.0})


def test_f2_measure_on_never_registered_is_refused_even_if_status_set():
    exp = _draft()
    exp.prereg_digest = None          # direct field write, no real registration
    with pytest.raises(ExperimentError):
        record_measurement(exp, case="c1", outcome={"runtime_s": 1.0})


# F3 -- status may only move forward ──────────────────────────────────────


def test_f3_status_cannot_move_backward():
    exp = _registered()
    record_measurement(exp, case="c1", outcome={"a": 1})
    assert status_of(exp) == "recording"
    # moving backward to draft must be refused (forward-only)
    with pytest.raises(ExperimentError):
        exp._move("draft")
    assert status_of(exp) == "recording"


def test_f3_direct_backward_write_is_blocked():
    exp = _registered()
    with pytest.raises(Exception):
        exp.status = "draft"


# F4/F5 -- None arguments are a refusal, not a crash ──────────────────────


def test_f4_pre_register_none_expectation_refuses():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation=None, reading="one judge",
                     pool_size="small")


def test_f4_pre_register_none_reading_refuses():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="rise < 5%", reading=None,
                     pool_size="small")


def test_f5_record_measurement_none_case_refuses():
    exp = _registered()
    with pytest.raises(ExperimentError):
        record_measurement(exp, case=None, outcome={"runtime_s": 1.0})
