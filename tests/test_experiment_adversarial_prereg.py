"""Package G, reviewer s6: adversarial tests for phase 2 (Experiment +
pre-registration). These try to break the gate: move the goalposts after
registration, bypass the state machine, or push the refusals off the
ExperimentError channel.
"""


import pytest

from delfin.agent.experiment import (
    Experiment,
    ExperimentError,
    _prereg_digest,
    pre_register,
    record_measurement,
)


def _draft(**over) -> Experiment:
    base = dict(
        id="adv_test",
        hypothesis="the change makes runs slower",
        why_chain=["one judge per comparison", "null rate from two runs"],
        switch="DELFIN_SOME_SWITCH",
    )
    base.update(over)
    return Experiment(**base)


def _registered(**over) -> Experiment:
    exp = _draft(**over)
    return pre_register(exp, expectation="runtime rises < 5%",
                        reading="wall-clock per run, one judge",
                        pool_size="small")


# ── A2.2 goalpost-moving after the first measurement must be detectable ──


def test_prereg_content_is_tamper_evident_after_measurement():
    """The recorded measurement must carry a digest of the pre-registered
    content, so a later mutation of expectation/reading is provable
    against the record instead of invisible.

    Reviewer s6 post-fix note: the fix implements the digest with its own
    construction (_prereg_digest: delimited ordered fields, 16-hex
    truncation).  The FINDING was the tamper-evidence PROPERTY, so this
    test pins the property -- record carries the experiment's digest,
    and the digest changes when pre-registered content is mutated --
    not one particular hash scheme.
    """
    exp = _registered()
    rec = record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})
    assert rec.prereg_digest == exp.prereg_digest
    assert rec.prereg_digest == _prereg_digest(exp)
    mutated = _registered()
    mutated.expectation = "MUTATED GOALPOST"
    assert _prereg_digest(mutated) != rec.prereg_digest
    # Moving the goalposts afterwards must be provable: the digest of the
    # mutated experiment no longer matches the recorded one.
    exp.expectation = "quietly redefined as: anything better than nothing"
    assert rec.prereg_digest != _prereg_digest(exp)


# ── A2.6 direct status writes must not bypass the pre-registration gate ──


def test_status_write_to_recording_bypasses_the_gate():
    """A caller (or the agent being protected against) that sets status
    directly must NOT be able to measure without pre-registration.

    Reviewer s6 post-fix note: the fix makes status a read-only property
    (experiment.py:125-128) -- a direct write raises AttributeError
    (Python's mechanism for a setter-less property), which IS the
    refusal.  The finding demanded 'cannot bypass', not a specific
    exception type, so this pins the write failing and the gate holding.
    """
    exp = _draft()
    with pytest.raises(AttributeError):
        exp.status = "recording"
    assert exp.status == "draft"          # nothing moved
    with pytest.raises(ExperimentError):  # the gate itself still holds
        record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})


def test_status_cannot_move_backwards():
    """The docstring promises the status 'only moves forward'
    (experiment.py module docstring, Experiment class); a backwards write
    must not succeed -- otherwise a 'measured' experiment can be rewound
    and re-registered with different expectations.

    Reviewer s6 post-fix note: _move now enforces monotonic order across
    _STATUSES (experiment.py:130-141).  A backwards _move raises
    ExperimentError; a direct assignment raises AttributeError because
    the property has no setter.  The finding demanded 'no backwards
    movement', so both refusals are pinned.
    """
    exp = _registered()
    record_measurement(exp, case="c1", outcome={})
    with pytest.raises(AttributeError):
        exp.status = "draft"              # direct write: no setter
    assert exp.status != "draft"          # still measured
    with pytest.raises(ExperimentError):  # the internal path refuses too
        exp._move("draft")


# ── A2.3 refusals must come through ExperimentError, not a crash ─────────


def test_pre_register_none_expectation_is_a_refusal_not_a_crash():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation=None, reading="r", pool_size="small")


def test_record_measurement_none_case_is_a_refusal_not_a_crash():
    exp = _registered()
    with pytest.raises(ExperimentError):
        record_measurement(exp, case=None, outcome={})


# ── Green side: refusals that DO hold get pinned too ─────────────────────


def test_whitespace_only_hypothesis_is_refused():
    with pytest.raises(ExperimentError):
        _registered(hypothesis="   ")


def test_reregister_after_measuring_is_refused():
    exp = _registered()
    record_measurement(exp, case="c1", outcome={})
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="different now",
                     reading="r", pool_size="small")


def test_refusal_message_names_pre_registration():
    exp = _draft()
    with pytest.raises(ExperimentError) as excinfo:
        record_measurement(exp, case="c1", outcome={})
    assert "pre-register" in str(excinfo.value)
