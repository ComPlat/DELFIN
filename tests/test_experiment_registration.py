"""Package G, Phase 2: the ``Experiment`` object and pre-registration.

An experiment the agent cannot fool itself with must refuse to *measure*
before its hypothesis, why-chain, switch, expectation and reading are
written down.  A measurement taken without pre-registration is exactly the
self-bias this package exists to stop, so the recorder is the gate: it
*raises* on a draft experiment instead of quietly recording.

The candidate change sits behind ONE switch, default off (additive by
construction); ``pool_size`` is small or large.
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
    base = dict(
        id="m2_test",
        hypothesis="raising StrictReporting makes the run slower",
        why_chain=["one judge per comparison",
                    "measuring twice with the same judge gives the null rate"],
        switch="DELFIN_STRICT_REPORTING",
    )
    base.update(over)
    return Experiment(**base)


# ── Switch behind one knob, default off ──────────────────────────────────


def test_candidate_change_sits_behind_one_switch_default_off():
    exp = _draft()
    assert exp.switch == "DELFIN_STRICT_REPORTING"
    assert status_of(exp) == "draft"


# ── Pre-registration refuses an incomplete experiment ────────────────────


def test_pre_register_refuses_empty_hypothesis():
    exp = _draft(hypothesis="")
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="runtime rises < 5%",
                     reading="wall-clock per run, one judge", pool_size="small")


def test_pre_register_refuses_empty_why_chain():
    exp = _draft(why_chain=[])
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="runtime rises < 5%",
                     reading="wall-clock per run, one judge", pool_size="small")


def test_pre_register_refuses_missing_expectation():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="",
                     reading="wall-clock per run, one judge", pool_size="small")


def test_pre_register_refuses_missing_reading():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="runtime rises < 5%",
                     reading="", pool_size="small")


def test_pre_register_refuses_unknown_pool_size():
    exp = _draft()
    with pytest.raises(ExperimentError):
        pre_register(exp, expectation="runtime rises < 5%",
                     reading="wall-clock per run", pool_size="medium")


@pytest.mark.parametrize("pool", ["small", "large"])
def test_pre_register_accepts_small_and_large(pool):
    exp = _draft()
    pre_register(exp, expectation="runtime rises < 5%",
                 reading="wall-clock per run, one judge", pool_size=pool)
    assert exp.pool_size == pool
    assert status_of(exp) == "registered"


# ── A measurement without pre-registration is refused ────────────────────


def test_measure_without_pre_registration_is_refused():
    exp = _draft()
    with pytest.raises(ExperimentError):
        record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})


def test_measure_after_pre_registration_records():
    exp = _draft()
    pre_register(exp, expectation="runtime rises < 5%",
                 reading="wall-clock per run, one judge", pool_size="small")
    record_measurement(exp, case="c1", outcome={"runtime_s": 12.0})
    assert len(exp.measurements) == 1
    assert status_of(exp) == "recording"


def test_pre_registered_fields_are_stored():
    exp = _draft()
    pre_register(exp, expectation="runtime rises < 5%",
                 reading="wall-clock per run, one judge", pool_size="large")
    assert exp.expectation == "runtime rises < 5%"
    assert exp.reading == "wall-clock per run, one judge"
