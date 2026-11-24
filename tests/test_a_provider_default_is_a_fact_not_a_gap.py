"""An unresolved field reads as an oversight; a stated default does not.

The run manifest records ``model.temperature`` and ``model.seed``, and a
survey of 26 autonomous-research systems counts released seeds and
traces as the field's largest gap (38% against 83% for code,
arXiv:2608.05179). Our manifest asked the engine for both; the engine has
neither, and ``temperature`` and ``seed`` appear NOWHERE in this
product. Nothing here has ever set them.

So both fields were ``unresolved`` in every manifest ever written, which
reads as "we forgot to record it" when the truth is "nobody asked for
one and the provider's default applied". Those are different facts, and
the second is worth recording: a run at an unstated default is NOT
reproducible, and saying so is more useful than a hole that looks like
an oversight.

WHAT IS ASSERTED, in order of what it protects:

  * **today's behaviour is unchanged.** With nothing configured the
    request carries no sampling parameters at all. A module that records
    what was asked must not change what is asked, or every existing run
    becomes a different experiment.
  * **the three cases are distinguishable**: provider_default,
    configured, refused_by_endpoint. The last reports the values as
    ASKED and never as applied.
  * **"reproducible" is not claimed lightly**: only a configured seed
    with temperature 0.
  * **a refusal does not end a run.** An endpoint that does not know
    ``seed`` answers 400 to the whole request, and a parameter added for
    the record must never cost an unattended run its turn.
"""
from __future__ import annotations

import pytest

from delfin.agent import run_manifest, sampling


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    """No environment leaking in, and the process flag cleared: a test
    that disabled sampling must not decide the next test's answer."""
    monkeypatch.delenv(sampling.ENV_TEMPERATURE, raising=False)
    monkeypatch.delenv(sampling.ENV_SEED, raising=False)
    sampling._reset_for_tests()
    yield
    sampling._reset_for_tests()


# --- today's behaviour must not change ---------------------------------

def test_nothing_configured_sends_nothing(monkeypatch):
    """The request carries no sampling parameters, exactly as before."""
    assert sampling.configured(settings={}) == {}


def test_nothing_configured_is_recorded_as_the_provider_default():
    rec = sampling.record({})
    assert rec["source"] == sampling.PROVIDER_DEFAULT
    assert rec["temperature"] is None and rec["seed"] is None
    assert rec["reproducible"] is False


# --- the three cases ---------------------------------------------------

def test_the_environment_configures_sampling_for_one_process(monkeypatch):
    """A benchmark arm fixes sampling without editing the user's
    settings file."""
    monkeypatch.setenv(sampling.ENV_TEMPERATURE, "0")
    monkeypatch.setenv(sampling.ENV_SEED, "11")
    params = sampling.configured(settings={})
    assert params == {"temperature": 0.0, "seed": 11}
    rec = sampling.record(params)
    assert rec["source"] == sampling.CONFIGURED
    assert rec["reproducible"] is True


def test_the_settings_file_configures_it_too():
    params = sampling.configured(
        settings={"agent": {"sampling_temperature": 0.7}})
    assert params == {"temperature": 0.7}
    assert sampling.record(params)["reproducible"] is False, (
        "a temperature without a seed is not reproducible")


def test_the_environment_wins_over_the_settings_file(monkeypatch):
    monkeypatch.setenv(sampling.ENV_TEMPERATURE, "0.2")
    params = sampling.configured(
        settings={"agent": {"sampling_temperature": 0.9}})
    assert params["temperature"] == pytest.approx(0.2)


def test_a_refusal_reports_what_was_asked_and_not_what_applied():
    rec = sampling.record({"temperature": 0.0, "seed": 5},
                          source=sampling.REFUSED)
    assert rec["source"] == sampling.REFUSED
    assert rec["temperature"] == 0.0 and rec["seed"] == 5
    assert rec["reproducible"] is False, (
        "a refused parameter never made the run reproducible")


# --- reproducible is a strong word -------------------------------------

@pytest.mark.parametrize("params", [
    {},                                   # nothing asked
    {"temperature": 0.0},                 # no seed
    {"seed": 3},                          # unstated temperature
    {"temperature": 0.7, "seed": 3},      # a seed does not fix sampling
])
def test_reproducible_is_false_unless_both_are_pinned(params):
    assert sampling.record(params)["reproducible"] is False


def test_reproducible_is_true_only_for_a_pinned_pair():
    assert sampling.record({"temperature": 0.0,
                            "seed": 3})["reproducible"] is True


# --- a bad value is not sent -------------------------------------------

@pytest.mark.parametrize("raw", ["-1", "3", "abc", ""])
def test_a_temperature_outside_the_accepted_range_is_not_sent(monkeypatch,
                                                              raw):
    """Sending it would be a 400 on every turn of the run."""
    monkeypatch.setenv(sampling.ENV_TEMPERATURE, raw)
    assert "temperature" not in sampling.configured(settings={})


# --- a refusal must not cost the rest of the run -----------------------

def test_a_refusal_stops_the_asking_for_the_process(monkeypatch):
    monkeypatch.setenv(sampling.ENV_SEED, "4")
    assert sampling.configured(settings={}) == {"seed": 4}
    sampling.disable()
    assert sampling.configured(settings={}) == {}, (
        "after a refusal the parameters must not be offered again")


def test_without_strips_only_the_sampling_parameters():
    payload = {"model": "m", "messages": [], "temperature": 0.0, "seed": 1,
               "stream": True}
    out = sampling.without(payload)
    assert "temperature" not in out and "seed" not in out
    assert out["model"] == "m" and out["stream"] is True


@pytest.mark.parametrize("text,expected", [
    ("Error code: 400 - unsupported parameter: 'seed'", True),
    ("Error code: 400 - invalid value for temperature", True),
    ("Error code: 400 - messages: too long", False),
    ("Error code: 400 - model not found", False),
    ("Error code: 500 - internal", False),
])
def test_only_a_refusal_naming_a_parameter_triggers_the_retry(text, expected):
    """Narrow on purpose: retrying a 400 about the messages without the
    sampling parameters would hide a real error behind a pointless
    second request."""
    assert sampling.looks_like_a_parameter_refusal(
        RuntimeError(text)) is expected


# --- what the manifest then says ---------------------------------------

def test_the_manifest_states_the_default_instead_of_leaving_a_gap(tmp_path,
                                                                 monkeypatch):
    monkeypatch.setenv("DELFIN_STATE", str(tmp_path / "state"))
    man = run_manifest.build(session_id="s1", workspace=tmp_path,
                             termination="done", role="solo_agent",
                             settings={"model": "m", "provider": "p",
                                       "sampling": sampling.record({})})
    assert man["model"]["sampling_source"] == sampling.PROVIDER_DEFAULT
    assert man["model"]["reproducible"] is False
    assert "model.temperature" not in man["unresolved"]
    assert "model.seed" not in man["unresolved"]


def test_the_manifest_carries_a_configured_pair(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_STATE", str(tmp_path / "state"))
    man = run_manifest.build(
        session_id="s1", workspace=tmp_path, termination="done",
        role="solo_agent",
        settings={"model": "m", "provider": "p",
                  "sampling": sampling.record({"temperature": 0.0,
                                               "seed": 7})})
    assert man["model"]["temperature"] == 0.0
    assert man["model"]["seed"] == 7
    assert man["model"]["sampling_source"] == sampling.CONFIGURED
    assert man["model"]["reproducible"] is True


def test_a_caller_that_passes_bare_values_still_works(tmp_path, monkeypatch):
    """Backwards compatible: the benchmark helpers pass temperature and
    seed directly."""
    monkeypatch.setenv("DELFIN_STATE", str(tmp_path / "state"))
    man = run_manifest.build(session_id="s1", workspace=tmp_path,
                             termination="done", role="solo_agent",
                             settings={"model": "m", "provider": "p",
                                       "temperature": 0.0, "seed": 2})
    assert man["model"]["seed"] == 2
    assert man["model"]["sampling_source"] == sampling.CONFIGURED
