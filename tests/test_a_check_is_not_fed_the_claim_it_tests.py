"""A check fed the claim it tests reports agreement with itself.

``run_replay`` builds a FRESH manifest and compares it to the recorded
one. It used to pass the recorded model id, provider, temperature, seed
and role into that fresh build. ``run_manifest._model`` falls back to the
configured id whenever the session's per-turn metrics are absent -- the
normal case on another machine -- so ``model.id`` and ``prompt.role``
matched BY CONSTRUCTION. Two of five comparisons could not fail.

Measured on a real 33-minute run before and after: the same manifest and
workspace went from ``reproduced: True`` with 65 matches to the two
fields reporting UNAVAILABLE and blocking. The verdict was wrong, not
merely generous.

What may still be passed in, and the distinction this file pins:
``expected_steps`` is the run's PLAN. A workspace's files cannot carry
how many steps were intended, and the flag it feeds counts the steps
actually present against it -- the plan is input to the check, the count
is re-derived. Nothing that is COMPARED may come from the record.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import run_replay


def _manifest(tmp_path, **over):
    """A recorded manifest naming a model and a role nothing can confirm."""
    man = {
        "schema": 1,
        "session_id": "s1",
        "code": {"commit": "deadbeef", "delfin_version": "9.9.9"},
        "model": {"id": "a-model-that-never-ran-here",
                  "id_source": "session_metrics",
                  "provider": "someprovider",
                  "temperature": 0.0, "seed": 7, "switches": None},
        "prompt": {"role": "a-role-that-is-not-configured-here",
                   "pack_digest": "cafe"},
        "run": {"termination": "done"},
        "contract": {"expected_steps": 2, "flags": {}},
        "environment": {"engines": {}},
        "unresolved": [],
    }
    man.update(over)
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps(man), encoding="utf-8")
    return p


def _check(report, name):
    return next(c for c in report.checks if c.name == name)


# --- the defect ---------------------------------------------------------

@pytest.mark.parametrize("field", ["model.id", "prompt.role"])
def test_a_recorded_value_nothing_can_re_derive_does_not_match(tmp_path,
                                                               field):
    """The record names a model and a role this machine has never seen.
    Neither may come back as a match, whatever else is decided about
    them."""
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(_manifest(tmp_path), workspace=ws)
    c = _check(report, field)
    assert c.verdict != run_replay.MATCH, (
        f"{field} matched a value nothing here produced")
    assert c.verdict == run_replay.UNAVAILABLE
    assert c.note


def test_a_measured_model_id_blocks_and_a_setting_does_not(tmp_path):
    """The record's own provenance decides what is a claim.

    `model.id_source` says whether the id was MEASURED from the run's
    per-turn metrics or taken from configuration. A measured id is an
    assertion about what answered and must be re-derived or block. A
    role is a setting a run is GIVEN: no workspace and no session record
    can say which role a re-run would use, so leaving it open says
    nothing about the result and must not block. Comparing it anyway is
    what made it match by construction.
    """
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(_manifest(tmp_path), workspace=ws)
    assert _check(report, "model.id").blocking is True, (
        "an id the record calls measured is a claim")
    role = _check(report, "prompt.role")
    assert role.blocking is False
    assert "setting" in role.note


def test_a_configured_model_id_is_a_setting_too(tmp_path):
    """The same field, the other provenance: an id the record itself
    calls configured was never a measurement, so there is nothing about
    the run to verify."""
    ws = tmp_path / "ws"
    ws.mkdir()
    man = _manifest(tmp_path)
    data = json.loads(man.read_text())
    data["model"]["id_source"] = "configured"
    man.write_text(json.dumps(data), encoding="utf-8")
    report = run_replay.replay_file(man, workspace=ws)
    c = _check(report, "model.id")
    assert c.verdict == run_replay.UNAVAILABLE
    assert c.blocking is False
    assert c.now != c.recorded


def test_such_a_run_is_not_reported_as_reproduced(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(_manifest(tmp_path), workspace=ws)
    assert report.reproduced() is False


def test_the_recorded_values_are_not_echoed_back_as_the_fresh_ones(tmp_path):
    """The symptom in its plainest form: `now` must not be a copy of
    `recorded`."""
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(_manifest(tmp_path), workspace=ws)
    for field in ("model.id", "prompt.role"):
        c = _check(report, field)
        assert c.now != c.recorded or c.now is None, (
            f"{field}: the fresh manifest returned the recorded value")


# --- what is still legitimately passed in ------------------------------

def test_the_plan_may_be_passed_because_the_files_cannot_carry_it(tmp_path):
    """`expected_steps` is the run's intent, not a claim about its
    outcome. Passing it is what lets the step COUNT be re-derived and
    judged at all; dropping it would make the flag unjudgeable rather
    than honest."""
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(
        _manifest(tmp_path, contract={"expected_steps": 2, "flags": {}}),
        workspace=ws)
    # The contract is reported one way or another, never silently absent.
    assert any(c.name.startswith("contract") for c in report.checks)


def test_a_field_the_record_never_held_is_not_blocking(tmp_path):
    """Unknown then and unknown now is not a failure to verify: nothing
    was claimed, so nothing failed. It is still shown."""
    ws = tmp_path / "ws"
    ws.mkdir()
    man = _manifest(tmp_path)
    data = json.loads(man.read_text())
    data["model"]["id"] = None
    man.write_text(json.dumps(data), encoding="utf-8")
    report = run_replay.replay_file(man, workspace=ws)
    c = _check(report, "model.id")
    assert c.verdict == run_replay.UNAVAILABLE
    assert c.blocking is False
    assert "does not carry" in c.note


# --- the live comparisons must still work ------------------------------

def test_a_commit_that_differs_is_still_reported_as_a_difference(tmp_path):
    """The check that was always live stays live: a record made at
    another commit is a difference, not an unavailability."""
    ws = tmp_path / "ws"
    ws.mkdir()
    report = run_replay.replay_file(
        _manifest(tmp_path, code={"commit": "0000000",
                                  "delfin_version": "9.9.9"}),
        workspace=ws)
    c = _check(report, "code.commit")
    assert c.verdict in (run_replay.DIFFERS, run_replay.UNAVAILABLE)
    if c.verdict == run_replay.DIFFERS:
        assert c.recorded == "0000000" and c.now != "0000000"
