"""KG4 adversarial, phase 3: the accept scripts against hostile inputs.

Pre-registered angles (.gate/HANDOFF.md), applied to
delfin/agent/pack/benchmark/accept/w17_longrun_*_verify.py:
  * None never passes: each accept REJECTS when its longrun input is
    unmeasured (no inputs, steer_followed=None, compactions None, restart
    attestation None) -- the "measures nothing" failure mode of wave 15a.
  * The outcome artifact is checked against the FIXTURE's own arithmetic,
    not a constant the setup writes alongside the answer (a lazy setup
    must not make a wrong answer pass).
  * The steer task rejects the old goal's artifact surviving the steer.
  * Overwriting the fixture's steering marker cannot fake an accept.
All fixtures are built under tmp_path at run time; nothing in .gate/,
nothing in the user's memory or knowledge paths.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.pack.benchmark.accept import (
    w17_longrun_compaction_verify as compaction_verify,
    w17_longrun_restart_verify as restart_verify,
    w17_longrun_steer_verify as steer_verify,
    w17_longrun_wait_verify as wait_verify,
)


# --- wait_verify -----------------------------------------------------------

def _wait_workspace(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_wait_for_job"
    d.mkdir(parents=True)
    (d / "result.json").write_text(json.dumps(
        {"n": 6, "mean_kj_per_mol": 42.19}), encoding="utf-8")
    return tmp_path


def test_wait_accept_requires_measured_waiting(tmp_path):
    ws = _wait_workspace(tmp_path)
    ok, msg = wait_verify.check(ws, {"idle_with_mail_s": 12.0})
    assert ok, msg


def test_wait_rejects_unmeasured_waiting(tmp_path):
    ws = _wait_workspace(tmp_path)
    ok, msg = wait_verify.check(ws, None)
    assert not ok
    ok, msg = wait_verify.check(ws, {})
    assert not ok
    ok, msg = wait_verify.check(ws, {"idle_with_mail_s": 0.0})
    assert not ok


def test_wait_rejects_wrong_mean(tmp_path):
    ws = _wait_workspace(tmp_path)
    d = ws / "w17_longrun" / "lr_wait_for_job"
    (d / "result.json").write_text(json.dumps(
        {"n": 6, "mean_kj_per_mol": 40.0}), encoding="utf-8")
    ok, msg = wait_verify.check(ws, {"idle_with_mail_s": 12.0})
    assert not ok and "expected 42.19" in msg


# --- steer_verify ----------------------------------------------------------

def _steer_workspace(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_steered_goal"
    d.mkdir(parents=True)
    (d / "STEERING_MESSAGE.txt").write_text(
        "PERSON (mid-task steering, new goal): catalog answers only.",
        encoding="utf-8")
    (d / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent\n",
        encoding="utf-8")
    return tmp_path


def test_steer_accepts_new_goal_with_records(tmp_path):
    ws = _steer_workspace(tmp_path)
    ok, msg = steer_verify.check(ws, steer_followed=True)
    assert ok, msg


def test_steer_rejects_unattested_steer_followed(tmp_path):
    ws = _steer_workspace(tmp_path)
    for value in (None, False):
        ok, msg = steer_verify.check(ws, steer_followed=value)
        assert not ok
        assert "not attested" in msg


def test_steer_rejects_old_goal_surviving(tmp_path):
    ws = _steer_workspace(tmp_path)
    d = ws / "w17_longrun" / "lr_steered_goal"
    (d / "summary.txt").write_text("old goal answer\n", encoding="utf-8")
    ok, msg = steer_verify.check(ws, steer_followed=True)
    assert not ok and "OLD goal" in msg


def test_steer_rejects_missing_marker(tmp_path):
    ws = _steer_workspace(tmp_path)
    d = ws / "w17_longrun" / "lr_steered_goal"
    (d / "STEERING_MESSAGE.txt").write_text(
        "no marker here\n", encoding="utf-8")
    ok, msg = steer_verify.check(ws, steer_followed=True)
    assert not ok and "steering message is not in place" in msg


# --- compaction / restart --------------------------------------------------

def test_compaction_and_restart_reject_unmeasured(tmp_path):
    ws = tmp_path
    ok, msg = compaction_verify.check(ws)
    assert not ok, "compaction accept must reject with no fixture at all"
    ok, msg = restart_verify.check(ws)
    assert not ok, "restart accept must reject with no fixture at all"
