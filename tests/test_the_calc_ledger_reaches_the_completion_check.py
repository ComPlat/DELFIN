"""The calc-ledger reaches the completion check, not just the schema.

Welle-10/s1 added the consumer: ``check_completion_claim(calcs=...)``
verdicts a calculation task on clean runs. It had no live caller --
``api_client._check_completion`` passed only ``observed`` + ``tests``,
so a calc task could never be verified (or honestly flagged) in a real
session; it always fell through to ``unchecked``.

Red on the previous commit: the wiring is missing, so every verdict
here comes back ``unchecked`` regardless of what the ledger carries.
The test goes through the PUBLIC path the same way the completion
tests do -- a real task in a real store, a client whose bound
``evidence_provider`` serves the ledgers, then
``client._check_completion(task_id, perms)``. Subject wording matches
what the ``_CALC_TASK_RE`` keying actually fires on ("berechne").
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

import delfin.agent.api_client as A
from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
)
from delfin.agent import agent_tasks

CALC_SUBJECT = "Berechne die Geometrie von Wasser"


@pytest.fixture
def calc_fixture(tmp_path):
    """A tiny calculation tree the indexer recognises: a run folder with
    an .inp (the file a calculation leaves behind) plus its .out, same
    shape as the corpus-names test."""
    calc = tmp_path / "my-calcs"
    (calc / "run_a").mkdir(parents=True)
    (calc / "run_a" / "run_a.inp").write_text(
        "! PBE0 def2-SVP\n* xyz 0 1\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n*\n",
        encoding="utf-8")
    (calc / "run_a" / "run_a.out").write_text(
        "ORCA\n! PBE0 def2-SVP\nFINAL SINGLE POINT ENERGY  -76.4\n",
        encoding="utf-8")
    return calc


@pytest.fixture
def agent_home(monkeypatch, tmp_path):
    """Redirect the change journal + task store into the tmp dir."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    agent_tasks._STORES.clear()
    return ws


def _client(agent_home):
    """The real document executor -- _check_completion lives on it and
    reads the task through its own task store."""
    return _doc_executor


def _perms(ws: Path, ledgers: dict) -> KitToolPermissions:
    """The permissions the completion check reads its ledgers from --
    the same attribute set_permissions binds, here hand-bound so the
    test pins the ledgers exactly."""
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id="sess-calc-wiring")
    perms.evidence_provider = lambda: dict(ledgers)
    return perms


def _create_and_start(ws: Path) -> int:
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id="sess-calc-wiring")
    out = json.loads(_doc_executor.execute(
        "task_create", {"subject": CALC_SUBJECT}, perms))
    tid = out["task"]["id"]
    _doc_executor.execute(
        "task_update", {"task_id": tid, "status": "in_progress"}, perms)
    return tid


def _ledgers(calcs):
    return {"observed": set(), "tests": [], "calcs": calcs}


def test_a_clean_run_verifies_a_calc_task(agent_home):
    tid = _create_and_start(agent_home)
    client = _client(agent_home)
    perms = _perms(agent_home, _ledgers([
        {"ts": time.time(), "folder": "calc/water",
         "outcome": "succeeded", "worst": "ok"}]))
    res = client._check_completion(tid, perms)
    assert res["verdict"] == "verified", res
    assert res["kind"] == "calc", res


def test_a_calc_task_without_a_ledger_is_unchecked_not_accused(agent_home):
    """Without a served "calcs" ledger the task falls through to the
    honest unknown -- a provider that predates the wiring must behave
    exactly as before it."""
    tid = _create_and_start(agent_home)
    client = _client(agent_home)
    perms = _perms(agent_home, {"observed": set(), "tests": []})
    res = client._check_completion(tid, perms)
    assert res["verdict"] == "unchecked", res


def test_a_failed_run_flags_the_calc_task(agent_home):
    tid = _create_and_start(agent_home)
    client = _client(agent_home)
    perms = _perms(agent_home, _ledgers([
        {"ts": time.time(), "folder": "calc/water",
         "outcome": "failed", "worst": ""}]))
    res = client._check_completion(tid, perms)
    assert res["verdict"] == "unmet", res
    assert res["kind"] == "calc_red", res


def test_a_run_outside_the_window_does_not_verify(agent_home):
    """The window is the task's start: a run recorded before the task
    began proves nothing about the work done since."""
    tid = _create_and_start(agent_home)
    client = _client(agent_home)
    perms = _perms(agent_home, _ledgers([
        {"ts": 0.0, "folder": "calc/water",
         "outcome": "succeeded", "worst": "ok"}]))
    res = client._check_completion(tid, perms)
    assert res["verdict"] == "unmet", res
    assert res["kind"] == "calc_none", res


# ---------------------------------------------------------------------------
# The producer: something actually FILLS ledgers["calcs"]
# ---------------------------------------------------------------------------

def _fresh_executor():
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    ex._calc_engine = None
    ex._calc_dirs = {}
    ex._calc_roots = {}
    ex._calc_evidence = []
    return ex


def test_get_calc_info_fills_the_calc_evidence_ledger(calc_fixture):
    """The wiring is only real when the ledger fills itself: a
    get_calc_info call through the executor records the run's outcome
    and the result-critic verdict into _calc_evidence with the exact
    keys the completion check's calc branch reads. Without this the
    provider serves an always-empty list and calcs is None in every
    real session -- the check never fires (operator gate, Welle 10)."""
    ex = _fresh_executor()
    ex._calc_dirs = {"calc": str(calc_fixture)}
    assert ex._ensure_calc_loaded()
    out = ex.execute("get_calc_info", {"calc_id": "run_a"},
                     A.KitToolPermissions(mode="default",
                                          workspace=str(calc_fixture)))
    info = json.loads(out)
    assert "error" not in info, info
    entries = ex._calc_evidence
    assert len(entries) == 1, entries
    e = entries[0]
    assert e["folder"].endswith("run_a"), e
    assert e["outcome"], e
    assert set(e) >= {"ts", "folder", "outcome", "worst"}, e


def test_get_calc_info_miss_still_leaves_the_ledger_untouched(calc_fixture):
    """A lookup that found nothing records nothing: the ledger says what
    the session observed, not what it tried."""
    ex = _fresh_executor()
    ex._calc_dirs = {"calc": str(calc_fixture)}
    assert ex._ensure_calc_loaded()
    out = ex.execute("get_calc_info", {"calc_id": "no_such_run"},
                     A.KitToolPermissions(mode="default",
                                          workspace=str(calc_fixture)))
    assert json.loads(out).get("error")
    assert ex._calc_evidence == []


def test_the_provider_serves_the_ledger_the_executor_fills():
    """End-to-end over the wiring seam: set_permissions binds a provider
    whose "calcs" key reads the executor's _calc_evidence, so a run the
    executor observed reaches the completion check without anyone
    copying it across. The executor is the module singleton -- exactly
    the instance get_calc_info runs on."""
    A._doc_executor._calc_evidence.append({
        "ts": 1500.0, "folder": "calc/water",
        "outcome": "succeeded", "worst": "ok"})
    try:
        perms = A.KitToolPermissions(mode="default", workspace="/tmp")
        A.OpenAIClient.set_permissions(
            A.OpenAIClient.__new__(A.OpenAIClient), perms)
        ledgers = perms.evidence_provider()
        assert ledgers["calcs"] == [{
            "ts": 1500.0, "folder": "calc/water",
            "outcome": "succeeded", "worst": "ok"}]
    finally:
        A._doc_executor._calc_evidence.clear()
