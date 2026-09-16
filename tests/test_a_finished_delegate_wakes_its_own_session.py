"""A finished background delegate wakes the session that started it, not every
session in the dashboard.

Review 2026-09-16: the wake path told markers apart by the process stamp
only, and two sessions of one dashboard are one process -- both woke on the
same finished delegate, each with its own set of already-woken ids.
"""
import inspect
import json
import pathlib

from delfin.agent import background_view as BV
from delfin.agent import subagents as SA


def _markers(tmp_path, monkeypatch, owners):
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path)
    monkeypatch.setattr(SA, "_entry_owned_by_us", lambda entry: True)
    monkeypatch.setattr(SA, "read_running", lambda **_: {})
    for sa_id, owner in owners.items():
        rec = {"sa_id": sa_id, "type": "explore", "description": f"run {sa_id}"}
        if owner:
            rec["owner_session"] = owner
        (tmp_path / f"{sa_id}.json").write_text(json.dumps(rec))


def test_each_session_is_woken_by_its_own_delegate(tmp_path, monkeypatch):
    _markers(tmp_path, monkeypatch, {"a1": "s1", "b2": "s2"})
    assert [w["job_id"] for w in BV.finished_background_agents(set(), session_id="s1")] == ["a1"]
    assert [w["job_id"] for w in BV.finished_background_agents(set(), session_id="s2")] == ["b2"]


def test_an_unowned_marker_wakes_no_session(tmp_path, monkeypatch):
    _markers(tmp_path, monkeypatch, {"a1": ""})
    assert BV.finished_background_agents(set(), session_id="s1") == []


def test_without_a_session_the_old_rule_holds(tmp_path, monkeypatch):
    _markers(tmp_path, monkeypatch, {"a1": "s1", "b2": ""})
    assert sorted(w["job_id"] for w in BV.finished_background_agents(set())) == ["a1", "b2"]


def test_the_marker_carries_the_reserving_session(tmp_path, monkeypatch):
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path)
    monkeypatch.setattr(SA, "reap_pending_reports", lambda *a, **k: None)
    monkeypatch.setattr(SA, "read_running", lambda **_: {})
    SA._note_pending_report("x9", subagent_type="explore", description="d", owner_session="s7")
    assert json.loads((tmp_path / "x9.json").read_text())["owner_session"] == "s7"


def test_the_reservation_hands_its_owner_to_the_marker():
    src = inspect.getsource(SA.reserve_running)
    assert "owner_session=owner_session" in src


def test_the_dashboard_wakes_with_its_own_session():
    src = pathlib.Path(inspect.getfile(__import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text()
    i = src.index("finished_background_agents(")
    assert "session_id=_background_owner()" in src[i:i + 200]
