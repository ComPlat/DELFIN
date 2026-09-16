"""A session's schedule tools reach its own entries only.

cron_list returned every entry, with the prompts other sessions had
scheduled, and cron_delete removed any entry by id, so one session could
silently cancel another's wake-up (review 2026-09-16). cron_create also
recorded the process's working directory instead of the session's
workspace, the fault schedule_wakeup had already been fixed for.
"""
import json
from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A
from delfin.agent import scheduler as S


@pytest.fixture
def sch(tmp_path, monkeypatch):
    s = S.Scheduler(path=tmp_path / "cron.json")
    monkeypatch.setattr(S, "get_scheduler", lambda path=None: s)
    return s


def _call(name, args, sid, ws):
    perms = SimpleNamespace(workspace=ws, task_session_id=sid)
    return json.loads(A._doc_executor._execute_scheduler(name, args, perms))


def test_a_session_sees_only_its_own_entries(sch, tmp_path):
    a = _call("cron_create", {"every_seconds": 600, "prompt": "A's check"}, "sess-a", tmp_path)
    _call("cron_create", {"every_seconds": 600, "prompt": "B's private plan"}, "sess-b", tmp_path)
    listed = _call("cron_list", {}, "sess-a", tmp_path)["entries"]
    assert [e["id"] for e in listed] == [a["id"]]
    assert "B's private plan" not in json.dumps(listed)


def test_a_session_cannot_delete_another_sessions_entry(sch, tmp_path):
    b = _call("cron_create", {"every_seconds": 600, "prompt": "B"}, "sess-b", tmp_path)
    out = _call("cron_delete", {"entry_id": b["id"]}, "sess-a", tmp_path)
    assert out["status"] == "not_found"
    assert any(e.id == b["id"] for e in sch.list_entries())
    assert _call("cron_delete", {"entry_id": b["id"]}, "sess-b", tmp_path)["status"] == "ok"


def test_a_caller_without_a_session_owns_the_unowned_entries(sch, tmp_path):
    mine = _call("cron_create", {"every_seconds": 600, "prompt": "cli"}, "", tmp_path)
    _call("cron_create", {"every_seconds": 600, "prompt": "dashboard"}, "sess-a", tmp_path)
    assert [e["id"] for e in _call("cron_list", {}, "", tmp_path)["entries"]] == [mine["id"]]


def test_a_recurring_entry_records_the_sessions_workspace(sch, tmp_path):
    ent = _call("cron_create", {"every_seconds": 600, "prompt": "x"}, "sess-a", tmp_path)
    stored = next(e for e in sch.list_entries() if e.id == ent["id"])
    assert stored.workspace == str(tmp_path)
