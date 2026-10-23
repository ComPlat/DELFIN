"""A CI watch, and a cluster job a background shell submitted, belong to the
session that started them.

Review 2026-09-16: watches carry a session since 81cc5cb1 and
check_agent_jobs(session_id=...) reads only that session's -- but
register_ci_watch stored no session and the SLURM ids found in a background
shell's output were registered without one. A session that pushed was never
told about its CI; the Background panel, which hides unowned work, did not
show the watch either. And _waiting_on_watched_jobs read every watch in the
workspace, so another session's cluster job held this session's
auto-continue.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as bj
from delfin.agent import job_monitor as jm

_SHA = "269429eb"


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    return tmp_path


def _entries(ws):
    return jm.load_watched(ws / ".delfin" / "agent_watched_jobs.json")["jobs"]


def test_a_ci_watch_carries_the_session_that_pushed(ws):
    jid = jm.register_ci_watch(ws, "ComPlat/DELFIN", _SHA, "main", session_id="s-1")
    assert _entries(ws)[jid]["session_id"] == "s-1"
    assert _entries(ws)[jid]["branch"] == "main"


def test_a_ci_watch_without_a_session_stays_unowned(ws):
    jid = jm.register_ci_watch(ws, "ComPlat/DELFIN", _SHA, "main")
    assert "session_id" not in _entries(ws)[jid]


def test_the_push_hands_its_session_to_the_watch():
    import inspect
    src = inspect.getsource(A._after_push)
    assert 'session_id=str(getattr(perms, "task_session_id", "") or "")' in src


def test_a_submitted_cluster_job_carries_the_shells_session(ws, monkeypatch):
    monkeypatch.setattr(bj, "_submitted_slurm_ids", lambda path: ["4976064"])
    got = bj._watch_submitted_jobs(ws, {"job_id": "bg-1", "stdout_path": "x", "session_id": "s-2"})
    assert got == ["4976064"]
    assert _entries(ws)["4976064"]["session_id"] == "s-2"


def _waiting(ws, sid):
    client = SimpleNamespace(_permissions=SimpleNamespace(workspace=ws, task_session_id=sid))
    return A.OpenAIClient._waiting_on_watched_jobs(client)


def test_another_sessions_watch_does_not_hold_this_turn(ws):
    jm.register_agent_job(ws, "4976064", "opt freq", extra={"session_id": "other"})
    assert _waiting(ws, "mine") is False


def test_an_unowned_watch_does_not_hold_a_session_either(ws):
    jm.register_agent_job(ws, "4976064", "opt freq")
    assert _waiting(ws, "mine") is False


def test_this_sessions_watch_still_holds_it(ws):
    jm.register_ci_watch(ws, "ComPlat/DELFIN", _SHA, "main", session_id="mine")
    assert _waiting(ws, "mine") is True


def test_a_client_without_a_session_waits_for_everything(ws):
    """The CLI and the daemon have no session; nothing changes for them."""
    jm.register_agent_job(ws, "4976064", "opt freq", extra={"session_id": "other"})
    assert _waiting(ws, "") is True
