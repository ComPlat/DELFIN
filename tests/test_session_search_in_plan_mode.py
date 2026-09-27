"""session_search is read-only research: plan mode may run it.

The plan-mode safe list is enumerated (a refusal list is wrong the moment
a tool is added). session_search reads the session index and changes
nothing, so it belongs beside history_search. The contrast case pins the
boundary: a writing tool (skill_propose_patch) stays refused in plan mode.
"""

from __future__ import annotations

import json


from delfin.agent import api_client as A
from delfin.agent import session_store


def _plan(tmp_path):
    return A.KitToolPermissions(workspace=tmp_path, mode="plan")


def _ex():
    return A._DocToolExecutor.__new__(A._DocToolExecutor)


def test_session_search_runs_in_plan_mode(tmp_path, monkeypatch):
    from delfin.agent import session_index

    monkeypatch.setattr(session_store, "_SESSIONS_DIR",
                        tmp_path / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: tmp_path / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: tmp_path / "session_index.sqlite")
    d = tmp_path / "agent_sessions"
    d.mkdir(parents=True)
    (d / "sess-plan.json").write_text(json.dumps({
        "session_id": "sess-plan", "title": "plan mode research",
        "created_at": 1760000000,
        "chat_messages": [
            {"role": "user", "content": "the iridium pincer complex debate"}],
    }), encoding="utf-8")
    session_index.index_session("sess-plan")

    out = _ex().execute("session_search", {"query": "iridium"},
                        _plan(tmp_path))
    assert '"sess-plan"' in out
    # Not an error payload (my fixture title contains the phrase, so
    # check the error key rather than the raw string).
    assert '"error"' not in out


def test_a_writing_tool_stays_refused_in_plan_mode(tmp_path):
    out = _ex().execute("skill_propose_patch",
                        {"name": "x", "old": "a", "new": "b",
                         "reason": "r", "evidence": [
                             {"kind": "test", "ref": "t"}]},
                        _plan(tmp_path))
    assert "plan mode" in out.lower() or "rejected" in out.lower()
