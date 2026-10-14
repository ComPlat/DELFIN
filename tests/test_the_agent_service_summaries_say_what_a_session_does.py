"""What the agent view's rows say comes from the session's own events.

A summary carries, besides the session's state, the last message sent to
it, what it does right now, the question it waits on, and the first line
of its last answer -- in a few words, written by no model. A session can
be renamed; a turn ended by a cancel says so.
"""

from __future__ import annotations

import pytest

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceError
from tui_service import ServiceProcess, turn
from tui_worker_client import WorkerClient


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _summary(c, sid: str) -> dict:
    return next(s for s in c.request("sessions.list")["sessions"] if s["sessionId"] == sid)


def test_a_summary_says_what_was_asked_what_runs_and_what_came_back(service):
    a = service.client()
    sid = service.session(a, permissionMode="default")
    turn(a, sid, "hello there")
    s = _summary(a, sid)
    assert s["lastPrompt"] == "hello there"
    assert s["answer"] == "Summary"                          # the first line that says something
    assert s["activity"] is None and s["waiting"] is None
    mark = len(a.events)
    a.request("session.prompt", {"sessionId": sid, "text": "please edit the parser"})
    a.wait_for(lambda e: e.get("type") == "permission.asked", after=mark)
    a.wait_for(lambda e: e.get("type") == "session.summary" and e["summary"]["waiting"], after=mark)
    s = _summary(a, sid)
    assert s["waiting"] == "edit_file delfin/xtb.py"
    assert s["lastPrompt"] == "please edit the parser"
    # While it ran a tool, a summary said which.
    tools = [e["summary"]["activity"] for e in a.events[mark:] if e.get("type") == "session.summary"
             and (e["summary"].get("activity") or "").startswith("read_file")]
    assert tools and tools[0] == "read_file  delfin/xtb.py"


def test_a_session_can_be_renamed(service):
    a, b = service.client(), service.client()
    sid = service.session(a)
    mark = len(b.events)
    assert a.request("session.setTitle", {"sessionId": sid, "title": "  parser   fix "}) == {"title": "parser fix"}
    b.wait_for(lambda e: e.get("type") == "session.summary" and e["summary"].get("title") == "parser fix",
               after=mark)
    with pytest.raises(ServiceError) as empty:
        a.request("session.setTitle", {"sessionId": sid, "title": "   "})
    assert empty.value.code == protocol.INVALID_PARAMS


def test_a_turn_ended_by_a_cancel_says_it_stopped(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("text.delta", after=start)
        w.request("session.cancel")
        ended, _ = w.wait_for("turn.ended", after=start)
        assert ended["stopped"] is True
        start, end = w.turn("hello")
        assert "stopped" not in w.events[end]
