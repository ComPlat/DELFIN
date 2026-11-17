"""The agent service keeps no session nobody asked anything (supervisor.py
copy pair, .gate/px/new_supervisor.py).

A session started new whose last window lets it go -- a detach, or the
window's connection closing (a hangup, a kill) -- ends after
UNUSED_GRACE_S when it was never asked anything: no idle worker per window
opened and closed. A session that was asked something, a resumed one and
one about an attention item (a brief) run on as before; a prompt sent
right behind the detach still reaches it.
"""
from __future__ import annotations

import time

import pytest

from delfin.agent.service import supervisor
from tui_service import ServiceProcess, turn

WAIT_S = supervisor.UNUSED_GRACE_S + 20.0


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _listed(c) -> set[str]:
    return {s["sessionId"] for s in c.request("sessions.list")["sessions"]}


def _gone(c, sid: str, timeout: float = WAIT_S) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if sid not in _listed(c):
            return True
        time.sleep(0.2)
    return False


def test_an_unused_session_ends_when_its_window_closes(service):
    window = service.client()
    sid = service.session(window)
    window.close()
    watcher = service.client()
    assert _gone(watcher, sid), "the unused session kept running after its window closed"
    hello = watcher.request("service.status")
    assert sid not in {s["sessionId"] for s in hello.get("dormant") or []}, "nor is it kept asleep"


def test_an_unused_session_ends_when_its_window_detaches(service):
    window = service.client()
    sid = service.session(window)
    window.request("session.detach", {"sessionId": sid})
    assert _gone(window, sid)


def test_a_session_that_was_asked_runs_on(service):
    window = service.client()
    sid = service.session(window)
    turn(window, sid, "hello")
    window.close()
    watcher = service.client()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert sid in _listed(watcher)


def test_a_prompt_right_behind_the_detach_still_reaches_it(service):
    window = service.client()
    sid = service.session(window)
    window.request("session.detach", {"sessionId": sid})
    window.request("session.prompt", {"sessionId": sid, "text": "just in time"})
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert sid in _listed(window)


def test_a_session_about_an_item_runs_on_unasked(service):
    window = service.client()
    sid = service.session(window, brief="job 4711 failed: OUT_OF_MEMORY", briefTitle="job 4711")
    window.close()
    watcher = service.client()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert sid in _listed(watcher)


def test_a_session_another_window_still_shows_runs_on(service):
    a = service.client()
    sid = service.session(a)
    b = service.client()
    b.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=60)
    a.close()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert sid in _listed(b)
