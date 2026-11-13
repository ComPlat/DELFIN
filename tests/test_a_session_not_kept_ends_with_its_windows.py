"""Keep session, in the agent service: one switch per session.

keep on (the default, as every service session ran before): the session
runs on after its windows close -- asleep after an idle hour, woken when
used. keep off: it ends (saved, resumable) once every window that had it
open with keep off has closed its connection. A window holds such a
session from creating it with ``keep: false``, switching it off
(``session.setKeep``) or attaching to it while it is off; a detach does
not let go (a dashboard detaches the sessions it hides), the window's
connection closing does. ``/keep [on|off]`` typed in a window is the same
switch. Every window's summary says ``keep``.
"""
from __future__ import annotations

import time

import pytest

from delfin.agent.service import protocol, supervisor
from tui_service import ServiceProcess, turn

WAIT_S = supervisor.UNUSED_GRACE_S + 30.0


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _summary(c, sid: str) -> dict:
    rows = c.request("sessions.list")["sessions"]
    return next((s for s in rows if s["sessionId"] == sid), {})


def _gone(c, sid: str, timeout: float = WAIT_S) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not _summary(c, sid):
            return True
        time.sleep(0.3)
    return False


def test_the_protocol_names_the_switch():
    assert protocol.SERVICE_REQUESTS["session.setKeep"] == {
        "sessionId": "string", "keep": "boolean"}
    assert protocol.SERVICE_REQUESTS["session.create"]["keep"] == "boolean?"
    assert "session.setKeep" in supervisor.ACTING


def test_a_session_is_kept_by_default_and_says_so(service):
    window = service.client()
    sid = service.session(window)
    turn(window, sid, "hello")
    assert _summary(window, sid)["keep"] is True
    window.close()
    watcher = service.client()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert _summary(watcher, sid), "a kept session ran on"


def test_switched_off_it_ends_when_its_window_closes(service):
    window = service.client()
    sid = service.session(window)
    turn(window, sid, "hello")
    assert window.request("session.setKeep", {"sessionId": sid, "keep": False}) == {
        "keep": False}
    watcher = service.client()
    assert _summary(watcher, sid)["keep"] is False
    window.close()
    assert _gone(watcher, sid), "a session switched off ran on after its window closed"
    status = watcher.request("service.status")
    assert sid not in {s.get("sessionId") for s in status.get("dormant") or []}, \
        "nor is it kept asleep"


def test_created_off_it_ends_with_its_window(service):
    window = service.client()
    sid = service.session(window, keep=False)
    turn(window, sid, "hello")
    watcher = service.client()
    window.close()
    assert _gone(watcher, sid)


def test_a_window_reopening_a_kept_session_does_not_unkeep_it(service):
    """keep false on session.create is for the session the window starts.
    One that already runs here -- the window resumes or continues it -- is
    attached as it is: a kept session stays kept."""
    a = service.client()
    sid = service.session(a, fakeSave=True)
    turn(a, sid, "hello")
    assert _summary(a, sid)["keep"] is True
    b = service.client()
    # --continue: the folder's latest session, which runs here already.
    got = b.request("session.create", {"cwd": str(service.workspace()), "continue": True,
                                       "fake": True, "keep": False}, timeout=60)
    assert got["sessionId"] == sid
    a.close()
    b.close()
    watcher = service.client()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert _summary(watcher, sid).get("keep") is True, "a kept session was ended"


def test_a_detach_does_not_end_it(service):
    """A dashboard detaches the sessions it hides; that is no closing."""
    window = service.client()
    sid = service.session(window, keep=False)
    turn(window, sid, "hello")
    window.request("session.detach", {"sessionId": sid})
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert _summary(window, sid), "a hidden view's session ended"


def test_every_window_that_opened_it_holds_it(service):
    a = service.client()
    sid = service.session(a, keep=False)
    turn(a, sid, "hello")
    b = service.client()
    b.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=60)
    a.close()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert _summary(b, sid), "ended while another window still had it open"
    watcher = service.client()
    b.close()
    assert _gone(watcher, sid)


def test_switched_on_again_it_runs_on(service):
    window = service.client()
    sid = service.session(window, keep=False)
    turn(window, sid, "hello")
    window.request("session.setKeep", {"sessionId": sid, "keep": True})
    window.close()
    watcher = service.client()
    time.sleep(supervisor.UNUSED_GRACE_S + 2.0)
    assert _summary(watcher, sid)["keep"] is True


def test_keep_typed_in_the_window_is_the_same_switch(service):
    window = service.client()
    sid = service.session(window)
    turn(window, sid, "hello")
    window.request("session.prompt", {"sessionId": sid, "text": "/keep off"})
    watcher = service.client()
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline and _summary(watcher, sid).get("keep") is not False:
        time.sleep(0.2)
    assert _summary(watcher, sid)["keep"] is False
    window.close()
    assert _gone(watcher, sid)

