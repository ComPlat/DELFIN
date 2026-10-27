"""`delfin-agent tui --continue` joins the workspace's session when it runs.

Input: a private agent service with the scripted engine. Measured on a
Pi: a window was closed, its session kept running in the service, and
``tui --continue`` in the same folder started a second worker for that
session, which refused at once ("is being written by another process").
The latest session of the folder was the one already running.

Output: ``session.create`` with ``continue`` attaches to that session --
same id, no second worker -- as a resume of a sleeping session does.
"""
from __future__ import annotations

import time

import pytest

from tui_service import ServiceProcess, turn


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def test_continue_attaches_to_the_running_session_of_the_folder(service):
    first = service.client()
    sid = service.session(first, fakeSave=True)
    turn(first, sid, "hello")
    first.close()
    window = service.client()
    # The worker saves after the turn's end event: wait until it is listed.
    deadline = time.monotonic() + 30
    while not window.request("sessions.past", {"workspace": str(service.workspace())})["sessions"]:
        assert time.monotonic() < deadline, "the session was never saved"
        time.sleep(0.05)
    joined = window.request("session.create", {
        "cwd": str(service.workspace()), "continue": True, "fake": True,
        "fakeDelay": 0.003}, timeout=60)
    assert joined["sessionId"] == sid
    listed = window.request("sessions.list")["sessions"]
    assert [s["sessionId"] for s in listed] == [sid]
