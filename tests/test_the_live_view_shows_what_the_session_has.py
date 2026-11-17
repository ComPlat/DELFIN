"""The dashboard's view of a service session shows what the session has.

Rule: the dashboard must not lie. Each control of the live view
(live_tab) shows the session's own value, from its session.info /
session.updated, never only the last click. Driven against a real private
service with the scripted engine; each check has a control (the model,
which was always reported) so a failure says which field is missing.
"""

from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent.service.client import ServiceClient          # noqa: E402
from delfin.dashboard import live_sessions as LS                # noqa: E402
from tui_service import ServiceProcess                          # noqa: E402


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


@pytest.fixture
def tui(service):
    c = ServiceClient(service.socket, client="delfin-tui")
    service.clients.append(c)
    return c


@pytest.fixture
def dash(service):
    c = LS.LiveClient(service.socket)
    yield c
    c.close()


def test_a_session_started_with_an_effort_shows_that_effort(service, tui, dash):
    from delfin.dashboard import live_tab

    sid = service.session(tui, effort="high", model="fake-large")
    attached = tui.request("session.attach", {"sessionId": sid})
    # Control: the model it started with is reported.
    assert attached["session"]["model"] == "fake-large"
    assert attached["session"].get("effort") == "high"
    _widget, refs = live_tab.build_live_tab(dash, sid)
    assert _until(lambda: refs["model"].value == "fake-large")
    assert _until(lambda: refs["effort"].value == "high"), refs["effort"].value
    refs["close"]()



def test_a_mode_switch_reaches_every_window_and_the_picker_offers_no_dashboard(service, tui, dash):
    from delfin.dashboard import live_tab

    sid = service.session(tui)
    tui.request("session.attach", {"sessionId": sid})
    _widget, refs = live_tab.build_live_tab(dash, sid)
    assert _until(lambda: refs["mode"].value == "solo")
    # A service session cannot drive a dashboard's tabs: not offered.
    assert _until(lambda: "research" in [v for _, v in refs["mode"].options])
    offered = [v for _, v in refs["mode"].options]
    assert "research" in offered and "dashboard" not in offered
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "/mode research"})
    # Control: the command ran (its notice), then the posture-like update.
    tui.wait_for(lambda e: e.get("type") == "notice" and "research" in str(e.get("text")), after=mark)
    tui.wait_for(lambda e: e.get("type") == "session.updated" and e.get("mode") == "research",
                 after=mark, timeout=10)
    assert _until(lambda: refs["mode"].value == "research")
    refs["close"]()
