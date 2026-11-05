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

