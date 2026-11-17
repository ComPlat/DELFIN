"""Where the dashboard starts its sessions, and what it does with a service.

Default: in the dashboard's kernel, shown in the agent tab's full chat. A
running agent service is still joined for its list (sessions started in the
terminal or on the phone), but it is never started for the list alone.
``agent.dashboard_sessions: "service"`` (or DELFIN_DASHBOARD_SESSIONS) starts
dashboard sessions in the service, which is then started when none runs.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import agent_sessions as AS                # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402


@pytest.fixture
def no_setting(monkeypatch):
    monkeypatch.delenv(AS.SESSIONS_ENV, raising=False)
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})


def test_without_a_setting_sessions_start_in_the_kernel(no_setting):
    assert AS.sessions_where() == "kernel"


@pytest.mark.parametrize("raw, where", [("service", "service"), ("SERVICE ", "service"),
                                        ("kernel", "kernel"), ("anything", "kernel")])
def test_the_variable_chooses(monkeypatch, raw, where):
    monkeypatch.setenv(AS.SESSIONS_ENV, raw)
    assert AS.sessions_where() == where


def test_the_setting_chooses_when_no_variable(monkeypatch):
    monkeypatch.delenv(AS.SESSIONS_ENV, raising=False)
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings",
                        lambda *a, **k: {"agent": {"dashboard_sessions": "service"}})
    assert AS.sessions_where() == "service"


def _connect_calls(monkeypatch):
    calls = []

    def fake_connect(*, start=True, path=None):
        calls.append(start)
        raise LS.ServiceUnavailable("no service in this test")
    monkeypatch.setattr(LS, "connect", fake_connect)
    return calls


def test_kernel_sessions_join_a_running_service_but_never_start_one(no_setting, monkeypatch):
    calls = _connect_calls(monkeypatch)
    connector = AS._service_connector(None, explicit_build=False)
    assert connector is not None                  # the list of a running service
    with pytest.raises(LS.ServiceUnavailable):
        connector()
    assert calls == [False]


def test_service_sessions_start_the_service(monkeypatch):
    monkeypatch.setenv(AS.SESSIONS_ENV, "service")
    calls = _connect_calls(monkeypatch)
    with pytest.raises(LS.ServiceUnavailable):
        AS._service_connector(None, explicit_build=False)()
    assert calls == [True]


def test_a_stand_in_build_has_no_service(no_setting):
    assert AS._service_connector(None, explicit_build=True) is None
    assert AS._service_connector(False, explicit_build=False) is None
