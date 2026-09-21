"""The doctor, the terminal and the dashboard say what protects this host.

The audit of 2026-09-16 found every degradation of the command isolation
silent outside the dashboard's panel, and a provider key exported in the
user's shell readable by every process of the user through /proc. On this
very host the settings had agent.bash_isolation = "off", which nothing
reported.
"""
import inspect

from delfin import doctor as D
from delfin.agent import process_guard as PG


def test_the_doctor_names_a_switched_off_isolation(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **k: {"agent": {"bash_isolation": "off"}})
    r = D.check_command_isolation()
    assert r.status == D.MISSING and "switched off" in r.detail
    assert "auto" in r.fix_hint


def test_the_doctor_reports_what_isolates_commands(monkeypatch):
    import delfin.agent.api_client as A
    import delfin.agent.socket_guard as SG
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **k: {"agent": {"bash_isolation": "auto"}})
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(SG, "available", lambda: True)
    assert D.check_command_isolation().status == D.OK
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    r = D.check_command_isolation()
    assert r.status == D.MISSING and "Landlock unavailable" in r.detail


def test_an_exported_key_is_named_never_shown(monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "value-that-must-not-appear-0123")
    r = D.check_no_exported_key()
    assert r.status == D.MISSING and "KIT_TOOLBOX_API_KEY" in r.detail
    assert "value-that-must-not-appear" not in r.detail + r.fix_hint
    assert "credentials set" in r.fix_hint
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY")
    for name in PG.PROVIDER_KEYS:
        monkeypatch.delenv(name, raising=False)
    assert D.check_no_exported_key().status == D.OK


def test_the_server_guard_is_a_loadable_server_extension():
    from delfin.dashboard import server_guard
    assert server_guard._jupyter_server_extension_points() == [
        {"module": "delfin.dashboard.server_guard"}]
    from delfin import cli_voila
    assert "'delfin.dashboard.server_guard': True" in inspect.getsource(cli_voila)


def test_terminal_launcher_and_kernel_warn_about_an_exported_key():
    from delfin import cli_voila
    from delfin.agent import cli
    from delfin.dashboard import tab_agent
    for src in (inspect.getsource(cli.main), inspect.getsource(cli_voila.main),
                inspect.getsource(tab_agent._register_process_exit_cleanup)):
        assert "exported_provider_keys()" in src
