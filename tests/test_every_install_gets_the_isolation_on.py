"""Every install gets the shell isolation on, once, and keeps a later choice.

A settings file with agent.bash_isolation = "off" ran locked and
unattended sessions with no isolation at all -- no filesystem limits, no
session-socket guard, no network proxy. On the host this wave was built
on, the file said "off" and nothing reported it. The update switches it to
"auto" once, says so, and leaves a later "off" alone.
"""
import json

from delfin import user_settings as US


def _write(path, data):
    path.write_text(json.dumps(data))


def test_off_becomes_auto_once_with_a_notice(tmp_path):
    path = tmp_path / "settings.json"
    _write(path, {"agent": {"bash_isolation": "off"}})
    settings = US.load_settings(path)
    assert settings["agent"]["bash_isolation"] == "auto"
    on_disk = json.loads(path.read_text())
    assert on_disk["agent"]["bash_isolation"] == "auto"
    assert US._ISOLATION_ON in on_disk["security_updates_applied"]
    notices = US.take_security_notices(path)
    assert len(notices) == 1 and "set it back" in notices[0]
    assert US.take_security_notices(path) == []            # shown once


def test_a_later_off_is_kept(tmp_path):
    path = tmp_path / "settings.json"
    _write(path, {"agent": {"bash_isolation": "off"}})
    US.load_settings(path)
    data = json.loads(path.read_text())
    data["agent"]["bash_isolation"] = "off"                # the user's choice now
    _write(path, data)
    assert US.load_settings(path)["agent"]["bash_isolation"] == "off"


def test_auto_is_recorded_without_a_notice(tmp_path):
    path = tmp_path / "settings.json"
    _write(path, {"agent": {"bash_isolation": "auto"}})
    US.load_settings(path)
    assert US.take_security_notices(path) == []
    assert US._ISOLATION_ON in json.loads(path.read_text())["security_updates_applied"]


def test_an_existing_file_gets_the_network_defaults(tmp_path):
    path = tmp_path / "settings.json"
    _write(path, {"agent": {"model": "kit.glm-5.3"}})
    net = US.load_settings(path)["agent"]["sandbox_network"]
    assert net["mode"] == "proxy" and net["allowed_ports"] == [80, 443]
    assert json.loads(path.read_text())["agent"]["sandbox_network"]["mode"] == "proxy"


def test_the_front_ends_show_the_notice():
    import inspect
    from delfin import cli_voila
    from delfin.agent import cli
    from delfin.dashboard import tab_agent
    assert "take_security_notices" in inspect.getsource(cli._show_security_notices)
    assert "take_security_notices" in inspect.getsource(cli_voila.main)
    assert "take_security_notices" in inspect.getsource(tab_agent._register_process_exit_cleanup)
