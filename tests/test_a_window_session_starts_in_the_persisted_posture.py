"""A session starts in the posture the flag or the user's settings name,
else in default -- in the window (the agent service's worker) as in chat,
which share cli._start_interactive_session.

Two starts in plan, both measured live in a window before this test:

- ``~/.delfin_settings.json`` ``agent.permission_mode`` (the user's own
  settings file, the one the dashboard writes) was never read: a kit file's
  ``default_mode`` -- written once by a permissions tool call, say -- decided
  the start posture alone, and the window opened in plan although the user's
  setting says default.
- the scripted (``--fake``) engine began in plan, so a test window never
  showed the posture a real one starts in.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import cli
from delfin.agent import kit_settings
from delfin.agent.service import worker
from delfin.agent.service.fake_engine import FakeEngine


@pytest.fixture
def home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    # The set-up changes into the workspace; undone after the test.
    monkeypatch.chdir(tmp_path)
    # Evaluated at import: the test's own home, not the machine's.
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH", home / ".delfin" / "settings.json")
    ws = tmp_path / "ws"
    ws.mkdir()
    return home, ws


def _user_settings(home: Path, mode: str) -> None:
    (home / ".delfin_settings.json").write_text(json.dumps({"agent": {"permission_mode": mode}}))


def _kit(home: Path, mode: str) -> None:
    (home / ".delfin" / "settings.json").write_text(json.dumps({"kit": {"default_mode": mode}}))


def test_the_users_permission_mode_setting_decides_the_start(home):
    home, ws = home
    _user_settings(home, "default")
    _kit(home, "plan")
    mode, source = cli._persisted_default_mode(ws)
    assert mode == "default"
    assert source.endswith(".delfin_settings.json")


def test_the_users_setting_is_read_when_no_kit_file_names_a_mode(home):
    home, ws = home
    _user_settings(home, "acceptEdits")
    assert cli._persisted_default_mode(ws)[0] == "acceptEdits"


@pytest.mark.parametrize("value", ["", "erlaubt", "ask_all", 3, None])
def test_a_value_that_is_no_mode_is_not_a_setting(home, value):
    home, ws = home
    (home / ".delfin_settings.json").write_text(json.dumps({"agent": {"permission_mode": value}}))
    _kit(home, "acceptEdits")
    assert cli._persisted_default_mode(ws)[0] == "acceptEdits"


def test_a_persisted_bypass_in_the_users_settings_is_still_not_inherited(home):
    home, ws = home
    _user_settings(home, "bypassPermissions")
    from delfin.agent import launch_guard
    persisted, source = cli._persisted_default_mode(ws)
    mode, _why = launch_guard.resolve_posture(persisted_mode=persisted, settings_path=source)
    assert mode == "default"


def test_a_window_session_with_default_settings_starts_in_default(home, monkeypatch):
    """The worker's real start (open_real_session -> chat's set-up), with
    the user's settings saying default and a kit file saying plan."""
    home, ws = home
    _user_settings(home, "default")
    _kit(home, "plan")

    class _Perms:
        mode = "plan"                    # what an engine may be built with

    class _Engine:
        session_id = "s1"
        kit_permissions = _Perms()

        def set_kit_permission_mode(self, mode, announce=True):
            self.kit_permissions.mode = mode
            return True

    monkeypatch.setattr(cli, "_build_engine", lambda ns: _Engine())
    monkeypatch.setattr(cli, "_open_session", lambda engine, ns, w: True)
    monkeypatch.setattr(cli, "_startup_banner", lambda *a, **k: "BANNER")
    monkeypatch.setattr(cli, "_bounding_notices", lambda ns, engine: [])
    args = worker.parse_args(["--cwd", str(ws)])
    engine, _ns, _banner, _notices = worker.open_real_session(args, ws)
    assert engine.kit_permissions.mode == "default"


def test_the_scripted_engine_starts_in_default():
    assert FakeEngine().kit_permissions.mode == "default"
