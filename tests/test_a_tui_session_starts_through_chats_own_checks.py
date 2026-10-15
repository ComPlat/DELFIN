"""A session the agent service starts is set up by the same function as
`delfin-agent chat`: cli._start_interactive_session. The administrator's
managed policy therefore refuses it exactly where it refuses chat -- a
mode named for the session that the policy disallows, a required
isolation layer the host lacks -- and with chat's text.

The worker used to repeat chat's set-up step by step and had left out the
managed-policy start check, the policy's mode coercion and the isolation
start notice (isolation_posture.apply_at_start).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin.agent import cli
from delfin.agent import managed_policy as MP
from delfin.agent.service import worker


@pytest.fixture
def policy(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    path = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(path),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)
    monkeypatch.chdir(tmp_path)
    MP.clear_cache()

    def write(body: dict) -> None:
        path.write_text(json.dumps(body))
        path.chmod(0o644)
        MP.clear_cache()

    ws = tmp_path / "ws"
    ws.mkdir()
    yield write, ws
    MP.clear_cache()


def _engine_must_not_be_built(_ns):
    raise AssertionError("the engine was built for a start the policy refuses")


def test_the_worker_refuses_a_mode_the_managed_policy_disallows(policy, monkeypatch):
    write, ws = policy
    write({"disallowed_modes": ["acceptEdits"]})
    monkeypatch.setattr(cli, "_build_engine", _engine_must_not_be_built)
    expected = MP.launch_decision("acceptEdits", explicit=True, attended=True)
    assert expected[2], expected
    args = worker.parse_args(["--cwd", str(ws), "--permission-mode", "acceptEdits"])
    with pytest.raises(worker.SessionRefused) as refused:
        worker.open_real_session(args, ws)
    assert expected[1] in str(refused.value)


def test_the_worker_refuses_when_required_isolation_is_missing(policy, monkeypatch):
    write, ws = policy
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": False)
    write({"require_isolation": True})
    monkeypatch.setattr(cli, "_build_engine", _engine_must_not_be_built)
    args = worker.parse_args(["--cwd", str(ws)])
    with pytest.raises(worker.SessionRefused) as refused:
        worker.open_real_session(args, ws)
    assert "an interactive session does not start" in str(refused.value)


def test_a_mode_from_settings_falls_back_as_in_chat(policy, monkeypatch):
    """Not named for the session (a settings file's default): the policy
    lowers it with a note instead of refusing, through coerce_mode."""
    write, ws = policy
    write({"disallowed_modes": ["acceptEdits"]})

    class _Perms:
        mode = "default"

    class _Engine:
        session_id = "s1"
        kit_permissions = _Perms()

        def set_kit_permission_mode(self, mode, announce=True):
            self.kit_permissions.mode = mode
            return True

    built: list = []
    monkeypatch.setattr(cli, "_build_engine", lambda ns: built.append(ns) or _Engine())
    monkeypatch.setattr(cli, "_open_session", lambda engine, ns, w: True)
    monkeypatch.setattr(cli, "_persisted_default_mode",
                        lambda w: ("acceptEdits", str(ws / "settings.json")))
    monkeypatch.setattr(cli, "_startup_banner", lambda *a, **k: "BANNER")
    monkeypatch.setattr(cli, "_bounding_notices", lambda ns, engine: [])
    args = worker.parse_args(["--cwd", str(ws)])
    engine, _ns, banner, _notices = worker.open_real_session(args, ws)
    assert built, "the engine was not built"
    assert engine.kit_permissions.mode != "acceptEdits"


def test_chat_and_the_worker_call_one_set_up(policy, monkeypatch):
    _write, ws = policy
    calls: list = []

    class _Stop(Exception):
        pass

    def start(ns, report, workspace, *, can_ask):
        calls.append((Path(workspace), can_ask))
        raise _Stop()

    monkeypatch.setattr(cli, "_start_interactive_session", start)
    args = worker.parse_args(["--cwd", str(ws)])
    with pytest.raises(_Stop):
        worker.open_real_session(args, ws)
    assert calls == [(ws.resolve(), False)]


def test_the_worker_ends_a_session_through_chats_session_end(monkeypatch, tmp_path):
    ended: list = []
    monkeypatch.setattr(cli, "_finish_interactive_session",
                        lambda engine, agent, workspace, *, title, t0, learn=True:
                        ended.append((title, learn)))
    worker.end_real_session(object(), tmp_path, "t", agent=None, t0=0.0, learn=False)
    assert ended == [("t", False)]
