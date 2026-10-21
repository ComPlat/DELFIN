"""An allow rule the user wrote answers before the caged route.

Input: ``permissions.repos`` in a redirected user settings file with an
allow rule for ``./run.sh``; ``./run.sh`` run through the real executor
(``_DocToolExecutor.execute``) in acceptEdits, with the caged route
(``_caged_command_may_run_unasked``) replaced by a recorder that would
admit everything.

Output, pinned here:

* with the rule, the caged route is never consulted, the command runs,
  nobody is asked, and the rule's own event (``permission_rule_allow``)
  is recorded -- the command runs as the rule allows it, not in the
  caged route's narrower read cage;
* without the rule, the caged route is consulted as before;
* a rule for another command does not keep the caged route away.

Why: since #303 the caged route came first. Its Landlock read list holds
only the workspace and the toolchain, so a program a rule names outside
them (the wave test gate under ``.delfin/gate-tools``) failed with EACCES,
exit 126, and the rule never answered (wave 16, 2026-10-16). The existing
allow-rule tests switch the caged route off, which hid it.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

from conftest import child_env
from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP
from delfin.agent import permission_rules as PR
from delfin.agent import security_events
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


class Asker:
    def __init__(self, answer: bool):
        self.answer = answer
        self.asked: list = []

    def __call__(self, name, args, preview):
        self.asked.append(name)
        return self.answer


def _git(tmp_path, where, *args):
    subprocess.run(["git", *args], cwd=str(where), check=True,
                   capture_output=True, env=child_env(tmp_path))


def _repo(tmp_path, name, remote=None):
    ws = tmp_path / name
    ws.mkdir()
    _git(tmp_path, ws, "init", "-q")
    if remote:
        _git(tmp_path, ws, "remote", "add", "origin", remote)
    run = ws / "run.sh"
    run.write_text("#!/bin/sh\necho RUN_MARKER\n")
    run.chmod(0o755)
    return ws


@pytest.fixture
def env(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    settings = home / US.SETTINGS_FILE_NAME
    monkeypatch.setattr(
        US, "get_settings_path",
        lambda base_path=None: (Path(base_path).expanduser() if base_path
                                else settings))
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    policy = etc / "managed_policy.json"
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(policy),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE",
                        frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)
    consulted = []

    def _caged(cmd, *a, **k):
        consulted.append(cmd)
        return True
    monkeypatch.setattr(A, "_caged_command_may_run_unasked", _caged)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    a = _repo(tmp_path, "a", "git@GitHub.com:me/a.git")
    b = _repo(tmp_path, "b", "https://github.com/me/b")
    local = _repo(tmp_path, "local")
    MP.clear_cache()
    PR._CACHE.clear()
    security_events.clear()
    yield {"settings": settings, "policy": policy, "a": a, "b": b,
           "local": local, "tmp": tmp_path, "consulted": consulted}
    MP.clear_cache()
    PR._CACHE.clear()


def _user(env, **block):
    env["settings"].write_text(json.dumps({"permissions": block}))
    PR._CACHE.clear()


def _user(env, **block):
    env["settings"].write_text(json.dumps({"permissions": block}))
    PR._CACHE.clear()


def _run(ws):
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits")
    perms.confirm_callback = Asker(False)
    out = _DocToolExecutor().execute("bash", {"command": "./run.sh"}, perms)
    return "RUN_MARKER" in out, bool(perms.confirm_callback.asked)


def _kinds():
    return [getattr(e, "kind", None) for e in security_events.recent(200)]


def test_a_rule_for_the_command_answers_before_the_caged_route(env):
    top = str(env["local"].resolve())
    _user(env, repos={top: {"allow": ["bash(./run.sh)"]}})
    ran, asked = _run(env["local"])
    assert env["consulted"] == []
    assert ran and not asked
    kinds = _kinds()
    assert "permission_rule_allow" in kinds
    assert "caged_unasked" not in kinds


def test_without_a_rule_the_caged_route_is_consulted(env):
    _user(env)
    _run(env["local"])
    assert env["consulted"] == ["./run.sh"]


def test_a_rule_for_another_command_keeps_the_caged_route(env):
    top = str(env["local"].resolve())
    _user(env, repos={top: {"allow": ["bash(./other.sh)"]}})
    _run(env["local"])
    assert env["consulted"] == ["./run.sh"]
