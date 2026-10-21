"""MCP servers and hooks start with the environment bash children get.

Input: a process environment holding harmless stand-ins for credentials
(``GH_TOKEN``, ``GITHUB_TOKEN``, ``AWS_SECRET_ACCESS_KEY``,
``AWS_SESSION_TOKEN``, a provider key, a ``*_TOKEN`` and a ``*_SECRET``
variable), and a stdio MCP server or a hook command that writes which of
them it can see into a file under ``tmp_path``. Output: those files.

Semantics: the child's environment is ``api_client._scrubbed_bash_env``
-- no credential the agent's own commands do not get -- plus what its
configuration names explicitly: an ``env`` list passes those variables
from the parent, an ``env`` object sets values, where ``${NAME}`` takes
NAME from the parent. Ordinary variables (PATH, HOME, a plain setting)
still arrive.
"""
from __future__ import annotations

import json
import sys
import time

import pytest

from delfin.agent import hooks as HK
from delfin.agent import mcp_client as MC

_SECRETS = {
    "GH_TOKEN": "fake-gh", "GITHUB_TOKEN": "fake-github",
    "AWS_SECRET_ACCESS_KEY": "fake-aws", "AWS_SESSION_TOKEN": "fake-aws-s",
    "KIT_TOOLBOX_API_KEY": "fake-kit", "SOME_SERVICE_TOKEN": "fake-svc",
    "LAB_DB_SECRET": "fake-db",
}
_NAMES = sorted(_SECRETS) + ["PLAIN_SETTING"]


@pytest.fixture
def secrets(monkeypatch):
    for k, v in _SECRETS.items():
        monkeypatch.setenv(k, v)
    monkeypatch.setenv("PLAIN_SETTING", "kept")


_PROBE = ("import json, os, sys; "
          "json.dump({k: os.environ.get(k) for k in sys.argv[2:]}, "
          "open(sys.argv[1], 'w'))")


def _spawn_mcp(tmp_path, env):
    out = tmp_path / "mcp-env.json"
    srv = MC._server_from_config("probe", {
        "command": sys.executable, "args": ["-c", _PROBE, str(out), *_NAMES],
        "env": env})
    try:
        srv.start()
        deadline = time.monotonic() + 20
        while not out.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        for _ in range(40):
            try:
                return json.loads(out.read_text())
            except (OSError, ValueError):
                time.sleep(0.05)
        raise AssertionError(srv.last_error or "probe wrote nothing")
    finally:
        srv.stop() if hasattr(srv, "stop") else None


def test_an_mcp_server_gets_no_credential_it_did_not_name(secrets, tmp_path):
    seen = _spawn_mcp(tmp_path, {})
    leaked = {k: v for k, v in seen.items() if k in _SECRETS and v}
    assert leaked == {}, leaked
    assert seen["PLAIN_SETTING"] == "kept"


def test_an_mcp_server_receives_what_its_env_list_names(secrets, tmp_path):
    seen = _spawn_mcp(tmp_path, ["GH_TOKEN"])
    assert seen["GH_TOKEN"] == "fake-gh"
    assert seen["GITHUB_TOKEN"] is None and seen["AWS_SECRET_ACCESS_KEY"] is None


def test_an_mcp_server_env_object_sets_and_references(secrets, tmp_path):
    seen = _spawn_mcp(tmp_path, {"LAB_DB_SECRET": "${LAB_DB_SECRET}",
                                 "SOME_SERVICE_TOKEN": "literal"})
    assert seen["LAB_DB_SECRET"] == "fake-db"
    assert seen["SOME_SERVICE_TOKEN"] == "literal"
    assert seen["GH_TOKEN"] is None


def _hook_env(tmp_path, entry_env=None):
    out = tmp_path / "hook-env.json"
    cmd = {"type": "command",
           "command": f"{sys.executable} -c \"{_PROBE}\" {out} " + " ".join(_NAMES)}
    if entry_env is not None:
        cmd["env"] = entry_env
    cfg = HK.HooksConfig()
    HK._merge_hooks(cfg, {"hooks": {"PreToolUse": [{"matcher": "", "hooks": [cmd]}]}},
                    source=str(tmp_path / "settings.json"))
    res = HK.run_hooks("PreToolUse", cfg, tool_name="bash", arguments={},
                       workspace=tmp_path)
    assert res and res[0].exit_code == 0, res
    return json.loads(out.read_text())


def test_a_hook_gets_no_credential_it_did_not_name(secrets, tmp_path):
    seen = _hook_env(tmp_path)
    leaked = {k: v for k, v in seen.items() if k in _SECRETS and v}
    assert leaked == {}, leaked
    assert seen["PLAIN_SETTING"] == "kept"


def test_a_hook_receives_what_its_env_list_names(secrets, tmp_path):
    seen = _hook_env(tmp_path, ["GITHUB_TOKEN"])
    assert seen["GITHUB_TOKEN"] == "fake-github"
    assert seen["GH_TOKEN"] is None


def test_the_hook_environment_still_carries_the_event(secrets):
    env = HK._build_env("PreToolUse", tool_name="bash")
    assert env["DELFIN_HOOK_EVENT"] == "PreToolUse"
    assert "GH_TOKEN" not in env and "AWS_SECRET_ACCESS_KEY" not in env
