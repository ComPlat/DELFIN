"""The external agent CLIs run their own tools, so DELFIN's gates never see
them. What DELFIN hands them decides what they may do (review 2026-09-16):

- DELFIN's "auto" started the Claude CLI with --dangerously-skip-permissions,
  every check off, although the CLI has an auto mode of its own.
- Codex ran "acceptEdits" and "auto" with danger-full-access: the whole disk,
  no sandbox.
- Neither got the secret deny list, and both inherited every key in the
  environment, the KIT and the other provider's included.
"""
import json
import subprocess

import pytest

from delfin.agent import api_client as A

KIT = "kit-secret-value-0123456789"


class _Captured(Exception):
    pass


@pytest.fixture
def launch(monkeypatch):
    seen = {}

    def fake_popen(cmd, **kw):
        seen["cmd"], seen["env"] = cmd, kw.get("env")
        raise _Captured()

    monkeypatch.setattr(subprocess, "Popen", fake_popen)
    # The host is supplied, not measured: the client looks for the claude
    # binary, which a CI runner does not have.
    real_which = A.shutil.which
    monkeypatch.setattr(A.shutil, "which", lambda name, *a, **k:
                        "/usr/bin/claude" if name == "claude" else real_which(name, *a, **k))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", KIT)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-own-key-for-the-cli")
    monkeypatch.setenv("OPENAI_API_KEY", "sk-openai-own-key-for-codex")
    return seen


def _claude(mode, launch):
    c = A.CLIClient.__new__(A.CLIClient)
    A.CLIClient.__init__(c, model="sonnet", claude_path="claude", permission_mode=mode)
    with pytest.raises(_Captured):
        c._ensure_proc("system")
    return launch["cmd"], launch["env"]


def test_auto_is_the_clis_own_mode_not_every_check_off(launch):
    cmd, _ = _claude("auto", launch)
    assert "--dangerously-skip-permissions" not in cmd
    assert cmd[cmd.index("--permission-mode") + 1] == "auto"
    cmd, _ = _claude("bypassPermissions", launch)
    assert "--dangerously-skip-permissions" in cmd


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "auto", "bypassPermissions"])
def test_the_claude_cli_always_gets_the_secret_deny_rules(launch, mode):
    cmd, env = _claude(mode, launch)
    rules = json.loads(cmd[cmd.index("--settings") + 1])["permissions"]["deny"]
    for rule in ("Read(**/.env)", "Read(~/.ssh/**)", "Edit(**/*.pem)",
                 "Read(~/.delfin/credentials.json)"):
        assert rule in rules
    assert KIT not in json.dumps(env) and "OPENAI_API_KEY" not in env
    assert env["ANTHROPIC_API_KEY"] == "sk-ant-own-key-for-the-cli"


def _codex(mode, launch, monkeypatch):
    monkeypatch.setattr(A.shutil, "which", lambda x: "/usr/bin/codex")
    c = A.CodexCLIClient(model="gpt-5.4", codex_path="codex", permission_mode=mode)
    with pytest.raises(_Captured):
        list(c.stream_message(system="s", messages=[{"role": "user", "content": "hi"}]))
    return launch["cmd"], launch["env"]


@pytest.mark.parametrize("mode", ["acceptEdits", "auto"])
def test_codex_accepting_edits_stays_in_the_workspace(launch, monkeypatch, mode):
    cmd, env = _codex(mode, launch, monkeypatch)
    assert cmd[cmd.index("--sandbox") + 1] == "workspace-write"
    assert "danger-full-access" not in cmd
    assert KIT not in json.dumps(env) and "ANTHROPIC_API_KEY" not in env
    assert env["OPENAI_API_KEY"] == "sk-openai-own-key-for-codex"


def test_only_bypass_gives_codex_the_whole_disk(launch, monkeypatch):
    cmd, _ = _codex("bypassPermissions", launch, monkeypatch)
    assert cmd[cmd.index("--sandbox") + 1] == "danger-full-access"
