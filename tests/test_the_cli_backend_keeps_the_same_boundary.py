"""The Claude CLI backend keeps the boundary every other backend keeps.

Its tools run inside the ``claude`` process, and under Bypass that
process was started with ``--dangerously-skip-permissions``: reads and
writes anywhere the account could, while the same chip on every other
backend asked before reading outside the working folders and refused
writing outside them. The user's rule (2026-10-09): "the CLI and the
dashboard must behave identically" -- and Bypass does not mean the agent
may walk into everything and read it.

``cli_gate`` is DELFIN's own gate as a PreToolUse hook. These tests drive
its decision function the way the hook does, and the relay that carries
its questions to the dialog.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from delfin.agent import cli_gate, kit_settings


@pytest.fixture
def scene(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "in.txt").write_text("inside\n")
    far = tmp_path / "far"
    far.mkdir()
    (far / "note.txt").write_text("OUTSIDE\n")

    def build(mode="bypassPermissions", can_ask=False):
        token = cli_gate.new_token()
        cli_gate.write_state(token, workspace=str(ws), mode=mode,
                             can_ask=can_ask)
        return token

    return build, ws, far


def _ev(tool, cwd, **tool_input):
    return {"tool_name": tool, "tool_input": tool_input, "cwd": str(cwd)}


def test_inside_the_workspace_nothing_is_asked(scene):
    build, ws, _far = scene
    token = build()
    assert cli_gate.decide(token, _ev("Read", ws, file_path=str(ws / "in.txt"))) is None
    assert cli_gate.decide(token, _ev("Bash", ws, command="cat in.txt")) is None
    assert cli_gate.decide(token, _ev("Write", ws, file_path="new.txt")) is None


@pytest.mark.parametrize("tool,args", [
    ("Read", {"file_path": "{far}/note.txt"}),
    ("Grep", {"path": "{far}", "pattern": "x"}),
    ("Glob", {"path": "{far}", "pattern": "*"}),
    ("Bash", {"command": "cat {far}/note.txt"}),
    ("Bash", {"command": "cd {far} && cat note.txt"}),
    ("Bash", {"command": "ls {far}"}),
])
def test_bypass_does_not_read_outside_unasked(scene, tool, args):
    build, ws, far = scene
    token = build(can_ask=False)
    args = {k: v.format(far=far) for k, v in args.items()}
    reason = cli_gate.decide(token, _ev(tool, ws, **args))
    assert reason, f"{tool} {args} read outside the workspace under Bypass"


@pytest.mark.parametrize("tool,args", [
    ("Write", {"file_path": "{far}/x.txt", "content": "x"}),
    ("Edit", {"file_path": "{far}/note.txt"}),
    ("Bash", {"command": "echo x > {far}/x.txt"}),
])
def test_bypass_does_not_write_outside(scene, tool, args):
    build, ws, far = scene
    token = build()
    args = {k: v.format(far=far) for k, v in args.items()}
    assert cli_gate.decide(token, _ev(tool, ws, **args))


class _Approver:
    def __init__(self, answer=True, expire=False):
        self.asked = []
        self.answer = answer
        self.expire = expire
        self.last_timed_out = False
        self.last_refusal_reason = ""

    def callback(self, tool, args, preview):
        self.asked.append((tool, args.get("path")))
        self.last_timed_out = self.expire
        return self.answer and not self.expire


def _decide_with_relay(token, event, approver):
    relay = cli_gate.Relay(token, approver.callback, poll_s=0.05).start()
    try:
        return cli_gate.decide(token, event)
    finally:
        relay.stop()


def test_outside_read_is_asked_through_the_same_dialog(scene, monkeypatch):
    build, ws, far = scene
    token = build(can_ask=True)
    approver = _Approver(answer=True)
    ev = _ev("Read", ws, file_path=str(far / "note.txt"))
    assert _decide_with_relay(token, ev, approver) is None
    assert approver.asked and approver.asked[0][0] == "read_file"
    # The grant opened the directory for reads: the next call, a new hook
    # process, does not ask again ...
    again = _Approver(answer=False)
    assert _decide_with_relay(token, ev, again) is None
    assert again.asked == []
    # ... and still does not let anything be written there.
    assert cli_gate.decide(token, _ev("Write", ws, file_path=str(far / "y")))


def test_a_refusal_holds_for_the_next_call(scene):
    build, ws, far = scene
    token = build(can_ask=True)
    ev = _ev("Read", ws, file_path=str(far / "note.txt"))
    assert _decide_with_relay(token, ev, _Approver(answer=False))
    later = _Approver(answer=True)
    assert _decide_with_relay(
        token, _ev("Bash", ws, command=f"cat {far}/note.txt"), later)
    assert later.asked == []


def test_nobody_answering_is_not_a_refusal(scene, monkeypatch):
    build, ws, far = scene
    monkeypatch.setattr(cli_gate, "HOOK_WAIT_S", 1.0)
    token = build(can_ask=True)
    ev = _ev("Read", ws, file_path=str(far / "note.txt"))
    reason = _decide_with_relay(token, ev, _Approver(expire=True))
    assert reason and "TIMED OUT" in reason
    state = json.loads(cli_gate._state_path(token).read_text())
    assert str((far / "note.txt").resolve()) not in state["denied_paths"]


def test_a_saved_read_dir_is_readable_never_writable(scene):
    build, ws, far = scene
    kit_settings.persist_read_dir(far)
    token = build(can_ask=False)
    assert cli_gate.decide(
        token, _ev("Read", ws, file_path=str(far / "note.txt"))) is None
    assert cli_gate.decide(
        token, _ev("Write", ws, file_path=str(far / "z.txt")))


def test_the_hook_prints_a_deny_decision(scene, monkeypatch, capsys):
    build, ws, far = scene
    token = build()
    monkeypatch.setattr("sys.stdin", io.StringIO(json.dumps(
        _ev("Read", ws, file_path=str(far / "note.txt")))))
    assert cli_gate.main(["cli_gate", token]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_a_missing_state_denies(scene, monkeypatch, capsys):
    _build, ws, _far = scene
    monkeypatch.setattr("sys.stdin", io.StringIO(json.dumps(
        _ev("Read", ws, file_path=str(ws / "in.txt")))))
    cli_gate.main(["cli_gate", "nosuchtoken"])
    out = json.loads(capsys.readouterr().out)
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_the_cli_process_is_started_with_the_hook(scene, monkeypatch):
    import sys

    from delfin.agent import api_client as ac
    _build, ws, _far = scene
    seen = {}

    class _Proc:
        def poll(self):
            return None

    def _popen(cmd, **kw):
        seen["cmd"] = cmd
        return _Proc()

    monkeypatch.setattr(ac.subprocess, "Popen", _popen)
    client = ac.CLIClient(claude_path=sys.executable, cwd=str(ws),
                          permission_mode="bypassPermissions")
    client._ensure_proc("system")
    settings = json.loads(seen["cmd"][seen["cmd"].index("--settings") + 1])
    hook = settings["hooks"]["PreToolUse"][0]["hooks"][0]["command"]
    assert "delfin.agent.cli_gate" in hook
    assert client._gate_token in hook
    assert cli_gate._state_path(client._gate_token).is_file()


def test_a_repository_cannot_add_read_dirs(scene, tmp_path):
    _build, ws, far = scene
    (ws / ".delfin").mkdir()
    (ws / ".delfin" / "settings.json").write_text(json.dumps(
        {"kit": {"read_dirs": [str(far)]}}))
    assert kit_settings.load(repo_dir=ws).read_dirs == []


def test_home_cannot_be_saved_as_a_read_dir(scene):
    with pytest.raises(ValueError):
        kit_settings.persist_read_dir(Path.home())
