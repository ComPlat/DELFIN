"""What DELFIN's gate lets through, the Claude CLI runs.

The gate is a PreToolUse hook in front of the CLI's tools. It asked the
user -- and then said nothing, so the CLI asked again. Under ``-p`` the CLI
has nobody to ask: `delfin-agent tui` on its default backend refused every
write the user had just approved with y (measured on a Raspberry Pi 5,
claude 2.1.309). The hook now answers "allow" for the reads, writes and
commands it decided; the CLI's own deny rules still hold over an allow.

Commands: the CLI's shell runs outside DELFIN's cage, so in default and
acceptEdits a command off the auto-allow list is asked about, as on every
other backend, instead of being left to a CLI that cannot ask.
"""

from __future__ import annotations

import io
import json

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

    def build(mode="default", can_ask=True):
        token = cli_gate.new_token()
        cli_gate.write_state(token, workspace=str(ws), mode=mode, can_ask=can_ask)
        return token

    return build, ws


def _ev(tool, cwd, **tool_input):
    return {"tool_name": tool, "tool_input": tool_input, "cwd": str(cwd)}


@pytest.fixture
def answers(monkeypatch):
    """The user's answers, and what each question showed."""
    from delfin.agent import api_client as ac
    asked: list = []
    replies: list = []

    def ask(perms, name, args, preview):
        asked.append((name, preview))
        return replies.pop(0) if replies else False

    monkeypatch.setattr(ac, "_ask_confirm", ask)
    return asked, replies


def _hook(monkeypatch, capsys, token, event):
    monkeypatch.setattr("sys.stdin", io.StringIO(json.dumps(event)))
    assert cli_gate.main(["cli_gate", token]) == 0
    out = capsys.readouterr().out.strip()
    return json.loads(out)["hookSpecificOutput"] if out else None


def test_a_read_the_gate_passed_is_allowed(scene, monkeypatch, capsys):
    build, ws = scene
    decision = _hook(monkeypatch, capsys, build(),
                     _ev("Read", ws, file_path=str(ws / "in.txt")))
    assert decision["permissionDecision"] == "allow"


def test_an_approved_command_is_allowed_and_a_denied_one_denied(
        scene, answers, monkeypatch, capsys):
    build, ws = scene
    asked, replies = answers
    token = build()
    replies.append(True)
    decision = _hook(monkeypatch, capsys, token, _ev("Bash", ws, command="make build"))
    assert decision["permissionDecision"] == "allow"
    assert asked and asked[-1] == ("bash", "$ make build")

    replies.append(False)
    decision = _hook(monkeypatch, capsys, token, _ev("Bash", ws, command="make build"))
    assert decision["permissionDecision"] == "deny"
    assert "denied" in decision["permissionDecisionReason"]


def test_an_auto_allowed_command_is_not_asked(scene, answers, monkeypatch, capsys):
    build, ws = scene
    asked, _ = answers
    decision = _hook(monkeypatch, capsys, build(), _ev("Bash", ws, command="ls"))
    assert decision["permissionDecision"] == "allow"
    assert asked == []


def test_with_nobody_to_ask_a_command_is_refused_with_its_reason(scene, answers):
    build, ws = scene
    reason = cli_gate.decide(build(can_ask=False), _ev("Bash", ws, command="make build"))
    assert reason and "nobody can be asked" in reason


@pytest.mark.parametrize("mode", ["plan", "bypassPermissions"])
def test_plan_and_bypass_do_not_ask_for_a_command(scene, answers, mode):
    build, ws = scene
    asked, _ = answers
    assert cli_gate.decide(build(mode=mode), _ev("Bash", ws, command="make build")) is None
    assert asked == []


def test_a_tool_the_gate_does_not_decide_is_left_to_the_cli(scene, monkeypatch, capsys):
    build, ws = scene
    assert _hook(monkeypatch, capsys, build(), _ev("TodoWrite", ws, todos=[])) is None


def test_the_write_question_shows_the_content(scene, monkeypatch):
    """The dialog previewed the path alone: a new file was "(no changes)"."""
    from delfin.agent import api_client as ac
    build, ws = scene
    seen: list = []

    def gate(self, path, perms, name, args):
        seen.append((name, dict(args)))
        return None

    monkeypatch.setattr(type(ac._doc_executor), "_gate_write_path", gate)
    token = build()
    cli_gate.decide(token, _ev("Write", ws, file_path="new.txt", content="hi"))
    cli_gate.decide(token, _ev("Edit", ws, file_path="in.txt",
                               old_string="inside", new_string="outside"))
    assert seen[0] == ("write_file", {"file_path": "new.txt", "content": "hi",
                                      "path": "new.txt"})
    assert seen[1][0] == "edit_file" and seen[1][1]["new_string"] == "outside"
