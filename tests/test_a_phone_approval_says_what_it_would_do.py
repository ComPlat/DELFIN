"""A phone approval says in plain words what it would do, then shows it.

Input: the approval broker's request description
(delfin/agent/service/approvals.py), the bridge's pending items
(delfin/agent/service/bridge.py) and the gateway's text
(delfin/messenger/gateway.py).

Before, the phone showed the tool name and the raw preview: readable for
someone who reads shell, not for anyone else. The bash tool already makes
the model state a one-line purpose (its required ``description``); it
never reached a window.

- ``permission.asked`` carries ``purpose``: the tool call's own
  ``description``, masked like the preview, one line, at most 200
  characters; absent when the call has none.
- The phone shows, in this order: what kind of action it is in plain
  words (a command on this machine, a file change, a read outside the
  folder), the purpose labelled as the agent's own unchecked statement,
  and then the full command or change -- the purpose never replaces it.
"""
from __future__ import annotations

import asyncio

from delfin.agent import terminal_confirm as tc
from delfin.agent.service import protocol
from delfin.agent.service.approvals import ProtocolBroker
from delfin.agent.service.bridge import Bridge
from delfin.messenger import activation
from delfin.messenger.gateway import Rooms


def _describe(args, tool="bash", preview="$ ls -la"):
    broker = ProtocolBroker.__new__(ProtocolBroker)
    broker.suggest_pattern = lambda command: ""
    broker.request_id = lambda req: "r1"
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool=tool, args=args, preview=preview)
    return broker.describe(req)


def test_the_tool_calls_purpose_travels_with_the_request():
    fields = _describe({"command": "ls -la", "description": "Listet die Dateien im Ordner auf"})
    assert fields["purpose"] == "Listet die Dateien im Ordner auf"
    assert "purpose" in protocol.EVENTS["permission.asked"]


def test_no_purpose_no_field_and_one_bounded_line():
    assert "purpose" not in _describe({"command": "ls"})
    long = _describe({"command": "ls", "description": "a\nb " + "x" * 400})
    assert "\n" not in long["purpose"] and len(long["purpose"]) <= 200


def test_a_secret_in_the_purpose_is_masked():
    secret = "sk-ant-api03-" + "A" * 40
    fields = _describe({"command": "true", "description": f"use key {secret}"})
    assert secret not in fields["purpose"]


def test_the_bridge_hands_purpose_and_kind_to_the_phone(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, "_DIR", tmp_path / "activation")
    activation.configure("s1")
    bridge = Bridge(allowed_uid=1, service_path=tmp_path / "none.sock")
    asked = {"requestId": "r1", "kind": "confirm", "tool": "bash",
             "preview": "$ rm -rf build", "purpose": "Räumt den Build-Ordner auf",
             "flags": {"shell": True},
             "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}

    async def request(method, params=None):
        if method == "sessions.list":
            return {"sessions": [{"sessionId": "s1", "workspace": str(tmp_path)}]}
        return {"pending": [asked]}

    monkeypatch.setattr(bridge, "_request", request)
    reply = asyncio.run(bridge._run_command(
        {"cmd": "session", "session": "s1", "action": "pending"}))
    item, = reply["items"]
    assert item["purpose"] == "Räumt den Build-Ordner auf"
    assert item["action"] == "shell"
    assert "rm -rf build" in item["text"]


def test_the_phone_text_explains_then_shows_the_command(monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LABEL", "Pi")
    text = Rooms.pending_text([{
        "kind": "confirm", "tool": "bash", "action": "shell",
        "purpose": "Räumt den Build-Ordner auf", "text": "$ rm -rf build",
        "phoneApproval": True,
        "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"}]}])
    what = text.index("Ein Befehl wird auf DELFIN-Pi ausgeführt")
    purpose = text.index("Angabe des Agenten (ungeprüft): Räumt den Build-Ordner auf")
    command = text.index("$ rm -rf build")
    assert what < purpose < command


def test_without_a_purpose_the_command_is_still_explained(monkeypatch):
    monkeypatch.delenv("DELFIN_MESSENGER_LABEL", raising=False)
    text = Rooms.pending_text([{
        "kind": "confirm", "tool": "write_file", "action": "write",
        "text": "notes.md: +3 lines", "phoneApproval": True, "options": []}])
    assert "Eine Datei wird geändert" in text
    assert "ungeprüft" not in text
    assert "notes.md: +3 lines" in text
