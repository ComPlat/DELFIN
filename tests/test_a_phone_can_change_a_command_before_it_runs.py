"""From the phone, a shell approval can be answered with a changed command.

Input: a reply (SimpleX reply gesture) to a shell approval that offers
key c, typed as ``ändern: NEUER BEFEHL`` (also ``change:`` or ``c:``),
optionally in backticks.
Output and semantics:

- the gateway sends ``decide`` with choice c and the command; the bridge
  passes it on as ``permission.reply`` (choice c, command), bound to the
  same live, unchanged request as every phone decision;
- the phone is shown the command exactly as it will run;
- a bare ``c`` or the option's label, an empty command, and a command
  with typographic quotes or dashes (what autocorrect makes of ``"`` and
  ``--``) change nothing and say how to answer;
- spoken, a change frees nothing: a transcript is never a command;
- the bridge refuses a command with any other choice, and choice c
  without a command.
"""
from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms

from test_messenger_scoped_approvals import code_for, scoped  # noqa: F401
from test_the_messenger_gateway_only_serves_the_paired_phone import Bridge, Daemon


SHELL = {"id": "c1", "sessionId": "s1", "kind": "confirm", "tool": "bash",
         "text": "$ rm -r results", "phoneApproval": True, "action": "shell",
         "options": [{"key": "y", "label": "yes"}, {"key": "n", "label": "no"},
                     {"key": "c", "label": "change the command",
                      "detail": "your version runs instead, once, exactly as typed"}]}


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    store.put_room("#10", "s1", "same title", owner="7", ready=True)
    daemon, bridge = Daemon(), Bridge()
    return Rooms(daemon, bridge), daemon, bridge


def _request(rooms, item=SHELL):
    async def run():
        await rooms.send_request("#10", store.rooms()["#10"], item)
        return next(mid for (chat, mid), it in rooms._reply_targets.items()
                    if it["id"] == item["id"])
    return asyncio.run(run())


def _reply(rooms, text, mid):
    return asyncio.run(rooms.reply(text, store.rooms()["#10"], "#10", mid))


def _decisions(bridge):
    return [c for c in bridge.calls if c.get("action") == "decide"]


def test_the_approval_says_how_to_change_the_command(rig):
    rooms, daemon, _ = rig
    _request(rooms)
    assert "ändern:" in daemon.sent[-1][1]


@pytest.mark.parametrize("text", ["ändern: rm -r results/old",
                                  "Ändern:rm -r results/old",
                                  "change: rm -r results/old",
                                  "c: `rm -r results/old`",
                                  "ändern: ```rm -r results/old```"])
def test_a_typed_change_is_sent_with_the_command(rig, text):
    rooms, _, bridge = rig
    mid = _request(rooms)
    out = _reply(rooms, text, mid)
    assert [(c["choice"], c["command"]) for c in _decisions(bridge)] == [
        ("c", "rm -r results/old")]
    assert "rm -r results/old" in out


@pytest.mark.parametrize("text", ["c", "change the command", "ändern:", "ändern:   ",
                                  "ändern: rm -r “results”", "ändern: ls —all",
                                  "ändern: ls –all"])
def test_a_change_that_is_no_usable_command_changes_nothing(rig, text):
    rooms, _, bridge = rig
    mid = _request(rooms)
    out = _reply(rooms, text, mid)
    assert _decisions(bridge) == []
    assert "ändern:" in out


def test_a_change_where_none_is_offered_changes_nothing(rig):
    rooms, _, bridge = rig
    item = {**SHELL, "options": [o for o in SHELL["options"] if o["key"] != "c"]}
    mid = _request(rooms, item)
    _reply(rooms, "ändern: ls", mid)
    assert _decisions(bridge) == []


def test_a_spoken_change_frees_nothing(rig):
    rooms, _, bridge = rig
    mid = _request(rooms)
    for spoken in ("ändern: rm -r results/old", "change the command"):
        asyncio.run(rooms.voice(spoken, store.rooms()["#10"], "#10", mid))
    assert _decisions(bridge) == []


def test_the_bridge_passes_the_change_to_the_live_request(scoped):  # noqa: F811
    bridge, requests, calls = scoped
    requests[0]["options"].append({"key": "c"})

    async def run():
        code = await code_for(bridge)
        base = {"session": "s1", "action": "decide", "id": code}
        assert "error" in await bridge._cmd_session({**base, "choice": "c"})
        assert "error" in await bridge._cmd_session({**base, "choice": "y", "command": "ls"})
        assert (await bridge._cmd_session({**base, "choice": "c", "command": "ls -la"}))["ok"]
    asyncio.run(run())
    assert [p for m, p in calls if m == "permission.reply"] == [
        {"sessionId": "s1", "requestId": "p123-1", "choice": "c", "command": "ls -la",
         "persist": False}]
