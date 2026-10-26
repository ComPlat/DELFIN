"""A session can be started from the phone, in a folder of its own.

Input: the bridge's ``create`` command (delfin/agent/service/bridge.py) and
the gateway's ``neu NAME`` (delfin/messenger/gateway.py).

- The folder is chosen by the host, never by the phone: one directory per
  name directly under ``DELFIN_MESSENGER_WORKSPACES``; a name is reduced to
  letters, digits, ``.``, ``_`` and ``-``, so ``../x`` or ``/etc`` cannot
  leave the root, and a symlink in the root is not followed.
- The session starts with the permission profile ``default`` (it asks
  before commands and edits), whatever the phone sends; no model, mode or
  folder field is accepted from the gateway.
- It is enabled for the messenger at once (the activation record), so its
  chat follows without a ``/sync`` at the computer.
- Without ``DELFIN_MESSENGER_WORKSPACES`` the host has not allowed it, and
  nothing is created.
- At most three starts per ten minutes.
"""
from __future__ import annotations

import asyncio
import os

import pytest

from delfin.agent.service import bridge as bridge_mod, protocol
from delfin.agent.service.bridge import Bridge
from delfin.messenger import activation


@pytest.fixture
def rig(tmp_path, monkeypatch):
    root = tmp_path / "phone-sessions"
    root.mkdir()
    monkeypatch.setenv("DELFIN_MESSENGER_WORKSPACES", str(root))
    monkeypatch.setattr(activation, "_DIR", tmp_path / "activation")
    bridge = Bridge(allowed_uid=1, service_path=tmp_path / "none.sock")
    calls = []

    async def request(method, params=None):
        calls.append((method, dict(params or {})))
        return {"sessionId": f"s{len(calls)}", "session": {"title": params.get("title")}}

    monkeypatch.setattr(bridge, "_request", request)
    return bridge, root, calls


def _create(bridge, **params):
    return asyncio.run(bridge._run_command({"cmd": "create", **params}))


def test_a_named_session_gets_its_own_folder_and_asks_first(rig):
    bridge, root, calls = rig
    reply = _create(bridge, name="Butan Konformere")
    assert reply.get("ok"), reply
    (method, params), = calls
    assert method == "session.create"
    assert params["cwd"] == str(root / "Butan-Konformere")
    assert (root / "Butan-Konformere").is_dir()
    assert params["permissionMode"] == "default"
    assert params["title"] == "Butan Konformere"
    assert set(params) <= set(protocol.SERVICE_REQUESTS["session.create"]) | {
        k.rstrip("?") for k in protocol.SERVICE_REQUESTS["session.create"]}
    assert activation.enabled(reply["sessionId"])


def test_a_name_cannot_leave_the_root(rig):
    bridge, root, calls = rig
    for name in ("../../etc", "/etc/passwd", "..", ".hidden"):
        _create(bridge, name=name)
    for _, params in calls:
        cwd = os.path.realpath(params["cwd"])
        assert os.path.dirname(cwd) == os.path.realpath(root), params["cwd"]
        assert not os.path.basename(cwd).startswith(".")


def test_a_symlink_in_the_root_is_not_followed(rig, tmp_path):
    bridge, root, calls = rig
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (root / "link").symlink_to(elsewhere)
    reply = _create(bridge, name="link")
    assert reply.get("error") and not calls


def test_an_unnamed_session_still_gets_a_fresh_folder(rig):
    bridge, root, calls = rig
    assert _create(bridge).get("ok")
    (_, params), = calls
    assert os.path.dirname(params["cwd"]) == str(root)
    assert os.path.basename(params["cwd"]).startswith("handy-")


def test_the_gateway_cannot_choose_folder_mode_or_model(rig):
    bridge, _, calls = rig
    for extra in ({"cwd": "/"}, {"permissionMode": "bypassPermissions"},
                  {"model": "x"}, {"mode": "office"}):
        reply = _create(bridge, name="a", **extra)
        assert reply.get("error"), extra
    assert calls == []


def test_without_a_configured_root_nothing_is_started(rig, monkeypatch):
    bridge, _, calls = rig
    monkeypatch.delenv("DELFIN_MESSENGER_WORKSPACES")
    reply = _create(bridge, name="a")
    assert reply.get("error") and calls == []


def test_starts_are_rate_limited(rig):
    bridge, _, calls = rig
    replies = [_create(bridge, name=f"n{i}") for i in range(5)]
    assert [bool(r.get("ok")) for r in replies] == [True, True, True, False, False]
    assert len(calls) == 3


def test_the_command_is_in_the_bridges_fixed_set():
    assert "create" in bridge_mod.COMMANDS


# --- gateway: `neu NAME` -----------------------------------------------------

from delfin.messenger import store  # noqa: E402
from delfin.messenger.gateway import Rooms  # noqa: E402


class _Bridge:
    def __init__(self):
        self.calls = []

    def ask(self, command):
        self.calls.append(command)
        if command.get("cmd") == "create":
            return {"ok": True, "sessionId": "s-new", "title": command.get("name") or "handy"}
        if command.get("cmd") == "sessions":
            return {"sessions": [{"id": "s-new", "name": "x", "messengerEnabled": True}]}
        if command.get("action") == "status":
            return {"ok": True, "title": "Butan", "messengerEnabled": True}
        return {"ok": True}


class _Daemon:
    def __init__(self):
        self.groups, self.invites = [], []

    async def create_group(self, name):
        self.groups.append(name)
        return "#20"

    async def invite(self, chat, owner):
        self.invites.append((chat, owner))

    async def rename_group(self, chat, name):
        pass

    async def send_chat(self, chat, text):
        return {}


def test_neu_starts_a_session_and_invites_to_its_chat(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    daemon, bridge = _Daemon(), _Bridge()
    rooms = Rooms(daemon, bridge)
    out = asyncio.run(rooms.dispatch("neu Butan", None, chat="@7"))
    assert {"cmd": "create", "name": "Butan"} in bridge.calls
    assert daemon.invites == [("#20", "7")]
    assert store.rooms()["#20"]["session_id"] == "s-new"
    assert "Butan" in out


def test_neu_with_a_task_sends_it_as_the_first_prompt(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    bridge = _Bridge()
    rooms = Rooms(_Daemon(), bridge)
    asyncio.run(rooms.dispatch("neu Butan: Optimiere beide Konformere.", None, chat="@7"))
    assert {"cmd": "create", "name": "Butan"} in bridge.calls
    prompts = [c for c in bridge.calls if c.get("action") == "prompt"]
    assert [(c["session"], c["text"]) for c in prompts] == [
        ("s-new", "Optimiere beide Konformere.")]
