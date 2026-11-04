"""`neu NAME [worktree word]` and `zusammenführen [pfad]` from the phone.

Input: the bridge's ``create`` and ``merge`` commands
(delfin/agent/service/bridge.py) and the gateway's ``neu``/``zusammenführen``
(delfin/messenger/gateway.py).

- ``neu NAME arbeitskopie`` (or ``neu NAME in arbeitskopie``) asks for the
  session to run in its own worktree of the new folder; the folder itself
  may first become a local project (initProject) -- the phone's consent is
  in the word it typed, so the bridge passes it on.
- Bare ``neu NAME`` stays as it was: no worktree, no initProject.
- ``zusammenführen PFAD`` merges a session worktree back (session.merge);
  without a path the bridge lists the worktrees (session.worktrees). A
  conflict is reported, nothing overwritten.
- The matchers are the only German part; every literal the phone is shown
  has an English entry (tests/test_the_phone_speaks_english_on_request.py).
"""
from __future__ import annotations

import asyncio

import pytest

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
        if method == "session.worktrees":
            return {"worktrees": [{"path": str(root / "a" / ".delfin" / "worktrees" / "wt1"),
                                   "summary": "the notes task"}]}
        return {"sessionId": f"s{len(calls)}", "session": {"title": params.get("title")}}

    monkeypatch.setattr(bridge, "_request", request)
    return bridge, root, calls


def _create(bridge, **params):
    return asyncio.run(bridge._run_command({"cmd": "create", **params}))


def test_a_worktree_word_makes_the_session_run_in_its_own_worktree(rig):
    bridge, root, calls = rig
    reply = _create(bridge, name="Notizen", worktree=True)
    assert reply.get("ok"), reply
    (method, params), = calls
    assert method == "session.create"
    assert params["worktree"] is True
    assert params["initProject"] is True
    assert params["cwd"] == str(root / "Notizen")


def test_bare_neu_stays_without_a_worktree(rig):
    bridge, root, calls = rig
    reply = _create(bridge, name="Notizen")
    assert reply.get("ok"), reply
    (method, params), = calls
    assert "worktree" not in params
    assert "initProject" not in params


def test_zusammenfuehren_without_a_path_lists_the_worktrees(rig):
    bridge, root, calls = rig
    reply = asyncio.run(bridge._run_command({"cmd": "merge", "path": ""}))
    assert reply.get("ok"), reply
    assert reply.get("worktrees"), reply
    (method, params), = calls
    assert method == "session.worktrees"


def test_zusammenfuehren_merges_a_path_back_and_reports_a_conflict(rig):
    bridge, root, calls = rig

    async def request(method, params=None):
        calls.append((method, dict(params or {})))
        if method == "session.merge":
            return {"status": "conflict", "files": ["notes.md"],
                    "message": "both changed notes.md"}

    monkeypatch = None  # noqa: F841  (the local override is the test)
    bridge._request = request
    reply = asyncio.run(bridge._run_command(
        {"cmd": "merge", "path": str(root / "a" / ".delfin" / "worktrees" / "wt1")}))
    assert reply.get("ok"), reply
    assert reply.get("status") == "conflict", reply
    (method, params), = calls
    assert method == "session.merge"
    assert params["path"].endswith("wt1")


def test_the_gateway_sends_the_worktree_word(tmp_path, monkeypatch):
    """The German matchers on the phone map to the bridge's fields."""
    from delfin.messenger import store
    from delfin.messenger.gateway import Rooms

    class _Bridge:
        def __init__(self):
            self.calls = []

        def ask(self, command):
            self.calls.append(command)
            if command.get("cmd") == "create":
                return {"ok": True, "sessionId": "s-new",
                        "title": command.get("name") or "handy",
                        "folder": command.get("name") or "handy"}
            if command.get("cmd") == "sessions":
                return {"sessions": [{"id": "s-new", "name": "x",
                                      "messengerEnabled": True}]}
            if command.get("action") == "status":
                return {"ok": True, "title": "Notizen", "messengerEnabled": True}
            return {"ok": True}

    class _Daemon:
        async def create_group(self, name):
            return "#20"

        async def invite(self, chat, owner):
            pass

        async def rename_group(self, chat, name):
            pass

        async def send_chat(self, chat, text):
            return {}

    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    bridge = _Bridge()
    out = asyncio.run(Rooms(_Daemon(), bridge).dispatch(
        "neu Notizen arbeitskopie", None, chat="@7"))
    assert {"cmd": "create", "name": "Notizen", "worktree": True} in bridge.calls
    assert "Notizen" in out
    asyncio.run(Rooms(_Daemon(), bridge).dispatch(
        "neu Notizen in arbeitskopie", None, chat="@7"))
    assert {"cmd": "create", "name": "Notizen", "worktree": True} in bridge.calls


def test_the_gateway_passes_the_merge_path(tmp_path, monkeypatch):
    """`zusammenführen PFAD` asks the bridge to merge that worktree back."""
    from delfin.messenger import store
    from delfin.messenger.gateway import Rooms

    calls = []

    class _Bridge:
        def ask(self, command):
            calls.append(command)
            if command.get("cmd") == "merge":
                return {"ok": True, "status": "ok", "target": "main",
                        "files": ["notes.md"]}
            return {"ok": True}

    class _Daemon:
        async def send_chat(self, chat, text):
            return {}

    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    rooms = Rooms(_Daemon(), _Bridge())
    out = asyncio.run(rooms.dispatch(
        "zusammenführen " + str(tmp_path / "wt1"), None, chat="@7"))
    assert calls and calls[0]["cmd"] == "merge", calls
    assert calls[0]["path"].endswith("wt1")
    assert "notes.md" in out, out
