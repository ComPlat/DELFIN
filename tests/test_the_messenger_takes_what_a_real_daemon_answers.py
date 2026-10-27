"""The gateway against what the pinned SimpleX daemon (v7.0.3) really answers.

Input: the replies measured on a Raspberry Pi with the real daemon and a
paired phone. Three of them broke the first live setup:

- A group display name longer than 50 characters is refused with
  ``chatCmdError`` / ``invalidDisplayName`` and a ``validName`` (the name cut
  to 50). A session without a title was named after its id, twice, behind
  the bot name: 63 characters, so no session chat could be created.
- That refusal reached the phone as "Nicht ausgeführt oder Verbindung
  unterbrochen" and nowhere else: neither the reason nor a log line.
- Pairing connected to the daemon the moment its unit was started, before
  it listened, and failed with ConnectionRefusedError.

Output: names within the limit that say machine and session; a refused
name retried once with the daemon's own ``validName``; the reason in the
log and, in words, on the phone; pairing that waits for the daemon, prints
the link last (a phone screen shows the end of the output), shows the
security code, says hello over the new contact and opens the forum for
all agents, so it is there before the first session is enabled.
"""
from __future__ import annotations

import asyncio
import json
import socket

import pytest

from delfin.messenger import simplex, store
from delfin.messenger.bridge_client import BridgeDown
from delfin.messenger.gateway import Rooms, group_name

LIMIT = 50


def test_a_session_chat_name_fits_the_daemons_limit(monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LABEL", "x" * 24)
    sid = "2fdfa16304e944b2b72fd41201685302"
    for title in ("", sid, "test", "t" * 120):
        name = group_name(title or sid, sid)
        assert 0 < len(name) <= LIMIT, name
    assert len(group_name("Forum", None)) <= LIMIT


def test_a_session_chat_says_machine_and_session_name(monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LABEL", "delfin-pi")
    name = group_name("test", "2fdfa16304e944b2b72fd41201685302")
    assert "delfin-pi" in name and "test" in name
    assert "2fdfa16304e944b2b72fd41201685302" not in name


class _Socket:
    """The daemon's WebSocket: each command is answered with the next reply."""

    def __init__(self, daemon, replies):
        self.daemon, self.replies, self.sent = daemon, list(replies), []

    async def send(self, raw):
        payload = json.loads(raw)
        self.sent.append(payload["cmd"])
        reply = self.replies.pop(0)
        self.daemon._pending[payload["corrId"]].set_result(reply)


def _daemon(replies):
    daemon = simplex.Daemon(port=1)
    daemon._ws = _Socket(daemon, replies)
    return daemon


USER = {"type": "activeUser", "user": {"userId": 1}}
REFUSED = {"type": "chatCmdError", "chatError": {"type": "error", "errorType": {
    "type": "invalidDisplayName", "displayName": "D" * 52, "validName": "D" * 50}}}


def test_a_refused_name_is_retried_with_the_daemons_valid_name():
    daemon = _daemon([USER, REFUSED, USER, {"type": "groupCreated",
                                            "groupInfo": {"groupId": 9}}])
    assert asyncio.run(daemon.create_group("D" * 52)) == "#9"
    assert json.loads(daemon._ws.sent[-1].split(" ", 2)[2])["displayName"] == "D" * 50


def test_a_refusal_names_its_reason():
    daemon = _daemon([{"type": "chatCmdError", "chatError": {
        "type": "errorAgent", "agentError": {"type": "SMP"}}}])
    with pytest.raises(simplex.SimplexError) as raised:
        asyncio.run(daemon.command("/user"))
    assert "errorAgent" in str(raised.value)


class _FailingDaemon:
    def __init__(self, error):
        self.error, self.sent = error, []

    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {}

    async def create_group(self, title):
        raise self.error

    async def invite(self, chat, owner):
        pass


class _Bridge:
    def __init__(self, title="", down=False):
        self.title, self.down = title, down

    def ask(self, command):
        if self.down:
            raise BridgeDown("/run/delfin-messenger-bridge/messenger.sock: refused")
        if command["cmd"] == "sessions":
            return {"sessions": [{"id": "s1", "name": "", "state": "idle",
                                  "messengerEnabled": True}]}
        return {"ok": True, "title": self.title, "folder": "test", "state": "idle",
                "messengerEnabled": True}


def _rooms(tmp_path, monkeypatch, daemon, bridge):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    store.put_room("#1", None, "Forum", owner="7", ready=True)
    return Rooms(daemon, bridge)


def _sync(rooms):
    asyncio.run(rooms.handle({"chat": "@7", "text": "/sync", "sender": "7", "id": "1"}))


def test_a_simplex_refusal_is_logged_and_told_in_words(tmp_path, monkeypatch, capsys):
    daemon = _FailingDaemon(simplex.SimplexError("invalidDisplayName"))
    _sync(_rooms(tmp_path, monkeypatch, daemon, _Bridge()))
    told = daemon.sent[-1][1]
    assert "SimpleX" in told and "invalidDisplayName" in told
    assert "invalidDisplayName" in capsys.readouterr().err


def test_an_unreachable_bridge_is_told_as_such(tmp_path, monkeypatch, capsys):
    daemon = _FailingDaemon(AssertionError("unused"))
    _sync(_rooms(tmp_path, monkeypatch, daemon, _Bridge(down=True)))
    told = daemon.sent[-1][1]
    assert "delfin-messenger-bridge" in told
    assert "/run/" not in told                          # no paths to the phone
    assert "BridgeDown" in capsys.readouterr().err


def test_an_untitled_session_is_named_after_its_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LABEL", "delfin-pi")

    class Daemon(_FailingDaemon):
        async def create_group(self, title):
            self.created = title
            return "#5"

    daemon = Daemon(None)
    rooms = _rooms(tmp_path, monkeypatch, daemon, _Bridge(title=""))
    asyncio.run(rooms.ensure_room("s1"))
    assert "test" in daemon.created and "s1" not in daemon.created.split()


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_the_daemon_connection_waits_for_a_daemon_still_starting():
    import websockets

    port = _free_port()

    async def main():
        async def late():
            await asyncio.sleep(0.6)
            return await websockets.serve(lambda ws: ws.wait_closed(), "127.0.0.1", port)

        starting = asyncio.create_task(late())
        async with simplex.Daemon(port=port, wait=10):
            pass
        (await starting).close()

    asyncio.run(main())


def test_pairing_waits_prints_the_link_last_shows_the_code_and_says_hello(
        tmp_path, monkeypatch, capsys):
    from delfin.messenger import __main__ as cli
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    monkeypatch.setattr(cli, "_show_qr", lambda link: print("QRCODE"))
    said = {}

    class Daemon:
        def __init__(self, *args, **kwargs):
            said["wait"] = kwargs.get("wait")
        async def __aenter__(self):
            return self
        async def __aexit__(self, *args):
            pass
        async def subscribe(self):
            pass
        async def pending_invitation(self):
            self.invitation_id = "conn-1"
            return "simplex:/invitation#/?v=x"
        async def events(self):
            yield {"type": "contactConnected",
                   "contact": {"contactId": 3, "activeConn": {"connId": "conn-1"}}}
        async def contact_code(self, contact_id):
            return "12345 67890"
        async def send_text(self, contact_id, text):
            said["hello"] = (contact_id, text)
        async def create_group(self, title):
            said.setdefault("groups", []).append(title)
            return "#4"
        async def invite(self, chat, contact_id):
            said["invited"] = (chat, contact_id)

    monkeypatch.setattr(cli.simplex, "Daemon", Daemon)
    assert cli._pair() == 0
    out = capsys.readouterr().out
    assert said["wait"] and said["wait"] >= 30
    assert out.index("QRCODE") < out.index("simplex:/invitation")
    assert "12345 67890" in out
    assert said["hello"][0] == "3"
    assert store.paired_contact() == "3"
    # The forum for all agents is there from the start, the phone invited.
    assert said["groups"] and said["groups"][0].startswith("DELFIN Forum")
    assert said["invited"] == ("#4", "3")
    assert [r["session_id"] for r in store.rooms().values()] == [None]
