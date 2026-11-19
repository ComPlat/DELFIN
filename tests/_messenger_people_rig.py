"""A people-mode account for the messenger tests, without a network.

The account home and two persons (anna, ben) in tmp_path, people mode on;
the REAL bridge (delfin.agent.service.bridge.Bridge) with one fake agent
service per person behind it; the REAL gateway (delfin.messenger.gateway.
Rooms) talking to that bridge in-process; a fake SimpleX daemon that
records what it is told. Nothing reaches a daemon, a socket of another
process, an external agent CLI or a network endpoint.
"""
from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path

from delfin.agent.service.bridge import Bridge
from delfin.agent.service.redact import Redactor
from delfin.messenger import activation, people, store
from delfin.messenger.gateway import Rooms

CONFIRM = {"kind": "confirm", "tool": "bash", "preview": "$ rm -r results",
           "flags": {"shell": True}, "options": [{"key": "y", "label": "Once"},
                                                  {"key": "n", "label": "No"}]}


class FakeService:
    """One person's agent service: their sessions, their open requests."""

    def __init__(self, person, sessions):
        self.person = person
        self.sessions = [dict(s) for s in sessions]
        self.pending: dict[str, list] = {}
        self.asked: list[tuple[str, dict]] = []
        self.closed = False

    def request(self, method, params=None, timeout=30.0):
        params = params or {}
        self.asked.append((method, params))
        if method == "sessions.list":
            return {"sessions": [dict(s) for s in self.sessions], "elsewhere": []}
        if method == "session.pending":
            return {"pending": list(self.pending.get(params.get("sessionId"), []))}
        if method == "permission.reply":
            reqs = self.pending.get(params.get("sessionId"), [])
            hit = [r for r in reqs if r["requestId"] == params.get("requestId")]
            for r in hit:
                reqs.remove(r)
            return {"applied": bool(hit), "message": ""}
        if method == "session.create":
            sid = f"{self.person or 'owner'}-new-{len(self.sessions)}"
            self.sessions.append({"sessionId": sid, "title": params.get("title"),
                                  "workspace": params.get("cwd"), "state": "idle"})
            return {"sessionId": sid}
        if method == "session.setTitle":
            for s in self.sessions:
                if s["sessionId"] == params.get("sessionId"):
                    s["title"] = params.get("title")
            return {"title": params.get("title")}
        return {"accepted": True}

    def calls(self, method):
        return [p for m, p in self.asked if m == method]

    def close(self):
        self.closed = True


class FakeDaemon:
    def __init__(self):
        self.sent: list[tuple[str, str]] = []
        self.groups: list[str] = []
        self.invited: list[tuple[str, str]] = []
        self.renamed: list[tuple[str, str]] = []
        self.removed: list[tuple[str, str]] = []
        self.deleted: list[str] = []
        self.dropped: list[str] = []
        self.links: list[str] = []
        self.updated: list[tuple[str, str, str]] = []
        self.quotes: list[tuple[str, str]] = []
        self.profile_commands: list = []
        self.contact_menus: list[tuple[str, list]] = []
        self.group_menus: list[tuple[str, list]] = []
        self._item = 0

    async def send_chat(self, chat, text, quote=None):
        self.sent.append((chat, text))
        if quote is not None:
            self.quotes.append((chat, str(quote)))
        self._item += 1
        return {"chatItems": [{"chatItem": {"meta": {"itemId": self._item}}}]}

    async def set_bot_profile(self, commands):
        self.profile_commands = list(commands)

    async def set_contact_menu(self, cid, commands):
        self.contact_menus.append((str(cid), list(commands)))

    async def enable_receipts(self):
        pass

    async def update_message(self, chat, item_id, text):
        self.updated.append((chat, str(item_id), text))

    async def create_group(self, title):
        self.groups.append(title)
        return "#" + str(500 + len(self.groups))

    async def invite(self, chat, cid):
        self.invited.append((chat, cid))

    async def rename_group(self, chat, title, commands=None):
        self.renamed.append((chat, title))
        if commands is not None:
            self.group_menus.append((chat, list(commands)))

    async def new_invitation(self):
        n = len(self.links) + 1
        link = f"simplex:/invitation#/?v=2&smp=SECRETKEYMATERIAL{n:04d}"
        self.links.append(link)
        return link, str(900 + n)

    async def contact_code(self, cid):
        return "12345 67890 11111"

    async def remove_member(self, chat, cid):
        self.removed.append((chat, cid))

    async def leave_group(self, chat):
        pass

    async def delete_contact(self, cid):
        self.deleted.append(str(cid))

    async def delete_pending(self, conn):
        self.dropped.append(str(conn))

    def to(self, chat):
        return [text for c, text in self.sent if c == chat]


class SyncBridge:
    """The gateway's BridgeClient, answered by the real bridge in-process."""

    def __init__(self, bridge):
        self.bridge = bridge
        self.commands: list[dict] = []

    def ask(self, command):
        self.commands.append(dict(command))
        return asyncio.run(self.bridge._run_command(json.loads(json.dumps(command))))


class Account:
    """The account home with people mode on, persons anna and ben."""

    def __init__(self, tmp_path: Path, monkeypatch):
        self.home = tmp_path / "account"
        (self.home / ".config" / "delfin").mkdir(parents=True)
        (self.home / ".config" / "delfin" / "launcher-people.conf").write_text("on\n")
        for name in ("anna", "ben"):
            (self.home / "people" / name / ".delfin").mkdir(parents=True, mode=0o700)
        handy = tmp_path / "DELFIN-Handy"
        handy.mkdir()
        monkeypatch.setattr(people, "account_home", lambda: self.home)
        # Each person's private runtime directory (their pairing doors, the
        # bridge's delivery check), short enough for a Unix socket path.
        import shutil
        import tempfile
        import weakref
        self.runtime = {name: Path(tempfile.mkdtemp(prefix="dpair-", dir="/tmp"))
                        for name in ("", "anna", "ben")}
        for path in self.runtime.values():
            weakref.finalize(self, shutil.rmtree, str(path), True)
        monkeypatch.setattr(people, "door_root", lambda name, account=None: self.runtime.get(name), raising=False)
        from delfin.messenger import pairing as _pairing
        pairing_dirs = {name: (self.home if not name else self.home / "people" / name)
                        / ".delfin" / "messenger_pairing" for name in self.runtime}
        monkeypatch.setattr(_pairing, "_door_dir", lambda directory=None: next(
            (self.runtime[n] for n, d in pairing_dirs.items() if directory is not None
             and Path(directory) == d), self.runtime[""]), raising=False)
        monkeypatch.setattr(activation, "_DIR", self.home / ".delfin" / "messenger_sessions")
        monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "gateway-state"))
        monkeypatch.setenv("DELFIN_MESSENGER_WORKSPACES", str(handy))
        monkeypatch.setattr(store, "_REFUSED_SEEN", set())
        redactor = Redactor()
        redactor._known = []
        redactor._loaded = float("inf")
        self.bridge = Bridge(allowed_uid=os.getuid() + 1, redactor=redactor)
        def ws(name, person=""):
            # A session's folder in its person's own tree (hol reaches
            # only those): the account's for the owner, ~/people/<p> else.
            base = self.home if not person else self.home / "people" / person
            path = base / "work" / name
            path.mkdir(parents=True, exist_ok=True)
            return str(path)

        self.services = {
            "": FakeService("", [{"sessionId": "s-own", "title": "Eigene", "workspace": ws("own"),
                                  "state": "idle"}]),
            "anna": FakeService("anna", [
                {"sessionId": "s-anna", "title": "Butan", "workspace": ws("butan", "anna"), "state": "idle"},
                {"sessionId": "s-anna2", "title": "Ethan", "workspace": ws("ethan", "anna"), "state": "idle"}]),
            "ben": FakeService("ben", [{"sessionId": "s-ben", "title": "Benzol",
                                        "workspace": ws("benzol", "ben"), "state": "idle"}]),
        }
        self.bridge._service = self.services[""]
        for name in ("anna", "ben"):
            self.bridge._person_bridge(name)._service = self.services[name]
        self.daemon = FakeDaemon()
        self.client = SyncBridge(self.bridge)
        self.rooms = Rooms(self.daemon, self.client)

    def activation_dir(self, person):
        return people.activation_dir(person) if person else activation._DIR

    def enable(self, person, sid, args=""):
        return activation.configure(sid, args, directory=self.activation_dir(person))

    def pairing_dir(self, person):
        return people.pairing_dir(person)

    def say(self, sender, chat, text, mid=None, quote=None):
        self._mid = getattr(self, "_mid", 0) + 1
        message = {"chat": chat, "sender": sender, "text": text, "id": mid or str(self._mid)}
        if quote is not None:
            message["quote_id"] = str(quote)
        asyncio.run(self.rooms.handle(message))

    def connect(self, cid, conn, name="Phone"):
        asyncio.run(self.rooms.contact_connected({
            "type": "contactConnected",
            "contact": {"contactId": int(cid), "activeConn": {"connId": int(conn)},
                        "profile": {"displayName": name}}}))

    def round(self):
        asyncio.run(self.rooms.pairing_round())
