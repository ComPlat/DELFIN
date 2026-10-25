"""The simplex-chat daemon connection: its local WebSocket API.

Only the daemon that runs on this host, on 127.0.0.1 (the provisioning
opens that port to the gateway's user alone, nftables meta skuid). The
gateway subscribes to chat events, sends agent messages through the
paired contact, and pairs with a one-time connection invitation.

The API surface used here is narrow on purpose; everything else the
daemon says is ignored. What the daemon says before it is subscribed is
not trusted either: a contact request or a message from anyone but the
paired contact is refused by the gateway, not by the daemon alone.
"""

from __future__ import annotations

import json
import os

import websockets

DEFAULT_PORT = 5225


def daemon_port() -> int:
    raw = os.environ.get("DELFIN_MESSENGER_PORT", "").strip()
    try:
        return int(raw) if raw else DEFAULT_PORT
    except ValueError:
        return DEFAULT_PORT


def ws_url(port: int | None = None) -> str:
    """The daemon's local API URL. SimpleX v7's websocket server speaks
    plain ws:// on 127.0.0.1 only."""
    return f"ws://127.0.0.1:{port or daemon_port()}"


class Daemon:
    """One connection to the daemon, as JSON over WebSocket."""

    def __init__(self, port: int | None = None):
        self.url = ws_url(port)
        self._ws = None
        self._corr_id = 0

    async def __aenter__(self) -> "Daemon":
        self._ws = await websockets.connect(self.url)
        return self

    async def __aexit__(self, *exc) -> None:
        if self._ws is not None:
            await self._ws.close()

    async def command(self, params: dict) -> dict:
        """One command, one response, matched by corrId."""
        self._corr_id += 1
        corr = f"delfin-{self._corr_id}"
        payload = {"corrId": corr, "cmd": json.dumps(params, ensure_ascii=False)}
        await self._ws.send(json.dumps(payload))
        while True:
            raw = await self._ws.recv()
            try:
                message = json.loads(raw)
            except (TypeError, ValueError):
                continue
            if message.get("corrId") != corr:
                continue          # an event or another client's response
            response = message.get("resp")
            if response is None:
                continue
            if isinstance(response, str):
                try:
                    response = json.loads(response)
                except ValueError:
                    return {"type": "unknown", "text": response}
            return response

    async def events(self):
        """Everything the daemon pushes that is not a command response;
        the caller filters by contact."""
        while True:
            raw = await self._ws.recv()
            try:
                message = json.loads(raw)
            except (TypeError, ValueError):
                continue
            if message.get("corrId") is None:
                response = message.get("resp")
                if isinstance(response, str):
                    try:
                        response = json.loads(response)
                    except ValueError:
                        continue
                yield response

    # -- the narrow API surface ----------------------------------------------------------

    async def subscribe(self) -> dict:
        return await self.command({"type": "subscribe"})

    async def contacts(self) -> dict:
        return await self.command({"type": "getContacts"})

    async def send_text(self, contact_id: str, text: str) -> dict:
        return await self.command(
            {"type": "sendMessages", "contactId": contact_id,
             "messageText": text})

    async def pending_invitation(self) -> str | None:
        """The connection invitation the daemon generated for pairing
        (chat.exploreContact INVITATION response), when there is one."""
        reply = await self.command({"type": "createLink"})
        if isinstance(reply, dict):
            conn = reply.get("connection") or {}
            return conn.get("connectionCode") or reply.get("connection")
        return None


def extract_contact_list(reply: dict) -> list[dict]:
    """The contacts out of a getContacts response, however the daemon
    wraps them (chat.contacts info array)."""
    if not isinstance(reply, dict):
        return []
    data = reply.get("chatItems") or reply.get("contacts") or []
    if isinstance(data, dict):
        data = data.get("contacts") or []
    return [c for c in data if isinstance(c, dict)]


def contact_id_of(contact: dict) -> str | None:
    return str(contact.get("contactId") or "") or None


def contact_name_of(contact: dict) -> str:
    profile = contact.get("localAlias") or contact.get("profile") or {}
    if isinstance(profile, dict):
        return str(profile.get("displayName") or contact.get("localAlias") or "")
    return str(contact.get("localAlias") or "")


def new_contact_id(event: dict) -> str | None:
    """The contact id in a chat.contact.connected / newContact event, so
    ``pair`` can catch the moment the phone scanned the link."""
    info = event.get("user") or event.get("contact") or {}
    if isinstance(info, dict):
        return str(info.get("contactId") or "") or None
    return None


def text_of(event: dict) -> tuple[str, str] | None:
    """(contactId, text) out of a chat.contact.switchTo / newChatItems
    event's message, when it carries one."""
    items = (event.get("chatItems")
             or (event.get("user") and [])
             or [])
    if not isinstance(items, list):
        return None
    for item in items:
        content = (item.get("chatItem") or item).get("msgContent") or {}
        if isinstance(content, dict) and content.get("type") == "text":
            cid = contact_id_of((item.get("chatItem") or item)
                                .get("contact") or {})
            text = str(content.get("text") or "")
            if cid and text:
                return cid, text
    return None
