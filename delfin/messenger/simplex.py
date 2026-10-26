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

import asyncio
import json
import os

import websockets

DEFAULT_PORT = 5227


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
        self._pending = {}
        self._events = asyncio.Queue()

    async def __aenter__(self) -> "Daemon":
        self._ws = await websockets.connect(self.url)
        self._reader = asyncio.create_task(self._read())
        return self

    async def __aexit__(self, *exc) -> None:
        if self._ws is not None:
            await self._ws.close()
            self._reader.cancel()
            await asyncio.gather(self._reader, return_exceptions=True)

    async def command(self, params: str) -> dict:
        """One command, one response, matched by corrId."""
        self._corr_id += 1
        corr = f"delfin-{self._corr_id}"
        payload = {"corrId": corr, "cmd": params}
        future = asyncio.get_running_loop().create_future()
        self._pending[corr] = future
        try:
            await self._ws.send(json.dumps(payload))
            reply = await asyncio.wait_for(future, 30)
            if reply.get("type") == "chatCmdError":
                raise RuntimeError("SimpleX rejected the command")
            return reply
        finally:
            self._pending.pop(corr, None)

    async def _read(self):
        error = ConnectionError("SimpleX connection closed")
        try:
            async for raw in self._ws:
                message = json.loads(raw)
                response = message.get("resp")
                if isinstance(response, str):
                    response = json.loads(response)
                if not isinstance(response, dict):
                    continue
                corr = message.get("corrId")
                if corr is None:
                    await self._events.put(response)
                else:
                    future = self._pending.get(corr)
                    if future is not None and not future.done():
                        future.set_result(response)
        except Exception as exc:
            error = exc
        finally:
            for future in self._pending.values():
                if not future.done():
                    future.set_exception(error)
            await self._events.put(error)

    async def events(self):
        while True:
            event = await self._events.get()
            if isinstance(event, Exception):
                raise event
            yield event

    # -- the narrow API surface ----------------------------------------------------------

    async def subscribe(self) -> dict:
        # CLI WebSockets push events automatically. /user checks readiness.
        return await self.command("/user")

    async def user_id(self) -> str:
        reply = await self.command("/user")
        return str(int(reply["user"]["userId"]))

    async def contacts(self) -> dict:
        return await self.command(f"/_contacts {await self.user_id()}")

    async def send_text(self, contact_id: str, text: str) -> dict:
        return await self.send_chat("@" + str(int(contact_id)), text)

    async def send_chat(self, chat: str, text: str) -> dict:
        if chat[:1] not in {"@", "#"} or not chat[1:].isdigit():
            raise ValueError("invalid chat reference")
        messages = [{"msgContent": {"type": "text", "text": text}, "mentions": {}}]
        return await self.command(f"/_send {chat} json " + json.dumps(messages, ensure_ascii=False))

    async def receive_file(self, file_id, path):
        return await self.command(f"/freceive {int(file_id)} encrypt=off {path}")

    async def send_file(self, chat, path, caption=''):
        if chat[:1] not in {'@', '#'} or not chat[1:].isdigit():
            raise ValueError('invalid chat reference')
        messages = [{'fileSource': {'filePath': str(path)},
                     'msgContent': {'type': 'file', 'text': caption}, 'mentions': {}}]
        return await self.command(f'/_send {chat} json ' + json.dumps(messages))

    async def create_group(self, title: str) -> str:
        profile = {"displayName": title, "fullName": title}
        reply = await self.command(f"/_group {await self.user_id()} " + json.dumps(profile))
        return "#" + str(int(reply["groupInfo"]["groupId"]))

    async def invite(self, chat: str, contact_id: str) -> None:
        await self.command(f"/_add #{int(chat[1:])} {int(contact_id)} member")

    async def rename_group(self, chat: str, title: str) -> None:
        await self.command(f"/_group_profile #{int(chat[1:])} " +
                           json.dumps({'displayName': title, 'fullName': title}))

    async def pending_invitation(self) -> str | None:
        reply = await self.command(f"/_connect {await self.user_id()}")
        self.invitation_id = (reply.get("connection") or {}).get("pccConnId")
        return (reply.get("connLinkInvitation") or {}).get("connFullLink")


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
    if event.get("type") != "contactConnected":
        return None
    return contact_id_of(event.get("contact") or {})


def messages_of(event: dict) -> list[dict]:
    """Received messages/file completions, preserving sender, chat and item ID."""
    complete = event.get('type') == 'rcvFileComplete'
    if event.get("type") not in {"newChatItems", "rcvFileDescrReady"} and not complete:
        return []
    messages = []
    wrappers = event.get('chatItems') or [] if event.get('type') == 'newChatItems' else [event.get('chatItem') or {}]
    for wrapped in wrappers:
        info = wrapped.get("chatInfo") or {}
        item = wrapped.get("chatItem") or {}
        direction = item.get("chatDir") or {}
        content = item.get("content") or {}
        msg = content.get("msgContent") or {}
        meta = item.get("meta") or {}
        if content.get("type") != "rcvMsgContent" or msg.get("type") not in {"text", "image", "file", "voice", "video"}:
            continue
        if meta.get("itemDeleted") or meta.get("itemId") is None:
            continue
        if info.get("type") == "direct" and direction.get("type") == "directRcv":
            sender = contact_id_of(info.get("contact") or {})
            chat = "@" + str(sender)
        elif info.get("type") == "group" and direction.get("type") == "groupRcv":
            member = direction.get("groupMember") or {}
            sender = str(member.get("memberContactId") or "")
            group_id = (info.get("groupInfo") or {}).get("groupId")
            if group_id is None:
                continue
            chat = "#" + str(group_id)
        else:
            continue
        if sender:
            message = {"chat": chat, "sender": sender, "text": str(msg.get("text") or ""),
                             "id": str(meta["itemId"])}
            if msg.get('type') != 'text' and item.get('file'):
                message['file'] = item['file']
                message['media_kind'] = msg['type']
                message['await_description'] = (item['file'].get('fileProtocol') == 'xftp' and
                                                event.get('type') == 'newChatItems')
                message['complete'] = complete or (item['file'].get('fileStatus') or {}).get('type') == 'rcvComplete'
            elif complete:
                continue
            quote = (item.get("quotedItem") or {}).get("itemId")
            if quote is not None:
                message["quote_id"] = str(quote)
            messages.append(message)
    return messages


def text_of(event: dict) -> tuple[str, str] | None:
    messages = messages_of(event)
    return (messages[0]["sender"], messages[0]["text"]) if messages else None
