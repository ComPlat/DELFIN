"""The messenger bridge: the agent service's narrow door for the gateway.

The gateway (``delfin.messenger``, user ``delfin-gw`` on the Pi) connects
over a Unix socket. A connection is accepted only when the kernel says the
peer's uid is the one configured (SO_PEERCRED, :func:`rundir.peer_uid`) --
not what a client claims, but what the process that connected runs as.

What the gateway may do is a fixed set of JSON-line commands; everything
that goes out runs through the service's
:class:`~delfin.agent.service.redact.Redactor` and is cut to a length limit,
so a credential never reaches the phone and an event never becomes a
megabyte. The bridge itself connects to the service's own socket
(``agent.sock``) as a window client (``client="delfin-tui"``-like name) --
a process that is not an agent child, so it may act; the supervisor's own
agent-detection and cage stay untouched. Scoped room commands use exact session IDs. A phone decision consumes a
short-lived code bound to the live request and goes through permission.reply;
only one-shot confirmation or a plan's default mode is allowed. No permanent
allow rule or accept-edits mode can be selected from the messenger.

Commands (one JSON object per line, both ways):
    {"cmd": "status"}                     -> the service's status (trimmed)
    {"cmd": "sessions"}                   -> the session roster (trimmed)
    {"cmd": "pending"}                    -> what needs the user (questions,
                                             approvals named but not granted)
    {"cmd": "answer", "id", "text"}       -> answer a parked question
    {"cmd": "prompt", "session", "text"}  -> send a task line to a session
                                             (creates the phone session
                                             on first use)
    {"cmd": "stop"}                       -> stop the whole agent service
Outgoing events (one JSON object per line):
    {"event": "session", ...}   a session's summary changed
    {"event": "attention", ...} what needs the user changed
    {"event": "service", ...}   the service stops or changed
"""

from __future__ import annotations

import asyncio
import json
import os
import secrets
import time
import sys
from pathlib import Path

from . import rundir

#: Who may connect. On the Pi the gateway runs as the user ``delfin-gw``;
#: overridable here so the tests run under the test runner's uid.
PEER_UID_ENV = "DELFIN_MESSENGER_PEER_UID"

#: Where the socket lives. Deliberately under the service's own runtime
#: directory (0700), not in /tmp at large; the provisioning script chgrps
#: it to ``delfin-gw`` and opens it with the bridge's helper below.
SOCKET_NAME = "messenger.sock"

#: A command or an event line is never larger than this.
MAX_LINE_BYTES = 8_192

#: Commands the gateway may send, and what each takes.
COMMANDS = {
    "session": {"session": "string", "action": "string", "text": "string?", "id": "string?"},
    "status": {},
    "sessions": {},
    "pending": {},
    "answer": {"id": "string", "text": "string"},
    "prompt": {"session": "string?", "text": "string"},
    "stop": {},
}


class PeerRefused(RuntimeError):
    """A connection whose peer uid is not the configured one."""


def configured_peer_uid() -> int | None:
    """The uid allowed to connect, from the environment; None means the
    check cannot pass (nobody may connect) unless testing disabled it."""
    raw = os.environ.get(PEER_UID_ENV, "").strip()
    if not raw:
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def peer_uid_allowed(uid: int | None, allowed: int | None) -> bool:
    """The single check a connection must pass. Never inverted: an unknown
    peer uid (platform without the answer) refuses, too."""
    return uid is not None and allowed is not None and uid == allowed


def socket_path() -> Path:
    return Path(os.environ.get("DELFIN_MESSENGER_SOCKET") or (rundir.runtime_dir() / SOCKET_NAME))


def _truncate(text: str, limit: int = 400) -> str:
    text = text if isinstance(text, str) else str(text)
    if len(text) <= limit:
        return text
    return text[:limit - 1] + "…"


def _pending_count(value) -> int:
    return value if isinstance(value, int) else len(value or [])


def _short_event(event: dict, redactor) -> dict:
    """The one line about an event that goes to the gateway: kind, which
    session, the title and state, a question's text -- through the
    redactor, cut."""
    event = redactor.value(event)
    etype = str(event.get("type") or "")
    kind = ("attention" if etype in ("attention.updated",) else
            "service" if etype in ("service.stopping", "service.status")
            else "session")
    if kind == "attention":
        items = event.get("items") or []
        trimmed = []
        for item in items[:8]:
            trimmed.append({
                "id": str(item.get("id") or ""),
                "kind": str(item.get("kind") or ""),
                "sessionId": str(item.get("sessionId") or item.get("session") or ""),
                "sessionTitle": _truncate(str(item.get("sessionTitle") or "")),
                "text": redactor.text(_truncate(str(item.get("text") or ""))),
                "detail": redactor.text(_truncate(str(item.get("detail") or ""))),
            })
        return {"event": kind, "ts": event.get("ts"), "items": trimmed}
    if kind == "session":
        summary = event.get("summary") or {}
        return {"event": kind, "ts": event.get("ts"), "session": event.get("session"),
                "state": str(summary.get("state") or ""),
                "title": redactor.text(_truncate(str(summary.get("title") or ""))),
                "pending": _pending_count(summary.get("pending")),
                "answer": _truncate(str(summary.get("answer") or ""), 1500)}
    return {"event": kind, "ts": event.get("ts"),
            "reason": _truncate(str(event.get("reason") or ""))}


class Bridge:
    """Serves the gateway; talks to the agent service through one
    long-lived ServiceClient (its window name says what it is)."""

    CLIENT_NAME = "delfin-messenger-bridge"

    def __init__(self, *, allowed_uid: int | None = None,
                 service_path: Path | None = None,
                 redactor=None, queue_depth: int = 256):
        import socket as socket_mod

        self._socket = socket_mod
        self.allowed_uid = configured_peer_uid() if allowed_uid is None else allowed_uid
        self.service_path = service_path
        self.redactor = redactor
        self.queue_depth = queue_depth
        self._subscribers: dict[asyncio.StreamWriter, object] = {}
        self._service: object | None = None
        self._pump_task: object | None = None
        self._phone_session: str | None = None
        self._issued: dict[str, tuple] = {}

    # -- the service side ---------------------------------------------------------------

    def _connect_service(self):
        from .client import ServiceClient

        client = ServiceClient(self.service_path, client=self.CLIENT_NAME)
        # Everything the gateway may ask needs may_act: prompt creates and
        # sends, stop ends the service. This process is not an agent child.
        self._service = client
        return client

    def _redactor(self):
        if self.redactor is None:
            from .redact import Redactor
            self.redactor = Redactor()
        return self.redactor

    # -- the gateway side ----------------------------------------------------------------

    async def serve(self, path: Path | None = None) -> None:
        """Bind the socket, enforce its rights, serve until cancelled."""
        path = path or socket_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() or path.is_symlink():
            try:
                path.unlink()
            except OSError:
                pass
        server = await asyncio.start_unix_server(self._client_task, path=str(path))
        os.chmod(path, 0o660)
        # The provisioned unit runs with the gateway group, so the socket
        # inherits that group without opening the private service directory.
        self._socket_path = path

        self._loop = asyncio.get_running_loop()
        serving = None
        try:
            self._connect_service()
            self._pump_task = self._loop.run_in_executor(None, self._pump_service)
            serving = asyncio.create_task(server.serve_forever())
            done, _ = await asyncio.wait(
                [serving, self._pump_task], return_when=asyncio.FIRST_COMPLETED)
            for task in done:
                task.result()
            raise ConnectionError("agent service connection closed")
        finally:
            server.close()
            await server.wait_closed()
            if serving is not None:
                serving.cancel()
                await asyncio.gather(serving, return_exceptions=True)
            if self._service is not None:
                self._service.close()
            if self._pump_task is not None:
                await self._pump_task

    def _pump_service(self):
        """Thread: shorten and fan out what the service sends. The client
        collects events in ``client.events`` (its reader thread runs on its
        own); this walks the list on, so nothing is sent twice."""
        index = 0
        while True:
            client = self._service
            if client is None or client.closed:
                break
            with client._cond:
                if index >= len(client.events):
                    client._cond.wait(0.5)
                    continue
            while index < len(client.events):
                event = client.events[index]
                index += 1
                if event.get("type") not in {"session.summary", "attention.updated",
                                              "service.stopping", "service.status"}:
                    continue
                short = _short_event(event, self._redactor())
                self._loop.call_soon_threadsafe(self._fan_out, short)

    def _fan_out(self, short: dict) -> None:
        for writer, queue in list(self._subscribers.items()):
            try:
                queue.put_nowait(short)
            except Exception:
                pass

    async def _client_task(self, reader, writer) -> None:
        sock = writer.get_extra_info("socket")
        uid = rundir.peer_uid(sock) if sock is not None else None
        if not peer_uid_allowed(uid, self.allowed_uid):
            writer.close()
            return
        queue: asyncio.Queue = asyncio.Queue(maxsize=self.queue_depth)
        self._subscribers[writer] = queue
        async def send_events():
            while True:
                await self._send(writer, await queue.get())

        sender = asyncio.create_task(send_events())
        try:
            # Events flow even while the gateway sends no commands.
            while True:
                line = await asyncio.wait_for(reader.readline(), timeout=None)
                if not line:
                    break
                if len(line) > MAX_LINE_BYTES:
                    await self._send(writer, {"error": "line too long"})
                    continue
                try:
                    request = json.loads(line)
                except ValueError:
                    await self._send(writer, {"error": "not JSON"})
                    continue
                if not isinstance(request, dict) or "cmd" not in request:
                    await self._send(writer, {"error": "send {\"cmd\": ...}"})
                    continue
                reply = await self._run_command(request)
                await self._send(writer, reply)
        except (ConnectionError, asyncio.IncompleteReadError, OSError):
            pass
        finally:
            sender.cancel()
            await asyncio.gather(sender, return_exceptions=True)
            self._subscribers.pop(writer, None)
            writer.close()

    async def _send(self, writer, payload: dict) -> None:
        payload = self._redactor().value(payload)
        data = (json.dumps(payload, ensure_ascii=False) + "\n").encode()
        if len(data) > MAX_LINE_BYTES:
            payload = {"error": "output too long",
                       "fallback": str(payload.get("cmd") or payload.get("event") or "")}
            data = (json.dumps(payload) + "\n").encode()
        writer.write(data)
        await writer.drain()

    # -- commands ------------------------------------------------------------------------

    async def _run_command(self, request: dict) -> dict:
        cmd = str(request.get("cmd") or "")
        params = {k: v for k, v in request.items() if k != "cmd"}
        try:
            handler = getattr(self, "_cmd_" + cmd, None)
            if handler is None:
                return {"error": f"no command {cmd!r}",
                        "commands": sorted(COMMANDS)}
            problem = self._params_problem(COMMANDS.get(cmd, {}), params)
            if problem:
                return {"error": problem}
            return await handler(params)
        except Exception as exc:  # noqa: BLE001
            return {"error": f"{type(exc).__name__}: {exc}"}

    @staticmethod
    def _params_problem(spec: dict, params: dict) -> str | None:
        for key, kind in spec.items():
            optional = kind.endswith("?")
            value = params.get(key)
            if value is None and optional:
                continue
            if value is None:
                return f"missing {key}"
            if kind.startswith("string") and not isinstance(value, str):
                return f"{key} must be a string"
        extra = set(params) - {k.split("?")[0] for k in spec}
        if extra:
            return f"unknown fields: {sorted(extra)}"
        return None

    def _client(self):
        if self._service is None:
            self._connect_service()
        return self._service

    async def _request(self, method: str, params: dict | None = None) -> dict:
        return await asyncio.to_thread(self._client().request, method, params or {})

    async def _cmd_session(self, params: dict) -> dict:
        """Room commands require an exact, local session ID; never fall back."""
        sid = params["session"]
        listing = await self._request("sessions.list")
        row = next((r for r in listing.get("sessions", [])
                    if str(r.get("sessionId") or r.get("id") or "") == sid), None)
        if row is None:
            return {"error": "Diese Sitzung ist nicht verfügbar. Kein anderer Agent wurde angesprochen."}
        action = params["action"]
        if action == "status":
            return {"ok": True, "state": row.get("state"), "title": row.get("title"), "sessionId": sid,
                    "activity": _truncate(self._redactor().text(str(row.get("activity") or ""))),
                    "answer": _truncate(self._redactor().text(str(row.get("answer") or "")), 1000)}
        if action == "pending":
            live = await self._request("session.pending", {"sessionId": sid})
            items = []
            now = time.monotonic()
            self._issued = {k: v for k, v in self._issued.items() if now - v[3] < 600}
            for req in live.get("pending", []):
                rid, kind = str(req.get("requestId") or ""), req.get("kind")
                preview = self._redactor().text(str(req.get("question") or req.get("preview") or ""))
                if kind == "ask" and req.get("choices"):
                    preview += "\n" + self._redactor().text(" / ".join(map(str, req["choices"])))
                fingerprint = json.dumps(req, sort_keys=True)
                code = next((k for k, v in self._issued.items() if v[:3] == (sid, rid, fingerprint)), None)
                if code is None:
                    code = secrets.token_hex(4)
                    self._issued[code] = (sid, rid, fingerprint, now)
                keys = {o.get("key") for o in req.get("options", [])}
                # Never approve a truncated action or a lasting permission change.
                allowed = len(preview) <= 1400 and ((kind == "confirm" and "y" in keys) or
                                                   (kind == "plan" and "d" in keys))
                items.append({"id": code, "kind": "question" if kind == "ask" else kind,
                              "sessionId": sid, "text": _truncate(preview, 1400),
                              "phoneApproval": allowed, "tool": req.get("tool", "")})
            return {"ok": True, "items": items}
        if action in {"answer", "approve", "deny"}:
            code = str(params.get("id") or "")
            issued = self._issued.get(code)
            if issued is None or issued[0] != sid or time.monotonic() - issued[3] >= 600:
                return {"error": "Anfrage abgelaufen oder gehört zu einer anderen Sitzung. pending erneut aufrufen."}
            live = await self._request("session.pending", {"sessionId": sid})
            req = next((r for r in live.get("pending", []) if r.get("requestId") == issued[1]), None)
            if req is None or json.dumps(req, sort_keys=True) != issued[2]:
                self._issued.pop(code, None)
                return {"error": "Anfrage nicht mehr offen oder geändert."}
            kind = req.get("kind")
            reply_params = {"sessionId": sid, "requestId": issued[1]}
            if action == "answer" and kind == "ask":
                reply_params["text"] = str(params.get("text") or "")
                if not reply_params["text"].strip():
                    return {"error": "Leere Antwort."}
            elif action in {"approve", "deny"} and kind in {"confirm", "plan"}:
                choice = "n" if action == "deny" else "d" if kind == "plan" else "y"
                if choice not in {o.get("key") for o in req.get("options", [])}:
                    return {"error": "Diese Entscheidung ist für die Anfrage nicht verfügbar."}
                preview = self._redactor().text(str(req.get("preview") or ""))
                if action == "approve" and len(preview) > 1400:
                    return {"error": "Aktion zu lang für die Handyansicht. Bitte am Dashboard freigeben."}
                reply_params.update(choice=choice, persist=False)
            else:
                return {"error": "Antwort und Freigabe dürfen nicht verwechselt werden."}
            # Consume before awaiting: concurrent replies cannot apply this code twice.
            if time.monotonic() - issued[3] >= 600 or self._issued.pop(code, None) != issued:
                return {"error": "Anfrage bereits beantwortet."}
            reply = await self._request("permission.reply", reply_params)
            return {"ok": bool(reply.get("applied")), "error": None if reply.get("applied") else "Anfrage bereits erledigt."}
        if action in {"prompt", "update"}:
            text = str(params.get("text") or "").strip()
            if not text:
                return {"error": "Leerer Auftrag."}
            reply = await self._request("session.prompt", {"sessionId": sid, "text": text, "how": "steer" if action == "update" else "queue"})
            return {"ok": True, "sessionId": sid, **reply}
        if action == "stop":
            reply = await self._request("session.cancel", {"sessionId": sid})
            return {"ok": True, **reply}
        return {"error": "Unbekannter Sitzungsbefehl."}

    async def _cmd_status(self, _params: dict) -> dict:
        reply = await self._request("service.status")
        sessions = []
        for row in reply.get("sessions") or []:
            sessions.append({
                "id": str(row.get("sessionId") or row.get("id") or ""),
                "title": self._redactor().text(_truncate(str(row.get("title") or ""))),
                "state": str(row.get("state") or ""),
                "pending": _pending_count(row.get("pending")),
            })
        return {"ok": True, "host": reply.get("host"), "pid": reply.get("pid"),
                "sessions": sessions}

    async def _cmd_sessions(self, _params: dict) -> dict:
        reply = await self._request("sessions.list")
        rows = []
        for row in (reply.get("sessions") or []) + (reply.get("elsewhere") or []):
            rows.append({
                "id": str(row.get("sessionId") or row.get("id") or ""),
                "name": self._redactor().text(_truncate(str(row.get("title") or row.get("name") or ""))),
                "state": str(row.get("state") or ""),
                "elsewhere": bool(row.get("host")) if row in (reply.get("elsewhere") or []) else False,
            })
        return {"ok": True, "sessions": rows}

    async def _cmd_pending(self, _params: dict) -> dict:
        reply = await self._request("attention.list")
        items = []
        for item in reply.get("items") or []:
            kind = str(item.get("kind") or "")
            text = self._redactor().text(_truncate(str(item.get("detail") or item.get("text") or "")))
            items.append({
                "id": str(item.get("id") or ""),
                "kind": kind,
                "session": self._redactor().text(_truncate(str(item.get("sessionTitle") or ""))),
                "sessionId": str(item.get("sessionId") or item.get("session") or ""),
                "text": text,
                # An approval is NAMED here, never granted: the phone is
                # told where to approve, the dashboard keeps the right.
                "approvable": kind in ("confirm", "plan"),
            })
        return {"ok": True, "items": items}

    async def _cmd_answer(self, params: dict) -> dict:
        # Questions only. attention.answer refuses a confirm/plan itself;
        # we refuse it here already, so the refusal is not an accident.
        item_id = str(params.get("id") or "")
        text = str(params.get("text") or "")
        if not item_id or not text:
            return {"error": "answer needs id and text"}
        listing = await self._request("attention.list")
        kinds = {str(it.get("id") or ""): str(it.get("kind") or "")
                 for it in listing.get("items") or []}
        kind = kinds.get(item_id)
        if kind is None:
            return {"error": f"no pending item {item_id!r}"}
        if kind != "question":
            return {"error": "approvals are never answered from the phone; "
                             "approve at the dashboard or terminal"}
        reply = await self._request("attention.answer", {"id": item_id, "text": text})
        if not reply.get("ok"):
            return {"error": str(reply.get("error") or "the answer was not taken")}
        return {"ok": True}

    async def _cmd_prompt(self, params: dict) -> dict:
        text = str(params.get("text") or "").strip()
        if not text:
            return {"error": "a task needs some text"}
        session = str(params.get("session") or "").strip()
        if session and session != "phone":
            known = await self._known_session(session)
            if not known:
                return {"error": f"no session {session!r}"}
            session = known
        else:
            session = await self._ensure_phone_session()
        reply = await self._request("session.prompt",
                                    {"sessionId": session, "text": text,
                                     "how": "queue"})
        return {"ok": True, "session": session,
                "accepted": bool(reply.get("accepted")),
                "queued": bool(reply.get("queued"))}

    async def _known_session(self, name_or_id: str) -> str | None:
        try:
            reply = await self._request("sessions.list")
        except Exception:
            return None
        rows = (reply.get("sessions") or []) + (reply.get("elsewhere") or [])
        for row in rows:
            if name_or_id in (str(row.get("sessionId") or row.get("id") or ""), str(row.get("title") or "")):
                return str(row.get("sessionId") or row.get("id") or "")
        return None

    async def _ensure_phone_session(self) -> str:
        if self._phone_session:
            # Still there? A stopped service or an ended session would
            # make prompts vanish; the list is the truth.
            if await self._known_session(self._phone_session):
                return self._phone_session
            self._phone_session = None
        # Named session, created once, reused. The cwd is the service's
        # HOME -- the gateway picks none.
        reply = await self._request(
            "session.create",
            {"cwd": str(Path.home()), "title": "phone", "name": "phone"})
        self._phone_session = str(reply.get("sessionId") or "")
        return self._phone_session

    async def _cmd_stop(self, _params: dict) -> dict:
        reply = await self._request("service.shutdown")
        return {"ok": True, "stopping": bool(reply.get("stopping")),
                "sessions": int(reply.get("sessions") or 0)}


def main(argv: list[str] | None = None) -> int:
    """`python -m delfin.agent.service.bridge`: serve until interrupted."""
    asyncio.run(Bridge().serve())
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
