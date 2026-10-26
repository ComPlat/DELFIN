"""The messenger bridge: the agent service's narrow door for the gateway.

The gateway (``delfin.messenger``, user ``delfin-gw`` on the Pi) connects
over a Unix socket. A connection is accepted only when the kernel says the
peer's uid is the one configured (SO_PEERCRED, :func:`rundir.peer_uid`) --
not what a client claims, but what the process that connected runs as.

What the gateway may do is a fixed set of JSON-line commands; everything
that goes out runs through the service's
:class:`~delfin.agent.service.redact.Redactor` and is cut to a length limit,
so summaries mask credentials and events stay bounded. Explicit file downloads
are binary artifacts from the session outbox; their contents are not text-redacted. The bridge itself connects to the service's own socket
(``agent.sock``) as a window client (``client="delfin-tui"``-like name) --
a process that is not an agent child, so it may act; the supervisor's own
agent-detection and cage stay untouched. Scoped room commands use exact session IDs. A phone decision consumes a
short-lived code bound to the live request and goes through permission.reply;
the offered local choices are preserved, with a second confirmation for lasting
rules. Sessions must be explicitly enabled locally; unscoped mutations are refused.

Commands (one JSON object per line, both ways):
    {"cmd": "status"}                     -> the service's status (trimmed)
    {"cmd": "sessions"}                   -> the session roster (trimmed)
    {"cmd": "pending"}                    -> what needs the user (questions,
                                             approvals named but not granted)
    {"cmd": "answer", "id", "text"}       -> refused; use a scoped reply
    {"cmd": "prompt", "session", "text"}  -> send a task line to a session
                                             (exact enabled ID required)
    {"cmd": "stop"}                       -> refused; use scoped stop
    {"cmd": "create", "name"?}            -> start a session in its own folder
                                             under DELFIN_MESSENGER_WORKSPACES,
                                             permission profile "default",
                                             enabled for the messenger
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
from delfin.messenger import activation

#: Who may connect. On the Pi the gateway runs as the user ``delfin-gw``;
#: overridable here so the tests run under the test runner's uid.
PEER_UID_ENV = "DELFIN_MESSENGER_PEER_UID"

#: Where the socket lives. Deliberately under the service's own runtime
#: directory (0700), not in /tmp at large; the provisioning script chgrps
#: it to ``delfin-gw`` and opens it with the bridge's helper below.
SOCKET_NAME = "messenger.sock"

#: A command or an event line is never larger than this.
MAX_LINE_BYTES = 8_192

#: The folder root of sessions started from the phone. Unset: the host does
#: not allow it. The phone names a session; the folder is root/<name>.
WORKSPACES_ENV = "DELFIN_MESSENGER_WORKSPACES"

#: Phone-started sessions per window of CREATE_WINDOW_S.
CREATE_LIMIT = 3
CREATE_WINDOW_S = 600.0

#: Commands the gateway may send, and what each takes.
COMMANDS = {
    "session": {"session": "string", "action": "string", "text": "string?", "id": "string?", "choice": "string?", "persist": "boolean?"},
    "status": {},
    "sessions": {},
    "pending": {},
    "answer": {"id": "string", "text": "string"},
    "prompt": {"session": "string?", "text": "string"},
    "stop": {},
    "create": {"name": "string?"},
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
        uid = int(raw)
    except ValueError:
        return None
    # The gateway is a user of its own. The service's own uid would admit
    # every process of this user -- the agent's shell commands among them,
    # which could then answer the agent's own approvals as "the phone".
    if uid == os.getuid():
        print(f"[messenger] {PEER_UID_ENV} is this service's own uid; "
              "the bridge admits nobody. Run the gateway as its own user "
              "(scripts/pi/pi_messenger.sh).", file=sys.stderr)
        return None
    return uid


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
        self._issued: dict[str, tuple] = {}
        self._created: list[float] = []

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
            if kind.startswith("boolean") and not isinstance(value, bool):
                return f"{key} must be a boolean"
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
        params = params or {}
        sid = params.get("sessionId")
        if sid and activation.record(sid).get("endpoint"):
            return await asyncio.to_thread(activation.ask_dashboard, sid, method, params)
        reply = await asyncio.to_thread(self._client().request, method, params)
        if method == "sessions.list":
            for record in activation.records():
                if not record.get("endpoint"):
                    continue
                try:
                    status = await asyncio.to_thread(activation.ask_dashboard, record["session_id"], "session.status", {})
                    if not status.get("error"):
                        reply.setdefault("sessions", []).append(status)
                except (OSError, ValueError):
                    pass
        return reply

    async def _cmd_session(self, params: dict) -> dict:
        """Room commands require an exact, local session ID; never fall back."""
        sid = params["session"]
        if not activation.enabled(sid):
            return {"error": "SimpleX ist für diese Sitzung ausgeschaltet. In TUI/Dashboard mit /sync aktivieren."}
        generation = activation.record(sid).get("generation")
        listing = await self._request("sessions.list")
        if not activation.enabled(sid) or activation.record(sid).get("generation") != generation:
            return {"error": "SimpleX-Aktivierung wurde geändert. Anfrage nicht ausgeführt."}
        row = next((r for r in listing.get("sessions", [])
                    if str(r.get("sessionId") or r.get("id") or "") == sid), None)
        if row is None:
            return {"error": "Diese Sitzung ist nicht verfügbar. Kein anderer Agent wurde angesprochen."}
        action = params["action"]
        if action in {'files', 'announcement'} or action.startswith(('upload.', 'download.')):
            from delfin.messenger.media import Exchange
            if not hasattr(self, '_media'):
                self._media = Exchange()
            return self._media.command(sid, activation.record(sid).get('generation'),
                row.get('workspace'), action, token=params.get('id', ''), text=params.get('text', ''))
        if action == "status":
            return {"ok": True, "state": row.get("state"), "title": row.get("title"), "sessionId": sid,
                    "messengerEnabled": True, "generation": activation.record(sid).get("generation"),
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
                    details = req.get("choiceDetails") or []
                    choices = [str(choice) + (": " + str(details[i]) if i < len(details) and details[i] else "")
                               for i, choice in enumerate(req["choices"])]
                    preview += "\n" + self._redactor().text("\n".join(choices))
                fingerprint = json.dumps([activation.record(sid).get("generation"), req], sort_keys=True)
                code = next((k for k, v in self._issued.items() if v[:3] == (sid, rid, fingerprint)), None)
                if code is None:
                    code = secrets.token_hex(4)
                    self._issued[code] = (sid, rid, fingerprint, now)
                keys = {o.get("key") for o in req.get("options", [])}
                # Never approve a truncated action; lasting rules require confirmation.
                allowed = len(preview) <= 1400 and ((kind == "confirm" and "y" in keys) or
                                                   (kind == "plan" and bool(keys & {"d", "e"})))
                items.append({"id": code, "kind": "question" if kind == "ask" else kind,
                              "sessionId": sid, "text": _truncate(preview, 1400),
                              "phoneApproval": allowed, "tool": req.get("tool", ""),
                              "options": self._redactor().value(req.get("options", [])),
                              "choices": req.get("choices", [])})
            return {"ok": True, "items": items}
        if action in {"answer", "approve", "deny", "decide"}:
            code = str(params.get("id") or "")
            issued = self._issued.get(code)
            if issued is None or issued[0] != sid or time.monotonic() - issued[3] >= 600:
                return {"error": "Anfrage abgelaufen oder gehört zu einer anderen Sitzung. pending erneut aufrufen."}
            live = await self._request("session.pending", {"sessionId": sid})
            req = next((r for r in live.get("pending", []) if r.get("requestId") == issued[1]), None)
            if req is None or json.dumps([activation.record(sid).get("generation"), req], sort_keys=True) != issued[2]:
                self._issued.pop(code, None)
                return {"error": "Anfrage nicht mehr offen oder geändert."}
            kind = req.get("kind")
            reply_params = {"sessionId": sid, "requestId": issued[1]}
            if action == "answer" and kind == "ask":
                reply_params["text"] = str(params.get("text") or "")
                if not reply_params["text"].strip():
                    return {"error": "Leere Antwort."}
            elif action in {"approve", "deny", "decide"} and kind in {"confirm", "plan"}:
                choice = params.get("choice", "") if action == "decide" else "n" if action == "deny" else "d" if kind == "plan" else "y"
                if choice not in {o.get("key") for o in req.get("options", [])}:
                    return {"error": "Diese Entscheidung ist für die Anfrage nicht verfügbar."}
                preview = self._redactor().text(str(req.get("preview") or ""))
                if choice != "n" and len(preview) > 1400:
                    return {"error": "Aktion zu lang für die Handyansicht. Bitte am Dashboard freigeben."}
                if choice == "?":
                    return {"error": "Bitte eine der angezeigten Optionen auswählen."}
                # Lasting rules need an explicit second reply to a new message.
                if choice in {"A", "k", "p"} and not params.get("persist"):
                    return {"ok": False, "confirmationRequired": True, "id": code, "choice": choice}
                reply_params.update(choice=choice, persist=bool(params.get("persist", False)))
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
            if row.get('workspace'):
                from delfin.messenger.media import directory
                with directory(row['workspace'], sid, 'outbox') as (_, outbox):
                    text += ('\n\n[Messenger: Dateien/Bilder nur auf ausdrücklichen Wunsch des Nutzers '
                             'für den Download in ' + str(outbox) + ' bereitstellen. '
                             'Keine Geheimnisse oder Zugangsdaten exportieren. Den Dateinamen nennen; '
                             'der Nutzer ruft ihn mit datei NAME ab. Kein kompletter Chat-Export. '
                             'Für eine ausdrückliche kurze Mitteilung an andere Agenten kannst du forum.txt '
                             'dort ablegen (max. 1400 UTF-8-Bytes). Nur bei vom Nutzer eingeschaltetem '
                             'Austausch wird sie weitergegeben. Keine automatischen Antworten auf fremde Mitteilungen.]')
            reply = await self._request("session.prompt", {"sessionId": sid, "text": text, "how": "steer" if action == "update" else "queue"})
            return {"ok": True, "sessionId": sid, **reply}
        if action == "stop":
            reply = await self._request("session.cancel", {"sessionId": sid})
            return {"ok": True, **reply}
        if action == "rename":
            from .protocol import clean_title
            title = clean_title(params.get("text"))
            if not title:
                return {"error": "Bitte einen Namen angeben."}
            reply = await self._request("session.setTitle", {"sessionId": sid, "title": title})
            return {"ok": not reply.get("error"), **reply}
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
                "messengerEnabled": activation.enabled(str(row.get("sessionId") or row.get("id") or "")),
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

    async def _cmd_answer(self, _params: dict) -> dict:
        return {"error": "Use a session-scoped question or approval reply."}

    async def _cmd_prompt(self, params: dict) -> dict:
        sid = str(params.get('session') or '')
        if not sid:
            return {'error': 'An explicitly activated session ID is required.'}
        return await self._cmd_session({'session': sid, 'action': 'prompt', 'text': params['text']})

    async def _cmd_stop(self, _params: dict) -> dict:
        return {'error': 'Use session-scoped stop; the phone cannot shut down the entire service.'}

    async def _cmd_create(self, params: dict) -> dict:
        """A new session in root/<name>, asking before it acts.

        The phone supplies a name only. The folder is one level under the
        configured root (a name keeps [A-Za-z0-9._-], no leading dot, at
        most 40 characters; an existing entry is reused only when it is a
        real directory, never a symlink). The permission profile is
        "default" whatever the user's configured one is, and the session is
        enabled for the messenger, so its chat follows without /sync.
        """
        import re
        from .protocol import clean_title
        root = os.environ.get(WORKSPACES_ENV, "").strip()
        if not root or not os.path.isabs(root) or not os.path.isdir(root):
            return {"error": "Neue Sitzungen vom Handy sind auf diesem Rechner nicht eingerichtet."}
        now = time.monotonic()
        self._created = [t for t in self._created if now - t < CREATE_WINDOW_S]
        if len(self._created) >= CREATE_LIMIT:
            return {"error": "Höchstens drei neue Sitzungen in zehn Minuten. Bitte später erneut."}
        title = clean_title(params.get("name")) or ""
        slug = re.sub(r"[^A-Za-z0-9._-]+", "-", title).strip("-.")[:40].strip("-.")
        if not slug:
            slug = "handy-" + time.strftime("%Y%m%d-%H%M%S")
        folder = os.path.join(root, slug)
        if os.path.lexists(folder):
            if os.path.islink(folder) or not os.path.isdir(folder):
                return {"error": f"{slug} gibt es dort schon und ist kein normaler Ordner."}
        else:
            os.mkdir(folder, 0o700)
        self._created.append(now)
        reply = await self._request("session.create", {
            "cwd": folder, "title": title or slug, "permissionMode": "default",
            "attach": False})
        sid = str(reply.get("sessionId") or "")
        if not sid:
            return {"error": "Sitzung wurde nicht gestartet."}
        activation.configure(sid)
        return {"ok": True, "sessionId": sid, "title": title or slug, "folder": slug}


def main(argv: list[str] | None = None) -> int:
    """`python -m delfin.agent.service.bridge`: serve until interrupted."""
    asyncio.run(Bridge().serve())
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
