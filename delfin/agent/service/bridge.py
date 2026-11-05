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
    {"cmd": "schedules", "after"?}        -> the user's schedules and the runs
                                            recorded after ``after``
    {"cmd": "approvals"}                  -> every pending approval of the
                                            service's sessions, numbered,
                                            with the check to answer with
    {"cmd": "answer", "n"?, "ok"?, "check"?, "id"?, "text"?}
                                          -> the numbered approval (n, ok,
                                             check) or a session question
                                             (id, text; scope unchanged)
    {"cmd": "schedule", "action", "id"?, "text"?, "session"?, "check"?}
                                          -> new (a proposal), ok (with the
                                             check it was shown), show, pause,
                                             resume, delete, run
    {"cmd": "attach", "session", "label"?} -> people mode's ``hol``: switch
                                            a running session of the person
                                            on for the messenger (not one
                                            switched off in the session)
    {"cmd": "pairing", "action": poll|answer|note, ...}
                                          -> the persons' pairing requests
                                             (delfin.messenger.pairing) and
                                             the gateway's answers
Every command may carry ``"person"`` (people mode, delfin.messenger.people):
it then runs in that person's own context -- their agent service socket
(``~/people/<name>/.run``), their activation records, their phone folder,
their issued approval codes -- and never reaches another person's
sessions. Same Unix user: an application-level boundary, not an OS one.
Outgoing events (one JSON object per line):
    {"event": "session", ...}   a session's summary changed
    {"event": "attention", ...} what needs the user changed
    {"event": "service", ...}   the service stops or changed
"""

from __future__ import annotations

import asyncio
import hashlib
import hmac
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
    "session": {"session": "string", "action": "string", "text": "string?", "id": "string?", "choice": "string?", "persist": "boolean?", "command": "string?"},
    "status": {},
    "sessions": {},
    "pending": {},
    "approvals": {},
    "answer": {"n": "integer?", "ok": "boolean?", "check": "string?",
               "id": "string?", "text": "string?"},
    "prompt": {"session": "string?", "text": "string"},
    "stop": {},
    "create": {"name": "string?"},
    "schedules": {"after": "string?"},
    "schedule": {"action": "string", "id": "string?", "text": "string?",
                 "session": "string?", "check": "string?"},
    "profile": {"session": "string?", "workspace": "string?",
                "text": "string"},
    # People mode: bring a running session of the person to the phone
    # (``hol``). Refused for a session switched off in the session itself.
    "attach": {"session": "string", "label": "string?"},
    # The pairing requests of every person (delfin.messenger.pairing):
    # poll -> the open requests; answer -> the gateway's answer to one;
    # note -> a line for the person (a new phone, a revocation).
    "pairing": {"action": "string", "person": "string?", "id": "string?",
                "result": "object?", "text": "string?"},
}

#: Every command may carry ``person`` (a person of people mode, see
#: delfin.messenger.people; absent or "" is the account owner). The command
#: then runs in that person's own context: their agent service, their
#: activation records, their phone folder -- never another person's.
PERSON_FIELD = "person"

#: Commands a person's context refuses: schedules and freedom profiles stay
#: the account owner's (they are one store per account today).
OWNER_ONLY = frozenset({"schedules", "schedule", "profile"})

#: Schedules made from the phone per window of CREATE_WINDOW_S.
SCHEDULE_LIMIT = 5

#: The folder (under DELFIN_MESSENGER_WORKSPACES) a schedule made in the
#: phone's schedules chat runs in.
SCHEDULE_FOLDER = "zeitplaene"

#: Why a schedule was refused, in the phone's words (schedules.ScheduleError
#: codes); the phone translates these whole.
SCHEDULE_ERRORS = {
    "when": "Zeit nicht verstanden. Beispiele: täglich 07:30 · mo-fr 18:00 · alle 2h · 30 7 * * 1-5",
    "time": "Uhrzeit fehlt oder ist ungültig. Beispiele: täglich 07:30 · mo-fr 18:00 · alle 2h · 30 7 * * 1-5",
    "interval": "Der kürzeste Abstand ist 5 Minuten.",
    "zone": "Unbekannte Zeitzone, z. B. tz=Europe/Berlin.",
    "never": "Dieser Zeitplan würde nie laufen.",
    "task": "Nach der Zeit fehlt die Aufgabe.",
    "workspace": "Der Ordner fehlt.",
    "unknown": "Zeitplan nicht gefunden. pläne zeigt die aktuellen.",
    "budget": "Ungültige Grenzen.",
    "command": "Ein Watcher führt nur lesende Befehle aus: squeue, sacct, sinfo, ls, tail, head, wc, grep, du, df und ähnliche.",
    "path": "Pfade müssen im Ordner des Zeitplans liegen.",
}


class PeerRefused(RuntimeError):
    """A connection whose peer uid is not the configured one."""


class ServiceDown(RuntimeError):
    """A person's agent service is not running; the text is for the phone."""


import contextvars

#: Set for one gateway request that must not start a person's service.
_NO_WAKE: contextvars.ContextVar = contextvars.ContextVar("bridge_no_wake", default=False)


class ServiceStarting(ServiceDown):
    """A person's agent service is being started (in a thread of its own);
    the command is answered at once and asked again once it is up."""


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

    #: The owner's context unless __init__ says otherwise (a bridge built
    #: without __init__ in a test is the owner's).
    person = ""
    activation_dir = None
    runtime_dir = None

    def __init__(self, *, allowed_uid: int | None = None,
                 service_path: Path | None = None,
                 redactor=None, queue_depth: int = 256, person: str = "",
                 activation_dir: Path | None = None, runtime_dir: Path | None = None):
        import socket as socket_mod

        self._socket = socket_mod
        self.allowed_uid = configured_peer_uid() if allowed_uid is None else allowed_uid
        self.service_path = service_path
        #: The person this bridge acts for ("" the account owner), their
        #: activation records (None: this HOME's) and the runtime directory
        #: a dashboard endpoint of theirs must lie in (None: unchecked).
        self.person = person
        self.activation_dir = activation_dir
        self.runtime_dir = runtime_dir
        self._people: dict[str, "Bridge"] = {}
        self.redactor = redactor
        self.queue_depth = queue_depth
        self._subscribers: dict[asyncio.StreamWriter, object] = {}
        self._service: object | None = None
        self._pump_task: object | None = None
        self._issued: dict[str, tuple] = {}
        self._approvals_numbered: dict[tuple, int] = {}
        self._approvals_issued: dict[tuple, float] = {}
        self._approvals_raw: dict[tuple, dict] = {}
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
        # The provisioned unit makes the socket's directory setgid delfin-gw
        # (pi_messenger.sh), so the socket takes that group while this
        # process keeps the service user's own: the service reads a window's
        # environment only with uid and gid equal, and lets it act only then.
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
        # wake: false -- the gateway's background look (notifications): it
        # never starts a person's stopped service; a phone's command does.
        token = _NO_WAKE.set(request.get("wake") is False)
        try:
            return await self._run_one({k: v for k, v in request.items() if k != "wake"})
        finally:
            _NO_WAKE.reset(token)

    async def _run_one(self, request: dict) -> dict:
        cmd = str(request.get("cmd") or "")
        person = request.get(PERSON_FIELD)
        if cmd != "pairing" and PERSON_FIELD in request:
            request = {k: v for k, v in request.items() if k != PERSON_FIELD}
            if person not in (None, ""):
                # Another person's context: their service and records only.
                if self.person:
                    return {"error": "Person bereits gewählt."}
                sub = self._person_bridge(person)
                if sub is None:
                    return {"error": "Unbekannte Person. Nichts ausgeführt."}
                return await sub._run_one(request)
        if self.person and cmd in OWNER_ONLY | {"pairing"}:
            return {"error": "Für diese Person vom Handy nicht verfügbar."}
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
        except ServiceStarting as exc:
            return {"error": str(exc), "starting": True}
        except ServiceDown as exc:
            return {"error": str(exc)}
        except Exception as exc:  # noqa: BLE001
            if self.person:
                # A person's paths stay off the phone.
                return {"error": f"Nicht ausgeführt ({type(exc).__name__})."}
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
            if kind.startswith("integer") and not isinstance(value, int):
                return f"{key} must be an integer"
            if kind.startswith("string") and not isinstance(value, str):
                return f"{key} must be a string"
            if kind.startswith("object") and not isinstance(value, dict):
                return f"{key} must be an object"
        extra = set(params) - {k.split("?")[0] for k in spec}
        if extra:
            return f"unknown fields: {sorted(extra)}"
        return None

    def _client(self):
        if self._service is None or (self.person and getattr(self._service, "closed", False)):
            try:
                self._connect_service()
            except OSError:
                if not self.person:
                    raise
                # A person's service exits when idle: started for the phone
                # the way their windows start it (people.request_start), in
                # a thread of its own -- never on the bridge's event loop.
                if _NO_WAKE.get():
                    raise ServiceDown("Der Agent-Dienst dieser Person läuft gerade nicht.") from None
                if not self._kick_person_start():
                    raise ServiceDown("Der Agent-Dienst dieser Person läuft gerade nicht. "
                                      "Am Rechner eine Sitzung öffnen.") from None
                raise ServiceStarting("Der Agent-Dienst dieser Person wird gestartet; "
                                      "die Antwort folgt gleich.") from None
        return self._service

    def _kick_person_start(self) -> bool:
        """A start of this person's service running in a thread: True when
        one runs (begun now or before), False when it is refused."""
        import threading
        running = self.__dict__.get("_start_thread")
        if running is not None and running.is_alive():
            return True
        if not self._may_start_person_service():
            return False
        thread = threading.Thread(target=self._start_person_service, daemon=True,
                                  name=f"person-start-{self.person}")
        self._start_thread = thread
        thread.start()
        return True

    #: Starts of one person's service per window, and the window.
    PERSON_STARTS = 1
    PERSON_START_WINDOW_S = 60.0

    def _may_start_person_service(self) -> bool:
        """The rate limit (PERSON_STARTS per PERSON_START_WINDOW_S); a
        refused start is audited."""
        from delfin.messenger import audit
        now = time.monotonic()
        starts = [t for t in self.__dict__.get("_person_starts", []) if now - t < self.PERSON_START_WINDOW_S]
        if len(starts) >= self.PERSON_STARTS:
            self._person_starts = starts
            audit.event("messenger_service_start_refused",
                        f"agent service of {self.person}: started less than "
                        f"{int(self.PERSON_START_WINDOW_S)} s ago", blocked=True)
            return False
        self._person_starts = [*starts, now]
        return True

    def _start_person_service(self) -> bool:
        """(In its own thread.) Start this person's agent service and
        connect to it; audited. True once connected."""
        from delfin.messenger import audit, people
        if not people.request_start(self.person):
            audit.event("messenger_service_start_refused",
                        f"agent service of {self.person}: start refused or failed", blocked=True)
            return False
        audit.event("messenger_service_started",
                    f"agent service of {self.person} started for the phone")
        deadline = time.monotonic() + people.START_WAIT_S
        while True:
            try:
                self._connect_service()
                return True
            except OSError:
                if time.monotonic() >= deadline:
                    return False
                time.sleep(0.25)

    # -- this context's activation records (the owner's: the module's own
    #    defaults, so a test that replaces a function keeps working) -----------

    def _where(self) -> dict:
        return {} if self.activation_dir is None else {"directory": self.activation_dir}

    def _enabled(self, sid) -> bool:
        return activation.enabled(sid, **self._where())

    def _record(self, sid) -> dict:
        return activation.record(sid, **self._where())

    def _records(self) -> list:
        return activation.records(**self._where())

    def _configure(self, sid, args="", **kwargs) -> str:
        return activation.configure(sid, args, **self._where(), **kwargs)

    def _switched_off(self, sid) -> bool:
        return activation.switched_off_in_session(sid, **self._where())

    def _ask_dashboard(self, sid: str, method: str, params: dict) -> dict:
        if self.activation_dir is None:
            return activation.ask_dashboard(sid, method, params)
        return activation.ask_dashboard(sid, method, params, directory=self.activation_dir,
                                        within=self.runtime_dir)

    def _person_bridge(self, name: str) -> "Bridge | None":
        """The context of person ``name`` (people mode), or None when the
        name is not a person of this account. One per person, kept: its
        issued approval codes and numbers are its own."""
        from delfin.messenger import people
        if not people.valid(name) or name == people.OWNER or people.home(name) is None:
            return None
        if "_people" not in self.__dict__:
            self._people = {}
        sub = self._people.get(name)
        if sub is None:
            sub = Bridge(allowed_uid=self.allowed_uid, service_path=people.runtime_socket(name),
                         redactor=self._redactor(), person=name,
                         activation_dir=people.activation_dir(name),
                         runtime_dir=people.runtime_dir(name))
            self._people[name] = sub
        return sub

    def _workspaces_root(self) -> str:
        """Where ``neu NAME`` makes folders for this context ("" : not allowed)."""
        from delfin.messenger import people
        owner_root = os.environ.get(WORKSPACES_ENV, "").strip()
        if not owner_root or not os.path.isabs(owner_root) or not os.path.isdir(owner_root):
            return ""
        if not self.person:
            return owner_root
        root = people.workspaces_root(self.person, owner_root)
        if root is None:
            return ""
        if not os.path.lexists(root):
            os.mkdir(root, 0o700)
        if os.path.islink(root) or not os.path.isdir(root):
            return ""
        return str(root)

    async def _request(self, method: str, params: dict | None = None,
                       timeout: float | None = None) -> dict:
        params = params or {}
        sid = params.get("sessionId")
        if sid and self._record(sid).get("endpoint"):
            return await asyncio.to_thread(self._ask_dashboard, sid, method, params)
        client = self._client()
        reply = await asyncio.to_thread(
            client.request, method, params, **({"timeout": timeout}
                                               if timeout is not None else {}))
        if method == "sessions.list":
            for record in self._records():
                if not record.get("endpoint"):
                    continue
                try:
                    status = await asyncio.to_thread(self._ask_dashboard, record["session_id"], "session.status", {})
                    if not status.get("error"):
                        reply.setdefault("sessions", []).append(status)
                except (OSError, ValueError):
                    pass
        return reply

    async def _cmd_session(self, params: dict) -> dict:
        """Room commands require an exact, local session ID; never fall back."""
        sid = params["session"]
        if not self._enabled(sid):
            return {"error": "SimpleX ist für diese Sitzung ausgeschaltet. In TUI/Dashboard mit /sync aktivieren."}
        generation = self._record(sid).get("generation")
        listing = await self._request("sessions.list")
        if not self._enabled(sid) or self._record(sid).get("generation") != generation:
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
            return self._media.command(sid, self._record(sid).get('generation'),
                row.get('workspace'), action, token=params.get('id', ''), text=params.get('text', ''))
        if action == "status":
            return {"ok": True, "state": row.get("state"), "title": row.get("title"), "sessionId": sid,
                    # The folder's name only, for naming an untitled session's chat.
                    "folder": _truncate(self._redactor().text(Path(str(row.get("workspace") or "")).name), 60),
                    "messengerEnabled": True, "generation": self._record(sid).get("generation"),
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
                fingerprint = json.dumps([self._record(sid).get("generation"), req], sort_keys=True)
                code = next((k for k, v in self._issued.items() if v[:3] == (sid, rid, fingerprint)), None)
                if code is None:
                    code = secrets.token_hex(4)
                    self._issued[code] = (sid, rid, fingerprint, now)
                keys = {o.get("key") for o in req.get("options", [])}
                # Never approve a truncated action; lasting rules require confirmation.
                allowed = len(preview) <= 1400 and ((kind == "confirm" and "y" in keys) or
                                                   (kind == "plan" and bool(keys & {"d", "e"})))
                flags = req.get("flags") or {}
                what = ("shell" if flags.get("shell") else "write" if flags.get("write")
                        else "read" if flags.get("outsideRead") else "")
                items.append({"id": code, "kind": "question" if kind == "ask" else kind,
                              "sessionId": sid, "text": _truncate(preview, 1400),
                              "action": what,
                              "purpose": _truncate(self._redactor().text(str(req.get("purpose") or "")), 200),
                              "phoneApproval": allowed, "tool": req.get("tool", ""),
                              "options": self._redactor().value(req.get("options", [])),
                              "choices": req.get("choices", [])})
            return {"ok": True, "items": items}
        if action == "branch":
            # Branch operations in this session's workspace: the dashboard's
            # rules (branch_ops), from the phone. A deletion rides
            # branch.cleanup with confirm=true only; branch.test can run
            # minutes, so it asks the service directly, not through a turn.
            ws = str(row.get("workspace") or "")
            if not ws:
                return {"error": "Diese Sitzung hat keinen Arbeitsordner."}
            sub = str(params.get("sub") or "")
            request = {"workspace": ws}
            if sub in {"list", "pull", "test", "cleanup"}:
                pass
            else:
                return {"error": "Unbekannter Branch-Befehl."}
            if sub == "list":
                request.update(view=str(params.get("view") or "mine"),
                               text=str(params.get("text") or ""))
                reply = await self._request("branch.list", request)
            elif sub == "pull":
                reply = await self._request("branch.pull", request)
            elif sub == "test":
                request.update(ref=str(params.get("ref") or ""))
                reply = await self._request("branch.test", request,
                                            timeout=1200.0)
            else:
                request.update(delete=params.get("delete") or [],
                               confirm=bool(params.get("confirm", False)))
                reply = await self._request("branch.cleanup", request,
                                            timeout=300.0)
            if reply.get("error"):
                return {"error": str(reply["error"])}
            return {"ok": bool(reply.get("ok")), "text": str(reply.get("text") or ""),
                    **({"rows": reply["rows"]} if "rows" in reply else {}),
                    **({"candidates": reply["candidates"]}
                       if "candidates" in reply else {})}
        if action in {"answer", "approve", "deny", "decide"}:
            code = str(params.get("id") or "")
            issued = self._issued.get(code)
            if issued is None or issued[0] != sid or time.monotonic() - issued[3] >= 600:
                return {"error": "Anfrage abgelaufen oder gehört zu einer anderen Sitzung. pending erneut aufrufen."}
            live = await self._request("session.pending", {"sessionId": sid})
            req = next((r for r in live.get("pending", []) if r.get("requestId") == issued[1]), None)
            if req is None or json.dumps([self._record(sid).get("generation"), req], sort_keys=True) != issued[2]:
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
                # A changed command goes with choice c only, and c only with one.
                command = str(params.get("command") or "")
                if (choice == "c") != bool(command.strip()):
                    return {"error": "Befehl ändern: mit „ändern: NEUER BEFEHL“ antworten."}
                reply_params.update(choice=choice, persist=bool(params.get("persist", False)),
                                    **({"command": command} if choice == "c" else {}))
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
            # The one display name: the service's session title, 50 characters.
            title = clean_title(params.get("text"))[:activation.NAME_MAX].strip()
            if not title:
                return {"error": "Bitte einen Namen angeben."}
            reply = await self._request("session.setTitle", {"sessionId": sid, "title": title})
            return {"ok": not reply.get("error"), **reply}
        if action == "detach":
            # weg / off in the session's chat: off for the phone, recorded as
            # the phone's switch, told in the session.
            self._configure(sid, "off", by="phone", phone=str(params.get("text") or ""))
            await self._tell_session(sid)
            return {"ok": True, "sessionId": sid}
        return {"error": "Unbekannter Sitzungsbefehl."}

    async def _tell_session(self, sid: str) -> None:
        """Show the session its messenger state as a notice (``/sync status``
        through the session's own handler). Best effort: a session that
        cannot be told is still switched."""
        try:
            await self._request("session.prompt", {"sessionId": sid, "text": "/sync status"})
        except Exception:  # noqa: BLE001
            pass

    async def _cmd_attach(self, params: dict) -> dict:
        """``hol`` on the phone: switch a running session of this person on
        for the messenger. Its exact id, from this person's own service; a
        session switched off in the session (``/sync off``) stays off."""
        sid = str(params.get("session") or "")
        label = str(params.get("label") or "")[:40]
        if self._switched_off(sid):
            return {"error": "In der Sitzung mit /sync off ausgeschaltet. Dort mit /sync wieder einschalten."}
        listing = await self._request("sessions.list")
        row = next((r for r in listing.get("sessions", [])
                    if str(r.get("sessionId") or r.get("id") or "") == sid), None)
        if row is None or not sid:
            return {"error": "Diese Sitzung ist nicht verfügbar. Kein anderer Agent wurde angesprochen."}
        already = self._enabled(sid)
        if not already:
            self._configure(sid, "on", by="phone", phone=label)
            await self._tell_session(sid)
        name = row.get("title") or ""
        return {"ok": True, "sessionId": sid, "already": already,
                "title": self._redactor().text(_truncate(str(name), 80))}

    async def _cmd_pairing(self, params: dict) -> dict:
        """The pairing requests of every person, and the gateway's answers
        (delfin.messenger.pairing). The person is the directory a request
        was found in, never a field the gateway sends unchecked."""
        from delfin.messenger import pairing, people
        action = str(params.get("action") or "")
        if action == "poll":
            found = []
            for name in people.names():
                path = people.pairing_dir(name)
                if path is None or not path.is_dir():
                    continue
                for item in pairing.collect(path):
                    found.append({**item, "person": name})
            return {"ok": True, "requests": found[:20]}
        name = str(params.get("person") or "")
        path = people.pairing_dir(name) if people.valid(name) else None
        if path is None:
            return {"error": "Unbekannte Person."}
        if action == "answer":
            return pairing.answer(path, str(params.get("id") or ""), dict(params.get("result") or {}),
                                  door_root=people.door_root(name))
        if action == "note":
            pairing.note(path, str(params.get("text") or ""))
            return {"ok": True}
        return {"error": "Unbekannte Pairing-Aktion."}

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
            sid = str(row.get("sessionId") or row.get("id") or "")
            rows.append({
                "id": sid,
                "name": self._redactor().text(_truncate(str(row.get("title") or row.get("name") or ""))),
                "folder": self._redactor().text(_truncate(Path(str(row.get("workspace") or "")).name, 40)),
                "switchedOff": self._switched_off(sid),
                "state": str(row.get("state") or ""),
                "messengerEnabled": self._enabled(sid),
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

    async def _cmd_approvals(self, _params: dict) -> dict:
        """Every pending approval of every session of this service --
        interactive, scheduled, unattended -- numbered, one short entry per
        request, with the HMAC check to answer with. Numbers stay stable
        while a request waits; a request that is gone frees its number."""
        listing = await self._request("sessions.list")
        rows = listing.get("sessions") or []
        items = []
        seen: dict[str, tuple] = {}
        for row in rows:
            sid = str(row.get("sessionId") or row.get("id") or "")
            title = self._redactor().text(_truncate(str(row.get("title") or "")))
            live = await self._request("session.pending", {"sessionId": sid})
            for req in live.get("pending", []):
                if str(req.get("kind")) not in {"confirm", "plan"}:
                    continue
                rid = str(req.get("requestId") or "")
                flags = req.get("flags") or {}
                what = ("shell" if flags.get("shell")
                        else "write" if flags.get("write")
                        else "read" if flags.get("outsideRead") else "")
                preview = self._redactor().text(str(req.get("preview") or ""))
                options = {o.get("key") for o in req.get("options", [])}
                # Always-explicit stays explicit: a request the approval
                # screen marks protected, or the existing outward
                # classifier (auto_guard) marks outward -- git push, gh pr
                # create/merge, forge API writes, mail, uploads -- is shown
                # but never phone-approved.
                outward = self._outward(req)
                approvable = ("y" in options and len(preview) <= 1400
                              and not flags.get("protected") and not outward)
                items.append({
                    "tool": str(req.get("tool") or ""),
                    "text": _truncate(preview, 300),
                    "risk": what,
                    "session": title,
                    "sessionId": sid,
                    "requestId": rid,
                    "approvable": approvable,
                    "_req": req,
                })
        if not items:
            return {"ok": True, "items": [], "check": ""}
        # The numbers: stable while a request waits, freed when it is gone.
        order = [(i["sessionId"], i["requestId"]) for i in items]
        now = time.monotonic()
        self._approvals_numbered = {
            key: num for key, num in self._approvals_numbered.items()
            if key in set(order) and now - self._approvals_issued.get(key, now) < 600}
        start = (max(self._approvals_numbered.values()) + 1
                 if self._approvals_numbered else 1)
        for i, key in enumerate(order, start):
            if key not in self._approvals_numbered:
                self._approvals_numbered[key] = i
                self._approvals_issued[key] = now
        for item in items:
            key = (item["sessionId"], item["requestId"])
            item["n"] = self._approvals_numbered[key]
            # The unchanged request, kept only inside this process: the
            # answer must apply to the same live request that was shown.
            self._approvals_raw[key] = item.pop("_req")
        return {"ok": True, "items": items, "check": self._approvals_check(items)}

    def _outward(self, req: dict) -> bool:
        """True when the existing classifier (auto_guard) marks this
        request outward: it publishes or sends something. One classifier,
        not a second list: the same rules the agent's own guard applies
        to git push, gh pr create/merge, forge API writes, mail, uploads
        and the outward MCP tools. Never raises (an error refuses)."""
        from delfin.agent import auto_guard
        name = str(req.get("tool") or "")
        text = str(req.get("preview") or "")
        # The preview holds the command as the screen shows it ("$ cmd");
        # auto_guard matches shell commands either way.
        verdict = auto_guard.classify_call(name, {"command": text})
        return bool(verdict and verdict.kind in {"outward", "protected"})

    def _approvals_check(self, items: list) -> str:
        """The code an answer must carry: an HMAC over the whole shown
        listing, so a message from another chat, a stale listing or a
        guessed number cannot approve. With no key (tests, headless) a
        plain fingerprint of the same payload -- never a constant."""
        from delfin.agent import trust_store
        payload = json.dumps(
            [[i["n"], i["sessionId"], i["requestId"], i["text"]] for i in items],
            sort_keys=True, ensure_ascii=False).encode()
        tag = trust_store.mac("user_approvals", payload)
        return tag if tag else "sha:" + hashlib.sha256(payload).hexdigest()

    def _approvals_check_ok(self, items: list, tag: str) -> bool:
        from delfin.agent import trust_store
        payload = json.dumps(
            [[i["n"], i["sessionId"], i["requestId"], i["text"]] for i in items],
            sort_keys=True, ensure_ascii=False).encode()
        good = trust_store.check_mac("user_approvals", payload, tag)
        if good is None:  # no key: only the fingerprint form can match
            good = hmac.compare_digest(
                "sha:" + hashlib.sha256(payload).hexdigest(), str(tag or ""))
        return bool(good)

    async def _cmd_answer(self, params: dict) -> dict:
        """A numbered approval (n, ok, check) or a session-scoped answer
        (id, text). The numbered path answers every session of the service,
        so it is bound to the HMAC check the listing showed."""
        n, ok, check = params.get("n"), params.get("ok"), params.get("check")
        if n is not None or check is not None or ok is not None:
            if not isinstance(n, int) or not isinstance(ok, bool) or not isinstance(check, str) or not check:
                return {"error": "Antwortform: ok <nummer> / nein <nummer>, mit der beigelegten Prüfung der Anfrage."}
            approvals = await self._cmd_approvals({})
            items = approvals.get("items", [])
            item = next((i for i in items if i["n"] == n), None)
            if item is None:
                return {"error": f"Keine offene Anfrage mit der Nummer {n}. Mit approvals die offenen ansehen."}
            if not self._approvals_check_ok(items, check):
                return {"error": "Prüfung stimmt nicht oder die Anfrage hat sich geändert. Erst approvals, dann antworten."}
            if not item["approvable"]:
                return {"error": "Diese Anfrage ist immer ausdrücklich: bitte am Dashboard freigeben."}
            key = (item["sessionId"], item["requestId"])
            shown = self._approvals_raw.get(key)
            live = await self._request("session.pending", {"sessionId": item["sessionId"]})
            req = next((r for r in live.get("pending", [])
                        if r.get("requestId") == item["requestId"]), None)
            if req is None or shown is None or json.dumps(req, sort_keys=True) != json.dumps(shown, sort_keys=True):
                return {"error": "Anfrage nicht mehr offen oder geändert. Mit approvals die offenen ansehen."}
            # Consume before awaiting: concurrent replies cannot approve twice;
            # the service applies a requestId once regardless.
            self._approvals_raw.pop(key, None)
            reply = await self._request("permission.reply", {
                "sessionId": item["sessionId"], "requestId": item["requestId"],
                "choice": "y" if ok else "n"})
            return {"ok": bool(reply.get("applied")),
                    "error": None if reply.get("applied") else "Anfrage bereits erledigt."}
        code = str(params.get("id") or "")
        text = str(params.get("text") or "")
        if not code or not text:
            return {"error": "Use a session-scoped question or approval reply."}
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
        root = self._workspaces_root()
        if not root:
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
        self._configure(sid)
        return {"ok": True, "sessionId": sid, "title": title or slug, "folder": slug}


    # -- the user's schedules ---------------------------------------------------------

    def _schedule_view(self, e: dict) -> dict:
        from delfin.agent import schedules
        red = self._redactor()
        if e.get("watch_done"):
            state = "done"
        elif e.get("host") != rundir.hostname():
            state = "elsewhere"
        elif not e["confirmed"]:
            state = "unconfirmed"
        elif e.get("paused"):
            state = "paused"
        elif e.get("running"):
            state = "running"
        else:
            state = "active"
        return {"id": e["id"], "when": _truncate(str(e["when"]), 80), "tz": str(e["tz"]),
                "next": schedules.describe_next(e),
                "task": red.text(_truncate(str(e["task"]), 300)),
                "folder": red.text(_truncate(Path(str(e["workspace"])).name, 60)),
                "host": _truncate(str(e.get("host") or ""), 60), "state": state,
                "lastStatus": str(e.get("last_status") or ""),
                "confirmed": bool(e["confirmed"]), "check": schedules.fingerprint(e),
                "runUsd": e.get("budget_run_usd"), "dayUsd": e.get("budget_day_usd"),
                "minutes": e.get("max_minutes")}

    async def _cmd_profile(self, params: dict) -> dict:
        """Set a workspace's freedom profile from the phone.

        ``profil <session|ordner> frei|standard`` -- a session's own
        workspace (exact, enabled ID, like every room command) or a
        named folder. ``frei`` is refused with the reason where the host
        cannot hold a shell command in the cage, and the word is
        validated before anything is stored. The store is the user's,
        masked like the schedules; the phone is the user.
        """
        from delfin.agent import workspace_profiles
        word = str(params.get("text") or "").strip().lower()
        if word not in {"frei", "free", "standard", "default"}:
            return {"error": "profil <sitzung|ordner> frei|standard"}
        sid = str(params.get("session") or "")
        workspace = str(params.get("workspace") or "")
        if not workspace and sid:
            if not self._enabled(sid):
                return {"error": "SimpleX ist für diese Sitzung "
                                 "ausgeschaltet. In TUI/Dashboard mit "
                                 "/sync aktivieren."}
            listing = await self._request("sessions.list")
            row = next((r for r in listing.get("sessions", [])
                        if str(r.get("sessionId") or r.get("id") or "") == sid),
                       None)
            workspace = str((row or {}).get("workspace") or "")
            if not workspace:
                return {"error": "Diese Sitzung ist nicht verfügbar. "
                                 "Kein anderer Agent wurde angesprochen."}
        if not workspace:
            return {"error": "profil <sitzung|ordner> frei|standard"}
        try:
            result = workspace_profiles.set_from_word(workspace, word)
        except PermissionError as exc:
            return {"error": str(exc)}
        except ValueError as exc:
            return {"error": str(exc)}
        return {"ok": True, **result}

    async def _cmd_schedules(self, params: dict) -> dict:
        """The schedules (cut) and the runs after ``after``, within one line."""
        from delfin.agent import schedules
        try:
            after = int(params.get("after") or 0)
        except ValueError:
            after = 0
        red = self._redactor()
        rows = [self._schedule_view(e) for e in schedules.entries()][:20]
        for row in rows:
            row["task"] = _truncate(row["task"], 100)
            row.pop("check", None)
        runs = [{"seq": r["seq"], "id": r["id"], "status": r["status"],
                 "text": red.text(_truncate(str(r.get("text") or ""), 1500)),
                 "costUsd": r.get("costUsd") or 0.0,
                 "task": red.text(_truncate(str(r.get("task") or ""), 120)),
                 "when": _truncate(str(r.get("when") or ""), 80),
                 "notify": r.get("notify") or "always"}
                for r in schedules.runs(after, limit=10)]
        head = schedules.last_seq()
        reply = {"ok": True, "host": rundir.hostname(), "schedules": rows,
                 "runs": runs, "head": head}
        # One line holds at most MAX_LINE_BYTES: later runs wait for the
        # next look, which asks after the last one sent.
        while len(json.dumps(reply, ensure_ascii=False).encode()) > MAX_LINE_BYTES - 512 and runs:
            runs.pop()
        while len(json.dumps(reply, ensure_ascii=False).encode()) > MAX_LINE_BYTES - 512 and rows:
            rows.pop()
        reply["last"] = runs[-1]["seq"] if runs else (after if after else head)
        return reply

    async def _cmd_schedule(self, params: dict) -> dict:
        """Make (as a proposal), confirm and manage a schedule by number.

        The folder is the host's: a session's own (exact enabled ID) or
        the phone folder; the phone never names one. ``ok`` confirms only
        the content whose check the phone was shown."""
        from delfin.agent import schedules
        action = str(params.get("action") or "")
        try:
            if action == "new":
                return await self._schedule_new(params)
            entry = schedules.get(str(params.get("id") or ""))
            if action == "show":
                return {"ok": True, **self._schedule_view(entry)}
            if action == "ok":
                check = str(params.get("check") or "")
                if entry["confirmed"]:
                    return {"ok": True, **self._schedule_view(entry)}
                if not check or check != schedules.fingerprint(entry):
                    return {"ok": False, "confirmationRequired": True,
                            **self._schedule_view(entry)}
                return {"ok": True, **self._schedule_view(
                    schedules.confirm(entry["id"], expect=check))}
            if action == "pause":
                return {"ok": True, **self._schedule_view(schedules.pause(entry["id"]))}
            if action == "resume":
                return {"ok": True, **self._schedule_view(schedules.resume(entry["id"]))}
            if action == "delete":
                schedules.delete(entry["id"])
                return {"ok": True, "id": entry["id"]}
            if action == "run":
                if entry.get("host") != rundir.hostname():
                    return {"error": "Dieser Zeitplan gehört zu einem anderen Rechner und läuft nur dort."}
                if not entry["confirmed"]:
                    return {"error": "Noch nicht bestätigt: zuerst ok mit seiner Nummer senden."}
                return {"ok": True, **self._schedule_view(schedules.request_run(entry["id"]))}
            return {"error": "Unbekannter Zeitplan-Befehl."}
        except schedules.ScheduleError as exc:
            return {"error": SCHEDULE_ERRORS.get(exc.code, SCHEDULE_ERRORS["when"]),
                    "code": exc.code}

    async def _schedule_new(self, params: dict) -> dict:
        from delfin.agent import schedules
        sid = str(params.get("session") or "")
        if sid:
            if not self._enabled(sid):
                return {"error": "SimpleX ist für diese Sitzung ausgeschaltet. In TUI/Dashboard mit /sync aktivieren."}
            listing = await self._request("sessions.list")
            row = next((r for r in listing.get("sessions", [])
                        if str(r.get("sessionId") or r.get("id") or "") == sid), None)
            folder = str((row or {}).get("workspace") or "")
            if not folder:
                return {"error": "Diese Sitzung ist nicht verfügbar. Kein anderer Agent wurde angesprochen."}
        else:
            root = os.environ.get(WORKSPACES_ENV, "").strip()
            if not root or not os.path.isabs(root) or not os.path.isdir(root):
                return {"error": "Zeitpläne aus diesem Chat brauchen einen Handy-Ordner auf dem Rechner. "
                                 "Im Sitzungs-Chat mit plan ZEIT AUFGABE anlegen oder in der TUI mit /schedule add."}
            folder = os.path.join(root, SCHEDULE_FOLDER)
            if os.path.lexists(folder):
                if os.path.islink(folder) or not os.path.isdir(folder):
                    return {"error": "Der Ordner für Zeitpläne ist kein normaler Ordner."}
            else:
                os.mkdir(folder, 0o700)
        now = time.monotonic()
        made = [t for t in getattr(self, "_schedules_made", []) if now - t < CREATE_WINDOW_S]
        if len(made) >= SCHEDULE_LIMIT:
            return {"error": "Höchstens fünf neue Zeitpläne in zehn Minuten. Bitte später erneut."}
        entry = schedules.add(str(params.get("text") or ""), workspace=folder,
                              created_by="phone", confirm=False)
        self._schedules_made = [*made, now]
        return {"ok": True, **self._schedule_view(entry)}


def main(argv: list[str] | None = None) -> int:
    """`python -m delfin.agent.service.bridge`: serve until interrupted."""
    asyncio.run(Bridge().serve())
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
