"""The dashboard as a window of the agent service.

The agent service (``delfin.agent.service.supervisor``) holds every session
of one person on one host; the TUI and the messenger bridge are its windows,
and so is the dashboard through this module. Everything goes over the
service's own Unix socket (``rundir.socket_path()``: the 0700 runtime
directory of this process's DELFIN_RUNTIME_DIR, so in the launcher's people
mode each person reaches their own service only) and the existing protocol
(``delfin.agent.service.protocol``). No listener is opened here.

* :class:`LiveClient` -- one connection: list, create, attach (recent turns
  replayed, then live), prompt, cancel, answer, rename, detach, stop. Events
  go to subscribers on the reader thread and are not kept: a dashboard
  kernel lives for days.
* :class:`Transcript` -- one session as the dashboard shows it, built from
  its replay and live events; durable events are applied once (``seq``).
* :class:`RunningList` -- every session of the service, kept current from
  ``session.summary`` / ``session.removed`` events (the service sends them
  to every window); nothing asks ``sessions.list`` again.

Rejected: polling sessions.list on the sidebar's tick. The service already
pushes every change, and the cluster's shared home makes the presence read
behind sessions.list the slow part.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path
from typing import Any, Callable

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceClient, ServiceError

__all__ = ["ServiceUnavailable", "ServiceError", "LiveClient", "connect", "Transcript",
           "RunningList", "state_word", "age_text", "client_words"]

#: What the dashboard calls itself in its hello: what permission.resolved
#: names when an answer came from here.
CLIENT_NAME = "delfin-dashboard"

#: Turns replayed when a session is opened (the TUI's replayTurns).
REPLAY_TURNS = 20

#: A session awake this long since anyone wrote to it (summary lastInput,
#: else its start) is marked "running N h without input" in the list.
IDLE_NOTICE_S = 3600.0


class ServiceUnavailable(RuntimeError):
    """No agent service can be reached or started here; the text says why."""


class LiveClient(ServiceClient):
    """The dashboard's connection to the person's agent service."""

    def __init__(self, path: Path | None = None, *, timeout: float = 10.0) -> None:
        # Set before the base starts its reader thread.
        self._subscribers: list[Callable[[dict], None]] = []
        self._closers: list[Callable[[], None]] = []
        self._sub_lock = threading.Lock()
        super().__init__(path, timeout=timeout, client=CLIENT_NAME)

    # -- the reader: events to subscribers, never kept --------------------------------
    def _read(self) -> None:
        try:
            for message in protocol.read_messages(self._file):
                if message.get("method") == "event":
                    event = message.get("params") or {}
                    with self._sub_lock:
                        subscribers = list(self._subscribers)
                    for fn in subscribers:
                        try:
                            fn(event)
                        except Exception:
                            pass             # one view's fault never stops the others
                    continue
                with self._cond:
                    self._responses[message.get("id")] = message
                    self._cond.notify_all()
        except (OSError, ValueError):
            pass
        finally:
            with self._cond:
                self.closed = True
                self._cond.notify_all()
            with self._sub_lock:
                closers = list(self._closers)
            for fn in closers:
                try:
                    fn()
                except Exception:
                    pass

    def subscribe(self, fn: Callable[[dict], None]) -> Callable[[], None]:
        """Call ``fn(event)`` for every event (on the reader thread; keep it
        short). Returns the unsubscribe."""
        with self._sub_lock:
            self._subscribers.append(fn)

        def unsubscribe() -> None:
            with self._sub_lock:
                if fn in self._subscribers:
                    self._subscribers.remove(fn)
        return unsubscribe

    def on_close(self, fn: Callable[[], None]) -> None:
        """Call ``fn()`` once the connection is gone (the service stopped)."""
        with self._sub_lock:
            self._closers.append(fn)

    # -- requests -------------------------------------------------------------------
    def sessions(self) -> dict:
        return self.request("sessions.list")

    def create(self, cwd: str, **settings: Any) -> dict:
        params = {k: v for k, v in settings.items() if v not in (None, "", [])}
        params["cwd"] = str(cwd)
        params.setdefault("attach", False)
        return self.request("session.create", params, timeout=200)

    def attach(self, sid: str, turns: int = REPLAY_TURNS) -> dict:
        return self.request("session.attach", {"sessionId": sid, "turns": int(turns)},
                            timeout=200)

    def detach(self, sid: str) -> None:
        if not self.closed:
            self.request("session.detach", {"sessionId": sid}, timeout=10)

    def prompt(self, sid: str, text: str, how: str = "queue") -> dict:
        return self.request("session.prompt", {"sessionId": sid, "text": text, "how": how})

    def cancel(self, sid: str) -> dict:
        return self.request("session.cancel", {"sessionId": sid})

    def reply(self, sid: str, request_id: str, **answer: Any) -> dict:
        params = {"sessionId": sid, "requestId": request_id}
        params.update({k: v for k, v in answer.items() if v not in (None, "", [])})
        return self.request("permission.reply", params)

    def rename(self, sid: str, title: str) -> str:
        return str(self.request("session.setTitle",
                                {"sessionId": sid, "title": str(title or "")})["title"])

    def posture(self, sid: str, mode: str, *, confirmed: bool = False) -> str:
        params: dict[str, Any] = {"sessionId": sid, "mode": mode}
        if confirmed:
            params["confirmed"] = True
        return str(self.request("session.setMode", params).get("approval") or "")

    def choices(self, sid: str, what: str) -> dict:
        return self.request("session.choices", {"sessionId": sid, "what": what}, timeout=60)

    def stop(self, sid: str) -> dict:
        """End the session (it is saved): the explicit stop, never a close."""
        return self.request("session.shutdown", {"sessionId": sid}, timeout=60)


def connect(*, start: bool = True, path: Path | None = None) -> LiveClient:
    """A connection to this person's agent service; with ``start`` the
    service is started when none runs (``launcher.ensure_service``, the
    TUI's own way). Raises :class:`ServiceUnavailable` with the reason."""
    try:
        from delfin.agent.service import client as _client
        from delfin.agent.service import launcher, rundir

        target = Path(path) if path is not None else rundir.socket_path()
        hello, why = _client.probe(target)
        if hello is None:
            if not start or path is not None:
                raise ServiceUnavailable(
                    f"no agent service answers here ({launcher.why_in_words(why)})")
            launcher.ensure_service()
        return LiveClient(target)
    except ServiceUnavailable:
        raise
    except Exception as exc:                                  # noqa: BLE001
        raise ServiceUnavailable(f"the agent service could not be reached: {exc}") from None


# ---------------------------------------------------------------------------
# What the dashboard shows
# ---------------------------------------------------------------------------

def state_word(summary: dict) -> str:
    """One of four words: working, needs you, done, asleep."""
    state = str(summary.get("state") or "")
    if state == "sleeping":
        return "asleep"
    if int(summary.get("pending") or 0) > 0 or state == "needs_input":
        return "needs you"
    if state in ("running", "starting"):
        return "working"
    return "done"


_STATE_ORDER = {"needs you": 0, "working": 1, "done": 2, "asleep": 3}


def age_text(seconds: float) -> str:
    """How long ago, in one short word: now, 4 min, 3 h, 2 d."""
    s = max(0.0, float(seconds or 0.0))
    if s < 60:
        return "now"
    if s < 3600:
        return f"{int(s // 60)} min"
    if s < 86400:
        return f"{int(s // 3600)} h"
    return f"{int(s // 86400)} d"


#: How the dashboard names the window an answer came from.
_CLIENT_WORDS = {"delfin-dashboard": "the dashboard", "delfin-tui": "the terminal",
                 "delfin-messenger-bridge": "the phone"}


def client_words(client: str, by: str = "client") -> str:
    """Who answered, in words: a window by its hello name, else how it ended."""
    if by == "client":
        return _CLIENT_WORDS.get(client, client or "another window")
    return {"outside": "outside the session", "expired": "the time limit",
            "abort": "the session's end"}.get(by, by or "someone")


def _folder(workspace: str) -> str:
    return Path(str(workspace or "")).name or str(workspace or "")


class RunningList:
    """Every session of the service, current through its events."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._summaries: dict[str, dict] = {}
        self.version = 0

    @classmethod
    def of(cls, listing: dict) -> "RunningList":
        out = cls()
        out.load(listing)
        return out

    def load(self, listing: dict) -> None:
        with self._lock:
            self._summaries = {str(s.get("sessionId")): dict(s)
                               for s in (listing or {}).get("sessions") or []
                               if s.get("sessionId")}
            self.version += 1

    def apply(self, event: dict) -> bool:
        kind = event.get("type")
        sid = str(event.get("session") or "")
        with self._lock:
            if kind == "session.summary" and isinstance(event.get("summary"), dict):
                summary = dict(event["summary"])
                sid = str(summary.get("sessionId") or sid)
                if not sid:
                    return False
                if summary.get("state") == "ended":
                    changed = self._summaries.pop(sid, None) is not None
                else:
                    changed = self._summaries.get(sid) != summary
                    self._summaries[sid] = summary
            elif kind == "session.removed":
                changed = self._summaries.pop(sid, None) is not None
            elif kind == "service.stopping":
                changed = bool(self._summaries)
                self._summaries.clear()
            else:
                return False
            if changed:
                self.version += 1
            return changed

    def summary(self, sid: str) -> dict | None:
        with self._lock:
            got = self._summaries.get(sid)
            return dict(got) if got else None

    def rows(self, now: float | None = None) -> list[dict]:
        """One row per session: needs you first, then working, done,
        asleep; within a state the newest session first (rows keep their
        place while sessions work)."""
        now = time.time() if now is None else now
        with self._lock:
            summaries = [dict(s) for s in self._summaries.values()]
        rows = []
        for s in summaries:
            sid = str(s["sessionId"])
            word = state_word(s)
            # A session runs on when its window is gone: one nobody wrote
            # to for an hour is said so, beside its Stop.
            idle = ""
            if word != "asleep":
                since = float(s.get("lastInput") or s.get("createdAt") or now)
                if now - since >= IDLE_NOTICE_S:
                    idle = f"running {age_text(now - since)} without input"
            rows.append({
                "sessionId": sid,
                "name": str(s.get("name") or protocol.display_name(
                    s.get("title"), s.get("workspace"), sid)),
                "shortId": sid[:8],
                "workspace": str(s.get("workspace") or ""),
                "folder": _folder(str(s.get("workspace") or "")),
                "model": str(s.get("model") or ""),
                "state": word,
                "needs": word == "needs you",
                "waiting": str(s.get("waiting") or ""),
                "activity": str(s.get("activity") or ""),
                "age": age_text(now - float(s.get("lastActivity") or now)),
                "idle": idle,
                "createdAt": float(s.get("createdAt") or 0.0),
            })
        rows.sort(key=lambda r: (_STATE_ORDER[r["state"]], -r["createdAt"]))
        return rows


class Transcript:
    """One session as the dashboard shows it: what was said, what runs,
    what waits for an answer, and its settings."""

    def __init__(self, sid: str) -> None:
        self.sid = sid
        self.lock = threading.RLock()
        self.entries: list[dict] = []
        self.pending: dict[str, dict] = {}
        self.info: dict = {}
        self.summary: dict = {}
        self.state = "starting"
        self.wait = ""
        self.ended = ""
        self.last_seq = 0
        self.version = 0
        self._parts: dict[str, dict] = {}
        self._ended_parts: set[str] = set()
        self._tools: dict[str, dict] = {}
        #: Live events that came before the attach's answer was loaded: the
        #: reader thread hands them over while the attaching thread is still
        #: waking up. Applied after the replay (durable ones once, by seq).
        self._early: list[dict] | None = []

    # -- what it is ---------------------------------------------------------------
    @property
    def name(self) -> str:
        with self.lock:
            if self.summary.get("name"):
                return str(self.summary["name"])
            return protocol.display_name(self.info.get("title"), self.info.get("workspace"),
                                         self.sid)

    @property
    def busy(self) -> bool:
        return self.state in ("running", "starting") and not self.pending

    @property
    def needs_you(self) -> bool:
        return bool(self.pending) or self.state == "needs_input"

    # -- filling it ---------------------------------------------------------------
    def load(self, attached: dict) -> None:
        """An attach's answer: the replay, then where the session stands."""
        with self.lock:
            self.entries = []
            self.pending = {}
            self._parts, self._ended_parts, self._tools = {}, set(), {}
            self.last_seq = 0
            for event in attached.get("events") or []:
                self._apply(event, replay=True)
            self.info.update(attached.get("session") or {})
            self.state = str(attached.get("state") or self.state)
            self.pending = {str(a.get("requestId")): a for a in attached.get("pending") or []
                            if a.get("requestId")}
            if attached.get("summary"):
                self.summary = dict(attached["summary"])
            early, self._early = self._early or [], None
            for event in early:
                self._apply(event, replay=False)
            self.version += 1

    def apply(self, event: dict) -> bool:
        """A live event; True when what is shown changed."""
        if str(event.get("session") or "") != self.sid:
            return False
        with self.lock:
            if self._early is not None:
                self._early.append(event)         # before load(): kept for it
                return False
            changed = self._apply(event, replay=False)
            if changed:
                self.version += 1
            return changed

    def _entry(self, role: str, text: str, **extra: Any) -> dict:
        entry = {"role": role, "text": text, **extra}
        self.entries.append(entry)
        return entry

    def _apply(self, event: dict, *, replay: bool) -> bool:
        kind = str(event.get("type") or "")
        seq = event.get("seq")
        if isinstance(seq, int) and seq > 0:
            if seq <= self.last_seq:
                return False                      # in the replay already
            self.last_seq = seq
        extra = {"seq": seq} if isinstance(seq, int) and seq > 0 else {}
        if kind in ("session.info", "session.updated"):
            self.info.update({k: v for k, v in event.items()
                              if k not in ("type", "session", "seq", "ts")})
            if "title" in event and self.summary:
                self.summary["title"] = event.get("title") or None
                self.summary["name"] = protocol.display_name(
                    event.get("title"), self.info.get("workspace"), self.sid)
            return True
        if kind == "session.summary":
            self.summary = dict(event.get("summary") or {})
            self.state = str(self.summary.get("state") or self.state)
            return True
        if kind == "session.state":
            self.state = str(event.get("state") or self.state)
            if self.state != "running":
                self.wait = ""
            return True
        if kind in ("session.ended", "session.removed"):
            self.state = "ended"
            self.ended = str(event.get("reason") or "ended")
            self.pending.clear()
            return True
        if kind == "user.message":
            self._entry("user", str(event.get("text") or ""),
                        source=str(event.get("source") or "user"), **extra)
            return True
        if kind == "text.started":
            part = str(event.get("partId") or "")
            if part and part not in self._parts:
                self._parts[part] = self._entry("agent", "", streaming=True, part=part)
            return True
        if kind == "text.delta":
            part = str(event.get("partId") or "")
            if part in self._ended_parts:
                return False                      # a part the replay already ended
            entry = self._parts.get(part)
            if entry is None:
                entry = self._parts[part] = self._entry("agent", "", streaming=True, part=part)
            entry["text"] += str(event.get("delta") or "")
            return True
        if kind == "text.ended":
            part = str(event.get("partId") or "")
            self._ended_parts.add(part)
            entry = self._parts.pop(part, None)
            if entry is None:
                entry = self._entry("agent", "", part=part)
            entry.update(text=str(event.get("text") or ""), streaming=False, **extra)
            return True
        if kind == "reasoning.ended":
            text = str(event.get("text") or "")
            if text.strip():
                self._entry("thinking", text, **extra)
            return bool(text.strip())
        if kind == "tool.started":
            call = str(event.get("callId") or "")
            entry = self._entry("tool", _tool_line(event), name=str(event.get("name") or ""),
                                done=False, ok=True, **extra)
            if call:
                self._tools[call] = entry
            return True
        if kind == "tool.ended":
            entry = self._tools.pop(str(event.get("callId") or ""), None)
            meta = event.get("meta") or {}
            if entry is None:
                entry = self._entry("tool", _tool_line(event), name=str(event.get("name") or ""))
            entry.update(done=True, ok=not (isinstance(meta, dict) and meta.get("error")),
                         output=str(event.get("output") or "")[:2000])
            return True
        if kind == "tool.denied":
            self._entry("tool", f"{event.get('name') or 'tool'}: refused", done=True, ok=False,
                        **extra)
            return True
        if kind in ("notice", "error"):
            self._entry(kind, str(event.get("text") or ""), **extra)
            return True
        if kind == "context":
            self._entry("context", str(event.get("text") or ""),
                        title=str(event.get("title") or ""), **extra)
            return True
        if kind == "wait":
            self.wait = str(event.get("text") or "")
            return True
        if kind == "turn.ended":
            self.wait = ""
            for entry in self._parts.values():
                entry["streaming"] = False
            self._parts.clear()
            if event.get("error"):
                self._entry("error", str(event["error"]), **extra)
            return True
        if kind == "permission.asked":
            rid = str(event.get("requestId") or "")
            if rid:
                self.pending[rid] = dict(event)
            return True
        if kind == "permission.resolved":
            rid = str(event.get("requestId") or "")
            asked = self.pending.pop(rid, None)
            by = str(event.get("by") or "")
            client = str(event.get("client") or "")
            what = _asked_line(asked) if asked else "the request"
            self._entry("resolved", f"{what}: {event.get('decision') or 'answered'} "
                        f"by {client_words(client, by)}",
                        requestId=rid, client=client, by=by,
                        decision=str(event.get("decision") or ""), **extra)
            return True
        return False


def _tool_line(event: dict) -> str:
    name = str(event.get("name") or "tool")
    args = event.get("input") or {}
    if isinstance(args, dict):
        for key in ("command", "path", "file_path", "pattern", "query", "url", "description"):
            value = args.get(key)
            if isinstance(value, str) and value.strip():
                first = " ".join(value.split())
                return f"{name}  {first[:120]}"
    return name


def _asked_line(asked: dict) -> str:
    if asked.get("kind") == "ask":
        return str(asked.get("question") or "a question")[:120]
    if asked.get("kind") == "plan":
        return "the plan"
    return f"{asked.get('tool') or 'a tool'}"
