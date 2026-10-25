"""The command loop: only the paired contact commands, with limits.

What the phone may say:
    status                    -> what the sessions are doing
    sitzungen / sessions      -> the session list
    stop                      -> the agent service stops (sessions are saved)
    hilfe / help              -> this list
    @name text                -> the task goes to the named session
    1 / 2 / text as a reply   -> the answer to the newest question, by its
                                 number from the last `pending` message
    anything else             -> a task for the "phone" session

What the phone never does: approve. A confirm/plan in `pending` is shown
with where to approve (dashboard or terminal), and an answer that names
an approval is refused by the bridge as well.
"""

from __future__ import annotations

import json
import time
from collections import deque

from . import simplex
from .bridge_client import BridgeClient, BridgeDown, default_socket_path

#: The message the phone may send, and how often: at most RATE_LIMIT
#: messages per RATE_WINDOW_S, at most SIZE_LIMIT characters.
RATE_LIMIT = 15
RATE_WINDOW_S = 60.0
SIZE_LIMIT = 4_096
MAX_BRIDGE_TEXT = 2_000

HELP_TEXT = (
    "Befehle:\n"
    "status - was die Sitzungen tun\n"
    "sessions / sitzungen - die Liste der Sitzungen\n"
    "pending - offene Fragen und Hinweise\n"
    "stop - alle Agenten anhalten (Sitzungen bleiben erhalten)\n"
    "hilfe / help - diese Liste\n"
    "@name text - Auftrag an die Sitzung »name«\n"
    "Antwort auf eine Frage: Nummer (1, 2, ...) oder Text schicken\n"
    "alles andere: Auftrag an die Handysitzung\n"
    "Freigaben gibt es nur am Dashboard/Terminal, nie per Handy")


class Refused(RuntimeError):
    """A message the gateway will not act on."""


def _trim(text: str, limit: int = MAX_BRIDGE_TEXT) -> str:
    text = str(text or "")
    if len(text) <= limit:
        return text
    return text[:limit - 1] + "…"


class PhoneSession:
    """One message the phone sent, with its id: processed once."""

    def __init__(self, contact_id: str, text: str, msg_id: str):
        self.contact_id = contact_id
        self.text = text
        self.msg_id = msg_id


class Gateway:
    def __init__(self, *, contact_id: str, bridge: BridgeClient | None = None,
                 send_text=None, now=time.monotonic):
        self.contact_id = contact_id
        self.bridge = bridge or BridgeClient(default_socket_path())
        self.send_text = send_text          # async (cid, text) -> None
        self.now = now
        self._seen_ids: deque[str] = deque(maxlen=2048)
        self._rate: deque[float] = deque()
        self.questions: dict[int, str] = {}      # number -> attention item id
        self._selected: str | None = None        # the item a bare number chose

    # -- limits ---------------------------------------------------------------------------

    def _check_limits(self, message: PhoneSession) -> None:
        now = self.now()
        self._rate.append(now)
        while self._rate and now - self._rate[0] > RATE_WINDOW_S:
            self._rate.popleft()
        if len(self._rate) > RATE_LIMIT:
            raise Refused("zu viele Nachrichten; warte kurz")
        if len(message.text) > SIZE_LIMIT:
            raise Refused(f"Nachricht zu lang ({len(message.text)} > {SIZE_LIMIT})")
        if message.msg_id and message.msg_id in self._seen_ids:
            raise Refused("dieselbe Nachricht noch einmal")
        if message.msg_id:
            self._seen_ids.append(message.msg_id)

    # -- dispatch --------------------------------------------------------------------------

    async def handle(self, message: PhoneSession) -> str:
        """One message from the paired contact, the reply text."""
        if message.contact_id != self.contact_id:
            raise Refused("nicht gekoppelt")
        self._check_limits(message)
        text = message.text.strip()
        if not text:
            return ""
        if text.lower() in ("hilfe", "help"):
            return HELP_TEXT
        if text.lower() in ("status", "sitzungen", "sessions", "pending"):
            return self._bridge_command("status" if text.lower() == "status"
                                        else "pending" if text.lower() == "pending"
                                        else "sessions")
        if text.lower() == "stop":
            return self._bridge_command("stop")
        # @name text: the task goes to the named session.
        if text.startswith("@"):
            name, _, rest = text.partition(" ")
            if not rest.strip():
                return "Gib nach dem @Namen den Auftrag an."
            return self._bridge_task(rest.strip(), session=name[1:])
        # A number chooses the item it names (a question or, visible as
        # such, an approval); the answer follows as the next message and
        # the bridge refuses an approval verbatim. An unknown number is
        # not an answer: it becomes a task like any other text.
        if text.isdigit():
            number = int(text)
            if number in self.questions:
                self._selected = self.questions[number]
                return "Ausgewählt — jetzt den Antworttext schicken."
            return self._bridge_task(text)
        # Free text: an answer to the chosen item, else the first open
        # question (the bridge refuses a non-question), else a task.
        if self._selected:
            item_id, self._selected = self._selected, None
            return self._bridge_answer(item_id, text)
        if self.questions:
            try:
                reply = self.bridge.ask({"cmd": "pending"})
            except (BridgeDown, OSError):
                reply = {}
            items = reply.get("items") or []
            if items and items[0].get("kind") == "question":
                return self._bridge_answer(str(items[0].get("id") or ""), text)
        return self._bridge_task(text)

    # -- bridge calls, each cut and with its errors turned into text -----

    def _bridge_command(self, cmd: str) -> str:
        try:
            reply = self.bridge.ask({"cmd": cmd})
        except BridgeDown as exc:
            return f"Der Agentendienst ist nicht erreichbar ({exc})."
        except OSError as exc:
            return f"Der Agentendienst ist nicht erreichbar ({exc})."
        if reply.get("error"):
            return f"Fehler: {reply['error']}"
        return self._render(cmd, reply)

    def _bridge_task(self, text: str, session: str | None = None) -> str:
        command = {"cmd": "prompt", "text": _trim(text)}
        if session:
            command["session"] = _trim(session, 64)
        try:
            reply = self.bridge.ask(command)
        except (BridgeDown, OSError) as exc:
            return f"Der Agentendienst ist nicht erreichbar ({exc})."
        if reply.get("error"):
            return f"Fehler: {reply['error']}"
        name = reply.get("session") or session or "Handy"
        return (f"Auftrag an »{name}« übergeben."
                if reply.get("accepted") or reply.get("queued")
                else f"Auftrag an »{name}« weitergereicht.")

    def _bridge_answer(self, item_id: str, text: str) -> str:
        if not item_id:
            return "Die Frage ist nicht mehr offen."
        try:
            reply = self.bridge.ask({"cmd": "answer", "id": item_id,
                                     "text": _trim(text)})
        except (BridgeDown, OSError) as exc:
            return f"Der Agentendienst ist nicht erreichbar ({exc})."
        if reply.get("error"):
            return reply["error"]      # the refusal of an approval, verbatim
        return "Antwort weitergegeben."

    # -- rendering -----------------------------------------------------------------------

    def _render(self, cmd: str, reply: dict) -> str:
        if cmd == "status":
            rows = reply.get("sessions") or []
            if not rows:
                return "Keine Sitzungen."
            lines = [f"{r.get('title') or r.get('id')}: {r.get('state')}"
                     for r in rows]
            out = "\n".join(lines)
            self._remember_questions_from_pending()
            return _trim(out, 2_000)
        if cmd == "sessions":
            rows = reply.get("sessions") or []
            if not rows:
                return "Keine Sitzungen."
            return _trim("\n".join(
                f"{r.get('name') or r.get('id')}: {r.get('state')}"
                for r in rows), 2_000)
        if cmd == "pending":
            return self._render_pending(reply)
        if cmd == "stop":
            return ("Alle Agenten werden angehalten "
                    f"({reply.get('sessions')} Sitzungen enden und werden gesichert).")
        return json.dumps(reply, ensure_ascii=False)[:1_000]

    def _render_pending(self, reply: dict) -> str:
        items = reply.get("items") or []
        if not items:
            return "Nichts offen."
        self._last_number = 0
        self.questions = {}
        lines: list[str] = []
        for item in items[:5]:
            kind = item.get("kind")
            title = item.get("session") or ""
            text = item.get("text") or item.get("detail") or ""
            self._last_number += 1
            self.questions[self._last_number] = str(item.get("id") or "")
            if kind == "question":
                lines.append(f"{self._last_number}. Frage ({title}): "
                             f"{_trim(text, 200)}")
            elif kind in ("confirm", "plan"):
                lines.append(f"{self._last_number}. Freigabe nötig ({title}): "
                             f"{_trim(text, 200)}\n"
                             "→ am Dashboard/Terminal genehmigen (Antwort hier wird abgelehnt)")
            else:
                lines.append(f"{kind} ({title}): {_trim(text, 200)}")
        return _trim("\n".join(lines), 2_000)

    def _remember_questions_from_pending(self) -> None:
        # The numbers the phone answers with stay valid: after a status
        # the question set is refreshed quietly.
        try:
            reply = self.bridge.ask({"cmd": "pending"})
        except (BridgeDown, OSError):
            return
        self._render_pending(reply)


class Listener:
    """The run loop: daemon events in, gateway decisions, bridge calls,
    answers back to the paired contact only."""

    def __init__(self, *, contact_id: str, gateway: Gateway):
        self.contact_id = contact_id
        self.gateway = gateway

    async def run(self, daemon: simplex.Daemon) -> None:
        await daemon.subscribe()
        async for event in daemon.events():
            cid_text = simplex.text_of(event)
            if cid_text is None:
                # A contact request from a stranger: refused, logged.
                if simplex.new_contact_id(event) and \
                        simplex.new_contact_id(event) != self.contact_id:
                    continue
                continue
            cid, text = cid_text
            if cid != self.contact_id:
                continue                     # only the paired contact
            message = PhoneSession(cid, text, str(event.get("itemId") or "")
                                   or str(time.time()))
            try:
                reply = await self.gateway.handle(message)
            except Refused as exc:
                reply = f"Aufgabe abgelehnt: {exc}"
            except Exception as exc:  # noqa: BLE001
                reply = f"Interner Fehler: {type(exc).__name__}: {exc}"
            if reply:
                await daemon.send_text(self.contact_id, _trim(reply, 2_000))
