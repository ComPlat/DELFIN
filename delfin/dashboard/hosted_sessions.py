"""An agent tab's session, hosted in the person's agent service.

The session runs where it always ran -- in this kernel, on the agent tab's
own turn code -- and the service holds what every window shares: its
name, its transcript, whether it runs, Keep session
(``delfin.agent.service.host``, ``supervisor.HostedSession``). Another
dashboard, the TUI and the phone list it, open it, send to it and answer
it; what they send is run here exactly as if typed in this tab.

    host = TabHost(adapter, sid=..., cwd=..., title=..., keep=...)
    host.start()          # False: no service, or it will not host it
    host.release(end=False)

Nothing in the agent tab is changed to be hosted. It is watched and
driven through ``adapter`` (``tab_agent.create_tab``'s ``host_adapter``):

* what it says -- ``state["chat_messages"]``, ``state["streaming"]`` and
  its open questions (``messenger.dashboard.DashboardSession.pending``,
  the phone's adapter) -- is read every POLL_S and said as protocol events
  (:class:`ChatFollower`); the TUI shows them as it shows a worker's;
* the chat as the tab renders it (``chat_html.value``) goes to the
  person's other dashboards whole, as ``dashboard.chat`` diffs: they show
  the same markup with the same page scripts -- streaming text, tool rows,
  cards, the 3D viewer;
* a request is answered by the tab's own callbacks: a prompt goes through
  its send (slash commands, queueing and steering as typed; the box's
  draft is kept), a cancel through its Stop, an answer through the
  phone's adapter, a permission change through its permission selector.

Rejected: events emitted from inside the turn code (fifty call sites in a
17,594-line closure, each a place to forget); a second renderer of the
chat in other dashboards (they would show a different chat).
"""

from __future__ import annotations

import re
import threading
import time
from html.parser import HTMLParser
from typing import Any, Callable, Optional

#: How often the tab is looked at.
POLL_S = 0.15

#: At most this often the chat HTML goes out while it changes (the tab
#: renders about eight times a second while streaming).
CHAT_EVERY_S = 0.12

#: The permission selector's profiles by the posture each one is
#: (agent_sessions._PROFILE_START, and Bypass).
PROFILE_OF_MODE = {"plan": "plan", "default": "ask_all", "acceptEdits": "repo_free",
                   "bypassPermissions": "all_free"}
MODE_OF_PROFILE = {v: k for k, v in PROFILE_OF_MODE.items()}

_KEEP_LINE = re.compile(r"^\s*/keep(?:\s+(on|off))?\s*$", re.IGNORECASE)


class _Text(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []
        self.name = ""
        self._in_name = False

    def handle_starttag(self, tag, attrs):
        if tag in ("br", "div", "p", "pre", "li", "tr"):
            self.parts.append("\n")
        if "tool-name" in (dict(attrs).get("class") or "").split():
            self._in_name = True

    def handle_endtag(self, tag):
        if self._in_name and tag == "span":
            self._in_name = False

    def handle_data(self, data):
        if self._in_name and not self.name:
            self.name = data.strip()
        self.parts.append(data)


def html_text(markup: str) -> tuple[str, str]:
    """``(tool name, text)`` of a chat line's HTML: the ``.tool-name`` span
    (empty when it has none) and all of its text, tags dropped."""
    parser = _Text()
    try:
        parser.feed(str(markup or ""))
        parser.close()
    except Exception:
        return "", re.sub(r"<[^>]+>", "", str(markup or ""))
    text = re.sub(r"\n{3,}", "\n\n", "".join(parser.parts)).strip()
    return parser.name, text


def chat_diff(old: str, new: str) -> tuple[int, str]:
    """``(keep, tail)``: ``new == old[:keep] + tail``, ``keep`` as long as
    the two agree. A binary search over slice comparisons -- each one a
    memory compare -- so a 2 MB chat costs some twenty compares, not two
    million steps."""
    lo, hi = 0, min(len(old), len(new))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if old[:mid] == new[:mid]:
            lo = mid
        else:
            hi = mid - 1
    return lo, new[lo:]


class ChatFollower:
    """The agent tab's chat, said as protocol events.

    ``step(messages, streaming, pending)`` is called with the tab as it is
    now; each change is said once, in order, through ``emit(type, **fields)``:
    a message once it is settled (an answer still streaming is said as
    ``text.delta``s and then ended; a live line that rewrites itself waits
    until something follows it or the turn is over), a turn's start and
    end, the state, open questions as they come and go.
    """

    def __init__(self, emit: Callable[..., Any], *, published: int = 0,
                 status: Optional[Callable[[], dict]] = None) -> None:
        self.emit = emit
        self.published = int(published)
        self.status = status or (lambda: {})
        self.streaming = False
        self.turn = 0
        self.turn_id = ""
        self.turn_started = 0.0
        self.turn_tools = 0
        self.turn_base: dict = {}
        self.part: Optional[dict] = None      # {"msg": id, "id": partId, "sent": n}
        self.parts = 0
        self.calls = 0
        self.asked: dict[str, dict] = {}
        self.state = "idle"
        self.repliers: dict[str, str] = {}

    # -- one look ---------------------------------------------------------------------
    def step(self, messages: list, streaming: bool, pending: list[dict]) -> None:
        messages = list(messages or ())
        streaming = bool(streaming)
        if streaming and not self.streaming:
            self._publish(messages, streaming=False, hold_last=True)
            self._turn_start(messages)
        if len(messages) < self.published:
            # Cleared (/clear, a new conversation): nothing to take back.
            self.published = len(messages)
            self.part = None
            self.emit("notice", text="The chat was cleared in the dashboard.")
        self._publish(messages, streaming=streaming, hold_last=streaming)
        if self.streaming and not streaming:
            self._turn_end()
        self.streaming = streaming
        self._questions(pending)
        state = "needs_input" if self.asked else ("running" if streaming else "idle")
        if state != self.state:
            self.state = state
            self.emit("session.state", state=state)

    def _turn_start(self, messages: list) -> None:
        self.turn += 1
        self.turn_id = f"dt-{int(time.time() * 1000):x}-{self.turn}"
        self.turn_started = time.monotonic()
        self.turn_tools = 0
        self.turn_base = dict(self.status() or {})
        last_user = next((m for m in reversed(messages) if m.get("role") == "user"), {})
        source = "harness" if last_user.get("origin") == "event" else "user"
        self.emit("turn.started", turnId=self.turn_id, source=source)

    def _turn_end(self) -> None:
        now = dict(self.status() or {})

        def grew(key: str, kind=int):
            try:
                return kind(max(0, kind(now.get(key) or 0) - kind(self.turn_base.get(key) or 0)))
            except (TypeError, ValueError):
                return kind(0)

        self.emit("turn.ended", turnId=self.turn_id or "dt",
                  durationMs=int((time.monotonic() - self.turn_started) * 1000),
                  inputTokens=grew("input_tokens"), outputTokens=grew("output_tokens"),
                  costUsd=float(grew("cost_usd", float)), toolCalls=self.turn_tools)

    # -- messages -----------------------------------------------------------------------
    def _publish(self, messages: list, *, streaming: bool, hold_last: bool) -> None:
        n = len(messages)
        while self.published < n:
            i = self.published
            msg = messages[i] if isinstance(messages[i], dict) else {}
            last = i == n - 1
            role = str(msg.get("role") or "")
            if (streaming and msg.get("_streaming")) or (last and hold_last and (
                    msg.get("_live_key") or role in ("assistant", "thinking"))):
                if streaming and role == "assistant":
                    self._stream(msg)
                return                       # said once it is settled
            self._say(msg)
            self.published += 1

    def _stream(self, msg: dict) -> None:
        text = str(msg.get("content") or "")
        if self.part is None or self.part["msg"] is not msg:
            self.parts += 1
            self.part = {"msg": msg, "id": f"dp-{self.parts}", "sent": 0}
            self.emit("text.started", turnId=self.turn_id, partId=self.part["id"])
        sent = self.part["sent"]
        if len(text) > sent:
            self.emit("text.delta", partId=self.part["id"], delta=text[sent:])
            self.part["sent"] = len(text)

    def _say(self, msg: dict) -> None:
        role = str(msg.get("role") or "")
        content = msg.get("content")
        text = content if isinstance(content, str) else ""
        if role == "user":
            self.emit("user.message", text=text,
                      source="harness" if msg.get("origin") == "event" else "user")
        elif role == "assistant":
            part = self.part if (self.part and self.part["msg"] is msg) else None
            self.part = None
            if part is None:
                if not text.strip():
                    return                   # a turn that answered only in tool calls
                self.parts += 1
                part = {"id": f"dp-{self.parts}", "sent": 0}
                self.emit("text.started", turnId=self.turn_id, partId=part["id"])
            elif len(text) > part["sent"]:
                self.emit("text.delta", partId=part["id"], delta=text[part["sent"]:])
            self.emit("text.ended", turnId=self.turn_id, partId=part["id"], text=text)
        elif role == "thinking":
            self.parts += 1
            self.emit("reasoning.ended", turnId=self.turn_id, partId=f"dp-{self.parts}",
                      text=text)
        elif role == "tool":
            name, line = html_text(text)
            self.calls += 1
            self.turn_tools += 1
            call = f"dc-{self.calls}"
            rest = line[len(name):].strip() if name and line.startswith(name) else line
            tool = "Bash" if name == "$" else (name or "tool")
            fields = {"command": rest} if tool == "Bash" else {"description": rest}
            self.emit("tool.started", turnId=self.turn_id, callId=call, name=tool, input=fields)
            self.emit("tool.ended", turnId=self.turn_id, callId=call, name=tool,
                      output="", meta={"ok": True, "chars": 0, "truncated": False})
        else:
            if role in ("approval", "gate") or "<" in text:
                text = html_text(text)[1] if "<" in text else text
            if not text.strip() and role == "delegate":
                text = str(msg.get("title") or msg.get("description") or "")
            if text.strip():
                self.emit("notice", text=text)

    # -- questions ----------------------------------------------------------------------
    def _questions(self, pending: list[dict]) -> None:
        now = {str(p.get("requestId")): p for p in pending or () if p.get("requestId")}
        for rid, asked in now.items():
            if rid in self.asked:
                continue
            fields = {k: v for k, v in asked.items() if k != "requestId"}
            fields.setdefault("tool", "dashboard")
            fields.setdefault("preview", "")
            fields.setdefault("options", [])
            fields["flags"] = dict(fields.get("flags") or {})
            self.emit("permission.asked", requestId=rid, **fields)
        for rid in [r for r in self.asked if r not in now]:
            who = self.repliers.pop(rid, "")
            self.emit("permission.resolved", requestId=rid, decision="answered",
                      by="client", client=who or "delfin-dashboard")
        self.asked = now


def shared_questions(pending: list[dict]) -> list[dict]:
    """What the tab has open that the other windows are asked.

    A question the tab only guessed from its answer's last line (a greeting
    that ends in "?") with no choices to pick is left out: in the dashboard
    it only focuses the box, but the terminal showed it as a question that
    took Enter and Esc -- the session could neither be answered nor stopped
    (2026-11-16). A guessed question with choices goes out marked
    ``detected`` (it does not hold the turn)."""
    return [p for p in pending or () if not (p.get("detected") and not p.get("choices"))]


class TabHost:
    """One agent tab's session hosted in the agent service.

    ``adapter``: the tab's ``host_adapter`` (tab_agent.create_tab) --
    ``state``, ``send``, ``stop``, ``answer``, ``plan``, ``approval``,
    ``question``, ``set_permission``, ``choices``, ``info``, ``chat_html``,
    ``status``. ``on_title(title)``: a rename from any window (the session
    list renames its row and the saved session); ``on_keep(on)``: ``/keep``
    typed in any window; ``on_end()``: ``session.shutdown`` from a window
    (End session: the list closes the tab, saved); ``on_lost()``: the link
    went away without this side letting go (the service stopped).
    """

    def __init__(self, adapter: dict, *, sid: str, cwd: str, title: str = "",
                 keep: bool = True, on_title: Callable[[str], None] = lambda _t: None,
                 on_keep: Callable[[bool], None] = lambda _k: None,
                 on_end: Callable[[], None] = lambda: None,
                 on_lost: Callable[[], None] = lambda: None,
                 path=None, fake: bool = False) -> None:
        from delfin.messenger.dashboard import DashboardSession

        self.adapter = adapter
        self.state = adapter["state"]
        self.sid = str(sid)
        self.cwd = str(cwd)
        self.title = str(title or "")
        self.keep = bool(keep)
        self.on_title, self.on_keep, self.on_end, self.on_lost = on_title, on_keep, on_end, on_lost
        self.path = path
        self.fake = fake
        self.link = None
        self.follower: Optional[ChatFollower] = None
        self.why = ""
        self._stop = threading.Event()
        self._released = False
        self._thread: Optional[threading.Thread] = None
        self._chat = {"sent": "", "v": 0, "at": 0.0}
        self._info: dict = {}
        self.session = DashboardSession(
            self.sid, self.state, send=adapter["send"], stop=adapter["stop"],
            answer=adapter["answer"], plan=adapter["plan"], approval=adapter["approval"],
            question=adapter["question"])

    # -- the link -----------------------------------------------------------------------
    def start(self) -> bool:
        """Host the session; False (``why`` says) when the service cannot."""
        from delfin.agent.service.host import HostLink, HostRefused

        info = self._now_info()
        try:
            self.link = HostLink(self.sid, self.cwd, handle=self.handle, path=self.path,
                                 title=self.title, keep=self.keep, info=info,
                                 fake=self.fake, on_close=self._closed)
        except HostRefused as exc:
            self.why = exc.message
            return False
        except Exception as exc:                       # noqa: BLE001
            self.why = f"{type(exc).__name__}: {exc}"
            return False
        self._info = info
        messages = list(self.state.get("chat_messages") or ())
        self.link.emit("session.info", sessionId=self.sid, workspace=self.cwd,
                       banner=[], notices=[], **info)
        # History once: a session that never ran in the service has none
        # in its log; one that did has it there already.
        self.follower = ChatFollower(self.link.emit,
                                     published=0 if self.link.log_was_empty else len(messages),
                                     status=self.adapter.get("status") or (lambda: {}))
        self._tick()
        self._thread = threading.Thread(target=self._run, daemon=True,
                                        name=f"tab-host-{self.sid[:8]}")
        self._thread.start()
        return True

    @property
    def hosted(self) -> bool:
        return self.link is not None and not self.link.closed

    def _closed(self) -> None:
        self._stop.set()
        if not self._released:
            try:
                self.on_lost()
            except Exception:
                pass

    def release(self, *, end: bool = False) -> None:
        """Let go on purpose (x in the list, after the tab saved)."""
        self._released = True
        self._stop.set()
        if self.link is not None:
            self.link.release(end=end)

    def close(self) -> None:
        """The link goes without a word (this kernel stops hosting)."""
        self._released = True
        self._stop.set()
        if self.link is not None:
            self.link.close()

    def set_keep(self, on: bool) -> None:
        self.keep = bool(on)
        if self.link is not None:
            self.link.set_keep(on)

    # -- watching the tab -------------------------------------------------------------
    def _run(self) -> None:
        while not self._stop.wait(POLL_S):
            try:
                self._tick()
            except Exception:
                pass

    def _now_info(self) -> dict:
        try:
            info = dict(self.adapter["info"]() or {})
        except Exception:
            info = {}
        if self.title:
            info["title"] = self.title
        info["suggestions"] = self._suggestions()
        return {k: v for k, v in info.items() if v not in (None,)}

    def _suggestions(self) -> list[str]:
        """What the tab offers in its empty message box (its placeholder,
        while it is an offer or an open task and not the resting hint):
        the other windows show the same, Tab takes it."""
        box = self.adapter.get("input")
        try:
            if box is None or "delfin-agent-input-suggests" not in tuple(box._dom_classes):
                return []
            line = " ".join(str(box.placeholder or "").split())
        except Exception:
            return []
        return [line[:300]] if line else []

    def _tick(self) -> None:
        link = self.link
        if link is None or link.closed or self.follower is None:
            return
        try:
            pending = self.session.pending()
        except Exception:
            pending = []
        pending = shared_questions(pending)
        self.follower.step(self.state.get("chat_messages") or (),
                           bool(self.state.get("streaming")), pending)
        info = self._now_info()
        changed = {k: v for k, v in info.items() if self._info.get(k) != v}
        if changed:
            self._info = info
            link.emit("session.updated", **changed)
        self._send_chat()

    def _send_chat(self, *, whole: bool = False) -> None:
        widget = self.adapter.get("chat_html")
        if widget is None or self.link is None:
            return
        value = str(getattr(widget, "value", "") or "")
        chat = self._chat
        if value == chat["sent"] and chat["v"]:
            return
        now = time.monotonic()
        if chat["v"] and now - chat["at"] < CHAT_EVERY_S:
            return
        version = chat["v"] + 1
        if whole or not chat["v"]:
            ok = self.link.emit_live("dashboard.chat", v=version, html=value)
        else:
            keep, tail = chat_diff(chat["sent"], value)
            ok = self.link.emit_live("dashboard.chat", v=version, base=chat["v"],
                                     keep=keep, tail=tail)
        if ok:
            chat.update(sent=value, v=version, at=now)

    # -- requests from any window -------------------------------------------------------
    def handle(self, method: str, params: dict) -> Any:
        from delfin.agent.service import protocol
        from delfin.agent.service.host import HostRefused

        streaming = bool(self.state.get("streaming"))
        if method == "session.prompt":
            text = str(params.get("text") or "")
            asked = _KEEP_LINE.match(text)
            if asked:
                if asked.group(1):
                    self.on_keep(asked.group(1).lower() == "on")
                return {"accepted": True}
            if params.get("via") == "messenger" and text.lstrip()[:1] in ("/", "!"):
                raise HostRefused(protocol.NOT_ALLOWED,
                                  "A command from the phone is not run in a dashboard.")
            if not text.strip():
                return {"accepted": False}
            accepted = bool(self.adapter["send"](text))
            return {"accepted": accepted, **({"queued": True} if streaming and accepted else {})}
        if method == "session.cancel":
            self.adapter["stop"]()
            return {"stopping": streaming, "refused": 0}
        if method == "permission.reply":
            rid = str(params.get("requestId") or "")
            if params.get("client") and self.follower is not None:
                self.follower.repliers[rid] = str(params["client"])
            got = self.session({"method": "permission.reply", **params})
            return {"applied": bool(got.get("applied")), "message": str(got.get("error") or "")}
        if method == "session.setTitle":
            title = protocol.clean_title(params.get("title"))
            got = self.session({"method": "session.setTitle", "title": title}) if title else {}
            if got.get("error"):
                raise HostRefused(protocol.INVALID_PARAMS, str(got["error"]))
            self.title = title
            self.on_title(title)
            return {"title": title}
        if method == "session.setMode":
            mode = str(params.get("mode") or "")
            if mode not in PROFILE_OF_MODE:
                raise HostRefused(protocol.INVALID_PARAMS, f"no posture {mode!r} here")
            approval = self.adapter["set_permission"](PROFILE_OF_MODE[mode],
                                                     bool(params.get("confirmed")))
            return {"approval": MODE_OF_PROFILE.get(approval, mode)}
        if method == "session.status":
            info = self._now_info()
            return {"state": "running" if streaming else "idle",
                    "approval": str(info.get("approval") or "default"),
                    "queued": len(self.state.get("message_queue") or ())}
        if method == "session.choices":
            return self.adapter["choices"](str(params.get("what") or ""))
        if method == "session.commands":
            from delfin.dashboard.tab_agent import _BUILTIN_SLASH_PREFIXES
            return {"commands": [{"name": name, "summary": "", "category": "dashboard",
                                  "takesArgs": True, "source": "builtin"}
                                 for name in sorted(_BUILTIN_SLASH_PREFIXES)]}
        if method == "session.completePath":
            from pathlib import Path

            from delfin.agent.service.choices import complete_path
            return {"paths": complete_path(str(params.get("fragment") or ""), Path(self.cwd))}
        if method == "session.toggleTasks":
            return {"showTasks": False}
        if method == "session.shutdown":
            busy = "a turn is running" if streaming else ""
            if params.get("dryRun") or (params.get("ifIdle") and busy):
                return {"closing": False, "busy": busy}
            self.link.after_reply(lambda: threading.Thread(
                target=self.on_end, daemon=True,
                name=f"tab-host-end-{self.sid[:8]}").start())
            return {"closing": True}
        raise HostRefused(protocol.METHOD_NOT_FOUND, f"no {method} for a dashboard session")
