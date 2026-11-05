"""One agent-service session, shown in the dashboard.

The session runs in the person's agent service, not in this kernel: the
TUI and the phone can have it open at the same time, and closing this view
leaves it running. What the view shows comes from the service's events
(``live_sessions.Transcript``); what it sends goes through
``live_sessions.LiveClient``.

    header    name (pencil: rename), folder, model, state; End session…
    controls  agent mode, model, effort, permission -- each shows what the
              session has (its session.info/updated), never only what was
              clicked; a change is sent as the session's own command
              (/mode, /model, /effort) or session.setMode
    chat      the transcript, in the agent tab's own markup and classes, so
              the page's chat scripts (morph, follow-the-end, scroll
              guards) treat it exactly as the in-kernel chat
    foot      open approvals and questions with their own buttons
    input     send (queued while a turn runs), into the running turn, stop
              the turn

Rendering is coalesced: an event marks the view changed and arms one timer
(RENDER_DELAY_S); nothing polls. Finished messages are rendered once and
cached; only the message still streaming is rendered again.
"""

from __future__ import annotations

import html
import threading
import time
from typing import Any, Callable, Optional

import ipywidgets as widgets

from delfin.dashboard import live_sessions as LS

#: At most this often the chat is rendered again while a session streams
#: (the in-kernel tab's rate, about 8 per second).
RENDER_DELAY_S = 0.12

#: The effort picker's words (the in-kernel tab's labels).
_EFFORT_LABELS = {"": "Default", "low": "Low", "medium": "Medium", "high": "High",
                  "xhigh": "Extra High"}

#: How long the second click of a lasting rule, of Bypass, waits.
CONFIRM_S = 10.0

#: The permission ladder, as the in-kernel tab labels it, by the posture
#: each rung is (tab_agent.PROFILE_TO_CLI_PERM's values).
PERMISSIONS = (("Plan", "plan"), ("Ask All", "default"),
               ("Accept Edits", "acceptEdits"), ("Bypass", "bypassPermissions"))

_LIVE_CSS = """<style>
.delfin-live-head { align-items: center; gap: 6px; padding: 2px 2px 4px; }
.delfin-live-name { font-size: 14px; font-weight: 600; color: #0f172a; }
.delfin-live-where { font-size: 12px; color: #64748b; }
.delfin-live-chip { display: inline-block; padding: 0 7px; border-radius: 9px;
    font-size: 11px; line-height: 18px; background: #e2e8f0; color: #334155; }
.delfin-live-chip.needs { background: #ffedd5; color: #9a3412; font-weight: 600; }
.delfin-live-chip.working { background: #dcfce7; color: #166534; }
.delfin-live-iconbtn { width: 26px !important; height: 24px !important; padding: 0 !important;
    background: transparent !important; border: 0; color: #475569; }
.delfin-live-iconbtn:hover { background: #e2e8f0 !important; }
.delfin-live-stop { margin-left: auto !important; width: auto !important; height: 24px !important;
    font-size: 12px; color: #b91c1c !important; background: transparent !important;
    border: 1px solid #fca5a5 !important; border-radius: 6px; }
.delfin-live-confirm { gap: 6px; align-items: center; padding: 4px 8px;
    background: #fef2f2; border: 1px solid #fecaca; border-radius: 6px; font-size: 12px; }
.delfin-live-controls { gap: 6px; flex-wrap: wrap; }
.delfin-live-controls .widget-dropdown { width: auto !important; }
.delfin-live-controls select { font-size: 12px; }
.delfin-live-notice { font-size: 12px; color: #b45309; }
.delfin-live-asked { gap: 4px; padding: 6px 8px; margin: 4px 6px; background: #fff7ed;
    border: 1px solid #fdba74; border-radius: 6px; }
.delfin-live-asked pre { max-height: 180px; overflow: auto; font-size: 11px; margin: 2px 0;
    white-space: pre; line-height: 1.35; font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace; }
.delfin-live-asked .jupyter-button { width: auto !important; height: 26px !important;
    font-size: 12px; border-radius: 6px; }
.delfin-live-input textarea { font-size: 13px; min-height: 44px; }
.delfin-live-send-row { gap: 6px; }
.delfin-live-send-row .jupyter-button { width: auto !important; height: 30px !important; }
.delfin-live-resolved { font-size: 11px; color: #475569; font-style: italic; padding: 0 10px; }
</style>"""

#: Enter sends from the live input; Esc in a rename field cancels it.
#: Listeners on the document, installed once per page; they act only on
#: this view's own classes (the agent tab's keys act on theirs).
LIVE_KEYS_JS = """
(function() {
    if (window.__delfinLiveKeys) return;
    window.__delfinLiveKeys = true;
    document.addEventListener('keydown', function(e) {
        var t = e.target;
        if (!t || !t.closest) return;
        if (e.key === 'Enter' && !e.shiftKey && !e.ctrlKey && !e.metaKey
                && t.tagName === 'TEXTAREA' && t.closest('.delfin-live-input')) {
            var root = t.closest('.delfin-live-root');
            var send = root && root.querySelector('.delfin-live-send');
            if (send) { e.preventDefault(); e.stopPropagation(); send.click(); }
            return;
        }
        if (e.key === 'Escape' && t.closest('.delfin-live-rename')) {
            var box = t.closest('.delfin-live-rename');
            var cancel = box.querySelector('.delfin-live-rename-cancel');
            if (cancel) { e.preventDefault(); e.stopPropagation(); cancel.click(); }
        }
    }, true);
})();
"""


def _classed(w, *names):
    for name in names:
        w.add_class(name)
    return w


def _render_entry(T, entry: dict) -> str:
    role = entry["role"]
    text = str(entry.get("text") or "")
    if role == "user":
        body = T._md_to_html(text)
        if entry.get("source") == "harness":
            return f'<div class="delfin-chat-msg delfin-chat-system">{body}</div>'
        return f'<div class="delfin-chat-msg delfin-chat-user">{body}</div>'
    if role == "agent":
        if not text.strip() and not entry.get("streaming"):
            return ""
        body = T._streaming_html(text) if entry.get("streaming") else T._md_to_html(text)
        return ('<div class="delfin-chat-msg delfin-chat-agent">'
                '<div class="delfin-chat-role">Agent</div>' + body + '</div>')
    if role == "thinking":
        preview = html.escape(" ".join(text.split())[:80])
        return ('<details class="delfin-chat-msg delfin-chat-thinking">'
                f'<summary>Thinking: {preview}</summary>'
                f'<div class="thinking-content">{html.escape(text)}</div></details>')
    if role == "tool":
        mark = "…" if not entry.get("done") else ("✓" if entry.get("ok") else "✗")
        out = str(entry.get("output") or "")
        detail = (f"<details><summary>{html.escape(text)}</summary>"
                  f"<pre>{html.escape(out)}</pre></details>" if out.strip()
                  else html.escape(text))
        return f'<div class="delfin-chat-msg delfin-chat-tool">{mark} {detail}</div>'
    if role == "resolved":
        return f'<div class="delfin-chat-msg delfin-live-resolved">{html.escape(text)}</div>'
    if role == "error":
        return ('<div class="delfin-chat-msg delfin-chat-system" style="color:#b91c1c">'
                + T._md_to_html(text) + '</div>')
    if role == "context":
        title = html.escape(str(entry.get("title") or "context"))
        return (f'<div class="delfin-chat-msg delfin-chat-system"><b>{title}</b>'
                + T._md_to_html(text) + '</div>')
    return f'<div class="delfin-chat-msg delfin-chat-system">{T._md_to_html(text)}</div>'


class _Renderer:
    """The transcript's HTML; finished messages rendered once."""

    def __init__(self) -> None:
        from delfin.dashboard import tab_agent as T

        self.T = T
        self._cache: dict[int, tuple[tuple, str]] = {}

    def __call__(self, transcript: LS.Transcript) -> str:
        parts = []
        with transcript.lock:
            entries = list(transcript.entries)
            keys = [(e["role"], len(str(e.get("text") or "")), bool(e.get("streaming")),
                     bool(e.get("done")), len(str(e.get("output") or ""))) for e in entries]
        seen = set()
        for entry, key in zip(entries, keys):
            ident = id(entry)
            seen.add(ident)
            cached = self._cache.get(ident)
            if cached is None or cached[0] != key:
                cached = (key, _render_entry(self.T, entry))
                self._cache[ident] = cached
            parts.append(cached[1])
        for ident in list(self._cache):
            if ident not in seen:
                del self._cache[ident]
        if not parts:
            parts = ['<div class="delfin-chat-msg delfin-chat-system">'
                     'Nothing said yet in this session.</div>']
        return self.T._chat_value('<div class="delfin-agent-chat">' + "".join(parts) + "</div>")


def _state_chip(transcript: LS.Transcript) -> str:
    word = LS.state_word({"state": transcript.state, "pending": len(transcript.pending)})
    if transcript.state == "ended":
        word = "ended"
    cls = {"needs you": "needs", "working": "working"}.get(word, "")
    return f'<span class="delfin-live-chip {cls}">{html.escape(word)}</span>'


def build_live_tab(client: LS.LiveClient, sid: str, *, attached: Optional[dict] = None,
                   on_change: Optional[Callable[[], None]] = None,
                   page_scripts: bool = True) -> tuple[widgets.Widget, dict]:
    """The view of service session ``sid``: ``(widget, refs)``.

    Attaches (``attached``: an attach answer already in hand, e.g. from
    session.create) and subscribes; ``refs["close"]`` unsubscribes and
    detaches -- the session keeps running. ``on_change`` is called after
    each render (the session list's dot and name follow it)."""
    from delfin.dashboard import tab_agent as T

    transcript = LS.Transcript(sid)
    render = _Renderer()
    view: dict[str, Any] = {"timer": None, "shown": -1, "syncing": False, "closed": False,
                            "asked_keys": (), "armed": {}, "bypass_at": 0.0,
                            "lock": threading.Lock(), "unsubscribe": None}

    css = widgets.HTML(T._AGENT_CSS + _LIVE_CSS)
    name_html = _classed(widgets.HTML(""), "delfin-live-name")
    rename_btn = _classed(widgets.Button(description="✎", tooltip="Rename this session"),
                          "delfin-live-iconbtn")
    rename_input = widgets.Text(placeholder="Name (empty: folder and short id)",
                                continuous_update=False, layout=widgets.Layout(width="260px"))
    rename_ok = widgets.Button(description="Save", layout=widgets.Layout(width="auto"))
    rename_cancel = _classed(widgets.Button(description="Cancel", layout=widgets.Layout(width="auto")),
                             "delfin-live-rename-cancel")
    rename_box = _classed(widgets.HBox([rename_input, rename_ok, rename_cancel],
                                       layout=widgets.Layout(display="none")), "delfin-live-rename")
    where_html = _classed(widgets.HTML(""), "delfin-live-where")
    stop_btn = _classed(widgets.Button(description="End session…", tooltip=(
        "End this session for every window (dashboard, terminal, phone). It is "
        "saved and can be resumed. Closing this view does not end it.")), "delfin-live-stop")
    stop_yes = widgets.Button(description="End it", button_style="danger",
                              layout=widgets.Layout(width="auto"))
    stop_no = widgets.Button(description="Keep it running", layout=widgets.Layout(width="auto"))
    stop_confirm = _classed(widgets.HBox([
        widgets.HTML("End this session for every window? It is saved and can be resumed."),
        stop_yes, stop_no], layout=widgets.Layout(display="none")), "delfin-live-confirm")
    head = _classed(widgets.HBox([name_html, rename_btn, rename_box, where_html, stop_btn]),
                    "delfin-live-head")

    mode_dd = widgets.Dropdown(description="Mode", options=[("…", "")], value="",
                               layout=widgets.Layout(width="170px"))
    model_dd = widgets.Dropdown(description="Model", options=[("…", "")], value="",
                                layout=widgets.Layout(width="300px"))
    effort_dd = widgets.Dropdown(description="Effort", options=[
        ("Default", ""), ("Low", "low"), ("Medium", "medium"), ("High", "high"),
        ("Extra High", "xhigh")], value="", layout=widgets.Layout(width="170px"))
    perm_dd = widgets.Dropdown(description="Permissions", options=list(PERMISSIONS),
                               value="default", layout=widgets.Layout(width="230px"))
    provider_html = _classed(widgets.HTML(""), "delfin-live-where")
    for dd in (mode_dd, model_dd, effort_dd, perm_dd):
        dd.style.description_width = "initial"
    controls = _classed(widgets.HBox([mode_dd, model_dd, effort_dd, perm_dd, provider_html]),
                        "delfin-live-controls")
    notice_html = _classed(widgets.HTML(""), "delfin-live-notice")

    chat = _classed(widgets.HTML(""), "delfin-agent-chat-host")
    approval_box = widgets.VBox([])
    wait_html = _classed(widgets.HTML(""), "delfin-live-where")
    frame = _classed(widgets.VBox([chat, approval_box, wait_html]), "delfin-agent-chat-frame")

    input_box = _classed(widgets.Textarea(placeholder=(
        "Message this session. Enter sends; while it works the message waits "
        "for the turn to end."), layout=widgets.Layout(width="100%")), "delfin-live-input")
    send_btn = _classed(widgets.Button(description="Send", button_style="primary"),
                        "delfin-live-send")
    steer_btn = widgets.Button(description="Send into the turn", tooltip=(
        "The model reads it at its next step of the running turn"),
        layout=widgets.Layout(display="none", width="auto"))
    cancel_btn = widgets.Button(description="Stop", tooltip="End the running turn (Esc)",
                                layout=widgets.Layout(display="none", width="auto"))
    send_row = _classed(widgets.HBox([send_btn, steer_btn, cancel_btn]), "delfin-live-send-row")

    scripts = widgets.Output()
    if page_scripts:
        with scripts:
            from IPython.display import Javascript, display
            display(Javascript(T._VISIBLE_LOOKUP_JS + T._AGENT_PAGE_JS + LIVE_KEYS_JS
                               + T._MEASURE_AGENT_TOP_JS))

    root = _classed(widgets.VBox([css, scripts, head, stop_confirm, controls, notice_html, frame,
                                  input_box, send_row]),
                    "delfin-agent-root", "delfin-live-root")

    # -- saying things ---------------------------------------------------------------
    def say(text: str) -> None:
        notice_html.value = html.escape(text) if text else ""

    def call(fn: Callable[[], Any], what: str) -> Any:
        try:
            return fn()
        except LS.ServiceError as exc:
            say(f"{what}: {exc.message}")
        except TimeoutError:
            say(f"{what}: the agent service did not answer in time")
        except OSError as exc:
            say(f"{what}: the agent service is not reachable ({exc})")
        return None

    # -- what is on screen -------------------------------------------------------------
    def _sync_controls() -> None:
        info = transcript.info
        view["syncing"] = True
        try:
            for dd, key in ((mode_dd, "mode"), (model_dd, "model"), (effort_dd, "effort")):
                value = str(info.get(key) or "")
                values = [v for _, v in dd.options]
                if value not in values:
                    dd.options = [*dd.options, (value or "Default", value)]
                if dd.value != value:
                    dd.value = value
            approval = str(info.get("approval") or "default")
            if approval in {v for _, v in PERMISSIONS} and perm_dd.value != approval:
                perm_dd.value = approval
            provider = str(info.get("provider") or "")
            provider_html.value = (f"provider {html.escape(provider)} (set when it started)"
                                   if provider else "")
        finally:
            view["syncing"] = False

    def _sync_head() -> None:
        name_html.value = html.escape(transcript.name)
        info = transcript.info
        bits = [html.escape(LS._folder(str(info.get("workspace") or ""))),
                html.escape(str(info.get("model") or "")), _state_chip(transcript)]
        background = [r for r in (info.get("background") or []) if isinstance(r, dict)]
        if background:
            agents = [r for r in background if str(r.get("kind") or "").startswith("agent")
                      or r.get("group") == "agents"]
            label = (f"{len(agents)} sub-agent{'s' * (len(agents) != 1)}, "
                     if agents else "")
            bits.append(f"{label}{len(background)} in the background")
        if transcript.ended:
            bits.append(f"ended: {html.escape(transcript.ended)}")
        where_html.value = " · ".join(b for b in bits if b)
        running = transcript.state == "running"
        cancel_btn.layout.display = "" if running else "none"
        steer_btn.layout.display = "" if running else "none"
        ended = transcript.state == "ended"
        send_btn.disabled = ended
        stop_btn.disabled = ended
        wait_html.value = html.escape(transcript.wait) if transcript.wait else ""

    def _asked_widget(asked: dict) -> widgets.Widget:
        rid = str(asked.get("requestId") or "")
        kind = str(asked.get("kind") or "confirm")
        head_bits = []
        if kind == "ask":
            head_bits.append(f"<b>{html.escape(str(asked.get('header') or 'Question'))}</b>: "
                             f"{html.escape(str(asked.get('question') or ''))}")
        else:
            what = "The plan" if kind == "plan" else f"{asked.get('tool') or 'A tool'} asks to run"
            head_bits.append(f"<b>{html.escape(what)}</b>")
            if asked.get("purpose"):
                head_bits.append(html.escape(str(asked["purpose"])))
            preview = str(asked.get("preview") or "")
            if preview:
                head_bits.append(f"<pre>{html.escape(preview[:6000])}</pre>")
        text = widgets.HTML("<br>".join(head_bits))
        if kind == "ask":
            labels = [str(c) for c in asked.get("choices") or []]
            details = [str(d) for d in asked.get("choiceDetails") or []]
            picks = []
            for i, label in enumerate(labels):
                b = widgets.Button(description=label[:60], tooltip=f"answer: {details[i] if i < len(details) else label}")
                b.on_click(lambda _b, _l=label: call(
                    lambda: client.reply(sid, rid, answers=[_l]), "Not answered"))
                picks.append(b)
            typed = widgets.Text(placeholder="or type the answer")
            send = widgets.Button(description="Answer")
            send.on_click(lambda _b: typed.value.strip() and call(
                lambda: client.reply(sid, rid, text=typed.value.strip()), "Not answered"))
            dismiss = widgets.Button(description="Dismiss")
            dismiss.on_click(lambda _b: call(lambda: client.reply(sid, rid, answers=[]),
                                             "Not dismissed"))
            return _classed(widgets.VBox([text, widgets.HBox(picks), widgets.HBox([typed, send, dismiss])]),
                            "delfin-live-asked")
        buttons = []
        for option in asked.get("options") or []:
            key = str(option.get("key") or "")
            label = str(option.get("label") or key)
            b = widgets.Button(description=label[:48], tooltip=f"{key}: {option.get('detail') or label}")

            def on_click(_b, _key=key, _btn=b, _label=label):
                lasting = _key in ("A", "k")
                if lasting:
                    armed = view["armed"].get((rid, _key), 0.0)
                    if time.monotonic() - armed > CONFIRM_S:
                        view["armed"][(rid, _key)] = time.monotonic()
                        _btn.description = "Click again: save this rule"
                        return
                call(lambda: client.reply(sid, rid, choice=_key, persist=lasting or None),
                     "Not answered")
            b.on_click(on_click)
            buttons.append(b)
        return _classed(widgets.VBox([text, widgets.HBox(buttons)]), "delfin-live-asked")

    def _sync_approvals() -> None:
        with transcript.lock:
            pending = list(transcript.pending.values())
        keys = tuple(str(a.get("requestId")) for a in pending)
        if keys == view["asked_keys"]:
            return
        view["asked_keys"] = keys
        approval_box.children = tuple(_asked_widget(a) for a in pending)

    def render_now() -> None:
        with view["lock"]:
            view["timer"] = None
            version = transcript.version
            if version != view["shown"]:
                view["shown"] = version
                value = render(transcript)
                if chat.value != value:
                    chat.value = value
            _sync_head()
            _sync_controls()
            _sync_approvals()
        if on_change is not None:
            try:
                on_change()
            except Exception:
                pass

    def schedule() -> None:
        with view["lock"]:
            if view["timer"] is not None or view["closed"]:
                return
            timer = threading.Timer(RENDER_DELAY_S, render_now)
            timer.daemon = True
            view["timer"] = timer
        timer.start()

    def on_event(event: dict) -> None:
        if transcript.apply(event):
            schedule()

    # -- acting ------------------------------------------------------------------------
    def send(_b=None, how: str = "queue") -> None:
        text = str(input_box.value or "").strip()
        if not text:
            return
        if call(lambda: client.prompt(sid, text, how=how), "Not sent") is not None:
            input_box.value = ""
            say("")

    def on_rename_open(_b=None) -> None:
        rename_input.value = str(transcript.info.get("title") or "")
        rename_box.layout.display = ""
        rename_btn.layout.display = "none"
        name_html.layout.display = "none"

    def on_rename_close(_b=None) -> None:
        rename_box.layout.display = "none"
        rename_btn.layout.display = ""
        name_html.layout.display = ""

    def on_rename_save(_b=None) -> None:
        got = call(lambda: client.rename(sid, rename_input.value), "Not renamed")
        if got is not None:
            on_rename_close()

    def on_rename_enter(change) -> None:
        # continuous_update=False: the value arrives on Enter (or on blur).
        if rename_box.layout.display != "none" and change.get("new") != change.get("old"):
            on_rename_save()

    def on_stop(_b=None) -> None:
        stop_confirm.layout.display = ""

    def on_stop_no(_b=None) -> None:
        stop_confirm.layout.display = "none"

    def on_stop_yes(_b=None) -> None:
        stop_confirm.layout.display = "none"
        call(lambda: client.stop(sid), "Not ended")

    def on_perm(change) -> None:
        if view["syncing"]:
            return
        mode = str(change.get("new") or "")
        confirmed = False
        if mode == "bypassPermissions":
            if time.monotonic() - view["bypass_at"] > CONFIRM_S:
                view["bypass_at"] = time.monotonic()
                say("Bypass asks nothing for work inside the allowed folders. "
                    "Choose Bypass again to switch to it.")
                _sync_controls()
                return
            confirmed = True
        if call(lambda: client.posture(sid, mode, confirmed=confirmed), "Not switched") is not None:
            say("")
        _sync_controls()

    def command_control(word: str):
        def on_change_(change) -> None:
            if view["syncing"]:
                return
            value = str(change.get("new") or "")
            if not value:
                _sync_controls()
                return
            answer = call(lambda: client.prompt(sid, f"/{word} {value}"), f"/{word} not sent")
            if answer and answer.get("queued"):
                say(f"/{word} {value} waits for the running turn to end")
        return on_change_

    send_btn.on_click(send)
    steer_btn.on_click(lambda _b: send(how="steer"))
    cancel_btn.on_click(lambda _b: call(lambda: client.cancel(sid), "Not stopped"))
    rename_btn.on_click(on_rename_open)
    rename_ok.on_click(on_rename_save)
    rename_cancel.on_click(on_rename_close)
    rename_input.observe(on_rename_enter, names="value")
    stop_btn.on_click(on_stop)
    stop_no.on_click(on_stop_no)
    stop_yes.on_click(on_stop_yes)
    perm_dd.observe(on_perm, names="value")
    mode_dd.observe(command_control("mode"), names="value")
    model_dd.observe(command_control("model"), names="value")
    effort_dd.observe(command_control("effort"), names="value")

    def load_choices() -> None:
        """The pickers' values, from the session (a model list may ask the
        provider: off the kernel's main thread)."""
        for dd, what in ((mode_dd, "mode"), (model_dd, "model"), (effort_dd, "effort")):
            try:
                got = client.choices(sid, what)
            except Exception:
                continue
            labels = _EFFORT_LABELS if what == "effort" else {}
            options = [(labels.get(str(c.get("value") or ""), str(c.get("value") or "")),
                        str(c.get("value") or "")) for c in got.get("choices") or []
                       if c.get("value") or what == "effort"]
            # The session's own value first (its info), then what the
            # picker says is current.
            current = str(transcript.info.get(what) or got.get("current") or "")
            if current and current not in {v for _, v in options}:
                options.insert(0, (current, current))
            if not options:
                continue
            view["syncing"] = True
            try:
                dd.options = options
                if current:
                    dd.value = current
            finally:
                view["syncing"] = False

    # -- attach ------------------------------------------------------------------------
    def attach() -> bool:
        transcript._early = []
        view["unsubscribe"] = client.subscribe(on_event)
        got = attached if attached is not None else call(lambda: client.attach(sid), "Not opened")
        if got is None:
            return False
        transcript.load(got)
        render_now()
        threading.Thread(target=load_choices, daemon=True, name="live-choices").start()
        return True

    def close() -> None:
        """This view goes; the session keeps running in the service."""
        view["closed"] = True
        unsubscribe = view.get("unsubscribe")
        if unsubscribe is not None:
            unsubscribe()
        timer = view.get("timer")
        if timer is not None:
            timer.cancel()
        try:
            client.detach(sid)
        except Exception:
            pass

    opened = attach()
    refs = {"transcript": transcript, "send": send, "input": input_box, "send_btn": send_btn,
            "steer_btn": steer_btn, "cancel_btn": cancel_btn, "approval_box": approval_box,
            "name_html": name_html, "rename_btn": rename_btn, "rename_box": rename_box,
            "rename_input": rename_input, "rename_ok": rename_ok, "rename_cancel": rename_cancel,
            "stop_btn": stop_btn, "stop_confirm": stop_confirm, "stop_yes": stop_yes,
            "stop_no": stop_no, "perm": perm_dd, "mode": mode_dd, "model": model_dd,
            "effort": effort_dd, "chat": chat, "notice": notice_html, "close": close,
            "render_now": render_now, "opened": opened, "say": say,
            "name": lambda: transcript.name, "busy": lambda: transcript.busy,
            "needs": lambda: transcript.needs_you, "sid": sid}
    return root, refs
