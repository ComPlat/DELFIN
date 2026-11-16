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
import json
import threading
import time
from typing import Any, Callable, Optional

import ipywidgets as widgets

from delfin.dashboard import live_sessions as LS

#: At most this often the chat is rendered again while a session streams
#: (the in-kernel tab's rate, about 8 per second).
RENDER_DELAY_S = 0.12

#: Agent modes a service session cannot run: they need this kernel's page.
_KERNEL_ONLY_MODES = frozenset({"dashboard"})

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
.delfin-live-transit-line { display: flex; align-items: center; gap: 8px;
  margin: 6px 0 2px; padding: 6px 10px; border-radius: 8px; font-size: 13px;
  background: #eff6ff; color: #1e3a8a; border: 1px solid #bfdbfe; }
.delfin-live-spin { width: 12px; height: 12px; flex: 0 0 auto; border-radius: 50%;
  border: 2px solid #93c5fd; border-top-color: #1d4ed8;
  animation: delfin-live-spin 0.9s linear infinite; }
@keyframes delfin-live-spin { to { transform: rotate(360deg); } }
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
.delfin-live-asked, .delfin-live-change { gap: 4px; padding: 6px 8px; margin: 4px 6px; background: #fff7ed;
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
        badge = ('<span class="delfin-chat-via" style="font-size:11px;opacity:.7">via phone</span>'
                 if entry.get("via") == "phone" else "")
        return f'<div class="delfin-chat-msg delfin-chat-user">{badge}{body}</div>'
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
                   page_scripts: bool = True, lazy: bool = False) -> tuple[widgets.Widget, dict]:
    """The view of service session ``sid``: ``(widget, refs)``.

    ``lazy``: not attached until ``refs["show"]()`` (the session list shows
    a view when it is the one on screen and ``refs["hide"]()``s the others).

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
    # Laid out as the agent tab's full chat is (tab_agent.create_tab): the
    # header row with the mode, the provider and Stop; the chat frame; the
    # composer with the message box, the permission, model and effort chips
    # and the round Send -- same widgets in the same places, same classes,
    # so the page's own scripts and stylesheet treat it as they treat the
    # full chat (the person, 2026-11-16: "identical, not a second layout").
    # The session's name lives in the session list, as for the full chat;
    # renaming and ending it for every window stay at the end of the row.
    name_html = _classed(widgets.HTML("", layout=widgets.Layout(display="none")),
                         "delfin-live-name")
    rename_btn = _classed(widgets.Button(description="✎", tooltip="Rename this session",
                                         layout=widgets.Layout(width="auto", flex="0 0 auto")),
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
        "saved and can be resumed. Closing this view does not end it."),
        layout=widgets.Layout(width="auto", flex="0 0 auto")), "delfin-live-stop")
    stop_yes = widgets.Button(description="End it", button_style="danger",
                              layout=widgets.Layout(width="auto"))
    stop_no = widgets.Button(description="Keep it running", layout=widgets.Layout(width="auto"))
    stop_confirm = _classed(widgets.HBox([
        widgets.HTML("End this session for every window? It is saved and can be resumed."),
        stop_yes, stop_no], layout=widgets.Layout(display="none")), "delfin-live-confirm")

    mode_dd = widgets.Dropdown(description="Mode:", options=[("…", "")], value="",
                               layout=widgets.Layout(width="160px"),
                               style={"description_width": "45px"})
    provider_dd = widgets.Dropdown(description="Provider:", options=[("…", "")], value="",
                                   disabled=True, layout=widgets.Layout(width="160px"),
                                   style={"description_width": "55px"})
    provider_html = _classed(widgets.HTML(""), "delfin-live-where")
    # Stop ends the running turn (Esc), where the full chat has it.
    cancel_btn = widgets.Button(description="Stop", tooltip="End the running turn (Esc)",
                                layout=widgets.Layout(width="auto", flex="0 0 auto"))
    cancel_btn.disabled = True
    # Auto mode, as the full chat's toggle: the session's own /auto.
    auto_btn = widgets.Button(description="Auto: off", layout=widgets.Layout(width="90px"),
                              tooltip=("Auto mode: when nobody types, start the next item of "
                                       "work in scope by itself. Never pushes, opens pull "
                                       "requests or merges; never widens the posture."))
    # The report beetle: a report is sent from the chat that runs the
    # session, so here it says where (the full chat takes over in a moment).
    bug_btn = widgets.Button(description="🐞", layout=widgets.Layout(width="44px", flex="0 0 auto"),
                             tooltip=("A bug report is sent from the chat that runs this "
                                      "session -- here once it has moved into this dashboard."))
    controls = _classed(widgets.HBox([mode_dd, provider_dd, auto_btn, cancel_btn, stop_btn,
                                      rename_btn, rename_box, bug_btn, name_html],
                                     layout=widgets.Layout(flex_flow="row wrap")),
                        "delfin-live-controls")
    head = _classed(widgets.VBox([controls], layout=widgets.Layout(margin="0 0 6px 0")),
                    "delfin-live-head")
    notice_html = _classed(widgets.HTML(""), "delfin-live-notice")
    # Said while the session is on its way into this dashboard's full chat
    # (agent_sessions sets it): what is happening and that typing works --
    # so a wait never reads as a hang.
    transit_html = _classed(widgets.HTML(""), "delfin-live-transit")

    model_dd = widgets.Dropdown(description="Model:", options=[("…", "")], value="",
                                layout=widgets.Layout(width="auto"),
                                style={"description_width": "45px"})
    effort_dd = widgets.Dropdown(description="Effort:", options=[
        ("Default", ""), ("Low", "low"), ("Medium", "medium"), ("High", "high"),
        ("Extra High", "xhigh")], value="", layout=widgets.Layout(width="auto"),
        style={"description_width": "42px"})
    perm_dd = widgets.Dropdown(description="Perms:", options=list(PERMISSIONS),
                               value="default", layout=widgets.Layout(width="auto"),
                               style={"description_width": "48px"})
    for dd, kind in ((perm_dd, "perm"), (model_dd, "model"), (effort_dd, "effort")):
        dd.add_class("delfin-agent-chip")
        dd.add_class(f"delfin-agent-chip-{kind}")

    chat = _classed(widgets.HTML(""), "delfin-agent-chat-host")
    approval_box = widgets.VBox([])
    # "change the command" (key c): the asked command, to edit and run once.
    change_input = widgets.Text(layout=widgets.Layout(width="100%"))
    change_run = widgets.Button(description="Run this instead", button_style="warning",
                                layout=widgets.Layout(width="auto"))
    change_cancel = widgets.Button(description="Cancel", layout=widgets.Layout(width="auto"))
    change_box = _classed(widgets.HBox([change_input, change_run, change_cancel],
                                       layout=widgets.Layout(display="none")), "delfin-live-change")
    wait_html = _classed(widgets.HTML(""), "delfin-live-where")
    frame = _classed(widgets.VBox([chat, transit_html, approval_box, change_box, wait_html]),
                     "delfin-agent-chat-frame")

    input_box = _classed(widgets.Textarea(placeholder=(
        "Message this session. Enter sends; while it works the message waits "
        "for the turn to end."), layout=widgets.Layout(width="100%", height="56px")),
        "delfin-agent-input")
    send_btn = _classed(widgets.Button(description="Send", button_style="primary",
                                       tooltip="Send (Enter)",
                                       layout=widgets.Layout(width="34px", height="34px")),
                        "delfin-agent-send", "delfin-live-send")
    # Enter sends what is in the box now, through this channel (the page
    # script __delfinAgentSubmit, as in the full chat).
    submit_channel = _classed(widgets.Text(value="", layout=widgets.Layout(display="none")),
                              "delfin-agent-submit")
    steer_btn = widgets.Button(description="Send into the turn", tooltip=(
        "The model reads it at its next step of the running turn"),
        layout=widgets.Layout(display="none", width="auto"))
    send_row = _classed(widgets.HBox([send_btn, submit_channel]),
                        "delfin-agent-send-row", "delfin-live-send-row")
    gap = _classed(widgets.HTML(""), "delfin-agent-composer-gap")
    composer_bar = _classed(widgets.HBox([perm_dd, model_dd, effort_dd, gap, steer_btn,
                                          send_row]), "delfin-agent-composer-bar")
    composer = _classed(widgets.VBox([input_box, composer_bar]), "delfin-agent-composer")
    input_area = _classed(widgets.VBox([composer]), "delfin-agent-input-area")

    scripts = widgets.Output()
    if page_scripts:
        with scripts:
            from IPython.display import Javascript, display
            display(Javascript(T._VISIBLE_LOOKUP_JS + T._AGENT_PAGE_JS + LIVE_KEYS_JS
                               + T._MEASURE_AGENT_TOP_JS))

    # The folded "Details" at the foot, as in the full chat: where the
    # session works, its model and state, what runs in the background.
    details_body = _classed(widgets.VBox([where_html, provider_html]),
                            "delfin-agent-below-body")
    details = _classed(widgets.Accordion(children=[details_body], titles=("Details",),
                                         selected_index=None), "delfin-agent-below")
    root = _classed(widgets.VBox([css, scripts, head, stop_confirm, notice_html, frame,
                                  input_area, details]),
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
            if provider and provider_dd.value != provider:
                provider_dd.options = [(provider, provider)]
                provider_dd.value = provider
            provider_dd.tooltip = "Set when the session started"
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
        auto_line = str(info.get("auto") or "").strip()
        auto_on = bool(auto_line) and not auto_line.lower().startswith(("off", "auto: off"))
        label = "Auto: on" if auto_on else "Auto: off"
        if auto_btn.description != label:
            auto_btn.description = label
        auto_btn.button_style = "info" if auto_on else ""
        running = transcript.state == "running"
        cancel_btn.disabled = not running
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
            question = str(asked.get("question") or "")
            with transcript.lock:
                said = next((e["text"] for e in reversed(transcript.entries)
                             if e["role"] == "agent"), "")
            if LS.ends_with_question(said, question):
                # The session's last answer ends with it: not said twice.
                question = "the question above"
            head_bits.append(f"<b>{html.escape(str(asked.get('header') or 'Question'))}</b>: "
                             f"{html.escape(question)}")
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
            from delfin.agent.approval_answers import answer_hint
            # The real range, how to pick several, free text allowed: the
            # same line every other surface shows.
            hint = widgets.HTML(f"<small>{html.escape(answer_hint(len(labels), bool(asked.get('multiSelect'))))}</small>")
            typed = widgets.Text(placeholder="or type the answer")
            send = widgets.Button(description="Answer")
            send.on_click(lambda _b: typed.value.strip() and call(
                lambda: client.reply(sid, rid, text=typed.value.strip()), "Not answered"))
            dismiss = widgets.Button(description="Dismiss")
            dismiss.on_click(lambda _b: call(lambda: client.reply(sid, rid, answers=[]),
                                             "Not dismissed"))
            return _classed(widgets.VBox([text, widgets.HBox(picks), hint, widgets.HBox([typed, send, dismiss])]),
                            "delfin-live-asked")
        buttons = []
        for option in asked.get("options") or []:
            key = str(option.get("key") or "")
            label = str(option.get("label") or key)
            if key == "c" and not asked.get("command"):
                # Nothing to edit (the command holds a secret the preview
                # masks): the key is not offered here.
                continue
            b = widgets.Button(description=label[:48], tooltip=f"{key}: {option.get('detail') or label}")

            def on_click(_b, _key=key, _btn=b, _label=label):
                if _key == "c":
                    # The command to change, in the edit box under the cards;
                    # "Run this instead" sends it (approved once, as edited).
                    view["change_rid"] = rid
                    change_input.value = str(asked.get("command") or "")
                    change_box.layout.display = ""
                    return
                lasting = _key in ("A", "k")
                if lasting:
                    armed = view["armed"].get((rid, _key), 0.0)
                    window = CONFIRM_S
                    if time.monotonic() - armed > window:
                        view["armed"][(rid, _key)] = time.monotonic()
                        _btn.description = "Click again: save this rule"

                        def disarm(_btn=_btn, _label=_label, _k=(rid, _key)):
                            # Past its window it says again what it does.
                            if time.monotonic() - view["armed"].get(_k, 0.0) >= window:
                                view["armed"].pop(_k, None)
                                _btn.description = _label[:48]
                        timer = threading.Timer(window, disarm)
                        timer.daemon = True
                        timer.start()
                        return
                    view["armed"].pop((rid, _key), None)
                got = call(lambda: client.reply(sid, rid, choice=_key, persist=lasting or None),
                           "Not answered")
                if got and not got.get("applied") and got.get("message"):
                    say(str(got["message"]))        # e.g. "?": what each key does
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
        if view.get("change_rid") not in keys:
            view["change_rid"] = ""
            change_box.layout.display = "none"
        approval_box.children = tuple(_asked_widget(a) for a in pending)

    def render_now() -> None:
        with view["lock"]:
            view["timer"] = None
            version = transcript.version
            if version != view["shown"]:
                view["shown"] = version
                with transcript.lock:
                    mirrored = transcript.chat_html
                # A session another dashboard hosts: that dashboard's chat,
                # as it shows it -- same markup, same page scripts.
                value = mirrored or render(transcript)
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

    def on_rename_close(_b=None) -> None:
        rename_box.layout.display = "none"
        rename_btn.layout.display = ""

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

    def on_change_run(_b=None) -> None:
        rid = str(view.get("change_rid") or "")
        command = str(change_input.value or "").strip()
        if not rid or not command:
            return
        if call(lambda: client.reply(sid, rid, choice="c", command=command),
                "Not changed") is not None:
            view["change_rid"] = ""
            change_box.layout.display = "none"

    def on_change_cancel(_b=None) -> None:
        view["change_rid"] = ""
        change_box.layout.display = "none"

    def on_auto(_b=None) -> None:
        on = auto_btn.description == "Auto: on"
        answer = call(lambda: client.prompt(sid, "/auto off" if on else "/auto on"),
                      "Auto mode not switched")
        if answer and answer.get("queued"):
            say("/auto waits for the running turn to end")

    def on_bug(_b=None) -> None:
        say("A bug report is sent from the chat that runs this session; "
            "it is here once the session has moved into this dashboard.")

    auto_btn.on_click(on_auto)
    bug_btn.on_click(on_bug)
    submitted = {"seq": 0}

    def on_submit_channel(change) -> None:
        try:
            data = json.loads(change.get("new") or "{}")
        except (TypeError, ValueError):
            return
        if not isinstance(data, dict):
            return
        seq = int(data.get("seq") or 0)
        if seq <= submitted["seq"]:
            return
        submitted["seq"] = seq
        text = str(data.get("text") or "")
        if text.strip():
            input_box.value = text
            send()

    submit_channel.observe(on_submit_channel, names="value")

    def transit(text: str = "") -> None:
        """The hand-over line: empty hides it."""
        value = (f'<div class="delfin-live-transit-line"><span class="delfin-live-spin">'
                 f'</span>{html.escape(text)}</div>') if text else ""
        if transit_html.value != value:
            transit_html.value = value

    change_run.on_click(on_change_run)
    change_cancel.on_click(on_change_cancel)
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
            if what == "mode":
                # Dashboard mode drives a dashboard's own tabs, which a
                # session in the agent service does not have: not offered.
                options = [o for o in options if o[1] not in _KERNEL_ONLY_MODES]
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
    # Attached only while shown (show/hide): the service lets a session
    # nobody watches sleep, and a reload must not wake every remembered one.
    def attach() -> bool:
        if view.get("attached"):
            return True
        transcript._early = []
        view["unsubscribe"] = client.subscribe(on_event)
        got = view.pop("first", None) or call(lambda: client.attach(sid), "Not opened")
        if got is None:
            unsubscribe = view.pop("unsubscribe", None)
            if unsubscribe is not None:
                unsubscribe()
            return False
        view["attached"] = True
        transcript.load(got)
        render_now()
        threading.Thread(target=load_choices, daemon=True, name="live-choices").start()
        return True

    def hide() -> None:
        """Off screen: detached, so its session may sleep; shown again, it
        is attached again (its recent turns, then live)."""
        if not view.get("attached"):
            return
        view["attached"] = False
        unsubscribe = view.pop("unsubscribe", None)
        if unsubscribe is not None:
            unsubscribe()
        timer = view.get("timer")
        if timer is not None:
            timer.cancel()
            view["timer"] = None
        try:
            client.detach(sid)
        except Exception:
            pass

    def close() -> None:
        """This view goes; the session keeps running in the service."""
        hide()
        view["closed"] = True

    view["first"] = attached
    opened = attach() if not lazy else True
    refs = {"transcript": transcript, "send": send, "input": input_box, "send_btn": send_btn,
            "steer_btn": steer_btn, "cancel_btn": cancel_btn, "approval_box": approval_box,
            "name_html": name_html, "rename_btn": rename_btn, "rename_box": rename_box,
            "rename_input": rename_input, "rename_ok": rename_ok, "rename_cancel": rename_cancel,
            "stop_btn": stop_btn, "stop_confirm": stop_confirm, "stop_yes": stop_yes,
            "stop_no": stop_no, "perm": perm_dd, "mode": mode_dd, "model": model_dd,
            "effort": effort_dd, "chat": chat, "notice": notice_html, "close": close,
            "render_now": render_now, "opened": opened, "say": say,
            "change_box": change_box, "change_input": change_input, "change_run": change_run,
            "change_cancel": change_cancel, "submit": submit_channel,
            "transit": transit, "transit_html": transit_html,
            "draft": lambda: str(input_box.value or ""),
            "provider": provider_dd, "auto_btn": auto_btn, "bug_btn": bug_btn,
            "details": details,
            "show": attach, "hide": hide, "attached": lambda: bool(view.get("attached")),
            "name": lambda: transcript.name, "busy": lambda: transcript.busy,
            "needs": lambda: transcript.needs_you, "sid": sid,
            # The control row: the session list puts its Keep session
            # switch here, as it does in the in-kernel tab's row.
            "controls": controls}
    return root, refs
