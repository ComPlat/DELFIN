"""The DELFIN Agent tab: sessions listed on the left, the open one on the right.

Each session is a whole agent tab (``tab_agent.create_tab``) with its own
engine, settings, chat and bug report. They run side by side and the list
shows one at a time. A session picks its working directory when it starts
(where delfin-voila was launched, by default), and the sessions that were
open come back when the dashboard starts again.

DELFIN.md and memory are files every session reads; the conversation is
each session's own.
"""

from __future__ import annotations

import html
import json
import os
import re
import subprocess
import threading
import uuid
from pathlib import Path
from typing import Any, Callable, Optional

import ipywidgets as widgets

_OPEN_SESSIONS_PATH = Path.home() / ".delfin" / "agent_open_sessions.json"
# Each open session is a full agent tab with its own engine and timers.
_MAX_OPEN = 8
_REFRESH_S = 3.0
# Reading the saved sessions parses every session file, so the resume list
# is rebuilt on open/close and at this interval, not on every refresh.
_RESUME_REFRESH_S = 60.0
# How long the first click on "Stop all agents" waits for the second.
_STOP_ARM_S = 8.0

#: Ids of contexts whose splitter init script is already registered. The
#: dashboard builds the agent tab once, but a settings change can rebuild a
#: tab on the same context; without the guard the drag script would be
#: emitted twice into the one keep_js() payload.
_REGISTERED_INIT_JS = set()


def _register_init_js(ctx: Any, script: str) -> None:
    """Add ``script`` to a context's init-js payload, at most once."""
    add = getattr(ctx, "add_init_js", None)
    if not callable(add):
        return
    key = id(ctx)
    if key in _REGISTERED_INIT_JS:
        return
    _REGISTERED_INIT_JS.add(key)
    try:
        add(script)
    except Exception:
        _REGISTERED_INIT_JS.discard(key)


# The session list, in the agent tab's own palette (tab_agent._AGENT_CSS):
# system font, slate greys, the chat's light blue for the session on screen.
# Everything is border-box and the column clips horizontally -- widgets set
# to 100% width plus their padding made it scroll sideways.
_SIDEBAR_CSS = """<style>
.delfin-session-shell { width: 100%; align-items: flex-start; overflow-x: clip; }
.delfin-sessions {
    /* Stays on screen while a long chat scrolls: sticky against the page.
       The shell clips with `clip`, not `hidden` -- `hidden` makes the
       shell a scroll box and a sticky child sticks to that, not the page. */
    position: sticky; top: 8px; align-self: flex-start;
    /* The drag writes --delfin-sessions-w on the shell and this rule
       reads it. It cannot write `flex` directly: the declaration below is
       `!important` (widget stylesheets override a plain one), and a
       stylesheet !important beats an inline declaration -- so the handler
       on this branch set sidebar.style.flex on every move and the column
       never changed width. A custom property is read by the !important
       rule itself, so the drag wins without the rule being weakened. */
    box-sizing: border-box;
    flex: 0 0 var(--delfin-sessions-w, 236px) !important;
    width: var(--delfin-sessions-w, 236px);
    margin: 0; padding: 8px 8px 10px; gap: 6px;
    background: #f8fafc; border: 1px solid #e5e7eb; border-radius: 10px;
    overflow: hidden;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
}
.delfin-sessions * { box-sizing: border-box; }
.delfin-sessions .widget-html, .delfin-sessions .widget-html-content {
    margin: 0; line-height: normal; }
.delfin-sessions .jupyter-button { margin: 0; box-shadow: none; }
.delfin-session-head { align-items: center; justify-content: space-between;
    padding: 0 2px 0 6px; }
.delfin-session-heading { font-size: 11px; font-weight: 600;
    letter-spacing: 0.06em; text-transform: uppercase; color: #6b7280; }
.delfin-session-add { width: 26px !important; height: 26px !important;
    padding: 0 !important; border: 0; border-radius: 7px;
    background: transparent !important; color: #475569; font-size: 18px;
    line-height: 26px; }
.delfin-session-add:hover { background: #e2e8f0 !important; color: #0f172a; }
.delfin-session-list { gap: 2px; max-height: calc(100vh - 320px);
    overflow-y: auto !important; overflow-x: hidden !important; }
.delfin-session-item { width: 100%; align-items: center; gap: 2px;
    padding: 5px 4px 5px 8px; border-radius: 8px; }
.delfin-session-item:hover { background: #eef2f7; }
.delfin-session-item.delfin-session-active { background: #dbeafe; }
.delfin-session-text { flex: 1 1 auto; min-width: 0; overflow: hidden; }
.delfin-session-title { width: 100% !important; height: 20px !important;
    padding: 0 !important; border: 0; background: transparent !important;
    text-align: left; font-size: 13px; font-weight: 500; color: #1f2937;
    line-height: 20px; white-space: nowrap; overflow: hidden;
    text-overflow: ellipsis; cursor: pointer; }
.delfin-session-title:focus-visible { outline: 2px solid #93c5fd;
    outline-offset: 1px; border-radius: 4px; }
.delfin-session-active .delfin-session-title { font-weight: 600; color: #0f172a; }
.delfin-session-path { font-size: 11px; color: #6b7280; line-height: 15px;
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.delfin-session-busy .delfin-session-title::before { content: "";
    display: inline-block; width: 7px; height: 7px; margin: 0 6px 1px 0;
    border-radius: 50%; background: #16a34a; vertical-align: middle;
    animation: delfin-session-pulse 1.4s ease-in-out infinite; }
@keyframes delfin-session-pulse { 50% { opacity: 0.35; } }
.delfin-session-needs .delfin-session-title::before { content: "";
    display: inline-block; width: 7px; height: 7px; margin: 0 6px 1px 0;
    border-radius: 50%; background: #f97316; vertical-align: middle; }
.delfin-session-unseen .delfin-session-title::before { content: "";
    display: inline-block; width: 7px; height: 7px; margin: 0 6px 1px 0;
    border-radius: 50%; background: #9ca3af; vertical-align: middle; }
@media (prefers-reduced-motion: reduce) {
    .delfin-session-busy .delfin-session-title::before { animation: none; } }
.delfin-session-close { flex: 0 0 22px; width: 22px !important;
    height: 22px !important; padding: 0 !important; border: 0;
    border-radius: 6px; background: transparent !important; color: #9ca3af;
    font-size: 15px; line-height: 22px; opacity: 0; }
.delfin-session-item:hover .delfin-session-close,
.delfin-session-active .delfin-session-close,
.delfin-session-close:focus-visible { opacity: 1; }
.delfin-session-close:hover { background: #cbd5e1 !important; color: #1f2937; }
.delfin-session-form { gap: 6px; padding: 8px; background: #ffffff;
    border: 1px solid #e5e7eb; border-radius: 8px; overflow: hidden; }
/* The Combobox renders as .widget-text; without the width rule it kept
   its default width, 4 px wider than the card, and the card (overflow
   auto by default) grew a sideways scrollbar (user, 2026-09-16). */
.delfin-session-form .widget-combobox, .delfin-session-form .widget-text,
.delfin-session-form .widget-checkbox {
    width: 100% !important; max-width: 100%; margin: 0; }
.delfin-session-form input[type="text"] { font-size: 12px; border-radius: 6px; }
.delfin-session-label { font-size: 11px; color: #6b7280; }
.delfin-session-hint { font-size: 11px; color: #475569; line-height: 1.35; }
.delfin-session-actions { gap: 6px; justify-content: flex-end; }
.delfin-session-actions .jupyter-button { width: auto !important;
    height: 26px !important; padding: 0 12px !important; border-radius: 6px;
    font-size: 12px; }
.delfin-session-start { background: #2563eb !important; color: #ffffff !important; }
.delfin-session-start:hover { background: #1d4ed8 !important; }
.delfin-session-cancel { background: transparent !important; color: #475569; }
.delfin-session-cancel:hover { background: #f1f5f9 !important; }
.delfin-session-notice { font-size: 11px; color: #b45309; padding: 0 6px; }
.delfin-session-resume { width: 100% !important; margin: 2px 0 0 !important; }
.delfin-session-resume select { font-size: 12px; color: #475569;
    border-radius: 6px; }
.delfin-session-stage { flex: 1 1 0% !important; min-width: 0 !important; }
.delfin-session-elsewhere { font-size: 11px; color: #475569; padding: 4px 6px 0; }
.delfin-session-stopall { width: 100% !important; margin: 8px 0 0 !important;
    background: transparent !important; color: #b91c1c !important;
    border: 1px solid #fca5a5 !important; border-radius: 6px; font-size: 12px; }
.delfin-session-stopall.delfin-armed { background: #b91c1c !important;
    color: #ffffff !important; }
/* Draggable divide between the sidebar and the chat, mirroring the calc
   browser's .calc-splitter. The shell aligns items at flex-start, so a
   height:100% strip would collapse to zero -- the splitter is sticky too,
   tracking the same 100vh column the sidebar sticks to. */
.delfin-session-shell { align-items: flex-start; }
/* The gap belongs to the SPLITTER, equally on both sides, and not to the
   sidebar. It was the sidebar's `margin-right: 12px`, which put the whole
   gap on one side: 12px between column and bar, nothing between bar and
   chat, so the handle sat flush against the chat with a wide empty strip
   behind it. One margin here is symmetric by construction -- there is no
   second number that can drift away from the first. */
.delfin-splitter-host { align-self: flex-start; display: flex;
    flex-direction: column; position: sticky; top: 8px;
    margin: 0 var(--delfin-splitter-gap, 10px); }
/* Height ends at the Send row, not at the bottom of the page. The JS
   below measures it; this is what shows until then and if the measure
   fails -- the chat's own height (tab_agent: 100vh - 460px) plus the
   composer, so the first paint is already close. */
.delfin-session-splitter { width: 8px;
    height: calc(100vh - 400px); min-height: 160px;
    cursor: col-resize; touch-action: none;
    background: linear-gradient(to right, #d6d6d6, #f2f2f2, #d6d6d6);
    border-radius: 4px; z-index: 10; pointer-events: auto; }
.delfin-session-splitter:hover { background: linear-gradient(
    to right, #b0b0b0, #e0e0e0, #b0b0b0); }
</style>"""

#: Drag the sessions column to a width. Runs once per shell (data-bound
#: guard) after the page is assembled, along every tab's init script in the
#: one keep_js() call from the dashboard. The widths it writes override the
#: 236px CSS default on the sidebar, so a drag sticks and a refresh keeps
#: the choice. Mirrors the calc browser's .calc-splitter handler.
_SPLITTER_INIT_JS = """\
(function () {
    // The agent tab's DOM may not exist when this runs. ipywidgets Tab
    // renders the SELECTED child only, and the dashboard sends every
    // tab's init script once the page is assembled -- so if the agent tab
    // is not the one in front, querySelectorAll finds no shell, the loop
    // body never executes, and nothing is ever bound. That is why the
    // handle could not be dragged even after the width was made
    // overridable: the handler was never attached in the first place.
    //
    // Every other tab with a splitter retries (tab_literature: 40 tries;
    // tab_calculations_browser: 400). This one ran once. It now retries
    // AND keeps watching, because a session opened minutes later brings
    // its own shell, and a bounded retry that has already expired would
    // leave that one dead.
    function bindAll() {
        var shells = document.querySelectorAll('.delfin-session-shell');
        var bound = 0;
        for (var i = 0; i < shells.length; i++) {
            if (bind(shells[i])) bound += 1;
        }
        return bound;
    }
    function bind(shell) {
            var bound = false;
            if (!shell || shell.dataset.delfinSplitterBound) return false;
            var splitter = shell.querySelector('.delfin-session-splitter');
            var sidebar = shell.querySelector('.delfin-sessions');
            var host = shell.querySelector('.delfin-splitter-host');
            // Not marked until the parts are really there. Marking first
            // meant a shell seen half-built was written off for good.
            if (!splitter || !sidebar || !host) return false;
            shell.dataset.delfinSplitterBound = '1';
            var MIN = 170, MAX = 420;
            // Read and write the one custom property the stylesheet's
            // !important rule consumes. Writing sidebar.style.flex cannot
            // work against that rule, and writing it anyway is what made
            // the drag look dead on the real page while a DOM test that
            // models no cascade reported it working.
            function current() {
                var v = shell.style.getPropertyValue('--delfin-sessions-w');
                var n = parseFloat(v);
                return isNaN(n) ? 236 : n;
            }
            function apply(w) {
                w = Math.max(MIN, Math.min(MAX, Math.round(w)));
                if (Math.abs(w - current()) < 1) return;
                shell.style.setProperty('--delfin-sessions-w', w + 'px');
            }
            function onMove(e) {
                var box = host.parentElement.getBoundingClientRect();
                apply(e.clientX - box.left);
            }
            function onUp(e) {
                if (e && e.pointerId != null) {
                    try { splitter.releasePointerCapture(e.pointerId); }
                    catch (err) { /* already released */ }
                }
                document.removeEventListener('pointermove', onMove);
                document.removeEventListener('pointerup', onUp);
                document.removeEventListener('pointercancel', onUp);
            }
            splitter.addEventListener('pointerdown', function (e) {
                e.preventDefault();
                try { splitter.setPointerCapture(e.pointerId); }
                catch (err) { /* document listeners carry it */ }
                document.addEventListener('pointermove', onMove);
                document.addEventListener('pointerup', onUp);
                document.addEventListener('pointercancel', onUp);
            });
            bound = true;

            // The bar reaches down to the Send row and stops there. Run to
            // the bottom of the column it looked like a page divider and
            // kept growing with the transcript; the composer is where the
            // conversation ends, so that is where the handle ends.
            //
            // Measured rather than computed from a CSS expression: the
            // composer's height depends on the textarea, which the user can
            // grow, and on whether the quick-cycle button is shown. The
            // stylesheet carries an approximation for the first paint.
            // Three guards, and each one is here because the first version
            // without it made the window jitter visibly -- reported from a
            // real session where the chat could not be read.
            //
            // The loop it closes: fit() writes the height of a CHILD of
            // `shell`, and the observer watched `shell` itself. Every write
            // could change the parent's measurement, which fired the
            // observer, which wrote again.
            var lastH = 0;
            var pending = false;
            function fit() {
                var send = shell.querySelector('.delfin-agent-send-row');
                if (!send) return;          // no composer yet: keep the CSS
                var top = shell.getBoundingClientRect().top;
                var bottom = send.getBoundingClientRect().bottom;
                var h = Math.round(bottom - top);
                if (h < 160) return;        // laid out but collapsed
                // 1. A tolerance, not equality. Equality stops a repeat of
                //    the SAME value but not an alternation between two that
                //    differ by a pixel, which is what sub-pixel layout
                //    produces and what oscillates forever.
                if (Math.abs(h - lastH) <= 2) return;
                lastH = h;
                splitter.style.height = h + 'px';
            }
            // 2. One measurement per frame. A burst of observer callbacks
            //    during a drag or a reflow collapses into a single write
            //    instead of a write per callback.
            function schedule() {
                if (pending) return;
                pending = true;
                var run = function () { pending = false; fit(); };
                if (typeof requestAnimationFrame === 'function') {
                    requestAnimationFrame(run);
                } else {
                    setTimeout(run, 16);
                }
            }
            fit();
            // The composer changes height when the textarea is dragged or a
            // long line wraps, and the window changes it too. A one-shot
            // measure was right for one layout and wrong after the first
            // resize.
            window.addEventListener('resize', schedule);
            try {
                // 3. Watch the ROW that decides the height, not the
                //    container the splitter lives in. The row's size does
                //    not depend on the splitter's, so writing the splitter
                //    cannot feed back into what is being measured.
                var watched = shell.querySelector('.delfin-agent-send-row');
                new ResizeObserver(schedule).observe(watched || shell);
            } catch (err) { /* older browser: resize alone */ }
            return bound;
    }

    // Try now, then keep trying: the tab may be rendered later, and a
    // session opened later brings a shell of its own. Bounded polling
    // first (fast, for the ordinary "tab not in front yet" case), then a
    // MutationObserver, which costs nothing while nothing changes and
    // catches a shell that appears long after the polling gave up.
    bindAll();
    var tries = 0;
    (function poll() {
        tries += 1;
        bindAll();
        if (tries < 60) setTimeout(poll, tries < 20 ? 100 : 500);
    })();
    try {
        var seen = new MutationObserver(function () { bindAll(); });
        seen.observe(document.body, {childList: true, subtree: true});
    } catch (err) { /* no observer: the polling above did what it could */ }
})();
"""



#: The row class for each dot; at most one is set at a time.
_DOT_CLASSES = {
    "needs": "delfin-session-needs",
    "busy": "delfin-session-busy",
    "unseen": "delfin-session-unseen",
}


#: Markup a title may open with that says nothing about the session.
_LABEL_STRIP = re.compile(r"^[\s>#*`\-=|)\]]+")
#: Only backticks and asterisks. NOT `_` or `~`: markdown uses them for
#: emphasis, and this domain uses them in names -- stripping `_` turned
#: the title "test_calc.py schlaegt fehl" into "testcalc.py", a file
#: that does not exist. An emphasis mark lost is cosmetic; a filename
#: altered is a wrong answer.
_LABEL_MARKS = re.compile(r"[`*]+")

#: How long a label may be before it is cut. Wide enough for a sentence
#: a person wrote, short enough that three of them fit the hint.
_LABEL_MAX = 40


def session_label(title: Any, session_id: Any = "") -> str:
    """A readable name for a session, from the text it opened with.

    Input: the stored ``title`` (the user's first message) and the
    session id. Output: a single line, at most ~40 characters plus an
    ellipsis, never empty.

    A title is a prompt, not a name. Of the 52 sessions on the machine
    this was written on, 38 were unusable as a label: multi-line,
    markdown headings, code fences, cut off mid-word. They are shown in
    the Resume dropdown and in the hint warning that another session
    already works in this repository -- so the control that exists to
    keep two sessions out of one checkout was unreadable exactly where
    it had to be read.

    Truncation alone does not help: the first 40 characters of a
    handover briefing are "# Übergabe: Sechs Arbeitszweige prüfen, z".
    The FIRST LINE, with its markup removed and cut on a word boundary,
    is a name. Falls back to a short form of the id, because a row with
    no label cannot be told from the row above it.
    """
    for raw in str(title or "").splitlines():
        line = _LABEL_MARKS.sub("", _LABEL_STRIP.sub("", raw)).strip()
        line = " ".join(line.split())
        if not line:
            continue
        if len(line) <= _LABEL_MAX:
            return line
        # Cut on a word boundary: a label that ends mid-word reads as a
        # different word, and three of them in one sentence read as noise.
        cut = line[:_LABEL_MAX].rsplit(" ", 1)[0].rstrip(" ,;:.")
        # ...unless that leaves almost nothing. A title that opens with a
        # long path ("Bau in tests/fixtures/user_project_workspace/ ein
        # kleines Modul ...") has no space inside the budget, and the
        # word boundary reduces it to "Bau in" -- which names no session.
        # Below half the budget the character cut says more.
        if len(cut) < _LABEL_MAX // 2:
            cut = line[:_LABEL_MAX].rstrip()
        return cut + "\u2026"
    sid = str(session_id or "").strip()
    return f"Session {sid[:8]}" if sid else "Untitled session"


def needs_you(state: dict) -> bool:
    """Whether a session is waiting for the user: an approval dialog, an
    open question, or a plan to accept. Never raises."""
    try:
        broker = state.get("_kit_confirm_broker")
        if broker is not None and getattr(broker, "_pending", None):
            return True
        ev = state.get("_ask_user_event")
        if ev is not None and not ev.is_set():
            return True
        return bool(state.get("_pending_plan_body"))
    except Exception:
        return False


def session_dot(*, busy: bool, needs: bool, unseen: bool) -> str:
    """The row class for one session: orange when it needs you, green while
    it works, grey when it finished while you were in another chat."""
    if needs:
        return _DOT_CLASSES["needs"]
    if busy:
        return _DOT_CLASSES["busy"]
    if unseen:
        return _DOT_CLASSES["unseen"]
    return ""

def _give_emergency_stop() -> None:
    """Give the emergency stop from a process of its own.

    The kernel this runs in ends with the stop within seconds; the stop's
    direct end of the other kernels on this machine must not end with it.
    """
    import subprocess
    import sys

    subprocess.Popen(
        [sys.executable, "-m", "delfin.agent.cli", "stop-all", "--yes",
         "--reason", "dashboard button"],
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, start_new_session=True)


class _SessionContext:
    """The dashboard context as one session sees it.

    Everything is the shared context except what belongs to the session: its
    working directory and the conversation it opens with, its engine and
    state (the shared ones are set to the session on screen, for the
    Activity tab), its status chip, and whether it adds page scripts (the
    first session does; the others' scripts are the same ones).
    """

    OWN = frozenset({
        "agent_workspace", "initial_session_id", "sessions_managed",
        "session_open_elsewhere", "agent_engine", "agent_state",
        "agent_status_html", "add_init_js", "presence_key",
    })

    def __init__(self, base: Any, **own: Any) -> None:
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_own", dict(own))

    def __getattr__(self, name: str) -> Any:
        own = object.__getattribute__(self, "_own")
        if name in own:
            return own[name]
        return getattr(object.__getattribute__(self, "_base"), name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name in self.OWN:
            object.__getattribute__(self, "_own")[name] = value
        else:
            setattr(object.__getattribute__(self, "_base"), name, value)


# ---------------------------------------------------------------------------
# What is remembered
# ---------------------------------------------------------------------------

def load_open_sessions() -> list[dict]:
    """The sessions that were open: ``[{"session_id", "workspace"}]``."""
    try:
        data = json.loads(_OPEN_SESSIONS_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    rows = data.get("sessions") if isinstance(data, dict) else None
    out: list[dict] = []
    here = _host()
    for row in rows or []:
        if not (isinstance(row, dict) and str(row.get("session_id") or "").strip()):
            continue
        # The home directory is shared by every login node. A session left
        # open on node A must not be reopened by a dashboard on node B --
        # two views would save over one conversation. Rows from before the
        # field carry no host and are restored as before.
        host = str(row.get("host") or "")
        if host and host != here:
            continue
        out.append({"session_id": str(row["session_id"]).strip(),
                    "workspace": str(row.get("workspace") or "")})
    return out[:_MAX_OPEN]


def _host() -> str:
    try:
        import socket
        return socket.gethostname()
    except Exception:
        return ""


def save_open_sessions(rows: list[dict]) -> None:
    """Remember which sessions are open, whole and with this host's name.

    Never raises. Written atomically: a torn write of this file is the
    whole list of open sessions gone on the next start."""
    try:
        from delfin.agent.state_paths import ensure_dir, write_text_atomic
        ensure_dir(_OPEN_SESSIONS_PATH.parent)
        here = _host()
        stamped = [{**row, "host": row.get("host") or here}
                   for row in rows if isinstance(row, dict)]
        write_text_atomic(_OPEN_SESSIONS_PATH,
                          json.dumps({"sessions": stamped}, indent=2) + "\n")
    except Exception:
        pass


def _saved_session(session_id: str) -> Optional[dict]:
    try:
        from delfin.agent.session_store import load_session
        return load_session(session_id)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Where a session works
# ---------------------------------------------------------------------------

def default_workspace(ctx: Any) -> str:
    """Where a new session works unless told otherwise.

    The directory delfin-voila was started in, by the rule the agent tab
    applies to it (inside a DELFIN tree, that tree), and the agent workspace
    when that is the home directory or above it.
    """
    from delfin.dashboard.tab_agent import _agent_workspace_from_launch
    repo = getattr(ctx, "repo_dir", None) or Path.cwd()
    workspace = str(_agent_workspace_from_launch(
        os.environ.get("DELFIN_LAUNCH_CWD", ""), repo))
    try:
        from delfin.agent.api_client import _is_forbidden_workspace_root
        agent_dir = getattr(ctx, "agent_dir", None)
        if agent_dir and _is_forbidden_workspace_root(workspace):
            return str(agent_dir)
    except Exception:
        pass
    return workspace


_PICKER_MAX = 200


def workspace_choices(ctx: Any, typed: str = "") -> list[str]:
    """Directories offered for a new session.

    First the four the dashboard knows -- the default, the agent workspace,
    the calculations and the DELFIN checkout -- then the folders a person
    can pick the way an editor's open dialog lets them: with nothing typed,
    the folders in the home directory; with a path typed, the folders in
    it (when it ends with a slash or is a directory) or those of its parent
    that start with what was typed. Hidden folders only when the typed
    name starts with a dot. Asked on every keystroke, so bounded and
    never raising.
    """
    out: list[str] = []
    for d in (default_workspace(ctx), getattr(ctx, "agent_dir", None),
              getattr(ctx, "calc_dir", None), getattr(ctx, "repo_dir", None)):
        text = str(d or "").strip()
        if text and text not in out:
            out.append(text)
    for text in _folders_like(typed):
        if text not in out:
            out.append(text)
        if len(out) >= _PICKER_MAX:
            break
    return out


def _folders_like(typed: str) -> list[str]:
    """The folders an editor's open dialog would show for ``typed``."""
    try:
        raw = str(typed or "").strip()
        if not raw:
            base, prefix = Path.home(), ""
        else:
            head, _sep, tail = raw.rpartition("/")
            cand = Path(raw).expanduser()
            # A trailing "." is a prefix for hidden folders, not the
            # directory itself (Path would fold it away).
            if raw.endswith("/") or (tail != "." and cand.is_dir()):
                base, prefix = cand, ""
            else:
                base, prefix = Path(head or "/").expanduser(), tail
        if not base.is_dir():
            return []
        shown = []
        for entry in sorted(base.iterdir(), key=lambda e: e.name.lower()):
            name = entry.name
            if not name.startswith(prefix):
                continue
            if name.startswith(".") and not prefix.startswith("."):
                continue
            try:
                if not entry.is_dir():
                    continue
            except OSError:
                continue
            shown.append(str(entry))
            if len(shown) >= _PICKER_MAX:
                break
        return shown
    except Exception:
        return []


def _exclude_locally(root: str, common_dir: str, pattern: str = ".delfin/") -> None:
    """Keep ``pattern`` out of ``git status`` in this clone only: an entry
    in its info/exclude, nothing committed. Skipped when already ignored."""
    try:
        probe = subprocess.run(
            ["git", "-C", root, "check-ignore", "-q", ".delfin/worktrees/probe"],
            capture_output=True, timeout=5)
        if probe.returncode == 0:
            return
        exclude = Path(common_dir) / "info" / "exclude"
        exclude.parent.mkdir(parents=True, exist_ok=True)
        text = exclude.read_text(encoding="utf-8") if exclude.exists() else ""
        if pattern not in text.splitlines():
            sep = "" if not text or text.endswith("\n") else "\n"
            exclude.write_text(f"{text}{sep}{pattern}\n", encoding="utf-8")
    except Exception:
        pass


def session_worktree(workspace: str) -> str:
    """A fresh worktree of the repository ``workspace`` is in, on a branch of
    its own, for a session that must not share a checkout with another.

    Under ``<repo>/.delfin/worktrees``, the way Claude Code keeps its own
    under ``.claude/worktrees``, and kept out of ``git status`` there.
    Returns the path at the same place inside the worktree as ``workspace``
    is inside the repository.
    """
    from delfin.agent import session_presence as _presence
    from delfin.agent import worktree as _wt
    repo = _presence.repository_of(workspace)
    if not repo["root"]:
        raise ValueError(f"{workspace} is not in a git repository")
    _exclude_locally(repo["root"], repo["common_dir"])
    info = _wt.enter_worktree(repo["root"], branch_prefix="session",
                              parent=Path(repo["root"]) / ".delfin" / "worktrees")
    # What a release will need, written INSIDE the tree so it goes away
    # with it. Nothing kept the WorktreeInfo: this function returned a path
    # string and dropped the rest on the floor, so no caller COULD release
    # the tree -- neither close_session nor the tab's shutdown had a
    # worktree step, and every "Own worktree" session left a checkout and
    # an orphan session/<hex> branch behind for good.
    #
    # pid and host are recorded so a later reader can tell a tree whose
    # session is gone from one that is still being worked in, without
    # judging another node's pid.
    _write_worktree_sidecar(info)
    try:
        inside = Path(workspace).resolve().relative_to(Path(repo["root"]))
    except ValueError:
        inside = Path(".")
    return str((info.path / inside).resolve())


#: Beside the worktree's own git metadata, not in a central registry: the
#: record and the thing it describes are then removed by the same act, so
#: there is no list that can outlive its entries.
_WORKTREE_SIDECAR = ".delfin/session_worktree.json"


def _write_worktree_sidecar(info) -> None:
    """Record repo, branch and owner for a session worktree. Never raises."""
    import json
    import os
    import socket
    try:
        path = Path(info.path) / _WORKTREE_SIDECAR
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({
            "repo_dir": str(info.repo_dir),
            "path": str(info.path),
            "branch": str(info.branch),
            "base_ref": str(info.base_ref),
            "created_at": float(info.created_at),
            "pid": os.getpid(),
            "host": (socket.gethostname() or "").strip(),
        }), encoding="utf-8")
    except Exception:
        pass


def read_worktree_sidecar(workspace) -> dict | None:
    """The sidecar for the session worktree ``workspace`` is in, or None.

    ``workspace`` may be a directory inside the worktree rather than its
    root -- session_worktree returns the path at the same place inside the
    tree as the original workspace was inside the repository -- so the
    search walks upwards. It stops at the first sidecar found: a worktree
    nested inside another would otherwise be released by its parent's
    record.
    """
    import json
    try:
        here = Path(workspace).expanduser().resolve()
    except Exception:
        return None
    for candidate in (here, *here.parents):
        side = candidate / _WORKTREE_SIDECAR
        try:
            if side.is_file():
                data = json.loads(side.read_text(encoding="utf-8"))
                return data if isinstance(data, dict) else None
        except (OSError, ValueError):
            return None
    return None


#: How many trees one sweep may remove. A pass that walks a whole
#: worktrees directory after a crash should not turn into an unbounded
#: amount of git work at startup; what is left over is picked up by the
#: next one.
_SWEEP_LIMIT = 12


def _sweep_candidates(worktrees_dir) -> list[dict]:
    """Sidecars under ``worktrees_dir`` that THIS host wrote, newest last.

    Only our own host: the sidecar records a pid, and a pid written on
    another login node names nothing here, so a tree belonging to one
    cannot be judged from this side. Each host reclaims its own.
    """
    import socket
    here = (socket.gethostname() or "").strip()
    out: list[dict] = []
    try:
        entries = sorted(Path(worktrees_dir).iterdir())
    except OSError:
        return out
    for tree in entries:
        if not tree.is_dir():
            continue
        side = read_worktree_sidecar(tree)
        if not side:
            continue
        if str(side.get("host") or "").strip() != here:
            continue
        out.append(side)
    return out


def _worktree_is_orphaned(side: dict, *, saved_workspaces=()) -> str:
    """"" if the tree in ``side`` may be reclaimed, else why it may not.

    Three questions this function owns, and it answers none that
    ``exit_worktree`` already answers -- uncommitted changes and running
    background jobs are its decision, asked afterwards by the caller, so
    there is one answer to "is this tree spare" rather than two.

      * is its owning process still alive on this host
      * is a live session working in it
      * would a saved session be reopened into it

    The third matters because a reopened session reuses its tree: the
    session list remembers the workspace, so removing it would turn a
    resume into a session that starts somewhere else without saying so.
    """
    import os
    import socket
    tree = str(side.get("path") or "")
    if not tree:
        return "the record names no path"
    here = (socket.gethostname() or "").strip()
    if str(side.get("host") or "").strip() != here:
        return "it belongs to another host"
    pid = int(side.get("pid") or 0)
    if pid > 0:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            pass
        except OSError:
            return f"its owning process ({pid}) may still be running"
        else:
            return f"its owning process ({pid}) is still running"
    holder = _live_session_in(tree)
    if holder:
        return f"a live session is working in it ({holder})"
    try:
        root = Path(tree).resolve()
    except Exception:
        return "its path cannot be resolved"
    for ws in saved_workspaces:
        try:
            candidate = Path(str(ws)).expanduser().resolve()
        except Exception:
            continue
        if candidate == root or root in candidate.parents:
            return "a saved session would be reopened into it"
    return ""


def _saved_session_workspaces() -> list[str]:
    """Workspaces of the sessions that can still be resumed."""
    try:
        from delfin.agent import session_store as _store
        return [str(row.get("workspace") or "")
                for row in (_store.list_sessions() or ())
                if row.get("workspace")]
    except Exception:
        return []


def reclaim_orphaned_worktrees(workspace, *, limit: int = _SWEEP_LIMIT) -> dict:
    """Remove session worktrees nobody is in. Returns a short report.

    ``{"released": [paths], "kept": {path: why}}``. Never raises: this
    runs while a dashboard is starting, and a directory that could not be
    tidied must not stop a session from opening.

    A session that crashes -- a killed kernel, a lost node -- never
    reaches close_session, so its tree and its session/<hex> branch stay.
    Nothing collected them: measured on this installation before the
    sidecar existed, there was no record to collect them BY.

    Five conditions, and a tree is removed only when all five hold: it
    carries a sidecar this host wrote, its owning process is gone, no live
    session is working in it, no saved session would be reopened into it,
    and -- asked by exit_worktree, not here -- it has no uncommitted
    changes and no background jobs running inside it.

    Deliberately NOT a general cleaner. It touches only directories whose
    sidecar says DELFIN created them for a session on this host, so a
    worktree a person made by hand, and anything else that happens to sit
    in that directory, is not its business.
    """
    report: dict = {"released": [], "kept": {}}
    try:
        from delfin.agent import session_presence as _presence
        repo = _presence.repository_of(str(workspace or ""))
        root = repo.get("root") or ""
        if not root:
            return report
        worktrees_dir = Path(root) / ".delfin" / "worktrees"
        if not worktrees_dir.is_dir():
            return report
        saved = _saved_session_workspaces()
        for side in _sweep_candidates(worktrees_dir)[:max(0, int(limit))]:
            tree = str(side.get("path") or "")
            why = _worktree_is_orphaned(side, saved_workspaces=saved)
            if why:
                report["kept"][tree] = why
                continue
            out = release_session_worktree(tree)
            if out.get("released"):
                report["released"].append(tree)
            else:
                report["kept"][tree] = out.get("kept") or "git kept it"
        # Administrative entries for trees whose directories are gone.
        if report["released"]:
            try:
                from delfin.agent import worktree as _wt
                _wt._run_git(Path(root), "worktree", "prune")
            except Exception:
                pass
    except Exception:
        return report
    return report


def _live_session_in(workspace) -> str:
    """A label for the live session whose workspace is ``workspace``, or "".

    Liveness is session_presence's answer, not a second one: its records
    carry pid and host and its reader already drops the stale ones. A
    record from another host is trusted as live -- a pid from another login
    node names nothing here, so it cannot be checked and must not be
    assumed dead.
    """
    try:
        from delfin.agent import session_presence as _presence
        here = Path(workspace).expanduser().resolve()
    except Exception:
        return ""
    try:
        rows = _presence.open_sessions()
    except Exception:
        return ""
    for row in rows or ():
        try:
            if Path(str(row.get("workspace") or "")).resolve() != here:
                continue
        except Exception:
            continue
        title = str(row.get("title") or "").strip()
        key = str(row.get("key") or "").strip()
        host = str(row.get("host") or "").strip()
        label = title or key or "unnamed"
        return f"{label} on {host}" if host else label
    return ""


def release_session_worktree(workspace) -> dict:
    """Release the session worktree ``workspace`` is in, if it is spare.

    Returns ``{"released": bool, "kept": str}`` -- ``kept`` says why, and
    is empty when there was nothing to release.

    The decision is ``worktree.exit_worktree(keep_if_changed=True)``, which
    already keeps a tree that has changes and one that background jobs are
    still running inside. That is the whole rule and it is not restated
    here: a second copy of "is this tree spare" would be a second answer
    that can disagree with the one the subagent teardown path uses.

    Never raises: closing a session must not fail because a directory
    could not be tidied.
    """
    from delfin.agent import worktree as _wt
    side = read_worktree_sidecar(workspace)
    if not side:
        return {"released": False, "kept": ""}
    try:
        info = _wt.WorktreeInfo(
            repo_dir=Path(side["repo_dir"]),
            path=Path(side["path"]),
            branch=str(side.get("branch") or ""),
            base_ref=str(side.get("base_ref") or ""),
            created_at=float(side.get("created_at") or 0.0),
        )
    except (KeyError, TypeError, ValueError):
        return {"released": False, "kept": "the worktree record is unreadable"}
    try:
        done = _wt.exit_worktree(info, keep_if_changed=True)
    except Exception as exc:
        return {"released": False, "kept": f"git refused the teardown: {exc}"}
    if done.held_by_jobs:
        return {"released": False,
                "kept": "background jobs are still running in it"}
    if done.had_changes:
        return {"released": False, "kept": "it has uncommitted changes"}
    return {"released": bool(done.cleaned_up), "kept": ""}


def _short_path(path: str) -> str:
    home = str(Path.home())
    return "~" + path[len(home):] if path.startswith(home) else path


def _display_path(path: str, limit: int = 34) -> str:
    """A working directory short enough for one line of the list: the home
    directory as ~, and a long path as its last two parts."""
    text = _short_path(path)
    if len(text) <= limit:
        return text
    parts = Path(text).parts
    return "…/" + "/".join(parts[-2:]) if len(parts) > 2 else text


# A message another session left arrives through the input box like a typed
# one. It is not the user's words: it must not name the session (the title
# then wrapped itself into every later message: "[Message from the session
# "[Message from ..." -- driven 2026-09-16), and it must not pin its language.
_DELIVERED_PREFIX = "[Message from the session"


def is_delivered_message(text: str) -> bool:
    return str(text or "").lstrip().startswith(_DELIVERED_PREFIX)


def _session_title(state: dict) -> str:
    for message in state.get("chat_messages") or []:
        if isinstance(message, dict) and message.get("role") == "user":
            text = " ".join(str(message.get("content") or "").split())
            if text and not is_delivered_message(text):
                return text[:38] + ("…" if len(text) > 38 else "")
    return "New session"


# ---------------------------------------------------------------------------
# The tab
# ---------------------------------------------------------------------------

def create_tab(ctx: Any, *, build: Optional[Callable] = None):
    """Build the session list and reopen the sessions that were open.

    ``build`` makes one session's tab from its context
    (``tab_agent.create_tab`` by default). Returns ``(widget, refs)``.
    """
    if build is None:
        from delfin.dashboard import tab_agent
        build = tab_agent.create_tab

    sessions: list[dict] = []
    view: dict[str, Any] = {"active": "", "restoring": True,
                            "scripts_added": False, "resume_at": 0.0}
    shared_status = getattr(ctx, "agent_status_html", None)

    def _classed(w, *names):
        for name in names:
            w.add_class(name)
        return w

    heading = _classed(widgets.HTML("Sessions"), "delfin-session-heading")
    new_btn = _classed(widgets.Button(description="+", tooltip="New session"),
                       "delfin-session-add")
    head = _classed(widgets.HBox([heading, new_btn]), "delfin-session-head")
    workdir_box = widgets.Combobox(
        options=workspace_choices(ctx), value=default_workspace(ctx),
        placeholder="Working directory", ensure_option=False,
        layout=widgets.Layout(width="100%"))
    start_btn = _classed(widgets.Button(description="Start"),
                         "delfin-session-start")
    cancel_btn = _classed(widgets.Button(description="Cancel"),
                          "delfin-session-cancel")
    own_worktree_box = widgets.Checkbox(
        value=False, description="Own worktree", indent=False,
        tooltip="Work on a branch of its own in a separate checkout",
        layout=widgets.Layout(display="none"))
    worktree_hint = _classed(widgets.HTML(""), "delfin-session-hint")
    new_form = _classed(widgets.VBox(
        [_classed(widgets.HTML("Working directory"), "delfin-session-label"),
         workdir_box, own_worktree_box, worktree_hint,
         _classed(widgets.HBox([cancel_btn, start_btn]),
                  "delfin-session-actions")],
        layout=widgets.Layout(display="none")), "delfin-session-form")
    notice = _classed(widgets.HTML(""), "delfin-session-notice")
    list_box = _classed(widgets.VBox(), "delfin-session-list")
    resume_dropdown = _classed(widgets.Dropdown(
        options=[("Resume a saved session…", "")], value=""),
        "delfin-session-resume")
    elsewhere_note = _classed(widgets.HTML(""), "delfin-session-elsewhere")
    stop_all_btn = _classed(widgets.Button(
        description="Stop all agents",
        tooltip=("Emergency stop: ends every agent of yours on every login "
                 "node -- dashboard kernels with their shells and MCP "
                 "servers, and the daemons. Schedules are disabled and "
                 "nothing starts on its own afterwards; cluster jobs keep "
                 "running.")), "delfin-session-stopall")
    sidebar = _classed(widgets.VBox(
        [widgets.HTML(_SIDEBAR_CSS), head, new_form, notice, list_box,
         resume_dropdown, elsewhere_note, stop_all_btn]), "delfin-sessions")
    stage = _classed(widgets.VBox(), "delfin-session-stage")
    splitter_host = _classed(widgets.HBox([widgets.HTML("")]),
                             "delfin-splitter-host")
    splitter = _classed(widgets.HTML(""), "delfin-session-splitter")
    splitter_host.children = (splitter,)
    widget = _classed(widgets.HBox([sidebar, splitter_host, stage]),
                      "delfin-session-shell")
    # Bind the drag once per shell when the whole page is assembled. The
    # init scripts are gathered by the dashboard into one keep_js() call, so
    # this runs after the shell is really on the page (the sidebar HTML
    # alone would render too early to find it).
    _register_init_js(ctx, _SPLITTER_INIT_JS)

    # -- lookup -------------------------------------------------------------

    def _find(key: str) -> Optional[dict]:
        return next((rec for rec in sessions if rec["key"] == key), None)

    def _state(rec: dict) -> dict:
        return rec["refs"].get("state") or {}

    def _mark(rec: dict) -> None:
        """One dot per row: what the session needs from you beats that it
        is working, which beats that it finished while you looked away."""
        state = _state(rec)
        busy = bool(state.get("streaming"))
        if rec.get("_was_busy") and not busy and view["active"] != rec["key"]:
            rec["unseen"] = True
        rec["_was_busy"] = busy
        wanted = session_dot(busy=busy, needs=needs_you(state),
                             unseen=bool(rec.get("unseen")))
        for cls in _DOT_CLASSES.values():
            if cls == wanted:
                rec["row"].add_class(cls)
            else:
                rec["row"].remove_class(cls)

    def _session_id(rec: dict) -> str:
        return str(_state(rec).get("active_session_id") or "")

    def _held_elsewhere(session_id: str, key: str) -> bool:
        return any(rec["key"] != key and _session_id(rec) == session_id
                   for rec in sessions)

    def _say(text: str) -> None:
        notice.value = html.escape(text) if text else ""

    # -- what is on screen ----------------------------------------------------

    def _render_list() -> None:
        rows = []
        for rec in sessions:
            if rec["key"] == view["active"]:
                rec["row"].add_class("delfin-session-active")
            else:
                rec["row"].remove_class("delfin-session-active")
            rows.append(rec["row"])
        list_box.children = tuple(rows)

    def _mirror_status(key: str, value: str) -> None:
        if key == view["active"] and shared_status is not None:
            try:
                shared_status.value = value
            except Exception:
                pass

    def activate(key: str) -> None:
        rec = _find(key)
        if rec is None:
            return
        view["active"] = key
        rec["unseen"] = False
        _mark(rec)
        for other in sessions:
            other["tab"].layout.display = "" if other["key"] == key else "none"
        state = _state(rec)
        try:
            ctx.agent_state = state
            ctx.agent_engine = state.get("engine")
        except Exception:
            pass
        _mirror_status(key, rec["ctx"].agent_status_html.value)
        _render_list()

    def _persist() -> None:
        if view["restoring"]:
            return
        save_open_sessions([
            {"session_id": _session_id(rec), "workspace": rec["workspace"]}
            for rec in sessions if _session_id(rec)])

    def _refresh_resume_options() -> None:
        try:
            from delfin.agent.session_store import list_sessions
            saved = list_sessions(limit=40)
        except Exception:
            saved = []
        open_ids = {_session_id(rec) for rec in sessions}
        options = [("Resume a saved session…", "")]
        for row in saved:
            sid = str(row.get("session_id") or "")
            if not sid or sid in open_ids:
                continue
            title = session_label(row.get("title"), sid)
            where = _short_path(str(row.get("workspace") or ""))
            options.append((f"{title} — {where}" if where else title, sid))
        resume_dropdown.options = options
        resume_dropdown.value = ""

    # -- open, close ----------------------------------------------------------

    def open_session(workspace: str = "", session_id: str = "") -> Optional[dict]:
        """Open a session working in ``workspace``; with ``session_id``, the
        saved conversation. A conversation already open is shown instead."""
        sid = str(session_id or "").strip()
        if sid:
            for rec in sessions:
                if _session_id(rec) == sid:
                    activate(rec["key"])
                    return rec
        if len(sessions) >= _MAX_OPEN:
            _say(f"{_MAX_OPEN} sessions are open — close one first.")
            return None
        _say("")
        key = uuid.uuid4().hex[:8]
        where = str(workspace or "").strip() or default_workspace(ctx)
        session_ctx = _SessionContext(
            ctx,
            agent_workspace=where,
            initial_session_id=sid,
            sessions_managed=True,
            session_open_elsewhere=(
                lambda other, _key=key: _held_elsewhere(other, _key)),
            agent_engine=None,
            agent_state=None,
            agent_status_html=widgets.HTML(""),
            add_init_js=(ctx.add_init_js if not view["scripts_added"]
                         else (lambda script: None)),
            presence_key=key,
        )
        tab, refs = build(session_ctx)
        view["scripts_added"] = True
        title_btn = _classed(widgets.Button(description="New session",
                                            tooltip=where),
                             "delfin-session-title")
        close_btn = _classed(widgets.Button(description="×", tooltip=(
            "Close this session — it stays saved and can be resumed")),
            "delfin-session-close")
        info = _classed(widgets.HTML(html.escape(_display_path(where))),
                        "delfin-session-path")
        row = _classed(widgets.HBox([
            _classed(widgets.VBox([title_btn, info]), "delfin-session-text"),
            close_btn]), "delfin-session-item")
        rec = {"key": key, "workspace": where, "tab": tab, "refs": refs or {},
               "ctx": session_ctx, "title_btn": title_btn, "info": info,
               "row": row}
        title_btn.on_click(lambda _b, _key=key: activate(_key))
        close_btn.on_click(lambda _b, _key=key: close_session(_key))
        session_ctx.agent_status_html.observe(
            lambda change, _key=key: _mirror_status(_key, change["new"]),
            names="value")
        sessions.append(rec)
        stage.children = tuple(r["tab"] for r in sessions)
        activate(key)
        refresh()
        _persist()
        return rec

    def close_session(key: str) -> None:
        """Close a session: it is saved, stops, and leaves the list."""
        rec = _find(key)
        if rec is None:
            return
        try:
            shutdown = rec["refs"].get("shutdown")
            if callable(shutdown):
                shutdown()
        except Exception:
            pass
        try:
            from delfin.agent import session_presence as _presence
            _presence.withdraw(key)
        except Exception:
            pass
        # Release the session's own worktree if it is spare. After the
        # shutdown above, so a background job that was holding the tree has
        # been asked to stop first; and after withdraw, so the presence
        # record is already gone when the decision is made. A tree with
        # changes or with jobs still in it is kept and said so -- the
        # decision is exit_worktree's, not this function's.
        try:
            _rel = release_session_worktree(rec.get("workspace") or "")
            if _rel.get("kept"):
                _say(f"The worktree was kept: {_rel['kept']}")
        except Exception:
            pass
        sessions.remove(rec)
        stage.children = tuple(r["tab"] for r in sessions)
        if not sessions:
            open_session()
            return
        if view["active"] == key:
            activate(sessions[-1]["key"])
        _render_list()
        _persist()
        _refresh_resume_options()

    # -- keeping the list current --------------------------------------------

    def refresh() -> None:
        """Titles, the working marker, the shared engine for the Activity
        tab, and the remembered list once a session has an id."""
        changed = False
        for rec in sessions:
            state = _state(rec)
            label = _session_title(state)
            if rec["title_btn"].description != label:
                rec["title_btn"].description = label
            _mark(rec)
            sid = _session_id(rec)
            if sid != rec.get("remembered_id"):
                rec["remembered_id"] = sid
                changed = True
            try:
                from delfin.agent import session_presence as _presence
                _presence.announce(rec["key"], session_id=sid,
                                   title=_session_title(state),
                                   workspace=rec["workspace"])
            except Exception:
                pass
            # Messages from other sessions, for every open session -- also
            # one that has not built an engine yet.
            try:
                deliver = rec["refs"].get("deliver_messages")
                if callable(deliver):
                    deliver()
            except Exception:
                pass
        active = _find(view["active"])
        if active is not None:
            try:
                ctx.agent_engine = _state(active).get("engine")
            except Exception:
                pass
        if changed:
            _persist()

    def _refresh_elsewhere() -> None:
        """Name the sessions open on other machines: the stop below reaches
        them, and nothing else on this page does."""
        try:
            import socket as _socket

            from delfin.agent import session_presence as _presence
            here = _socket.gethostname()
            counts: dict[str, int] = {}
            for record in _presence.open_sessions():
                host = str(record.get("host") or "")
                if host and host != here:
                    short = host.split(".")[0]
                    counts[short] = counts.get(short, 0) + 1
        except Exception:
            counts = {}
        elsewhere_note.value = (
            "Also open on other login nodes: " + html.escape(", ".join(
                f"{n} on {h}" for h, n in sorted(counts.items())))
            if counts else "")

    def _tick() -> None:
        try:
            refresh()
            import time as _time
            if _time.monotonic() >= view["resume_at"]:
                view["resume_at"] = _time.monotonic() + _RESUME_REFRESH_S
                _refresh_resume_options()
                _refresh_elsewhere()
        except Exception:
            pass
        finally:
            timer = threading.Timer(_REFRESH_S, _tick)
            timer.daemon = True
            timer.start()
            view["timer"] = timer

    # -- controls ---------------------------------------------------------------

    def _suggest_worktree(_change=None) -> None:
        """Offer a worktree inside a git repository, and pre-select it when
        another open session already works in that repository."""
        where = str(workdir_box.value or "").strip()
        try:
            from delfin.agent import session_presence as _presence
            in_repo = bool(where) and bool(_presence.repository_of(where)["root"])
            busy = _presence.in_same_repository(where) if in_repo else []
        except Exception:
            in_repo, busy = False, []
        own_worktree_box.layout.display = "" if in_repo else "none"
        own_worktree_box.value = bool(busy)
        if busy:
            names = ", ".join(
                session_label(r.get("title"), r.get("session_id"))
                for r in busy[:3])
            worktree_hint.value = (
                "<div style='font-size:10px;color:#546e7a'>Also working in "
                f"this repository: {html.escape(names)}. A worktree of its "
                "own keeps the two sessions apart.</div>")
        else:
            worktree_hint.value = ""

    def _on_new(_btn=None) -> None:
        workdir_box.options = workspace_choices(ctx)
        new_form.layout.display = ""
        new_btn.layout.display = "none"
        _suggest_worktree()

    def _on_cancel(_btn=None) -> None:
        new_form.layout.display = "none"
        new_btn.layout.display = ""

    def _on_start(_btn=None) -> None:
        where = str(workdir_box.value or "").strip()
        if where and not Path(where).expanduser().is_dir():
            _say(f"Not a directory: {where}")
            return
        target = str(Path(where).expanduser()) if where else ""
        # A directory another LIVE session is working in is refused, not
        # silently shared. The picker is a plain filesystem listing, so a
        # worktree belonging to a running session was offered and typeable;
        # starting there makes repository_of() resolve the worktree as the
        # root and nests a worktree-of-a-worktree under it. The hint
        # already existed -- _suggest_worktree says "also working in this
        # repository" and pre-ticks the box -- but nothing refused.
        #
        # Only when no fresh worktree is being made: ticking the box means
        # a tree of one's own, which is the answer to this, not the problem.
        _fresh = (own_worktree_box.value
                  and own_worktree_box.layout.display != "none")
        if target and not _fresh:
            _holder = _live_session_in(target)
            if _holder:
                _say(f"{_short_path(target)} is already being worked in by "
                     f"an open session ({_holder}). Tick \u201cOwn "
                     f"worktree\u201d to get a checkout of your own, or "
                     f"pick another directory.")
                return
        if own_worktree_box.value and own_worktree_box.layout.display != "none":
            try:
                target = session_worktree(target or default_workspace(ctx))
            except Exception as exc:
                _say(f"No worktree was created: {exc}")
                return
        _on_cancel()
        open_session(target)

    def _on_resume(change) -> None:
        sid = str(change.get("new") or "")
        if not sid:
            return
        data = _saved_session(sid) or {}
        open_session(str(data.get("workspace") or ""), sid)
        _refresh_resume_options()

    workdir_box.observe(_suggest_worktree, names="value")

    def _refresh_folder_choices(change=None) -> None:
        """The list under the box follows what is typed, like an editor's
        open dialog: the folders of the typed path, or of its parent."""
        try:
            typed = str((change or {}).get("new") if change else workdir_box.value or "")
            choices = workspace_choices(ctx, typed)
            if list(workdir_box.options) != choices:
                workdir_box.options = choices
        except Exception:
            pass

    workdir_box.observe(_refresh_folder_choices, names="value")
    def _on_stop_all(_btn=None) -> None:
        """Two clicks: the first arms the button for a few seconds."""
        import time as _time
        if view.get("stop_armed_until", 0.0) < _time.monotonic():
            view["stop_armed_until"] = _time.monotonic() + _STOP_ARM_S
            stop_all_btn.description = "Click again to stop everything"
            stop_all_btn.add_class("delfin-armed")

            def _disarm() -> None:
                if view.get("stop_armed_until", 0.0) <= _time.monotonic():
                    stop_all_btn.description = "Stop all agents"
                    stop_all_btn.remove_class("delfin-armed")

            timer = threading.Timer(_STOP_ARM_S + 0.5, _disarm)
            timer.daemon = True
            timer.start()
            return
        view["stop_armed_until"] = 0.0
        stop_all_btn.disabled = True
        stop_all_btn.description = "Stopping every agent…"
        try:
            _give_emergency_stop()
        except Exception as exc:
            stop_all_btn.disabled = False
            stop_all_btn.description = "Stop all agents"
            stop_all_btn.remove_class("delfin-armed")
            _say(f"The stop could not be given: {exc}. In a terminal: "
                 "delfin-agent stop-all")
            return
        _say("Emergency stop given: every agent on every login node ends "
             "within seconds, this dashboard's too. Reload the page to "
             "start again; nothing starts on its own until you send a "
             "message yourself.")

    stop_all_btn.on_click(_on_stop_all)
    new_btn.on_click(_on_new)
    cancel_btn.on_click(_on_cancel)
    start_btn.on_click(_on_start)
    resume_dropdown.observe(_on_resume, names="value")

    # -- start ------------------------------------------------------------------

    rows = load_open_sessions()
    wanted = (str(getattr(ctx, "initial_session_id", "") or "").strip()
              or os.environ.get("DELFIN_RESUME_SESSION", "").strip())
    if wanted == "latest":
        try:
            from delfin.agent.session_store import resume_latest
            wanted = str((resume_latest() or {}).get("session_id") or "")
        except Exception:
            wanted = ""
    if wanted and all(row["session_id"] != wanted for row in rows):
        rows.insert(0, {"session_id": wanted, "workspace": ""})
    for row in rows[:_MAX_OPEN]:
        data = _saved_session(row["session_id"])
        if not data:
            continue                    # deleted since: nothing to reopen
        try:
            open_session(row["workspace"] or str(data.get("workspace") or ""),
                         row["session_id"])
        except Exception:
            continue
    if not sessions:
        open_session()
    if wanted:
        for rec in sessions:
            if _session_id(rec) == wanted:
                activate(rec["key"])
    view["restoring"] = False
    _persist()
    # Reclaim session worktrees nobody is in. AFTER the open sessions are
    # restored, so their presence records exist and a tree one of them is
    # working in reads as live -- the saved-workspace check would catch it
    # anyway, and two independent reasons to keep a tree is the right
    # number for a pass that removes directories.
    #
    # Said out loud rather than done quietly: this removes a checkout, and
    # a reader who sees a folder disappear should be able to find out from
    # the session list why.
    try:
        _swept = reclaim_orphaned_worktrees(
            (sessions[0]["workspace"] if sessions else "")
            or default_workspace(ctx))
        if _swept["released"]:
            _say(f"Reclaimed {len(_swept['released'])} worktree(s) left by "
                 f"sessions that are gone.")
    except Exception:
        pass
    _tick()

    return widget, {
        "open": open_session,
        "close": close_session,
        "activate": activate,
        "sessions": lambda: list(sessions),
        "active": lambda: _find(view["active"]),
        "refresh": refresh,
        "form": {"new": new_btn, "workdir": workdir_box,
                 "own_worktree": own_worktree_box, "start": start_btn},
        # Exported so the reclamation can be driven in a test; a pass that
        # removes directories should not be reachable only at startup.
        "reclaim_worktrees": reclaim_orphaned_worktrees,
        "notice": notice,
    }
