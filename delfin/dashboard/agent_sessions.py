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
# How often, at most, a dashboard whose agent service went away tries to
# reach it again (on the refresh tick; no loop of its own).
_RECONNECT_S = 10.0

#: Ids of contexts whose splitter init script is already registered. The
#: dashboard builds the agent tab once, but a settings change can rebuild a
#: tab on the same context; without the guard the drag script would be
#: emitted twice into the one keep_js() payload.
_REGISTERED_INIT_JS = set()

#: The stop of every running refresh loop, keyed by ``id(ctx)`` of the
#: context whose session list started it. A list rebuilt on the same
#: context stops the previous loop; ``stop_refresh_loops`` stops them all.
_REFRESH_LOOPS: dict[int, Callable[[], None]] = {}
_REFRESH_LOOPS_LOCK = threading.Lock()


def stop_refresh_loops() -> int:
    """Stop every session list's refresh loop in this process.

    Output: how many loops were stopped. Each stop cancels the pending
    timer and waits for a tick already running to finish, so nothing a
    loop writes lands after this returns. Used at the end of a test, where
    a loop started by one test would otherwise keep writing presence
    records and the open-sessions file into the next test's redirected
    home; a dashboard kernel ends its loops by ending."""
    with _REFRESH_LOOPS_LOCK:
        stops = list(_REFRESH_LOOPS.values())
        _REFRESH_LOOPS.clear()
    for stop in stops:
        try:
            stop()
        except Exception:
            pass
    return len(stops)


def _register_init_js(ctx: Any, script: str) -> None:
    """Add ``script`` to a context's init-js payload, at most once per
    script. (Keyed on the context alone, the first script registered was
    the only one: the rename, folder and suggestion keys never arrived.)"""
    add = getattr(ctx, "add_init_js", None)
    if not callable(add):
        return
    key = (id(ctx), hash(script))
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
   browser's .calc-splitter. The shell aligns items at flex-start (the
   sidebar keeps its own height); the divider alone stretches over the
   whole column, from the control row to the Details bar. */
.delfin-session-shell { align-items: flex-start; }
/* The gap belongs to the SPLITTER, equally on both sides, and not to the
   sidebar. It was the sidebar's `margin-right: 12px`, which put the whole
   gap on one side: 12px between column and bar, nothing between bar and
   chat, so the handle sat flush against the chat with a wide empty strip
   behind it. One margin here is symmetric by construction -- there is no
   second number that can drift away from the first. */
.delfin-splitter-host { align-self: stretch; display: flex;
    flex-direction: column;
    margin: 0 var(--delfin-splitter-gap, 10px); }
/* The length is the column's, in CSS. It was measured by a script from the
   shell's top to the Send row's bottom and re-measured only when the Send
   row or the window changed size: its top sat 10 px below the column's,
   and anything that moved the Send row without resizing it left the
   length stale. The column fits the window since the agent tab's height
   script, so the full column is no longer the full transcript. */
.delfin-session-splitter { width: 8px;
    flex: 1 1 auto; min-height: 160px;
    cursor: col-resize; touch-action: none;
    background: linear-gradient(to right, #d6d6d6, #f2f2f2, #d6d6d6);
    border-radius: 4px; z-index: 10; pointer-events: auto; }
.delfin-session-splitter:hover { background: linear-gradient(
    to right, #b0b0b0, #e0e0e0, #b0b0b0); }
/* The working-directory list (the box's suggestions). */
.delfin-path-suggest { gap: 0; max-height: 220px; overflow-y: auto !important;
    overflow-x: hidden !important; border: 1px solid #e5e7eb;
    border-radius: 6px; background: #ffffff; }
.delfin-path-suggest:empty { display: none; }
.delfin-path-option { width: 100% !important; height: 24px !important;
    margin: 0 !important; padding: 0 8px !important; border: 0;
    background: transparent !important; text-align: left; font-size: 12px;
    color: #1f2937; white-space: nowrap; overflow: hidden;
    text-overflow: ellipsis; font-family: 'SF Mono', 'Consolas', monospace; }
.delfin-path-option:hover, .delfin-path-option.delfin-path-active {
    background: #dbeafe !important; color: #0f172a; }
.delfin-session-rename-input { width: 100% !important; }
/* Suggestions above the message box (_COMPLETE_JS). */
.delfin-complete-pop { position: fixed; z-index: 10000; max-height: 260px;
    overflow-y: auto; background: #ffffff; border: 1px solid #e5e7eb;
    border-radius: 10px; box-shadow: 0 8px 24px rgba(15, 23, 42, 0.14);
    padding: 4px; font-family: 'SF Mono', 'Consolas', monospace; font-size: 12px; }
.delfin-complete-item { display: flex; justify-content: space-between; gap: 12px;
    padding: 4px 8px; border-radius: 6px; cursor: pointer; color: #1f2937; }
.delfin-complete-item:hover, .delfin-complete-active { background: #dbeafe; color: #0f172a; }
.delfin-complete-label { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.delfin-complete-kind { color: #94a3b8; font-size: 11px; flex: 0 0 auto; }
/* "Newest" at the foot of the visible chat, also in a short chat. The
   button is the transcript's last child; sticky keeps it at the foot only
   once the chat scrolls, so in a new session it sat under the last
   message. margin-top:auto pushes it to the end of the (min 100 % high)
   column, the negative bottom margin into the 84 px foot reserve; with a
   working line it stays above that line (the sticky offset of the agent
   tab's own rule). */
.delfin-agent-chat-host > .delfin-agent-chat > .delfin-chat-jump {
    margin-top: auto !important; margin-bottom: -78px !important; }
.delfin-agent-chat-frame:has(> .delfin-agent-working-slot .delfin-agent-working)
    .delfin-agent-chat > .delfin-chat-jump { margin-bottom: 0 !important; }
</style>"""

#: Drag the sessions column to a width. Runs once per shell (data-bound
#: guard) after the page is assembled, along every tab's init script in the
#: one keep_js() call from the dashboard. The widths it writes override the
#: 236px CSS default on the sidebar, so a drag sticks and a refresh keeps
#: the choice. Mirrors the calc browser's .calc-splitter handler.
#: Esc in a rename box cancels it, in the session list as in the running
#: list (live_tab binds the same key for its own view; bound once here so
#: the list does not depend on a service view having been opened).
#: Folders shown under the working-directory box at once (the list
#: scrolls; workspace_choices bounds what is read at _PICKER_MAX).
_SUGGEST_SHOWN = 40

#: Keys of the working-directory box: Tab, Up/Down, Enter, Esc.
_PATH_KEYS_JS = """\
(function () {
    if (window.__delfinPathKeys) return;
    window.__delfinPathKeys = true;
    function opts(form) {
        return Array.prototype.slice.call(
            form.querySelectorAll('.delfin-path-suggest .delfin-path-option'));
    }
    function mark(list, i) {
        list.forEach(function (b, j) {
            b.classList.toggle('delfin-path-active', j === i);
            if (j === i && b.scrollIntoView) b.scrollIntoView({block: 'nearest'});
        });
    }
    document.addEventListener('keydown', function (e) {
        var t = e.target;
        if (!t || !t.closest || t.tagName !== 'INPUT') return;
        var host = t.closest('.delfin-path-input');
        var form = host && t.closest('.delfin-session-form');
        if (!form) return;
        var list = opts(form);
        var at = list.findIndex(function (b) {
            return b.classList.contains('delfin-path-active'); });
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
            if (!list.length) return;
            e.preventDefault();
            at = e.key === 'ArrowDown' ? Math.min(list.length - 1, at + 1)
                                       : Math.max(0, at - 1);
            mark(list, at);
        } else if (e.key === 'Tab' && !e.shiftKey && list.length) {
            e.preventDefault();
            (list[at >= 0 ? at : 0]).click();
            setTimeout(function () { t.focus(); }, 0);
        } else if (e.key === 'Enter') {
            e.preventDefault();
            if (at >= 0) { list[at].click(); setTimeout(function () { t.focus(); }, 0); return; }
            var start = form.querySelector('.delfin-session-start');
            // Let the typed value reach the kernel before Start reads it.
            if (start) setTimeout(function () { start.click(); }, 120);
        } else if (e.key === 'Escape') {
            var cancel = form.querySelector('.delfin-session-cancel');
            if (cancel) { e.preventDefault(); cancel.click(); }
        }
    }, true);
    document.addEventListener('click', function (e) {
        var b = e.target && e.target.closest && e.target.closest('.delfin-path-option');
        var form = b && b.closest('.delfin-session-form');
        var input = form && form.querySelector('.delfin-path-input input');
        if (input) setTimeout(function () { input.focus(); }, 0);
    }, true);
})();
"""

#: Suggestions in the message box (delfin.dashboard.completion): for a
#: `!command` and for a word that starts like a path. The box sends what is
#: typed up to the caret to a hidden field; the kernel answers in a hidden
#: HTML; a popup above the box shows the answer. Tab takes the highlighted
#: (or first) suggestion, Up/Down move, Enter takes one only when one was
#: highlighted (otherwise it sends, as always), Esc closes. Keys are taken
#: on `window` in the capture phase, ahead of the box's own Enter-to-send.
_COMPLETE_JS = """\
(function () {
    if (window.__delfinComplete) return;
    window.__delfinComplete = true;
    var st = {seq: 0, items: [], start: 0, at: -1, ta: null, timer: null, pop: null};
    function pop() {
        if (!st.pop) {
            st.pop = document.createElement('div');
            st.pop.className = 'delfin-complete-pop';
            st.pop.style.display = 'none';
            document.body.appendChild(st.pop);
            st.pop.addEventListener('mousedown', function (e) {
                var row = e.target.closest('.delfin-complete-item');
                if (!row) return;
                e.preventDefault();
                take(parseInt(row.getAttribute('data-i'), 10));
            });
        }
        return st.pop;
    }
    function hide() { st.items = []; st.at = -1; if (st.pop) st.pop.style.display = 'none'; }
    function wanted(before) {
        var line = before.slice(before.lastIndexOf('\\n') + 1);
        if (before.indexOf('\\n') < 0 && line.charAt(0) === '!' && line.length > 1) return true;
        var word = line.split(/\\s/).pop();
        return /^(\\/|~\\/|~$|\\.\\/|\\.\\.\\/)/.test(word);
    }
    function ask(ta) {
        var before = ta.value.slice(0, ta.selectionStart);
        if (!wanted(before)) { hide(); return; }
        var q = document.querySelector('.delfin-complete-query input');
        if (!q) return;
        st.ta = ta;
        st.seq += 1;
        q.value = JSON.stringify({text: before, seq: st.seq});
        q.dispatchEvent(new Event('input', {bubbles: true}));
    }
    function render() {
        var p = pop();
        if (!st.items.length || !st.ta) { hide(); return; }
        var r = st.ta.getBoundingClientRect();
        p.innerHTML = '';
        st.items.forEach(function (it, i) {
            var row = document.createElement('div');
            row.className = 'delfin-complete-item' + (i === st.at ? ' delfin-complete-active' : '');
            row.setAttribute('data-i', String(i));
            var lab = document.createElement('span');
            lab.className = 'delfin-complete-label';
            lab.textContent = it.label;
            var kind = document.createElement('span');
            kind.className = 'delfin-complete-kind';
            kind.textContent = it.kind;
            row.appendChild(lab); row.appendChild(kind);
            p.appendChild(row);
        });
        p.style.left = Math.max(8, r.left) + 'px';
        p.style.width = Math.min(Math.max(320, r.width * 0.6), 640) + 'px';
        p.style.bottom = (window.innerHeight - r.top + 6) + 'px';
        p.style.display = 'block';
        var act = p.querySelector('.delfin-complete-active');
        if (act && act.scrollIntoView) act.scrollIntoView({block: 'nearest'});
    }
    function take(i) {
        var it = st.items[i], ta = st.ta;
        if (!it || !ta) return;
        var caret = ta.selectionStart;
        var v = ta.value.slice(0, st.start) + it.insert + ta.value.slice(caret);
        ta.value = v;
        var pos = st.start + it.insert.length;
        ta.setSelectionRange(pos, pos);
        ta.dispatchEvent(new Event('input', {bubbles: true}));
        ta.focus();
        hide();
        if (/\\/$/.test(it.insert)) ask(ta);   // a folder: go on into it
    }
    document.addEventListener('input', function (e) {
        var t = e.target;
        if (!t || t.tagName !== 'TEXTAREA' || !t.closest || !t.closest('.delfin-agent-input')) return;
        clearTimeout(st.timer);
        st.timer = setTimeout(function () { ask(t); }, 110);
    }, true);
    window.addEventListener('keydown', function (e) {
        if (!st.items.length || !st.ta || e.target !== st.ta) return;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
            e.preventDefault(); e.stopImmediatePropagation();
            var n = st.items.length;
            st.at = e.key === 'ArrowDown' ? (st.at + 1) % n : (st.at <= 0 ? n - 1 : st.at - 1);
            render();
        } else if (e.key === 'Tab' && !e.shiftKey) {
            e.preventDefault(); e.stopImmediatePropagation();
            take(st.at >= 0 ? st.at : 0);
        } else if (e.key === 'Enter' && !e.shiftKey && st.at >= 0) {
            e.preventDefault(); e.stopImmediatePropagation();
            take(st.at);
        } else if (e.key === 'Escape') {
            e.preventDefault(); e.stopImmediatePropagation();
            hide();
        }
    }, true);
    document.addEventListener('focusout', function (e) {
        if (e.target === st.ta) setTimeout(hide, 150);
    }, true);
    new MutationObserver(function () {
        var el = document.querySelector('.delfin-complete-result [data-json]');
        if (!el) return;
        var data;
        try { data = JSON.parse(el.getAttribute('data-json')); } catch (_e) { return; }
        if (!data || data.seq !== st.seq) return;
        st.items = data.items || []; st.start = data.start || 0; st.at = -1;
        render();
    }).observe(document.documentElement, {childList: true, subtree: true, characterData: true});
})();
"""

_RENAME_KEYS_JS = """\
(function () {
    if (window.__delfinRenameKeys) return;
    window.__delfinRenameKeys = true;
    document.addEventListener('keydown', function (e) {
        var t = e.target;
        if (!t || !t.closest || e.key !== 'Escape') return;
        var box = t.closest('.delfin-live-rename');
        var cancel = box && box.querySelector('.delfin-live-rename-cancel');
        if (cancel) { e.preventDefault(); e.stopPropagation(); cancel.click(); }
    }, true);
})();
"""

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
        item = {"session_id": str(row["session_id"]).strip(),
                "workspace": str(row.get("workspace") or "")}
        if row.get("live"):
            # Ran in the agent service: reopened as a view on it.
            item["live"] = True
        out.append(item)
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

    A directory inside a session worktree that is held -- see
    ``_worktree_hold_reason`` -- is NOT offered: starting there nests a
    worktree-of-a-worktree. A worktree whose session is gone stays offered.
    When the typed path is a held tree, ``_suggest_worktree`` says why it
    is missing from the list.

    Cost: the presence directory is read at most once per call, and the
    decision is made once per sidecar -- the subfolders of one worktree all
    find the same record (measured before: 199 presence reads for 200
    entries, 200 ms per keystroke).
    """
    rows: list = []
    decided: dict = {}

    def _offer(text: str) -> bool:
        text = str(text or "").strip()
        if not text:
            return False
        try:
            side = read_worktree_sidecar(text)
        except Exception:
            return True
        if side is None:
            return True
        key = str(side.get("path") or text)
        if key not in decided:
            if not rows:
                rows.append(_open_session_rows())
            try:
                decided[key] = not _worktree_hold_reason(
                    key, side, sessions=rows[0])
            except Exception:
                decided[key] = True
        return decided[key]

    out: list[str] = []
    for d in (default_workspace(ctx), getattr(ctx, "agent_dir", None),
              getattr(ctx, "calc_dir", None), getattr(ctx, "repo_dir", None)):
        text = str(d or "").strip()
        if _offer(text) and text not in out:
            out.append(text)
    for text in _folders_like(typed):
        if _offer(text) and text not in out:
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


#: The sidecar a session worktree carries (delfin.agent.session_worktree,
#: which makes the tree; read below by read_worktree_sidecar).
from delfin.agent.session_worktree import WORKTREE_SIDECAR as _WORKTREE_SIDECAR  # noqa: E402


def _write_worktree_sidecar(info) -> None:
    """Record repo, branch and owner for a session worktree. Never raises.
    (delfin.agent.session_worktree.write_worktree_sidecar.)"""
    from delfin.agent import session_worktree as _sw
    _sw.write_worktree_sidecar(info)


def session_worktree(workspace: str) -> str:
    """A fresh worktree of the repository ``workspace`` is in, on a branch of
    its own, for a session that must not share a checkout with another.

    ``delfin.agent.session_worktree.create`` -- the service's
    ``session.create`` with ``worktree`` makes the same tree. Returns the
    path at the same place inside the worktree as ``workspace`` is inside
    the repository. Raises ValueError outside a repository.
    """
    from delfin.agent import session_worktree as _sw
    return _sw.create(workspace)[0]


#: The start form's label for the box inside a repository, where an own
#: worktree costs nothing else. In a plain folder the label says what the
#: start does there (a local project) -- OWN_WORKTREE_PLAIN_LABEL below.
OWN_WORKTREE_LABEL_IN_REPO = "Own worktree"
OWN_WORKTREE_PLAIN_LABEL = ("Own worktree (makes this folder a local project)")


def plain_workspace(workspace: str) -> bool:
    """Whether ``workspace`` is a directory outside any git repository
    (the case the own-worktree offer needs consent in)."""
    try:
        from delfin.agent import session_presence as _presence
        return bool(workspace) and not bool(
            _presence.repository_of(workspace)["root"])
    except Exception:
        return False


def start_session_worktree(workspace: str, *, ask=None) -> dict:
    """A session worktree in a PLAIN folder: consent first, then the
    engine's ``session_worktree.start``.

    ``ask`` receives the consent question and returns True to proceed;
    None (nobody to ask) is a no. The folder is not touched before a yes;
    a refusal (home, system paths, inside another repository) raises
    ValueError with the reason. Inside a repository no consent is needed:
    the worktree costs nothing there, and ``session_worktree`` above
    applies.
    """
    from delfin.agent import session_worktree as _sw

    ws = str(Path(str(workspace).strip() or ".").expanduser())
    if not plain_workspace(ws):
        # A repository already: the plain consent path does not apply.
        return {"path": session_worktree(ws)}
    from delfin.agent import local_project as _lp
    if not (callable(ask) and ask(_lp.consent_question(_lp.plan(ws)))):
        raise PermissionError("no consent: the folder was not changed")
    try:
        made = _sw.start(ws, init_project=True)
    except _sw.NeedsLocalProject as exc:          # pragma: no cover
        raise PermissionError(exc.question) from None
    except Exception as exc:
        raise ValueError(str(exc)) from None
    # The folder is a repository now: the presence layer's git cache still
    # says plain for it (a TTL cache), which would send the next start back
    # through the consent path.
    try:
        from delfin.agent import session_presence as _presence
        _presence._git_cache.pop(ws, None)
        _presence._git_cache.pop(str(Path(ws).resolve()), None)
    except Exception:
        pass
    return made


def merge_session_worktree(workspace: str) -> dict:
    """Merge the session worktree ``workspace`` is in back into its main
    checkout: the engine's ``session_worktree.merge``. Returns
    ``{"status": "ok"|"conflict"|"error", "files", "message", "target"}``;
    a conflict is reported, the main checkout is left as it was."""
    from delfin.agent import session_worktree as _sw
    return _sw.merge(workspace)


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
    import socket
    tree = str(side.get("path") or "")
    if not tree:
        return "the record names no path"
    here = (socket.gethostname() or "").strip()
    if str(side.get("host") or "").strip() != here:
        return "it belongs to another host"
    owner = _worktree_owner_alive(side)
    if owner is None:
        return f"its owning process ({side.get('pid')}) may still be running"
    if owner:
        return f"its owning process ({side.get('pid')}) is still running"
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


def _worktree_owner_alive(side: dict):
    """Whether the process that made the worktree in ``side`` still runs.

    True / False, or None when it cannot be asked. ``proc_identity.alive``
    with the recorded pid, start time and host: a pid written on another
    host is None (it names nothing here), and a pid whose start time
    differs is a later process, so False. A sidecar written before the
    start time was recorded degrades to the plain pid check.
    """
    from delfin.agent import proc_identity as _proc_identity
    try:
        return _proc_identity.alive(
            int(side.get("pid") or 0), str(side.get("proc_start") or ""),
            str(side.get("host") or "").strip())
    except (TypeError, ValueError):
        return False


def _open_session_rows() -> list:
    """session_presence's open sessions, or [] when they cannot be read."""
    try:
        from delfin.agent import session_presence as _presence
        return list(_presence.open_sessions() or ())
    except Exception:
        return []


def _worktree_hold_reason(workspace, side: dict | None = None, *,
                          sessions=None) -> str:
    """Why the session worktree ``workspace`` is in may not be started in,
    or "" when it may (or when ``workspace`` is not in a session worktree).

    Two holders: a live session working in the tree (session_presence),
    and the process that made it, still running on this host. The second
    covers a live session whose presence record went stale. An owner on
    another host cannot be asked and does not hide the tree; a live session
    there is found by the first check, which trusts other hosts' records.

    ``sessions`` is a prefetched ``open_sessions()`` list so a caller that
    asks about many paths reads the presence directory once.
    """
    if side is None:
        side = read_worktree_sidecar(workspace)
    if side is None:
        return ""
    tree = str(side.get("path") or workspace)
    # inside=True: session_worktree hands a session the path at the same
    # place inside the tree as its workspace was inside the repository, so
    # a session working in tree/sub holds the tree.
    holder = _live_session_in(tree, rows=sessions, inside=True)
    if holder:
        return f"a live session is working in it ({holder})"
    if _worktree_owner_alive(side) is True:
        return ("the dashboard process that created it is still running "
                f"(pid {side.get('pid')})")
    return ""


def _live_session_in(workspace, rows=None, *, inside: bool = False) -> str:
    """A label for the live session whose workspace is ``workspace``, or "".

    ``inside=True`` also counts a session whose workspace lies below
    ``workspace``. ``rows`` is a prefetched ``open_sessions()`` list.

    Liveness is session_presence's answer, not a second one: its records
    carry pid and host and its reader already drops the stale ones. A
    record from another host is trusted as live -- a pid from another login
    node names nothing here, so it cannot be checked and must not be
    assumed dead.
    """
    try:
        here = Path(workspace).expanduser().resolve()
    except Exception:
        return ""
    if rows is None:
        rows = _open_session_rows()
    for row in rows or ():
        try:
            ws = Path(str(row.get("workspace") or "")).resolve()
            if ws != here and not (inside and here in ws.parents):
                continue
        except Exception:
            continue
        title = str(row.get("title") or "").strip()
        key = str(row.get("key") or "").strip()
        host = str(row.get("host") or "").strip()
        label = title or key or "unnamed"
        return f"{label} on {host}" if host else label
    return ""


def session_worktree_state(workspace, *, hold_fn=None) -> dict:
    """Classify a directory as a session worktree: its holding state.

    Returns ``{"is_worktree": bool, "holder": str, "why": str,
    "releasable": bool}``.

      - ``is_worktree`` — the path carries a ``.delfin/session_worktree.json``
        sidecar, so it is a session worktree on disk (or a session is started
        inside one).
      - ``holder`` — the label of a LIVE session working in it; empty when no
        live session is in it (its owning session has gone, or the path is
        not a worktree at all).
      - ``why`` — ``_worktree_hold_reason``: why the picker does not offer
        it, "" when it does. Set also when no session holds it but the
        process that created it still runs, so the start form can say why
        the entry is missing.
      - ``releasable`` — a worktree spare enough to hand back: it is a
        session worktree and ``_worktree_is_orphaned`` says it may be
        reclaimed (owning process dead on this host, no live session, no
        saved session that would reopen it). Releasing never deletes a
        checkout; that decision is ``exit_worktree``'s. A worktree whose
        owning process is still alive is never releasable, even when its
        session has left the presence file.

    ``hold_fn`` is the liveness probe (default module ``_live_session_in``).
    It is injectable so the picker -- and tests -- can drive the held/gone
    decision without touching the presence files. Never raises: a directory
    that cannot be classified is simply "not a worktree".
    """
    if hold_fn is None:
        hold_fn = _live_session_in
    side = read_worktree_sidecar(workspace)
    is_worktree = side is not None
    holder = ""
    why = ""
    releasable = False
    if is_worktree:
        try:
            holder = str(hold_fn(workspace) or "")
        except Exception:
            holder = ""
        try:
            why = (f"a live session is working in it ({holder})" if holder
                   else _worktree_hold_reason(workspace, side))
        except Exception:
            why = ""
        # Releasing must agree with the orphan check: an alive owning process
        # is not releasable even when the session has gone from the presence
        # file, and a worktree a SAVED session would reopen is not releasable
        # either (the same saved-workspaces the reclaim sweep passes, so the
        # spare answer never disagrees with the release path). ``_worktree_is_
        # orphaned`` owns all of this.
        try:
            releasable = (_worktree_is_orphaned(
                side, saved_workspaces=_saved_session_workspaces()) == "")
        except Exception:
            releasable = False
    return {
        "is_worktree": is_worktree,
        "holder": holder,
        "why": why,
        "releasable": releasable,
    }


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


def _session_title(state: dict, session_id: str = "") -> str:
    """The name in the list: the one the person gave it, else the start of
    its first message. A chosen name lives in ``_messenger_titles`` -- the
    map the agent tab's own save writes as the saved session's title, so
    it survives a reload and shows in Resume."""
    chosen = (state.get("_messenger_titles") or {}).get(session_id) \
        if session_id else ""
    if chosen:
        return str(chosen)
    for message in state.get("chat_messages") or []:
        if isinstance(message, dict) and message.get("role") == "user":
            text = " ".join(str(message.get("content") or "").split())
            if text and not is_delivered_message(text):
                return text[:38] + ("…" if len(text) > 38 else "")
    return "New session"


# ---------------------------------------------------------------------------
# Sessions in the agent service
# ---------------------------------------------------------------------------

#: Where the dashboard starts its sessions: "kernel" (the default: this
#: dashboard's kernel, shown in the agent tab's full chat) or "service" (the
#: person's agent service, so the terminal and the phone open them too; shown
#: in live_tab's view, which does not yet have the chat's controls). Either
#: way the sidebar lists the sessions a running service has and opens them.
#: The environment variable wins over agent.dashboard_sessions in the settings.
SESSIONS_ENV = "DELFIN_DASHBOARD_SESSIONS"

#: The permission profiles the agent tab's ladder stores, and the posture
#: each starts a service session in. Bypass cannot be a start posture
#: (the worker's --permission-mode stops at acceptEdits); it is switched to
#: after the start, confirmed, and refused where policy or isolation do not
#: allow it -- the page then says so.
_PROFILE_START = {"plan": "plan", "ask_all": "default", "repo_free": "acceptEdits",
                  "all_free": "acceptEdits"}
_PROFILE_OLD_NAMES = {"default": "ask_all", "erlaubt": "repo_free", "full": "all_free",
                      "acceptEdits": "repo_free", "auto": "all_free",
                      "bypassPermissions": "all_free"}


def sessions_where() -> str:
    """``"kernel"`` (default) or ``"service"``: where this dashboard starts sessions."""
    raw = os.environ.get(SESSIONS_ENV, "").strip().lower()
    if not raw:
        try:
            from delfin.user_settings import load_settings
            raw = str(((load_settings() or {}).get("agent") or {})
                      .get("dashboard_sessions") or "").strip().lower()
        except Exception:
            raw = ""
    return "service" if raw == "service" else "kernel"


def dashboard_session_settings() -> dict:
    """The create settings of a new service session: the agent tab's saved
    choices (provider, model, effort, agent mode, permission profile), so a
    session started here runs as the in-kernel tab would have run it.
    ``_bypass``: the profile is Bypass (switched to after the start)."""
    try:
        from delfin.user_settings import load_settings
        saved = (load_settings() or {}).get("agent", {}) or {}
    except Exception:
        saved = {}
    out: dict = {}
    provider = str(saved.get("provider") or "")
    if provider in ("claude", "openai", "kit"):
        out["provider"] = provider
    for key in ("model", "effort", "mode"):
        if str(saved.get(key) or "").strip():
            out[key] = str(saved[key]).strip()
    profile = str(saved.get("permission_profile", saved.get("permission_mode", "")) or "")
    profile = _PROFILE_OLD_NAMES.get(profile, profile)
    if profile in _PROFILE_START:
        out["permissionMode"] = _PROFILE_START[profile]
        out["_bypass"] = profile == "all_free"
    return out


def _service_connector(service: Any, explicit_build: bool):
    """What gives the session list its agent-service client, or None for
    sessions in the kernel only. ``service``: False (none), a callable
    returning a ``live_sessions.LiveClient``, a client, or None -- the
    person's service unless a stand-in ``build`` was given. With sessions in
    the kernel (the default) a service is only joined when one runs: its
    sessions (started in the terminal or on the phone) are listed, and
    nothing is started for the list alone."""
    if service is False:
        return None
    if callable(service):
        return service
    if service is not None:
        return lambda: service
    if explicit_build:
        return None
    from delfin.dashboard import live_sessions
    start = sessions_where() == "service"
    return lambda: live_sessions.connect(start=start)


_RUNNING_CSS = """<style>
.delfin-running-head { align-items: center; justify-content: space-between; padding: 6px 2px 0 6px; }
.delfin-running-count { font-size: 11px; color: #9a3412; font-weight: 600; }
.delfin-running-list { gap: 2px; max-height: calc(50vh - 120px); overflow-y: auto !important;
    overflow-x: hidden !important; }
.delfin-running-item { width: 100%; align-items: center; gap: 2px; padding: 4px 4px 4px 8px;
    border-radius: 8px; }
.delfin-running-item:hover { background: #eef2f7; }
.delfin-running-item.delfin-running-here { box-shadow: inset 3px 0 0 #93c5fd; }
.delfin-running-item.delfin-running-needs { background: #fff7ed; }
.delfin-running-item.delfin-running-needs .delfin-running-state { color: #c2410c; font-weight: 700; }
.delfin-running-state { font-weight: 600; }
.delfin-running-idle { color: #b45309; font-weight: 600; }
.delfin-running-item .delfin-running-stop { opacity: 1; color: #b91c1c; }
.delfin-running-rename { gap: 4px; }
.delfin-running-rename .widget-text { width: 100% !important; margin: 0; }
.delfin-running-rename .jupyter-button { width: auto !important; height: 24px !important;
    padding: 0 6px !important; font-size: 11px; }
.delfin-running-empty { font-size: 11px; color: #6b7280; padding: 2px 6px; }
</style>"""


class RunningPanel:
    """The sidebar's "Running sessions": every session of the person's
    agent service -- all folders, started anywhere -- with its name, folder,
    model, state and time since its last activity; "needs you" first and
    marked. One click opens it here; the pencil renames it (Enter saves,
    Esc cancels). Kept current by the service's summary events: an event
    marks it changed and arms one short timer, nothing polls."""

    RENDER_DELAY_S = 0.25

    def __init__(self, client: Any, on_open: Callable[[str], None],
                 is_open: Callable[[str], bool]) -> None:
        from delfin.dashboard import live_sessions as LS

        self.client = client
        self.on_open = on_open
        self.is_open = is_open
        self.running = LS.RunningList()
        self.buttons: dict[str, widgets.Button] = {}
        #: Each row's Stop, its confirmation and the two answers: nothing
        #: runs on forgotten (a session outlives every window).
        self.stops: dict[str, widgets.Button] = {}
        self.confirms: dict[str, widgets.HBox] = {}
        self.confirm_yes: dict[str, widgets.Button] = {}
        self.confirm_no: dict[str, widgets.Button] = {}
        self._rows: dict[str, dict] = {}
        self._timer: Optional[threading.Timer] = None
        #: Held for a whole render, and for the timer: the event timer and
        #: the 3 s tick both render, and two at once built two rows for one
        #: new session (re-entrant: render clears the timer under it).
        self._lock = threading.RLock()
        self._shown = -1
        self._unsubscribe: Optional[Callable[[], None]] = None
        heading = widgets.HTML("Running sessions")
        heading.add_class("delfin-session-heading")
        self.count = widgets.HTML("")
        self.count.add_class("delfin-running-count")
        head = widgets.HBox([heading, self.count])
        head.add_class("delfin-running-head")
        self.list_box = widgets.VBox()
        self.list_box.add_class("delfin-running-list")
        self.empty = widgets.HTML("No session runs in the agent service.")
        self.empty.add_class("delfin-running-empty")
        self.widget = widgets.VBox([widgets.HTML(_RUNNING_CSS), head, self.empty, self.list_box])

    # -- feeding it ---------------------------------------------------------------
    def attach(self, client: Any) -> None:
        """Follow ``client``'s service: its list now, its events from here."""
        if self._unsubscribe is not None:
            self._unsubscribe()
        self.client = client
        self._unsubscribe = client.subscribe(self._on_event)
        try:
            self.running.load(client.sessions())
        except Exception:
            pass
        self.render()

    def detach(self) -> None:
        if self._unsubscribe is not None:
            self._unsubscribe()
            self._unsubscribe = None
        with self._lock:
            if self._timer is not None:
                self._timer.cancel()
                self._timer = None

    def _on_event(self, event: dict) -> None:
        if self.running.apply(event):
            with self._lock:
                if self._timer is not None:
                    return
                self._timer = threading.Timer(self.RENDER_DELAY_S, self.render)
                self._timer.daemon = True
                timer = self._timer
            timer.start()

    def rows(self) -> list[dict]:
        return self.running.rows()

    # -- what is shown ------------------------------------------------------------
    def _row(self, sid: str) -> dict:
        row = self._rows.get(sid)
        if row is not None:
            return row
        title = widgets.Button(description="", tooltip="Open this session here")
        title.add_class("delfin-session-title")
        info = widgets.HTML("")
        info.add_class("delfin-session-path")
        pencil = widgets.Button(description="✎", tooltip="Rename (Enter saves, Esc cancels)")
        pencil.add_class("delfin-session-close")
        edit = widgets.Text(continuous_update=False,
                            placeholder="Name (empty: folder and short id)")
        save = widgets.Button(description="Save")
        cancel = widgets.Button(description="Cancel")
        cancel.add_class("delfin-live-rename-cancel")
        rename = widgets.HBox([edit, save, cancel], layout=widgets.Layout(display="none"))
        rename.add_class("delfin-running-rename")
        rename.add_class("delfin-live-rename")
        stop = widgets.Button(description="■", tooltip=(
            "Stop this session for every window (it is saved and can be resumed)"))
        stop.add_class("delfin-session-close")
        stop.add_class("delfin-running-stop")
        yes = widgets.Button(description="End it", button_style="danger")
        no = widgets.Button(description="Keep it")
        confirm = widgets.HBox([widgets.HTML("End for every window?"), yes, no],
                               layout=widgets.Layout(display="none"))
        confirm.add_class("delfin-running-rename")
        text = widgets.VBox([title, info, rename, confirm])
        text.add_class("delfin-session-text")
        box = widgets.HBox([text, pencil, stop])
        box.add_class("delfin-running-item")
        row = {"box": box, "title": title, "info": info, "pencil": pencil, "edit": edit,
               "rename": rename, "confirm": confirm, "classes": set()}

        def ask_stop(_b=None, _row=row):
            _row["confirm"].layout.display = ""

        def keep(_b=None, _row=row):
            _row["confirm"].layout.display = "none"

        def end(_b=None, _row=row, _sid=sid):
            _row["confirm"].layout.display = "none"
            try:
                self.client.stop(_sid)
            except Exception as exc:
                _row["info"].value = html.escape(f"not stopped: {getattr(exc, 'message', exc)}")

        stop.on_click(ask_stop)
        no.on_click(keep)
        yes.on_click(end)
        self.stops[sid], self.confirms[sid] = stop, confirm
        self.confirm_yes[sid], self.confirm_no[sid] = yes, no

        def open_it(_b=None, _sid=sid):
            self.on_open(_sid)

        def start_rename(_b=None, _row=row, _sid=sid):
            summary = self.running.summary(_sid) or {}
            _row["edit"].value = str(summary.get("title") or "")
            _row["rename"].layout.display = ""
            _row["title"].layout.display = "none"

        def end_rename(_b=None, _row=row):
            _row["rename"].layout.display = "none"
            _row["title"].layout.display = ""

        def save_rename(_b=None, _row=row, _sid=sid):
            try:
                self.client.rename(_sid, _row["edit"].value)
            except Exception as exc:
                _row["info"].value = html.escape(f"not renamed: {getattr(exc, 'message', exc)}")
                return
            end_rename()

        def on_enter(change, _row=row):
            if _row["rename"].layout.display != "none" and change.get("new") != change.get("old"):
                save_rename()

        title.on_click(open_it)
        pencil.on_click(start_rename)
        save.on_click(save_rename)
        cancel.on_click(end_rename)
        edit.observe(on_enter, names="value")
        self._rows[sid] = row
        self.buttons[sid] = title
        return row

    def render(self) -> None:
        """The rows as the list stands; widgets are reused, only text and
        order change."""
        with self._lock:
            self._timer = None
            self._render()

    def _render(self) -> None:
        rows = self.running.rows()
        children = []
        needs = 0
        for r in rows:
            sid = r["sessionId"]
            row = self._row(sid)
            if row["title"].description != r["name"]:
                row["title"].description = r["name"]
            row["title"].tooltip = f"{r['name']} ({r['shortId']}): open it here"
            state = (f'<span class="delfin-running-state">{html.escape(r["state"])}</span>')
            bits = [html.escape(r["folder"]), html.escape(r["model"]), state, html.escape(r["age"])]
            if r["needs"] and r["waiting"]:
                bits.append(html.escape(r["waiting"][:60]))
            if r["idle"]:
                bits.append(f'<span class="delfin-running-idle">{html.escape(r["idle"])}</span>')
            value = " · ".join(b for b in bits if b)
            if row["info"].value != value:
                row["info"].value = value
            wanted = {"delfin-running-needs"} if r["needs"] else set()
            if self.is_open(sid):
                wanted.add("delfin-running-here")
            for cls in ("delfin-running-needs", "delfin-running-here"):
                if cls in wanted and cls not in row["classes"]:
                    row["box"].add_class(cls)
                elif cls not in wanted and cls in row["classes"]:
                    row["box"].remove_class(cls)
            row["classes"] = wanted
            needs += r["needs"]
            children.append(row["box"])
        for sid in [s for s in self._rows if s not in {r["sessionId"] for r in rows}]:
            self._rows.pop(sid, None)
            for widgets_of in (self.buttons, self.stops, self.confirms, self.confirm_yes,
                               self.confirm_no):
                widgets_of.pop(sid, None)
        if tuple(self.list_box.children) != tuple(children):
            self.list_box.children = tuple(children)
        self.empty.layout.display = "none" if rows else ""
        count = f"{needs} need{'s' * (needs == 1)} you" if needs else ""
        if self.count.value != count:
            self.count.value = count


# ---------------------------------------------------------------------------
# The tab
# ---------------------------------------------------------------------------

def create_tab(ctx: Any, *, build: Optional[Callable] = None, service: Any = None,
               session_defaults: Optional[Callable[[], dict]] = None):
    """Build the session list and reopen the sessions that were open.

    ``build`` makes one in-kernel session's tab from its context
    (``tab_agent.create_tab`` by default). ``service``: where sessions run
    (see ``_service_connector``; by default the person's agent service,
    started if needed, and the kernel when there is none -- said on the
    page). ``session_defaults``: the create settings of a new service
    session (``dashboard_session_settings`` by default). Returns
    ``(widget, refs)``.
    """
    connector = _service_connector(service, explicit_build=build is not None)
    if build is None:
        from delfin.dashboard import tab_agent
        build = tab_agent.create_tab
    if session_defaults is None:
        session_defaults = dashboard_session_settings

    sessions: list[dict] = []
    view: dict[str, Any] = {"active": "", "restoring": True,
                            "scripts_added": False, "resume_at": 0.0,
                            # A stand-in service (tests, explicit callers) starts
                            # sessions there; the person's own only when set so.
                            "start_in_service": connector is not None and (
                                service is not None or sessions_where() == "service")}
    shared_status = getattr(ctx, "agent_status_html", None)

    def _classed(w, *names):
        for name in names:
            w.add_class(name)
        return w

    heading = _classed(widgets.HTML("Sessions"), "delfin-session-heading")
    new_btn = _classed(widgets.Button(description="+", tooltip="New session"),
                       "delfin-session-add")
    head = _classed(widgets.HBox([heading, new_btn]), "delfin-session-head")
    # The working directory is picked the way a terminal or an editor's
    # quick-open does it: type a path, the folders that fit are listed
    # under the box; Tab takes the highlighted (or first) one, Up/Down move
    # the highlight, a click goes into a folder, Enter starts (or takes the
    # highlighted folder first), Esc closes the card. The browser's own
    # datalist could do none of that.
    workdir_box = widgets.Text(
        value=default_workspace(ctx), continuous_update=True,
        placeholder="Working directory -- type a path, Tab completes",
        layout=widgets.Layout(width="100%"))
    workdir_box.add_class("delfin-path-input")
    path_suggest = _classed(widgets.VBox(), "delfin-path-suggest")
    start_btn = _classed(widgets.Button(description="Start"),
                         "delfin-session-start")
    cancel_btn = _classed(widgets.Button(description="Cancel"),
                          "delfin-session-cancel")
    own_worktree_box = widgets.Checkbox(
        value=False, description=OWN_WORKTREE_LABEL_IN_REPO, indent=False,
        tooltip="Work on a branch of its own in a separate checkout",
        layout=widgets.Layout(display="none"))
    worktree_hint = _classed(widgets.HTML(""), "delfin-session-hint")
    # The consent row for an own worktree in a plain folder (hidden until
    # a start in one asks): the question, Continue and Cancel.
    consent_box = _classed(widgets.VBox(), "delfin-session-consent")
    consent_yes = _classed(widgets.Button(description="Continue",
                                          tooltip=("Make this folder a local "
                                                   "git project and start in "
                                                   "a worktree of it")),
                           "delfin-session-consent-yes")
    consent_no = _classed(widgets.Button(description="Cancel",
                                         tooltip="Leave the folder as it is"),
                          "delfin-session-consent-no")
    consent_yes.on_click(lambda _b: _on_consent_yes())
    consent_no.on_click(lambda _b: _hide_consent())
    new_form = _classed(widgets.VBox(
        [_classed(widgets.HTML("Working directory"), "delfin-session-label"),
         workdir_box, path_suggest, own_worktree_box, worktree_hint,
         consent_box,
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
    # The agent service: every session started here runs in it, and the
    # running list shows every session it holds. Without one, the kernel.
    live: dict[str, Any] = {"client": None, "why": "", "retry_at": 0.0, "scripts": False}
    running_panel = RunningPanel(None, on_open=lambda sid: open_live(sid),
                                 is_open=lambda sid: any(r.get("session_id") == sid
                                                         for r in sessions))
    running_panel.widget.layout.display = "none"
    # The message box's suggestions (_COMPLETE_JS): what is typed up to the
    # caret comes in here, the answer goes out in complete_result.
    complete_query = _classed(widgets.Text(
        continuous_update=True, layout=widgets.Layout(display="none")),
        "delfin-complete-query")
    complete_result = _classed(widgets.HTML(
        "", layout=widgets.Layout(display="none")), "delfin-complete-result")

    def _on_complete(change) -> None:
        try:
            ask = json.loads(change.get("new") or "{}")
        except ValueError:
            return
        if not isinstance(ask, dict):
            return
        try:
            from delfin.dashboard import completion as _completion
            rec = next((r for r in sessions if r["key"] == view.get("active")), None)
            cwd = (rec or {}).get("workspace") or default_workspace(ctx)
            own = []
            if rec is not None and not rec.get("live"):
                for m in _state(rec).get("chat_messages") or []:
                    text = str((m or {}).get("content") or "") if isinstance(m, dict) else ""
                    if (m.get("role") == "user" and text.startswith("!")
                            and "\n" not in text):
                        own.append(text)
            answer = _completion.suggest(str(ask.get("text") or ""), cwd,
                                         own_commands=own)
        except Exception:
            answer = {"start": 0, "items": []}
        answer["seq"] = ask.get("seq")
        complete_result.value = ('<div data-json="%s"></div>'
                                 % html.escape(json.dumps(answer), quote=True))

    complete_query.observe(_on_complete, names="value")

    sidebar = _classed(widgets.VBox(
        [widgets.HTML(_SIDEBAR_CSS), head, new_form, notice, list_box,
         running_panel.widget, resume_dropdown, elsewhere_note, stop_all_btn,
         complete_query, complete_result]),
        "delfin-sessions")
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
    _register_init_js(ctx, _RENAME_KEYS_JS)
    _register_init_js(ctx, _PATH_KEYS_JS)
    _register_init_js(ctx, _COMPLETE_JS)

    # -- lookup -------------------------------------------------------------

    def _find(key: str) -> Optional[dict]:
        return next((rec for rec in sessions if rec["key"] == key), None)

    def _state(rec: dict) -> dict:
        return rec["refs"].get("state") or {}

    def _mark(rec: dict) -> None:
        """One dot per row: what the session needs from you beats that it
        is working, which beats that it finished while you looked away."""
        state = _state(rec)
        if rec.get("live"):
            # A view on a service session: what its events (shown) or the
            # service's summary (hidden) say.
            busy, needs = _live_flags(rec)
        else:
            busy, needs = bool(state.get("streaming")), needs_you(state)
        if rec.get("_was_busy") and not busy and view["active"] != rec["key"]:
            rec["unseen"] = True
        rec["_was_busy"] = busy
        wanted = session_dot(busy=busy, needs=needs,
                             unseen=bool(rec.get("unseen")))
        for cls in _DOT_CLASSES.values():
            if cls == wanted:
                rec["row"].add_class(cls)
            else:
                rec["row"].remove_class(cls)

    def _address_for(session_id: str) -> str:
        """The key this session is known by, to a person and to a peer.

        Input: the session id being resumed, or "". Output: an eight-hex
        address. Semantics: the SAME address every time a given saved
        conversation is opened, and a fresh one for a conversation that
        does not have an id yet.

        It was `uuid.uuid4().hex[:8]` unconditionally -- a new address per
        TAB, not per session. A peer that had been told an address kept
        writing to it after the session was re-opened, and was answered
        "no other open session": observed between three sessions working
        together, which recovered by reading the roster out of the
        refusal and sending again, a round each time. Anything still
        unread in the old mailbox stays there.

        The head of the session id is what `cli._presence_key_for` already
        uses for a terminal session, so the two surfaces now mean the same
        session by the same word instead of holding two schemes for one
        idea.

        A collision with a session already open takes a random address
        rather than two tabs sharing one inbox; eight hex digits make that
        vanishingly rare, and being wrong in the other direction would
        deliver one session's messages to another.
        """
        head = str(session_id or "").strip()[:8]
        if head and not any(rec["key"] == head for rec in sessions):
            return head
        return uuid.uuid4().hex[:8]

    def _session_id(rec: dict) -> str:
        return str(_state(rec).get("active_session_id") or "")

    def rename_session(key: str, name: str) -> str:
        """Give an open session the name *name* (empty: back to its first
        message). Kept with the session: its next save writes it as the
        saved title; a session no process holds is renamed on disk now.
        Returns the name shown."""
        name = " ".join(str(name or "").split())[:80]
        rec = next((r for r in sessions if r["key"] == key), None)
        if rec is None or rec.get("live"):
            return ""
        state = _state(rec)
        sid = _session_id(rec)
        if sid:
            titles = state.setdefault("_messenger_titles", {})
            if name:
                titles[sid] = name
            else:
                titles.pop(sid, None)
            try:
                from delfin.agent import session_store as _store
                _store.rename_session(sid, name)
            except Exception:
                pass
        else:
            rec["chosen_title"] = name
        refresh()
        _refresh_resume_options()
        return rec["title_btn"].description

    view["rename_session"] = rename_session

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
        # A service session's view is attached only while it is on screen:
        # the hidden ones let their sessions sleep.
        for other in sessions:
            if other.get("live") and other["key"] != key:
                other["refs"]["hide"]()
        if rec.get("live") and not rec["refs"]["show"]():
            _say(f"Session {rec['session_id'][:8]} could not be opened: "
                 f"{html.unescape(rec['refs']['notice'].value) or 'the service refused it'}")
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
            {"session_id": sid, "workspace": rec["workspace"],
             **({"live": True} if rec.get("live") else {})}
            for rec in sessions for sid in [_remembered_id(rec)] if sid])

    def _remembered_id(rec: dict) -> str:
        """The conversation this session stands for in the open list.

        Its own id once it has one. Before that, the conversation it was
        asked to open, while a live writer still holds it: that one is
        open in another view, and leaving it out of the list meant the
        next start did not reopen it anywhere."""
        sid = _session_id(rec)
        if sid:
            return sid
        wanted = str(getattr(rec["ctx"], "initial_session_id", "") or "")
        if not wanted:
            return ""
        try:
            from delfin.agent.session_store import session_lock_holder
            return wanted if session_lock_holder(wanted) else ""
        except Exception:
            return ""

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

    # -- the agent service --------------------------------------------------

    def _connect() -> bool:
        """Reach the service (``connector``); True when there is a client.
        A failure is said on the page, never swallowed."""
        if connector is None:
            return False
        import time as _time
        from delfin.dashboard import live_sessions as LS
        live["retry_at"] = _time.monotonic() + _RECONNECT_S
        try:
            client = connector()
        except LS.ServiceUnavailable as exc:
            live["client"], live["why"] = None, str(exc)
            return False
        except Exception as exc:                          # noqa: BLE001
            live["client"], live["why"] = None, f"the agent service could not be reached: {exc}"
            return False
        live["client"], live["why"] = client, ""
        live["had"], live["stopped"] = True, False
        client.on_close(_service_lost)

        def _heard(event: dict) -> None:
            if event.get("type") == "service.stopping" and not event.get("restart"):
                live["stopped"] = True          # stopped on purpose: not restarted from here
        client.subscribe(_heard)
        running_panel.attach(client)
        running_panel.widget.layout.display = ""
        return True

    def _kernel_notice() -> str:
        return (f"The agent service is not available ({live['why']}). Sessions run "
                "in this dashboard's kernel: the terminal and the phone cannot open "
                "them, and they end with this dashboard.")

    def _service_lost() -> None:
        """The service stopped (a restart, an update, a stop): the views
        say so; the refresh tick connects again and reopens them."""
        live["client"] = None
        running_panel.detach()
        if live.get("stopped"):
            live["why"] = "it was stopped"
            text = ("The agent service was stopped; its sessions are saved. Reload "
                    "the page to start it again.")
        else:
            live["why"] = "it went away; connecting again"
            text = "The agent service went away. Connecting again…"
        _say(text)
        for rec in sessions:
            if rec.get("live"):
                rec["refs"]["say"](text)

    def _reconnect() -> None:
        """After a restart or a service that vanished (never after a stop
        given on purpose, nor when there never was one): at most every
        _RECONNECT_S, on the refresh tick."""
        import time as _time
        if connector is None or live["client"] is not None or not live.get("had") \
                or live.get("stopped") or _time.monotonic() < live["retry_at"]:
            return
        if not _connect():
            return
        _say("")
        for rec in [r for r in sessions if r.get("live")]:
            _rebuild_live(rec)

    def _live_tab(sid: str, key: str, attached: Optional[dict] = None):
        from delfin.dashboard import live_tab
        scripts = not live["scripts"]
        live["scripts"] = True
        # Lazy: attached only while it is the view on screen (activate).
        return live_tab.build_live_tab(
            live["client"], sid, attached=attached, page_scripts=scripts, lazy=True,
            on_change=lambda _key=key: _on_live_change(_key))

    def _live_summary(rec: dict) -> dict:
        return running_panel.running.summary(rec.get("session_id") or "") or {}

    def _live_name(rec: dict) -> str:
        if rec["refs"]["attached"]():
            return rec["refs"]["name"]()
        summary = _live_summary(rec)
        return str(summary.get("name") or rec["refs"]["name"]())

    def _live_flags(rec: dict) -> tuple[bool, bool]:
        """(working, needs you) of a view's session: from its events while
        it is shown, from the service's summary while it is hidden."""
        if rec["refs"]["attached"]():
            return bool(rec["refs"]["busy"]()), bool(rec["refs"]["needs"]())
        from delfin.dashboard import live_sessions as LS
        summary = _live_summary(rec)
        word = LS.state_word(summary) if summary else "done"
        return word == "working", word == "needs you"

    def _on_live_change(key: str) -> None:
        rec = _find(key)
        if rec is None:
            return
        name = _live_name(rec)
        if rec["title_btn"].description != name:
            rec["title_btn"].description = name
        _mark(rec)

    def _rebuild_live(rec: dict) -> None:
        """A view on the new connection, in the old one's place."""
        tab, refs = _live_tab(rec["session_id"], rec["key"])
        refs["state"] = rec["refs"]["state"]
        refs["shutdown"] = refs["close"]
        rec["tab"], rec["refs"] = tab, refs
        shown = view["active"] == rec["key"]
        tab.layout.display = "" if shown else "none"
        stage.children = tuple(r["tab"] for r in sessions)
        if shown:
            refs["show"]()

    def open_live(sid: str, workspace: str = "", *, show: bool = True) -> Optional[dict]:
        """A view on service session ``sid`` (one click in the running
        list): its recent turns, then live. Shown again when it is open.
        ``show`` False: listed only (a restored view off screen), not
        attached -- its session is not woken."""
        for rec in sessions:
            if rec.get("session_id") == sid or _session_id(rec) == sid:
                if show:
                    activate(rec["key"])
                return rec
        if live["client"] is None:
            _say("The agent service is not reachable now; try again in a moment.")
            return None
        if len(sessions) >= _MAX_OPEN:
            _say(f"{_MAX_OPEN} sessions are open — close one first.")
            return None
        summary = running_panel.running.summary(sid) or {}
        where = str(workspace or summary.get("workspace") or "")
        key = _address_for(sid)
        tab, refs = _live_tab(sid, key)
        if show and not refs["show"]():
            _say(f"Session {sid[:8]} could not be opened: "
                 f"{html.unescape(refs['notice'].value) or 'the service refused it'}")
            refs["close"]()
            return None
        where = where or str(refs["transcript"].info.get("workspace") or "")
        refs["state"] = {"active_session_id": sid, "chat_messages": [],
                         "streaming": False, "engine": None}
        refs["shutdown"] = refs["close"]
        status = widgets.HTML("")
        title_btn = _classed(widgets.Button(
            description=str(summary.get("name") or refs["name"]()), tooltip=where),
            "delfin-session-title")
        close_btn = _classed(widgets.Button(description="×", tooltip=(
            "Close this view -- the session keeps running in the agent service; "
            "the terminal and the phone still have it")), "delfin-session-close")
        info = _classed(widgets.HTML(html.escape(_display_path(where))), "delfin-session-path")
        row = _classed(widgets.HBox([
            _classed(widgets.VBox([title_btn, info]), "delfin-session-text"),
            close_btn]), "delfin-session-item")
        rec = {"key": key, "workspace": where, "tab": tab, "refs": refs, "live": True,
               "session_id": sid,
               "ctx": _SessionContext(ctx, agent_workspace=where, initial_session_id=sid,
                                      sessions_managed=True, session_open_elsewhere=None,
                                      agent_engine=None, agent_state=refs["state"],
                                      agent_status_html=status,
                                      add_init_js=lambda script: None, presence_key=key),
               "title_btn": title_btn, "info": info, "row": row}
        title_btn.on_click(lambda _b, _key=key: activate(_key))
        close_btn.on_click(lambda _b, _key=key: close_session(_key))
        sessions.append(rec)
        stage.children = tuple(r["tab"] for r in sessions)
        if show:
            activate(key)
        else:
            tab.layout.display = "none"
            _render_list()
        running_panel.render()
        _persist()
        return rec

    def _open_in_service(where: str, resume: str = "") -> Optional[dict]:
        """A new session in the service (or ``resume``: a saved one, woken
        there), then a view on it. None when the service refused it -- the
        caller then runs it in the kernel and the page says why."""
        from delfin.dashboard import live_sessions as LS
        client = live["client"]
        if resume and running_panel.running.summary(resume):
            return open_live(resume, where)        # it runs there already
        settings = dict(session_defaults() or {})
        bypass = bool(settings.pop("_bypass", False))
        if str(settings.get("mode") or "") == "dashboard":
            live["note"] = ("Dashboard mode works on this page's own tabs: this session "
                            "runs in this dashboard's kernel, and the terminal and the "
                            "phone cannot open it.")
            return None
        try:
            started = client.create(where, resume=resume or None, **settings)
        except (LS.ServiceError, TimeoutError, OSError) as exc:
            live["note"] = (f"The agent service did not start the session "
                            f"({getattr(exc, 'message', exc)}); it runs in this "
                            "dashboard's kernel instead.")
            return None
        sid = str(started.get("sessionId") or "")
        rec = open_live(sid, where)
        if rec is not None and bypass:
            try:
                client.posture(sid, "bypassPermissions", confirmed=True)
            except LS.ServiceError as exc:
                rec["refs"]["say"](f"Bypass is not available here ({exc.message}); "
                                   "this session asks as Accept Edits.")
        return rec

    # -- open, close ----------------------------------------------------------

    def open_session(workspace: str = "", session_id: str = "", *,
                     in_kernel: bool = False) -> Optional[dict]:
        """Open a session working in ``workspace``; with ``session_id``, the
        saved conversation. A conversation already open is shown instead.
        In the agent service when there is one (``in_kernel``: never)."""
        sid = str(session_id or "").strip()
        if sid:
            for rec in sessions:
                if _session_id(rec) == sid:
                    activate(rec["key"])
                    return rec
        if len(sessions) >= _MAX_OPEN:
            _say(f"{_MAX_OPEN} sessions are open — close one first.")
            return None
        in_kernel = in_kernel or not view["start_in_service"]
        if live["client"] is not None and not in_kernel:
            where = str(workspace or "").strip() or default_workspace(ctx)
            live["note"] = ""
            rec = _open_in_service(where, resume=sid)
            if rec is not None:
                return rec
            _say(live["note"] or _kernel_notice())
        elif connector is not None and not in_kernel:
            _say(_kernel_notice())
        else:
            _say("")
        key = _address_for(sid)
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
        # A session in its own worktree gets a "Merge back" action: its
        # committed and uncommitted changes go to the main checkout,
        # uncommitted, for review. A conflict is said, not overwritten --
        # the worktree stays so the work is not lost.
        merge_btn = None
        if read_worktree_sidecar(where) is not None:
            merge_btn = _classed(widgets.Button(
                description="Merge back",
                tooltip=("Bring this session's changes back to the main "
                         "checkout, uncommitted, for review; a conflict is "
                         "shown, nothing overwritten")),
                "delfin-session-merge")
            merge_btn.on_click(
                lambda _b, _where=where: _merge_back(_where))
        # Rename in place: the pencil opens a box over the name; Enter or
        # Save keeps it, Esc or Cancel leaves it (the running rows do the
        # same, so both lists read alike).
        pencil = _classed(widgets.Button(
            description="✎", tooltip="Rename (Enter saves, Esc cancels)"),
            "delfin-session-close")
        pencil.add_class("delfin-session-rename-btn")
        name_edit = widgets.Text(continuous_update=False,
                                 placeholder="Name (empty: its first message)")
        name_edit.add_class("delfin-session-rename-input")
        name_save = widgets.Button(description="Save")
        name_cancel = _classed(widgets.Button(description="Cancel"),
                               "delfin-live-rename-cancel")
        rename_box = _classed(widgets.HBox(
            [name_edit, name_save, name_cancel],
            layout=widgets.Layout(display="none")), "delfin-running-rename")
        rename_box.add_class("delfin-live-rename")
        row = _classed(widgets.HBox([
            _classed(widgets.VBox([title_btn, info, rename_box]),
                     "delfin-session-text"),
            *([merge_btn] if merge_btn else []), pencil, close_btn]),
            "delfin-session-item")
        rec = {"key": key, "workspace": where, "tab": tab, "refs": refs or {},
               "ctx": session_ctx, "title_btn": title_btn, "info": info,
               "row": row, "merge_btn": merge_btn}

        def _open_rename(_b=None):
            name_edit.value = title_btn.description \
                if title_btn.description != "New session" else ""
            rename_box.layout.display = ""

        def _close_rename(_b=None):
            rename_box.layout.display = "none"

        def _save_rename(_b=None, _key=key):
            rename_session(_key, name_edit.value)
            _close_rename()

        pencil.on_click(_open_rename)
        name_save.on_click(_save_rename)
        name_cancel.on_click(_close_rename)
        name_edit.on_submit(_save_rename)
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

    def _merge_back(where: str) -> None:
        """Merge the session worktree ``where`` is in back into its main
        checkout: said in the notice, conflicts included; a conflict leaves
        the main checkout as it was and keeps the worktree."""
        try:
            out = merge_session_worktree(where)
        except Exception as exc:
            _say(f"Merge back failed: {exc}")
            return
        if out["status"] == "ok":
            files = ", ".join(out.get("files") or []) or "no file changes"
            _say(f"Merged back to {out.get('target')}: {files}. "
                 "Uncommitted, for review.")
        elif out["status"] == "conflict":
            files = ", ".join(out.get("files") or []) or "the changed files"
            _say(f"Conflict, nothing overwritten: {files}. "
                 f"{out.get('message') or ''} The worktree is kept; "
                 "resolve it there and merge again.")
        else:
            _say(f"Merge back did not apply: {out.get('message') or 'error'}")

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
        if rec.get("live"):
            # Only this view goes: the session runs on in the service, its
            # worktree and its presence record are the worker's.
            sessions.remove(rec)
            stage.children = tuple(r["tab"] for r in sessions)
            if sessions and view["active"] == key:
                activate(sessions[-1]["key"])
            elif not sessions:
                view["active"] = ""
            _render_list()
            running_panel.render()
            _persist()
            _refresh_resume_options()
            return
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
        _reconnect()
        if live["client"] is not None:
            running_panel.render()        # the ages; the rows follow events
        for rec in sessions:
            if rec.get("live"):
                # The service's session: named, announced and messaged by
                # its worker; its view follows the service's events.
                _on_live_change(rec["key"])
                continue
            state = _state(rec)
            sid_now = _session_id(rec)
            if rec.get("chosen_title") and sid_now:
                # A name given before the session had an id goes where the
                # save reads it as soon as there is one.
                state.setdefault("_messenger_titles", {})[sid_now] = \
                    rec.pop("chosen_title")
            label = rec.get("chosen_title") or _session_title(state, sid_now)
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
                                   title=_session_title(state, sid),
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

    # The loop's state. ``tick_lock`` is held for a whole tick and by the
    # stop, so a stop either precedes a tick (which then does nothing) or
    # waits for it, and no tick reschedules after a stop. Reentrant: a stop
    # called from inside a tick (a refresh that closes the panel) must not
    # wait for itself.
    loop = {"stopped": False, "timer": None,
            "tick_lock": threading.RLock()}

    def _tick() -> None:
        with loop["tick_lock"]:
            if loop["stopped"]:
                return
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
                if not loop["stopped"]:
                    timer = threading.Timer(_REFRESH_S, _tick)
                    timer.daemon = True
                    timer.start()
                    loop["timer"] = view["timer"] = timer

    def stop_refresh() -> None:
        """End this list's refresh loop: no tick runs after it returns.
        Idempotent; called from inside a tick it marks the stop and the
        tick does not reschedule."""
        loop["stopped"] = True
        running_panel.detach()
        timer = loop["timer"]
        if timer is not None:
            timer.cancel()
        with loop["tick_lock"]:
            pass
        # A timer that fired before the cancel returns at once (stopped);
        # joined so no thread of this loop is left when the stop returns.
        timer = loop["timer"]
        if timer is not None and timer is not threading.current_thread():
            timer.join(timeout=5.0)
        with _REFRESH_LOOPS_LOCK:
            if _REFRESH_LOOPS.get(id(ctx)) is stop_refresh:
                del _REFRESH_LOOPS[id(ctx)]

    # -- controls ---------------------------------------------------------------

    def _suggest_worktree(_change=None) -> None:
        """Offer a worktree inside a git repository, and pre-select it when
        another open session already works in that repository."""
        where = str(workdir_box.value or "").strip()
        # A held session worktree is left out of workspace_choices; typed
        # here, it is named with the reason, so the missing entry is not a
        # silent one.
        try:
            held = session_worktree_state(where)
        except Exception:
            held = {"is_worktree": False, "holder": "", "why": ""}
        if held.get("is_worktree") and held.get("why"):
            own_worktree_box.layout.display = "none"
            own_worktree_box.value = False
            worktree_hint.value = (
                "<div style='font-size:10px;color:#546e7a'>Not offered as a "
                "start folder: this session worktree is held -- "
                f"{html.escape(held['why'])}. Choose another directory.</div>")
            return
        try:
            from delfin.agent import session_presence as _presence
            in_repo = bool(where) and bool(_presence.repository_of(where)["root"])
            busy = _presence.in_same_repository(where) if in_repo else []
        except Exception:
            in_repo, busy = False, []
        own_worktree_box.layout.display = "" if (in_repo or where) else "none"
        own_worktree_box.value = bool(busy)
        own_worktree_box.description = (OWN_WORKTREE_LABEL_IN_REPO if in_repo
                                        else OWN_WORKTREE_PLAIN_LABEL)
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
        _refresh_folder_choices()
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
        plain = plain_workspace(target or default_workspace(ctx)) \
            if (own_worktree_box.value
                and own_worktree_box.layout.display != "none") else False
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
            if plain:
                # A plain folder: making it a repository needs consent.
                # ipywidgets has no blocking dialog, so the start stops at
                # a visible confirm row in the form (question, Continue,
                # Cancel); the Start continues only after the yes, so
                # nothing is done before the consent -- the same question
                # the engine's own window asks, asked once per start.
                _say("")
                from delfin.agent import local_project as _LP
                consent_box.children = (widgets.HTML(
                    "<div style='font-size:11px'>"
                    f"{html.escape(_LP.consent_question(_LP.plan(target)))}"
                    "</div>"), consent_yes, consent_no)
                consent_box.layout.display = ""
                start_btn.layout.display = "none"
                return
            try:
                target = session_worktree(target or default_workspace(ctx))
            except Exception as exc:
                _say(f"No worktree was created: {exc}")
                return
        _on_cancel()
        open_session(target)

    def _on_consent_yes(_btn=None) -> None:
        """Continue after the consent row: hide it, make the folder a
        local project and start in a worktree of it."""
        _hide_consent()
        where = str(workdir_box.value or "").strip()
        target = str(Path(where).expanduser()) if where \
            else default_workspace(ctx)
        try:
            made = start_session_worktree(target, ask=lambda _q: True)
            target = made["path"]
        except Exception as exc:
            _say(f"No worktree was created: {exc}")
            return
        _on_cancel()
        open_session(target)

    def _hide_consent(_btn=None) -> None:
        consent_box.children = ()
        consent_box.layout.display = "none"
        start_btn.layout.display = ""

    def _on_resume(change) -> None:
        sid = str(change.get("new") or "")
        if not sid:
            return
        data = _saved_session(sid) or {}
        open_session(str(data.get("workspace") or ""), sid)
        _refresh_resume_options()

    workdir_box.observe(_suggest_worktree, names="value")

    def _pick_folder(path: str) -> None:
        """Go into *path*: the box shows it with a trailing slash, so the
        list below becomes its sub-folders (as Tab does in a shell)."""
        workdir_box.value = str(path).rstrip("/") + "/"

    def _show_suggestions(choices: list, typed: str) -> None:
        typed_dir = typed if typed.endswith("/") else \
            (typed.rsplit("/", 1)[0] + "/" if "/" in typed else "")
        typed_dir = os.path.expanduser(typed_dir) if typed_dir else ""
        # What fits the typed path first: Tab takes the first entry, and
        # it must complete what was typed, not jump to the default folder.
        full = os.path.expanduser(typed) if typed else ""
        fits = [c for c in choices if full and str(c).startswith(full)]
        choices = fits + [c for c in choices if c not in fits]
        buttons = []
        for path in choices[:_SUGGEST_SHOWN]:
            if typed_dir and str(path).startswith(typed_dir):
                label = str(path)[len(typed_dir):] + "/"
            else:
                label = _display_path(str(path))
            b = _classed(widgets.Button(description=label, tooltip=str(path)),
                         "delfin-path-option")
            b.on_click(lambda _b, _p=path: _pick_folder(_p))
            buttons.append(b)
        path_suggest.children = tuple(buttons)

    def _refresh_folder_choices(change=None) -> None:
        """The list under the box follows what is typed, like a shell's Tab:
        the folders of the typed path, or of its parent."""
        try:
            typed = str((change or {}).get("new") if change else workdir_box.value or "")
            choices = workspace_choices(ctx, typed)
            key = (typed, tuple(choices))
            if view.get("suggested") != key:
                view["suggested"] = key
                _show_suggestions(choices, typed)
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
    _connect()
    for row in rows[:_MAX_OPEN]:
        if live["client"] is not None and running_panel.running.summary(row["session_id"]):
            # The service runs it (or keeps it asleep): a view on it,
            # attached only once it is the one on screen (below) -- a
            # reload must not wake every remembered session.
            try:
                open_live(row["session_id"], row["workspace"], show=False)
            except Exception:
                pass
            continue
        data = _saved_session(row["session_id"])
        if not data:
            continue                    # deleted since: nothing to reopen
        try:
            # Saved: resumed in the service when there is one, else here.
            open_session(row["workspace"] or str(data.get("workspace") or ""),
                         row["session_id"])
        except Exception:
            continue
    if not sessions:
        open_session()
    if sessions and _find(view["active"]) is None:
        activate(sessions[-1]["key"])
    if view["start_in_service"] and live["client"] is None:
        _say(_kernel_notice())
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
    with _REFRESH_LOOPS_LOCK:
        previous = _REFRESH_LOOPS.pop(id(ctx), None)
    if previous is not None:
        previous()
    with _REFRESH_LOOPS_LOCK:
        _REFRESH_LOOPS[id(ctx)] = stop_refresh
    _tick()

    return widget, {
        "open": open_session,
        "close": close_session,
        "activate": activate,
        "sessions": lambda: list(sessions),
        "active": lambda: _find(view["active"]),
        "refresh": refresh,
        "stop_refresh": stop_refresh,
        "form": {"new": new_btn, "workdir": workdir_box,
                 "own_worktree": own_worktree_box, "start": start_btn,
                 "consent_yes": consent_yes, "consent_no": consent_no},
        # Exported so the reclamation can be driven in a test; a pass that
        # removes directories should not be reachable only at startup.
        "reclaim_worktrees": reclaim_orphaned_worktrees,
        "notice": notice,
        "service": lambda: live["client"],
        "running": running_panel.running,
        "running_panel": running_panel,
        "open_live": open_live,
    }
