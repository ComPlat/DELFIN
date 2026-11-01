"""Shared pieces of the dashboard header: messages, per-user state, style, keys.

Input: text for the person (and optional details and actions), the header's
widgets. Output:

* :class:`Toasts` -- one message area (bottom right, fixed; the top right
  is where the bar's own controls end): the latest
  message only, a dismiss action on it, success and info messages expire
  after ``auto_s`` seconds, errors and questions wait; earlier messages on
  demand (the last ``keep``). Details (paths, ports, commands) sit in a
  collapsed ``<details>`` under the text.
* :func:`load_prefs` / :func:`save_prefs` -- the header's per-user state in
  ``~/.delfin/dashboard_header.json`` (more row open, status hint compact,
  last view); best effort: an unreadable or unwritable file means
  defaults, never an error.
* :data:`HEADER_CSS` -- the style of every header element, scoped to the
  ``delfin-hdr`` class and written against the JupyterLab theme variables
  (so the Voila dark theme applies) with light fallbacks. Visibility is
  never set here: widgets show and hide through ``layout.display``, which
  a rule with ``!important`` would override.
* :func:`keys_js` -- Esc presses the Cancel of an open confirmation,
  else closes the open menu, ``/`` or Ctrl+K opens the
  branch menu and focuses its search, a click outside closes it, opening
  it focuses the search.

Rejected: messages printed into an Output widget under the bar (they
accumulated until a reload, the field complaint this replaces); a
JavaScript dependency for the menu (the ipywidgets stack plus one inline
script is enough).
"""

from __future__ import annotations

import html
import json
import threading
import time
from pathlib import Path
from typing import Callable, Optional

# -- per-user state ------------------------------------------------------------


def _prefs_path() -> Path:
    """Resolved per call, so a redirected home moves it."""
    return Path.home() / ".delfin" / "dashboard_header.json"


def load_prefs() -> dict:
    try:
        data = json.loads(_prefs_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def save_prefs(**changes) -> None:
    try:
        from delfin.agent.state_paths import write_text_atomic
        data = load_prefs()
        data.update(changes)
        write_text_atomic(_prefs_path(), json.dumps(data, indent=1, sort_keys=True))
    except Exception:
        pass


def _esc(text) -> str:
    return html.escape(str(text or ""))


# -- messages ------------------------------------------------------------------

_ICON = {"ok": "✓", "info": "ℹ", "warn": "!", "error": "✕", "busy": "…",
         "ask": "?"}


class Toasts:
    """The header's one message area (see the module docstring)."""

    def __init__(self, *, auto_s: float = 6.0, keep: int = 20, timer: bool = True):
        import ipywidgets as W
        self.W = W
        self.auto_s = float(auto_s)
        self.keep = int(keep)
        self.timer = timer
        self.history: list[dict] = []
        self._current: Optional[dict] = None
        self._seq = 0
        self._lock = threading.RLock()
        self.box = W.VBox([])
        self.box.add_class("delfin-toasts")
        self.earlier = W.VBox([])
        self.earlier.add_class("delfin-toast-earlier")
        self.earlier.layout.display = "none"

    # the latest message, as data (what the tests and the page read)
    def latest(self) -> Optional[dict]:
        return self._current

    def show(self, text: str, *, level: str = "info", details: str = "",
             actions: Optional[list] = None, key: str = "") -> int:
        """Show *text* as the latest message. ``actions``: ``[(label, fn)]``.
        A message with the same *key* as the shown one replaces it."""
        with self._lock:
            self._seq += 1
            msg = {"id": self._seq, "text": str(text or ""), "level": level,
                   "details": str(details or ""), "actions": list(actions or []),
                   "key": key, "at": time.time()}
            if (self._current is not None and key and self._current.get("key") == key
                    and self.history and self.history[-1] is self._current):
                self.history[-1] = msg
            else:
                self.history.append(msg)
            del self.history[:-self.keep]
            self._current = msg
            self._render()
        if level in ("ok", "info") and self.timer and self.auto_s > 0:
            t = threading.Timer(self.auto_s, self.expire, kwargs={"msg_id": msg["id"]})
            t.daemon = True
            t.start()
        return msg["id"]

    def dismiss(self, msg_id: Optional[int] = None) -> None:
        with self._lock:
            if self._current is None:
                return
            if msg_id is not None and self._current["id"] != msg_id:
                return
            self._current = None
            self._render()

    def expire(self, msg_id: Optional[int] = None) -> None:
        """Drop the shown message when it is a success or an info one."""
        with self._lock:
            cur = self._current
            if cur is None or cur["level"] not in ("ok", "info"):
                return
            if msg_id is not None and cur["id"] != msg_id:
                return
            self._current = None
            self._render()

    def run_action(self, index: int) -> None:
        cur = self._current
        if cur is None or index >= len(cur["actions"]):
            return
        _label, fn = cur["actions"][index]
        self.dismiss(cur["id"])
        try:
            fn()
        except Exception as exc:
            self.show(f"That did not work: {exc}", level="error")

    def toggle_earlier(self, *_a) -> None:
        open_ = self.earlier.layout.display == "none"
        self.earlier.layout.display = "" if open_ else "none"
        self._render()

    def _render(self) -> None:
        W = self.W
        cur = self._current
        if cur is None:
            self.box.children = ((self._earlier_btn(), self.earlier)
                                 if self.earlier.layout.display != "none" else ())
            return
        lvl = cur["level"]
        body = (f'<span class="delfin-toast-icon">{_ICON.get(lvl, "")}</span>'
                f'<span class="delfin-toast-text">{_esc(cur["text"])}</span>')
        if cur["details"]:
            body += ('<details class="delfin-toast-details"><summary>Details</summary>'
                     f'<pre>{_esc(cur["details"])}</pre></details>')
        text = W.HTML(body)
        text.add_class("delfin-toast-body")
        buttons = []
        for i, (label, _fn) in enumerate(cur["actions"]):
            b = W.Button(description=str(label))
            b.add_class("delfin-btn")
            b.add_class("delfin-toast-action")
            b.on_click(lambda _b, i=i: self.run_action(i))
            buttons.append(b)
        close = W.Button(description="✕", tooltip="Dismiss")
        close.add_class("delfin-icon-btn")
        close.add_class("delfin-toast-close")
        close.on_click(lambda _b, mid=cur["id"]: self.dismiss(mid))
        row = W.HBox([text, *buttons, close])
        row.add_class("delfin-toast")
        row.add_class(f"is-{lvl}")
        kids = [row]
        n = len(self.history) - 1
        if n > 0:
            kids.append(self._earlier_btn())
        if self.earlier.layout.display != "none":
            kids.append(self.earlier)
        self.box.children = tuple(kids)

    def _earlier_btn(self):
        W = self.W
        n = max(0, len(self.history) - (1 if self._current else 0))
        open_ = self.earlier.layout.display != "none"
        b = W.Button(description=("Hide earlier" if open_ else f"Earlier ({n})"),
                     tooltip="Earlier messages")
        b.add_class("delfin-link-btn")
        b.on_click(self.toggle_earlier)
        if open_:
            items = []
            for m in reversed(self.history):
                if m is self._current:
                    continue
                when = time.strftime("%H:%M", time.localtime(m["at"]))
                items.append(f'<li class="is-{_esc(m["level"])}"><span>{when}</span> '
                             f'{_esc(m["text"])}</li>')
            self.earlier.children = (W.HTML("<ul>" + "".join(items) + "</ul>"),)
        return b


# -- style ---------------------------------------------------------------------

HEADER_CSS = """
<style>
.delfin-hdr {
  --dh-bg: var(--jp-layout-color1, #ffffff);
  --dh-bg2: var(--jp-layout-color2, #f3f5f8);
  --dh-fg: var(--jp-ui-font-color1, #1f2933);
  --dh-muted: var(--jp-ui-font-color2, #5f6b7a);
  --dh-border: var(--jp-border-color1, #cfd6df);
  --dh-accent: #3d74b6;
  --dh-accent-fg: #ffffff;
  --dh-ok: #2f8a4c; --dh-warn: #b7791f; --dh-err: #c43c3c;
  --dh-h: 30px; --dh-r: 8px;
  overflow: visible;
}
.delfin-hdr .widget-box, .delfin-hdr .widget-hbox, .delfin-hdr .widget-vbox {
  overflow: visible;
}
/* The branch menu hangs below the bar: the cell and output areas around the
   header must not clip it, and the header's cell stacks above the tabs. */
.jp-Cell:has(.delfin-hdr), .jp-Cell-outputWrapper:has(.delfin-hdr),
.jp-OutputArea:has(.delfin-hdr), .jp-OutputArea-child:has(.delfin-hdr),
.jp-OutputArea-output:has(.delfin-hdr), .widget-output:has(.delfin-hdr) {
  overflow: visible;
}
/* Above page content (in-page panels use up to 60; the molecule viewer's
   buttons 30), below full-screen overlays (9990 and up). */
.jp-Cell:has(.delfin-hdr), .jp-OutputArea-child:has(.delfin-hdr) {
  position: relative; z-index: 1000;
}
.delfin-topbar { flex-wrap: wrap; align-items: center; gap: 8px; row-gap: 6px; }
.delfin-hdr .delfin-group { align-items: center; gap: 6px; flex-wrap: wrap; }
/* buttons: one height, rounded, no wrapped labels, auto width */
.delfin-hdr .jupyter-button.delfin-btn,
.delfin-hdr .jupyter-button.delfin-icon-btn,
.delfin-hdr .jupyter-button.delfin-link-btn {
  height: var(--dh-h); line-height: calc(var(--dh-h) - 2px);
  width: auto; min-width: 0; padding: 0 12px; margin: 0;
  border-radius: var(--dh-r); border: 1px solid var(--dh-border);
  background: var(--dh-bg); color: var(--dh-fg);
  white-space: nowrap; box-shadow: none; font-size: 13px;
  flex: 0 0 auto;
}
.delfin-hdr .jupyter-button.delfin-btn:hover,
.delfin-hdr .jupyter-button.delfin-icon-btn:hover { background: var(--dh-bg2); }
.delfin-hdr .jupyter-button.delfin-btn.is-primary {
  background: var(--dh-accent); border-color: var(--dh-accent); color: var(--dh-accent-fg);
  font-weight: 600;
}
.delfin-hdr .jupyter-button.delfin-btn.is-primary:hover { filter: brightness(1.08); }
.delfin-hdr .jupyter-button.delfin-btn.is-danger { color: var(--dh-err); border-color: var(--dh-err); }
.delfin-hdr .jupyter-button.delfin-btn.is-danger.is-filled { background: var(--dh-err); color: #fff; }
.delfin-hdr .jupyter-button.delfin-icon-btn { padding: 0 8px; border-color: transparent; background: transparent; }
.delfin-hdr .jupyter-button.delfin-link-btn {
  border-color: transparent; background: transparent; color: var(--dh-accent); padding: 0 6px;
}
.delfin-hdr .jupyter-button.mod-active.delfin-btn { background: var(--dh-bg2); border-color: var(--dh-accent); }
.delfin-hdr .jupyter-button:disabled { opacity: .5; }
.delfin-hdr .jupyter-button:focus-visible { outline: 2px solid var(--dh-accent); outline-offset: 1px; }
/* the branch switcher: a pill that opens a menu */
.delfin-branch-wrap { position: relative; flex: 0 0 auto; }
.delfin-hdr .jupyter-button.delfin-branch-btn {
  font-family: var(--jp-code-font-family, monospace); max-width: 280px;
  overflow: hidden; text-overflow: ellipsis; border-radius: 999px;
}
.delfin-hdr .delfin-branch-pop {
  position: absolute; top: calc(var(--dh-h) + 6px); left: 0; z-index: 60;
  width: min(560px, calc(100vw - 32px)); max-height: 70vh;
  background: var(--dh-bg); color: var(--dh-fg);
  border: 1px solid var(--dh-border); border-radius: 12px;
  box-shadow: 0 10px 30px rgba(0,0,0,.18); padding: 8px; gap: 6px;
  display: flex; flex-direction: column; overflow: hidden;
}
.delfin-hdr .delfin-pop-head { gap: 6px; align-items: center; flex-wrap: wrap; }
.delfin-hdr .delfin-branch-pop > :not(.delfin-branch-rows) { flex: 0 0 auto; }
.delfin-hdr .delfin-branch-rows > * { flex-shrink: 0; }
.delfin-hdr .delfin-branch-search { flex: 1 1 200px; min-width: 160px; width: auto; margin: 0; }
.delfin-hdr .widget-text input, .delfin-hdr .widget-text input[type="text"] {
  height: var(--dh-h); border-radius: var(--dh-r); border: 1px solid var(--dh-border);
  background: var(--dh-bg); color: var(--dh-fg); padding: 0 10px;
}
.delfin-hdr .widget-text input:focus { border-color: var(--dh-accent); outline: none; box-shadow: 0 0 0 2px rgba(61,116,182,.25); }
.delfin-hdr .delfin-views { width: auto; margin: 0; }
.delfin-hdr .delfin-views .widget-toggle-button {
  width: auto; height: 26px; line-height: 24px; padding: 0 10px; margin: 0 4px 0 0;
  border-radius: 999px; border: 1px solid var(--dh-border); background: transparent;
  color: var(--dh-muted); font-size: 12px; box-shadow: none;
}
.delfin-hdr .delfin-views .widget-toggle-button.mod-active {
  background: var(--dh-accent); border-color: var(--dh-accent); color: var(--dh-accent-fg);
}
.delfin-hdr .delfin-branch-rows { overflow-y: auto; max-height: calc(70vh - 110px);
  min-height: 0; flex: 1 1 auto; gap: 1px; }
.delfin-hdr .delfin-brow { align-items: center; gap: 4px; border-radius: var(--dh-r); padding: 0 4px; min-height: 32px; }
.delfin-hdr .delfin-brow:hover, .delfin-hdr .delfin-brow.is-selected { background: var(--dh-bg2); }
.delfin-hdr .jupyter-button.delfin-brow-main {
  flex: 1 1 auto; text-align: left; justify-content: flex-start; border-color: transparent;
  background: transparent; font-family: var(--jp-code-font-family, monospace);
  overflow: hidden; text-overflow: ellipsis; padding: 0 6px;
}
.delfin-hdr .delfin-brow-meta { color: var(--dh-muted); font-size: 11px; white-space: nowrap; flex: 0 1 auto; overflow: hidden; }
.delfin-hdr .delfin-brow-meta .ok { color: var(--dh-ok); } .delfin-hdr .delfin-brow-meta .bad { color: var(--dh-err); }
.delfin-hdr .delfin-brow .delfin-brow-act { visibility: hidden; }
.delfin-hdr .delfin-brow:hover .delfin-brow-act, .delfin-hdr .delfin-brow.is-selected .delfin-brow-act,
.delfin-hdr .delfin-brow .delfin-brow-act:focus { visibility: visible; }
.delfin-hdr .delfin-pop-note { color: var(--dh-muted); font-size: 11px; }
.delfin-hdr .delfin-confirm { gap: 6px; align-items: center; flex-wrap: wrap; padding: 6px; border-radius: var(--dh-r); background: var(--dh-bg2); }
.delfin-hdr .delfin-confirm-wide { max-width: 900px; margin: 4px 0; border-left: 4px solid var(--dh-warn); }
.delfin-hdr .delfin-confirm-wide .widget-html { flex: 1 1 100%; }
.delfin-hdr .delfin-confirm-wide .widget-html-content { white-space: normal; line-height: 1.5; }
/* chips: preview, status, home */
.delfin-hdr .delfin-chip {
  align-items: center; gap: 2px; height: var(--dh-h); padding: 0 4px 0 10px;
  border-radius: 999px; border: 1px solid var(--dh-border); background: var(--dh-bg2);
  color: var(--dh-fg); font-size: 12px; white-space: nowrap; flex: 0 0 auto;
}
.delfin-hdr .delfin-chip a { color: var(--dh-accent); font-weight: 600; text-decoration: none; margin-left: 6px; }
.delfin-hdr .delfin-chip .dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; background: var(--dh-ok); }
.delfin-hdr .delfin-chip.is-busy .dot { background: var(--dh-warn); animation: delfinPulse 1s ease-in-out infinite; }
.delfin-hdr .delfin-chip.is-error .dot { background: var(--dh-err); }
@keyframes delfinPulse { 50% { opacity: .3; } }
.delfin-hdr .delfin-status { color: var(--dh-muted); font-size: 12px; white-space: nowrap; }
.delfin-hdr .jupyter-button.delfin-hint {
  border-radius: 999px; font-size: 12px; color: var(--dh-muted); background: transparent;
}
.delfin-hdr .jupyter-button.delfin-hint.is-behind { color: var(--dh-accent); border-color: var(--dh-accent); font-weight: 600; }
.delfin-hdr .jupyter-button.delfin-hint.is-warn { color: var(--dh-err); border-color: var(--dh-err); }
.delfin-hdr .delfin-home { font-size: 12px; color: var(--dh-muted); white-space: nowrap; }
.delfin-hdr .delfin-home.is-warn { color: var(--dh-err); font-weight: 600; }
/* switch (Keep session) */
.delfin-hdr .jupyter-button.delfin-switch {
  height: var(--dh-h); width: auto; padding: 0 10px 0 40px; position: relative;
  border-radius: 999px; border: 1px solid var(--dh-border); background: var(--dh-bg);
  color: var(--dh-fg); white-space: nowrap; box-shadow: none; font-size: 13px; flex: 0 0 auto;
}
.delfin-hdr .jupyter-button.delfin-switch::before {
  content: ""; position: absolute; left: 8px; top: 50%; width: 26px; height: 14px;
  margin-top: -7px; border-radius: 999px; background: var(--dh-border); transition: background .15s;
}
.delfin-hdr .jupyter-button.delfin-switch::after {
  content: ""; position: absolute; left: 10px; top: 50%; width: 10px; height: 10px;
  margin-top: -5px; border-radius: 50%; background: #fff; transition: left .15s;
}
.delfin-hdr .jupyter-button.delfin-switch.mod-active { background: var(--dh-bg); }
.delfin-hdr .jupyter-button.delfin-switch.mod-active::before { background: var(--dh-ok); }
.delfin-hdr .jupyter-button.delfin-switch.mod-active::after { left: 22px; }
.delfin-hdr .jupyter-button.delfin-switch i { display: none; }
/* the more row and panels under the bar */
.delfin-hdr .delfin-more-row { gap: 6px; align-items: center; flex-wrap: wrap; padding: 6px 0 2px 0; }
.delfin-hdr .delfin-panel {
  border: 1px solid var(--dh-border); border-radius: 12px; padding: 10px; gap: 8px;
  background: var(--dh-bg); color: var(--dh-fg); margin: 6px 0; max-width: 900px;
}
.delfin-hdr .delfin-panel-title { font-weight: 600; }
.delfin-hdr .delfin-cleanup-rows { max-height: 40vh; overflow-y: auto; gap: 0; }
.delfin-hdr .delfin-cleanup-rows .widget-checkbox { width: auto; }
.delfin-hdr .delfin-cleanup-rows .widget-checkbox label { white-space: normal; }
.delfin-hdr .delfin-count { color: var(--dh-muted); font-size: 12px; white-space: nowrap; }
.delfin-hdr .delfin-panel .widget-inline-hbox input[type="number"] { min-width: 56px; }
/* the message area */
.delfin-toasts {
  position: fixed; bottom: 16px; right: 16px; z-index: 1000; gap: 4px;
  max-width: min(460px, calc(100vw - 32px)); align-items: flex-end;
  pointer-events: none;
}
.delfin-toasts > * { pointer-events: auto; }
.delfin-hdr .delfin-toast {
  align-items: flex-start; gap: 6px; padding: 8px 6px 8px 12px; border-radius: 10px;
  background: var(--dh-bg); color: var(--dh-fg); border: 1px solid var(--dh-border);
  border-left: 4px solid var(--dh-accent); box-shadow: 0 8px 24px rgba(0,0,0,.18);
  max-width: 100%;
}
.delfin-hdr .delfin-toast.is-ok { border-left-color: var(--dh-ok); animation: delfinFade .4s ease 5.6s forwards; }
.delfin-hdr .delfin-toast.is-info { animation: delfinFade .4s ease 5.6s forwards; }
.delfin-hdr .delfin-toast.is-error { border-left-color: var(--dh-err); }
.delfin-hdr .delfin-toast.is-warn, .delfin-hdr .delfin-toast.is-ask { border-left-color: var(--dh-warn); }
@keyframes delfinFade { to { opacity: 0; } }
.delfin-hdr .delfin-toast-body { flex: 1 1 auto; min-width: 0; }
.delfin-hdr .delfin-toast-body .widget-html-content { white-space: normal; line-height: 1.4; }
.delfin-hdr .delfin-toast-icon { font-weight: 700; margin-right: 6px; }
.delfin-hdr .delfin-toast.is-ok .delfin-toast-icon { color: var(--dh-ok); }
.delfin-hdr .delfin-toast.is-error .delfin-toast-icon { color: var(--dh-err); }
.delfin-hdr .delfin-toast-details summary { cursor: pointer; color: var(--dh-muted); font-size: 11px; }
.delfin-hdr .delfin-toast-details pre { white-space: pre-wrap; font-size: 11px; max-height: 30vh; overflow: auto; margin: 4px 0 0 0; }
.delfin-hdr .delfin-toast-earlier { background: var(--dh-bg); border: 1px solid var(--dh-border); border-radius: 10px; padding: 4px 8px; }
.delfin-hdr .delfin-toast-earlier ul { margin: 0; padding-left: 16px; font-size: 12px; max-height: 30vh; overflow: auto; }
.delfin-hdr .delfin-toast-earlier li span { color: var(--dh-muted); }
/* the preview banner */
.delfin-hdr .delfin-preview-bar { align-items: center; gap: 10px; padding: 6px 10px; border-radius: 10px; background: #e65100; color: #fff; flex-wrap: wrap; }
.delfin-hdr .delfin-preview-bar.is-wrong { background: #b71c1c; }
.delfin-hdr .delfin-preview-bar .jupyter-button.delfin-btn { background: #fff; color: #b33c00; border-color: #fff; font-weight: 600; }
.delfin-hdr .delfin-preview-bar details { font-size: 11px; }
@media (max-width: 640px) {
  .delfin-hdr .jupyter-button.delfin-branch-btn { max-width: 60vw; }
  .delfin-hdr .delfin-brow-meta { display: none; }
}
</style>
"""


def keys_js() -> str:
    """Keyboard and click-outside behaviour of the branch menu (inline)."""
    return r"""
(function() {
  if (window.__delfinHeaderKeys) return;
  window.__delfinHeaderKeys = true;
  function pop() { return document.querySelector('.delfin-branch-pop'); }
  function open() { var p = pop(); return !!(p && p.offsetParent !== null); }
  function toggle() { var b = document.querySelector('.delfin-branch-btn'); if (b) b.click(); }
  function focusSearch() {
    var tries = 0;
    (function again() {
      var i = document.querySelector('.delfin-branch-search input');
      if (i && open()) { i.focus(); i.select(); return; }
      if (++tries < 40) setTimeout(again, 50);
    })();
  }
  document.addEventListener('keydown', function(ev) {
    if (ev.defaultPrevented) return;
    var t = ev.target, tag = (t && t.tagName || '').toUpperCase();
    var typing = tag === 'INPUT' || tag === 'TEXTAREA' || (t && t.isContentEditable);
    // Esc presses the Cancel of an open confirmation first; no key
    // presses its confirm button (rollback.Confirm).
    if (ev.key === 'Escape') {
      var c = Array.prototype.find.call(
        document.querySelectorAll('.delfin-confirm .delfin-confirm-cancel'),
        function(b) { return b.offsetParent !== null; });
      if (c) { ev.preventDefault(); ev.stopImmediatePropagation(); c.click(); return; }
    }
    if (ev.key === 'Escape' && open()) {
      ev.preventDefault(); ev.stopImmediatePropagation(); toggle(); return;
    }
    if ((ev.key.toLowerCase() === 'k' && (ev.ctrlKey || ev.metaKey)) || (ev.key === '/' && !typing)) {
      ev.preventDefault();
      ev.stopImmediatePropagation();
      if (!open()) toggle();
      focusSearch();
    }
  }, true);
  document.addEventListener('click', function(ev) {
    var t = ev.target;
    if (t && t.closest && t.closest('.delfin-branch-btn')) { setTimeout(focusSearch, 0); return; }
    if (!open()) return;
    if (t && t.closest && t.closest('.delfin-branch-wrap')) return;
    if (!document.body.contains(t)) return;
    toggle();
  }, true);
})();
"""


def css_widget():
    import ipywidgets as W
    return W.HTML(HEADER_CSS)


def classed(widget, *names):
    for n in names:
        widget.add_class(n)
    return widget


def icon_button(description: str, tooltip: str = "", *names: str,
                on_click: Optional[Callable] = None):
    import ipywidgets as W
    b = W.Button(description=description, tooltip=tooltip)
    classed(b, *(names or ("delfin-btn",)))
    if on_click is not None:
        b.on_click(lambda _b: on_click())
    return b
