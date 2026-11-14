"""The header's Dashboards menu: the kept dashboards, Open and End, and New.

Input: the kept-session records (``session.live_records``: the registry a
resume follows), the server's kernel list for windows and last activity
(``server_record.kernels``, asked with the token every kernel of the
server has), and the two addresses the launcher hands every kernel
(``session.RESUME_PATH_ENV``, ``session.DASHBOARD_PATH_ENV``).

Output: a "Dashboards" button in the header that opens a small menu:

* one row per kept dashboard: its name, open windows and last activity;
  **Open** follows the resume address in a new tab -- relative to the
  page, so it works under whatever host and port the browser used (an
  ssh forward, a tunnel) and carries no token (the browser holds the
  login cookie); **End** asks once, then drops the record and asks the
  server to shut that kernel down, which also hands its kept agent
  sessions to the agent service (resume_server.shutdown_kernel);
* the dashboard this window shows is marked and has neither (the Keep
  dashboard switch is its control);
* one kept on another login node is named, with nothing to click: its
  kernel is served there;
* **New dashboard** opens a fresh dashboard in a new tab.

Rejected: a second registry of dashboards (the kept-session records
already are one, and the server's kernel list is the truth about
windows); ending a kernel by signalling its process (the server would
restart it, and its records and hand-over would be skipped).
"""

from __future__ import annotations

import html
import os
import time
from typing import Optional
from urllib.parse import quote

from delfin.dashboard import server_record as _sr
from delfin.dashboard import session as _session


def _server(env: Optional[dict] = None) -> Optional[dict]:
    return _sr.from_kernel_environment(env)


def entries(env: Optional[dict] = None) -> list[dict]:
    """One row per kept dashboard, newest first.

    Keys: name, kernel_id, host, here, current, windows, last_activity
    (seconds, or None), started_at.
    """
    server = _server(env)
    kernel_list = _sr.kernels(server, timeout=3.0) if server else None
    by_id = {str(k.get("id") or ""): k for k in (kernel_list or [])
             if isinstance(k, dict)}
    me = _session.kernel_id()
    here_host = _session._hostname()
    rows = []
    for rec in _session.live_records():
        kid = str(rec.get("kernel_id") or "")
        host = str(rec.get("host") or "")
        kernel = by_id.get(kid)
        windows = last = None
        if kernel is not None:
            try:
                windows = int(kernel.get("connections") or 0)
            except (TypeError, ValueError):
                windows = None
            last = _sr._seconds_since(str(kernel.get("last_activity") or ""))
        rows.append({
            "name": str(rec.get("session_name") or ""),
            "kernel_id": kid,
            "host": host,
            "here": not host or host == here_host,
            "current": bool(me) and kid == me,
            "windows": windows,
            "last_activity": last,
            "started_at": rec.get("started_at"),
        })
    return rows


def open_href(name: str, env: Optional[dict] = None) -> str:
    """The resume address of *name*, relative to the page, or ""."""
    e = os.environ if env is None else env
    path = str(e.get(_session.RESUME_PATH_ENV) or "").split("?")[0]
    if not path or not name:
        return ""
    return f"{path}?{_session.SESSION_QUERY_KEY}={quote(name)}"


def new_dashboard_href(env: Optional[dict] = None) -> str:
    e = os.environ if env is None else env
    return str(e.get(_session.DASHBOARD_PATH_ENV) or "")


def end(name: str, *, env: Optional[dict] = None) -> tuple[bool, str]:
    """End the kept dashboard *name* on this host: (ended, what to say).

    The record goes first, or the server keeps the kernel (a kept one
    ignores a shutdown request); it is put back as it was when the server
    cannot be asked, so nothing is half-ended.
    """
    server = _server(env)
    if server is None:
        return False, ("This dashboard cannot reach its server to end "
                       f"{name}; `delfin-voila stop` in a terminal ends the "
                       "whole server.")
    record = next((r for r in _session.list_records()
                   if r.get("session_name") == name), None)
    if record is None:
        return True, f"{name} has ended already."
    if str(record.get("host") or "") not in ("", _session._hostname()):
        return False, (f"{name} runs on {record.get('host')}; end it from "
                       "a dashboard there.")
    kid = str(record.get("kernel_id") or "")
    path = _session.record_path(name)
    try:
        with open(path, encoding="utf-8") as handle:
            saved = handle.read()
    except OSError:
        saved = ""
    for _attempt in range(2):
        _session.drop_record(name, reason="ended from the Dashboards menu")
        if not _sr.end_kernel(server, kid):
            if saved and not os.path.exists(path):
                _restore(path, saved)
            return False, (f"The server did not end {name}; it is still "
                           "kept. `delfin-voila stop` ends the whole server.")
        if not os.path.exists(path):
            return True, (f"Ended {name}. Its agent sessions with Keep "
                          "session on go on in the agent service.")
        # The record is back: its heartbeat was writing it as it was
        # dropped, and the server kept the kernel. Once more.
        time.sleep(0.2)
    return False, f"{name} is still running; try End once more."


def _restore(path: str, text: str) -> None:
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except OSError:
        return
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(text)


def _row_markup(row: dict) -> str:
    name = html.escape(row["name"])
    bits = []
    if row["current"]:
        bits.append("this window")
    if not row["here"]:
        bits.append(f"runs on <code>{html.escape(row['host'])}</code>")
    elif row["windows"] is not None:
        n = row["windows"]
        bits.append("no window" if n == 0 else f"{n} window{'s' * (n != 1)}")
    if row["here"] and row["last_activity"] is not None:
        bits.append(f"active {_sr.ago(row['last_activity'])}")
    meta = " &middot; ".join(bits)
    return (f'<span class="delfin-dash-name">{name}</span>'
            f'<span class="delfin-dash-meta">{meta}</span>')


def _link(href: str, label: str, title: str) -> str:
    return (f'<a class="delfin-dash-link" href="{html.escape(href, quote=True)}" '
            f'target="_blank" rel="noopener" '
            f'title="{html.escape(title, quote=True)}">{html.escape(label)}</a>')


class DashboardsMenu:
    """The button and its menu (see the module docstring)."""

    def __init__(self, *, env: Optional[dict] = None, say=None):
        import ipywidgets as W

        self.W = W
        self.env = env
        self.say = say or (lambda text, level="info": None)
        self.button = W.Button(description="Dashboards",
                               tooltip="Kept dashboards: open, end, or start a new one")
        self.button.add_class("delfin-btn")
        self.button.add_class("delfin-dash-btn")
        self.button.on_click(lambda _b: self.toggle())
        self.rows = W.VBox([])
        self.rows.add_class("delfin-dash-rows")
        self.note = W.HTML("")
        self.footer = W.HTML("")
        self.pop = W.VBox([self.rows, self.note, self.footer])
        self.pop.add_class("delfin-dash-pop")
        self.pop.layout.display = "none"
        self.wrap = W.VBox([self.button, self.pop])
        self.wrap.add_class("delfin-dash-wrap")
        self._ends: dict = {}
        self._confirm = None
        self._entries: list[dict] = []

    # -- what tests and the page read --------------------------------------
    def all_widgets(self) -> list:
        out = []

        def walk(widget):
            out.append(widget)
            for child in getattr(widget, "children", ()) or ():
                walk(child)
        walk(self.pop)
        return out

    def end_buttons(self) -> dict:
        return dict(self._ends)

    def confirm_button(self):
        return self._confirm

    # -- behaviour ------------------------------------------------------------
    def toggle(self) -> None:
        opening = self.pop.layout.display == "none"
        self.pop.layout.display = "" if opening else "none"
        if opening:
            self.refresh()

    def refresh(self) -> None:
        W = self.W
        try:
            self._entries = entries(self.env)
        except Exception as exc:                     # noqa: BLE001
            self._entries = []
            self.note.value = (f'<span class="delfin-pop-note">The kept '
                               f'dashboards could not be read: {html.escape(str(exc))}</span>')
        self._ends = {}
        self._confirm = None
        children = []
        for row in self._entries:
            parts = [W.HTML(_row_markup(row))]
            parts[0].add_class("delfin-dash-label")
            if row["here"] and not row["current"]:
                href = open_href(row["name"], self.env)
                if href:
                    parts.append(W.HTML(_link(href, "Open",
                                              f"Open {row['name']} in a new tab")))
                end_btn = W.Button(description="End",
                                   tooltip=f"End {row['name']}: its tabs, forms "
                                           "and agent sessions without Keep "
                                           "session end")
                end_btn.add_class("delfin-btn")
                end_btn.add_class("is-danger")
                end_btn.on_click(lambda _b, name=row["name"]: self._ask_end(name))
                self._ends[row["name"]] = end_btn
                parts.append(end_btn)
            line = W.HBox(parts)
            line.add_class("delfin-dash-row")
            children.append(line)
        if not self._entries:
            children.append(W.HTML(
                '<span class="delfin-pop-note">No dashboard is kept. Switch '
                'on Keep dashboard to find this one here again.</span>'))
        self.rows.children = tuple(children)
        href = new_dashboard_href(self.env)
        self.footer.value = (_link(href, "New dashboard",
                                   "Start another dashboard in a new tab")
                             if href else "")
        self.note.value = ""
        self.button.description = (f"Dashboards ({len(self._entries)})"
                                   if self._entries else "Dashboards")

    def _ask_end(self, name: str) -> None:
        W = self.W
        yes = W.Button(description=f"End {name}")
        yes.add_class("delfin-btn")
        yes.add_class("is-danger")
        yes.add_class("is-filled")
        no = W.Button(description="Cancel")
        no.add_class("delfin-btn")
        no.add_class("delfin-confirm-cancel")
        text = W.HTML(
            f'<span>End <b>{html.escape(name)}</b>? Its tabs and forms '
            'end; agent sessions with Keep session on go on in the agent '
            'service.</span>')
        box = W.HBox([text, yes, no])
        box.add_class("delfin-confirm")
        self._confirm = yes

        def _do(_b):
            ok, said = end(name, env=self.env)
            self.say(said, "ok" if ok else "error")
            self.refresh()
            self.note.value = f'<span class="delfin-pop-note">{html.escape(said)}</span>'

        def _cancel(_b):
            self._confirm = None
            self.note.value = ""
            self.rows.children = tuple(c for c in self.rows.children if c is not box)

        yes.on_click(_do)
        no.on_click(_cancel)
        self.rows.children = tuple(list(self.rows.children) + [box])
