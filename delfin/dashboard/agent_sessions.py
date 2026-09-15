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
        "agent_status_html", "add_init_js",
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
    for row in rows or []:
        if isinstance(row, dict) and str(row.get("session_id") or "").strip():
            out.append({"session_id": str(row["session_id"]).strip(),
                        "workspace": str(row.get("workspace") or "")})
    return out[:_MAX_OPEN]


def save_open_sessions(rows: list[dict]) -> None:
    """Remember which sessions are open. Never raises."""
    try:
        from delfin.agent.state_paths import ensure_dir, write_text
        ensure_dir(_OPEN_SESSIONS_PATH.parent)
        write_text(_OPEN_SESSIONS_PATH,
                   json.dumps({"sessions": rows}, indent=2) + "\n")
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


def workspace_choices(ctx: Any) -> list[str]:
    """Directories offered for a new session: the default, the agent
    workspace, the calculations and the DELFIN checkout."""
    out: list[str] = []
    for d in (default_workspace(ctx), getattr(ctx, "agent_dir", None),
              getattr(ctx, "calc_dir", None), getattr(ctx, "repo_dir", None)):
        text = str(d or "").strip()
        if text and text not in out:
            out.append(text)
    return out


def _short_path(path: str) -> str:
    home = str(Path.home())
    return "~" + path[len(home):] if path.startswith(home) else path


def _session_title(state: dict) -> str:
    for message in state.get("chat_messages") or []:
        if isinstance(message, dict) and message.get("role") == "user":
            text = " ".join(str(message.get("content") or "").split())
            if text:
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

    heading = widgets.HTML("<b>Sessions</b>")
    new_btn = widgets.Button(description="+ New session",
                             tooltip="Start another session",
                             layout=widgets.Layout(width="100%"))
    workdir_box = widgets.Combobox(
        options=workspace_choices(ctx), value=default_workspace(ctx),
        placeholder="Working directory", ensure_option=False,
        layout=widgets.Layout(width="100%"))
    start_btn = widgets.Button(description="Start", button_style="primary",
                               layout=widgets.Layout(width="50%"))
    cancel_btn = widgets.Button(description="Cancel",
                                layout=widgets.Layout(width="50%"))
    new_form = widgets.VBox(
        [widgets.HTML("<span style='font-size:11px;color:#546e7a'>"
                      "Working directory</span>"),
         workdir_box, widgets.HBox([start_btn, cancel_btn])],
        layout=widgets.Layout(display="none"))
    notice = widgets.HTML("")
    list_box = widgets.VBox()
    resume_dropdown = widgets.Dropdown(
        options=[("Resume a saved session…", "")], value="",
        layout=widgets.Layout(width="100%"))
    sidebar = widgets.VBox(
        [heading, new_btn, new_form, notice, list_box, resume_dropdown],
        layout=widgets.Layout(width="240px", min_width="200px",
                              flex="0 0 240px", margin="0 12px 0 0"))
    stage = widgets.VBox(layout=widgets.Layout(flex="1 1 auto", min_width="0"))
    widget = widgets.HBox([sidebar, stage],
                          layout=widgets.Layout(width="100%",
                                                align_items="flex-start"))

    # -- lookup -------------------------------------------------------------

    def _find(key: str) -> Optional[dict]:
        return next((rec for rec in sessions if rec["key"] == key), None)

    def _state(rec: dict) -> dict:
        return rec["refs"].get("state") or {}

    def _session_id(rec: dict) -> str:
        return str(_state(rec).get("active_session_id") or "")

    def _held_elsewhere(session_id: str, key: str) -> bool:
        return any(rec["key"] != key and _session_id(rec) == session_id
                   for rec in sessions)

    def _say(text: str) -> None:
        notice.value = (f"<span style='font-size:11px;color:#b45309'>"
                        f"{html.escape(text)}</span>" if text else "")

    # -- what is on screen ----------------------------------------------------

    def _render_list() -> None:
        rows = []
        for rec in sessions:
            rec["title_btn"].button_style = (
                "info" if rec["key"] == view["active"] else "")
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
            title = str(row.get("title") or "Untitled")[:40]
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
        )
        tab, refs = build(session_ctx)
        view["scripts_added"] = True
        title_btn = widgets.Button(description="New session", tooltip=where,
                                   layout=widgets.Layout(flex="1 1 auto",
                                                         width="auto"))
        close_btn = widgets.Button(description="×", tooltip=(
            "Close this session — it stays saved and can be resumed"),
            layout=widgets.Layout(width="32px"))
        info = widgets.HTML(
            f"<div style='font-size:10px;color:#78909c;overflow-wrap:anywhere;"
            f"margin:-2px 0 4px 2px'>{html.escape(_short_path(where))}</div>")
        rec = {"key": key, "workspace": where, "tab": tab, "refs": refs or {},
               "ctx": session_ctx, "title_btn": title_btn, "info": info,
               "row": widgets.VBox([widgets.HBox([title_btn, close_btn]), info])}
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
            label = (("● " if state.get("streaming") else "")
                     + _session_title(state))
            if rec["title_btn"].description != label:
                rec["title_btn"].description = label
            sid = _session_id(rec)
            if sid != rec.get("remembered_id"):
                rec["remembered_id"] = sid
                changed = True
        active = _find(view["active"])
        if active is not None:
            try:
                ctx.agent_engine = _state(active).get("engine")
            except Exception:
                pass
        if changed:
            _persist()

    def _tick() -> None:
        try:
            refresh()
            import time as _time
            if _time.monotonic() >= view["resume_at"]:
                view["resume_at"] = _time.monotonic() + _RESUME_REFRESH_S
                _refresh_resume_options()
        except Exception:
            pass
        finally:
            timer = threading.Timer(_REFRESH_S, _tick)
            timer.daemon = True
            timer.start()
            view["timer"] = timer

    # -- controls ---------------------------------------------------------------

    def _on_new(_btn=None) -> None:
        workdir_box.options = workspace_choices(ctx)
        new_form.layout.display = ""
        new_btn.layout.display = "none"

    def _on_cancel(_btn=None) -> None:
        new_form.layout.display = "none"
        new_btn.layout.display = ""

    def _on_start(_btn=None) -> None:
        where = str(workdir_box.value or "").strip()
        if where and not Path(where).expanduser().is_dir():
            _say(f"Not a directory: {where}")
            return
        _on_cancel()
        open_session(str(Path(where).expanduser()) if where else "")

    def _on_resume(change) -> None:
        sid = str(change.get("new") or "")
        if not sid:
            return
        data = _saved_session(sid) or {}
        open_session(str(data.get("workspace") or ""), sid)
        _refresh_resume_options()

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
    _tick()

    return widget, {
        "open": open_session,
        "close": close_session,
        "activate": activate,
        "sessions": lambda: list(sessions),
        "active": lambda: _find(view["active"]),
        "refresh": refresh,
    }
