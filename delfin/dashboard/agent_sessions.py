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

# The session list, in the agent tab's own palette (tab_agent._AGENT_CSS):
# system font, slate greys, the chat's light blue for the session on screen.
# Everything is border-box and the column clips horizontally -- widgets set
# to 100% width plus their padding made it scroll sideways.
_SIDEBAR_CSS = """<style>
.delfin-session-shell { width: 100%; align-items: flex-start; overflow-x: hidden; }
.delfin-sessions {
    box-sizing: border-box; flex: 0 0 236px !important; width: 236px;
    max-width: 236px; margin: 0 12px 0 0; padding: 8px 8px 10px; gap: 6px;
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
    border: 1px solid #e5e7eb; border-radius: 8px; }
.delfin-session-form .widget-combobox, .delfin-session-form .widget-checkbox {
    width: 100% !important; margin: 0; }
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
</style>"""


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
    try:
        inside = Path(workspace).resolve().relative_to(Path(repo["root"]))
    except ValueError:
        inside = Path(".")
    return str((info.path / inside).resolve())


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
    sidebar = _classed(widgets.VBox(
        [widgets.HTML(_SIDEBAR_CSS), head, new_form, notice, list_box,
         resume_dropdown]), "delfin-sessions")
    stage = _classed(widgets.VBox(), "delfin-session-stage")
    widget = _classed(widgets.HBox([sidebar, stage]), "delfin-session-shell")

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
            if state.get("streaming"):
                rec["row"].add_class("delfin-session-busy")
            else:
                rec["row"].remove_class("delfin-session-busy")
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
            names = ", ".join(str(r.get("title") or "a session") for r in busy[:3])
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
        "form": {"new": new_btn, "workdir": workdir_box,
                 "own_worktree": own_worktree_box, "start": start_btn},
    }
