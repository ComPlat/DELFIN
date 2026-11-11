"""Swarm tab (package SW2a): build and manage swarms in the dashboard.

A VIEW and a set of buttons over delfin/agent/swarm (SW1) functions --
spec, plan, launch, status, triage, merge. No orchestration logic lives
here: every action calls one SW1 library function, and everything the
user must decide (a swarm start, a merge) goes through an explicit
click; nothing starts or merges without it.

Testable logic lives in the module-level helpers below (plain dicts in,
plain dicts/strings out) so tests run without ipywidgets, like
tests/test_tab_agent_activity.py. ``create_tab`` is the ipywidgets half.

Pause/Resume: the user asked for it, but SW1 does not ship pause/resume
verbs yet (spec/launch/status/triage/merge only). The tab fails closed:
it probes the SW1 API surface and renders pause/resume ONLY when the
library exposes them (QS ruling nacht-s11, 2026-11-02). A button that
exists but does nothing would be a lie; a missing button is the honest UI.
"""
from __future__ import annotations

import html as _html
import threading
from pathlib import Path
from typing import Any

try:  # ipywidgets is optional (headless test runs, plain Python sessions)
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:  # pragma: no cover - exercised on widget-less hosts
    widgets = None
    HAS_WIDGETS = False

from delfin.agent.swarm import launch as _launch
from delfin.agent.swarm import merge as _merge
from delfin.agent.swarm import status as _status
from delfin.agent.swarm import triage as _triage
from delfin.agent.swarm.spec import SpecError, load_spec


def available_actions(_report: dict[str, Any] | None = None) -> tuple[str, ...]:
    """The actions this tab offers, derived from SW1's real API surface.

    Fail-closed capability probe (QS ruling): a verb is listed only when
    the SW1 library actually exposes it. 'start' backs the guided flow's
    Start button (launch.start_members); 'merge_plan' backs the Pre-PR
    check display (merge.merge_plan, which prints commands and never
    executes); 'pause'/'resume' appear only if SW1 grows the verbs.
    Probed here: delfin.agent.swarm (the package root), so a future
    top-level pause_all/resume_all is picked up; the probe never imports
    a name that does not exist (getattr with a default, not hasattr on a
    submodule import that would raise).
    """
    verbs = ["start", "merge_plan"]
    import delfin.agent.swarm as _swarm_pkg

    for verb in ("pause_all", "resume_all"):
        if getattr(_swarm_pkg, verb, None) is not None:
            verbs.append(verb)
    return tuple(verbs)


def triage_decide(item: dict[str, Any], decision: str, *,
                  reason: str = "") -> dict[str, Any]:
    """Approve/Deny one triage item -- through SW1's own machinery ONLY.

    Approve: SW1's triage_answer (triage.py:108), which re-derives the
    verdict itself and refuses escalations with TriageError -- the tab
    passes NO verdict of its own (T4/D13: a caller's verdict is untrusted
    data). Deny: SW1 ships no deny verb, so the item is left for the user
    (fail-closed; the tab never answers on the user's behalf without SW1).
    The dialog text and the auto-answered log come from SW1's triage_list
    / the log entry triage_answer writes -- no tab-side log exists.
    """
    if decision == "approve":
        try:
            out = _triage.triage_answer(item, None, answer=reason)
        except _triage.TriageError as exc:
            # SW1's refusal is the answer to show, never swallowed.
            return {"answered": False, "item_id": item.get("id"),
                    "error": str(exc)}
        return dict(out, answered=True)
    # deny (or anything else): fail closed.
    return {"answered": False, "item_id": item.get("id"),
            "left_for_user": True}


# ---------------------------------------------------------------------------
# Board: a status report -> card rows (pure data)
# ---------------------------------------------------------------------------


def board_cards(report: dict[str, Any]) -> list[dict[str, Any]]:
    """Turn a swarm status report into one card row per package.

    The view does not rename or drop anything from SW1's schema: phase
    and member health/waiting pass through verbatim, so the dashboard
    can never show a state the status contract does not define. A
    document that is not a status report raises -- an empty board would
    hide the mismatch from the user.
    """
    if not isinstance(report, dict) or "packages" not in report:
        raise ValueError(
            "not a swarm status report: no 'packages' key "
            "(collect_status output expected)")
    cards = []
    for pkg in report.get("packages", []):
        members = [
            {"name": m["name"], "role": m.get("role", ""),
             "health": m.get("health", ""),
             "waits_on": m.get("waits_on")}
            for m in pkg.get("members", [])
        ]
        cards.append({
            "name": pkg["name"],
            "phase": pkg.get("phase", ""),
            "members": members,
        })
    return cards


def waiting_panel(report: dict[str, Any]) -> list[dict[str, Any]]:
    """The 'waits on you' rows: members waiting on the user."""
    if not isinstance(report, dict):
        raise ValueError("not a swarm status report (dict expected)")
    return list(report.get("waiting_for_user") or [])


# ---------------------------------------------------------------------------
# Guided flow: spec text -> validation -> plan preview (pure data)
# ---------------------------------------------------------------------------


def validate_spec_text(spec_text: str) -> dict[str, Any]:
    """Validate swarm.yaml text through SW1, for the guided flow.

    Returns {'ok': bool, 'error': prose, 'packages': [...]} -- a
    validation failure is a message the flow shows inline, never a
    traceback. The spec is loaded to a temp file so SW1's load_spec (the
    function the CLI calls) stays the single validation path.
    """
    import tempfile
    from pathlib import Path

    result: dict[str, Any] = {"ok": False, "error": "", "packages": []}
    if not isinstance(spec_text, str) or not spec_text.strip():
        result["error"] = "The swarm.yaml text is empty."
        return result
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "swarm.yaml"
        path.write_text(spec_text, encoding="utf-8")
        try:
            spec = load_spec(path)
        except SpecError as exc:
            result["error"] = str(exc)
            return result
        # The CLI only validates through plan_launch (launch.py:54 calls
        # validate_spec inside plan_members), so the guided flow uses the
        # same call: rules enforced identically, one path, no drift.
        try:
            _launch.plan_launch(spec)
        except SpecError as exc:
            result["error"] = str(exc)
            return result
        # A load error that is not a SpecError (YAML crash etc.) is not
        # swallowed: let it raise so the user sees the real failure.
    result["ok"] = True
    result["packages"] = [
        {"name": p.get("name", ""),
         "owned": list(p.get("owned") or []),
         "roles": dict(p.get("roles") or {"builder": True})}
        for p in spec.get("packages", [])
    ]
    return result


def plan_preview(spec_text: str) -> dict[str, Any]:
    """The dry-run preview of exactly what would start (guided flow).

    Calls SW1's plan_launch + render_dry_run -- the same functions and
    the same text the CLI shows -- and starts NOTHING: start_members
    is only reachable from the explicit Start button handler.
    """
    ok = validate_spec_text(spec_text)
    if not ok["ok"]:
        _clear_previewed_spec()
        return {"ok": False, "error": ok["error"], "preview": ""}
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "swarm.yaml"
        path.write_text(spec_text, encoding="utf-8")
        spec = load_spec(path)  # validated above
        plan = _launch.plan_launch(spec)
    # F1 (reviewer adversarial round 2): record exactly what was previewed,
    # so the Start handler launches THAT text and not a re-read of the
    # textarea the user may have edited since the preview (TOCTOU).
    global last_previewed_spec_text
    last_previewed_spec_text = spec_text
    return {"ok": True, "error": "", "preview": _launch.render_dry_run(plan),
            "goal": plan.get("goal", "")}


# ---------------------------------------------------------------------------
# Triage queue -> panel rows (pure data)
# ---------------------------------------------------------------------------


def triage_rows(queue: list[dict[str, Any]], *,
                roster: set[str] | frozenset[str] = frozenset()) \
        -> list[dict[str, Any]]:
    """The triage queue as panel rows, sorted by SW1's triage_list.

    Escalations come first (the user's queue), auto-answerables last.
    Each row carries the item's real dialog text -- the user decides on
    what was actually said, not a summary -- and SW1's verdict.
    """
    if not isinstance(queue, list):
        raise ValueError("queue must be a list of item dicts")
    roster_fs = frozenset(roster or ())
    prepared = [dict(item, roster=roster_fs) for item in queue]
    sorted_items = _triage.triage_list(prepared)
    rows = []
    for item in sorted_items:
        rows.append({
            "id": item.get("id", ""),
            "member": item.get("member", ""),
            "kind": item.get("kind", ""),
            "from": item.get("from", ""),
            "text": item.get("text", ""),
            "action": item.get("action", ""),
            "reason": item.get("reason", ""),
        })
    return rows


# ---------------------------------------------------------------------------
# SW1 identity pins (reviewer contract): the tab re-exports SW1's status
# functions instead of re-implementing any of them.
# ---------------------------------------------------------------------------

_validate_status = _status.validate_report
_collect_status = _status.collect_status

STATUS_SCHEMA_VERSION = _status.STATUS_SCHEMA_VERSION

# F1: the one record of what a successful preview showed. Set by
# plan_preview on success, cleared by _clear_previewed_spec on failure;
# _on_start launches exactly this text or refuses.
last_previewed_spec_text: str | None = None


def _clear_previewed_spec() -> None:
    """Forget the last successfully previewed spec (F1 part 2).

    A failed preview must not leave a stale record Start could launch:
    the user may have edited the textarea between preview and click.
    """
    global last_previewed_spec_text
    last_previewed_spec_text = None


# ---------------------------------------------------------------------------
# Save / load of swarm.yaml + Pre-PR display (phase 2, all pure data)
# ---------------------------------------------------------------------------


def save_spec(spec_text: str, target: Path | str, *,
              force: bool = False) -> Path:
    """Persist the guided-flow spec through SW1's own writer/loader.

    load_spec is the consumer: the file must round-trip through the same
    loader the CLI and start_members use, so nothing is written by hand.
    A non-empty spec is validated (validate_spec, the CLI's own rules)
    before anything lands on disk; an empty one falls back to SW1's
    commented example via write_example_spec, whose overwrite protection
    (SpecError via the loader path / LaunchError in launch.py) stays
    unchanged.
    """
    target_path = Path(target)
    if target_path.exists() and not force:
        # Same rule as SW1's write_example_spec: refuse to overwrite.
        raise SpecError(f"{target_path} already exists -- not overwriting")
    if spec_text and spec_text.strip():
        # The CLI's own validator (validate_spec_text -> validate_spec);
        # anything invalid never reaches the disk.
        result = validate_spec_text(spec_text)
        if not result["ok"]:
            raise SpecError(result["error"])
        target_path.write_text(spec_text, encoding="utf-8")
        return target_path
    # Empty/whitespace spec: write the commented example instead.
    return _launch.write_example_spec(target_path, force=force)


def load_spec_text(target: Path | str) -> str:
    """Read a swarm.yaml back as text for the guided flow textarea."""
    return Path(target).read_text(encoding="utf-8")


def pre_pr_display(package: str, *, base: str = "main") -> dict[str, Any]:
    """The Pre-PR check as the COMMAND for the user (never executed).

    Delegates to SW1's merge_plan: commands as text, ``executed`` always
    False, check state ``not_run`` until the user actually ran it.
    """
    return _merge.merge_plan(package, base=base)


# ---------------------------------------------------------------------------
# ipywidgets half
# ---------------------------------------------------------------------------


def create_tab(ctx: Any):
    """Build the Swarm tab widget. Read-only until the user clicks.

    Refresh reuses the dashboard's throttled timer pattern
    (tab_agent_activity.py:787-802): one collect_status call per interval
    through _collect_status (SW1's own function), cooperative stop, no
    second scheduler loop of its own.
    """
    if not HAS_WIDGETS:
        raise ImportError("ipywidgets is required for the Swarm tab.")

    state: dict[str, Any] = {"_timer": None, "_stopped": False}

    title = widgets.HTML("<h3>Swarm</h3>")

    # -- board -------------------------------------------------------------
    board_box = widgets.VBox([])
    waits_box = widgets.VBox([])
    status_html = widgets.HTML("")

    def _refresh() -> None:
        spec_path = getattr(ctx, "swarm_spec_path", None)
        try:
            if not spec_path or not Path(spec_path).exists():
                status_html.value = (
                    "<em>No swarm.yaml yet. Use 'New swarm' below.</em>")
                board_box.children = ()
                waits_box.children = ()
                return
            spec = load_spec(Path(spec_path))
            report = _collect_status(spec)
            cards = board_cards(report)
            rows = []
            for card in cards:
                healths = ", ".join(
                    f"{m['name']} ({m['health']})" for m in card["members"])
                rows.append(widgets.HTML(
                    f"<b>{_html.escape(card['name'])}</b> — "
                    f"{_html.escape(card['phase'])}<br>"
                    f"<small>{_html.escape(healths)}</small>"))
            board_box.children = tuple(rows) or (
                widgets.HTML("<em>No packages.</em>"),)
            wait_rows = [widgets.HTML(
                f"<b>{_html.escape(w['member'])}</b> waits on: "
                f"{_html.escape(str(w.get('waits_on', '')))}")
                for w in waiting_panel(report)]
            waits_box.children = tuple(wait_rows) or (
                widgets.HTML("<em>Nothing waits on you.</em>"),)
            status_html.value = ""
        except Exception as exc:  # a broken spec/report must not kill the tab
            status_html.value = (
                f"<span style='color:#b91c1c'>{_html.escape(str(exc))}</span>")

    def _tick() -> None:
        if state["_stopped"]:
            return
        _refresh()
        timer = threading.Timer(5.0, _tick)
        timer.daemon = True
        state["_timer"] = timer
        timer.start()

    refresh_btn = widgets.Button(description="Refresh", icon="refresh",
                                 button_style="info")
    refresh_btn.on_click(lambda _b: _refresh())

    board = widgets.VBox([widgets.HTML("<b>Packages</b>"), board_box])
    waits = widgets.VBox([widgets.HTML("<b>Waits on you</b>"), waits_box])
    top = widgets.HBox([board, waits])

    # -- guided flow: goal -> editable spec -> preview -> Start -------------
    spec_area = widgets.Textarea(
        value="", placeholder="Edit swarm.yaml here…",
        layout=widgets.Layout(width="100%", height="240px"))
    validate_btn = widgets.Button(description="Validate", icon="check")
    preview_btn = widgets.Button(description="Dry-run preview",
                                 icon="eye", button_style="info")
    start_btn = widgets.Button(
        description="Start", icon="play", button_style="danger")
    start_btn.disabled = True  # nothing starts without preview + click
    flow_status = widgets.HTML("")
    preview_out = widgets.HTML("")

    spec_state: dict[str, Any] = {"preview_ok": False}

    def _on_validate(_b=None) -> None:
        result = validate_spec_text(spec_area.value or "")
        if result["ok"]:
            flow_status.value = "<span style='color:#15803d'>Valid.</span>"
        else:
            flow_status.value = (
                f"<span style='color:#b91c1c'>"
                f"{_html.escape(result['error'])}</span>")
            spec_state["preview_ok"] = False
            start_btn.disabled = True
        _on_preview()

    def _on_preview(_b=None) -> None:
        result = plan_preview(spec_area.value or "")
        if result["ok"]:
            preview_out.value = "<pre>" + _html.escape(result["preview"]) + \
                "</pre>"
            spec_state["preview_ok"] = True
            start_btn.disabled = False
        else:
            preview_out.value = ""
            spec_state["preview_ok"] = False
            start_btn.disabled = True
            flow_status.value = (
                f"<span style='color:#b91c1c'>"
                f"{_html.escape(result['error'])}</span>")

    def _on_start(_b) -> None:
        """The ONLY path that starts sessions; launches the PREVIEWED text.

        F1: never re-reads the textarea -- the spec the user saw in the
        preview is the spec that starts. Without a fresh preview record
        nothing starts.
        """
        if not spec_state["preview_ok"] or not last_previewed_spec_text:
            return  # nothing starts without the user having seen the plan
        import tempfile
        from pathlib import Path as _P
        with tempfile.TemporaryDirectory() as tmp:
            path = _P(tmp) / "swarm.yaml"
            path.write_text(last_previewed_spec_text, encoding="utf-8")
            spec = load_spec(path)
        try:
            result = _launch.start_members(spec, execute=True)
        except Exception as exc:
            flow_status.value = (
                f"<span style='color:#b91c1c'>"
                f"{_html.escape(str(exc))}</span>")
            return
        started = ", ".join(m["name"] for m in result.get("started", []))
        flow_status.value = (
            f"<span style='color:#15803d'>Started: "
            f"{_html.escape(started)}</span>")

    validate_btn.on_click(_on_validate)
    preview_btn.on_click(_on_preview)
    start_btn.on_click(_on_start)

    actions = available_actions()
    pause_btn = widgets.Button(description="Pause all", icon="pause") \
        if "pause" in actions else None
    resume_btn = widgets.Button(description="Resume", icon="forward") \
        if "resume" in actions else None

    buttons = [validate_btn, preview_btn, start_btn]
    if pause_btn:
        buttons.append(pause_btn)
    if resume_btn:
        buttons.append(resume_btn)

    flow = widgets.VBox([
        widgets.HTML("<b>New swarm</b> — goal → plan → preview → Start"),
        spec_area, widgets.HBox(buttons), flow_status, preview_out,
    ])

    tab_widget = widgets.VBox([title, top, flow, status_html])
    _tick()
    return tab_widget, {"board_box": board_box, "waits_box": waits_box,
                        "flow_status": flow_status, "state": state}


# ---------------------------------------------------------------------------
# SWV watch pane (read-only): observe a swarm and one member session.
# Pure, HTML-able formatters over delfin.agent.swarm.observe — the SAME data
# API the TUI and an SW3 lead call. The widgets half is create_watch_tab,
# registered as its own tab so the SW2a manage flow above stays untouched.
# ---------------------------------------------------------------------------

from delfin.agent.swarm.observe import member_feed as watch_feed
from delfin.agent.swarm.observe import swarm_overview as watch_overview


def watched_rows_html(rows: list[dict[str, Any]]) -> list[str]:
    """One HTML string per member of an ``swarm_overview`` result.

    Pure and write-free: every field is HTML-escaped so hostile input stays
    inert, a member awaiting approval is visually flagged, and the rendered
    HTML carries no control (no approve button, no input) — observing a
    member never acts on it.
    """
    out: list[str] = []
    for row in rows or []:
        name = _html.escape(str(row.get("name") or ""))
        package = _html.escape(str(row.get("package") or ""))
        role = _html.escape(str(row.get("role") or ""))
        model = _html.escape(str(row.get("model") or ""))
        state = _html.escape(str(row.get("state") or ""))
        last = row.get("last_tool_call") or {}
        last_tool = _html.escape(
            f"{last.get('tool') or ''} {last.get('arg') or ''}".strip())
        n_commits = int(row.get("commits_on_branch") or 0)
        n_approvals = int(row.get("pending_approvals_count") or 0)
        row_cls = (" color:#b91c1c;" if state == "NEEDS APPROVAL"
                   else "")
        cells = [
            f"<b>{name}</b>", package, role, model,
            f"<span style='{row_cls}'>{state}</span>",
            _html.escape(last_tool)[:120],
            f"{n_commits} commit(s)",
            (f"{n_approvals} approval(s)" if n_approvals else "no approvals"),
        ]
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    return out


def _ago(seconds: Any) -> str:
    try:
        s = float(seconds)
    except (TypeError, ValueError):
        return "-"
    if s < 90:
        return f"{int(s)}s"
    if s < 5400:
        return f"{int(s // 60)}m"
    return f"{s / 3600:.1f}h"


#: A session whose presence heartbeat (about once a minute) is older than
#: this is not listed: its kernel or its job has ended.
HEARTBEAT_STALE_S = 300.0

#: The all-sessions view is read at most this often (the tab ticks every
#: 5 s; the read is a directory of small files, kept rare on a login node).
RUNNING_REFRESH_S = 30.0


def running_session_rows(records: list[dict] | None = None, *,
                         now: float | None = None,
                         settings: dict | None = None,
                         beats: dict | None = None) -> list[dict]:
    """One row per open agent session of this user, on ANY host.

    The Swarm Watch tab without a roster: a wave or swarm started by a
    script (a SLURM job on a compute node) has no roster the dashboard
    knows, and the tab showed nothing. The rows come from the presence
    records in the shared HOME (``session_health.collect``), so sessions
    on other nodes are listed with their host. Read-only.
    """
    import time as _time

    from delfin.agent import session_health
    if records is None:
        try:
            records = session_health.collect(with_git=False)
        except Exception:
            records = []
    now = _time.time() if now is None else now
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = (load_settings() or {}).get("agent") or {}
        except Exception:
            settings = {}
    if beats is None:
        try:
            from delfin.agent import session_presence
            beats = {str(r.get("key") or ""): float(r.get("updated_at") or 0)
                     for r in session_presence.open_sessions()}
        except Exception:
            beats = {}
    rows = []
    for rec in records or []:
        if not isinstance(rec, dict):
            continue
        if not rec.get("host") and not rec.get("session_id"):
            continue          # a pause flag with no session behind it
        beat = beats.get(str(rec.get("name") or ""))
        if beats and (beat is None or now - beat > HEARTBEAT_STALE_S):
            continue          # its node stopped answering (a job that ended)
        try:
            verdict = session_health.classify(rec, now, settings).get(
                "verdict") or ""
        except Exception:
            verdict = ""
        last = rec.get("last_turn") if isinstance(rec.get("last_turn"), dict) \
            else {}
        rows.append({
            "name": str(rec.get("name") or rec.get("session_id") or ""),
            "host": str(rec.get("host") or "").split(".")[0],
            "branch": str(rec.get("branch") or ""),
            "model": str(last.get("model") or ""),
            "state": verdict,
            "idle": _ago(rec.get("idle_s")),
            "approvals": len(rec.get("open_approvals") or []),
            "inbox": int(rec.get("inbox") or 0),
        })
    rows.sort(key=lambda r: (r["host"], r["name"]))
    return rows


RUNNING_HEADER = ("<tr><td><b>session</b></td><td>host</td><td>branch</td>"
                  "<td>model</td><td>state</td><td>last turn</td>"
                  "<td>approvals</td><td>inbox</td></tr>")


def running_rows_html(rows: list[dict[str, Any]]) -> list[str]:
    """HTML rows for ``running_session_rows``; escaped, no controls."""
    out: list[str] = []
    for row in rows or []:
        n_appr = int(row.get("approvals") or 0)
        state = _html.escape(str(row.get("state") or ""))
        flag = " color:#b91c1c;" if n_appr or "waiting" in state else ""
        cells = [
            f"<b>{_html.escape(str(row.get('name') or ''))}</b>",
            _html.escape(str(row.get("host") or "")),
            _html.escape(str(row.get("branch") or ""))[:60],
            _html.escape(str(row.get("model") or "")),
            f"<span style='{flag}'>{state}</span>",
            _html.escape(str(row.get("idle") or "")) + " ago",
            (f"<span style='color:#b91c1c'>{n_appr} open</span>" if n_appr
             else "-"),
            str(int(row.get("inbox") or 0)),
        ]
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    return out


def session_view_html(name: str, *, records: list[dict] | None = None,
                      load: Any = None, feed: Any = None,
                      last_n: int = 12) -> str:
    """What one running session said and did lately, as escaped HTML.

    Read-only, from the shared HOME, so a session on another node is
    shown the same way: its saved conversation (``session_store``, written
    at each turn's end) and its latest tool steps (``member_feed``, the
    audit log). Nothing here acts on the session.
    """
    from delfin.agent import session_store
    load = load or session_store.load_session
    feed = feed or watch_feed
    if records is None:
        try:
            from delfin.agent import session_presence
            records = session_presence.open_sessions()
        except Exception:
            records = []
    sid = ""
    for rec in records or []:
        if str(rec.get("key") or rec.get("name") or "") == name:
            sid = str(rec.get("session_id") or "")
            break
    parts = [f"<h4>{_html.escape(name)}</h4>"]
    data = load(sid) if sid else None
    msgs = (data or {}).get("chat_messages") or []
    if msgs:
        parts.append("<div><b>Conversation</b> (as of its last turn end)</div>")
        for m in msgs[-last_n:]:
            if not isinstance(m, dict):
                continue
            role = str(m.get("role") or "")
            text = str(m.get("content") or "")
            who = ("message" if role == "user" and m.get("_harness")
                   else role or "?")
            colour = {"assistant": "#1d4ed8", "user": "#065f46"}.get(role, "#6b7280")
            parts.append(
                f"<div style='margin:4px 0; padding:4px 6px; border-left:3px solid {colour}; "
                f"white-space:pre-wrap'><b>{_html.escape(who)}</b>: "
                f"{_html.escape(text[:1500])}{' ...' if len(text) > 1500 else ''}</div>")
    else:
        parts.append("<div><em>No saved conversation yet.</em></div>")
    try:
        steps = feed(name, since=0) or []
    except Exception:
        steps = []
    if steps:
        parts.append("<div style='margin-top:6px'><b>Latest steps</b></div>")
        parts.extend(f"<div>{h}</div>" for h in watched_feed_html(steps[-15:]))
    return "".join(parts)


def send_to_session(name: str, text: str, *, send: Any = None) -> str:
    """Leave *text* in session *name*'s inbox (from the operator, never as
    the user -- the session reads it as a session message). Returns what to
    say."""
    text = str(text or "").strip()
    if not name or not text:
        return "Pick a session and write a message first."
    if send is None:
        from delfin.agent import ops_interface
        send = ops_interface.send
    try:
        out = send(name, text) or {}
    except Exception as exc:
        return f"Not sent: {exc}"
    return (f"Sent to {name} -- it reads it at its next step "
            f"(id {str(out.get('id') or '')[:8]}).")


def watched_feed_html(feed: list[dict[str, Any]]) -> list[str]:
    """One HTML string per ``member_feed`` entry, newest last.

    Refusals are visually flagged; the feed renders what happened, it never
    offers to act on it (no approve/answer/input control).
    """
    out: list[str] = []
    for entry in feed or []:
        kind = str(entry.get("kind") or "tool")
        text = _html.escape(str(entry.get("text") or ""))
        at = _html.escape(str(entry.get("at") or ""))
        if kind == "refusal":
            line = (f"<span style='color:#b91c1c'>[refusal]</span> "
                    f"{at} — {text}")
        elif kind == "commit":
            line = f"[commit] {text}"
        else:
            line = f"[tool] {at} — {text}"
        out.append(line)
    return out


#: Rules of the agent tab's own stylesheets that give Swarm Watch the same
#: look (chat bubbles, the composer and its round Send, the session list).
_LOOK_KEYS = ("delfin-chat-msg", "delfin-chat-user", "delfin-chat-agent",
              "delfin-chat-system", "delfin-agent-composer", "delfin-agent-send",
              "delfin-agent-input", "delfin-sessions", "delfin-session-item",
              "delfin-session-title", "delfin-session-path",
              "delfin-session-heading", "delfin-session-head",
              "delfin-session-list", "delfin-session-text",
              "delfin-session-active")


def _strip_blocks(css: str) -> str:
    """CSS without comments and without at-rule blocks (@media, @keyframes,
    @supports), which a flat rule parser would cut in half."""
    import re as _re
    css = _re.sub(r"/\*.*?\*/", "", css, flags=_re.S)
    out, i = [], 0
    while i < len(css):
        at = css.find("@", i)
        if at < 0:
            out.append(css[i:])
            break
        out.append(css[i:at])
        brace = css.find("{", at)
        semi = css.find(";", at)
        if brace < 0 or (0 <= semi < brace):
            i = (semi + 1) if semi >= 0 else len(css)
            continue
        depth, j = 0, brace
        while j < len(css):
            if css[j] == "{":
                depth += 1
            elif css[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        i = j + 1
    return "".join(out)


def agent_look_css() -> str:
    """The agent tab's look for Swarm Watch, taken from its own stylesheets
    (tab_agent._AGENT_CSS, agent_sessions._SIDEBAR_CSS) at import time, with
    every class renamed into Swarm Watch's own namespace (``delfin-swv-``).

    One source: a change to the agent tab's design reaches this tab too, so
    no second, diverging Send button or bubble can grow here. Renamed
    because the agent tab's page scripts find its controls by class across
    the whole page (its Enter-to-send clicks the first
    ``.delfin-agent-send-row button``) and must never find these.
    """
    import re as _re
    try:
        from delfin.dashboard import agent_sessions as _as
        from delfin.dashboard import tab_agent as _ta
        sources = [_ta._AGENT_CSS, _as._SIDEBAR_CSS]
    except Exception:
        return ""
    rules = []
    for src in sources:
        body = src.replace("<style>", "").replace("</style>", "")
        for sel, decl in _re.findall(r"([^{}]+)\{([^{}]*)\}", _strip_blocks(body)):
            sel = " ".join(sel.split())
            if not any(k in sel for k in _LOOK_KEYS):
                continue
            sel = sel.replace(".delfin-agent-root ", ".delfin-swv ")
            sel = _re.sub(r"\.delfin-(agent|chat|session)", r".delfin-swv-\1", sel)
            rules.append(f"{sel} {{{decl.strip()}}}")
    return "<style>\n" + "\n".join(rules) + "\n" + _SWV_OWN_CSS + "\n</style>"


#: What Swarm Watch adds around the borrowed look: the two columns, the
#: scrolling chat, the header line.
_SWV_OWN_CSS = """
.delfin-swv-shell { width: 100%; align-items: flex-start; gap: 12px; }
.delfin-swv-main { flex: 1 1 auto; min-width: 0; gap: 8px; }
.delfin-swv-head { font-size: 13px; color: #334155; padding: 2px 4px; }
.delfin-swv-head b { color: #0f172a; font-size: 14px; }
.delfin-swv-chat { max-height: calc(100vh - 380px); min-height: 240px;
    overflow-y: auto !important; padding: 8px; border: 1px solid #e5e7eb;
    border-radius: 10px; background: #fafafa; }
.delfin-swv-chat .widget-html-content { display: flex; flex-direction: column; }
.delfin-swv-note { font-size: 11px; color: #64748b; }
.delfin-swv-steps { font-size: 11px; color: #475569; }
.delfin-swv-session-title.delfin-swv-needs::before { content: "";
    display: inline-block; width: 7px; height: 7px; margin: 0 6px 1px 0;
    border-radius: 50%; background: #f97316; vertical-align: middle; }
"""

#: Enter sends from the Swarm Watch box (Shift+Enter: a new line).
_SWV_KEYS_JS = """
(function () {
    if (window.__delfinSwvKeys) return;
    window.__delfinSwvKeys = true;
    document.addEventListener('keydown', function (e) {
        var t = e.target;
        if (!t || t.tagName !== 'TEXTAREA' || !t.closest) return;
        if (e.key !== 'Enter' || e.shiftKey || e.ctrlKey || e.metaKey) return;
        var box = t.closest('.delfin-swv-agent-composer');
        var send = box && box.querySelector('.delfin-swv-agent-send');
        if (send) { e.preventDefault(); setTimeout(function () { send.click(); }, 80); }
    }, true);
})();
"""


def chat_bubbles_html(name: str, *, records: list[dict] | None = None,
                      load: Any = None, last_n: int = 40) -> str:
    """One session's saved conversation as the agent tab's bubbles
    (renamed classes): the person's messages right, the agent's left with
    its markdown, messages from sessions and the harness as system lines.
    Escaped throughout (chat_markdown emits only its own markup)."""
    from delfin.agent import session_store
    from delfin.dashboard import chat_markdown
    load = load or session_store.load_session
    if records is None:
        try:
            from delfin.agent import session_presence
            records = session_presence.open_sessions()
        except Exception:
            records = []
    sid = next((str(r.get("session_id") or "") for r in records or []
                if str(r.get("key") or r.get("name") or "") == name), "")
    data = load(sid) if sid else None
    msgs = [m for m in (data or {}).get("chat_messages") or [] if isinstance(m, dict)]
    if not msgs:
        return "<i>No saved conversation yet -- it is written at each turn's end.</i>"

    def code(lang: str, src: str) -> str:
        return "<pre><code>" + _html.escape(src) + "</code></pre>"

    out = []
    for m in msgs[-last_n:]:
        role = str(m.get("role") or "")
        text = str(m.get("content") or "")
        if role == "assistant":
            body = chat_markdown.render(text[:20000], code) or _html.escape(text[:2000])
            out.append(f'<div class="delfin-swv-chat-msg delfin-swv-chat-agent">{body}</div>')
        elif role == "user" and not m.get("_harness"):
            out.append('<div class="delfin-swv-chat-msg delfin-swv-chat-user">'
                       + _html.escape(text[:4000]) + "</div>")
        else:
            out.append('<div class="delfin-swv-chat-msg delfin-swv-chat-system">'
                       + _html.escape(text[:600]) + ("…" if len(text) > 600 else "")
                       + "</div>")
    return "".join(out)


def create_watch_tab(ctx: Any):
    """Swarm Watch, laid out like the agent tab: the running sessions of
    every node on the left, the picked one's conversation as the agent
    tab's bubbles on the right, and the same composer below it (it sends a
    message from the operator to the session's inbox, never as the user).
    Read-only otherwise: no approval or stop control here (RX1/SV1).

    With a roster or spec on the context, a table of its members follows
    the conversation (the observe API, as before).
    """
    if not HAS_WIDGETS:
        raise ImportError("ipywidgets is required for the Swarm Watch tab.")
    import time as _time

    state: dict[str, Any] = {"_timer": None, "_stopped": False,
                             "picked": "", "rows": [], "_list_at": 0.0,
                             "_chat_at": 0.0}

    try:
        add_js = getattr(ctx, "add_init_js", None)
        if callable(add_js):
            add_js(_SWV_KEYS_JS)
    except Exception:
        pass

    heading = widgets.HTML("Running sessions")
    heading.add_class("delfin-swv-session-heading")
    refresh_btn = widgets.Button(description="↻", tooltip="Read the sessions now",
                                 layout=widgets.Layout(width="28px", height="26px"))
    head = widgets.HBox([heading, refresh_btn])
    head.add_class("delfin-swv-session-head")
    list_box = widgets.VBox([])
    list_box.add_class("delfin-swv-session-list")
    count_html = widgets.HTML("")
    count_html.add_class("delfin-swv-note")
    sidebar = widgets.VBox([head, list_box, count_html])
    sidebar.add_class("delfin-swv-sessions")

    head_html = widgets.HTML("<i>Pick a session on the left.</i>")
    head_html.add_class("delfin-swv-head")
    chat_html = widgets.HTML("")
    chat_html.add_class("delfin-swv-chat")
    steps_html = widgets.HTML("")
    steps_html.add_class("delfin-swv-steps")
    msg_box = widgets.Textarea(
        placeholder="Message to this session (Enter sends; it reads it as a "
                    "message from the operator, not as the user)",
        layout=widgets.Layout(width="100%"))
    msg_box.add_class("delfin-swv-agent-input")
    send_btn = widgets.Button(description="Send", tooltip="Send (Enter)",
                              layout=widgets.Layout(width="34px", height="34px"))
    send_btn.add_class("delfin-swv-agent-send")
    send_row = widgets.HBox([send_btn])
    send_row.add_class("delfin-swv-agent-send-row")
    send_note = widgets.HTML("")
    send_note.add_class("delfin-swv-note")
    gap = widgets.HTML("")
    gap.add_class("delfin-swv-agent-composer-gap")
    bar = widgets.HBox([send_note, gap, send_row])
    bar.add_class("delfin-swv-agent-composer-bar")
    composer = widgets.VBox([msg_box, bar])
    composer.add_class("delfin-swv-agent-composer")
    table_html = widgets.HTML("")
    main = widgets.VBox([head_html, chat_html, steps_html, composer, table_html])
    main.add_class("delfin-swv-main")
    shell = widgets.HBox([sidebar, main])
    shell.add_class("delfin-swv-shell")
    root = widgets.VBox([widgets.HTML(agent_look_css()), shell])
    root.add_class("delfin-swv")

    def _pick(name: str) -> None:
        state["picked"] = name
        state["_chat_at"] = 0.0
        _render_list()
        _render_chat()

    def _render_list() -> None:
        items = []
        for row in state["rows"]:
            name = row["name"]
            title = widgets.Button(description=name,
                                   tooltip=f"{name} on {row['host']}")
            title.add_class("delfin-swv-session-title")
            if row.get("approvals"):
                title.add_class("delfin-swv-needs")
            title.on_click(lambda _b, _n=name: _pick(_n))
            path = widgets.HTML(_html.escape(
                f"{row['host']} · {row['state'] or '?'} · {row['branch']}"))
            path.add_class("delfin-swv-session-path")
            text = widgets.VBox([title, path])
            text.add_class("delfin-swv-session-text")
            item = widgets.HBox([text])
            item.add_class("delfin-swv-session-item")
            if name == state["picked"]:
                item.add_class("delfin-swv-session-active")
            items.append(item)
        list_box.children = tuple(items)
        count_html.value = (f"{len(items)} running on "
                            f"{len({r['host'] for r in state['rows']})} node(s)"
                            if items else "No session runs.")

    def _render_chat() -> None:
        name = state["picked"]
        if not name:
            return
        row = next((r for r in state["rows"] if r["name"] == name), None)
        if row:
            appr = (f" · <span style='color:#c2410c'>{row['approvals']} approval(s) open</span>"
                    if row["approvals"] else "")
            head_html.value = (f"<b>{_html.escape(name)}</b> · {_html.escape(row['host'])}"
                               f" · {_html.escape(row['branch'])} · {_html.escape(row['model'])}"
                               f" · {_html.escape(row['state'])}{appr}")
        chat_html.value = chat_bubbles_html(name)
        try:
            steps = watch_feed(name, since=0) or []
        except Exception:
            steps = []
        steps_html.value = ("<b>Latest steps</b><br>" + "<br>".join(
            watched_feed_html(steps[-8:]))) if steps else ""
        state["_chat_at"] = _time.time()

    def _send(_b=None) -> None:
        said = send_to_session(state["picked"], msg_box.value)
        send_note.value = _html.escape(said)
        if said.startswith("Sent"):
            msg_box.value = ""

    send_btn.on_click(_send)

    def _refresh(force: bool = False) -> None:
        now = _time.time()
        try:
            if force or now - state["_list_at"] >= RUNNING_REFRESH_S:
                state["rows"] = running_session_rows()
                state["_list_at"] = now
                if not state["picked"] and state["rows"]:
                    state["picked"] = state["rows"][0]["name"]
                _render_list()
                _render_chat()
            elif state["picked"] and now - state["_chat_at"] >= 10.0:
                _render_chat()
            roster = getattr(ctx, "roster_path", None)
            spec_path = getattr(ctx, "swarm_spec_path", None)
            rows = []
            if roster:
                rows = watched_rows_html(watch_overview(str(roster)))
            elif spec_path and Path(spec_path).exists():
                rows = watched_rows_html(watch_overview(load_spec(Path(spec_path))))
            table_html.value = ("<table><tr><td><b>member</b></td><td>pkg</td><td>role</td>"
                                "<td>model</td><td>state</td><td>last tool</td>"
                                "<td>commits</td><td>approvals</td></tr>"
                                + "".join(rows) + "</table>") if rows else ""
        except Exception as exc:  # a broken source must not kill the tab
            send_note.value = f"<span style='color:#b91c1c'>{_html.escape(str(exc))}</span>"

    refresh_btn.on_click(lambda _b: _refresh(force=True))

    def _shown() -> bool:
        """Whether Swarm Watch is the tab on screen: reading the sessions
        for a tab nobody looks at costs the login node for nothing."""
        try:
            tabs = getattr(ctx, "tabs_widget", None)
            title_of = getattr(ctx, "tab_indices_title_of", None)
            if tabs is None or not callable(title_of):
                return True
            return title_of(tabs.selected_index) == "Swarm Watch"
        except Exception:
            return True

    def _tick() -> None:
        if state["_stopped"]:
            return
        if _shown():
            _refresh()
        timer = threading.Timer(5.0, _tick)
        timer.daemon = True
        state["_timer"] = timer
        timer.start()

    _refresh(force=True)
    _tick()
    return root


# Self-register AFTER create_tab is defined (module end, the tab_review.py
# pattern, tab_review.py:1393-1397). An earlier attempt registered at line 56,
# before create_tab existed: the NameError was swallowed by the except and the
# tab silently never appeared in the registry — the exact defect the registry
# pin test caught.
try:  # pragma: no cover - import-side effect
    from delfin.dashboard.tab_registry import register_tab
    register_tab("swarm", "Swarm", create_tab, order=57)
except Exception:  # registration must never break the import
    pass

try:  # pragma: no cover - import-side effect
    from delfin.dashboard.tab_registry import register_tab
    register_tab("swarm-watch", "Swarm Watch", create_watch_tab, order=58)
except Exception:  # a watch-tab problem must never hide the manage tab
    pass
