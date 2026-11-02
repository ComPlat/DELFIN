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
from delfin.agent.swarm.spec import SpecError, load_spec, validate_spec


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
