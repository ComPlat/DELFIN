"""swarm status: the JSON contract SW2 (dashboard/TUI) renders.

STATUS_SCHEMA_VERSION is the contract version: SW2 pins against it and a
schema change is a release event, not a silent shape change.

All functions here are read-only; they never start, stop or answer
anything. Side-effecting swarm actions live in their own modules.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Callable

STATUS_SCHEMA_VERSION = 1

#: The JSON schema of ``swarm status --json``. Documented in place so a
#: consumer (SW2, SW3 lead) can validate what it gets.
STATUS_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["schema_version", "goal", "packages"],
    "properties": {
        "schema_version": {"type": "integer",
                           "const": STATUS_SCHEMA_VERSION},
        "goal": {"type": "string"},
        "packages": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name", "phase", "members"],
                "properties": {
                    "name": {"type": "string"},
                    # Phase of the package as a whole, from the members'
                    # states: starting -> building -> in_review -> qs_ok
                    # -> pr -> merged.
                    "phase": {"type": "string",
                              "enum": ["starting", "building", "in_review",
                                       "qs_ok", "pr", "merged", "blocked",
                                       "stopped"]},
                    "members": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "required": ["name", "role", "health"],
                            "properties": {
                                "name": {"type": "string"},
                                "role": {"type": "string",
                                         "enum": ["builder", "reviewer",
                                                  "qs"]},
                                # Health from delfin.agent.session_health.
                                "health": {"type": "string"},
                                "waits_on": {"type": ["string", "null"]},
                            },
                        },
                    },
                },
            },
        },
        # What the user must answer, aggregated across members.
        "waiting_for_user": {"type": "array", "items": {"type": "object"}},
    },
}

#: An example document, valid against STATUS_SCHEMA. Used by tests and as
#: documentation for SW2.
STATUS_EXAMPLE: dict[str, Any] = {
    "schema_version": STATUS_SCHEMA_VERSION,
    "goal": "Ship the swarm feature",
    "packages": [
        {
            "name": "core",
            "phase": "in_review",
            "members": [
                {"name": "s6", "role": "builder", "health": "active",
                 "waits_on": None},
                {"name": "s12", "role": "reviewer", "health": "active",
                 "waits_on": "builder commit"},
            ],
        },
    ],
    "waiting_for_user": [],
}

_EVENT_KINDS = frozenset({
    "phase_change", "approval_waiting", "member_stalled", "member_ended",
    "ci_red", "ci_green",
})

#: Subscribers receive every event dict; SW3/SW2 subscribe, nothing in
#: DELFIN depends on there being one.
_event_subscribers: list[Callable[[dict], None]] = []


def emit_event(event: dict[str, Any]) -> None:
    """Send one swarm event to every subscriber. Never raises."""
    event = dict(event)
    event.setdefault("kind", "")
    if event["kind"] not in _EVENT_KINDS:
        # An unknown kind is a programming error upstream; deliver it
        # anyway but tagged, so a subscriber can reject it.
        event["kind_invalid"] = True
    for sub in list(_event_subscribers):
        try:
            sub(event)
        except Exception:  # a subscriber must not break the swarm
            continue


def subscribe(callback: Callable[[dict], None]) -> Callable[[], None]:
    """Register an event subscriber; returns the unsubscribe callable."""
    _event_subscribers.append(callback)

    def _off() -> None:
        try:
            _event_subscribers.remove(callback)
        except ValueError:
            pass

    return _off


#: A member the job started this long ago with no session record of its own
#: has ended (a session announces itself within seconds of its start).
PRESENCE_GRACE_S = 180.0

#: Session-health verdicts (delfin.agent.session_health.classify) as the
#: member health this report shows.
_VERDICT_HEALTH = {
    "working": "active", "waiting_on_model": "active",
    "waiting_on_person": "waiting", "stalled": "stalled",
    "erroring": "stalled", "paused": "stalled", "stopped": "idle",
    "gone": "ended",
}


def _member_names(spec: dict[str, Any]) -> set[str]:
    return {f"{pkg['name']}-{role}" for pkg in spec.get("packages", [])
            for role in (pkg.get("roles") or {"builder": True})}


def running_job_for(spec: dict[str, Any],
                    jobs: list[dict[str, Any]] | None) -> dict | None:
    """The RUNNING swarm job (a row of ``job.jobs()``) that runs *spec*:
    same goal, a member of the spec among its members. The newest wins.
    A pending, unknown or ended job is none: its members have not started
    or their records are history (D8)."""
    names = _member_names(spec)
    goal = str(spec.get("goal") or "")
    for row in reversed(list(jobs or [])):
        if str(row.get("state") or "") != "RUNNING":
            continue
        if str(row.get("goal") or "") != goal:
            continue
        if names & set(row.get("members") or []):
            return row
    return None


def _same_path(a: str, b: Path) -> bool:
    try:
        return Path(a).resolve() == b.resolve()
    except (OSError, ValueError):
        return False


def _health_of(rec: dict[str, Any], now: float) -> tuple[str, Any]:
    """(health, waits_on) of one session-health record."""
    if rec.get("state"):  # a record that already carries its health
        return str(rec["state"]), rec.get("waits_on")
    from .. import session_health as _sh
    verdict = _sh.classify(rec, now, {})["verdict"]
    health = _VERDICT_HEALTH.get(verdict, "active")
    waits_on = None
    if health == "waiting":
        ids = [str(i) for i in rec.get("open_approvals") or []]
        waits_on = "approval " + ", ".join(ids) if ids else "approval"
    return health, waits_on


def _in_job(name: str, worktree: str, row: dict[str, Any],
            records: dict[str, dict], now: float) -> tuple[str, Any]:
    """(health, waits_on) of member *name* of the running job *row*.

    The job's progress file says whether it started the member; the
    member's own session record says what it does now, and counts only
    when it is on the job's node in the member's worktree -- a same-named
    session elsewhere is not this member (D8)."""
    prog = row.get("progress") or {}
    if name not in (prog.get("started") or []):
        return "starting", None
    if prog.get("phase") in ("ended", "failed"):
        return "ended", None
    node = str(prog.get("node") or row.get("node") or "").split(".")[0]
    rec = records.get(name)
    if rec is not None and node and \
            str(rec.get("host") or "").split(".")[0] == node and \
            _same_path(str(rec.get("workspace") or ""),
                       Path(str(row.get("repo") or "")) / worktree):
        return _health_of(rec, now)
    started_at = (prog.get("started_at") or {}).get(name) or \
        prog.get("at") or 0
    try:
        age = now - float(started_at)
    except (TypeError, ValueError):
        age = 0.0
    return ("started", None) if age < PRESENCE_GRACE_S else ("ended", None)


def collect_status(spec: dict[str, Any], *,
                   jobs: list[dict[str, Any]] | None = None,
                   now: float | None = None) -> dict[str, Any]:
    """Build a status report for a spec: phases + member health.

    Read-only. A member of a RUNNING swarm job for this spec (*jobs*: the
    rows ``delfin.agent.swarm.job.jobs()`` gives; None asks it) is
    'starting' until the job starts it, 'started' until its session
    announces itself, then its session health ('active', 'waiting',
    'idle', 'stalled'), and 'ended' when its session is gone. The job's
    run id (``swarm-<job>``) is this run's id. Outside a job, a member
    counts only with a record in the current swarm_memory run
    (DELFIN_SWARM_RUN). D8 stays fail-closed: with neither, every member
    reports 'starting' -- SW2 must not read stale runs as healthy.
    """
    from .. import session_health as _sh
    from .. import swarm_memory

    now = time.time() if now is None else now
    if jobs is None:
        try:
            from . import job as _job
            jobs = _job.jobs()
        except Exception:
            jobs = []
    row = running_job_for(spec, jobs)
    # Scope health to THIS swarm run: a same-named member of an earlier
    # run must not leak its final state into a fresh one (D8).
    run_id = swarm_memory.current_run_id()
    if not run_id and row is not None:
        run_id = str((row.get("progress") or {}).get("run") or "")
    members_by_pkg: dict[str, list[dict[str, Any]]] = {}
    waiting: list[dict[str, Any]] = []
    records: dict[str, dict] = {}
    if run_id or row is not None:
        try:
            for r in _sh.collect(with_git=False):
                key = r.get("name") or r.get("session")
                if key:
                    records[str(key)] = r
        except Exception:
            records = {}
    entries = []
    if run_id:
        try:
            entries = swarm_memory.read(run_id)
        except Exception:
            entries = []
    from .launch import WORKTREE_DIR
    for pkg in spec.get("packages", []):
        roles = pkg.get("roles") or {"builder": True}
        members = []
        for role in roles:
            name = f"{pkg['name']}-{role}"
            rec = records.get(name)
            if row is not None and name in (row.get("members") or []):
                health, waits_on = _in_job(
                    name, f"{WORKTREE_DIR}/{name}", row, records, now)
            # D8: outside a job a health record only counts for THIS run.
            # swarm_memory keeps the run id each entry was written under;
            # a record from an earlier run (or no record) means 'starting'.
            elif rec is not None and run_id and any(
                    getattr(e, "author", "") == name for e in entries):
                health, waits_on = _health_of(rec, now)
            else:
                health, waits_on = "starting", None
            if waits_on:
                waiting.append({"member": name, "waits_on": waits_on})
            members.append({"name": name, "role": role, "health": health,
                            "waits_on": waits_on})
        members_by_pkg[pkg["name"]] = members
    report = {
        "schema_version": STATUS_SCHEMA_VERSION,
        "goal": spec.get("goal", ""),
        "packages": [],
        "waiting_for_user": waiting,
    }
    for pkg in spec.get("packages", []):
        members = members_by_pkg[pkg["name"]]
        report["packages"].append({
            "name": pkg["name"],
            "phase": _phase_of(members),
            "members": members,
        })
    return report


def _phase_of(members: list[dict[str, Any]]) -> str:
    """The package phase from member health states."""
    states = {m["health"] for m in members}
    if "stalled" in states:
        return "blocked"
    if "waiting" in states:
        return "in_review"
    if "ended" in states:
        return "merged"
    if states & {"active", "idle"}:
        return "building"
    if states & {"starting", "started"}:
        return "starting"
    return "stopped"


def to_json(report: dict[str, Any]) -> str:
    """Serialise a status report for ``swarm status --json``."""
    return json.dumps(report, indent=2, sort_keys=True)


def validate_report(report: dict[str, Any]) -> list[str]:
    """Structural self-check of a status report against STATUS_SCHEMA.

    Deliberately not a full JSON-Schema engine: a closed set of checks
    over the fields SW2 actually renders, so the contract cannot drift
    silently. Returns a list of violations (empty = valid).
    """
    problems: list[str] = []
    for key in STATUS_SCHEMA["required"]:
        if key not in report:
            problems.append(f"missing required key: {key}")
    if report.get("schema_version") != STATUS_SCHEMA_VERSION:
        problems.append(
            f"schema_version must be {STATUS_SCHEMA_VERSION}")
    for pkg in report.get("packages", []):
        for key in ("name", "phase", "members"):
            if key not in pkg:
                problems.append(f"package missing key: {key}")
        if pkg.get("phase") not in \
                STATUS_SCHEMA["properties"]["packages"]["items"][
                    "properties"]["phase"]["enum"]:
            problems.append(f"unknown phase: {pkg.get('phase')!r}")
        for m in pkg.get("members", []):
            for key in ("name", "role", "health"):
                if key not in m:
                    problems.append(f"member missing key: {key}")
    return problems
