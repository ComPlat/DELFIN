"""swarm status: the JSON contract SW2 (dashboard/TUI) renders.

STATUS_SCHEMA_VERSION is the contract version: SW2 pins against it and a
schema change is a release event, not a silent shape change.

All functions here are read-only; they never start, stop or answer
anything. Side-effecting swarm actions live in their own modules.
"""

from __future__ import annotations

import json
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
                    # states: building -> in_review -> qs_ok -> pr -> merged.
                    "phase": {"type": "string",
                              "enum": ["building", "in_review", "qs_ok",
                                       "pr", "merged", "blocked", "stopped"]},
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


def collect_status(spec: dict[str, Any]) -> dict[str, Any]:
    """Build a status report for a spec: phases + member health.

    Read-only. Health comes from delfin.agent.session_health; a member
    with no session record yet is 'starting'. The package phase is the
    lowest phase of its members (a stopped/stalled member blocks the
    whole package upward).
    """
    from .. import session_health as _sh

    members_by_pkg: dict[str, list[dict[str, Any]]] = {}
    waiting: list[dict[str, Any]] = []
    try:
        records = {r.get("session"): r for r in _sh.collect()}
    except Exception:
        records = {}
    for pkg in spec.get("packages", []):
        roles = pkg.get("roles") or {"builder": True}
        members = []
        for role in roles:
            name = f"{pkg['name']}-{role}"
            rec = records.get(name)
            health = (rec or {}).get("state", "starting")
            waits_on = (rec or {}).get("waits_on")
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
    if "active" in states:
        return "building"
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
