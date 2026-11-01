"""Service-level swarm RPC methods (package SW2b).

The methods the TUI's /swarm view calls. Every function here is a THIN
wrapper over SW1 (delfin/agent/swarm/) -- JSON-able params in, JSON-able
result out, no orchestration logic -- so the CLI (cli_swarm.py), the
dashboard tab (tab_swarm.py) and this service surface call the same
library code and cannot drift apart.

Nothing here starts a session on its own: swarm.start runs the dry-run
plan by default and refuses execute without the TUI's explicit
confirmation flag (SW1's start_members itself refuses without
execute=True, launch.py:116-120; the double gate means params alone can
never flip the side effect on). swarm.pause/resume answer with the
structured refusal until SW1 ships those verbs (tab_swarm.py:13-15
renders the buttons only when they appear).

SW1's validation errors (SpecError, LaunchError, TriageError -- all
RuntimeError subclasses or ValueError) pass through: the service layer
turns them into protocol error answers, never tracebacks.
"""

from __future__ import annotations

import json
from typing import Any

from ..swarm import launch as _launch
from ..swarm import spec as _spec
from ..swarm import status as _status
from ..swarm import triage as _triage


class SwarmRpcError(ValueError):
    """A swarm request cannot be served, with the reason in prose."""


def _spec_of(params: dict[str, Any]) -> dict[str, Any]:
    """The spec a request carries: a dict under "spec", or YAML/JSON
    text under "specYaml". A missing spec is an error, never a silent
    default -- a start with no spec would plan nothing and say nothing."""
    spec = params.get("spec")
    if spec is None:
        text = params.get("specYaml")
        if not text or not isinstance(text, str):
            raise SwarmRpcError("missing spec: pass spec (object) or "
                                "specYaml (swarm.yaml text)")
        try:
            spec = json.loads(text)
        except json.JSONDecodeError:
            import yaml
            try:
                spec = yaml.safe_load(text)
            except yaml.YAMLError as exc:
                raise SwarmRpcError(f"specYaml is not valid YAML: {exc}")
    if not isinstance(spec, dict):
        raise SwarmRpcError("spec must be an object")
    return spec


def swarm_status(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.status: the read-only report SW1's collect_status builds."""
    spec = _spec_of(params)
    try:
        _spec.validate_spec(spec)
    except _spec.SpecError as exc:
        raise SwarmRpcError(str(exc)) from None
    report = _status.collect_status(spec)
    problems = _status.validate_report(report)
    if problems:
        # SW1's own self-check failed: surface it, never serve a
        # silently drifted report (the CLI does the same, cli_swarm.py:165).
        raise SwarmRpcError(
            "status report failed its schema check: " + "; ".join(problems))
    return report


def swarm_plan(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.plan: the launch plan + the dry-run text a user confirms."""
    spec = _spec_of(params)
    plan = _launch.plan_launch(spec)
    plan["dry_run_text"] = _launch.render_dry_run(plan)
    return plan


def swarm_validate(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.validate: is this spec startable? Errors, never raised:
    the TUI shows the reasons next to the editor."""
    spec = _spec_of(params)
    try:
        _spec.validate_spec(spec)
    except _spec.SpecError as exc:
        return {"valid": False, "errors": [str(exc)]}
    return {"valid": True, "errors": []}


def swarm_start(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.start: dry-run plan by default; the real start only with
    the explicit confirmation the TUI sends after the user approves the
    plan (double gate: this flag AND SW1's own execute=True)."""
    spec = _spec_of(params)
    dry_run = params.get("dryRun", True)
    if dry_run is None:
        dry_run = True
    confirmed = bool(params.get("confirm", False))
    if dry_run or not confirmed:
        plan = _launch.plan_launch(spec)
        return {
            "dry_run": True,
            "goal": plan["goal"],
            "members": plan["members"],
            "dry_run_text": _launch.render_dry_run(plan),
            "started": [],
            "refused": None if dry_run else
            "execute needs confirm=true: the user has not approved "
            "this plan in the UI",
        }
    try:
        return _launch.start_members(spec, execute=True)
    except _launch.LaunchError as exc:
        raise SwarmRpcError(str(exc)) from None


def swarm_pause(_params: dict[str, Any]) -> dict[str, Any]:
    """swarm.pause: SW1 ships no pause verb; structured refusal (SW2a
    precedent) -- never a crash, never an invented side effect."""
    return {"available": False, "reason": "SW1 ships no pause verb; "
            "the swarm tab renders pause only when SW1 grows it "
            "(tab_swarm.py probes for it)"}


def swarm_resume(_params: dict[str, Any]) -> dict[str, Any]:
    """swarm.resume: see swarm_pause."""
    return {"available": False, "reason": "SW1 ships no resume verb; "
            "the swarm tab renders resume only when SW1 grows it "
            "(tab_swarm.py probes for it)"}


def swarm_triage_list(params: dict[str, Any]) -> list[dict[str, Any]]:
    """swarm.triage.list: the queue sorted, escalations first."""
    queue = params.get("queue")
    if not isinstance(queue, list):
        raise SwarmRpcError("queue must be a list of items")
    return _triage.triage_list(queue)


def swarm_triage_answer(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.triage.answer: answer one item on the user's behalf.

    NO verdict param reaches SW1: triage_answer re-derives the verdict
    from the item (triage.py:108-121, finding T4/D13) and fails closed
    on an escalation item, so a crafted verdict dict cannot widen the
    answer -- the service surface passes only item + answer.
    """
    item = params.get("item")
    if not isinstance(item, dict):
        raise SwarmRpcError("item must be an object")
    answer = str(params.get("answer") or "")
    if not answer:
        raise SwarmRpcError("answer needs text")
    try:
        return _triage.triage_answer(item, {}, answer=answer)
    except _triage.TriageError as exc:
        raise SwarmRpcError(str(exc)) from None
