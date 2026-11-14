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
from pathlib import Path
from typing import Any

from ..swarm import launch as _launch
from ..swarm import observe as _observe
from ..swarm import spec as _spec
from ..swarm import status as _status
from ..swarm import triage as _triage
from .. import session_messages as _session_messages


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
    # A strict boolean: a truthy string or int must not become a
    # confirmation (F1). The protocol gate enforces boolean too; this is
    # the wrapper's own gate, so the module keeps its double-gate
    # promise independently of the boundary.
    confirmed = params.get("confirm") is True
    # Any confirm that is PRESENT and not strictly true is an invalid
    # confirmation attempt: refuse it with a reason, even on a dry run
    # (F1 -- the TUI shows what was wrong, the plan is still returned).
    bad_confirm = params.get("confirm") is not None and not confirmed
    reason = None
    if bad_confirm:
        reason = ("confirm must be a boolean true (got "
                  f"{params.get('confirm')!r})")
    elif not dry_run:
        reason = ("execute needs confirm=true: the user has not "
                  "approved this plan in the UI")
    if dry_run or not confirmed:
        plan = _launch.plan_launch(spec)
        return {
            "dry_run": True,
            "goal": plan["goal"],
            "members": plan["members"],
            "dry_run_text": _launch.render_dry_run(plan),
            "started": [],
            "refused": reason,
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


def swarm_message(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.message: send one plain partner message to a member.

    Delivers through session_messages (no shell, no approval surface);
    the lead's swarm_message tool dispatches here. The lead session
    confirms state-changing sends with the user before dispatch; this
    wrapper is transport only.
    """
    to_key = str(params.get("to") or "")
    text = str(params.get("text") or "")
    if not to_key:
        raise SwarmRpcError("to is required")
    if not text:
        raise SwarmRpcError("text is required")
    sent = _session_messages.send(to_key, text)
    return {"sent": bool(sent), "to": to_key}


def swarm_prcheck(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.prcheck: read-only PR status; runs no tests in-process.

    Policy (operator, 2026-11-12, after the login-node suspension
    #640): the pre-PR check runs tests. On a login node the in-process
    route would bypass the compute-node rule hpc_guard.MOVED_COMMANDS
    enforces for the CLI (test-on-slurm, cli.py:6205). So this handler
    NEVER runs the check itself: on a login node it answers with the
    structured unavailability and names the compute-node route
    (delfin prcheck on a compute node, or delfin test-on-slurm);
    outside a login node it delegates to the read-only git-read report
    pre_pr_check.run (no tests are executed by it).

    The checked tree is the WORKTREE the request names (params
    "worktree", required) -- inside the agent service the process cwd
    is the service's runtime dir, not the swarm repo, so a '.' default
    would check the wrong tree. Refused when the path does not exist:
    an absent worktree is a wrong parameter, not an empty report.
    """
    pr = params.get("pr")
    if pr is None:
        raise SwarmRpcError("pr is required")
    if not isinstance(pr, int) or isinstance(pr, bool):
        raise SwarmRpcError("pr cannot be bool; pass an integer")
    if pr <= 0:
        raise SwarmRpcError("pr must be a positive integer")
    worktree = str(params.get("worktree") or "")
    if not worktree:
        raise SwarmRpcError(
            "worktree is required: the repository worktree to check "
            "(the service process cwd is not the swarm repo)")
    if not Path(worktree).expanduser().is_dir():
        raise SwarmRpcError(f"worktree does not exist: {worktree!r}")
    from .. import pre_pr_check  # delfin.agent.pre_pr_check
    from delfin.agent import calc_policy as _calc_policy
    host = _calc_policy.host_kind()
    if host == "login_node":
        return {"pr": pr, "checks": [], "available": False,
                "reason": "login node: the pre-PR check runs tests, so it "
                          "must not run here in-process; use 'delfin "
                          "prcheck' on a compute node or 'delfin "
                          "test-on-slurm' (slurm_tests.py)"}
    # Read-only git reads only (pre_pr_check.run, pre_pr_check.py:314):
    # selects the tests, never executes them. Nothing is written to the
    # checked tree.
    report = pre_pr_check.run(worktree, base="origin/main")
    return {"pr": pr,
            "checks": [{"blocker": b.kind, "detail": b.detail}
                       for b in report.blockers],
            "available": report.verdict == "READY",
            "reason": (f"pre-PR check {report.verdict} on host {host} "
                       f"({len(report.changed_files)} changed files, "
                       f"{len(report.selected_tests)} tests selected)")}


def swarm_observe_overview(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.observe.overview: one row per member (read-only watch).

    ``source`` may be an SW1 spec object or a roster-style dict; a
    ``sourcePath`` names a roster file (name -> {package, role}) so an
    operator wave with no spec shows too. Pure read from observe.py --
    nothing here starts, pauses, answers or prompts a session.
    """
    source = params.get("source")
    if source is None and params.get("sourcePath"):
        source = str(params.get("sourcePath"))
    if source is None:
        raise SwarmRpcError("missing source: pass source (spec/roster dict) "
                            "or sourcePath (roster file path)")
    if not isinstance(source, (dict, str)):
        raise SwarmRpcError("source must be a spec/roster dict or a path")
    return {"rows": _observe.swarm_overview(source)}


def swarm_observe_feed(params: dict[str, Any]) -> dict[str, Any]:
    """swarm.observe.feed: one member's recent steps, newest last.

    Ownership -- a feed is shown only to the user who owns the session --
    is enforced before dispatch, in the supervisor (an agent's own
    connection gets NOT_ALLOWED). This thin wrapper only reads
    observe.member_feed and re-labels the list under ``entries``.
    """
    session = str(params.get("session") or "")
    if not session:
        raise SwarmRpcError("session is required")
    since = params.get("since")
    since = int(since) if since is not None else 0
    return {"entries": _observe.member_feed(session, since)}
