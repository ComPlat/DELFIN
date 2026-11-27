"""The lab tools an agent calls: single writes, and a whole planned run.

Input: a tool name, its arguments as the model sent them, and the
session's ``KitToolPermissions``. Output: a JSON string.

=================  ==========================================================
lab_devices        devices, capabilities, levels, fresh state, lock, latch
lab_read           one device's readings with their age
lab_set            one write through ``Kernel.set``; a refusal carries its
                   ``code`` and exact ``reason`` (and ``approval_request``)
lab_plan_check     a list of steps judged like ``lab_set`` with no write
lab_estop          global or one module; clearing is not a tool
lab_run            a whole workflow document: judged by ``lab.plan``, then
                   run by the ONE coordinator; nothing is written by this
                   tool itself
=================  ==========================================================

THE TWO HALVES, AND WHY THEY ARE TOOLS RATHER THAN ONE. ``lab_plan_check``
answers what the kernel would say about a list of *writes*. A workflow is
more than its writes -- it transfers a sample between modules and measures
it -- and a run cannot be judged by the write check alone (see
``delfin/lab/plan.py`` for the measured reason). So ``lab_run`` with
``dry_run`` true is the planner's check, and ``lab_run`` without it runs
what that check passed. The validation happens in BOTH cases and the run
refuses a document that did not pass, because an advisory check on the path
to a hotplate is a comment.

The actor is the person (``store.person``) and the session
(``permissions.task_session_id``, else this process). Every call is a
heartbeat for that actor. Arguments the schema does not name are ignored:
nothing in a call (``override``, ``max``, ``force``, a claimed approval)
reaches a limit, because the kernel reads limits from the pinned manifest
only. The permission mode is not passed on. ``lab_run`` does not bypass the
kernel either: the coordinator reaches every device through an edge whose
own kernel verifies the envelope, the manifest pin and its replay window.
"""

from __future__ import annotations

import json
import os

from . import store as _store

TOOL_NAMES = frozenset({"lab_devices", "lab_read", "lab_set",
                        "lab_plan_check", "lab_estop", "lab_run"})


def _actor(permissions):
    from .safety import Actor
    sid = str(getattr(permissions, "task_session_id", "") or "") \
        or f"engine-{os.getpid()}"
    return Actor.here(sid)


def _lab_run(kernel, actor: dict, args: dict) -> dict:
    """Judge a workflow document, and run it when ``dry_run`` is false."""
    from . import plan as _plan
    doc = args.get("workflow")
    if isinstance(doc, str):
        try:
            doc = json.loads(doc)
        except ValueError:
            return {"error": "workflow must be a JSON object (or the path to "
                             "one, as `path`)"}
    if doc is None and args.get("path"):
        try:
            doc = _plan.load(str(args["path"]))
        except ValueError as exc:
            return {"error": str(exc)}
    if not isinstance(doc, dict):
        return {"error": "lab_run needs a workflow document: an object with "
                         "`id`, `samples` and `steps`"}
    if kernel.manifest is None:
        return {"error": kernel.load_error or "no pinned manifest"}

    proposal = _plan.validate(kernel.manifest, doc, kernel)
    out = proposal.as_dict()
    if args.get("dry_run", True):
        out["ran"] = False
        out["note"] = ("nothing was written: this is the plan check. Call "
                       "again with dry_run=false to run it.")
        return out
    if not proposal.ok:
        out["ran"] = False
        out["note"] = "the plan was refused, so nothing ran"
        return out

    run_id = str(args.get("run_id") or "").strip() \
        or f"run-{os.urandom(4).hex()}"
    out["run_id"] = run_id
    try:
        cp = _plan.run(kernel.manifest, proposal, run_id,
                       person=actor["person"], coordinator=actor["session"])
    except ValueError as exc:                    # the refusal just described
        out["ran"] = False
        out["error"] = str(exc)
        return out
    out["ran"] = True
    out["status"] = cp.get("status")
    out["done"] = cp.get("done")
    out["deviations"] = cp.get("deviations")
    out["measurements"] = cp.get("measurements")
    # Where the run's own record is. `record.json` and the Chemotion export
    # are written when a run ENDS (completed, failed, abandoned) -- a run
    # stopped for a person has a checkpoint and nothing else, which is what
    # `checkpoint.json` says rather than the record being silently absent.
    from .workflow import run_dir
    here = run_dir(run_id)
    out["run_dir"] = str(here)
    out["record_written"] = (here / "record.json").is_file()
    # A run that stopped for a person is reported as stopped, with the id
    # that person needs. Saying "done" here would be the one lie that
    # matters: a hazardous step has NOT happened.
    pend = cp.get("pending") or {}
    if cp.get("status") == "waiting_approval" and pend.get("approval"):
        out["approval_request"] = pend["approval"]
        out["waiting_step"] = pend.get("step")
        out["ask"] = (f"step {pend.get('step')!r} needs a person: they run "
                      f"`delfin-agent lab approve {pend['approval']}` and then "
                      "the run is resumed")
        out["note"] = ("NOT finished: the hazardous step has not happened, so "
                       "no run record exists yet. Only the checkpoint does.")
    return out


def execute(name: str, arguments: dict, permissions=None) -> str:
    from .safety import Kernel
    args = arguments if isinstance(arguments, dict) else {}
    if not _store.lab_configured():
        return json.dumps({"error": "no lab manifest on this installation "
                           f"({_store.manifest_path().name} is missing)"})
    kernel = Kernel.load()
    actor = _actor(permissions)
    if name == "lab_estop":
        from . import estop
        out = estop.trigger(module=args.get("module") or None,
                            reason=str(args.get("reason") or "agent e-stop"),
                            by=actor.person, kernel=kernel)
        out["note"] = ("writes are refused until a person runs "
                       "`delfin-agent lab estop --clear`")
        return json.dumps(out)
    if kernel.manifest is None:
        return json.dumps({"error": kernel.load_error})
    kernel.heartbeat(actor)
    if name == "lab_devices":
        from . import estop
        return json.dumps({"devices": kernel.devices(),
                           "estop": estop.status()}, default=str)
    if name == "lab_read":
        return json.dumps(kernel.read(str(args.get("device") or ""), actor),
                          default=str)
    if name == "lab_set":
        out = kernel.set(actor, str(args.get("device") or ""),
                         str(args.get("channel") or ""), args.get("value"),
                         args.get("unit"), args.get("approval") or None)
        if not out.get("ok"):
            out = {"error": f"refused ({out['code']}): {out['reason']}", **out}
        return json.dumps(out, default=str)
    if name == "lab_plan_check":
        steps = args.get("steps")
        if isinstance(steps, str):
            try:
                steps = json.loads(steps)
            except ValueError:
                return json.dumps({"error": "steps must be a list of "
                                   "{device, channel, value, unit}"})
        return json.dumps(kernel.plan_check(actor, steps or []), default=str)
    if name == "lab_run":
        return json.dumps(_lab_run(kernel, actor.as_dict(), args), default=str)
    return json.dumps({"error": f"unknown lab tool {name!r}"})

