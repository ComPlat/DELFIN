"""The five agent tools, executed through the safety kernel.

Input: a tool name, its arguments as the model sent them, and the
session's ``KitToolPermissions``. Output: a JSON string.

=================  ==========================================================
lab_devices        devices, capabilities, levels, fresh state, lock, latch
lab_read           one device's readings with their age
lab_set            one write through ``Kernel.set``; a refusal carries its
                   ``code`` and exact ``reason`` (and ``approval_request``)
lab_plan_check     a list of steps judged like ``lab_set`` with no write
lab_estop          global or one module; clearing is not a tool
=================  ==========================================================

The actor is the person (``store.person``) and the session
(``permissions.task_session_id``, else this process). Every call is a
heartbeat for that actor. Arguments the schema does not name are ignored:
nothing in a call (``override``, ``max``, ``force``, a claimed approval)
reaches a limit, because the kernel reads limits from the pinned manifest
only. The permission mode is not passed on.
"""

from __future__ import annotations

import json
import os

from . import store as _store

TOOL_NAMES = frozenset({"lab_devices", "lab_read", "lab_set",
                        "lab_plan_check", "lab_estop"})


def _actor(permissions):
    from .safety import Actor
    sid = str(getattr(permissions, "task_session_id", "") or "") \
        or f"engine-{os.getpid()}"
    return Actor.here(sid)


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
    return json.dumps({"error": f"unknown lab tool {name!r}"})
