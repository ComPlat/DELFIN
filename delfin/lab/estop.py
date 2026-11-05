"""Emergency stop: global or per module, cleared only by a person.

``trigger`` is callable from anywhere -- the CLI, the ``lab_estop`` tool,
a dashboard/TUI/phone handler -- and does its work in the CALLING process:
it does not ask the agent loop, which may be the thing that hangs.

Order: (1) the stop flag is written under ``flock``, so every write that
checks after this instant is refused; (2) every device in scope is driven
to its safe state, continuing past a device that fails; (3) locks and
dead-man obligations in scope are revoked and pending approvals cancelled;
(4) audit and a security event. The e-stop daemon re-asserts the safe
state while the flag stands.

``clear`` lifts the flag and nothing else: devices stay in their safe
state (a reset permits restarting; it does not restart). It is refused in
a process an agent session started.
"""

from __future__ import annotations

from typing import Optional

from . import audit as _audit
from . import store as _store


class EstopError(PermissionError):
    pass


def trigger(*, module: Optional[str] = None, reason: str = "",
            by: Optional[str] = None, kernel=None) -> dict:
    from .safety import Kernel
    by = by or _store.person()
    now = _store.now()
    entry = {"by": by, "reason": reason or "e-stop", "ts": now}
    kernel = kernel or Kernel.load()
    # A module stop must name a module of the pinned manifest: a misspelled
    # scope stopped nothing and still answered ok. Refused before any flag is
    # written, for every caller (CLI, lab_estop tool, coordinator envelope).
    if module:
        known = sorted((kernel.manifest.modules if kernel.manifest else {}) or {})
        if module not in known:
            return {"ok": False, "code": "unknown_module", "scope": module,
                    "known": known, "stopped": [], "not_confirmed": [],
                    "error": f"no module {module!r} in the lab manifest"
                             + (f"; known: {', '.join(known)}" if known else "")
                             + "; nothing was stopped (a global e-stop needs no module)"}
    with _store.state() as doc:
        if module:
            doc["estop"]["modules"][module] = entry
        else:
            doc["estop"]["global"] = entry
    stopped, failed = [], []
    if kernel.manifest is not None:
        in_scope = [d for d in kernel.manifest.devices.values()
                    if (module is None or d.module == module)
                    and (kernel.node is None or d.node == kernel.node)]
        for d in in_scope:
            if not d.safe_state:
                continue
            try:
                res = kernel.safe(d.id, f"e-stop: {entry['reason']}",
                                  latch=False)
                (stopped if res.get("ok") else failed).append(d.id)
            except Exception as exc:
                failed.append(f"{d.id}: {exc}")
        ids = {d.id for d in in_scope}
        with _store.state() as doc:
            for dev in list(doc["locks"]):
                if dev in ids:
                    doc["locks"].pop(dev, None)
            for dev in list(doc["active"]):
                if dev in ids:
                    doc["active"].pop(dev, None)
            for mod in list(doc["module_locks"]):
                if module is None or mod == module:
                    doc["module_locks"].pop(mod, None)
            for rec in doc["approvals"].values():
                if rec["status"] in ("pending", "granted") and rec["device"] in ids:
                    rec["status"] = "cancelled"
    _audit.record("estop", module=module, by=by, reason=entry["reason"],
                  stopped=stopped, failed=failed or None)
    _audit.alert("lab_estop", f"{'module ' + module if module else 'global'} "
                 f"e-stop by {by}: {entry['reason']}; safe: {len(stopped)}, "
                 f"not confirmed: {len(failed)}")
    return {"ok": not failed, "scope": module or "global", "stopped": stopped,
            "not_confirmed": failed}


def clear(*, module: Optional[str] = None, by: Optional[str] = None) -> dict:
    from delfin.agent.approval_guard import inside_agent_session
    if inside_agent_session():
        _audit.alert("lab_estop_cleared", "clear refused: agent session")
        raise EstopError("an e-stop is cleared by a person, not by a process "
                         "an agent session started")
    by = by or _store.person()
    with _store.state() as doc:
        if module:
            had = doc["estop"]["modules"].pop(module, None)
        else:
            had = doc["estop"].get("global")
            doc["estop"]["global"] = None
    _audit.record("estop_cleared", module=module, by=by)
    _audit.alert("lab_estop_cleared", f"{'module ' + module if module else 'global'}"
                 f" e-stop cleared by {by}; devices stay in their safe state",
                 blocked=False)
    return {"ok": True, "cleared": bool(had), "scope": module or "global"}


def clear_latch(device: str, *, by: Optional[str] = None) -> dict:
    from delfin.agent.approval_guard import inside_agent_session
    if inside_agent_session():
        raise EstopError("a latched device is cleared by a person")
    by = by or _store.person()
    with _store.state() as doc:
        had = doc["latched"].pop(device, None)
    _audit.record("cleared", device=device, by=by)
    return {"ok": True, "cleared": bool(had), "device": device}


def status() -> dict:
    snap = _store.snapshot()
    return {"global": snap["estop"].get("global"),
            "modules": snap["estop"].get("modules") or {},
            "latched": snap["latched"]}
