"""Lab audit: DELFIN's hash-chained log, one per installation.

Input: an event name and its fields. Output: one record in
``<lab_dir>/audit.log`` through ``delfin.agent.audit_log.append`` (hash
chain, keyed MAC where the account key is available, rotation, the
administrator's copy), and for the alarming events one security event,
which reaches the containment panel and the attention inbox.

A write is three records in order: ``intent`` (before the driver is
called), ``result`` (what the driver said), ``readback`` (what the device
reports afterwards). A crash between intent and result leaves an intent
without a result, which is what a reader needs to see.

Never raises: an audit failure must not change the kernel's decision. The
audit module of the agent records its own ``audit_write_failed`` event.
"""

from __future__ import annotations

from typing import Any

from . import store as _store

#: Security-event kinds the lab records; labelled in
#: ``delfin.agent.security_events._KINDS`` and pinned by a test.
LAB_EVENT_KINDS = frozenset({
    "lab_refused", "lab_safe_state", "lab_estop", "lab_estop_cleared",
    "lab_deadman", "lab_readback_mismatch", "lab_approval_refused",
    "lab_manifest_refused",
})

#: Audit events (the ``event`` field of a record with ``tool: lab``).
LAB_AUDIT_EVENTS = frozenset({
    "intent", "result", "readback", "refused", "safe_state", "estop",
    "estop_cleared", "deadman", "approval_requested", "approval_granted",
    "approval_consumed", "latched", "cleared", "released", "manifest_pinned",
    "run_step",
})


def record(event: str, **fields: Any) -> None:
    try:
        from delfin.agent import audit_log
        rec = {"tool": "lab", "event": event, "pid": __import__("os").getpid()}
        rec.update({k: v for k, v in fields.items() if v is not None})
        audit_log.append(rec, log_path=_store.ensure_lab_dir() / "audit.log",
                         session_id=str(fields.get("session") or ""))
    except Exception:
        pass


def alert(kind: str, detail: str, *, blocked: bool = True) -> None:
    try:
        from delfin.agent import security_events
        security_events.record(kind, "lab", str(detail)[:300], blocked=blocked)
    except Exception:
        pass


def read_all() -> list[dict]:
    """Every lab audit record (active log), oldest first."""
    import json
    out = []
    try:
        path = _store.lab_dir() / "audit.log"
        for line in path.read_text(encoding="utf-8").splitlines():
            try:
                out.append(json.loads(line))
            except ValueError:
                continue
    except OSError:
        pass
    return out
