"""Visible containment surface: a thread-safe record of every safety decision
the agent's permission gate makes (a denied command, a blocked secret/script
payload, a self-modification attempt, an out-of-workspace access).

Blocks are otherwise only visible as a tool error to the model and a line in
the 0600 audit log. This module keeps a small in-process ring buffer so the
dashboard can show a prominent "🛡 Containment" panel — so an attempt to break
out is immediately VISIBLE to the user, not buried.

The ring buffer holds 200 events and lives in memory only, which means a
dashboard that is open sees everything and everyone else sees nothing: the
scheduler daemon and the job monitor record containment events inside child
processes no panel is attached to, and every one of them — including a
refusal the agent tried to walk around, the loudest signal the system has —
died with the process. A BLOCKED event therefore also goes to the durable
attention inbox, deduped per (kind, tool) so a loop cannot flood it. Kinds
that are the user's own decision or a routine grant stay panel-only.

Pure data + rendering otherwise; no I/O on the hot path beyond a bounded
deque append and (at most once per kind/tool) one inbox append.
"""

from __future__ import annotations

import threading
from collections import deque
from dataclasses import dataclass
from html import escape
from typing import Deque, Optional

# Severity-ish category. "block" = the gate refused; "ask" = it routed to a
# user confirmation; "info" = noteworthy but allowed.
#
# Every kind the gate records needs an entry here. An unknown kind renders
# as an unlabelled bullet, and the kinds nobody had added were the loudest
# ones: a refusal the agent tried to walk around, and every escape from a
# locked scope, appeared blanker in the panel than an ordinary denied
# command. ``tests/test_the_containment_panel_labels_what_it_records.py``
# fails when a new kind is recorded without a label.
_KINDS = {
    "deny_pattern":   ("⛔", "Denied command"),
    "secret_path":    ("🔑", "Secret-path access blocked"),
    "script_payload": ("📜", "Script payload blocked"),
    "self_mod":       ("🛡", "Self-modification guarded"),
    "calc_edit":      ("🧪", "Stored-calculation edit"),
    "read_only_write": ("📕", "Write into a read-only location"),
    "outside_ws":     ("📁", "Outside-workspace access"),
    "denied_by_user": ("✋", "User denied"),
    "denied_again":   ("🚫", "Refusal circumvented"),
    "denied_path_via_bash": ("🚫", "Refused path reached via shell"),
    "read_grant":     ("👁", "Directory opened for reading"),
    "approval_timeout": ("⌛", "Approval window expired"),
    "isolation":      ("🔒", "Filesystem isolation active"),
    "egress":         ("🌐", "Outbound data transfer"),
    "auto_verify":    ("🔁", "Auto-verify caught a problem"),
    "auto_verify_exhausted": ("🔁", "Auto-verify gave up"),
    "test_tamper":    ("🧪", "Failing test edited"),
    # Locked scope: the folder IS the promise, so each way out is named.
    "locked_scope_read":    ("🔒", "Locked scope: read outside the folder"),
    "locked_scope_bash":    ("🔒", "Locked scope: command outside the folder"),
    "locked_scope_symlink": ("🔒", "Locked scope: link out of the folder"),
    "locked_scope_exec":    ("🔒", "Locked scope: execution outside the folder"),
    "locked_scope_widen":   ("🔒", "Locked scope: attempt to widen it"),
    "locked_scope_parse":   ("🔒", "Locked scope: command could not be checked"),
    # MCP dispatch bypasses the native executor, so its refusals are their
    # own kinds.
    "plan_mode_mcp":     ("📝", "Plan mode: MCP tool refused"),
    "role_denied_mcp":   ("🎭", "Role: MCP tool refused"),
    "mcp_bash_no_perms": ("🔌", "MCP shell without permissions"),
    "mcp_write_no_perms": ("🔌", "MCP write without permissions"),
    "mcp_side_effect":   ("🔌", "MCP tool that is not read-only"),
    "mcp_side_effect_no_perms": ("🔌", "MCP side effect without permissions"),
    "mcp_module_name":   ("🔌", "MCP module name outside the adapters dir"),
    "mcp_budget":        ("🔌", "MCP tool surface over budget"),
}


# Kinds that stay in the panel even when blocked: the user's own refusal,
# an approval window that expired (the confirm broker already parks that
# one in the inbox itself), and the informational grants. Everything else
# is durable by default — a containment kind added later is loud until
# someone decides otherwise, which is the safer direction for this list.
_PANEL_ONLY_KINDS = frozenset({
    "denied_by_user", "approval_timeout", "read_grant", "isolation",
    "auto_verify", "auto_verify_exhausted",
})


@dataclass(frozen=True)
class SecurityEvent:
    seq: int
    kind: str
    tool: str
    detail: str
    blocked: bool


_LOCK = threading.Lock()
_EVENTS: Deque[SecurityEvent] = deque(maxlen=200)
_SEQ = 0
_ON_RECORD = None  # optional callback(SecurityEvent) for live UI refresh
# (kind, tool) pairs already raised in this process — the dedup key.
_ALERTED: set[tuple[str, str]] = set()


def set_listener(callback) -> None:
    """Register a single best-effort callback fired after each record()."""
    global _ON_RECORD
    _ON_RECORD = callback


def _should_alert(ev: SecurityEvent) -> bool:
    """First BLOCKED event of its (kind, tool) that is not panel-only."""
    if not ev.blocked or ev.kind in _PANEL_ONLY_KINDS:
        return False
    key = (ev.kind, ev.tool)
    with _LOCK:
        if key in _ALERTED:
            return False
        _ALERTED.add(key)
    return True


def _raise_alert(ev: SecurityEvent) -> None:
    """Durable inbox event for a containment block. Never raises.

    Outside the lock: this writes a file and may notify, and the gate
    that called ``record`` is holding the tool call open.
    """
    try:
        if not _should_alert(ev):
            return
        _, label = _KINDS.get(ev.kind, ("•", ev.kind))
        from .attention import emit_attention
        emit_attention(
            "security_alert",
            title=f"{label} ({ev.tool})" if ev.tool else label,
            detail=f"{ev.kind}: {ev.detail}"[:400],
        )
    except Exception:
        pass


def record(kind: str, tool: str, detail: str, *, blocked: bool = True) -> None:
    """Record one containment event. Never raises."""
    global _SEQ
    try:
        with _LOCK:
            _SEQ += 1
            ev = SecurityEvent(_SEQ, kind, tool, str(detail)[:300], blocked)
            _EVENTS.append(ev)
        cb = _ON_RECORD
        if cb is not None:
            try:
                cb(ev)
            except Exception:
                pass
        _raise_alert(ev)
    except Exception:
        pass


def recent(limit: int = 20) -> list[SecurityEvent]:
    with _LOCK:
        evs = list(_EVENTS)
    return evs[-limit:][::-1]            # newest first


def counts() -> dict:
    with _LOCK:
        evs = list(_EVENTS)
    return {"total": len(evs), "blocked": sum(1 for e in evs if e.blocked)}


def clear() -> None:
    with _LOCK:
        _EVENTS.clear()
        _ALERTED.clear()


def known_kinds() -> frozenset:
    """Event kinds this panel can label. Recorded kinds must be a subset."""
    return frozenset(_KINDS)


def format_panel_html(limit: int = 12) -> str:
    """Newest-first HTML feed for the dashboard containment panel."""
    evs = recent(limit)
    if not evs:
        return (
            "<div style='color:#888;font-size:12px;font-style:italic;'>"
            "🛡 Containment: no blocked or flagged actions this session.</div>"
        )
    n = counts()
    rows = []
    for e in evs:
        glyph, label = _KINDS.get(e.kind, ("•", e.kind))
        colour = "#d32f2f" if e.blocked else "#f59e0b"
        rows.append(
            f"<div style='font-family:monospace;font-size:12px;color:{colour};"
            f"margin:1px 0;'>{glyph} <b>{escape(label)}</b> "
            f"<span style='color:#888'>[{escape(e.tool)}]</span> "
            f"{escape(e.detail)[:160]}</div>"
        )
    head = (
        f"<div style='font-size:11px;color:#aaa;margin-bottom:3px;'>"
        f"🛡 Containment &nbsp; {n['blocked']} blocked / {n['total']} flagged"
        f"</div>"
    )
    return (
        "<div style='border-left:3px solid #d32f2f;padding:6px 10px;"
        "background:#d32f2f11;border-radius:4px;'>"
        + head
        + "<div style='max-height:160px;overflow-y:auto;'>"
        + "".join(rows)
        + "</div></div>"
    )


__all__ = [
    "SecurityEvent", "record", "recent", "counts", "clear",
    "known_kinds", "set_listener", "format_panel_html",
]
