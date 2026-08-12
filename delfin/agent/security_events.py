"""Visible containment surface: a thread-safe record of every safety decision
the agent's permission gate makes (a denied command, a blocked secret/script
payload, a self-modification attempt, an out-of-workspace access).

Blocks are otherwise only visible as a tool error to the model and a line in
the 0600 audit log. This module keeps a small in-process ring buffer so the
dashboard can show a prominent "🛡 Containment" panel — so an attempt to break
out is immediately VISIBLE to the user, not buried.

Pure data + rendering; no I/O on the hot path beyond a bounded deque append.
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
    "mcp_budget":        ("🔌", "MCP tool surface over budget"),
}


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


def set_listener(callback) -> None:
    """Register a single best-effort callback fired after each record()."""
    global _ON_RECORD
    _ON_RECORD = callback


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
