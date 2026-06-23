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
_KINDS = {
    "deny_pattern":   ("⛔", "Denied command"),
    "secret_path":    ("🔑", "Secret-path access blocked"),
    "script_payload": ("📜", "Script payload blocked"),
    "self_mod":       ("🛡", "Self-modification guarded"),
    "outside_ws":     ("📁", "Outside-workspace access"),
    "denied_by_user": ("✋", "User denied"),
    "isolation":      ("🔒", "Filesystem isolation active"),
    "egress":         ("🌐", "Outbound data transfer"),
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
    "set_listener", "format_panel_html",
]
