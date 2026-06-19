"""HTML / plaintext renderer for the live task ticker.

Implements TodoWrite display: a compact panel that
shows pending / in-progress / completed tasks with a colour-coded
glyph each. Designed to be re-rendered on every relevant tool
result (task_create / task_update) so the dashboard reflects
state immediately.

Usage::

    from delfin.agent.task_ticker import render_html, render_text

    html = render_html(workspace)   # for the dashboard widget
    txt  = render_text(workspace)   # for headless logs
"""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Iterable

from .agent_tasks import get_store


_GLYPHS = {
    "pending":     "&#9744;",      # ☐
    "in_progress": "&#9658;",      # ▶
    "completed":   "&#9745;",      # ☑
    "deleted":     "&#9747;",      # ☓
}
_TEXT_GLYPHS = {
    "pending":     "[ ]",
    "in_progress": "[>]",
    "completed":   "[x]",
    "deleted":     "[-]",
}
_COLOURS = {
    "pending":     "#888",
    "in_progress": "#0a84ff",
    "completed":   "#28a745",
    "deleted":     "#aa0000",
}
_ORDER = ("in_progress", "pending", "completed", "deleted")


def _sorted(tasks: Iterable[dict]) -> list[dict]:
    by_status: dict[str, list[dict]] = {s: [] for s in _ORDER}
    for t in tasks:
        s = t.get("status", "pending")
        by_status.setdefault(s, []).append(t)
    out: list[dict] = []
    for s in _ORDER:
        out.extend(sorted(by_status.get(s, []), key=lambda t: t.get("id", 0)))
    return out


def render_html(
    workspace: Path | str,
    *,
    session_id: str | None = None,
    show_completed: bool = True,
    max_rows: int = 30,
) -> str:
    """Return an HTML fragment listing tasks, ready for ipywidgets HTML."""
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False, session_id=session_id, with_seq=True)
    if not show_completed:
        raw = [t for t in raw if t.get("status") != "completed"]
    if not raw:
        return (
            "<div style='color:#888;font-size:12px;font-style:italic;'>"
            "No tasks yet — call task_create to start a plan.</div>"
        )
    rows: list[str] = []
    counts: dict[str, int] = {s: 0 for s in _ORDER}
    for t in _sorted(raw)[:max_rows]:
        status = str(t.get("status", "pending"))
        counts[status] = counts.get(status, 0) + 1
        glyph = _GLYPHS.get(status, _GLYPHS["pending"])
        colour = _COLOURS.get(status, "#444")
        subject = escape(str(t.get("subject", "")))
        active = escape(str(t.get("active_form", "")))
        label = active if status == "in_progress" and active else subject
        decorate = "text-decoration:line-through;" if status == "completed" else ""
        num = t.get("seq") if t.get("seq") is not None else t.get("id")
        rows.append(
            f"<div style='font-family:monospace;font-size:13px;"
            f"color:{colour};{decorate}'>{glyph} #{num} {label}</div>"
        )
    summary = (
        f"&#9658; {counts.get('in_progress', 0)} &nbsp; "
        f"&#9744; {counts.get('pending', 0)} &nbsp; "
        f"&#9745; {counts.get('completed', 0)}"
    )
    # Bounded, scrollable list: a long backlog stays a compact panel with a
    # side scrollbar instead of growing into a giant wall of tasks that pushes
    # the chat off-screen (Jerome 2026-06-13). The summary header stays fixed
    # above the scroll area so the counts are always visible.
    scroll_open = "<div style='max-height:220px;overflow-y:auto;'>"
    return (
        "<div style='border-left:3px solid #444;padding:6px 10px;"
        "background:#0001;border-radius:4px;'>"
        f"<div style='font-size:11px;color:#aaa;margin-bottom:4px;'>"
        f"Tasks &nbsp; {summary}</div>"
        + scroll_open
        + "".join(rows)
        + "</div>"
        + "</div>"
    )


def render_text(
    workspace: Path | str, *, session_id: str | None = None, max_rows: int = 30
) -> str:
    """Plain-text rendering for logs / headless contexts."""
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False, session_id=session_id, with_seq=True)
    if not raw:
        return "(no tasks)"
    lines: list[str] = []
    for t in _sorted(raw)[:max_rows]:
        status = str(t.get("status", "pending"))
        glyph = _TEXT_GLYPHS.get(status, "[?]")
        subject = str(t.get("subject", ""))
        num = t.get("seq") if t.get("seq") is not None else t.get("id")
        lines.append(f"{glyph} #{num} {subject}")
    return "\n".join(lines)


__all__ = ["render_html", "render_text"]
