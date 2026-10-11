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

from .agent_tasks import get_store, resolve_session_scope


_GLYPHS = {
    "pending":     "&#9744;",      # ☐
    "in_progress": "&#9658;",      # ▶
    "blocked":     "&#9940;",      # ⛔
    "completed":   "&#9745;",      # ☑
    "deleted":     "&#9747;",      # ☓
}
_TEXT_GLYPHS = {
    "pending":     "[ ]",
    "in_progress": "[>]",
    "blocked":     "[!]",
    "completed":   "[x]",
    "deleted":     "[-]",
}
_COLOURS = {
    "pending":     "#888",
    "in_progress": "#0a84ff",
    "blocked":     "#d97706",
    "completed":   "#28a745",
    "deleted":     "#aa0000",
}
# Blocked sits directly under in_progress: it is open work the user may
# have to unblock, so it must not be filed away below the done pile.
_ORDER = ("in_progress", "blocked", "pending", "completed", "deleted")


def _sorted(tasks: Iterable[dict]) -> list[dict]:
    by_status: dict[str, list[dict]] = {s: [] for s in _ORDER}
    for t in tasks:
        s = t.get("status", "pending")
        by_status.setdefault(s, []).append(t)
    out: list[dict] = []
    for s in _ORDER:
        out.extend(sorted(by_status.get(s, []), key=lambda t: t.get("id", 0)))
    return out


#: A block set before the session's last message may no longer hold: the
#: person has since said something, so the "waiting on" note may be about a
#: reason that no longer applies. Marked so on the panel rather than trusted
#: silently (U7, phase 2).
_STALE_MARKER = " (stale)"


def _block_epoch(task: dict) -> float | None:
    """Epoch seconds of when ``task`` was last updated (i.e. when it was
    blocked), or None when the timestamp is missing or unparseable. Never
    raises: a ticker is a courtesy that must not take the panel down over a
    malformed timestamp."""
    ts = str(task.get("updated_at") or task.get("created_at") or "").strip()
    if not ts:
        return None
    try:
        import datetime as _dt
        naive = _dt.datetime.fromisoformat(ts.rstrip("Z"))
    except ValueError:
        return None
    if naive.tzinfo is None:
        naive = naive.replace(tzinfo=_dt.timezone.utc)
    return naive.timestamp()


def _stale_suffix(task: dict, last_message_at: float | None) -> str:
    """``_STALE_MARKER`` when *task* is a blocked task whose reason was set
    strictly before ``last_message_at`` — i.e. the block outdates the last
    message. "" otherwise. A blocked task with no reason cannot reach the
    store (an invariant), and one whose timestamp is missing or not older
    is current, so neither is flagged.

    Never raises: a reference that cannot parse as a number is treated as
    "nothing stale" (same as ``None``), so a malformed value arriving over a
    wire cannot take the panel down with it (reviewer nacht-s12).
    """
    if last_message_at is None:
        return ""
    try:
        ref = float(last_message_at)
    except (TypeError, ValueError):
        return ""
    if str(task.get("status")) != "blocked":
        return ""
    if not str(task.get("blocked_reason") or "").strip():
        return ""
    block = _block_epoch(task)
    if block is None or not (block < ref):
        return ""
    return _STALE_MARKER


def session_last_activity_at(
    workspace: Path | str, *, session_id: str | None = None,
) -> float | None:
    """Epoch seconds of the session's last recorded activity — the proxy for
    "the last message" a stale blocked reason is judged against — or None
    when there is no such session or its timestamp cannot be read.

    The reference a call site passes as ``last_message_at`` to
    :func:`rows` / :func:`render_text` / :func:`render_html`. Sessions are
    stored globally keyed by session id (``session_store``), so the lookup
    is by ``session_id``; ``workspace`` is accepted for signature symmetry
    with the render functions but does not drive the lookup.

    Never raises: a session record that is absent, not a dict, or carries a
    missing / non-numeric ``updated_at`` yields None (nothing stale), the
    same guarantee the render functions already keep (a broken reference
    must not take the panel down).
    """
    sid = str(session_id or "").strip()
    if not sid:
        return None
    try:
        from . import session_store
        rec = session_store.load_session(sid)
    except Exception:
        return None
    if not isinstance(rec, dict):
        return None
    try:
        value = float(rec.get("updated_at") or 0)
    except (TypeError, ValueError):
        return None
    return value or None


def render_html(
    workspace: Path | str,
    *,
    session_id: str | None = None,
    show_completed: bool = True,
    max_rows: int = 30,
    bare: bool = False,
    last_message_at: float | None = None,
) -> str:
    """Return an HTML fragment listing tasks, ready for ipywidgets HTML.

    ``bare=True`` returns the rows alone, without the framed header: the
    body of a fold whose title already carries the counts.

    ``last_message_at`` marks a blocked reason set before it with a stale
    suffix in a dimmer colour (None marks nothing).

    The session filter goes through the shared resolver: an empty id is
    UNSCOPED, the same as everywhere else. Reading it as "no session,
    therefore no tasks" is what made this panel print "No tasks yet"
    while the model's own reminder listed the workspace's open work.
    """
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False,
                     session_id=resolve_session_scope(session_id),
                     with_seq=True)
    if not show_completed:
        raw = [t for t in raw if t.get("status") != "completed"]
    if not raw:
        # Nothing: a panel announcing that there is nothing took a line
        # from the chat on every session that had not planned anything.
        return ""
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
        if status == "blocked":
            reason = escape(str(t.get("blocked_reason", ""))[:60])
            label += f" — waiting on {reason}" if reason else ""
            stale = _stale_suffix(t, last_message_at)
            if stale:
                # A block that outdates the last message is a stale note, not
                # a fact: dimmer than the reason itself, so it reads as a
                # caution rather than as "waiting" all the same.
                label += (f" <span style='color:#999'>{escape(stale)}</span>"
                          .replace("&", "&amp;"))
        elif status == "completed" and t.get("verified") == "unmet":
            # A completion the evidence did not support stays visible as
            # such; the alternative is a green tick over an open hole.
            label += " (unverified)"
        decorate = "text-decoration:line-through;" if status == "completed" else ""
        num = t.get("seq") if t.get("seq") is not None else t.get("id")
        rows.append(
            f"<div style='font-family:monospace;font-size:13px;"
            f"color:{colour};{decorate}'>{glyph} #{num} {label}</div>"
        )
    summary = (
        f"&#9658; {counts.get('in_progress', 0)} &nbsp; "
        f"&#9744; {counts.get('pending', 0)} &nbsp; "
        + (f"&#9940; {counts.get('blocked', 0)} &nbsp; "
           if counts.get("blocked", 0) else "")
        + f"&#9745; {counts.get('completed', 0)}"
    )
    # Bounded, scrollable list: a long backlog stays a compact panel with a
    # side scrollbar instead of growing into a giant wall of tasks that pushes
    # the chat off-screen (Jerome 2026-06-13). The summary header stays fixed
    # above the scroll area so the counts are always visible.
    scroll_open = "<div style='max-height:220px;overflow-y:auto;'>"
    if bare:
        return scroll_open + "".join(rows) + "</div>"
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


def counts(
    workspace: Path | str,
    *,
    session_id: str | None = None,
) -> tuple[int, int]:
    """(open, completed) tasks of the session, for one-line summaries."""
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False,
                     session_id=resolve_session_scope(session_id))
    done = sum(1 for t in raw if t.get("status") == "completed")
    return len(raw) - done, done


_PLAIN_GLYPHS = {
    "pending": "\u2610", "in_progress": "\u25b6", "blocked": "\u26d4",
    "completed": "\u2611", "deleted": "\u2613",
}


def rows(
    workspace: Path | str,
    *,
    session_id: str | None = None,
    max_rows: int = 30,
    last_message_at: float | None = None,
) -> list[dict]:
    """The task list as data, in display order, for a surface that builds
    its own rows (the dashboard's clickable task list).

    Each row: ``status``, ``glyph`` (plain text), ``num`` (the session
    number the user reads), ``label`` (what the row shows: the active form
    while in progress, the blocking reason appended) and ``subject`` (what
    a click puts into the message box).

    ``last_message_at`` is the epoch-second time of the session's last
    message; a blocked task whose reason was set before it reads as stale
    (see :func:`_stale_suffix`). None marks nothing.
    """
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False,
                     session_id=resolve_session_scope(session_id),
                     with_seq=True)
    out: list[dict] = []
    for t in _sorted(raw)[:max_rows]:
        status = str(t.get("status", "pending"))
        subject = str(t.get("subject", ""))
        active = str(t.get("active_form", ""))
        label = active if status == "in_progress" and active else subject
        if status == "blocked":
            reason = str(t.get("blocked_reason", ""))[:60]
            label += f" \u2014 waiting on {reason}" if reason else ""
            label += _stale_suffix(t, last_message_at)
        elif status == "completed" and t.get("verified") == "unmet":
            label += " (unverified)"
        out.append({
            "status": status,
            "glyph": _PLAIN_GLYPHS.get(status, _PLAIN_GLYPHS["pending"]),
            "num": t.get("seq") if t.get("seq") is not None else t.get("id"),
            "label": label,
            "subject": subject,
        })
    return out


def rows_from_todos(todos: Iterable[dict]) -> list[dict]:
    """The same rows as :func:`rows`, from a TodoWrite payload.

    A CLI backend (the Anthropic CLI) keeps its plan in its own TodoWrite
    list, not in this task store; the dashboard's task line reads it
    through here so the line shows the plan whichever backend runs.
    Entries: ``content``, ``status``, ``activeForm``.
    """
    order = {s: i for i, s in enumerate(_ORDER)}
    items = []
    for n, t in enumerate(todos or (), 1):
        if not isinstance(t, dict):
            continue
        status = str(t.get("status") or "pending")
        if status not in _PLAIN_GLYPHS:
            status = "pending"
        subject = " ".join(str(t.get("content") or "").split())
        if not subject:
            continue
        active = " ".join(str(t.get("activeForm") or "").split())
        items.append((order.get(status, 0), n, {
            "status": status,
            "glyph": _PLAIN_GLYPHS[status],
            "num": n,
            "label": active if status == "in_progress" and active else subject,
            "subject": subject,
        }))
    return [r for _, _, r in sorted(items, key=lambda x: (x[0], x[1]))]


def title_for(rows_: list[dict]) -> str:
    """Plain-text title of the task fold for *rows_*: the counts and the
    task in hand, e.g. ``Tasks  ▶ 1  ☐ 2  ☑ 0 · Parsing``. Empty when no
    row is open, so the fold costs no height unless it has something to
    say."""
    counts = {s: 0 for s in _ORDER}
    current = ""
    for r in rows_ or ():
        counts[r["status"]] = counts.get(r["status"], 0) + 1
        if r["status"] == "in_progress" and not current:
            current = r["label"]
    open_count = counts["in_progress"] + counts["pending"] + counts["blocked"]
    if not open_count:
        return ""
    parts = [f"\u25b6 {counts['in_progress']}", f"\u2610 {counts['pending']}"]
    if counts["blocked"]:
        parts.append(f"\u26d4 {counts['blocked']}")
    parts.append(f"\u2611 {counts['completed']}")
    title = "Tasks  " + "  ".join(parts)
    if current:
        title += " \u00b7 " + current[:90]
    return title


def render_title(
    workspace: Path | str,
    *,
    session_id: str | None = None,
) -> str:
    """:func:`title_for` over the session's tasks in the store. Plain text
    because a fold title is not HTML."""
    return title_for(rows(workspace, session_id=session_id, max_rows=10_000))


def next_steps(
    workspace: Path | str, *, session_id: str | None = None,
    limit: int = 3,
) -> list[str]:
    """The open work, shortest first, as things a user could ask for next.

    Not invented: these are the tasks the agent itself left open, so a
    suggestion is never a guess about what the user wants and never
    costs a token to produce. Pending before blocked — a blocked item is
    waiting on something, and offering it as the next step would be an
    invitation into a wall.

    Never raises: a suggestion is a courtesy, and one that takes the
    prompt with it is not.
    """
    try:
        store = get_store(Path(workspace))
        raw = store.list(include_deleted=False,
                         session_id=resolve_session_scope(session_id),
                         with_seq=True)
    except Exception:
        return []
    out: list[str] = []
    for want in ("pending", "in_progress", "blocked"):
        for t in _sorted(raw or []):
            if str(t.get("status", "pending")) != want:
                continue
            subject = str(t.get("subject", "")).strip()
            if subject and subject not in out:
                out.append(subject)
            if len(out) >= max(0, int(limit)):
                return out
    return out


def render_text(
    workspace: Path | str, *, session_id: str | None = None, max_rows: int = 30,
    last_message_at: float | None = None,
) -> str:
    """Plain-text rendering for logs / headless contexts."""
    store = get_store(Path(workspace))
    raw = store.list(include_deleted=False,
                     session_id=resolve_session_scope(session_id),
                     with_seq=True)
    if not raw:
        return "(no tasks)"
    lines: list[str] = []
    for t in _sorted(raw)[:max_rows]:
        status = str(t.get("status", "pending"))
        glyph = _TEXT_GLYPHS.get(status, "[?]")
        subject = str(t.get("subject", ""))
        if status == "blocked" and t.get("blocked_reason"):
            subject += f" — waiting on {str(t['blocked_reason'])[:60]}"
            subject += _stale_suffix(t, last_message_at)
        num = t.get("seq") if t.get("seq") is not None else t.get("id")
        lines.append(f"{glyph} #{num} {subject}")
    return "\n".join(lines)


__all__ = ["render_html", "render_text", "next_steps"]
