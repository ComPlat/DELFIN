"""Episodic session memory for the DELFIN agent.

Sessions and handoff briefs used to be write-only: a finished session
left a full state dump under ``~/.delfin/agent_sessions/`` but nothing
ever recalled it, so a future session could not answer "have I worked on
something like this before?". This module closes that loop with one
COMPACT episode record per session, stored inside the per-project memory
namespace and recalled BM25-ranked against the current task text.

Design constraints (mirroring the typed memory store):
- LLM-free and deterministic: fields are derived extractively from the
  saved session state, never generated.
- Bounded: one small markdown file per session, newest 100 kept.
- Best-effort: every caller wraps the hook in try/except; nothing here
  may break a session save.

Files live under ``<project memory dir>/episodes/<date>_<sid8>.md`` with
scalar frontmatter (session_id, date, verdict, cost) and four body
sections (Goal / Outcome / Decisions / Open items).
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Any

from .memory_store import (
    _atomic_write,
    _delfin_memory_dir,
    _parse_frontmatter,
    bm25_scores,
)

# Newest episodes kept after every save; older ones are deleted.
_MAX_EPISODES = 100

# Compactness caps: an episode is a recall hook, not a transcript.
_GOAL_CAP = 500
_OUTCOME_CAP = 500
_DECISION_CAP = 300
_OPEN_ITEM_CAP = 200
_MAX_DECISIONS = 5
_MAX_OPEN_ITEMS = 8


def _episodes_dir(repo_root: Path | str) -> Path:
    """Episode store directory inside the per-project memory namespace."""
    return _delfin_memory_dir(Path(repo_root)) / "episodes"


def _one_line(text: str, cap: int) -> str:
    """Collapse whitespace to a single line and cap the length."""
    return " ".join(str(text or "").split())[:cap].strip()


def _safe_sid(session_id: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]", "", session_id or "")[:8] or "session"


def save_episode(
    session_id: str,
    *,
    repo_root: Path | str,
    goal: str,
    outcome: str,
    decisions: list[str],
    open_items: list[str],
    cost_usd: float = 0.0,
    verdict: str = "",
) -> Path:
    """Write one compact episode record for a finished (or auto-saved)
    session to ``<project memory dir>/episodes/<date>_<sid8>.md``.

    Saving the same session again on the same day overwrites the file in
    place, so repeated auto-saves refresh one record instead of piling
    up. After the write the store is pruned to the newest
    ``_MAX_EPISODES`` entries. Returns the written path.
    """
    d = _episodes_dir(repo_root)
    d.mkdir(parents=True, exist_ok=True)

    date = time.strftime("%Y-%m-%d")
    sid8 = _safe_sid(session_id)
    fpath = d / f"{date}_{sid8}.md"

    goal_line = _one_line(goal, _GOAL_CAP) or "(no goal recorded)"
    outcome_line = _one_line(outcome, _OUTCOME_CAP) or "(no outcome recorded)"
    decision_lines = [
        _one_line(item, _DECISION_CAP)
        for item in (decisions or [])[:_MAX_DECISIONS]
        if _one_line(item, _DECISION_CAP)
    ]
    open_lines = [
        _one_line(item, _OPEN_ITEM_CAP)
        for item in (open_items or [])[:_MAX_OPEN_ITEMS]
        if _one_line(item, _OPEN_ITEM_CAP)
    ]

    try:
        cost = float(cost_usd or 0.0)
    except (TypeError, ValueError):
        cost = 0.0

    lines = [
        "---",
        f"session_id: {_one_line(session_id, 80)}",
        f"date: {date}",
        f"verdict: {_one_line(verdict, 40)}",
        f"cost: {cost:.4f}",
        "---",
        "",
        "## Goal",
        goal_line,
        "",
        "## Outcome",
        outcome_line,
        "",
        "## Decisions",
    ]
    if decision_lines:
        lines += [f"- {item}" for item in decision_lines]
    else:
        lines.append("- (none recorded)")
    lines += ["", "## Open items"]
    if open_lines:
        lines += [f"- {item}" for item in open_lines]
    else:
        lines.append("- (none)")
    lines.append("")

    _atomic_write(fpath, "\n".join(lines))
    _prune_episodes(repo_root, keep=_MAX_EPISODES)
    return fpath


def _message_text(content: Any) -> str:
    """Extract plain text from a chat/engine message content field.

    Engine messages may carry a list of content blocks (text blocks,
    tool_use, tool_result); only ``text`` blocks contribute.
    """
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts: list[str] = []
        for block in content:
            if isinstance(block, dict) and block.get("type") == "text":
                parts.append(str(block.get("text", "")))
            elif isinstance(block, str):
                parts.append(block)
        return "\n".join(parts)
    return ""


def build_episode_from_state(
    engine_state: dict, chat_messages: list
) -> dict:
    """Derive episode fields deterministically from saved session state.

    - goal: first user message (chat first; engine messages as fallback
      when the chat list is empty, e.g. headless CLI saves)
    - outcome: tail (last paragraph) of the last assistant message
    - decisions: first paragraph of each of the last 3 assistant messages
    - open_items: pending / in-progress / blocked entries from the
      persisted task list (``todo_payload``), falling back to the task
      store the state points at
    - cost_usd: taken straight from the state

    Extractive and LLM-free by design; mirrors the handoff-brief
    derivation in ``session_store.build_handoff_brief``. Returns a dict
    matching ``save_episode``'s keyword parameters (minus verdict).
    """
    state = engine_state or {}
    messages: list[dict] = [
        m for m in (chat_messages or []) if isinstance(m, dict)
    ]
    if not messages:
        messages = [
            m for m in (state.get("engine_messages") or [])
            if isinstance(m, dict)
        ]

    goal = ""
    for m in messages:
        if m.get("role") == "user":
            goal = _message_text(m.get("content", "")).strip()
            if goal:
                break

    assistant_texts = [
        _message_text(m.get("content", "")).strip()
        for m in messages
        if m.get("role") == "assistant"
        and _message_text(m.get("content", "")).strip()
    ]

    outcome = ""
    if assistant_texts:
        outcome = assistant_texts[-1].split("\n\n")[-1].strip()

    decisions: list[str] = []
    for text in assistant_texts[-3:]:
        first_para = text.split("\n\n")[0].strip()
        if first_para:
            decisions.append(first_para)

    # The task list, from the state when a caller put it there and from
    # the real store otherwise. The engine's exported state carries no
    # task field, so on the headless path this was always empty and
    # every episode recorded "no open items" for work still in progress.
    todos = state.get("todo_payload")
    if not todos:
        try:
            from .session_store import tasks_as_todo_payload
            todos = tasks_as_todo_payload(
                state.get("workspace") or state.get("project_dir") or "",
                state.get("session_id") or "")
        except Exception:
            todos = []
    open_items: list[str] = []
    for t in todos or []:
        if not isinstance(t, dict):
            continue
        status = t.get("status", "")
        if status in ("pending", "in_progress", "blocked"):
            reason = str(t.get("blocked_reason", "") or "")
            open_items.append(
                f"#{t.get('id', '?')} {t.get('subject', '')} ({status}"
                + (f": {reason}" if status == "blocked" and reason else "")
                + ")"
            )

    try:
        cost = float(state.get("cost_usd", 0.0) or 0.0)
    except (TypeError, ValueError):
        cost = 0.0

    return {
        "goal": goal,
        "outcome": outcome,
        "decisions": decisions,
        "open_items": open_items,
        "cost_usd": cost,
    }


def _extract_section(body: str, header: str) -> str:
    m = re.search(
        rf"^## {re.escape(header)}\s*\n(.*?)(?=^## |\Z)",
        body or "", re.M | re.S,
    )
    return m.group(1).strip() if m else ""


def _section_items(section_text: str) -> list[str]:
    items = [
        line[2:].strip()
        for line in (section_text or "").splitlines()
        if line.startswith("- ")
    ]
    return [i for i in items if i and not i.startswith("(none")]


def list_episodes(repo_root: Path | str) -> list[dict]:
    """One record per stored episode, newest first.

    Ordering key is (date frontmatter, file mtime, filename) descending
    so same-day episodes still order deterministically. Unreadable files
    are skipped. Each record carries: file, path, session_id, date,
    verdict, cost, goal, outcome, decisions, open_items, mtime.
    """
    d = _episodes_dir(repo_root)
    if not d.is_dir():
        return []
    out: list[dict] = []
    for p in d.glob("*.md"):
        try:
            text = p.read_text(encoding="utf-8")
        except OSError:
            continue
        meta, body = _parse_frontmatter(text)
        try:
            mtime = int(p.stat().st_mtime)
        except OSError:
            mtime = 0
        try:
            cost = float(meta.get("cost", "0") or 0.0)
        except (TypeError, ValueError):
            cost = 0.0
        out.append({
            "file": p.name,
            "path": str(p),
            "session_id": meta.get("session_id", ""),
            "date": meta.get("date", ""),
            "verdict": meta.get("verdict", ""),
            "cost": cost,
            "goal": _extract_section(body, "Goal"),
            "outcome": _extract_section(body, "Outcome"),
            "decisions": _section_items(_extract_section(body, "Decisions")),
            "open_items": _section_items(
                _extract_section(body, "Open items")),
            "mtime": mtime,
        })
    out.sort(
        key=lambda r: (r["date"], r["mtime"], r["file"]), reverse=True)
    return out


def _prune_episodes(repo_root: Path | str, *, keep: int = _MAX_EPISODES) -> list[str]:
    """Delete all but the newest ``keep`` episodes; return deleted names."""
    deleted: list[str] = []
    for rec in list_episodes(repo_root)[max(keep, 0):]:
        try:
            Path(rec["path"]).unlink()
            deleted.append(rec["file"])
        except OSError:
            continue
    return deleted


def _placeholder(text: str) -> bool:
    return not text or text.startswith("(no ")


def recall_episodes(
    repo_root: Path | str,
    task_text: str,
    *,
    max_entries: int = 2,
    max_chars: int = 1200,
) -> str:
    """Compact prompt block of past sessions similar to ``task_text``.

    Episodes are BM25-ranked (goal + outcome + open items as the
    document) against the task text; only positive-scoring entries are
    returned, newest first on ties. Empty store, empty task text or no
    match all yield ``""`` so the prompt only grows when there is a real
    hook. Output is capped at ``max_entries`` lines / ``max_chars``.
    """
    episodes = list_episodes(repo_root)
    if not episodes or not (task_text or "").strip():
        return ""

    docs = [
        "\n".join([
            ep["goal"], ep["outcome"], " ".join(ep["open_items"]),
        ])
        for ep in episodes
    ]
    scores = bm25_scores(task_text, docs)
    # Episodes arrive newest-first, so index order breaks score ties in
    # favour of the most recent session.
    ranked = sorted(
        (i for i, s in enumerate(scores) if s > 0),
        key=lambda i: (-scores[i], i),
    )
    if not ranked:
        return ""

    header = "# Similar past sessions"
    lines = [header]
    used = len(header)
    for i in ranked[:max(max_entries, 0)]:
        ep = episodes[i]
        goal = _one_line(ep["goal"], 110) or "(no goal recorded)"
        outcome = _one_line(ep["outcome"], 150) or "(no outcome recorded)"
        verdict_part = f" ({ep['verdict']})" if ep["verdict"] else ""
        line = f"- {ep['date']}{verdict_part}: {goal} -> {outcome}"
        opens = [o for o in ep["open_items"] if not _placeholder(o)]
        if opens:
            line += f" [open: {_one_line('; '.join(opens[:3]), 120)}]"
        if used + 1 + len(line) > max_chars:
            remaining = max_chars - used - 1
            if len(lines) == 1 and remaining > 40:
                lines.append(line[:remaining])
            break
        lines.append(line)
        used += 1 + len(line)

    if len(lines) == 1:
        return ""
    return "\n".join(lines)
