"""Deterministic working-state block for compaction.

Compaction replaces history with a summary. The summary paths keep the
user's GOAL and the agent's own conclusions, but the working state --
files changed, last test outcome, open tasks, operator instructions,
denial reasons -- only survived by accident, inside whatever assistant
prose happened to mention it. This module rebuilds that state from the
session's own records on every compaction and hands it back as one
bounded block, verbatim into the new context.

Design rules:

* DETERMINISTIC — built from tool results and on-disk stores (change
  journal, task store), never from model text, never via a model call.
  Identical history and stores produce an identical block.
* BOUNDED — hard character ceiling; a recap that can itself blow the
  window it is meant to save is worse than none.
* NO SECRETS — every line is passed through ``scrub_secrets`` before it
  enters the block (the same redactor the transcript archive uses).
* NEWEST FIRST — a returning reader wants the latest state first.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterable

# Same prefixes the engine's _is_machine_turn uses, duplicated so this
# module never imports the engine (the engine imports this module).
_MACHINE_TURN_PREFIXES = (
    "[Command results]",
    "[Verify]",
    "[System]",
    "[Message from the session",
    "[scheduled",
    "[watch]",
)

# Operator/harness guidance is rare and load-bearing; keep a few, short.
_MAX_INSTRUCTIONS = 4
_MAX_INSTRUCTION_CHARS = 220
# Denials carry their reason on the error line; newest first.
_MAX_DENIALS = 4
_MAX_DENIAL_CHARS = 200
# Journal-backed file changes.
_MAX_FILES = 12
# Structured operator refusals (refusal_memory store), newest first.
_MAX_REFUSALS = 4
_MAX_REFUSAL_CHARS = 160
# Last test outcome per test file, newest wins.
_MAX_TESTS = 5
_MAX_TEST_LINE_CHARS = 160
# Open tasks come from the task store, capped by open_task_summary.
_MAX_TASKS = 6

# The whole block may never exceed this, whatever the history holds.
_MAX_BLOCK_CHARS = 2200

_HEADER = "[Working state — session recap, rebuilt at each compaction]"


def _message_text(msg: Any) -> str:
    content = msg.get("content", "") if isinstance(msg, dict) else ""
    if isinstance(content, str):
        return content
    try:
        return " ".join(
            str(b.get("text", "")) for b in content
            if isinstance(b, dict))
    except Exception:
        return ""


def _machine_texts(messages: Iterable[dict]) -> list[str]:
    """Text of every synthetic user turn, oldest first."""
    out: list[str] = []
    for msg in messages or []:
        if not isinstance(msg, dict) or msg.get("role") != "user":
            continue
        text = _message_text(msg)
        if text.lstrip().startswith(_MACHINE_TURN_PREFIXES):
            out.append(text)
    return out


def _clip(text: str, limit: int) -> str:
    text = " ".join(str(text).split())
    return text if len(text) <= limit else text[: limit - 3] + "..."


_INSTRUCTION_PREFIXES = ("[System]", "[Message from the session")


def _instructions(messages: Iterable[dict]) -> list[str]:
    """Operator / harness messages, newest first, deduped.

    Only turns that ARE guidance ([System], a message from another
    session): a keyword match on arbitrary machine output ("never",
    "operator") turned a line of test output into a standing order.
    """
    found: list[str] = []
    for msg in reversed(list(messages or [])):
        if not isinstance(msg, dict) or msg.get("role") != "user":
            continue
        text = _message_text(msg).strip()
        if not text.startswith(_INSTRUCTION_PREFIXES):
            continue
        short = _clip(text, _MAX_INSTRUCTION_CHARS)
        if short not in found:
            found.append(short)
        if len(found) >= _MAX_INSTRUCTIONS:
            break
    return found


_DENIAL_RE = re.compile(r"\b(?:denied|refused|not approved)\b", re.I)


def _denials(messages: Iterable[dict]) -> list[str]:
    """Newest-first denial lines (tool error / permission refusal), each
    with whatever reason the gate wrote next to it."""
    found: list[str] = []
    for text in reversed(_machine_texts(messages)):
        for line in text.splitlines():
            line = line.strip()
            # A refusal, not any error: "do not retry" on a plain failure
            # would stop the agent from fixing and re-running it.
            if not _DENIAL_RE.search(line):
                continue
            short = _clip(line, _MAX_DENIAL_CHARS)
            if short not in found:
                found.append(short)
            if len(found) >= _MAX_DENIALS:
                return found
    return found


# A pytest verdict line: "12 passed", "1 failed, 3 passed", "no tests ran".
_TEST_VERDICT_RE = re.compile(
    r"\b\d+ (?:passed|failed|errors?)\b|\bno tests ran\b", re.I)


def _test_outcomes(messages: Iterable[dict]) -> list[str]:
    """Last outcome line per test invocation, newest first.

    Test outcomes arrive as command results (``gate tests/... -> N
    passed``) or raw pytest output. The LAST matching line per run is
    the summary; keeping only those keeps the block short while every
    test file still shows its newest verdict.
    """
    lines: list[str] = []
    seen: set[str] = set()
    for text in reversed(_machine_texts(messages)):
        for line in text.splitlines():
            line = line.strip()
            if not _TEST_VERDICT_RE.search(line):
                continue
            short = _clip(line, _MAX_TEST_LINE_CHARS)
            if short not in seen:
                seen.add(short)
                lines.append(short)
            if len(lines) >= _MAX_TESTS:
                return lines
    return lines


_PATH_TOKEN_CHARS = set(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    "._-/")


def _touched_names(messages: Iterable[dict]) -> list[str]:
    """File (and file-adjacent) names from recent machine turns, newest
    first, deduped. This is the `where was I working` pointer; the change
    journal below says which of those files were actually modified."""
    names: list[str] = []
    for text in reversed(_machine_texts(messages)):
        for token in text.split():
            token = token.strip("`:;(),'\"")
            if not (token.endswith((".py", ".md", ".txt", ".json",
                                     ".yaml", ".yml", ".toml"))
                    or token.startswith(("delfin/", "tests/"))):
                continue
            if not set(token) <= _PATH_TOKEN_CHARS or len(token) < 6:
                continue
            if token not in names:
                names.append(token)
            if len(names) >= _MAX_FILES:
                return names
    return names


def _journal_files(session_id: str, limit: int = _MAX_FILES) -> list[str]:
    """Newest-first list of files this session changed, from the change
    journal (the store behind ``list_changes_made`` / undo)."""
    try:
        from . import change_journal as _cj
        records = _cj.list_changes(session_id) or []
    except Exception:
        return []
    out: list[str] = []
    for rec in reversed(records):
        path = str(rec.get("path", "") or "")
        if not path:
            continue
        if path not in out:
            out.append(path)
        if len(out) >= limit:
            break
    return out


def _open_tasks(workspace: Path) -> list[str]:
    """Open task-tool work for this workspace, via the existing summary
    surface (tri-state, capped, never raises)."""
    try:
        from . import agent_tasks as _at
        summary = _at.open_task_summary(Path(workspace), None, cap=_MAX_TASKS)
        if summary.get("state") != "open":
            return []
        lines: list[str] = []
        for status in ("in_progress", "pending", "blocked"):
            for t in summary.get(status, []) or []:
                subject = str(t.get("subject", "") or "").strip()
                if not subject:
                    continue
                why = str(t.get("blocked_reason", "") or "").strip()
                line = f"[{status}] #{t.get('id')} {subject}"
                if why:
                    line += f" — blocked: {_clip(why, 100)}"
                lines.append(_clip(line, 180))
        return lines[:_MAX_TASKS]
    except Exception:
        return []


def _refusals(workspace: Path) -> list[str]:
    """Structured operator refusals from the refusal_memory store, newest
    first, one line each. The message-based ``_denials`` above catches
    what the gate wrote this session; this store is what survives after
    compaction, recorded by the dialog itself."""
    try:
        from .refusal_memory import RefusalMemory
        mem = RefusalMemory.load(Path(workspace) / ".delfin" / "refusals.json")
    except Exception:
        return []
    lines: list[str] = []
    for r in reversed(mem.entries):
        what = r.target + ("/" if r.is_dir else "")
        line = f"{r.tool} {what} refused at {r.time or '?'}: {r.reason or '?'}"
        lines.append(_clip(line, _MAX_REFUSAL_CHARS))
        if len(lines) >= _MAX_REFUSALS:
            break
    return lines


def _scrub(text: str) -> str:
    """Remove credential material via the house redactor. Never raises."""
    try:
        from .output_guard import scrub_secrets
        return scrub_secrets(text)
    except Exception:
        return text


def build_working_state_block(
    messages: Iterable[dict],
    *,
    session_id: str = "",
    workspace: Path | None = None,
) -> str:
    """One bounded, redacted, deterministic working-state block.

    Returns ``""`` when there is nothing to say, so the caller adds no
    empty section to the context.
    """
    msgs = [m for m in (messages or []) if isinstance(m, dict)]

    tasks = _open_tasks(workspace) if workspace is not None else []
    files = _journal_files(session_id)
    names = _touched_names(msgs)
    tests = _test_outcomes(msgs)
    denials = _denials(msgs)
    instr = _instructions(msgs)
    refusals = _refusals(workspace) if workspace is not None else []

    if not (tasks or files or names or tests or denials or instr
            or refusals):
        return ""

    # Build by priority: tasks, instructions, denials, refusals, tests,
    # names, files. When the ceiling is hit, the LOWEST-priority
    # sections drop first; a hard cut with a marker guarantees the
    # ceiling regardless.
    pieces = [
        ("Open tasks (task tool):", tasks),
        ("Standing instructions:", instr),
        ("Recent denials (do not retry):", denials),
        ("Operator refusals (do not ask again):", refusals),
        ("Last test outcomes:", tests),
        ("Recently worked on:", names),
        ("Files changed this session (change journal):", files),
    ]
    chunks: list[str] = []
    used = len(_HEADER) + 60
    for header, data in pieces:
        if not data:
            continue
        chunk = f"{header}\n" + "\n".join(f"  {d}" for d in data)
        if used + len(chunk) + 2 > _MAX_BLOCK_CHARS:
            continue
        chunks.append(chunk)
        used += len(chunk) + 2
    block = _HEADER + "\n" + "\n\n".join(chunks) + "\n"
    if len(block) > _MAX_BLOCK_CHARS:
        block = block[:_MAX_BLOCK_CHARS] + "\n... [working state trimmed]\n"
    return _scrub(block)
