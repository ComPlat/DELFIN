"""Persistent session storage for the DELFIN Agent.

Saves and loads agent sessions so conversations survive dashboard restarts.
Each session is a JSON file in ``~/.delfin/agent_sessions/`` (per-user,
per-machine, 0600 permissions).  Sessions MUST NOT be committed to the repo.
"""

from __future__ import annotations

import json
import os
import re
import time
import zipfile
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

# The version stamped into every session file this build writes. The files
# had carried twenty-two keys and no version at all, so a reader had no way
# to tell an old file missing half its fields from a complete one, and a
# file from a newer build was read with today's loaders — silently dropping
# whatever the newer build had added.
#
#   v1  unversioned: no run clock, no outcome-cost baseline
#   v2  run_elapsed_s + last_outcome_cost carried across a resume
SESSION_SCHEMA_VERSION = 2


class SessionSchemaError(Exception):
    """A session file this build must not read.

    Raised for a schema version NEWER than ``SESSION_SCHEMA_VERSION``.
    Refusing is the point: loading it anyway would drop the fields this
    build does not know about and then present the result as a complete
    restore.
    """


def migrate_session_data(
    data: dict[str, Any], from_version: int,
) -> tuple[dict[str, Any], tuple[str, ...]]:
    """Bring an older session dict up to the current schema.

    Returns the upgraded copy (the input is never mutated) and the names of
    the migrations that ran, so the restore report can say which ones did.
    """
    out = dict(data)
    applied: list[str] = []
    if from_version < 2:
        # v1 files predate the run clock and the outcome-cost baseline.
        # Neither can be reconstructed: the elapsed time of the earlier run
        # was never recorded. The resumed run therefore starts its
        # wall-clock budget from zero — stated here rather than left to be
        # discovered as a budget that never expires.
        out.setdefault("run_elapsed_s", 0.0)
        out.setdefault("last_outcome_cost", float(out.get("cost_usd") or 0.0))
        applied.append("v1_to_v2_run_clock_and_cost_baseline")
    out["schema_version"] = SESSION_SCHEMA_VERSION
    return out, tuple(applied)


def _is_json_safe(value: Any) -> bool:
    try:
        json.dumps(value)
    except (TypeError, ValueError):
        return False
    return True


def _chmod_user_only(path: Path) -> None:
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def tasks_as_todo_payload(workspace: Any, session_id: Any) -> list[dict]:
    """The real task store, in the shape ``todo_payload`` consumers read.

    ``todo_payload`` had exactly one producer: the dashboard branch that
    handles an external CLI's todo tool. On every other path -- the
    open-weights backend, the headless CLI -- it stayed empty, so each
    saved episode and every handoff brief told a fresh agent that
    nothing was outstanding while the store still held in_progress
    work. The task store is the list; this is how it reaches the brief.

    Returns [] on any failure — a save must never break on bookkeeping.
    """
    if not workspace:
        return []
    try:
        from .agent_tasks import get_store, resolve_session_scope
        tasks = get_store(Path(workspace)).list(
            session_id=resolve_session_scope(session_id), with_seq=True)
    except Exception:
        return []
    out: list[dict] = []
    for t in tasks or []:
        try:
            out.append({
                "id": t.get("id"),
                "seq": t.get("seq"),
                "subject": str(t.get("subject", "")),
                "status": str(t.get("status", "")),
                "blocked_reason": str(t.get("blocked_reason", "") or ""),
            })
        except Exception:
            continue
    return out


def _atomic_write_text(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` atomically (temp file + ``os.replace``).

    Same pattern as ``memory_store._atomic_write``: a plain ``write_text``
    truncates the file first, so a crash mid-write (or a reader racing the
    writer — two dashboards can share one sessions dir) leaves a torn or
    empty session file. With replace, the file on disk is always either
    the complete old version or the complete new one.

    The transcript archive (``*.jsonl``) intentionally stays append-mode:
    appends can't be replaced atomically without rewriting the whole file,
    and its reader already skips a torn last line.
    """
    import tempfile
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
        _chmod_user_only(path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass

# Legacy mode name migration
# Kept identical to the engine's table on purpose: when the two disagreed,
# the same saved session migrated differently depending on which one read
# it -- the engine knew `pipeline`, this one did not. See the engine's copy
# for why the multi-role pipelines were retired.
_LEGACY_MODE_MAP = {
    "default": "solo",
    "high_risk": "solo",
    "runtime_cluster": "solo",
    "release": "solo",
    "quick": "solo",
    "reviewed": "solo",
    "tdd": "solo",
    "cluster": "solo",
    "full": "solo",
    "pipeline": "solo",
}


def _migrate_mode(mode: str) -> str:
    return _LEGACY_MODE_MAP.get(mode, mode)


_SESSIONS_DIR = Path.home() / ".delfin" / "agent_sessions"


def _ensure_dir() -> Path:
    _SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
    return _SESSIONS_DIR


# ---------------------------------------------------------------------------
# Single-writer guard — one live process per session id
# ---------------------------------------------------------------------------
#
# Two dashboards (or a dashboard + CLI) resuming the SAME session id used
# to silently clobber each other's saves turn-by-turn: last writer wins,
# the other's turns vanish on its next reload. A tiny per-session lock
# file (pid + timestamp) makes the second writer's save fail loudly
# instead of corrupting the first writer's history.

_LOCK_MAX_AGE_S = 3600.0   # locks older than 1h are stale regardless of pid


class SessionLockedError(RuntimeError):
    """Another live process holds the writer lock for this session.

    Raised by :func:`save_session` (via :func:`acquire_session_lock`)
    instead of overwriting a session another dashboard/CLI is actively
    writing. Callers that wrap saves in try/except surface it as a
    save-failed warning; the session file itself is left untouched.
    """

    def __init__(self, session_id: str, holder_pid: int) -> None:
        super().__init__(
            f"session '{session_id}' is locked by live pid {holder_pid}; "
            f"save skipped so two writers don't clobber each other")
        self.session_id = session_id
        self.holder_pid = holder_pid


def _lock_path(session_id: str) -> Path:
    return _SESSIONS_DIR / f"{session_id}.lock"


def _pid_alive(pid: int) -> bool:
    """True when ``pid`` is a live process (EPERM counts as alive)."""
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True   # exists, owned by someone else
    except OSError:
        return False
    return True


def acquire_session_lock(session_id: str) -> Path:
    """Acquire or refresh the per-session single-writer lock.

    The lock file holds ``{"pid": ..., "ts": ...}``. Our own pid just
    refreshes the timestamp. A DIFFERENT pid that is still alive and
    whose lock is younger than :data:`_LOCK_MAX_AGE_S` raises
    :class:`SessionLockedError`. Stale locks (dead pid, or older than
    1h) are broken silently.
    """
    d = _ensure_dir()
    p = d / f"{session_id}.lock"
    me = os.getpid()
    holder, ts = 0, 0.0
    try:
        info = json.loads(p.read_text(encoding="utf-8"))
        holder = int(info.get("pid", 0) or 0)
        ts = float(info.get("ts", 0) or 0)
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        pass
    if holder and holder != me:
        if (time.time() - ts) < _LOCK_MAX_AGE_S and _pid_alive(holder):
            raise SessionLockedError(session_id, holder)
        # Stale lock (dead pid or >1h old): break silently.
    _atomic_write_text(p, json.dumps({"pid": me, "ts": time.time()}))
    return p


def release_session_lock(session_id: str) -> None:
    """Drop the lock if THIS process holds it. Best-effort, never raises."""
    p = _lock_path(session_id)
    try:
        info = json.loads(p.read_text(encoding="utf-8"))
        if int(info.get("pid", 0) or 0) == os.getpid():
            p.unlink()
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        pass


# ---------------------------------------------------------------------------
# Mid-turn checkpoints — crash insurance inside a single long turn
# ---------------------------------------------------------------------------
#
# Full session saves happen only at turn boundaries, but one turn can run
# hundreds of tool rounds over hours. The engine writes this lightweight
# checkpoint (throttled, best-effort) during the tool loop and clears it
# at normal turn end — so a ``<sid>.turn.json`` left on disk is proof of
# an unclean process death, and feeds a recovery note on the next load.


def _turn_checkpoint_path(session_id: str) -> Path:
    return _SESSIONS_DIR / f"{session_id}.turn.json"


def save_turn_checkpoint(session_id: str, payload: dict[str, Any]) -> Path | None:
    """Atomically persist a lightweight mid-turn checkpoint.

    ``payload`` should stay cheap: the turn's user message, the partial
    response text so far, the tool-call count and a timestamp. Returns
    the checkpoint path, or ``None`` when there is no session id or the
    write failed (callers treat this as best-effort).
    """
    if not session_id:
        return None
    _ensure_dir()
    data = dict(payload or {})
    data.setdefault("ts", time.time())
    data["session_id"] = session_id
    p = _turn_checkpoint_path(session_id)
    try:
        _atomic_write_text(p, json.dumps(data, ensure_ascii=False))
    except OSError:
        return None
    return p


def load_turn_checkpoint(session_id: str) -> dict[str, Any] | None:
    """Return the mid-turn checkpoint dict, or ``None`` if absent/corrupt."""
    if not session_id:
        return None
    p = _turn_checkpoint_path(session_id)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def clear_turn_checkpoint(session_id: str) -> bool:
    """Remove the mid-turn checkpoint. True if a file was removed."""
    if not session_id:
        return False
    p = _turn_checkpoint_path(session_id)
    try:
        if p.exists():
            p.unlink()
            return True
    except OSError:
        pass
    return False


def consume_crash_recovery_note(session_id: str) -> str | None:
    """One-shot: turn a surviving mid-turn checkpoint into a recovery note.

    Called on session load. If a checkpoint exists AND is newer than the
    session's last full save, the previous process died mid-turn — build
    the "[recovered] ..." note (interrupted task, tool-call count, the
    partial answer) for injection into the restored conversation. The
    checkpoint is cleared in every case (stale ones silently), so the
    note fires exactly once. Returns ``None`` when there is nothing to
    recover.
    """
    ck = load_turn_checkpoint(session_id)
    if not ck:
        return None
    clear_turn_checkpoint(session_id)
    saved_at = 0.0
    saved = load_session(session_id)
    if saved:
        try:
            saved_at = float(saved.get("updated_at", 0) or 0)
        except (TypeError, ValueError):
            saved_at = 0.0
    try:
        ck_ts = float(ck.get("ts", 0) or 0)
    except (TypeError, ValueError):
        ck_ts = 0.0
    if ck_ts <= saved_at:
        return None   # stale — a full save landed after this checkpoint
    user_msg = str(ck.get("user_message", "") or "").strip()
    n_tools = int(ck.get("tool_calls", 0) or 0)
    partial = str(ck.get("partial_response", "") or "").strip()
    note = (
        "[recovered] The previous session ended mid-turn while working on: "
        f"{user_msg or '(unknown task)'}. Partial progress ({n_tools} tool "
        "calls, partial answer below) was not committed — verify state "
        "before continuing."
    )
    if partial:
        note += "\n\n" + partial
    return note


def _transcript_archive_dir() -> Path:
    """Directory where full pre-compaction transcripts are archived so
    long sessions never truly lose information — the user can scroll
    back even after the in-memory ``messages`` list was compacted."""
    p = Path.home() / ".delfin" / "transcript_archive"
    p.mkdir(parents=True, exist_ok=True)
    return p


def archive_pre_compaction_transcript(
    session_id: str,
    messages: list[dict[str, Any]],
    *,
    info: dict[str, Any] | None = None,
) -> Path:
    """Persist the full pre-compaction conversation to JSONL so the user
    can browse it later (e.g. via /session archive ls). Returns the
    archive file path.

    The file is append-only: every compaction in a session adds a new
    record block at the end with ``info`` (msg count, tokens saved,
    timestamp). Cheap to read tail-first; safe to leave for months.
    """
    if not session_id:
        return Path(os.devnull)
    d = _transcript_archive_dir()
    p = d / f"{session_id}.jsonl"
    rec = {
        "compacted_at": time.time(),
        "n_messages": len(messages),
        "info": info or {},
        "messages": [
            {
                "role": m.get("role", ""),
                "content": (m.get("content", "") if isinstance(m.get("content"), str)
                            else json.dumps(m.get("content"), ensure_ascii=False)),
            }
            for m in messages
        ],
    }
    try:
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        _chmod_user_only(p)
    except OSError:
        return Path(os.devnull)
    return p


def list_transcript_archives() -> list[dict[str, Any]]:
    """Return one record per archived session: session_id, size, mtime,
    n_compaction_records (one line per compaction event)."""
    d = _transcript_archive_dir()
    out: list[dict[str, Any]] = []
    for p in d.glob("*.jsonl"):
        try:
            stat = p.stat()
            n_lines = sum(1 for _ in p.open("r", encoding="utf-8"))
        except OSError:
            continue
        out.append({
            "session_id": p.stem,
            "path": str(p),
            "bytes": stat.st_size,
            "mtime": stat.st_mtime,
            "compactions": n_lines,
        })
    out.sort(key=lambda r: r["mtime"], reverse=True)
    return out


def load_transcript_archive(session_id: str) -> list[dict[str, Any]]:
    """Return every compaction record in order (oldest → newest) for the
    given session. Each record contains the messages snapshot at that
    compaction. Empty list when no archive exists."""
    if not session_id:
        return []
    p = _transcript_archive_dir() / f"{session_id}.jsonl"
    if not p.exists():
        return []
    out: list[dict[str, Any]] = []
    try:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    except OSError:
        return []
    return out


# ---------------------------------------------------------------------------
# Lossless elision store — full originals of in-place-trimmed message bodies
# ---------------------------------------------------------------------------
#
# The compaction ladder's sliding-window trim and hard-clear stages mutate
# message content IN PLACE; unlike the summary path they never archived the
# original, so the dropped middle was gone for good. Every destructive
# in-place edit now first appends the full original body here, keyed by a
# short ref id that the replacement marker carries, so history_search /
# history_get('elided:<ref>') can page the text back in.

# Documented cap: at most this many elision records are kept per session;
# once exceeded the OLDEST records are dropped (newest always survive).
_ELIDED_CAP = 2000


def elided_store_path(session_id: str) -> Path:
    """Path of the per-session elided-content JSONL (may not exist yet)."""
    return _SESSIONS_DIR / f"{session_id}.elided.jsonl"


def _enforce_elided_cap(path: Path, cap: int | None = None) -> None:
    """Drop the oldest records once the store exceeds the cap.

    Single streaming pass with a bounded tail buffer; the rewrite is
    atomic (temp file + replace) so a crash never tears the store.
    ``cap`` defaults to the module-level :data:`_ELIDED_CAP` at call time
    so tests can shrink it.
    """
    if cap is None:
        cap = _ELIDED_CAP
    try:
        from collections import deque
        tail: deque[str] = deque(maxlen=max(1, int(cap)))
        total = 0
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                total += 1
                tail.append(line)
        if total > len(tail):
            text = "".join(tail)
            if text and not text.endswith("\n"):
                text += "\n"
            _atomic_write_text(path, text)
    except OSError:
        pass


def append_elided_record(
    session_id: str,
    *,
    index: int,
    role: str,
    content: str,
    reason: str = "",
) -> str | None:
    """Append one lossless-elision record; returns its short ref id.

    Record fields: ``ref`` (short id the trim marker carries), ``ts``,
    ``index`` (live message index at elision time), ``role``, ``reason``
    (which compaction stage elided it), ``content`` (the FULL original).

    Same append-mode pattern as the transcript archive: single-line
    appends, and the readers skip torn lines. Best-effort by contract —
    any failure returns ``None`` and the caller keeps a plain marker, so
    a broken store can never break compaction.
    """
    if not session_id:
        return None
    import uuid
    ref = uuid.uuid4().hex[:8]
    rec = {
        "ref": ref,
        "ts": time.time(),
        "index": int(index),
        "role": str(role or ""),
        "reason": str(reason or ""),
        "content": content if isinstance(content, str) else str(content),
    }
    _ensure_dir()
    p = elided_store_path(session_id)
    try:
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        _chmod_user_only(p)
        _enforce_elided_cap(p)
    except OSError:
        return None
    return ref


def load_elided_record(session_id: str, ref: str) -> dict[str, Any] | None:
    """Resolve one elision record by its ref id (newest wins on duplicate).

    Returns the record dict or ``None`` when the store or ref is missing;
    torn/corrupt lines are skipped. Never raises on OS errors.
    """
    if not session_id or not ref:
        return None
    p = elided_store_path(session_id)
    if not p.is_file():
        return None
    found: dict[str, Any] | None = None
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(rec, dict) and rec.get("ref") == ref:
                    found = rec
    except OSError:
        return None
    return found


def save_session(
    session_id: str,
    *,
    mode: str = "quick",
    role_index: int = 0,
    route: list[str] | None = None,
    role_outputs: dict[str, str] | None = None,
    chat_messages: list[dict[str, Any]] | None = None,
    cycle_history: list[dict[str, Any]] | None = None,
    engine_messages: list[dict[str, Any]] | None = None,
    token_usage: dict[str, int] | None = None,
    cost_usd: float = 0.0,
    title: str = "",
    # Long-session extras (all optional, ignored if backend doesn't pass)
    perm_profile: str = "",
    provider: str = "",
    model: str = "",
    effort: str = "",
    active_gate: dict[str, Any] | None = None,
    last_compaction_info: dict[str, Any] | None = None,
    subagent_calls: list[dict[str, Any]] | None = None,
    pending_plan_body: str = "",
    todo_payload: list[dict[str, Any]] | None = None,
    transcript_archive_path: str = "",
    workspace: str = "",
    # Engine state the exporter already produces. Without these two a
    # resumed session forgets the directory it was pinned to, and the
    # context estimator loses its floor — so the window bar reads far
    # emptier than the conversation actually is, until the next turn
    # re-establishes it from the provider.
    project_dir: str = "",
    last_input_tokens: int = 0,
    # Size of the system prompt that drove the last turn. The prompt text
    # itself stays in the process (it carries the injected memory); the
    # estimate only needs its length, and without it a resumed session
    # counts the prompt as nothing.
    system_prompt_chars: int = 0,
    # What stands in for the history a long session has already trimmed
    # away. Dropping these on save means a resumed session keeps the
    # messages it still has and loses the condensed memory of the ones it
    # does not — the part that cannot be reconstructed from anything.
    compaction_summaries: dict[str, Any] | None = None,
    # The ledgers the grounding guards judge against. Without them a
    # resumed session runs the guards in their enforcing branch against an
    # empty record, and corrects the agent for restating what it verified
    # before the save.
    evidence: dict[str, Any] | None = None,
    # Everything else the engine's exporter produces, persisted verbatim.
    # The headless path forwards ``export_state()`` wholesale; without a
    # catch-all here, a field added to the exporter would raise TypeError
    # inside a best-effort save and the whole session would silently not be
    # written. Persisting it unrecognised is strictly better than losing it,
    # and the restore side reads the file by key.
    **extra: Any,
) -> Path:
    """Save a session to disk.

    Parameters
    ----------
    session_id : str
        The CLI session ID (used as filename).
    mode : str
        Active DELFIN_AGENT_LITE mode.
    role_index : int
        Current role index in the route.
    route : list[str]
        The mode's role route.
    role_outputs : dict
        Collected outputs from completed roles.
    chat_messages : list
        The UI chat messages (user/assistant/system with role_label).
    cycle_history : list
        Recent gate / handoff / retry events for the Cycle Inspector.
    engine_messages : list
        The engine's conversation messages (role + content).
    token_usage : dict
        ``{"input": N, "output": N}``.
    cost_usd : float
        Accumulated cost.
    title : str
        Short description (auto-generated from first user message if empty).

    Returns
    -------
    Path
        Path to the saved session file.

    Raises
    ------
    SessionLockedError
        When another live process holds this session's writer lock —
        the save is skipped entirely so two writers on the same session
        id can't silently clobber each other.
    """
    if not session_id:
        return Path(os.devnull)

    d = _ensure_dir()

    # Single-writer guard: acquire/refresh our lock, or bail loudly if a
    # different live process is writing this session id right now.
    acquire_session_lock(session_id)

    # Auto-title from first user message
    if not title and chat_messages:
        for msg in chat_messages:
            if msg.get("role") == "user":
                text = msg.get("content", "")
                title = text[:80].replace("\n", " ").strip()
                if len(text) > 80:
                    title += "..."
                break

    data = {
        # Unrecognised exporter keys first, so a named parameter always
        # wins over a same-named passthrough. Anything that would not
        # survive json.dumps is dropped here rather than taking the whole
        # save down with it.
        **{str(k): v for k, v in (extra or {}).items()
           if not str(k).startswith("_") and _is_json_safe(v)},
        "schema_version": SESSION_SCHEMA_VERSION,
        "session_id": session_id,
        "mode": mode,
        "role_index": role_index,
        "route": route or [],
        "role_outputs": role_outputs or {},
        "chat_messages": chat_messages or [],
        "cycle_history": cycle_history or [],
        "engine_messages": engine_messages or [],
        "token_usage": token_usage or {"input": 0, "output": 0},
        "cost_usd": cost_usd,
        "title": title,
        "updated_at": time.time(),
        # The directory this session worked in. Sessions are listed and
        # resumed per workspace: a session started elsewhere must not be
        # offered — let alone loaded — in another project's dashboard.
        "workspace": str(workspace or ""),
        # Long-session state (only persisted when callers pass it)
        "perm_profile": perm_profile or "",
        "provider": provider or "",
        "model": model or "",
        "effort": effort or "",
        "active_gate": active_gate or None,
        "last_compaction_info": last_compaction_info or None,
        "subagent_calls": subagent_calls or [],
        "pending_plan_body": pending_plan_body or "",
        # Fall back to the real task store when the caller passes nothing:
        # only one of the several save paths ever had a payload to pass,
        # and the two consumers ("Open items" in the handoff brief and in
        # the episode) read this key and nothing else.
        "todo_payload": (todo_payload
                         or tasks_as_todo_payload(workspace, session_id)),
        "transcript_archive_path": transcript_archive_path or "",
        "project_dir": str(project_dir or ""),
        "last_input_tokens": int(last_input_tokens or 0),
        "system_prompt_chars": int(system_prompt_chars or 0),
        "evidence": dict(evidence or {}),
        "compaction_summaries": compaction_summaries or {},
    }

    # Set created_at only on first save
    filepath = d / f"{session_id}.json"
    if filepath.exists():
        try:
            existing = json.loads(filepath.read_text())
            data["created_at"] = existing.get("created_at", data["updated_at"])
        except (json.JSONDecodeError, OSError):
            data["created_at"] = data["updated_at"]
    else:
        data["created_at"] = data["updated_at"]

    _atomic_write_text(filepath, json.dumps(data, ensure_ascii=False, indent=1))
    return filepath


def load_session(session_id: str) -> dict[str, Any] | None:
    """Load a session from disk.

    Returns None if the session file doesn't exist or is corrupt.
    """
    filepath = _SESSIONS_DIR / f"{session_id}.json"
    if not filepath.exists():
        return None
    try:
        return json.loads(filepath.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def list_sessions(
    limit: int = 50, *, workspace: str | Path | None = None,
) -> list[dict[str, Any]]:
    """List saved sessions, newest first.

    Returns a list of lightweight session summaries (no full chat history).

    ``workspace`` scopes the result to the directory the session worked
    in: a project's dashboard offers that project's history, not every
    session on the machine. Sessions saved before the field existed carry
    no workspace and stay visible everywhere, so nothing becomes
    unreachable.
    """
    d = _SESSIONS_DIR
    if not d.exists():
        return []

    want = ""
    if workspace:
        try:
            want = str(Path(workspace).expanduser().resolve())
        except Exception:
            want = str(workspace)

    sessions = []
    for f in d.glob("*.json"):
        if f.name.endswith(".turn.json"):
            continue   # mid-turn crash checkpoint, not a session
        try:
            data = json.loads(f.read_text())
            if want:
                own = str(data.get("workspace") or "")
                if own:
                    try:
                        own = str(Path(own).expanduser().resolve())
                    except Exception:
                        pass
                    if own != want:
                        continue
            sessions.append({
                "workspace": data.get("workspace", ""),
                "session_id": data.get("session_id", f.stem),
                "title": data.get("title", "Untitled"),
                "mode": _migrate_mode(data.get("mode", "quick")),
                "role_index": data.get("role_index", 0),
                "route": data.get("route", []),
                "cost_usd": data.get("cost_usd", 0.0),
                "token_usage": data.get("token_usage", {}),
                "message_count": len(data.get("chat_messages", [])),
                "created_at": data.get("created_at", 0),
                "updated_at": data.get("updated_at", 0),
            })
        except (json.JSONDecodeError, OSError):
            continue

    sessions.sort(key=lambda s: s.get("updated_at", 0), reverse=True)
    return sessions[:limit]


def latest_session(workspace: str | Path | None = None) -> dict[str, Any] | None:
    """Most recently updated session of ``workspace`` (None if there is
    none). The resume target: continuing a project means continuing THAT
    project's last conversation, never whatever ran last on the machine.
    """
    rows = list_sessions(limit=1, workspace=workspace)
    return rows[0] if rows else None


def delete_session(session_id: str) -> bool:
    """Delete a session file. Returns True if deleted."""
    filepath = _SESSIONS_DIR / f"{session_id}.json"
    # Companion files (writer lock, mid-turn checkpoint, elided store) go
    # with it so a deleted session can't leave orphan sidecars behind.
    for side in (_lock_path(session_id), _turn_checkpoint_path(session_id),
                 elided_store_path(session_id)):
        try:
            side.unlink()
        except OSError:
            pass
    if filepath.exists():
        filepath.unlink()
        return True
    return False


def fork_session(
    source_id: str,
    new_id: str | None = None,
    *,
    title_suffix: str = " (fork)",
    sessions_dir: Path | None = None,
) -> str | None:
    """Duplicate an existing session under a new ID.

    The fork inherits all chat / engine / cycle history but gets a fresh
    ``session_id``, fresh timestamps, and an appended title suffix.  The
    ``cost_usd`` and ``token_usage`` counters stay attached so the fork's
    cumulative state is honest.

    Parameters
    ----------
    source_id : str
        The session ID to fork.  Must exist on disk.
    new_id : str, optional
        New session ID.  If omitted, a timestamp-based ID is generated.
    title_suffix : str
        Appended to the title to disambiguate from the source.
    sessions_dir : Path, optional
        Override the default ``~/.delfin/agent_sessions/`` path (for tests).

    Returns
    -------
    str | None
        The new session ID, or ``None`` if the source could not be loaded.
    """
    base_dir = sessions_dir or _SESSIONS_DIR
    src_path = base_dir / f"{source_id}.json"
    if not src_path.exists():
        return None
    try:
        data = json.loads(src_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None

    if not new_id:
        new_id = f"fork-{int(time.time())}-{source_id[-6:] if source_id else 'src'}"

    # Mutate to make this a fresh standalone session
    data["session_id"] = new_id
    src_title = str(data.get("title") or "Untitled")
    if title_suffix and title_suffix not in src_title:
        data["title"] = f"{src_title}{title_suffix}"
    now = time.time()
    data["created_at"] = now
    data["updated_at"] = now
    data["forked_from"] = source_id

    out_path = base_dir / f"{new_id}.json"
    base_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write_text(out_path, json.dumps(data, ensure_ascii=False, indent=1))
    return new_id


def session_lineage(session_id: str) -> list[dict[str, Any]]:
    """Walk the ``forked_from`` chain back to the root and return one
    summary record per ancestor, newest (the given session) first.

    Used by ``/session tree`` to show the branch lineage so the user can
    see where a fork came from. Stops at the first session without a
    ``forked_from`` field or when a parent is missing on disk.
    """
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    current = session_id
    while current and current not in seen:
        seen.add(current)
        data = load_session(current)
        if not data:
            break
        out.append({
            "session_id": current,
            "title": data.get("title", "Untitled"),
            "forked_from": data.get("forked_from", ""),
            "updated_at": data.get("updated_at", 0),
            "message_count": len(data.get("chat_messages") or []),
            "cost_usd": data.get("cost_usd", 0.0),
        })
        parent = data.get("forked_from") or ""
        current = parent
    return out


def session_children(session_id: str) -> list[dict[str, Any]]:
    """Return immediate child sessions (i.e. sessions whose
    ``forked_from`` points at ``session_id``). Newest first."""
    out: list[dict[str, Any]] = []
    for summary in list_sessions(limit=200):
        sid = summary.get("session_id", "")
        if not sid:
            continue
        data = load_session(sid)
        if not data:
            continue
        if data.get("forked_from") == session_id:
            out.append({
                "session_id": sid,
                "title": data.get("title", "Untitled"),
                "updated_at": data.get("updated_at", 0),
                "message_count": len(data.get("chat_messages") or []),
            })
    out.sort(key=lambda r: r.get("updated_at", 0), reverse=True)
    return out


# ---------------------------------------------------------------------------
# Session handoff — pass a session to a fresh agent
# ---------------------------------------------------------------------------
#
# The handoff brief is the .delfin equivalent of the compaction summary
# a long-running agent produces before continuing in a fresh context: a
# structured, self-contained recap that lets a NEW agent (new engine,
# new session id, zero conversation history) pick up exactly where the
# old one left off.

_FILE_PATH_RE = re.compile(
    r"(?:(?<=\s)|(?<=^)|(?<=`))"
    r"((?:[\w./\-]+/[\w./\-]+|[\w./\-]+\.[a-zA-Z]{1,5}))"
    r"(?::\d+)?"
)


def _extract_files_touched(engine_messages: list[dict]) -> list[str]:
    """Best-effort scan of engine messages for file paths the agent
    read / wrote / edited. Deduplicated, capped at 20."""
    seen: list[str] = []
    for msg in engine_messages or []:
        content = msg.get("content", "")
        if not isinstance(content, str):
            continue
        # Only assistant / tool messages mention files in a useful way
        for m in _FILE_PATH_RE.finditer(content):
            path = m.group(1)
            if "://" in path or path.count("/") == 0 and "." not in path:
                continue
            if len(path) < 4 or len(path) > 120:
                continue
            if path not in seen:
                seen.append(path)
            if len(seen) >= 20:
                return seen
    return seen


def build_handoff_brief(data: dict[str, Any]) -> str:
    """Produce a structured, self-contained markdown recap of a session
    suitable for briefing a fresh agent.

    Sections: Goal / Key decisions & state / Files touched / Open items /
    Recommended next step. Pure-extractive (no LLM call) so it's cheap
    and deterministic — the format mirrors what a long-running agent
    emits at compaction time.
    """
    chat = data.get("chat_messages") or []
    engine_msgs = data.get("engine_messages") or []
    title = data.get("title") or data.get("session_id", "untitled")[:40]

    # Goal: first user message
    goal = ""
    for m in chat:
        if m.get("role") == "user":
            goal = str(m.get("content", "")).strip()
            break
    goal = goal[:500] or "(no explicit goal found in this session)"

    # Key decisions & state: last 3 assistant messages, trimmed
    assistant_msgs = [
        str(m.get("content", "")).strip()
        for m in chat if m.get("role") == "assistant" and m.get("content")
    ]
    decisions: list[str] = []
    for a in assistant_msgs[-3:]:
        first_para = a.split("\n\n")[0].replace("\n", " ").strip()
        if first_para:
            decisions.append(first_para[:280])

    # Files touched
    files = _extract_files_touched(engine_msgs)

    # Open items: pending / in_progress / blocked tasks from the persisted
    # payload — or, when the session was saved without one, from the task
    # store itself. A brief whose "Open items" is empty because nobody
    # filled the field reads exactly like a brief for finished work.
    todos = data.get("todo_payload") or tasks_as_todo_payload(
        data.get("workspace") or data.get("project_dir") or "",
        data.get("session_id") or "")
    open_items: list[str] = []
    for t in todos:
        status = t.get("status", "")
        if status in ("pending", "in_progress", "blocked"):
            mark = {"in_progress": "▶", "blocked": "⛔"}.get(status, "○")
            reason = str(t.get("blocked_reason", "") or "")
            open_items.append(
                f"{mark} #{t.get('id','?')} {t.get('subject','')}"
                + (f" (waiting on {reason})" if status == "blocked" and reason
                   else ""))
    # Active gate is also an open item
    gate = data.get("active_gate")
    if isinstance(gate, dict) and gate.get("type"):
        open_items.append(
            f"⏸ paused on {gate.get('type')} gate: {gate.get('title','')}"
        )
    # Pending plan
    if (data.get("pending_plan_body") or "").strip():
        open_items.append("📋 a plan-mode plan is awaiting approval")

    # Recommended next step: gist of the last assistant message
    next_step = ""
    if assistant_msgs:
        last = assistant_msgs[-1]
        next_step = last.split("\n\n")[-1].replace("\n", " ").strip()[:300]
    if not next_step:
        next_step = "(continue from the open items above)"

    cost = data.get("cost_usd", 0.0)
    tok = data.get("token_usage", {}) or {}
    n_msgs = len(chat)
    sid = data.get("session_id", "")

    lines = [
        f"# Session Handoff — {title}",
        "",
        "_Self-contained recap for a fresh agent. You have NO prior "
        "conversation history — everything you need to continue is below._",
        "",
        "## Goal",
        "",
        goal,
        "",
        "## Key decisions & current state",
        "",
    ]
    if decisions:
        for d in decisions:
            lines.append(f"- {d}")
    else:
        lines.append("- (no assistant turns recorded yet)")
    lines += ["", "## Files touched", ""]
    if files:
        for f in files:
            lines.append(f"- `{f}`")
    else:
        lines.append("- (none detected)")
    lines += ["", "## Open items", ""]
    if open_items:
        for o in open_items:
            lines.append(f"- {o}")
    else:
        lines.append("- (nothing explicitly pending)")
    lines += [
        "",
        "## Recommended next step",
        "",
        next_step,
        "",
        "---",
        f"_Generated {time.strftime('%Y-%m-%d %H:%M')} from session "
        f"`{sid}` — {n_msgs} messages, "
        f"{tok.get('input', 0):,} in / {tok.get('output', 0):,} out tokens, "
        f"${cost:.4f}._",
    ]
    return "\n".join(lines)


def _handoffs_dir() -> Path:
    p = Path.home() / ".delfin" / "handoffs"
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_handoff_brief(session_id: str, brief: str) -> Path:
    """Persist a handoff brief to ``~/.delfin/handoffs/<sid>_<ts>.md``."""
    d = _handoffs_dir()
    ts = time.strftime("%Y%m%d_%H%M%S")
    safe_sid = re.sub(r"[^a-zA-Z0-9_-]", "_", session_id or "session")[:40]
    p = d / f"{safe_sid}_{ts}.md"
    _atomic_write_text(p, brief)
    return p


# ---------------------------------------------------------------------------
# Portable session bundles — hand a session to another machine / agent
# ---------------------------------------------------------------------------


# Bundle envelope format — independent of the session schema inside it.
BUNDLE_KIND = "delfin-session-bundle"
BUNDLE_FORMAT_VERSION = 1


def _bundles_dir() -> Path:
    p = Path.home() / ".delfin" / "bundles"
    p.mkdir(parents=True, exist_ok=True)
    return p


def export_bundle(session_id: str) -> Path | None:
    """Pack a session into a single portable ``.delfin-session`` zip.

    Contents:
      - ``session.json``        — the full saved-session state
      - ``transcript.jsonl``    — the pre-compaction transcript archive
                                  (if one exists)
      - ``handoff.md``          — a freshly-generated handoff brief

    Returns the bundle path, or ``None`` if the session doesn't exist.
    """
    data = load_session(session_id)
    if not data:
        return None
    bundle_path = _bundles_dir() / f"{session_id}.delfin-session"
    brief = build_handoff_brief(data)
    archive = _transcript_archive_dir() / f"{session_id}.jsonl"
    with zipfile.ZipFile(bundle_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("session.json", json.dumps(data, ensure_ascii=False, indent=1))
        z.writestr("handoff.md", brief)
        if archive.is_file():
            try:
                z.writestr("transcript.jsonl", archive.read_text(encoding="utf-8"))
            except OSError:
                pass
        # A tiny manifest so importers can sanity-check the bundle. The
        # importer now actually opens it — for its first two releases it
        # did not, so the version written here was decoration.
        z.writestr("manifest.json", json.dumps({
            "kind": BUNDLE_KIND,
            "version": BUNDLE_FORMAT_VERSION,
            "schema_version": SESSION_SCHEMA_VERSION,
            "session_id": session_id,
            "exported_at": time.time(),
        }))
    _chmod_user_only(bundle_path)
    return bundle_path


def import_bundle(bundle_path: Path | str, *, new_id: str | None = None) -> str | None:
    """Ingest a ``.delfin-session`` bundle exported by ``export_bundle``.

    The session lands under a FRESH id (so it never clobbers an existing
    local session) and the transcript archive is restored alongside it.
    Returns the new session id, or ``None`` on failure.

    The manifest the exporter writes is READ here. It carried a version
    from the first release and nothing ever opened it, so a bundle from a
    newer build — or a zip that merely happened to contain a session.json —
    was ingested as if it were ours. A newer envelope or a newer session
    schema is a stated refusal (``SessionSchemaError``), not a silent
    partial import.
    """
    p = Path(bundle_path)
    if not p.is_file():
        return None
    try:
        with zipfile.ZipFile(p, "r") as z:
            names = set(z.namelist())
            if "session.json" not in names:
                return None
            manifest: dict[str, Any] = {}
            if "manifest.json" in names:
                try:
                    manifest = json.loads(
                        z.read("manifest.json").decode("utf-8"))
                except (json.JSONDecodeError, UnicodeDecodeError):
                    manifest = {}
            data = json.loads(z.read("session.json").decode("utf-8"))
            transcript = (
                z.read("transcript.jsonl").decode("utf-8")
                if "transcript.jsonl" in names else ""
            )
    except (zipfile.BadZipFile, json.JSONDecodeError, OSError, KeyError):
        return None

    if isinstance(manifest, dict) and manifest:
        kind = str(manifest.get("kind") or "")
        if kind and kind != BUNDLE_KIND:
            raise SessionSchemaError(
                f"not a session bundle: manifest kind {kind!r}")
        try:
            fmt = int(manifest.get("version") or BUNDLE_FORMAT_VERSION)
        except (TypeError, ValueError):
            fmt = BUNDLE_FORMAT_VERSION
        if fmt > BUNDLE_FORMAT_VERSION:
            raise SessionSchemaError(
                f"bundle envelope v{fmt}; this build understands "
                f"v{BUNDLE_FORMAT_VERSION}. Refusing to import it.")
    try:
        inner = int(data.get("schema_version") or 1)
    except (TypeError, ValueError):
        inner = 1
    if inner > SESSION_SCHEMA_VERSION:
        raise SessionSchemaError(
            f"bundle carries a session of schema v{inner}; this build "
            f"understands v{SESSION_SCHEMA_VERSION}. Refusing to import it.")
    if inner < SESSION_SCHEMA_VERSION:
        data, _applied = migrate_session_data(data, inner)

    sid = new_id or f"import-{int(time.time())}-{(data.get('session_id') or 'x')[-6:]}"
    data["session_id"] = sid
    data["imported_at"] = time.time()
    src_title = str(data.get("title") or "Untitled")
    if "(imported)" not in src_title:
        data["title"] = f"{src_title} (imported)"

    d = _ensure_dir()
    _atomic_write_text(
        d / f"{sid}.json", json.dumps(data, ensure_ascii=False, indent=1)
    )
    if transcript.strip():
        arch = _transcript_archive_dir() / f"{sid}.jsonl"
        _atomic_write_text(arch, transcript)
    return sid


def resume_latest(*, max_age_s: int | None = None) -> dict[str, Any] | None:
    """Return the most recently updated session, or None.

    Parameters
    ----------
    max_age_s : int, optional
        If set, ignore sessions whose ``updated_at`` is older than
        this many seconds. Useful for "resume only if I was just
        working on it" UX.
    """
    summaries = list_sessions(limit=1)
    if not summaries:
        return None
    head = summaries[0]
    if max_age_s is not None and head.get("updated_at"):
        if time.time() - float(head["updated_at"]) > max_age_s:
            return None
    return load_session(str(head["session_id"]))


def cleanup_old_sessions(*, max_age_days: int = 30) -> int:
    """Prune session files older than ``max_age_days``. Returns count."""
    if not _SESSIONS_DIR.exists():
        return 0
    cutoff = time.time() - max_age_days * 86_400
    removed = 0
    for f in _SESSIONS_DIR.glob("*.json"):
        try:
            mtime = f.stat().st_mtime
        except OSError:
            continue
        if mtime < cutoff:
            try:
                f.unlink()
                removed += 1
            except OSError:
                pass
    return removed


def find_sessions(query: str, *, limit: int = 20) -> list[dict[str, Any]]:
    """Substring search over session titles + first user message."""
    needle = query.strip().lower()
    if not needle:
        return list_sessions(limit=limit)
    out: list[dict[str, Any]] = []
    for s in list_sessions(limit=200):
        if needle in str(s.get("title", "")).lower():
            out.append(s)
            continue
        # Cheap dive into the file for first-message match
        full = load_session(str(s["session_id"]))
        if not full:
            continue
        msgs = full.get("chat_messages") or []
        first_user = ""
        for m in msgs:
            if m.get("role") == "user":
                first_user = str(m.get("content", ""))
                break
        if needle in first_user.lower():
            out.append(s)
        if len(out) >= limit:
            break
    return out
