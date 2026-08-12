"""Persistent task list for the KIT-Toolbox coding agent.

Implements TaskCreate/Update/List planning tools so a
multi-day project (Jerome's Bayesian-opt integration, an iterative
refactor across files) keeps a structured to-do that survives session
restarts. State lives in ``<workspace>/.delfin/session_tasks.json`` by
default — staying with the project rather than the user account so
two simultaneous projects don't share task lists.

Tasks are intentionally simple records:

* ``id`` — integer, monotonically increasing within a workspace
* ``session_id`` — owning chat/session id within that workspace
* ``subject`` — short imperative title shown in the dashboard list
* ``description`` — what needs to be done (multiline OK)
* ``status`` — ``pending`` / ``in_progress`` / ``blocked`` /
  ``completed`` / ``deleted``
* ``blocked_reason`` — required while ``blocked``: what it waits on
* ``started_at`` — when the task last entered ``in_progress``; the start
  of the window a completion claim is checked against
* ``verified`` / ``verify_note`` — what the completion check found
* ``created_at`` / ``updated_at`` — ISO-8601 timestamps

No owners, no priorities. Those add machinery the agent rarely needs and
the user rarely reads.
"""

from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


_VALID_STATUSES = {
    "pending", "in_progress", "blocked", "completed", "deleted",
}

# Statuses that mean "still outstanding". ``blocked`` is open work too:
# waiting on a user answer, a missing credential or a failed dependency
# is a different fact from "not started" and from "done", and before it
# existed the three were indistinguishable in the list the user reads.
OPEN_STATUSES: tuple[str, ...] = ("in_progress", "pending", "blocked")

# What the agent can advance by itself. ``blocked`` is deliberately not
# here: auto-continue must not spin on work that waits on someone else.
ACTIONABLE_STATUSES: tuple[str, ...] = ("in_progress", "pending")


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def resolve_session_scope(session_id: Any) -> Optional[str]:
    """The value every listing surface must pass to :meth:`TaskStore.list`.

    One resolver because the surfaces disagreed in OPPOSITE directions
    for the same empty id: the model-facing reminder passed ``None`` and
    saw every session's leftovers, while the user's ticker passed ``""``
    and printed "No tasks yet" about the very tasks the prompt was
    listing. The CLI backend mints no session id, so that state is real.

    Empty means UNSCOPED — the whole workspace — which is what
    ``task_list``, ``task_adopt`` and the reminder already document.
    """
    sid = str(session_id or "").strip()
    return sid or None


class TaskStore:
    """Thread-safe JSON-backed task list.

    Constructor takes a directory; the store keeps state in
    ``<dir>/session_tasks.json``. Concurrent calls are serialised
    through a per-instance lock; the on-disk format is the
    authoritative source of truth (every mutation reads-modify-writes
    so two processes operating on the same file don't lose data
    silently — they just race for the last write, which is acceptable
    for a single-user agent).
    """

    def __init__(self, base_dir: Path) -> None:
        self.base_dir = Path(base_dir)
        self.path = self.base_dir / ".delfin" / "session_tasks.json"
        self._lock = threading.Lock()
        # Why the last load produced nothing, when that is not the same
        # fact as "no tasks". A store that cannot be read used to be
        # indistinguishable from a finished plan at every consumer.
        self.last_load_error: str = ""

    # -- internal helpers --------------------------------------------------

    def _load(self) -> dict:
        if not self.path.exists():
            self.last_load_error = ""
            return {"next_id": 1, "tasks": []}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
            self.last_load_error = ""
            return data
        except Exception as exc:
            # Corrupt file — back it up and start fresh rather than crash
            # the tool. The user can recover from the backup if needed.
            # The reason is kept: an empty list from here is a READ
            # FAILURE, and a caller that reports it as "nothing
            # outstanding" is asserting something it does not know.
            self.last_load_error = f"{type(exc).__name__}: {exc}"[:200]
            try:
                bak = self.path.with_suffix(".json.bak")
                bak.write_text(self.path.read_text(encoding="utf-8"))
            except Exception:
                pass
            return {"next_id": 1, "tasks": []}

    def _save(self, data: dict) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False),
                       encoding="utf-8")
        tmp.replace(self.path)

    # -- public API --------------------------------------------------------

    def create(self, subject: str, description: str = "",
               active_form: str = "", session_id: str = "",
               blocked_by: list[int] | None = None) -> dict:
        if not subject.strip():
            raise ValueError("subject must be non-empty")
        with self._lock:
            data = self._load()
            tid = int(data.get("next_id", 1))
            data["next_id"] = tid + 1
            now = _now_iso()
            task = {
                "id": tid,
                "session_id": str(session_id or "").strip(),
                "subject": subject.strip(),
                "description": description,
                "active_form": active_form or "",
                "status": "pending",
                "created_at": now,
                "updated_at": now,
                # Set when the task enters in_progress; the start of the
                # window a completion claim is checked against.
                "started_at": "",
                # Required while blocked: what the task is waiting on.
                "blocked_reason": "",
                # What the completion check found: "verified" / "unmet" /
                # "unchecked". Empty until the task is completed.
                "verified": "",
                "verify_note": "",
                # DAG dependency edges. ``blocked_by`` = predecessors that
                # must reach ``completed`` before this task can leave
                # ``pending``. ``blocks`` = downstream task IDs we keep
                # in sync as a reverse-index for cheap traversal.
                "blocked_by": list(blocked_by or []),
                "blocks": [],
            }
            data.setdefault("tasks", []).append(task)
            # Maintain the reverse index: each blocker's ``blocks`` list
            # gets the new task ID appended.
            for parent_id in (blocked_by or []):
                for t in data["tasks"]:
                    if int(t.get("id", 0)) == int(parent_id):
                        t.setdefault("blocks", []).append(tid)
                        if tid not in t["blocks"]:
                            t["blocks"] = list(dict.fromkeys(t["blocks"]))
                        break
            self._save(data)
            return task

    def update(self, task_id: int, **fields: Any) -> dict:
        if not fields:
            raise ValueError("at least one field to update is required")
        if "status" in fields and fields["status"] not in _VALID_STATUSES:
            raise ValueError(
                f"status must be one of {sorted(_VALID_STATUSES)}, "
                f"got {fields['status']!r}"
            )
        allowed = {
            "status", "subject", "description", "active_form",
            "add_blocked_by", "remove_blocked_by",
            # What a blocked task waits on, and what the completion check
            # found — both written by the task tools, both part of the
            # record the user reads.
            "blocked_reason", "verified", "verify_note",
            # Adopt path: a fresh session takes over a task left open by a
            # previous one (task_adopt) by rewriting its owning session id.
            "session_id",
        }
        unknown = set(fields) - allowed
        if unknown:
            raise ValueError(f"unknown field(s): {sorted(unknown)}")
        add_blockers = fields.pop("add_blocked_by", None) or []
        rem_blockers = fields.pop("remove_blocked_by", None) or []
        with self._lock:
            data = self._load()
            for t in data.get("tasks", []):
                if int(t["id"]) == int(task_id):
                    new_status = fields.get("status")
                    # Guard: a task cannot LEAVE the queue past an
                    # unfinished predecessor. The guard used to cover
                    # in_progress only, so `update(child, "completed")`
                    # walked straight around the DAG without ever being
                    # started -- the one transition the whole dependency
                    # edge exists to order.
                    if new_status in ("in_progress", "completed"):
                        blockers = list(t.get("blocked_by") or [])
                        unmet = [
                            b for b in blockers
                            if any(
                                int(o.get("id", 0)) == int(b)
                                and o.get("status") not in ("completed", "deleted")
                                for o in data["tasks"]
                            )
                        ]
                        if unmet:
                            raise ValueError(
                                f"task #{task_id} is blocked by unfinished "
                                f"task(s): {unmet}"
                            )
                    if new_status == "blocked":
                        reason = str(
                            fields.get("blocked_reason")
                            or t.get("blocked_reason") or ""
                        ).strip()
                        if not reason:
                            raise ValueError(
                                f"task #{task_id}: status='blocked' needs "
                                "blocked_reason — name what it waits on "
                                "(a user answer, a missing credential, a "
                                "failed dependency). A blocked task with "
                                "no reason is a pending one nobody can act "
                                "on."
                            )
                    if new_status == "in_progress":
                        # One task at a time, per session: the ticker and
                        # the per-turn reminder both claim to show what
                        # the agent is ON, and a second parallel
                        # in_progress makes that claim false -- and the
                        # completion window ambiguous.
                        sid = str(t.get("session_id", "") or "")
                        other = next(
                            (o for o in data["tasks"]
                             if int(o.get("id", 0)) != int(task_id)
                             and str(o.get("session_id", "") or "") == sid
                             and o.get("status") == "in_progress"),
                            None,
                        )
                        if other is not None:
                            raise ValueError(
                                f"task #{other.get('id')} "
                                f"({str(other.get('subject', ''))[:60]}) is "
                                "already in_progress — finish it "
                                "(status='completed') or park it "
                                "(status='blocked' with blocked_reason) "
                                f"before starting #{task_id}"
                            )
                        if t.get("status") != "in_progress":
                            t["started_at"] = _now_iso()
                    if new_status == "completed" and t.get("status") not in (
                            "in_progress", "completed"):
                        # No silent pending -> completed. The step is what
                        # makes the work window exist at all; without it a
                        # completion claim has nothing to be checked
                        # against, and the list never showed the user what
                        # was being worked on.
                        raise ValueError(
                            f"task #{task_id} is '{t.get('status')}' — mark "
                            "it in_progress before completed, so the list "
                            "shows what you are on and the work has a "
                            "recorded window"
                        )
                    t.update({k: v for k, v in fields.items() if v is not None})
                    # Apply dependency edits + keep reverse index in sync
                    if add_blockers or rem_blockers:
                        existing = list(t.get("blocked_by") or [])
                        for b in rem_blockers:
                            try:
                                existing.remove(int(b))
                            except ValueError:
                                pass
                            for parent in data["tasks"]:
                                if int(parent.get("id", 0)) == int(b):
                                    blk = list(parent.get("blocks") or [])
                                    if task_id in blk:
                                        blk.remove(task_id)
                                        parent["blocks"] = blk
                                    break
                        for b in add_blockers:
                            if int(b) in existing or int(b) == task_id:
                                continue
                            existing.append(int(b))
                            for parent in data["tasks"]:
                                if int(parent.get("id", 0)) == int(b):
                                    blk = list(parent.get("blocks") or [])
                                    if task_id not in blk:
                                        blk.append(task_id)
                                    parent["blocks"] = blk
                                    break
                        t["blocked_by"] = existing
                    t["updated_at"] = _now_iso()
                    self._save(data)
                    return t
            raise KeyError(f"task #{task_id} not found")

    def get(self, task_id: int) -> Optional[dict]:
        with self._lock:
            data = self._load()
            for t in data.get("tasks", []):
                if int(t["id"]) == int(task_id):
                    return t
        return None

    def find_stalled(
        self,
        *,
        max_age_s: float = 600.0,
        session_id: str | None = None,
    ) -> list[dict]:
        """Return tasks that are ``in_progress`` but haven't been updated
        for ``max_age_s`` seconds.  Used by the dashboard's
        ``_record_turn_outcome`` to surface a one-shot suggestion to
        switch into plan mode when a step seems stuck.

        Returned tasks have an extra ``age_s`` field so the caller can
        log the worst offender without recomputing.
        """
        import datetime as _dt
        now = _dt.datetime.utcnow()
        out: list[dict] = []
        for t in self.list(session_id=session_id):
            if t.get("status") != "in_progress":
                continue
            ts = t.get("updated_at") or t.get("created_at") or ""
            try:
                parsed = _dt.datetime.fromisoformat(ts.rstrip("Z"))
            except Exception:
                continue
            age = (now - parsed).total_seconds()
            if age >= max_age_s:
                out.append({**t, "age_s": int(age)})
        out.sort(key=lambda r: r.get("age_s", 0), reverse=True)
        return out

    def list(
        self, *, include_deleted: bool = False, session_id: str | None = None,
        with_seq: bool = False,
    ) -> list[dict]:
        with self._lock:
            data = self._load()
        tasks = data.get("tasks", []) or []
        # Empty and None mean the same thing — UNSCOPED — through the one
        # resolver every caller shares. They used to be opposites here,
        # which is how the prompt and the panel described the same store
        # in contradictory ways.
        scope = resolve_session_scope(session_id)
        if scope is not None:
            tasks = [t for t in tasks
                     if str(t.get("session_id", "") or "") == scope]
        if with_seq:
            # Session-relative 1-based ordinal in creation order (id asc), so
            # the user sees a small "task 3" instead of the global, ever-rising
            # id ("task 90" — bug 20260619-172400). The global id stays the
            # reference for task_update/task_get/blocked_by; seq is display-only
            # and stable as tasks complete (completed tasks keep their slot).
            ordered = sorted(tasks, key=lambda t: int(t.get("id", 0) or 0))
            seq_by_id = {
                int(t.get("id", 0) or 0): i for i, t in enumerate(ordered, 1)
            }
            tasks = [
                {**t, "seq": seq_by_id.get(int(t.get("id", 0) or 0))}
                for t in tasks
            ]
        if not include_deleted:
            tasks = [t for t in tasks if t.get("status") != "deleted"]
        return tasks


# A simple registry keyed by base_dir so a single dashboard session can
# hold per-workspace stores without re-loading from disk on every call.
_STORES: dict[str, TaskStore] = {}
_STORES_LOCK = threading.Lock()


def get_store(base_dir: Path) -> TaskStore:
    """Return the cached TaskStore for ``base_dir``, creating if needed."""
    key = str(Path(base_dir).resolve())
    with _STORES_LOCK:
        store = _STORES.get(key)
        if store is None:
            store = TaskStore(Path(base_dir))
            _STORES[key] = store
        return store


def open_task_summary(
    base_dir: Path, session_id: Any = None, *, cap: int = 6,
) -> dict:
    """Tri-state view of one session's outstanding work.

    ``state`` is ``"open"`` / ``"none"`` / ``"unknown"``. The third value
    is the point: the boolean this replaces failed CLOSED, so a task
    store that could not be read reported exactly what a finished plan
    reports — nothing outstanding — and the turn ended on it.

    Returns ``{"state", "in_progress", "pending", "blocked", "error"}``
    where the three lists hold at most ``cap`` short summaries each
    (``{"id", "seq", "subject", "blocked_reason"}``) and ``counts`` holds
    the full totals. Never raises.
    """
    out: dict = {
        "state": "unknown", "in_progress": [], "pending": [], "blocked": [],
        "counts": {"in_progress": 0, "pending": 0, "blocked": 0},
        "error": "",
    }
    try:
        store = get_store(Path(base_dir))
        tasks = store.list(
            session_id=resolve_session_scope(session_id), with_seq=True)
        err = str(getattr(store, "last_load_error", "") or "")
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"[:200]
        return out
    if err:
        out["error"] = err
        return out
    for t in tasks or []:
        status = str(t.get("status", "") or "")
        if status not in OPEN_STATUSES:
            continue
        out["counts"][status] = out["counts"].get(status, 0) + 1
        if len(out[status]) < max(1, int(cap)):
            out[status].append({
                "id": t.get("id"),
                "seq": t.get("seq"),
                "subject": str(t.get("subject", ""))[:80],
                "blocked_reason": str(t.get("blocked_reason", ""))[:80],
            })
    out["state"] = "open" if any(out["counts"].values()) else "none"
    return out


def _task_label(t: dict) -> str:
    """``task 2 (id 7)`` — the session-relative number the user reads,
    with the global id task_update/task_get key off."""
    seq = t.get("seq")
    return (f"task {seq} (id {t.get('id')})" if seq is not None
            else f"task #{t.get('id')}")


def format_open_tasks_notice(summary: dict) -> str:
    """One user-visible block for a turn that is ending on open work.

    Empty string when there is demonstrably nothing outstanding. The
    turn used to end silently in exactly the two cases this covers —
    open tasks the agent stopped short of, and a task ledger nobody
    could read — so "no line" meant the same as "done".
    """
    try:
        state = str((summary or {}).get("state", "") or "")
        if state == "none":
            return ""
        if state != "open":
            err = str((summary or {}).get("error", "") or "")
            return (
                "⚠ The task list could not be read"
                + (f" ({err})" if err else "")
                + " — I cannot tell whether work is still open. Check "
                  "the task panel before treating this as finished."
            )
        counts = (summary or {}).get("counts") or {}
        n_open = int(counts.get("in_progress", 0)) + int(counts.get("pending", 0))
        n_blocked = int(counts.get("blocked", 0))
        lines = [
            f"⚠ This turn ended with {n_open} open task(s)"
            + (f" and {n_blocked} blocked" if n_blocked else "")
            + " — not everything on the list is done:"
        ]
        for t in (summary.get("in_progress") or []):
            lines.append(f"  ▶ {_task_label(t)} {t.get('subject', '')}")
        for t in (summary.get("pending") or []):
            lines.append(f"  ☐ {_task_label(t)} {t.get('subject', '')}")
        for t in (summary.get("blocked") or []):
            reason = t.get("blocked_reason") or "no reason recorded"
            lines.append(
                f"  ⛔ {_task_label(t)} {t.get('subject', '')} — waiting on "
                f"{reason}")
        shown = (len(summary.get("in_progress") or [])
                 + len(summary.get("pending") or [])
                 + len(summary.get("blocked") or []))
        remainder = n_open + n_blocked - shown
        if remainder > 0:
            lines.append(f"  … +{remainder} more")
        return "\n".join(lines)
    except Exception:
        return ""


def open_foreign_tasks(
    base_dir: Path, current_session_id: str, *, cap: int = 5,
) -> dict:
    """Summary of open (pending / in_progress / blocked) tasks owned by
    OTHER sessions of this workspace.

    The task store is workspace-scoped and outlives sessions, but every
    listing surface filters to the CURRENT session id — so a fresh
    session (fresh id) never sees the work a paused one left behind.
    This helper feeds the one-shot first-turn prompt block that surfaces
    that leftover work for explicit adoption (task_adopt). It never
    merges the lists automatically: session scoping is a deliberate
    leak fix (bug 20260616-183359) and stays authoritative.

    Returns ``{"count", "oldest_age_days", "tasks"}`` where ``tasks``
    holds at most ``cap`` summaries — oldest first, each
    ``{"id", "subject", "status", "session_id", "age_days"}`` — while
    ``count``/``oldest_age_days`` describe ALL foreign open tasks. Age
    is whole days since the task's last update (creation when never
    updated). An empty ``current_session_id`` means the session already
    sees every workspace task, so the summary is empty. Never raises;
    cheap (one JSON read).
    """
    empty: dict = {"count": 0, "oldest_age_days": 0, "tasks": []}
    sid = str(current_session_id or "").strip()
    if not sid:
        return empty
    try:
        tasks = get_store(Path(base_dir)).list(session_id=None)
    except Exception:
        return empty
    now = datetime.now(timezone.utc).replace(tzinfo=None)
    foreign: list[dict] = []
    for t in tasks or []:
        try:
            if t.get("status") not in OPEN_STATUSES:
                continue
            if str(t.get("session_id", "") or "") == sid:
                continue
            ts = str(t.get("updated_at") or t.get("created_at") or "")
            try:
                parsed = datetime.fromisoformat(ts.rstrip("Z"))
                age_days = max(0, int((now - parsed).total_seconds() // 86400))
            except Exception:
                age_days = 0
            foreign.append({
                "id": t.get("id"),
                "subject": str(t.get("subject", ""))[:100],
                "status": t.get("status"),
                "session_id": str(t.get("session_id", "") or ""),
                "age_days": age_days,
            })
        except Exception:
            continue
    if not foreign:
        return empty
    foreign.sort(key=lambda r: (-r["age_days"], int(r["id"] or 0)))
    return {
        "count": len(foreign),
        "oldest_age_days": foreign[0]["age_days"],
        "tasks": foreign[: max(1, int(cap))],
    }
