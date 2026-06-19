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
* ``status`` — ``pending`` / ``in_progress`` / ``completed`` / ``deleted``
* ``created_at`` / ``updated_at`` — ISO-8601 timestamps

No dependency / blocks tracking, no owners, no priorities. Those add
machinery the agent rarely needs and the user rarely reads.
"""

from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


_VALID_STATUSES = {"pending", "in_progress", "completed", "deleted"}


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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

    # -- internal helpers --------------------------------------------------

    def _load(self) -> dict:
        if not self.path.exists():
            return {"next_id": 1, "tasks": []}
        try:
            return json.loads(self.path.read_text(encoding="utf-8"))
        except Exception:
            # Corrupt file — back it up and start fresh rather than crash
            # the tool. The user can recover from the backup if needed.
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
                    # Guard: blocked tasks cannot transition to
                    # in_progress until every predecessor is completed
                    # or deleted. Tasks already in progress can be
                    # marked completed regardless.
                    if fields.get("status") == "in_progress":
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
        if session_id is not None:
            if session_id == "":
                tasks = []
            else:
                tasks = [t for t in tasks if t.get("session_id", "") == session_id]
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
