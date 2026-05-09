"""Persistent task list for the KIT-Toolbox coding agent.

Mirrors Claude Code's TaskCreate/Update/List planning tools so a
multi-day project (Jerome's Bayesian-opt integration, an iterative
refactor across files) keeps a structured to-do that survives session
restarts. State lives in ``<workspace>/.delfin/session_tasks.json`` by
default — staying with the project rather than the user account so
two simultaneous projects don't share task lists.

Tasks are intentionally simple records:

* ``id`` — integer, monotonically increasing within a workspace
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
               active_form: str = "") -> dict:
        if not subject.strip():
            raise ValueError("subject must be non-empty")
        with self._lock:
            data = self._load()
            tid = int(data.get("next_id", 1))
            data["next_id"] = tid + 1
            now = _now_iso()
            task = {
                "id": tid,
                "subject": subject.strip(),
                "description": description,
                "active_form": active_form or "",
                "status": "pending",
                "created_at": now,
                "updated_at": now,
            }
            data.setdefault("tasks", []).append(task)
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
        allowed = {"status", "subject", "description", "active_form"}
        unknown = set(fields) - allowed
        if unknown:
            raise ValueError(f"unknown field(s): {sorted(unknown)}")
        with self._lock:
            data = self._load()
            for t in data.get("tasks", []):
                if int(t["id"]) == int(task_id):
                    t.update({k: v for k, v in fields.items() if v is not None})
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

    def list(self, *, include_deleted: bool = False) -> list[dict]:
        with self._lock:
            data = self._load()
        tasks = data.get("tasks", []) or []
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
