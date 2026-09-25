"""The watchpost's own scheduler entry in ~/.delfin/cron.json.

enable/disable touch only this one entry -- other entries in the file
stay byte-identical. Nothing runs on its own: without `enable` there
is no schedule at all.
"""
from __future__ import annotations

import json
from pathlib import Path

from delfin.agent.scheduler import get_scheduler

PROMPT_MARKER = "watchpost run"


def _cron_path() -> Path:
    return Path.home() / ".delfin" / "cron.json"


def enable(every_minutes: int) -> dict:
    """Add (or replace) the watchpost interval entry."""
    every = int(every_minutes) * 60
    if every < 60:
        raise ValueError("every_minutes must be >= 1")
    sched = get_scheduler(_cron_path())
    # remove a previous watchpost entry first, so enable is idempotent
    disable(sched=sched)
    entry = sched.schedule_interval(
        every_seconds=every,
        prompt=f"watchpost run: scan the account for intrusion traces "
               f"and report (read-only)",
        reason="watchpost look-out, enabled by the user")
    return {"id": entry.id, "every_seconds": every}


def disable(sched=None) -> bool:
    """Remove the watchpost entry. True when one was removed."""
    sched = sched or get_scheduler(_cron_path())
    removed = False
    for ent in list(sched.list_entries()):
        if PROMPT_MARKER in ent.prompt:
            sched.delete(ent.id)
            removed = True
    return removed
