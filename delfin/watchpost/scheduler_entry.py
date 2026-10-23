"""The watchpost's own scheduler entry in ~/.delfin/cron.json.

enable/disable touch only this one entry -- other entries in the file
stay byte-identical. Nothing runs on its own: without `enable` there
is no schedule at all.
"""
from __future__ import annotations

import json
from pathlib import Path

from delfin.agent.scheduler import Scheduler

# The daemon runs entries with this prefix as a deterministic scan
# (scheduler_daemon.run_watchpost_entry), never as a model turn.
PROMPT_MARKER = "[watchpost]"


def _cron_path() -> Path:
    return Path.home() / ".delfin" / "cron.json"


def _schedule_file() -> Scheduler:
    """The schedule file itself, without a scheduler thread.

    get_scheduler() is the process-wide singleton: it ignores the path
    after its first call and STARTS the scheduler loop. Adding or removing
    one entry from the CLI needs neither.
    """
    return Scheduler(path=_cron_path())


def enable(every_minutes: int) -> dict:
    """Add (or replace) the watchpost interval entry."""
    every = int(every_minutes) * 60
    if every < 60:
        raise ValueError("every_minutes must be >= 1")
    sched = _schedule_file()
    # remove a previous watchpost entry first, so enable is idempotent
    disable(sched=sched)
    entry = sched.schedule_interval(
        every_seconds=every,
        prompt=f"{PROMPT_MARKER} diff scan of the account, report only",
        reason="watchpost look-out, enabled by the user")
    return {"id": entry.id, "every_seconds": every}


def disable(sched=None) -> bool:
    """Remove the watchpost entry. True when one was removed."""
    sched = sched or _schedule_file()
    removed = False
    for ent in list(sched.list_entries()):
        if PROMPT_MARKER in ent.prompt:
            sched.delete(ent.id)
            removed = True
    return removed
