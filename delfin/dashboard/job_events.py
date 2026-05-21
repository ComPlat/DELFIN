"""Diff + format helpers for proactive job-event notifications (D5).

The dashboard polls the active backend's ``list_jobs()`` periodically.
This module turns those snapshots into *events* (status transitions)
and decides which ones deserve a system message in the agent chat
("Job X just failed — want me to investigate?").

Pure helpers — no widgets, no threads — so they can be tested without
a notebook kernel or live backend.

Status taxonomy (matches what backend_local + backend_slurm emit):
    PENDING      job is queued, not running yet
    RUNNING      job is actively executing
    COMPLETED    job finished successfully
    FAILED       job exited with error / non-zero
    CANCELLED    user (or scheduler) cancelled the job
    TIMEOUT      job hit its wall-clock limit

A *terminal* status is anything except PENDING/RUNNING.  A *notable*
transition is any move INTO a terminal state — that's when the agent
should react.  Already-terminal jobs in the first snapshot are NOT
events (we don't notify on jobs we never saw running).
"""
from __future__ import annotations

from dataclasses import dataclass


_TERMINAL_STATES: frozenset[str] = frozenset({
    "COMPLETED",
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
})

_RUNNING_STATES: frozenset[str] = frozenset({"PENDING", "RUNNING"})


@dataclass(frozen=True)
class JobEvent:
    """One status transition for a single job."""

    job_id: str
    name: str
    job_dir: str
    from_status: str   # "" if first time we see the job
    to_status: str     # current status
    is_terminal: bool  # True if to_status is in _TERMINAL_STATES

    @property
    def transition(self) -> str:
        if self.from_status:
            return f"{self.from_status}→{self.to_status}"
        return self.to_status


def diff_job_states(
    before: dict, after: dict,
) -> list[JobEvent]:
    """Diff two ``{job_id: JobInfo}``-style maps into events.

    Both ``before`` and ``after`` accept anything with ``status``,
    ``name``, ``job_dir`` attributes (real ``JobInfo``) OR plain dicts
    with the same keys — handy for tests.

    Returns events for:
      - new jobs (from_status == "")
      - jobs whose status changed since ``before``
    Jobs absent from ``after`` are silently dropped (the user can see
    that in the job-status tab; not chat-worthy).
    """
    out: list[JobEvent] = []
    for job_id, info in (after or {}).items():
        new_status = _attr(info, "status") or ""
        if not new_status:
            continue
        prev_info = (before or {}).get(job_id)
        old_status = _attr(prev_info, "status") if prev_info is not None else ""
        if old_status == new_status:
            continue
        out.append(JobEvent(
            job_id=str(job_id),
            name=str(_attr(info, "name") or ""),
            job_dir=str(_attr(info, "job_dir") or ""),
            from_status=str(old_status or ""),
            to_status=str(new_status),
            is_terminal=new_status in _TERMINAL_STATES,
        ))
    return out


def _attr(obj, name: str):
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj.get(name)
    return getattr(obj, name, None)


def is_notable_event(event: JobEvent) -> bool:
    """Filter for events worth a chat message.

    Notable = transition INTO a terminal state from a non-terminal one.
    First-time appearances in PENDING/RUNNING are *not* notable (the
    user just submitted them — no news).  First-time appearances
    that are already terminal are also not notable (we missed the run
    — surfacing them as "events" is misleading).
    """
    if not event.is_terminal:
        return False
    # Skip if we never observed the job in a non-terminal state
    if event.from_status not in _RUNNING_STATES:
        return False
    return True


def format_job_event_message(event: JobEvent) -> str:
    """Compose a short German chat message describing the event.

    Includes a concrete next-step suggestion appropriate to the
    terminal status, so the user can either reply "ja" / approve via
    UI or ignore the notification.
    """
    name = event.name or event.job_id
    label = f"Job '{name}' ({event.job_id})"
    folder = f" in {event.job_dir}" if event.job_dir else ""
    if event.to_status == "COMPLETED":
        action = (
            "Soll ich die Energien analysieren oder einen Recalc-Check "
            "durchführen?"
        )
        return f"✅ {label} ist erfolgreich beendet{folder}.  {action}"
    if event.to_status == "FAILED":
        action = (
            "Soll ich die Output-Datei lesen und einen Recalc mit "
            "angepassten Parametern vorschlagen?"
        )
        return f"❌ {label} fehlgeschlagen{folder}.  {action}"
    if event.to_status == "TIMEOUT":
        action = "Soll ich einen Recalc mit höherem time_limit vorschlagen?"
        return f"⏱️ {label} hat das Zeitlimit erreicht{folder}.  {action}"
    if event.to_status == "CANCELLED":
        # Less actionable — user usually triggered it themselves
        return f"🚫 {label} wurde abgebrochen{folder}."
    return f"ℹ️ {label}: {event.transition}{folder}."


def snapshot_jobs(jobs) -> dict:
    """Helper: turn a ``list[JobInfo]`` into ``{job_id: JobInfo}``."""
    if not jobs:
        return {}
    out: dict = {}
    for job in jobs:
        jid = _attr(job, "job_id")
        if jid is None:
            continue
        out[str(jid)] = job
    return out


__all__ = [
    "JobEvent",
    "diff_job_states",
    "is_notable_event",
    "format_job_event_message",
    "snapshot_jobs",
]
