"""The cluster's jobs, reachable from a terminal and from a session.

``job_monitor`` already owned the whole mechanism — the shared and
per-workspace watch lists, the three-state scheduler query
(:func:`delfin.agent.job_monitor.query_job_states_detailed`), the daemon
with its PID lock — but only the dashboard had a window into it. The
``delfin-agent`` CLI had no subcommand, and a session's chat had no
slash command (background bash jobs had ``/bash``; SLURM jobs had
nothing). This module is the missing SURFACE, not a second mechanism:
every state below comes from the same query and files the daemon uses.

Two rules inherited from job_monitor:

* **Read-only by default.** Listing answers when asked — there is no
  polling loop in here. The only mutation is the watch flag
  (:func:`watch_set`), and that writes the existing
  ``agent.job_monitor.enabled`` setting and says so; it does not start
  or stop anything. Cancelling a job stays where it already was, behind
  a confirmation (``/cancel <id>``).
* **Universal.** A host without squeue/sacct gets one plain line and a
  carried-on exit — the same three-state honesty the daemon applies to
  one unreachable job, applied to the whole queue. No hostnames, no
  absolute paths baked in.

Everything is injectable (``run_fn``, paths, ``now``) so tests drive
real squeue/sacct output shapes without a scheduler.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Callable, Optional

from . import job_monitor as _jm

# A scheduler that answers nothing: every job's state is unknown, and the
# one line the user sees explains why (mirrors _default_run's contract
# of "None means it could not be run at all", but with the reason kept).
SchedulerUnavailable = FileNotFoundError


def _safe_run(cmd: list[str]) -> Optional[str]:
    """job_monitor._default_run, except a missing binary is NOT hidden.

    FileNotFoundError from subprocess.run is already swallowed into None
    by the daemon's runner — right there, where a degraded daemon keeps
    polling. Here the caller is a human asking once; they need the one
    line that says the host has no scheduler, so the exception carries up.
    """
    import subprocess
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return out.stdout if out.returncode == 0 else None
    except FileNotFoundError:
        raise
    except Exception:
        return None


def collect_job_rows(
    *,
    run_fn: Callable[[list[str]], Optional[str]] = _safe_run,
    shared_watched_path: Optional[Path] = None,
    agent_workspaces: Optional[list[str]] = None,
    now: Callable[[], float] = time.time,
) -> list[dict]:
    """One row per watched SLURM job, running/queued/finished alike.

    Sources are the same two watch lists the daemon polls: the shared
    ``~/.delfin/watched_jobs.json`` (dashboard-written) and each
    workspace's ``<ws>/.delfin/agent_watched_jobs.json``. Only SLURM
    entries are collected — bash jobs have ``/bash``, CI watches have
    their own report. Each row: ``job_id``, ``state`` (a SLURM state,
    ``UNKNOWN`` when the scheduler does not know the id, or
    ``UNAVAILABLE`` when it could not be asked), ``elapsed_s``,
    ``description``, ``workspace``, ``folder``.

    Raises :class:`SchedulerUnavailable` only when the host has no
    scheduler at all AND at least one job would have been listed — an
    empty watch list is not degraded, it is empty.
    """
    rows: list[dict] = []
    wanted = False     # any SLURM id we would have asked about

    shared = _jm.load_watched(
        shared_watched_path) if shared_watched_path is not None \
        else _jm.load_watched()
    for jid, info in (shared.get("jobs") or {}).items():
        info = info or {}
        if not (info.get("kind") in (None, "", "slurm") or jid.isdigit()):
            continue
        wanted = True
        added = info.get("added_at")
        rows.append({
            "job_id": jid,
            "state": None,       # filled by the batch query below
            "elapsed_s": (now() - float(added)
                          if added is not None else None),
            "description": str(info.get("description") or ""),
            "workspace": "",
            "folder": str(info.get("folder") or ""),
        })

    workspaces = (agent_workspaces if agent_workspaces is not None
                  else _jm.agent_watch_workspaces())
    for ws in workspaces:
        data = _jm.load_watched(_jm._agent_watch_path(ws))
        for jid, info in (data.get("jobs") or {}).items():
            info = info or {}
            if info.get("kind") != "slurm":
                continue    # bash/ci watches belong to their own surfaces
            wanted = True
            added = info.get("added_at")
            rows.append({
                "job_id": jid,
                "state": None,
                "elapsed_s": (now() - float(added)
                              if added is not None else None),
                "description": str(info.get("description") or ""),
                "workspace": str(ws),
                "folder": str(info.get("folder") or ""),
            })

    if not wanted:
        return []

    # run_fn raising SchedulerUnavailable (no squeue/sacct binary at all)
    # propagates to the caller, which prints scheduler_note — one line,
    # never a table of invented states, never a traceback.
    states = _jm.query_job_states_detailed(
        [r["job_id"] for r in rows], run_fn=run_fn)
    for r in rows:
        s = states.get(r["job_id"], _jm.STATE_UNAVAILABLE)
        # The daemon's tri-state, rendered for a human: a state stays as
        # the scheduler spelled it; "not known" and "not asked" both get
        # a word instead of an empty string.
        r["state"] = (s if s and s != _jm.STATE_UNAVAILABLE
                      else "UNKNOWN" if s == "" else _jm.STATE_UNAVAILABLE)

    def _rank(r: dict) -> tuple[int, str]:
        s = r["state"] or ""
        if s in ("RUNNING", "CONFIGURING", "COMPLETING"):
            return (0, r["job_id"])
        if s in ("PENDING", "REQUEUED", "REQUEUE_HOLD", "SUSPENDED"):
            return (1, r["job_id"])
        return (2, r["job_id"])

    rows.sort(key=_rank)
    return rows


def _fmt_elapsed(seconds: Optional[float]) -> str:
    if seconds is None or seconds < 0:
        return "—"
    s = int(seconds)
    if s < 60:
        return f"{s}s"
    m, s = divmod(s, 60)
    if m < 60:
        return f"{m}m{s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m"


def render_jobs(rows: list[dict]) -> str:
    """The one table every surface prints — CLI, /jobs, dashboard.

    Pure function over :func:`collect_job_rows` output so the three
    surfaces cannot drift apart.
    """
    if not rows:
        return ("no jobs are being watched — add one with "
                "/watch add <id> [folder] or the dashboard's watch list")
    lines = [f"{'JOB':<12} {'STATE':<14} {'ELAPSED':>10}  WORKSPACE / DESCRIPTION"]
    for r in rows:
        where = r.get("folder") or r.get("workspace") or ""
        desc = r.get("description") or ""
        tail = f"{where} — {desc}" if (where and desc) else (where or desc)
        lines.append(f"{r['job_id']:<12} {str(r['state']):<14} "
                     f"{_fmt_elapsed(r.get('elapsed_s')):>10}  {tail}")
    return "\n".join(lines)


def scheduler_note(n_jobs: int, exc: BaseException) -> str:
    """The one line a host without squeue gets. Never a traceback."""
    return (f"squeue/sacct are not available on this host "
            f"({type(exc).__name__}) — the states of {n_jobs} watched "
            f"job(s) are unknown, not green.")


# ---------------------------------------------------------------------------
# watch: the daemon's flag and status, without touching the daemon
# ---------------------------------------------------------------------------

def _settings_path() -> Path:
    from delfin.user_settings import get_settings_path
    return get_settings_path()


def watch_report(
    *,
    pid_path: Optional[Path] = None,
    watched_path: Optional[Path] = None,
) -> str:
    """`jobs watch --status`: daemon, setting, watched count, cost note."""
    st = _jm.monitor_status(pid_path)
    cfg = _jm.monitor_settings()
    watched = _jm.load_watched(watched_path).get("jobs", {})
    daemon = (f"running (PID {st['pid']})" if st["running"] else "off")
    setting = ("enabled" if cfg["enabled"]
               else "disabled (agent.job_monitor.enabled=false)")
    lines = [
        f"job monitor daemon: {daemon}",
        f"setting: {setting}",
        f"watched jobs: {len(watched)}",
        f"diagnosis: {'on (costs tokens, switchable: auto_diagnose)' if cfg['auto_diagnose'] else 'off (0 tokens)'}",
    ]
    if not cfg["enabled"]:
        lines.append(
            "turn it on with `delfin-agent jobs watch --on` — the watch "
            "loop itself is LLM-free; only failure diagnosis costs tokens "
            "(auto_diagnose, separately switchable)")
    return "\n".join(lines)


def watch_set(
    on_off: str,
    *,
    pid_path: Optional[Path] = None,
) -> str:
    """`jobs watch --on|--off`: persist agent.job_monitor.enabled.

    Only the setting. It does not launch or signal the daemon — the
    daemon is started where it already was (``/watch start`` in the
    dashboard) and exits by itself when it next reads a disabled
    setting (job_monitor.run_loop's token safety), so `--off` cannot
    strand a billed turn.
    """
    on_off = on_off.strip().lower()
    if on_off not in ("on", "off"):
        return "usage: jobs watch --on | --off | --status"
    from delfin.user_settings import load_settings, save_settings
    settings = load_settings(_settings_path())
    cfg = settings.setdefault("agent", {}).setdefault("job_monitor", {})
    cfg["enabled"] = (on_off == "on")
    save_settings(settings, _settings_path())
    if on_off == "on":
        st = _jm.monitor_status(pid_path)
        extra = ("a daemon is already running" if st["running"] else
                 "start it with `/watch start` in the dashboard")
        return (f"job monitoring enabled (agent.job_monitor.enabled=true). "
                f"The watch loop is LLM-free; {extra}.")
    return ("job monitoring disabled (agent.job_monitor.enabled=false). "
            "A running daemon exits on its next pass — no diagnosis "
            "tokens are spent from now on.")
