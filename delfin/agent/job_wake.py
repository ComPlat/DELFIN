"""What finished while nobody was looking — in the words both surfaces use.

The dashboard grew the wake-up first, and it grew it in two halves that
lived apart: the producer as a closure inside ``tab_agent``, the renderer
as a module function beside it. They disagreed about the names of the
keys, each half had a passing test, and six live wake-ups read

    [watch] A job you were watching has finished:
    - shell None [?]

That is what a second copy costs. One module owns both halves now, and
the terminal reads the same one rather than growing a third.

Nothing here polls, waits or starts anything: it reports what has already
finished, and the caller decides what that is worth.
"""

from __future__ import annotations

from typing import Any, Iterable, Optional


def finished_shells(seen: set, *, registry: Any = None,
                    session_id: Optional[str] = None) -> list[dict]:
    """Background shells that finished since the last look.

    ``seen`` is the caller's memory across looks; a job is reported once.
    Reported, not drained: the turn this wakes reads the output through
    ``bash_output`` like any other, and an event consumed here would be
    one the agent never sees.

    ``session_id`` is whose shells to report. Without it every session is
    told about every session's jobs, which is what happened on 2026-09-28:
    one benchmark shell finishing woke three sessions. Two of them spent a
    turn working out that the job was not theirs, one wrote to the other
    to ask, and the owner had to answer "please ignore". The two siblings
    of this call in the dashboard's wake tick -- watched jobs and
    background agents -- already take a session; this one did not, and the
    event that leaked was a shell.

    A job with no session recorded is reported to everyone, as before: an
    unowned job is better announced twice than lost.

    Never raises — a wake-up that throws is worse than one that misses.
    """
    out: list[dict] = []
    try:
        if registry is None:
            from . import bash_jobs as _bj
            registry = _bj.get_registry()
        from .bash_jobs import _signal_name
        want = str(session_id or "").strip()
        for job in registry.list_jobs(include_finished=True):
            code = job.poll()
            if code is None or job.job_id in seen:
                continue
            owner = str(getattr(job, "session_id", "") or "")
            if want and owner and owner != want:
                continue
            seen.add(job.job_id)
            if code == 0:
                state = "ok"
            elif code < 0:
                # Not the command's own status: something outside ended
                # it. "exit -9" reads as an ordinary failure.
                state = f"killed by {_signal_name(-code)}"
            else:
                state = f"exit {code}"
            out.append({
                "kind": "shell",
                "job_id": job.job_id,
                "state": state,
                "ok": code == 0,
                "description": str(getattr(job, "command", ""))[:80],
            })
    except Exception:
        return []
    return out


def finished_watched_jobs(workspace: str, seen: set, *,
                          session_id: Optional[str] = None,
                          run_fn: Any = None) -> list[dict]:
    """Agent-registered jobs that reached a terminal state since the last look.

    The terminal's counterpart to the dashboard tick's
    ``check_agent_jobs`` call. ``seen`` is the caller's memory across
    looks (the same set the shells use, ``_wake_seen``), so a job is
    reported once even though ``check_agent_jobs`` with
    ``consume=False`` never removes the entry. ``session_id`` scopes the
    report to this session's own submissions, exactly as the dashboard
    scopes it.

    Job state comes from ``job_monitor.query_job_states_detailed`` via
    ``check_agent_jobs`` — one tri-state scheduler query for all ids,
    never a private squeue loop. The look is pulled by the idle prompt
    behind ``_WAKE_EVERY_S`` (repl), so the throttle lives at the caller.

    Never raises — a wake-up that throws is worse than one that misses.
    """
    try:
        from .job_monitor import check_agent_jobs
        return check_agent_jobs(workspace, run_fn=run_fn,
                                consume=False, marker="wake_notified",
                                session_id=session_id or None)
    except Exception:
        return []


def wake_prompt(done: Iterable[dict]) -> str:
    """The message a finished watched job sends an idle agent; "" for none.

    It says who is speaking. A turn nobody typed arrives through the same
    input the user types into, so without a word to the contrary the
    model reads it as the user asking — and answers them for something
    they never said.
    """
    lines = []
    for ev in done or []:
        line = (f"- {ev.get('kind', 'job')} {ev.get('job_id')} "
                f"[{ev.get('state', '?')}]")
        if ev.get("description"):
            line += f" {str(ev['description'])[:80]}"
        if ev.get("signatures"):
            line += " — " + ", ".join(str(s) for s in ev["signatures"])
        if ev.get("degraded"):
            line += f" — {ev['degraded']}"
        if ev.get("url"):
            line += f" — {ev['url']}"
        lines.append(line)
    if not lines:
        return ""
    return ("[watch — a system event, not the user] A job you were "
            "watching has finished:\n"
            + "\n".join(lines)
            + "\n\nNobody typed this; the job watcher put it here because "
              "you asked to be told. Say what the result means for the "
              "work it was waiting on. If it failed, name the cause from "
              "the evidence before proposing a fix.")


def wake_enabled(settings: Optional[dict] = None) -> bool:
    """``agent.wake_on_job_end``; on unless the user turned it off."""
    try:
        if settings is None:
            from delfin.user_settings import load_settings
            settings = load_settings()
        return bool(((settings or {}).get("agent") or {}).get(
            "wake_on_job_end", True))
    except Exception:
        return True
