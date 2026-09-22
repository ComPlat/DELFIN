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


def finished_shells(seen: set, *, registry: Any = None) -> list[dict]:
    """Background shells that finished since the last look.

    ``seen`` is the caller's memory across looks; a job is reported once.
    Reported, not drained: the turn this wakes reads the output through
    ``bash_output`` like any other, and an event consumed here would be
    one the agent never sees.

    Never raises — a wake-up that throws is worse than one that misses.
    """
    out: list[dict] = []
    try:
        if registry is None:
            from . import bash_jobs as _bj
            registry = _bj.get_registry()
        from .bash_jobs import _signal_name
        for job in registry.list_jobs(include_finished=True):
            code = job.poll()
            if code is None or job.job_id in seen:
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
