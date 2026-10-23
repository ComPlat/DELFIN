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


_QUESTION_TAG = "QUESTION:"

#: Words a denial is spoken with. The turn-end note distinguishes a turn
#: that stopped because a tool was refused from one that stopped because
#: it asked something; the refusal phrases live here. Substrings, lower
#: case, matched against the last part of the turn's text.
_DENIAL_PHRASES = (
    "permission denied",
    "was denied",
    "was refused",
    "not on the auto-allow list",
    "refusing to overwrite",
    "could not be executed",
)


def _ends_with_open_question(text: str) -> bool:
    """Whether the turn's text ends with a question to the user.

    Two shapes: the explicit ``QUESTION:`` tag the role prompt prescribes,
    or a last line that ends in ``?``. A ``?`` anywhere but the end is a
    mention, not a question — the answer stands and nothing is pending.
    """
    t = (text or "").strip()
    if not t:
        return False
    if _QUESTION_TAG in t[-300:]:
        return True
    return t.endswith("?")


def _carries_denial(text: str) -> bool:
    """Whether the turn's text reports a refusal it could not work around."""
    t = (text or "").lower()
    tail = t[-400:]
    return any(p in tail for p in _DENIAL_PHRASES)


def blocked_note(text: str, open_tasks: Iterable[dict],
                 denied: bool = False) -> str:
    """The note a blocked turn leaves at the next prompt; "" for none.

    Wave-10 finding: a session sat an hour waiting at a question it had
    asked, and the next turn started as if nothing were pending. When a
    turn ENDS with an open question or a denial AND open tasks remain,
    the note says so in one line — "blocked on X; open: Y" — so the next
    turn reads the state, not just the task list.

    ``open_tasks`` are the caller's (agent_tasks OPEN_STATUSES entries:
    dicts with a ``subject``); they arrive pre-filtered, so this function
    stays a pure renderer over what it is given and never touches the
    task store itself.
    """
    tasks = list(open_tasks or [])
    if not tasks:
        return ""
    if denied or _carries_denial(text or ""):
        why = "a denied action"
    elif _ends_with_open_question(text or ""):
        why = "an unanswered question"
    else:
        return ""
    names = ", ".join(str(t.get("subject", "") or "")[:60]
                      for t in tasks[:3])
    more = len(tasks) - 3
    if more > 0:
        names += f" … +{more} more"
    return (f"[blocked — a system note, not the user] The turn you just "
            f"finished ended with {why}. You are blocked on that answer; "
            f"open tasks: {len(tasks)} — {names}. Ask again plainly if "
            f"the question is stale, or work a task that does not need it.")


def note_turn_blocked(text: str, open_tasks: Iterable[dict],
                      denied: bool = False) -> str:
    """Turn-end side of :func:`blocked_note`, for the terminal to call.

    Same answer, named for the call site: the terminal records the note
    after a turn and reads it at the next idle prompt. Never raises.
    """
    try:
        return blocked_note(text, open_tasks, denied=denied)
    except Exception:
        return ""


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
