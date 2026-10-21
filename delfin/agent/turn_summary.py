"""The end of a turn, told like Claude Code tells it.

``format_turn_summary`` builds the dimmed line the person sees after
every turn:

    Surfaced after 31s · done 11:11 · 1 shell still running

It is PURE: no clock, no registry, no I/O. The caller (the CLI footer
in :mod:`delfin.agent.cli_stream`, and the dashboard hook) supplies the
duration, the end time, the turn's outcome and how many background
jobs are still running. Only counts ever reach the string — a job
name, a command or a hostile label can never appear in the line.
"""

# The rotating verb set behind ``agent.ui.turn_verbs = "rotate"``.
_ROTATING_VERBS = ("Dived", "Surfaced", "Distilled", "Crystallized",
                   "Echolocated")
# A turn that did not end well never claims a success verb — "Stopped
# after 12s" says what happened instead.
_SUCCESS_VERB = "Surfaced"
_FAILURE_VERB = "Stopped"

_OK_OUTCOMES = frozenset({"ok", "success", "done", ""})


def format_duration(duration_s) -> str:
    """``31`` -> ``31s``, ``133`` -> ``2m 13s``, ``3840`` -> ``1h 4m``.

    A non-finite duration (inf/nan) renders nothing rather than
    crashing: OverflowError on round() is caught like TypeError.
    """
    try:
        seconds = int(round(float(duration_s)))
    except (TypeError, ValueError, OverflowError):
        return ""
    if seconds < 0:
        seconds = 0
    if seconds < 60:
        return f"{seconds}s"
    hours, rem = divmod(seconds, 3600)
    if hours:
        return f"{hours}h {rem // 60}m"
    return f"{rem // 60}m {rem % 60}s"


def _verb(outcome: str, verb_mode: str, pick: int) -> str:
    if str(outcome or "").strip().lower() not in _OK_OUTCOMES:
        return _FAILURE_VERB
    if str(verb_mode or "").strip().lower() == "rotate":
        return _ROTATING_VERBS[pick % len(_ROTATING_VERBS)]
    return _SUCCESS_VERB


def format_turn_summary(duration_s, ended_at, outcome,
                        running_jobs=None, verb_mode: str = "fixed",
                        pick: int = 0) -> str:
    """One line summarising how a turn ended.

    ``duration_s``: float seconds or None (omitted).
    ``ended_at``: local ``datetime`` of the turn's end, or None.
    ``outcome``: ``"ok"`` for a turn that ended without an error;
    anything else (``"error"``, ``"stopped"``, unknown, hostile) gets
    the failure verb, never a success claim.
    ``running_jobs``: count of background jobs of THIS session still
    running; 0 or None omits the part. A string by mistake is counted
    as one job and never echoed as text.
    ``verb_mode``: ``"fixed"`` (default) or ``"rotate"`` behind
    ``agent.ui.turn_verbs``.
    ``pick``: the rotating verb index for ``verb_mode="rotate"`` — the
    caller derives it, so this stays a pure function.
    """
    parts: list[str] = []

    verb = _verb(outcome, verb_mode, pick)
    duration = format_duration(duration_s)
    if duration:
        parts.append(f"{verb} after {duration}")

    if ended_at is not None:
        try:
            hhmm = ended_at.strftime("%H:%M")
        except AttributeError:
            hhmm = ""
        if hhmm:
            parts.append(f"done {hhmm}")

    # Only the COUNT of jobs is rendered, never a name or a command:
    # hostile text in a job label has no path into this line.
    try:
        jobs = int(running_jobs)
    except (TypeError, ValueError):
        # A non-numeric value (e.g. a label passed by mistake) counts
        # as one running job rather than silently disappearing.
        jobs = 1 if running_jobs else 0
    if jobs > 0:
        noun = "shell" if jobs == 1 else "jobs"
        parts.append(f"{jobs} {noun} still running")

    if not parts:
        return "Turn ended"
    return " · ".join(parts)
