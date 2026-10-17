"""UX2 session health: classify an open session's record into one verdict.

Pure module: records in, ``verdict`` + ``reason`` + ``next`` out. ``now`` is a
parameter, settings are passed in by the caller (never read from disk here),
and every threshold that is 0, negative, text, NaN or infinite falls back to
its default instead of switching a check off.

Wave 17, measured: a stale pause stopped 25 sessions and nobody saw it for
100 min; the watcher's anomaly brake paused a whole wave for 75 min unseen.
The verdicts here exist so the person is told.
"""

from __future__ import annotations

import math

VERDICTS = ("working", "waiting_on_person", "waiting_on_model", "stalled",
            "stopped", "paused", "gone", "erroring")

#: Default thresholds, minutes unless named otherwise. ``agent.health.*``
#: overrides these; the caller resolves them, classify only bounds them.
DEFAULT_WAIT_ON_PERSON_MIN = 15
DEFAULT_WAIT_ON_MODEL_MIN = 10
DEFAULT_STALLED_MIN = 20
DEFAULT_API_ERRORS = 3


def _threshold(value, default, minimum):
    """*value* bounded into [minimum, ...); anything unusable -> *default*."""
    if isinstance(value, bool):
        return default
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(f) or math.isinf(f):
        return default
    if f < minimum:
        return default
    return f


def _num(value):
    """*value* as a float, 0.0 when it is not a number (a malformed record
    field must never make classify raise)."""
    if isinstance(value, bool):
        return 0.0
    try:
        f = float(value)
    except (TypeError, ValueError):
        return 0.0
    if math.isnan(f) or math.isinf(f):
        return 0.0
    return f


#: Characters kept of a sanitised name/id in a suggested command.
_NAME_MAX = 80


def _one_line(value):
    """*value* as ONE safe line: every terminal control character and line
    separator (\n, \r, \v, \f, U+0085, U+2028, U+2029, tab) removed, capped
    at _NAME_MAX. A hostile session name or id must not be able to inject a
    second line into a suggested command or an attention title."""
    s = str(value or "")
    bad = "\n\r\v\f\t\x1b\x85\u2028\u2029"
    return "".join(c for c in s if c.isprintable() and c not in bad)[:_NAME_MAX]


def _threshold_int(value, default, minimum):
    """Like _threshold, but integer-valued (error counts)."""
    f = _threshold(value, default, minimum)
    if f is default:
        return default
    return int(f)


def _now_minus(now, seconds):
    return now - seconds


def classify(record, now, settings):
    """One verdict for one session record.

    *record* is a dict as ``ops_interface.sessions()`` produces, enriched
    with: ``paused`` (bool), ``pid_alive`` (bool), ``approval_opened_at``
    (epoch s of the oldest open approval), ``api_errors_in_a_row`` (int),
    ``no_token_s`` (s since the last token, when it is out), ``open_work``
    (bool: a turn ended with work it said it would continue) and
    ``resume_planned`` (bool: something will wake this session).
    *now* is epoch seconds; *settings* the ``agent`` dict from settings.
    Returns ``{"verdict", "reason", "next"}``.
    """
    health = settings.get("health", {}) if isinstance(settings, dict) \
        else {}
    if not isinstance(health, dict):
        health = {}
    wait_person_min = _threshold(
        health.get("waiting_on_person_min"), DEFAULT_WAIT_ON_PERSON_MIN, 1.0)
    wait_model_min = _threshold(
        health.get("wait_on_model_min"), DEFAULT_WAIT_ON_MODEL_MIN, 1.0)
    stalled_min = _threshold(
        health.get("stalled_min"), DEFAULT_STALLED_MIN, 1.0)
    api_errors = _threshold_int(
        health.get("api_errors"), DEFAULT_API_ERRORS, 1)

    name = _one_line(record.get("name") or record.get("session_id"))
    sid = _one_line(record.get("session_id"))
    if record.get("paused"):
        return {"verdict": "paused",
                "reason": "a pause flag is set for this session",
                "next": f"resume it with: delfin-agent resume {name}"}

    if record.get("pid_alive") is False:
        return {"verdict": "gone",
                "reason": "the process is dead and the session record is left",
                "next": "nothing to do; ignore it or clean the record"}

    errors = _num(record.get("api_errors_in_a_row"))
    if errors >= api_errors:
        return {"verdict": "erroring",
                "reason": f"{errors} API errors in a row",
                "next": "check the provider; the session retries on its own"}

    approvals = record.get("open_approvals") or []
    if approvals:
        opened = _num(record.get("approval_opened_at"))
        age_s = max(0.0, now - opened) if opened else 0.0
        if age_s >= wait_person_min * 60:
            return {"verdict": "waiting_on_person",
                    "reason": f"open approval for "
                              f"{int(age_s / 60)} min"
                              + (f" ({approvals[0]})" if approvals[0] else ""),
                    "next": "answer the approval or escalate it"}
        return {"verdict": "waiting_on_person",
                "reason": "open approval",
                "next": "answer the approval or escalate it"}

    no_token = _num(record.get("no_token_s"))
    if no_token >= wait_model_min * 60:
        return {"verdict": "waiting_on_model",
                "reason": f"no model token for {int(no_token / 60)} min",
                "next": "it continues on its own when the provider answers"}

    turn = record.get("last_turn") or {}
    ended = turn.get("stopped") is True
    if ended:
        if record.get("open_work"):
            idle = _num(record.get("idle_s"))
            if idle >= stalled_min * 60:
                return {"verdict": "stalled",
                        "reason": f"turn ended {int(idle / 60)} min ago "
                                  "with open work and nothing waiting",
                        "next": "look at it and nudge it yourself"}
        if not record.get("open_work") and not record.get("resume_planned"):
            idle = _num(record.get("idle_s"))
            if idle >= stalled_min * 60:
                return {"verdict": "stopped",
                        "reason": f"the last turn stopped and nothing "
                                  f"resumes it ({int(idle / 60)} min idle)",
                        "next": "resume it yourself: delfin-agent --resume "
                                f"{sid}"}
    return {"verdict": "working",
            "reason": "recent turn activity, nothing waiting",
            "next": "no action needed"}


# ---------------------------------------------------------------------------
# collect + watch (the emitting half; classify above stays pure)
# ---------------------------------------------------------------------------

#: Verdicts whose ENTERING is worth an attention event. waiting_on_person
#: is emitted only once it is OLD (the watch's own grace), so a question
#: answered within minutes does not page anyone.
_EMIT_VERDICTS = ("stalled", "stopped", "paused", "gone", "erroring",
                  "waiting_on_person")

#: A waiting_on_person verdict must be at least this old (minutes) before
#: the watch emits for it.
WAITING_EMIT_MIN = 15


def pause_watch(seconds: int) -> None:
    """Sleep between watch ticks, interruptibly: a sleep that a
    KeyboardInterrupt wakes must not swallow the interrupt (a person's
    Ctrl-C stops the loop). Tests pass 0 to keep runs fast."""
    import time as _t
    try:
        _t.sleep(max(0, int(seconds)))
    except KeyboardInterrupt:
        raise


def collect(settings=None):
    """One health record per open session, ready for classify.

    Reads only: presence records (session_presence.open_sessions), pause
    flags (session_pause.is_paused), pending approvals
    (ops_interface.pending_approvals), last turns (ops_interface.sessions,
    which enriches each row). Never writes, never sends.
    """
    from . import ops_interface
    from . import session_pause

    approvals = ops_interface.pending_approvals()
    paused = set(session_pause.paused_keys())
    oldest = {}
    for a in approvals:
        s = str(a.get("session") or "")
        opened = a.get("opened_at") or 0.0
        if opened and (s not in oldest or opened < oldest[s]):
            oldest[s] = opened
    out = []
    seen_keys = set()
    for row in ops_interface.sessions():
        key = str(row.get("name") or "")
        sid = str(row.get("session_id") or "")
        seen_keys.add(key)
        mine = [a for a in approvals
                if a["session"] in (key, sid) or a["session_id"] == sid]
        rec = dict(row)
        rec["paused"] = key in paused or bool(
            key and session_pause.is_paused(key))
        rec["pid_alive"] = True
        rec["open_approvals"] = [a["id"] for a in mine]
        rec["approval_opened_at"] = min(
            (a["opened_at"] for a in mine if a.get("opened_at")),
            default=None)
        rec["title"] = row.get("title") or ""
        rec.setdefault("api_errors_in_a_row", 0)
        rec.setdefault("no_token_s", 0)
        rec.setdefault("open_work", False)
        rec.setdefault("resume_planned", False)
        rec.setdefault("started_at", None)
        out.append(rec)
    # A paused session whose process died leaves only the pause flag: it is
    # exactly the wave-17 case nobody saw. Include it even though no
    # presence record announces it any more.
    for key in sorted(paused - seen_keys):
        out.append({"name": key, "session_id": "", "title": "",
                    "paused": True, "pid_alive": False,
                    "open_approvals": [], "approval_opened_at": None,
                    "last_turn_end": None, "last_turn": {}, "idle_s": None,
                    "api_errors_in_a_row": 0, "no_token_s": 0,
                    "open_work": False, "resume_planned": False,
                    "started_at": None})
    return out


def watch_events(records, now, state, *, settings=None, stop_all=None,
                 stop_all_seen=None, emit=True, _clock=None):
    """One watch tick: emit ONE attention event per (session, verdict)
    entering an emitting state, deduped until the session leaves it.

    *state* is a persistent dict ({"seen": {...}}); the caller owns it and
    passes the same dict on every tick. ``stop_all`` is the stop_all record
    (or None) and ``stop_all_seen`` the watch's dict remembering which stop
    was already announced. Returns the events emitted this tick. With
    ``emit=False`` it plans the events without sending them (the --json
    path uses this to show what WOULD be told).
    """
    health = settings.get("health", {}) if isinstance(settings, dict) else {}
    if not isinstance(health, dict):
        health = {}
    from . import attention
    emitted = []
    seen = state.setdefault("seen", {})
    stop_seen = stop_all_seen if stop_all_seen is not None else {}

    if stop_all:
        at = stop_all.get("at")
        mark = ("stop-all", _num(at))
        if mark not in stop_seen:
            stop_seen[mark] = True
            emitted.append(("stop-all", f"stop-all is SET: every agent of "
                            "yours ends itself (reason: "
                            f"{stop_all.get('reason') or 'not given'}). "
                            "Restart a session yourself when you are back."))
    else:
        stop_seen.clear() if isinstance(stop_seen, dict) else None

    fresh = False
    for rec in records:
        started = _num(rec.get("started_at"))
        at = _num((stop_all or {}).get("at"))
        if at:
            # a fresh session started after the stop: the stop is over
            if started > at and rec.get("pid_alive") and not fresh:
                fresh = True
                if ("cleared", at) not in stop_seen:
                    stop_seen[("cleared", at)] = True
                    emitted.append(("cleared", "stop-all is over: a fresh "
                                    "session started after it. New agents "
                                    "are allowed again."))

        verdict = classify(rec, now, health)
        title = _one_line(rec.get("title") or "")
        name = _one_line(rec.get("name") or rec.get("session_id"))
        grace = WAITING_EMIT_MIN if verdict == "waiting_on_person" else 0
        if verdict not in _EMIT_VERDICTS:
            seen.pop(name, None)
            continue
        if verdict == "waiting_on_person":
            opened = _num(rec.get("approval_opened_at"))
            if opened and now - opened < grace * 60:
                continue
        key = (name, verdict)
        if seen.get(key):
            continue
        seen[key] = True
        emitted.append((name, f"{name} {verdict}: "
                        f"{_one_line(rec.get('reason')) or verdict}. "
                        f"{_one_line(rec.get('next'))}"))
        _ = title
    for name, text in emitted:
        if emit:
            try:
                attention.emit_attention("session_state", session_id=name,
                                         title=text)
            except Exception:
                pass
    _ = _clock
    return emitted
