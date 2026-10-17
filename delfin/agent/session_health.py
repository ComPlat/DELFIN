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
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(f) or math.isinf(f):
        return default
    if f < minimum:
        return default
    return f


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
    health = (settings.get("health", {}) if isinstance(settings, dict)
              else {})
    wait_person_min = _threshold(
        health.get("waiting_on_person_min"), DEFAULT_WAIT_ON_PERSON_MIN, 1.0)
    wait_model_min = _threshold(
        health.get("wait_on_model_min"), DEFAULT_WAIT_ON_MODEL_MIN, 1.0)
    stalled_min = _threshold(
        health.get("stalled_min"), DEFAULT_STALLED_MIN, 1.0)
    api_errors = _threshold_int(
        health.get("api_errors"), DEFAULT_API_ERRORS, 1)

    if record.get("paused"):
        return {"verdict": "paused",
                "reason": "a pause flag is set for this session",
                "next": f"resume it with: delfin-agent resume "
                        f"{record.get('name') or record.get('session_id')}"}

    if record.get("pid_alive") is False:
        return {"verdict": "gone",
                "reason": "the process is dead and the session record is left",
                "next": "nothing to do; ignore it or clean the record"}

    errors = record.get("api_errors_in_a_row") or 0
    if isinstance(errors, int) and errors >= api_errors:
        return {"verdict": "erroring",
                "reason": f"{errors} API errors in a row",
                "next": "check the provider; the session retries on its own"}

    approvals = record.get("open_approvals") or []
    if approvals:
        opened = record.get("approval_opened_at")
        age_s = (now - opened) if opened else 0.0
        if age_s >= wait_person_min * 60:
            return {"verdict": "waiting_on_person",
                    "reason": f"open approval for "
                              f"{int(age_s / 60)} min"
                              + (f" ({approvals[0]})" if approvals[0] else ""),
                    "next": "answer the approval or escalate it"}
        return {"verdict": "waiting_on_person",
                "reason": "open approval",
                "next": "answer the approval or escalate it"}

    no_token = record.get("no_token_s") or 0
    if no_token >= wait_model_min * 60:
        return {"verdict": "waiting_on_model",
                "reason": f"no model token for {int(no_token / 60)} min",
                "next": "it continues on its own when the provider answers"}

    turn = record.get("last_turn") or {}
    ended = turn.get("stopped") is True
    if ended:
        if record.get("open_work"):
            idle = record.get("idle_s") or 0
            if idle >= stalled_min * 60:
                return {"verdict": "stalled",
                        "reason": f"turn ended {int(idle / 60)} min ago "
                                  "with open work and nothing waiting",
                        "next": "look at it and nudge it yourself"}
        if not record.get("open_work") and not record.get("resume_planned"):
            idle = record.get("idle_s") or 0
            if idle >= stalled_min * 60:
                return {"verdict": "stopped",
                        "reason": f"the last turn stopped and nothing "
                                  f"resumes it ({int(idle / 60)} min idle)",
                        "next": "resume it yourself: delfin-agent --resume "
                                f"{record.get('session_id')}"}
    return {"verdict": "working",
            "reason": "recent turn activity, nothing waiting",
            "next": "no action needed"}
