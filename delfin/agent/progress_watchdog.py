"""The progress watchdog: a stall's class, its age, and the safe action.

Input: one observation per tick (see :class:`ProgressWatchdog`), on a
clock the caller drives — no sleeps, no thread of its own. Output: the
classification (:func:`classify`), one stream-json event plus one status
line per entered stall (:func:`ProgressWatchdog.tick`), and, past a
configurable limit, the safe action the class allows.

The stall classes and what may be done about each:

``approval-wait``
    A question is parked in the approval room. A person must answer
    it; the session KEEPS WAITING and reports. Never retried, never
    killed, never answered by the machine.
``model-wait``
    A request is out and no byte has come back. The request deadline
    (api_client._request_timeout_s, DELFIN_REQUEST_TIMEOUT_S) already
    retries this; the watchdog only REPORTS it and, past its limit,
    leaves the retry to the existing path. It does not cancel.
``gate-queue``
    The endpoint gate has no free slot (wait_notice is set, the request
    has not started). The wait itself is the mechanism — the gate wait
    replaces waiting inside the endpoint — so the session KEEPS
    WAITING and reports; past the limit it reports louder. Never
    killed: the slot is flock'd and frees itself.
``tool-running``
    A tool call is executing. It may be the user's own long job; the
    watchdog NEVER kills a process it did not start (stall_guard's
    neighbourhood, ops: never kill). It reports and keeps waiting.
``idle``
    Nothing is running at all (no turn, no wait). Reported as the
    session's state, not as a stall.

The limit is ``agent.watchdog.<class>_s`` per class (0 or invalid: no
limit, report only) plus ``agent.watchdog.tick_s`` for the reporting
cadence. Every emitted event carries ``{"type": "stall", "class": ...,
"duration_s": ..., "action": ...}`` so a stream-json reader can alert
without parsing prose.
"""

from __future__ import annotations

import time

# The classes, with the safe action each may take past its limit.
KEEP_WAITING = "keep-waiting"     # approval-wait, gate-queue, tool-running
REPORT_MODEL_RETRY = "report-model-retry"   # model-wait: the deadline retries
NONE = "none"                     # idle: not a stall

CLASSES = ("approval-wait", "model-wait", "gate-queue",
           "tool-running", "idle")

#: Default limits, seconds, per class. 0 = report only, never act.
DEFAULT_LIMITS_S = {
    "approval-wait": 0.0,
    "model-wait": 0.0,
    "gate-queue": 0.0,
    "tool-running": 0.0,
    "idle": 0.0,
}

#: The safe action per class, past the limit (or, with no limit, the
#: action a limit of 0 would allow: reporting alone).
SAFE_ACTION = {
    "approval-wait": KEEP_WAITING,
    "model-wait": REPORT_MODEL_RETRY,
    "gate-queue": KEEP_WAITING,
    "tool-running": KEEP_WAITING,
    "idle": NONE,
}


def _positive(value, default: float) -> float:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    return f if f > 0 else default


def resolve_limits(agent_cfg: dict | None) -> dict[str, float]:
    """``agent.watchdog`` limits over :data:`DEFAULT_LIMITS_S`; malformed
    values keep the default; 0 means report-only."""
    cfg = (agent_cfg or {}).get("watchdog")
    out = dict(DEFAULT_LIMITS_S)
    if isinstance(cfg, dict):
        for key in out:
            if f"{key}_s" not in cfg:
                continue
            try:
                value = float(cfg[f"{key}_s"])
            except (TypeError, ValueError):
                continue
            if isinstance(cfg[f"{key}_s"], bool):
                continue
            out[key] = value if value >= 0 else DEFAULT_LIMITS_S[key]
    return out


def resolve_tick_s(agent_cfg: dict | None, default: float = 30.0) -> float:
    """``agent.watchdog.tick_s``; malformed or non-positive -> *default*."""
    cfg = (agent_cfg or {}).get("watchdog")
    if isinstance(cfg, dict):
        return _positive(cfg.get("tick_s"), default)
    return default


def classify(snapshot, gate_notice: str = "") -> dict:
    """One snapshot (``activity.State`` or a stand-in) -> one stall row.

    Returns ``{"stall": bool, "class": str, "duration_s": float,
    "action": str}``. ``gate_notice`` non-empty says the endpoint gate
    is queueing this request; it outranks the plain model-wait because
    the request has not reached the endpoint yet and the class decides
    what may be done.
    """
    phase = str(getattr(snapshot, "phase", "") or "")
    duration = max(0.0, float(getattr(snapshot, "phase_s", 0.0) or 0.0))
    if phase == "person":
        return {"stall": True, "class": "approval-wait",
                "duration_s": duration, "action": KEEP_WAITING}
    if phase == "model_wait":
        if str(gate_notice or "").strip():
            return {"stall": True, "class": "gate-queue",
                    "duration_s": duration, "action": KEEP_WAITING}
        return {"stall": True, "class": "model-wait",
                "duration_s": duration, "action": REPORT_MODEL_RETRY}
    if phase == "tool":
        return {"stall": True, "class": "tool-running",
                "duration_s": duration, "action": KEEP_WAITING}
    return {"stall": False, "class": "idle", "duration_s": duration,
            "action": NONE}


def status_line(row: dict) -> str:
    """The one status line for a classification row."""
    if not row.get("stall"):
        return ""
    from .activity import fmt_s
    what = {"approval-wait": "waiting for your approval",
            "model-wait": "waiting for the model",
            "gate-queue": "queued at the endpoint gate",
            "tool-running": "a tool is running"}.get(row["class"],
                                                     row["class"])
    return (f"stall: {what} for {fmt_s(row.get('duration_s') or 0)} "
            f"[{row['class']}]")


def event(row: dict) -> dict:
    """The stream-json event for a stall row (its own ``type``: a reader
    replaces it, like the ``wait`` tick, instead of appending)."""
    return {"type": "stall", "class": row["class"],
            "duration_s": row["duration_s"],
            "action": row["action"], "stall": bool(row["stall"])}


class ProgressWatchdog:
    """Drive the classification per tick; emit once per ENTERED stall.

    :meth:`tick` takes the current snapshot (+ gate notice), classifies,
    and — when the stall just began, its class changed, or the next
    report is due — returns ``(event_dict, status_line_str)``. A tick
    that changes nothing returns ``(None, "")``: an unattended session
    must not flood its stream with a row a second. A per-class limit
    past its duration marks :attr:`limit_crossed`, which the caller
    acts on; the watchdog itself never acts beyond reporting.
    """

    def __init__(self, *, now: "callable" = time.monotonic,
                 limits: dict[str, float] | None = None,
                 tick_s: float = 30.0) -> None:
        self._now = now
        self.limits = dict(limits or DEFAULT_LIMITS_S)
        self.tick_s = max(1.0, float(tick_s or 30.0))
        self._current: dict | None = None       # the open stall's row
        self._last_emit = 0.0
        self.limit_crossed = False              # read by the caller

    def tick(self, snapshot, gate_notice: str = "") -> "tuple[dict, str]":
        row = classify(snapshot, gate_notice)
        if not row["stall"]:
            self._current = None
            self._last_emit = 0.0
            self.limit_crossed = False
            return None, ""
        now = self._now()
        if self._current is not None and self._current["class"] == row["class"]:
            row["duration_s"] = max(row["duration_s"],
                                    self._current["duration_s"])
            if now - self._last_emit < self.tick_s:
                return None, ""
        self._last_emit = now
        limit = _positive(self.limits.get(row["class"], 0.0), 0.0)
        self.limit_crossed = limit > 0 and row["duration_s"] >= limit
        self._current = dict(row)
        return event(row), status_line(row)


# --- cross-session view -------------------------------------------------
# The in-process watchdog classifies the live phase (activity.py State).
# Across sessions only the persisted traces exist (ops_interface.sessions
# enriched by session_health), so classify_record maps those rows onto
# the same stall classes. Same classes, weaker evidence: tool-running
# cannot be seen from outside, so it is never reported here.

IDLE_STALL_MIN_S = 3600.0    # idle counts as a stall after an hour


def stall_class(record, now: float) -> dict:
    """One session row (``session_health.collect()``-shaped) -> stall row.

    Returns ``{"stall": bool, "class": str, "age_s": float, "action":
    str}``. Read-only: nothing here mutates. Precedence: open approval
    (approval-wait, age from ``approval_opened_at``) > gate queue
    (``gate_wait_s``) > no tokens back (model-wait, ``no_token_s``) >
    a running phase (tool-running, ``phase``/``phase_s``). States where
    nothing safe can be decided from outside — paused, dead pid,
    repeated api errors — are reported as ``unclassified`` rather than
    guessed; a quiet session (nothing waiting, turn over) is ``idle``
    and NOT a stall.
    """
    if not isinstance(record, dict):
        return {"stall": False, "class": "unclassified", "age_s": 0.0,
                "action": ""}
    paused = bool(record.get("paused"))
    pid_alive = record.get("pid_alive")
    errors = record.get("api_errors_in_a_row") or 0
    if paused or pid_alive is False or errors:
        return {"stall": False, "class": "unclassified", "age_s": 0.0,
                "action": ""}
    approvals = record.get("open_approvals") or []
    if approvals:
        opened = _num(record.get("approval_opened_at"))
        return {"stall": True, "class": "approval-wait",
                "age_s": max(0.0, now - opened) if opened else 0.0,
                "action": SAFE_ACTION["approval-wait"]}
    gate = _num(record.get("gate_wait_s"))
    if gate is not None and gate > 0:
        return {"stall": True, "class": "gate-queue", "age_s": gate,
                "action": SAFE_ACTION["gate-queue"]}
    no_token = _num(record.get("no_token_s"))
    if no_token is not None and no_token > 0:
        return {"stall": True, "class": "model-wait", "age_s": no_token,
                "action": SAFE_ACTION["model-wait"]}
    phase = str(record.get("phase") or "")
    if phase and phase not in ("idle", "done"):
        phase_s = _num(record.get("phase_s")) or 0.0
        return {"stall": True, "class": "tool-running", "age_s": phase_s,
                "action": SAFE_ACTION["tool-running"]}
    return {"stall": False, "class": "idle", "age_s": 0.0, "action": NONE}


#: The cross-session view is the same classes on weaker evidence;
#: :func:`classify_record` is kept as the first name it shipped under.
classify_record = stall_class


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
