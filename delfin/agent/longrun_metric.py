"""Long-run quality metrics for one session, reduced from its records.

Package KG4, phase 2. A pure reducer in the shape of autonomy_metric
(delfin/agent/autonomy_metric.py): every input is a record a REAL writer
produced, the module does no I/O, and anything the records cannot attest is
reported as ``None`` — never guessed, never silently a zero.

The questions it answers (the brief):
  * idle time with mail or a due event waiting,
  * turns without progress,
  * nudges delivered while the session was waiting,
  * tool calls written as text,
  * work redone after a compaction or restart,
  * drift (actions outside the stated goal),
  * how often the person steered and whether the steer was followed.

Writers of the inputs (one line each, as in autonomy_metric's docstring):
  * ``turns``       -- turn_metrics.record / turn_metrics.read (real JSONL).
  * ``next_turn_ts``-- the caller's turn-start observation (the timestamp of
                       the turn that FOLLOWED the idle window; None = no
                       observed next turn, so idle is unmeasurable -> None).
  * ``mail_waiting``/``due_event_waiting`` -- the caller's observation at the
                       idle window (session_messages.deliverable / its own
                       scheduler); boolean facts, not guesses.
  * ``nudges``      -- the caller's nudge records {ts, kind}.
  * ``text_calls``  -- the caller's probe of turns where the model emitted a
                       tool call as prose instead of a tool call.
  * ``compactions`` -- compaction_log.record_compaction / read_compactions.
  * ``commands``    -- audit_log.build_changes_report()["commands"] entries
                       extended with ``ts`` by the caller.
  * ``actions``     -- the caller's per-action goal verdicts {ts, in_goal};
                       the reducer never judges in_goal itself.
  * ``steers``      -- the person's steering messages {ts} (writer: the
                       session's own transcript, handed over by the caller).

Progress uses the real writer's own definition: turn_metrics
``_delivered_nothing`` (turn_metrics.py:216), imported directly so the
metric cannot drift from its definition.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from delfin.agent import turn_metrics as tm

from dataclasses import dataclass, field

__all__ = [
    "LongrunInputs",
    "LongrunResult",
    "summarize",
]


@dataclass
class LongrunInputs:
    """One session's raw records. Carried unchanged; nothing is written
    back to the sources. ``None`` means "not observed", not zero."""

    turns: List[Dict[str, Any]] = field(default_factory=list)
    next_turn_ts: Optional[float] = None
    mail_waiting: bool = False
    due_event_waiting: bool = False
    nudges: List[Dict[str, Any]] = field(default_factory=list)
    text_calls: List[Dict[str, Any]] = field(default_factory=list)
    compactions: List[Dict[str, Any]] = field(default_factory=list)
    commands: List[Dict[str, Any]] = field(default_factory=list)
    actions: List[Dict[str, Any]] = field(default_factory=list)
    steers: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class LongrunResult:
    """Per-session picture. ``None`` under a name marks a value the
    session's records cannot attest (unmeasurable), never zero."""

    turns: int = 0
    turns_without_progress: int = 0
    idle_wait_s: Optional[float] = None
    idle_with_mail_s: float = 0.0
    idle_with_due_event_s: float = 0.0
    nudges_while_waiting: int = 0
    text_calls: int = 0
    compactions_seen: int = 0
    redone_after_compaction: int = 0
    drift_actions: int = 0
    drift_fraction: Optional[float] = None
    steers: int = 0
    steers_followed: int = 0


def _ts(record: Dict[str, Any]) -> Optional[float]:
    """Record timestamp as float; None when the record carries none."""
    ts = record.get("ts")
    if ts is None:
        return None
    try:
        return float(ts)
    except (TypeError, ValueError):
        return None


def _delivered_nothing(entry: Dict[str, Any]) -> bool:
    """turn_metrics._delivered_nothing (turn_metrics.py:216) — the real
    one, imported so the metric cannot drift from its definition:
    no text, no tool call, no output token, and no crash to explain it."""
    return tm._delivered_nothing(entry)


def _norm_command(text: Any) -> str:
    return " ".join(str(text or "").split())


def summarize(inputs: LongrunInputs) -> LongrunResult:
    """Reduce one session's records to the long-run picture.

    Pure: no I/O; every fact is either an input or derived from one.
    """
    turns = list(inputs.turns or [])
    res = LongrunResult(turns=len(turns))
    no_progress = sum(1 for e in turns if _delivered_nothing(e))
    res.turns_without_progress = no_progress
    res.text_calls = len(inputs.text_calls or [])

    # --- Idle window: last turn's end -> the next turn's start. The idle
    # time is measured only when BOTH ends are observed; the waiting
    # attribution uses the caller's facts about the window.
    last_ts = _ts(turns[-1]) if turns else None
    if last_ts is not None and inputs.next_turn_ts is not None:
        idle = max(0.0, float(inputs.next_turn_ts) - last_ts)
        res.idle_wait_s = idle
        if inputs.mail_waiting:
            res.idle_with_mail_s = idle
        if inputs.due_event_waiting:
            res.idle_with_due_event_s = idle
        # Nudges delivered INSIDE the idle window, while waiting. Malformed
        # records (missing/non-numeric ts) are pre-filtered like comps and
        # steers below: the records cannot attest them, so they are not
        # counted and must not crash (reviewer nacht-s19, adversarial pass
        # on 1ebf8a8e).
        if idle > 0:
            nudged = (n for n in (inputs.nudges or []) if _ts(n) is not None)
            res.nudges_while_waiting = sum(
                1 for n in nudged
                if last_ts <= _ts(n) <= float(inputs.next_turn_ts))

    # --- Work redone after a compaction or restart: a command that ran
    # before a compaction and the SAME command text that runs again after
    # it. The boundary is real (compaction_log's own records); the match
    # is whitespace-normalized exact text, deliberately narrow so an
    # unrelated similar command is not counted as redone work.
    comps = sorted((c for c in (inputs.compactions or []) if _ts(c) is not None),
                   key=_ts)
    res.compactions_seen = len(comps)
    if comps:
        boundary = _ts(comps[0])
        before = {}
        for c in (inputs.commands or []):
            ts = _ts(c)
            if ts is not None and ts < boundary:
                key = _norm_command(c.get("command"))
                if key:
                    before[key] = ts
        for c in (inputs.commands or []):
            ts = _ts(c)
            if ts is not None and ts > boundary:
                key = _norm_command(c.get("command"))
                if key in before:
                    res.redone_after_compaction += 1

    # --- Drift and steers. The in_goal verdicts come from the caller
    # (the stated goal is its knowledge, not the reducer's). Drift fraction
    # is over actions the caller judged at all.
    actions = [a for a in (inputs.actions or []) if _ts(a) is not None]
    judged = [a for a in actions if "in_goal" in a and a["in_goal"] is not None]
    if judged:
        drift = sum(1 for a in judged if not a["in_goal"])
        res.drift_actions = drift
        res.drift_fraction = drift / len(judged)
    res.steers = len(inputs.steers or [])
    steers_sorted = sorted((s for s in (inputs.steers or []) if _ts(s) is not None),
                           key=_ts)
    followed = 0
    for s in steers_sorted:
        t = _ts(s)
        # Followed: some in-goal action AFTER the steer (the very next
        # action the caller judged in_goal=True).
        later = [a for a in judged if _ts(a) > t]
        if later and any(a["in_goal"] for a in later):
            followed += 1
    res.steers_followed = followed
    return res
