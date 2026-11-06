"""Deterministic replay of a long-horizon session transcript -> metrics.

F2 mode (a): a scripted transcript is played through the REAL record
writers (turn_metrics.record, tool_trace.record,
compaction_log.record_compaction) into a store the caller redirects (a
scratch dir in tests under tmp_path; the run's own home in live mode). The
records are then read back and reduced by :mod:`metrics` — so the metric
proof in CI is over the same real writers a live run uses, with no model
call anywhere.

A transcript is a list of events, each one a record the matching real
writer accepts:

    {"kind": "turn",  ...turn_metrics.record kwargs}
    {"kind": "tool",  tool=..., input=..., output=..., ...}
    {"kind": "compaction", kind=..., messages_compacted=..., ...}

``replay`` returns the machine-readable JSON string of the metrics, plus
the parsed records (for an accept script to re-check numbers).

The order of events is the order of the session; the real writers stamp
wall-clock ts, so the compaction boundary lands where the transcript places
it. Two transcripts that differ only in which records carry a planted fact
late (recalling vs forgetting) are the red control: the first must score
higher on late recall.
"""
from __future__ import annotations

from typing import Any, Dict, List

from delfin.agent.longhorizon import metrics


def _emit(kind: str, ev: Dict[str, Any]) -> None:
    """Write ONE event through its real writer. Never raises (a garbage
    event is skipped, and the run it would have recorded is unmeasured,
    matching the writers' own never-raises contract)."""
    ev = dict(ev)
    ev.pop("kind", None)
    if kind == "turn":
        from delfin.agent import turn_metrics
        sess = ev.pop("session", "") or "replay"
        turn_metrics.record(sess, **ev)
    elif kind == "tool":
        from delfin.agent import tool_trace
        sess = ev.pop("session", "") or "replay"
        tool_trace.record(
            sess,
            tool=ev.get("tool", ""),
            tool_input=ev.get("input", ""),
            output=ev.get("output", ""),
            duration_ms=ev.get("duration_ms", 0),
            ok=bool(ev.get("ok", True)),
            error=ev.get("error", ""),
        )
    elif kind == "compaction":
        from delfin.agent import compaction_log
        sess = ev.pop("session", "") or "replay"
        # the compaction KIND is a distinct key in the transcript so it does
        # not collide with the event-dispatch ``kind``.
        ck = ev.pop("compaction_kind", "auto")
        ev["kind"] = ck
        compaction_log.record_compaction(sess, **ev)
    # any other kind is ignored: a transcript is only as useful as its
    # recognised record kinds.


def replay(
    transcript: List[Dict[str, Any]],
    *,
    session: str = "replay",
    planted_facts: List[Dict[str, Any]] | None = None,
    decisions: List[Dict[str, Any]] | None = None,
) -> Dict[str, Any]:
    """Replay ``transcript`` into records and reduce them to the metric
    payload. Return {"metrics": {...}, "records": {"turns": [...],
    "tools": [...], "compactions": [...]}}.

    Uses one ``session`` id for every event so the reads find them all.
    Pure apart from the writers it drives (which are the point)."""
    session = (session or "replay")
    for ev in transcript:
        kind = str(ev.get("kind", ""))
        ev = dict(ev)
        if "session" not in ev:
            ev["session"] = session
        _emit(kind, ev)

    from delfin.agent import compaction_log, tool_trace, turn_metrics
    turns = turn_metrics.read(session)
    tools = tool_trace.read(session)
    compactions = compaction_log.read_compactions(session)

    res = metrics.measure(turns, tools, compactions,
                          planted_facts=planted_facts, decisions=decisions)
    return {
        "metrics": _parse(metrics.to_json(res)),
        "records": {"turns": turns, "tools": tools,
                    "compactions": compactions},
    }


def _parse(payload: str) -> Dict[str, Any]:
    import json
    return json.loads(payload)
