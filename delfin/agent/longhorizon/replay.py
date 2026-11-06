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
    {"kind": "compaction", compaction_kind=..., messages_compacted=...,}

``replay`` returns {"metrics": {...}, "records": {...}, "skipped_events":
[...]}. An event the writers cannot accept (a field a writer does not have,
a typo'd ``kind``) is SKIPPED and reported in ``skipped_events`` — it never
crashes the run and it never vanishes silently, matching the reducer's "a
value the records cannot attest is reported, not guessed" standard.

The order of events is the order of the session; the real writers stamp
wall-clock ts, so the compaction boundary lands where the transcript places
it. Two transcripts that differ only in which records carry a planted fact
late (recalling vs forgetting) are the red control: the first must score
higher on late recall.
"""
from __future__ import annotations

import inspect
from typing import Any, Dict, List

from delfin.agent.longhorizon import metrics


def _writer_kwargs(writer, ev: Dict[str, Any], key_map: Dict[str, str]) -> Dict[str, Any]:
    """Forward only the kwargs the writer's ``record`` signature accepts,
    after renaming transcript keys via ``key_map``. A field a real writer
    keeps (e.g. ``round`` on tool_trace.record) survives; a key the writer
    does not have is dropped here AND reported as a skipped event, so a
    transcript that names a field the writer never records is diagnosable
    instead of silently lost."""
    try:
        params = inspect.signature(writer).parameters
    except (TypeError, ValueError):
        params = {}
    out: Dict[str, Any] = {}
    for key, value in ev.items():
        if key in key_map:
            key = key_map[key]
        if key in params and params[key].kind in (
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
                inspect.Parameter.KEYWORD_ONLY):
            out[key] = value
    return out


def _emit(kind: str, ev: Dict[str, Any]):
    """Write ONE event through its real writer; return None on success or a
    short reason string when the event cannot be recorded."""
    if kind == "turn":
        from delfin.agent import turn_metrics
        sess = ev.pop("session", "") or "replay"
        kw = _writer_kwargs(turn_metrics.record, ev, {})
        turn_metrics.record(sess, **kw)
    elif kind == "tool":
        from delfin.agent import tool_trace
        sess = ev.pop("session", "") or "replay"
        kw = _writer_kwargs(tool_trace.record, ev, {"input": "tool_input"})
        tool_trace.record(sess, **kw)
    elif kind == "compaction":
        from delfin.agent import compaction_log
        sess = ev.pop("session", "") or "replay"
        # the compaction KIND is a distinct transcript key so it does not
        # collide with the event-dispatch ``kind``.
        ck = ev.pop("compaction_kind", "auto")
        ev["kind"] = ck
        kw = _writer_kwargs(compaction_log.record_compaction, ev, {})
        compaction_log.record_compaction(sess, **kw)
    else:
        return f"unrecognised kind {kind!r}"
    return None


def replay(
    transcript: List[Dict[str, Any]],
    *,
    session: str = "replay",
    planted_facts: List[Dict[str, Any]] | None = None,
    decisions: List[Dict[str, Any]] | None = None,
) -> Dict[str, Any]:
    """Replay ``transcript`` into records and reduce them to the metric
    payload. Return {"metrics": {...}, "records": {"turns": [...],
    "tools": [...], "compactions": [...]}, "skipped_events": [...]}.

    Each event is written through its real writer inside a narrow
    try/except: a malformed event is SKIPPED and reported, never a crash
    and never silent — so transcript-authoring bugs surface in
    ``skipped_events`` instead of aborting or quietly zeroing a metric.
    Uses one ``session`` id for every event so the reads find them all."""
    session = (session or "replay")
    skipped: List[Dict[str, Any]] = []
    for original in transcript:
        kind = str(original.get("kind", ""))
        ev = dict(original)
        if "session" not in ev:
            ev["session"] = session
        try:
            reason = _emit(kind, ev)
        except (TypeError, ValueError) as exc:
            reason = f"{type(exc).__name__}: {exc}"
        if reason:
            skipped.append({"event": original, "reason": reason})

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
        "skipped_events": skipped,
    }


def _parse(payload: str) -> Dict[str, Any]:
    import json
    return json.loads(payload)
