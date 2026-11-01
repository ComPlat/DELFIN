"""Build an OpenTelemetry-style span tree for one DELFIN agent session.

Every diagnostic of a wave session used to mean reading tmux panes and three
JSONL stores by hand. One trace per run — spans for turns, model calls and
tool calls — is what LangSmith, the OpenAI Agents tracing and Claude Code's
OpenTelemetry export give.

Sources are ONLY the stores that already exist; nothing new is written:
- ``tool_trace``: one JSON line per tool call (ts, tool, duration_ms, ok,
  error, input, output) — ``tool_trace.record`` is the real writer.
- ``turn_metrics``: one JSON line per turn (ts, provider, model, total_ms,
  ttft_ms, input/output/cached tokens) — ``turn_metrics.record`` is the
  real writer.
- ``session_store.load_session``: title/provider/mode for the session span.

A store that is missing or empty gives a PARTIAL trace with a note, never an
error: the export exists to answer "what happened", and "nothing recorded"
is an answer, not a crash.

Redaction: tool input/output/error text goes through
``output_guard.scrub_secrets`` (exact held credentials first, then the shape
scan — the same path ``turn_metrics._clean_error`` uses) and outputs are cut
to 500 characters. No secret, key, token or approval code reaches a span.

The held-credential values are resolved ONCE per export (OT2), not once per
span field: ``spans()`` snapshots them at entry and passes them to
``scrub_secrets(exact_values=...)``. A credential that enters the store
during a single export is caught by the shape scan, not by the exact pass
of that export — a documented, tested snapshot choice, see the comment on
``_held_exact_values``.

No network. Sending to a collector is a separate, person-gated step
(:func:`write_otlp_json` writes a file; nothing here phones home).
"""

from __future__ import annotations

import json
from typing import Any

#: Tool output kept in a span, after redaction (brief says 500).
OUTPUT_CAP = 500
#: Tool input is capped tighter than output: it is a replay aid, not content.
INPUT_CAP = 200

#: Redact failures since the last spans() call — drained into the notes so
#: a broken redact path is VISIBLE ("note"), never a silent raw-text export.
_REDACT_FAILURES: list = []

#: Snapshot semantics (OT2, QS ruling): the held-credential values are
#: resolved ONCE per spans() call, at export start. A credential that enters
#: the credential store DURING a single export call is therefore not caught
#: by the exact pass in that export (the pattern scan still applies); the
#: next export picks it up. This is a documented, tested choice, not an
#: accidental window: the export is one atomic snapshot of "what DELFIN
#: holds", exactly like every other attribute of the trace is a snapshot of
#: what happened.
_HELD_EXACT: list | None = None
_HELD_EXACT_FAILED = False


def _reset_held_exact() -> None:
    """Clear the once-per-export snapshot (called at spans() entry)."""
    global _HELD_EXACT, _HELD_EXACT_FAILED
    _HELD_EXACT = None
    _HELD_EXACT_FAILED = False


def _held_exact_values() -> "list | tuple":
    """The held-credential values, resolved at most once per spans() call.

    Mirrors the per-call behaviour of scrub_secrets: a raising lookup
    degrades to the empty exact pass (the shape scan still runs); today
    that exception is swallowed inside scrub_secrets (output_guard.py).
    The failure is remembered so every field of this export degrades the
    same way, not just the first one.
    """
    global _HELD_EXACT, _HELD_EXACT_FAILED
    if _HELD_EXACT is not None:
        return _HELD_EXACT
    try:
        from .web_tools import _known_secret_values
        _HELD_EXACT = list(_known_secret_values())
    except Exception:
        _HELD_EXACT = []
        _HELD_EXACT_FAILED = True
    return _HELD_EXACT


def _redact(text: Any, cap: int) -> str:
    """Redact then cut, or DROP the content on a redact failure.

    QS ruling on F1: scrub_secrets unavailable or raising must never yield
    raw text in a span. The field content is dropped ("[redaction
    unavailable]") and the caller records a visible note — the failure is
    loud, the raw text never leaves.

    OT2: uses the once-per-export held-credential snapshot (see
    _held_exact_values) instead of a per-field lookup.
    """
    if not isinstance(text, str):
        try:
            text = json.dumps(text, ensure_ascii=False) if isinstance(
                text, (dict, list)) else str(text or "")
        except Exception:
            text = str(text or "")
    try:
        from .output_guard import scrub_secrets
        return scrub_secrets(text, exact_values=_held_exact_values())[:cap]
    except Exception:
        _REDACT_FAILURES.append(1)
        return "[redaction unavailable]"[:cap]


def _id(session: str, part: str) -> str:
    """Deterministic span id: stable exports, consistent parent ids."""
    h = 0
    for ch in f"{session}:{part}":
        h = (h * 131 + ord(ch)) & 0xFFFFFFFFFFFF
    return f"{h:012x}"


def _tool_spans(session: str, turns: list[dict], notes: list[str]) -> list[dict]:
    """One span per recorded tool call, parented by its containing turn.

    Correlation: a tool call is bound to the turn whose [ts, ts+total_ms]
    window contains the call's ts (tool_trace entries carry no turn index;
    `round` is not written by every recorder, ts is). Calls outside every
    window parent to the session span instead.
    """
    from . import tool_trace
    entries = tool_trace.read(session)
    if not entries:
        notes.append("tool_trace: no entries recorded for this session")
        return []
    spans = []
    for i, e in enumerate(entries):
        ts = float(e.get("ts", 0) or 0)
        parent = _id(session, "session")
        for j, t in enumerate(turns):
            if not ts:
                break
            t_ts = float(t.get("ts", 0) or 0)
            if t_ts <= ts <= t_ts + float(t.get("total_ms", 0) or 0) / 1000.0:
                parent = _id(session, f"turn{j}")
                break
        attrs = {
            "delfin.span_type": "tool_call",
            "tool.name": str(e.get("tool", "")),
            "tool.ok": bool(e.get("ok", True)),
            "duration_ms": int(e.get("duration_ms", 0) or 0),
            "tool.error": _redact(e.get("error", ""), 300),
        }
        if "input" in e:
            attrs["tool.input"] = _redact(e.get("input", ""), INPUT_CAP)
        if "output" in e:
            attrs["tool.output"] = _redact(e.get("output", ""), OUTPUT_CAP)
        spans.append({
            "trace_id": _id(session, "trace"),
            "span_id": _id(session, f"tool{i}"),
            "parent_span_id": parent,
            "name": "tool_call",
            "start_time_unix_nano": int(ts * 1e9),
            "end_time_unix_nano": int((ts + int(
                e.get("duration_ms", 0) or 0) / 1000.0) * 1e9),
            "attributes": attrs,
        })
    return spans


def _turn_spans(session: str, notes: list[str]) -> tuple[list[dict], list[dict]]:
    from . import turn_metrics
    entries = turn_metrics.read(session)
    if not entries:
        notes.append("turn_metrics: no entries recorded for this session")
        return [], []
    spans = []
    for i, t in enumerate(entries):
        ts = float(t.get("ts", 0) or 0)
        total_s = float(t.get("total_ms", 0) or 0) / 1000.0
        attrs = {
            "delfin.span_type": "turn",
            "provider": str(t.get("provider", "")),
            "model": str(t.get("model", "")),
            "role": str(t.get("role", "")),
            "total_ms": int(t.get("total_ms", 0) or 0),
            "ttft_ms": t.get("ttft_ms"),
            "output_chars": int(t.get("output_chars", 0) or 0),
            "stopped": bool(t.get("stopped", False)),
        }
        # Tokens: an omitted key means "the recorder had no field", a 0
        # means the endpoint reported none -- kept apart, as in turn_metrics.
        for k in ("input_tokens", "output_tokens", "cached_tokens"):
            if k in t:
                attrs[k] = int(t.get(k) or 0)
        if t.get("error"):
            attrs["error"] = _redact(t.get("error", ""), 300)
        spans.append({
            "trace_id": _id(session, "trace"),
            "span_id": _id(session, f"turn{i}"),
            "parent_span_id": _id(session, "session"),
            "name": "model_call",
            "start_time_unix_nano": int(ts * 1e9),
            "end_time_unix_nano": int((ts + total_s) * 1e9),
            "attributes": attrs,
        })
    return spans, entries


def spans(session: str, *, notes: bool = False):
    """Span tree for ``session``: session > turn > model_call / tool_call.

    Returns a list of span dicts (and, with ``notes=True``, a second list of
    human-readable notes about what was missing). Builds from the existing
    stores only; a missing store is a note, not an error.
    """
    note_list: list[str] = []
    _REDACT_FAILURES.clear()
    _reset_held_exact()
    session_span = {
        "trace_id": _id(session, "trace"),
        "span_id": _id(session, "session"),
        "parent_span_id": "",
        "name": "session",
        "start_time_unix_nano": 0,
        "end_time_unix_nano": 0,
        "attributes": {"delfin.span_type": "session"},
    }
    try:
        from . import session_store
        rec = session_store.load_session(session)
        if rec:
            for key, attr in (("title", "session.title"),
                              ("provider", "session.provider"),
                              ("model", "session.model"),
                              ("mode", "session.mode")):
                if rec.get(key):
                    session_span["attributes"][attr] = str(rec[key])
        else:
            note_list.append("session_store: no saved session record")
    except Exception:
        note_list.append("session_store: not readable here")

    turn_spans, turn_entries = _turn_spans(session, note_list)
    tool_spans = _tool_spans(session, turn_entries, note_list)

    session_start = min(
        [s["start_time_unix_nano"] for s in turn_spans + tool_spans
         if s["start_time_unix_nano"]] or [0])
    session_end = max(
        [s["end_time_unix_nano"] for s in turn_spans + tool_spans] or [0])
    session_span["start_time_unix_nano"] = session_start
    session_span["end_time_unix_nano"] = max(session_end, session_start)
    session_span["attributes"]["tool_calls"] = len(tool_spans)
    session_span["attributes"]["turns"] = len(turn_spans)

    # Parent-fix: tool calls recorded outside any turn window parent to the
    # session span directly, so every span has a resolvable parent.
    all_spans = [session_span] + turn_spans + tool_spans
    ids = {s["span_id"] for s in all_spans}
    for s in all_spans:
        if s["parent_span_id"] and s["parent_span_id"] not in ids:
            s["parent_span_id"] = _id(session, "session")
    out = list(all_spans)
    if _REDACT_FAILURES:
        note_list.append(
            f"redaction: {len(_REDACT_FAILURES)} field(s) DROPPED — the "
            "redact path failed; raw text was never exported")
        _REDACT_FAILURES.clear()
    return (out, note_list) if notes else out


# ---------------------------------------------------------------------------
# OTLP/JSON export and text timeline
# ---------------------------------------------------------------------------

_ATTR_KEY_MAP = {
    "delfin.span_type": "delfin.span_type",
}


def _otlp_attr(key: str, value) -> dict:
    """One OTLP KeyValue. Bytes/arrays are not produced here."""
    if isinstance(value, bool):
        return {"key": key, "value": {"boolValue": value}}
    if isinstance(value, int):
        return {"key": key, "value": {"intValue": str(value)}}
    if isinstance(value, float):
        return {"key": key, "value": {"doubleValue": value}}
    return {"key": key, "value": {"stringValue": str(value)}}


def _otlp_span(s: dict) -> dict:
    """One span dict -> OTLP/JSON span (Jaeger / Grafana Tempo import shape)."""
    out = {
        "traceId": s["trace_id"],
        "spanId": s["span_id"],
        "name": s["name"],
        "kind": "SPAN_KIND_INTERNAL",
        "startTimeUnixNano": str(s["start_time_unix_nano"]),
        "endTimeUnixNano": str(s["end_time_unix_nano"]),
        "attributes": [_otlp_attr(k, v)
                       for k, v in s["attributes"].items()],
    }
    parent = s.get("parent_span_id") or ""
    if parent:
        out["parentSpanId"] = parent
    return out


def write_otlp_json(session: str, path) -> "str | None":
    """Write the span tree as an OTLP/JSON file a collector can import.

    No network: the file is written locally. Sending it to a collector is a
    separate, person-gated decision (``telemetry.otlp_endpoint`` in the
    dashboard, through the egress proxy rules) — this function never phones
    home, whatever that setting says.

    Returns the written path, or None when writing failed.
    """
    from pathlib import Path
    import time as _time
    tree, note_list = spans(session, notes=True)
    doc = {
        "resourceSpans": [{
            "resource": {
                "attributes": [
                    _otlp_attr("service.name", "delfin-agent"),
                    _otlp_attr("delfin.session", session),
                ],
            },
            "scopeSpans": [{
                "scope": {"name": "delfin.agent.otel_export",
                          "version": "1"},
                "spans": [_otlp_span(s) for s in tree],
            }],
        }],
    }
    p = Path(path)
    try:
        p.write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                     encoding="utf-8")
        return str(p)
    except OSError:
        return None


def timeline(session: str) -> str:
    """One line per span, indented, with durations; the longest 3 marked.

    Text timeline for ``delfin-agent trace <session>`` (cli.py). Notes are
    rendered so a terminal user sees a redaction failure too (s4, phase-3
    review): a dropped field must never look like an ordinary timeline.
    """
    tree, note_list = spans(session, notes=True)
    if not tree:
        return "(no spans — nothing recorded for this session)"
    durations = sorted((s["end_time_unix_nano"] - s["start_time_unix_nano"]
                        for s in tree), reverse=True)
    top3 = set(durations[:3])

    def _depth(s: dict) -> int:
        return {"session": 0, "model_call": 1, "turn": 1,
                "tool_call": 2}.get(s["name"], 1)

    lines = []
    for s in tree:
        dur_ms = (s["end_time_unix_nano"] - s["start_time_unix_nano"]) / 1e6
        ind = "  " * _depth(s)
        attrs = s["attributes"]
        label = attrs.get("tool.name") or attrs.get("model") or ""
        mark = "*" if (s["end_time_unix_nano"]
                       - s["start_time_unix_nano"]) in top3 else " "
        lines.append(
            f"{mark} {ind}{s['name']}"
            + (f" {label}" if label else "")
            + f"  ({dur_ms:.1f} ms)")
    for n in note_list:
        lines.append(f"! {n}")
    return "\n".join(lines)

