"""One report over every agent session that ran since a point in time.

After each supervised round the operator used to tally by hand what the
sessions did: tool calls and failures, dialogues, denials, calls over
five minutes, commits, tokens, ttft, endpoint errors — one session at a
time, plus a sum. DELFIN already holds every one of those numbers:

* ``tool_trace``          — per session: every call, ok/error, duration_ms
* ``turn_metrics``        — per session: per-turn ttft, tokens, cache, errors
* ``report_denials``      — the refusal categories (its ``_categorize``)

This module only re-reads those stores and aggregates; it never counts
anything a source does not already hold. Where a number is not available
(a session without a turn log) the field is ``None`` and the rendering
says ``n/a`` — a zero would read as measured.

Read-only: nothing here writes to any store. The tests run against
fabricated traces under ``tmp_path``; the first real run is the
operator's (CLI: ``delfin-agent report --since <time> [--name PREFIX]``).
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from statistics import median
from typing import Any

from . import tool_trace, turn_metrics
from .report_denials import _NOT_A_REFUSAL, _categorize

#: Calls at or above this duration count as "over five minutes".
FIVE_MIN_MS = 5 * 60 * 1000

_DIALOGUE_TOOLS = re.compile(r"ask_user|ask_user_question", re.I)
_COMMIT = re.compile(r"\bgit\s+commit\b")


# --------------------------------------------------------------------- since

_REL_SINCE = re.compile(r"^(\d+(?:\.\d+)?)\s*([smhd])$", re.I)
_REL_UNITS = {"s": 1.0, "m": 60.0, "h": 3600.0, "d": 86_400.0}


def parse_since(spec: str) -> float:
    """A ``--since`` value -> unix cutoff timestamp.

    Accepts a relative span (``90m``, ``12h``, ``3d``) or an ISO-ish
    timestamp (``2026-09-18``, ``2026-09-18T00:00:00``). Raises ValueError
    on anything else — a silently-reinterpreted cutoff would report the
    wrong sessions with a straight face.
    """
    spec = (spec or "").strip()
    m = _REL_SINCE.match(spec)
    if m:
        return time.time() - float(m.group(1)) * _REL_UNITS[m.group(2).lower()]
    for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            import calendar
            st = time.strptime(spec, fmt)
            return float(calendar.timegm(st))
        except ValueError:
            continue
    raise ValueError(f"cannot parse --since value: {spec!r}")


# ------------------------------------------------------------------ audit

def _audit_facts(audit_path: "Path | None") -> tuple[dict, dict]:
    """(session -> label, session -> dialogs) from the audit log.

    The label is the session's -n name where a dialog record carries it,
    else the name of its workspace directory (every audit record names
    the workspace). Dialogs are the ``event: dialog`` records the
    confirmation broker writes for each question a person answered.
    """
    labels: dict[str, str] = {}
    dialogs: dict[str, int] = {}
    try:
        import json as _json
        from . import audit_log as _audit
        path = Path(audit_path) if audit_path else _audit._default_log_path()
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    rec = _json.loads(line)
                except ValueError:
                    continue
                sid = str(rec.get("session_id") or "")
                if not sid:
                    continue
                if rec.get("event") == "dialog":
                    dialogs[sid] = dialogs.get(sid, 0) + 1
                    if rec.get("session_key"):
                        labels[sid] = str(rec["session_key"])
                if sid not in labels and rec.get("workspace"):
                    labels[sid] = Path(str(rec["workspace"])).name
    except Exception:
        pass
    return labels, dialogs


# ------------------------------------------------------------------ sessions

def _session_ids(trace_root: Path, name: str, since: float,
                 labels: "dict | None" = None) -> list[tuple[str, float]]:
    """(session, last activity ts) for sessions with activity >= since,
    newest first and filtered by name prefix. A session whose whole trace
    predates the cutoff is not part of the round."""
    out: list[tuple[str, float]] = []
    for sid in tool_trace.sessions(root=trace_root):
        entries = tool_trace.read(sid, root=trace_root)
        ts_list = [float(e.get("ts") or 0.0) for e in entries
                   if isinstance(e, dict)]
        last = max(ts_list, default=0.0)
        if last < since:
            continue
        label = (labels or {}).get(sid, "")
        if name and not (sid.startswith(name) or label.startswith(name)):
            continue
        out.append((sid, last))
    out.sort(key=lambda kv: kv[1], reverse=True)
    return out


def _commit_count(entries: list[dict]) -> int:
    """Successful git commits, read from the traces: a bash call whose
    command contains ``git commit`` and whose ``ok`` is true. The command
    comes from the trace's own reading (``command_of``), never re-parsed
    here."""
    n = 0
    for e in entries:
        if e.get("ok", True) is False:
            continue
        cmd = tool_trace.command_of(e)
        if cmd and _COMMIT.search(cmd):
            n += 1
    return n


def _top_tools(entries: list[dict], limit: int = 5) -> list[dict]:
    per: dict[str, dict] = {}
    for e in entries:
        name = str(e.get("tool") or "unknown")
        row = per.setdefault(name, {"name": name, "count": 0, "failed": 0})
        row["count"] += 1
        if e.get("ok", True) is False:
            row["failed"] += 1
    return sorted(per.values(),
                  key=lambda r: (-r["count"], r["name"]))[:limit]


def _trace_stats(entries: list[dict]) -> dict[str, Any]:
    total = len(entries)
    failed = sum(1 for e in entries if e.get("ok", True) is False)
    over = [e for e in entries
            if float(e.get("duration_ms") or 0) >= FIVE_MIN_MS]
    denials: dict[str, int] = {}
    dialogues = 0
    for e in entries:
        err = str(e.get("error") or "")
        if (e.get("ok", True) is False and err
                and not _NOT_A_REFUSAL.search(err)):
            cat = _categorize(err)
            if cat:
                denials[cat] = denials.get(cat, 0) + 1
        if _DIALOGUE_TOOLS.search(str(e.get("tool") or "")):
            dialogues += 1
    return {
        "tool_calls": total,
        "tool_calls_failed": failed,
        "failed_pct": (100.0 * failed / total) if total else 0.0,
        "top_tools": _top_tools(entries),
        "dialogues": dialogues,
        "denials": sum(denials.values()),
        "denial_reasons": denials,
        "calls_over_5min": len(over),
        "sum_over_5min_min": (
            sum(float(e.get("duration_ms") or 0) for e in over) / 60_000),
        "commits": _commit_count(entries),
    }


def _turn_stats(session_id: str, since: float) -> dict[str, Any]:
    """Token/cache/ttft/endpoint-error numbers from the session's turn log.

    None (rendered ``n/a``) when the session has no turn log or none in
    the window — not zero, which would read as measured.
    """
    try:
        rows = [r for r in turn_metrics.read(session_id)
                if float(r.get("ts") or 0.0) >= since]
    except Exception:
        rows = []
    if not rows:
        return {"input_tokens": None, "output_tokens": None,
                "cached_pct": None, "ttft_median_ms": None,
                "ttft_max_ms": None, "endpoint_errors": None,
                "turns": 0}
    in_tok = sum(int(r.get("input_tokens") or 0) for r in rows)
    out_tok = sum(int(r.get("output_tokens") or 0) for r in rows)
    cached = sum(int(r.get("cached_tokens") or 0) for r in rows)
    ttfts = [int(r["ttft_ms"]) for r in rows
             if r.get("ttft_ms") is not None]
    errors = sum(1 for r in rows
                 if str(r.get("error") or "").strip())
    return {
        "input_tokens": in_tok,
        "output_tokens": out_tok,
        "cached_pct": (100.0 * cached / in_tok) if in_tok else None,
        "ttft_median_ms": int(median(ttfts)) if ttfts else None,
        "ttft_max_ms": max(ttfts) if ttfts else None,
        "endpoint_errors": errors,
        "turns": len(rows),
    }


# -------------------------------------------------------------------- public

def collect(*, since_s: float, trace_root: "Path | str | None" = None,
            name: str = "", audit_path: "Path | str | None" = None) -> dict:
    """Aggregate every session with activity since ``since_s`` (unix ts).

    Returns ``{"since_s": …, "sessions": [per-session dicts, newest
    first], "totals": {…}}``. Never raises on missing stores; a number a
    source cannot give is None.
    """
    root = Path(trace_root) if trace_root else tool_trace._DIR
    labels, dialogs = _audit_facts(Path(audit_path) if audit_path else None)
    sessions: list[dict] = []
    try:
        for sid, last in _session_ids(root, name, since_s, labels):
            entries = tool_trace.read(sid, root=root)
            row = {"session_id": sid, "label": labels.get(sid, ""),
                   "last_activity": last}
            row.update(_trace_stats(entries))
            # A person answering a permission dialog is the operator
            # load the round is about; ask_user calls in the trace are
            # only the questions the model put itself.
            row["dialogues"] = row["dialogues"] + dialogs.get(sid, 0)
            row.update(_turn_stats(sid, since_s))
            sessions.append(row)
    except Exception:
        sessions = sessions or []

    totals: dict[str, Any] = {"sessions": len(sessions)}
    for key in ("tool_calls", "tool_calls_failed", "dialogues", "denials",
                "calls_over_5min", "commits", "endpoint_errors",
                "input_tokens", "output_tokens"):
        vals = [s.get(key) for s in sessions if s.get(key) is not None]
        totals[key] = sum(vals) if vals else (0 if key in (
            "tool_calls", "tool_calls_failed", "dialogues", "denials",
            "calls_over_5min", "commits") else None)
    calls = totals.get("tool_calls") or 0
    totals["failed_pct"] = (100.0 * totals["tool_calls_failed"] / calls
                            if calls else 0.0)
    if totals.get("input_tokens"):
        cached_sum = 0.0
        for s in sessions:
            cp = s.get("cached_pct")
            if cp is not None and s.get("input_tokens"):
                cached_sum += float(cp) * s["input_tokens"] / 100.0
        totals["cached_pct"] = (100.0 * cached_sum / totals["input_tokens"]
                                if totals["input_tokens"] else None)
    else:
        totals["cached_pct"] = None
    return {"since_s": since_s, "sessions": sessions, "totals": totals}


def _n(v: Any, suffix: str = "") -> str:
    if v is None:
        return "n/a"
    return f"{v:.1f}{suffix}" if isinstance(v, float) else f"{v}{suffix}"


def render_text(data: dict) -> str:
    """Terminal rendering: one block per session, then the sum."""
    lines: list[str] = []
    totals = data.get("totals") or {}
    sessions = data.get("sessions") or []
    lines.append(
        f"Round report — {len(sessions)} session(s) since "
        f"{time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(data.get('since_s') or 0))}"
    )
    lines.append("")
    if not sessions:
        lines.append("No sessions in the window.")
        return "\n".join(lines)
    for s in sessions:
        label = f" {s['label']}" if s.get("label") else ""
        lines.append(f"[{s['session_id'][:12]}]{label}")
        lines.append(
            f"  tools: {_n(s.get('tool_calls'))} calls, "
            f"{_n(s.get('tool_calls_failed'))} failed "
            f"({s.get('failed_pct', 0.0):.1f}%)"
        )
        top = ", ".join(f"{t['name']}×{t['count']}"
                        for t in (s.get("top_tools") or []))
        lines.append(f"  top tools: {top or 'n/a'}")
        lines.append(
            f"  dialogues: {s.get('dialogues', 0)}, "
            f"denials: {s.get('denials', 0)} "
            f"({', '.join(f'{k}:{v}' for k, v in sorted((s.get('denial_reasons') or {}).items())) or 'none'})"
        )
        lines.append(
            f"  >5min calls: {s.get('calls_over_5min', 0)} "
            f"(sum {_n(s.get('sum_over_5min_min'))} min), "
            f"commits: {s.get('commits', 0)}"
        )
        lines.append(
            f"  tokens: in {_n(s.get('input_tokens'))} / "
            f"out {_n(s.get('output_tokens'))}, "
            f"cached {_n(s.get('cached_pct'), '%')}, "
            f"ttft median {_n(s.get('ttft_median_ms'), ' ms')} "
            f"max {_n(s.get('ttft_max_ms'), ' ms')}"
        )
        lines.append(
            f"  turns: {s.get('turns', 0)}, "
            f"endpoint errors: {_n(s.get('endpoint_errors'))}"
        )
        lines.append("")
    t = totals
    lines.append("TOTAL")
    lines.append(
        f"  sessions: {t.get('sessions')}, "
        f"tools: {t.get('tool_calls')} calls, "
        f"{t.get('tool_calls_failed')} failed "
        f"({t.get('failed_pct', 0.0):.1f}%)"
    )
    lines.append(
        f"  dialogues: {t.get('dialogues')}, denials: {t.get('denials')}, "
        f">5min calls: {t.get('calls_over_5min')}, "
        f"commits: {t.get('commits')}"
    )
    lines.append(
        f"  tokens: in {_n(t.get('input_tokens'))} / "
        f"out {_n(t.get('output_tokens'))}, "
        f"cached {_n(t.get('cached_pct'), '%')}, "
        f"endpoint errors: {_n(t.get('endpoint_errors'))}"
    )
    return "\n".join(lines)
