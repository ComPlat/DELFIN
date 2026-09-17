"""Session report: aggregate one agent session's activity from real sources.

Sources (all optional, none may raise):
- delfin.agent.tool_trace.read            -> tool call entries
- delfin.agent.change_journal.list_changes -> file change records
- delfin.agent.security_events.recent     -> security/denial events
- delfin.agent.agent_metrics.read_turns   -> per-turn metrics (model, cost, tokens)

Shared contract: the ``SessionReport`` dataclass. Two other modules consume
this field list; keep the field names stable (see tests/test_session_report.py).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field


@dataclass
class SessionReport:
    """Everything known about one agent session, from real local sources.

    Shared contract — do NOT rename fields without notifying downstream
    consumers (Sessions B and C).
    """

    session_id: str = ""
    model: str = ""
    started_at: float = 0.0          # unix ts of first tool call (0.0 if unknown)
    ended_at: float = 0.0            # unix ts of last tool call (0.0 if unknown)
    tool_calls: list[dict] = field(default_factory=list)
    # each: {"name": str, "count": int, "ok": int, "failed": int}
    files_changed: list[dict] = field(default_factory=list)
    # each: {"path": str, "change": "created" | "modified" | "deleted"}
    commands_run: list[str] = field(default_factory=list)
    tests_run: list[dict] = field(default_factory=list)
    # each: {"target": str, "status": str, "passed": int, "failed": int}
    denials: list[dict] = field(default_factory=list)
    # each: {"kind": str, "detail": str}
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0


# --------------------------------------------------------------------------
# helpers: each wraps ONE source; never raises, never returns None.
# --------------------------------------------------------------------------

_TEST_TOOL_HINTS = ("run_tests", "pytest")


def _tool_entries(session_id: str) -> list[dict]:
    """Tool trace entries for this session; [] on any failure."""
    try:
        from delfin.agent import tool_trace

        return tool_trace.read(session_id) or []
    except Exception:
        return []


def _change_records(session_id: str) -> list[dict]:
    try:
        from delfin.agent import change_journal

        return change_journal.list_changes(session_id) or []
    except Exception:
        return []


def _security_events() -> list:
    try:
        from delfin.agent import security_events

        return security_events.recent() or []
    except Exception:
        return []


def _turn_rows(session_id: str) -> list[dict]:
    """Turn metric rows belonging to this session; [] on any failure."""
    try:
        from delfin.agent import agent_metrics

        rows = agent_metrics.read_turns() or []
    except Exception:
        return []
    out = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        if str(row.get("session_id", "")) == str(session_id):
            out.append(row)
    return out


# --------------------------------------------------------------------------
# field builders
# --------------------------------------------------------------------------

def _tool_call_summary(entries: list[dict]) -> list[dict]:
    """Aggregate trace entries per tool name -> {name, count, ok, failed}."""
    per: dict[str, dict] = {}
    for e in entries:
        if not isinstance(e, dict):
            continue
        name = str(e.get("tool", "") or "unknown")
        row = per.setdefault(name, {"name": name, "count": 0, "ok": 0, "failed": 0})
        row["count"] += 1
        if e.get("ok", True):
            row["ok"] += 1
        else:
            row["failed"] += 1
    return list(per.values())


def _files_changed(records: list[dict]) -> list[dict]:
    """Map journal records to {path, change} in journal order (oldest first)."""
    out = []
    seen: dict[str, dict] = {}
    for rec in records:
        if not isinstance(rec, dict):
            continue
        path = str(rec.get("path", "") or "")
        if not path:
            continue
        extra = rec.get("extra") or {}
        if isinstance(rec.get("created"), bool) and rec.get("created"):
            change = "created"
        elif rec.get("deleted") or (isinstance(extra, dict) and extra.get("deleted")):
            # _write_record merges `extra` into the record top level
            # (change_journal.py:478), so `deleted` sits there, not nested.
            change = "deleted"
        else:
            change = "modified"
        if path in seen:
            seen[path]["change"] = change  # latest change wins
        else:
            entry = {"path": path, "change": change}
            seen[path] = entry
            out.append(entry)
    return out


def _commands_run(entries: list[dict]) -> list[str]:
    """First line of every bash-like tool input (the auditable command)."""
    out = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        tool = str(e.get("tool", ""))
        if "bash" not in tool.lower():
            continue
        raw = e.get("input")
        if isinstance(raw, str):
            cmd = raw.strip().splitlines()
            if cmd:
                out.append(cmd[0])
        elif raw is not None:
            out.append(str(raw))
    return out


_PASS_FAIL = re.compile(r"(\d+) passed", re.I)
_FAILED = re.compile(r"(\d+) failed", re.I)


def _tests_run(entries: list[dict]) -> list[dict]:
    """Best-effort: pytest-style tool calls -> {target, status, passed, failed}.

    Counts are parsed from the recorded output; when the output is not
    parseable the status is "unknown" and counts stay 0 — a missing parse is
    a limitation of the trace, not an error.
    """
    out = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        tool = str(e.get("tool", ""))
        if not any(h in tool.lower() for h in _TEST_TOOL_HINTS):
            continue
        raw = e.get("input")
        if isinstance(raw, str):
            try:
                payload = json.loads(raw)
            except (json.JSONDecodeError, ValueError):
                payload = {}
            args = payload.get("pytest_args")
            target = payload.get("target") or payload.get("nodeid")
            if not target and isinstance(args, list) and args:
                target = str(args[0])
            target = str(target or "")
        elif raw is not None:
            target = str(raw)
        else:
            target = ""
        output = e.get("output")
        out_text = output if isinstance(output, str) else ""
        m_pass = _PASS_FAIL.search(out_text)
        m_fail = _FAILED.search(out_text)
        passed = int(m_pass.group(1)) if m_pass else 0
        failed = int(m_fail.group(1)) if m_fail else 0
        if e.get("ok", True) is False:
            status = "failed"
        elif m_pass or m_fail:
            status = "passed" if failed == 0 else "failed"
        else:
            status = "unknown"
        out.append({"target": target, "status": status, "passed": passed, "failed": failed})
    return out


def _denials() -> list[dict]:
    """Blocked security events -> {kind, detail}.

    Limitation: SecurityEvent carries no session_id, so these are
    process-global, not session-scoped.
    """
    out = []
    for ev in _security_events():
        try:
            if not getattr(ev, "blocked", True):
                continue
            out.append({"kind": str(ev.kind), "detail": str(ev.detail)})
        except Exception:
            continue
    return out


def _numbers(session_id: str, entries: list[dict]) -> tuple[str, float, float, float, int, int]:
    """model, started_at, ended_at, cost_usd, input_tokens, output_tokens."""
    rows = _turn_rows(session_id)
    model = ""
    for row in rows:  # last row with a non-empty model wins
        if row.get("model"):
            model = str(row["model"])
    cost = 0.0
    in_tok = 0
    out_tok = 0
    for row in rows:
        try:
            cost += float(row.get("cost_usd") or 0.0)
            in_tok += int(row.get("input_tokens") or 0)
            out_tok += int(row.get("output_tokens") or 0)
        except (TypeError, ValueError):
            continue
    timestamps = [float(e["ts"]) for e in entries
                  if isinstance(e, dict) and isinstance(e.get("ts"), (int, float))]
    started = min(timestamps) if timestamps else 0.0
    ended = max(timestamps) if timestamps else 0.0
    if rows and not timestamps:
        try:
            turn_ts = [float(r.get("ts") or 0.0) for r in rows]
            turn_ts = [t for t in turn_ts if t]
            if turn_ts:
                started, ended = min(turn_ts), max(turn_ts)
        except (TypeError, ValueError):
            pass
    return model, started, ended, cost, in_tok, out_tok


# --------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------

def collect_session_report(session_id: str) -> SessionReport:
    """Build a SessionReport for one session from the real local sources.

    Never raises: a missing or misbehaving source contributes an empty
    value for its fields.
    """
    report = SessionReport(session_id=str(session_id))
    try:
        entries = _tool_entries(session_id)
        records = _change_records(session_id)

        report.tool_calls = _tool_call_summary(entries)
        report.files_changed = _files_changed(records)
        report.commands_run = _commands_run(entries)
        report.tests_run = _tests_run(entries)
        report.denials = _denials()

        (model, started, ended,
         cost, in_tok, out_tok) = _numbers(session_id, entries)
        report.model = model
        report.started_at = started
        report.ended_at = ended
        report.cost_usd = cost
        report.input_tokens = in_tok
        report.output_tokens = out_tok
    except Exception:
        # Contract: never raise. Anything already filled stays filled.
        pass
    return report
