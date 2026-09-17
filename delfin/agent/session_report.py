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
        elif isinstance(extra, dict) and extra.get("deleted"):
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


# --------------------------------------------------------------------------
# rendering (Session B). Pure functions: no I/O, no colour, no imports of
# the collector's sources — they only read the dataclass.
# --------------------------------------------------------------------------

def _fmt_duration(started_at: float, ended_at: float) -> str:
    """Human duration between two unix timestamps; '-' when unknown."""
    try:
        seconds = float(ended_at) - float(started_at)
    except (TypeError, ValueError):
        return "-"
    if seconds < 0 or started_at <= 0.0 or ended_at <= 0.0:
        return "-"
    m, s = divmod(int(round(seconds)), 60)
    h, m = divmod(m, 60)
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def _fmt_time(ts: float) -> str:
    """UTC HH:MM:SS for a unix timestamp; '-' when unknown (0.0/negative)."""
    try:
        ts = float(ts)
    except (TypeError, ValueError):
        return "-"
    if ts <= 0.0:
        return "-"
    import datetime as _dt

    return _dt.datetime.fromtimestamp(ts, tz=_dt.timezone.utc).strftime(
        "%Y-%m-%d %H:%M:%S UTC")


def _fmt_cost(cost_usd: float) -> str:
    try:
        return f"${float(cost_usd):.2f}"
    except (TypeError, ValueError):
        return "-"


def _fmt_int(value) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return "-"


def render_markdown(report: SessionReport) -> str:
    """Render a SessionReport as a clean Markdown document.

    Sections: header (session id, model, duration), Tool Calls, Files
    Changed, Commands Run, Tests Run, Denials, Cost & Tokens.
    """
    lines: list[str] = []
    add = lines.append

    add(f"# Session Report: {report.session_id or '(unknown session)'}")
    add("")
    duration = _fmt_duration(report.started_at, report.ended_at)
    add(f"- **Model:** {report.model or '-'}")
    add(f"- **Started:** {_fmt_time(report.started_at)}")
    add(f"- **Ended:** {_fmt_time(report.ended_at)}")
    add(f"- **Duration:** {duration}")
    add("")

    add("## Tool Calls")
    add("")
    if report.tool_calls:
        add("| Tool | Calls | OK | Failed |")
        add("|---|---:|---:|---:|")
        for row in report.tool_calls:
            add("| {name} | {count} | {ok} | {failed} |".format(
                name=str(row.get("name", "") or "-").replace("|", "\\|"),
                count=_fmt_int(row.get("count", 0)),
                ok=_fmt_int(row.get("ok", 0)),
                failed=_fmt_int(row.get("failed", 0)),
            ))
    else:
        add("_(none)_")
    add("")

    add("## Files Changed")
    add("")
    if report.files_changed:
        add("| File | Change |")
        add("|---|---|")
        for row in report.files_changed:
            add("| {path} | {change} |".format(
                path=str(row.get("path", "") or "-").replace("|", "\\|"),
                change=str(row.get("change", "") or "-").replace("|", "\\|"),
            ))
    else:
        add("_(none)_")
    add("")

    add("## Commands Run")
    add("")
    if report.commands_run:
        add("```")
        for cmd in report.commands_run:
            add(str(cmd))
        add("```")
    else:
        add("_(none)_")
    add("")

    add("## Tests Run")
    add("")
    if report.tests_run:
        add("| Target | Status | Passed | Failed |")
        add("|---|---|---:|---:|")
        for row in report.tests_run:
            add("| {target} | {status} | {passed} | {failed} |".format(
                target=str(row.get("target", "") or "-").replace("|", "\\|"),
                status=str(row.get("status", "") or "-").replace("|", "\\|"),
                passed=_fmt_int(row.get("passed", 0)),
                failed=_fmt_int(row.get("failed", 0)),
            ))
    else:
        add("_(none)_")
    add("")

    add("## Denials")
    add("")
    if report.denials:
        for row in report.denials:
            add("- **{kind}**: {detail}".format(
                kind=str(row.get("kind", "") or "-"),
                detail=str(row.get("detail", "") or "-"),
            ))
    else:
        add("_(none)_")
    add("")

    add("## Cost & Tokens")
    add("")
    add(f"- **Cost:** {_fmt_cost(report.cost_usd)}")
    add(f"- **Input tokens:** {_fmt_int(report.input_tokens)}")
    add(f"- **Output tokens:** {_fmt_int(report.output_tokens)}")
    add("")

    return "\n".join(lines)


def render_terminal(report: SessionReport) -> str:
    """Render a SessionReport as a compact plain-text summary.

    Plain ASCII, no colour codes — safe for any terminal, log or pager.
    """
    total_calls = 0
    total_failed = 0
    for row in report.tool_calls or []:
        try:
            total_calls += int(row.get("count", 0) or 0)
            total_failed += int(row.get("failed", 0) or 0)
        except (TypeError, ValueError):
            pass
    n_files = len(report.files_changed or [])
    n_tests = len(report.tests_run or [])
    n_denials = len(report.denials or [])
    tests_passed = sum(
        1 for row in report.tests_run or []
        if str(row.get("status", "")) == "passed")
    tests_failed = n_tests - tests_passed

    lines = [
        f"Session {report.session_id or '(unknown)'}"
        f" | model {report.model or '-'}"
        f" | duration {_fmt_duration(report.started_at, report.ended_at)}",
        f"Tool calls: {total_calls} ({total_failed} failed)"
        f" | files changed: {n_files}"
        f" | tests: {tests_passed} passed, {tests_failed} failed"
        f" | denials: {n_denials}",
        f"Cost {_fmt_cost(report.cost_usd)}"
        f" | tokens in {_fmt_int(report.input_tokens)}"
        f" out {_fmt_int(report.output_tokens)}",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# shutdown hook: best-effort write of the Markdown report into
# ~/.delfin/session_reports/<safe-session-id>.md (per-session-file pattern
# shared with change_journal / pending_changes).
# --------------------------------------------------------------------------

def _report_dir() -> "Path":
    from pathlib import Path

    return Path.home() / ".delfin" / "session_reports"


def _safe_session_id(session_id: str) -> str:
    """Same sanitization as change_journal._safe_session_id (no traversal)."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", str(session_id or "") or "session")[:40]


def write_session_report(session_id: str) -> "Path | None":
    """Collect and write the Markdown session report; never raises.

    Returns the written path, or None when nothing could be done (empty
    session id, or every step failed — including the collection itself,
    which by contract also never raises).
    """
    try:
        sid = str(session_id or "").strip()
        if not sid:
            return None
        report = collect_session_report(sid)
        text = render_markdown(report)
        path = _report_dir() / (_safe_session_id(sid) + ".md")
        path.parent.mkdir(parents=True, exist_ok=True)
        # Atomic-ish: write a temp file next to the target, then replace,
        # so a reader never sees a half-written report.
        tmp = path.with_suffix(".md.tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(path)
        return path
    except Exception:
        # Best-effort by contract: a failed report must never break the
        # shutdown path that called us.
        return None
