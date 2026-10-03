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
import os
import re

from . import tool_trace
from dataclasses import dataclass, field
from pathlib import Path


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
    #: How many of this session's turns nobody could price, and how many
    #: ran where the provider charges nothing at all. Without them a cost
    #: of 0.00 says "free" for a session on a model with no published
    #: rate, which is the one thing it does not know.
    unpriced_turns: int = 0
    non_billing_turns: int = 0


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
    """The command of every bash-like tool call, as it was run.

    The trace stores a tool's input as the JSON it was called with --
    ``{"command": ..., "description": ...}`` -- so reading the first line
    of it put the JSON in the report where the command belonged. The
    payload is parsed, and raw text is still accepted for a recorder that
    stored the command itself.
    """
    out = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        tool = str(e.get("tool", ""))
        if "bash" not in tool.lower():
            continue
        raw = e.get("input")
        if raw is None:
            continue
        # The trace's own reading of an entry: the command is a field of
        # the recorded call, not the first line of its JSON.
        command = tool_trace.command_of(e)
        if not command and not tool_trace.call_args(e):
            command = str(raw)          # an older recorder stored the text
        first = command.strip().splitlines()
        if first and first[0].strip():
            out.append(first[0])
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


def _stored_record(session_id: str) -> dict:
    """The stored session record, or {} on any failure.

    Solo sessions never write agent_metrics rows (``record_turn`` is
    dashboard-only), but the record ``session_store.save_session`` keeps
    has carried model, token_usage, cost_usd and created/updated
    timestamps all along — the report just never looked at it (measured
    2026-09-19: a 13-minute session reported ``duration 1m 4s``,
    ``model -``, ``tokens 0``). This is the fallback, never the source
    of record: metrics rows keep precedence where they exist.
    """
    try:
        from delfin.agent import session_store

        data = session_store.load_session(session_id)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _numbers(session_id: str, entries: list[dict]) -> tuple:
    """model, started_at, ended_at, cost_usd, input/output tokens, and how
    many turns were unpriced or non-billing.

    Metrics rows are the source of record where they exist. Solo sessions
    write none, and for them the stored session record is the fallback:
    model, token_usage, and the created_at..updated_at span widened by
    run_elapsed_s. Duration is the UNION of trace and record — a trace
    spans the tool calls, not the session, and two calls at the end of a
    5-minute run must not report it as 30 ms.
    """
    rows = _turn_rows(session_id)
    record = _stored_record(session_id) if not rows else {}
    model = ""
    for row in rows:  # last row with a non-empty model wins
        if row.get("model"):
            model = str(row["model"])
    if not model:
        model = str(record.get("model", "") or "")
    cost = 0.0
    in_tok = 0
    out_tok = 0
    unpriced = 0
    non_billing = 0
    for row in rows:
        state = str(row.get("price_state") or "")
        if state == "unknown":
            unpriced += 1
        elif state == "non_billing":
            non_billing += 1
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
    if record:
        # The record's created/updated span — no metrics rows exist here,
        # so this is the session's own lifetime as the store saw it.
        try:
            created = float(record.get("created_at") or 0.0)
            updated = float(record.get("updated_at") or 0.0)
        except (TypeError, ValueError):
            created = updated = 0.0
        if created:
            rec_started = created
            # A solo session saves once at the end, so created == updated
            # and the span is zero. run_elapsed_s is the run clock the
            # engine persists for the resume budget — the real wall-clock
            # lifetime. It is a sum across resumes, never a span, so a
            # real updated_at still wins: never stretch a measured span.
            if updated and updated > created:
                rec_ended = updated
            else:
                try:
                    elapsed = float(record.get("run_elapsed_s") or 0.0)
                except (TypeError, ValueError):
                    elapsed = 0.0
                rec_ended = created + elapsed if elapsed > 0.0 else created
            started = min(started, rec_started) if started else rec_started
            ended = max(ended, rec_ended) if ended else rec_ended
    if not in_tok and not out_tok and record:
        usage = record.get("token_usage")
        if isinstance(usage, dict):
            try:
                in_tok = int(usage.get("input") or 0)
                out_tok = int(usage.get("output") or 0)
            except (TypeError, ValueError):
                pass
    return (model, started, ended, cost, in_tok, out_tok,
            unpriced, non_billing)


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

        (model, started, ended, cost, in_tok, out_tok,
         unpriced, non_billing) = _numbers(session_id, entries)
        report.model = model
        report.started_at = started
        report.ended_at = ended
        report.cost_usd = cost
        report.unpriced_turns = unpriced
        report.non_billing_turns = non_billing
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


def _fmt_cost(cost_usd: float, unpriced: int = 0, non_billing: int = 0) -> str:
    """The cost, or what its absence means.

    A zero is three different statements -- measured, charge-free, or
    never priced at all -- and printing "$0.00" for the third says the
    session was free when nobody knows what it cost.
    """
    try:
        value = float(cost_usd)
    except (TypeError, ValueError):
        return "-"
    if value > 0:
        return f"${value:.2f}"
    if int(unpriced or 0) > 0:
        return f"not measured ({int(unpriced)} turn(s) with no rate)"
    if int(non_billing or 0) > 0:
        return "no charge"
    return f"${value:.2f}"


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
    add(f"- **Cost:** {_fmt_cost(report.cost_usd, report.unpriced_turns, report.non_billing_turns)}")
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
        f"Cost {_fmt_cost(report.cost_usd, report.unpriced_turns, report.non_billing_turns)}"
        f" | tokens in {_fmt_int(report.input_tokens)}"
        f" out {_fmt_int(report.output_tokens)}",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# friction: lost time and stumbles, computed from the recorded trace only.
#
# The four wave-10 patterns -- lost time waiting on an approval, a denied
# action reworked in another form, one file read in many slices, and idle
# time after a job -- are all reconstructible from the per-call trace record
# (``tool_trace.record``: ts, tool, input, duration_ms, ok, error) with NO
# new logging inside the security code. Everything here reads the trace only;
# nothing reads security_events (an in-process ring buffer a CLI report
# always sees empty) and nothing depends on wall-clock time.
# --------------------------------------------------------------------------

_IDLE_MIN_S = 300.0          # consecutive-call gap this long is idle
_WAIT_MAX_S = 14400.0        # deny -> successful identical retry window (4h)
_REWORK_MAX_S = 600.0        # deny -> reworked-form window (a deliberate retry)
_READ_MIN_COUNT = 3          # same file sliced at least this often -> friction
_TOP_DEFAULT = 3

# Slice commands that read PART of a file -- the wave-10 pattern of one file
# read in many sed -n slices, each paying for an approval. grep deliberately
# excluded: grepping a big file is the intended lookup, not an inefficient
# partial read, so counting it would over-report ordinary use.
_EXTRACT_RE = re.compile(r"^\s*(?:sed\b|head\b|tail\b|awk\b)", re.I)


def _norm_cmd(command: str) -> str:
    """Normalise a bash command for same-action comparison (whitespace only)."""
    return re.sub(r"\s+", " ", (command or "").strip())


def _command(e: dict, is_bash: bool) -> str:
    """The command of a bash entry (via tool_trace.command_of), else ''."""
    if not is_bash:
        return ""
    try:
        return tool_trace.command_of(e) or ""
    except Exception:
        return ""


def _object_tokens(command: str) -> set[str]:
    """Path/file-like tokens a command touches, plus each one's basename.

    Used to decide whether a later call performs the same action in another
    form: the two commands must act on the SAME object (a shared path or file
    basename), so an unrelated later command is never a "rework".
    """
    tokens: set[str] = set()
    for t in (command or "").split():
        t = t.strip(",;|&")
        if not t or t.startswith("-") or t.startswith("'") or t.startswith('"'):
            continue
        if "/" in t or "." in t:
            tokens.add(t)
            tokens.add(t.rsplit("/", 1)[-1])  # basename, so paths/locals collide
    return tokens


def _is_denial(e: dict) -> bool:
    """A blocked call (ok is False, and the error is a refusal), not a timeout."""
    if e.get("ok", True):
        return False
    try:
        from .report_denials import _categorize
        return _categorize(str(e.get("error") or "")) is not None
    except Exception:
        return False


def _extract_target(command: str) -> str:
    """The file a slice command reads, or '' when none is identifiable."""
    toks = (command or "").strip().split()
    for t in reversed(toks):
        if t.startswith("-") or t.startswith("'") or t.startswith('"'):
            continue
        if "/" in t or "." in t:
            return t
    return ""


def _refusal_entries(entries: list[dict]) -> list[dict]:
    """The trace entries that are refusals, with the command of a bash one."""
    out = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        if not _is_denial(e):
            continue
        is_bash = "bash" in str(e.get("tool", "")).lower()
        cmd = _command(e, is_bash)
        out.append({"entry": e, "tool": str(e.get("tool") or ""),
                    "cmd": cmd, "ts": _ts(e)})
    return out


def _ts(e: dict) -> float:
    try:
        v = float(e.get("ts") or 0)
    except (TypeError, ValueError):
        return 0.0
    return v


def _idle_gaps(entries: list[dict], *, idle_min_s: float,
               wait_spans: list) -> list[dict]:
    """Consecutive-call gaps at or above the threshold, minus the spans already
    counted as approval waits (a wait and an idle gap over the same stretch
    would double-count the same lost time)."""
    gaps = []
    entries = [e for e in entries if isinstance(e, dict)]
    for i in range(1, len(entries)):
        prev, nxt = entries[i - 1], entries[i]
        p_ts, n_ts = _ts(prev), _ts(nxt)
        gap = n_ts - p_ts if (p_ts > 0 and n_ts > 0) else 0.0
        if gap < idle_min_s:
            continue
        if any(ws[0] - 1e-9 <= p_ts and n_ts <= ws[1] + 1e-9
               for ws in wait_spans):
            continue
        is_bash = "bash" in str(prev.get("tool", "")).lower()
        waiting = _command(prev, is_bash) or str(prev.get("tool") or "")
        gaps.append({"kind": "idle_gap", "magnitude": round(gap),
                     "unit": "s", "detail": str(prev.get("tool") or ""),
                     "command": waiting})
    return gaps


def compute_friction(
    entries: list[dict],
    *,
    idle_min_s: float = _IDLE_MIN_S,
    wait_max_s: float = _WAIT_MAX_S,
    rework_max_s: float = _REWORK_MAX_S,
    read_min_count: int = _READ_MIN_COUNT,
) -> dict:
    """Every friction event in one session's recorded trace.

    Returns a dict never raising on malformed entries, with one key per
    pattern, each a list of item dicts ``{kind, magnitude, unit, detail,
    command}``:
      approval_wait    a refusal followed later by the identical successful
                       call; magnitude = wait seconds, command = the denied cmd
      denial_reworked  a refusal followed later by a DIFFERENT-FORM bash call
                       on the same object (same path/basename); magnitude =
                       count, command = the reworked call
      repeated_reads   the same file sliced by sed/head/tail/awk at least
                       read_min_count times; magnitude = count, command = one
                       example slice
      idle_gap         a consecutive-call gap >= idle_min_s NOT already counted
                       as an approval wait; magnitude = gap seconds, command =
                       the call that was waiting
    """
    items = {"approval_wait": [], "denial_reworked": [],
             "repeated_reads": [], "idle_gap": []}
    waits = items["approval_wait"]
    wait_spans = []

    for d in _refusal_entries(entries):
        d_ts = d["ts"]
        matched = False
        for e in entries:
            if not isinstance(e, dict):
                continue
            if e is d["entry"]:
                continue
            if e.get("ok", True) is False:
                continue
            e_ts = _ts(e)
            if e_ts <= d_ts or e_ts - d_ts > wait_max_s:
                continue
            if str(e.get("tool") or "") != d["tool"]:
                continue
            is_bash = "bash" in str(e.get("tool", "")).lower()
            cmd = _command(e, is_bash)
            if _norm_cmd(cmd) == _norm_cmd(d["cmd"]):
                gap = round(e_ts - d_ts)
                waits.append({"kind": "approval_wait", "magnitude": gap,
                              "unit": "s", "detail": d["cmd"],
                              "command": d["cmd"]})
                wait_spans.append((d_ts, e_ts))
                matched = True
                break
        if matched:
            continue
        # no identical retry -> look for the same object reached another way
        if not d["cmd"]:
            continue
        denied_objs = _object_tokens(d["cmd"])
        for e in entries:
            if not isinstance(e, dict):
                continue
            if e is d["entry"] or e.get("ok", True) is False:
                continue
            if "bash" not in str(e.get("tool", "")).lower():
                continue
            e_ts = _ts(e)
            if e_ts <= d_ts or e_ts - d_ts > rework_max_s:
                continue
            cmd = _command(e, True)
            if not cmd or _norm_cmd(cmd) == _norm_cmd(d["cmd"]):
                continue
            if denied_objs & _object_tokens(cmd):
                items["denial_reworked"].append({
                    "kind": "denial_reworked", "magnitude": 1, "unit": "n",
                    "detail": cmd, "command": cmd})
                break

    # repeated partial reads, grouped per file
    by_file: dict[str, list[str]] = {}
    for e in entries:
        if not isinstance(e, dict) or "bash" not in str(e.get("tool", "")).lower():
            continue
        cmd = _command(e, True)
        if not _EXTRACT_RE.match(cmd):
            continue
        target = _extract_target(cmd)
        if target:
            by_file.setdefault(target, []).append(cmd)
    for path, cmds in by_file.items():
        if len(cmds) >= read_min_count:
            items["repeated_reads"].append({
                "kind": "repeated_reads", "magnitude": len(cmds),
                "unit": "reads", "detail": path, "command": cmds[0]})

    items["idle_gap"] = _idle_gaps(entries, idle_min_s=idle_min_s,
                                   wait_spans=wait_spans)
    return items


def friction_summary(
    entries: list[dict],
    *,
    top: int = _TOP_DEFAULT,
    idle_min_s: float = _IDLE_MIN_S,
    read_min_count: int = _READ_MIN_COUNT,
) -> list[dict]:
    """The top ``top`` frictions of a session, ranked by magnitude, with the
    command that caused each. Pure input for the REIBUNG report sessions write
    by hand. A clean trace returns [].
    """
    try:
        data = compute_friction(entries, idle_min_s=idle_min_s,
                                read_min_count=read_min_count)
    except Exception:
        return []
    ranked = []
    for kind, lst in data.items():
        ranked.extend(lst)
    ranked.sort(key=lambda r: (r.get("magnitude") or 0), reverse=True)
    return [dict(r) for r in ranked[:max(0, int(top))]]


# --------------------------------------------------------------------------
# shutdown hook: best-effort write of the Markdown report into
# ~/.delfin/session_reports/<safe-session-id>.md (per-session-file pattern
# shared with change_journal / pending_changes).
# --------------------------------------------------------------------------

def _report_dir() -> "Path":
    return Path.home() / ".delfin" / "session_reports"


def _own_it(path: "Path", mode: int) -> None:
    """Keep a report to the user it is about; never raises.

    This is the one file that collects every command of a session in one
    place, and it was the only one of these stores written with whatever
    the umask allowed -- tool_trace creates its file 0600, change_journal
    chmods 0600/0700. On a shared home the difference is who can read the
    session back.
    """
    try:
        os.chmod(path, mode)
    except OSError:
        pass


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
        path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
        _own_it(path.parent, 0o700)
        # Atomic-ish: write a temp file next to the target, then replace,
        # so a reader never sees a half-written report.
        tmp = path.with_suffix(".md.tmp")
        tmp.write_text(text, encoding="utf-8")
        # Narrow it BEFORE it takes the final name: a chmod afterwards
        # leaves a window in which the report is readable by others.
        _own_it(tmp, 0o600)
        tmp.replace(path)
        return path
    except Exception:
        # Best-effort by contract: a failed report must never break the
        # shutdown path that called us.
        return None
