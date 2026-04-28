"""Agent Activity tab: outcome history + cost/performance overview.

Reads the JSONL written by ``delfin.agent.outcome_tracker`` and renders a
timeline + aggregate stats. Read-only — no actions, just visibility.
"""
from __future__ import annotations

import html as _html
from collections import Counter
from datetime import datetime
from pathlib import Path

import ipywidgets as widgets

from delfin.agent.outcome_tracker import CycleOutcome, load_outcomes
from delfin.dashboard.claude_hooks import discover_hooks, render_hooks_html
from delfin.dashboard.claude_settings import (
    discover_settings,
    render_settings_html,
)
from delfin.dashboard.git_worktrees import (
    list_worktrees,
    render_worktrees_html,
)
from delfin.dashboard.schedules import (
    discover_schedules,
    render_schedules_html,
)
from delfin.dashboard.mcp_overview import (
    discover_mcp_servers,
    render_overview_html as render_mcp_overview_html,
)


_VERDICT_COLORS = {
    "PASS": "#10b981",
    "FAIL": "#ef4444",
    "PARTIAL": "#f59e0b",
    "": "#9ca3af",
}


# ---------------------------------------------------------------------------
# A1 — Live state pane
# ---------------------------------------------------------------------------

_LIVE_STATUS_BADGES = {
    "streaming": ("#3b82f6", "Streaming"),
    "idle":      ("#10b981", "Idle"),
    "stopped":   ("#9ca3af", "Stopped"),
    "missing":   ("#9ca3af", "Not started"),
}


def _format_live_state(snapshot: dict) -> str:
    """Render the current engine + jobs snapshot as a compact pane.

    ``snapshot`` keys (all optional, missing → block omitted):
    - ``engine_status``: ``"streaming"`` | ``"idle"`` | ``"stopped"`` | ``"missing"``
    - ``current_mode``: e.g. ``"dashboard"``
    - ``current_model``: e.g. ``"sonnet"``
    - ``current_provider``: e.g. ``"claude"``
    - ``session_cost_usd``: float
    - ``session_turns``: int
    - ``active_jobs``: dict[status_str, count]
    - ``active_session_id``: str (truncated for display)
    - ``perm_profile``: e.g. ``"repo_free"``

    Always returns HTML — even an empty snapshot renders a small
    "agent not started" placeholder so the user knows the pane is alive.
    """
    status = (snapshot.get("engine_status") or "missing").lower()
    badge_color, badge_label = _LIVE_STATUS_BADGES.get(status, _LIVE_STATUS_BADGES["missing"])

    status_html = (
        f'<span style="background:{badge_color};color:#fff;padding:2px 10px;'
        f'border-radius:12px;font-size:11px;font-weight:700;'
        f'letter-spacing:0.3px;">{_html.escape(badge_label)}</span>'
    )

    parts = []
    mode = snapshot.get("current_mode") or ""
    model = snapshot.get("current_model") or ""
    provider = snapshot.get("current_provider") or ""
    if mode or model:
        provider_str = f"{_html.escape(provider)}/{_html.escape(model)}" if provider else _html.escape(model)
        mode_str = _html.escape(mode) if mode else "—"
        parts.append(
            f'<span style="color:#6b7280;font-size:11px;">'
            f'mode <strong style="color:#111827;">{mode_str}</strong> · '
            f'<strong style="color:#111827;">{provider_str}</strong>'
            f'</span>'
        )

    cost = float(snapshot.get("session_cost_usd") or 0.0)
    turns = int(snapshot.get("session_turns") or 0)
    if cost > 0 or turns > 0:
        parts.append(
            f'<span style="color:#6b7280;font-size:11px;">'
            f'session <strong style="color:#6366f1;">{_fmt_cost(cost)}</strong>'
            f'{f" · {turns} turn(s)" if turns else ""}'
            f'</span>'
        )

    jobs = snapshot.get("active_jobs") or {}
    if jobs:
        job_summary = ", ".join(
            f'{_html.escape(str(s))}:<strong style="color:#111827;">{int(c)}</strong>'
            for s, c in jobs.items() if int(c or 0) > 0
        )
        if job_summary:
            parts.append(
                f'<span style="color:#6b7280;font-size:11px;">jobs {job_summary}</span>'
            )

    perm = snapshot.get("perm_profile") or ""
    if perm:
        parts.append(
            f'<span style="color:#6b7280;font-size:11px;">'
            f'perms <strong style="color:#111827;">{_html.escape(perm)}</strong>'
            f'</span>'
        )

    sid = snapshot.get("active_session_id") or ""
    if sid:
        sid_short = sid[:8] + "…" if len(sid) > 8 else sid
        parts.append(
            f'<span style="color:#9ca3af;font-size:10px;font-family:monospace;">'
            f'sid {_html.escape(sid_short)}</span>'
        )

    body = (
        '<span style="color:#9ca3af;font-size:11px;font-style:italic;">'
        'No active engine. Start a conversation in the Agent tab.</span>'
        if not parts and status == "missing"
        else " · ".join(parts) or
             '<span style="color:#9ca3af;font-size:11px;">'
             '(engine present but idle)</span>'
    )

    return (
        '<div style="display:flex;align-items:center;gap:10px;padding:10px 14px;'
        'background:#f9fafb;border:1px solid #e5e7eb;border-radius:8px;'
        'margin-bottom:12px;">'
        f'{status_html}<div style="flex:1;display:flex;flex-wrap:wrap;gap:14px;">'
        f'{body}</div></div>'
    )


def _fmt_ts(ts: str) -> str:
    """Render an ISO timestamp as ``YYYY-MM-DD HH:MM``."""
    if not ts:
        return ""
    try:
        dt = datetime.fromisoformat(ts)
        return dt.strftime("%Y-%m-%d %H:%M")
    except ValueError:
        return ts[:16]


def _fmt_cost(usd: float) -> str:
    if usd < 0.01:
        return f"${usd * 100:.2f}¢"
    return f"${usd:.3f}"


def _fmt_duration(seconds: float) -> str:
    if seconds <= 0:
        return ""
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}m"
    return f"{seconds / 3600:.1f}h"


def _truncate(text: str, n: int) -> str:
    if len(text) <= n:
        return text
    return text[: n - 1] + "…"


def _filter_outcomes(
    outcomes: list[CycleOutcome],
    *,
    provider: str = "",
    mode: str = "",
    verdict: str = "",
    task_class: str = "",
) -> list[CycleOutcome]:
    """Apply UI filters to the loaded outcomes."""
    out = outcomes
    if provider:
        out = [o for o in out if o.provider == provider]
    if mode:
        out = [o for o in out if o.mode == mode]
    if verdict:
        out = [o for o in out if o.verdict == verdict]
    if task_class:
        out = [o for o in out if o.task_class == task_class]
    return out


def _outcome_cost_delta(outcome: CycleOutcome) -> float:
    """Return the per-cycle cost Δ for an outcome.

    Newer outcomes (A6+) carry an explicit ``cost_usd_delta``. Legacy
    entries written before that field don't, but we still record the
    cumulative ``cost_usd``; the caller can recover Δ by diffing against
    the previous outcome's cumulative total when both belong to the same
    session (heuristic: timestamps within 4 h, same provider+mode).
    """
    delta = float(getattr(outcome, "cost_usd_delta", 0.0) or 0.0)
    if delta > 0:
        return delta
    # No delta on legacy entry — caller should fall back via _outcomes_total_delta.
    return 0.0


def _outcomes_total_delta(outcomes: list[CycleOutcome]) -> float:
    """Sum honest per-cycle costs across a list of outcomes.

    Uses ``cost_usd_delta`` when present; otherwise falls back to a
    "diff against previous same-session entry" heuristic, and finally
    uses ``cost_usd`` itself for the very first session entry.
    """
    if not outcomes:
        return 0.0
    total = 0.0
    # Group by (provider, mode) timeline so legacy session-Δ stays sane
    last_session_total: dict[tuple[str, str], float] = {}
    last_ts: dict[tuple[str, str], str] = {}
    SESSION_WINDOW_S = 4 * 3600
    for o in outcomes:
        delta = _outcome_cost_delta(o)
        if delta > 0:
            total += delta
            continue
        # Legacy entry — compute against previous same-session entry.
        key = (o.provider or "", o.mode or "")
        prev = last_session_total.get(key, 0.0)
        prev_ts = last_ts.get(key, "")
        same_session = bool(prev_ts) and _ts_close(prev_ts, o.timestamp, SESSION_WINDOW_S)
        if same_session and float(o.cost_usd) >= prev:
            total += float(o.cost_usd) - prev
        else:
            # First entry of a (provider, mode) timeline or far-apart turn:
            # treat its cumulative cost as its own delta — better than 0.
            total += float(o.cost_usd)
        last_session_total[key] = float(o.cost_usd)
        last_ts[key] = o.timestamp
    return total


def _ts_close(ts_a: str, ts_b: str, window_s: float) -> bool:
    """True iff two ISO timestamps are within ``window_s`` of each other."""
    try:
        a = datetime.fromisoformat(ts_a)
        b = datetime.fromisoformat(ts_b)
    except (ValueError, TypeError):
        return False
    return abs((b - a).total_seconds()) <= window_s


def _aggregate_costs_by_period(
    outcomes: list[CycleOutcome],
    *,
    now: datetime | None = None,
) -> dict:
    """Bucket honest per-cycle deltas into today / week / month / all.

    Boundaries are calendar-based in the local timezone:
    - "today" = same calendar date as ``now``
    - "week"  = within the same ISO week (Mon-Sun)
    - "month" = within the same calendar month
    - "all"   = sum across every outcome we have

    Outcomes with unparseable timestamps fall through into "all" but
    are NOT counted toward today/week/month (we can't place them).
    """
    if not outcomes:
        return {"today": 0.0, "week": 0.0, "month": 0.0, "all": 0.0}

    now = now or datetime.now()
    today = now.date()
    iso_year, iso_week, _ = today.isocalendar()
    cur_month = (today.year, today.month)

    today_total = 0.0
    week_total = 0.0
    month_total = 0.0
    all_total = _outcomes_total_delta(outcomes)

    # We need the per-outcome delta in the same way _outcomes_total_delta
    # computes it, including the legacy heuristic. Reuse it row-by-row.
    last_cum: dict[tuple[str, str], float] = {}
    last_ts: dict[tuple[str, str], str] = {}
    SESSION_WINDOW_S = 4 * 3600

    for o in outcomes:
        delta = _outcome_cost_delta(o)
        if delta <= 0:
            key = (o.provider or "", o.mode or "")
            prev = last_cum.get(key, 0.0)
            prev_ts = last_ts.get(key, "")
            same_session = bool(prev_ts) and _ts_close(prev_ts, o.timestamp, SESSION_WINDOW_S)
            if same_session and float(o.cost_usd) >= prev:
                delta = float(o.cost_usd) - prev
            else:
                delta = float(o.cost_usd)
            last_cum[key] = float(o.cost_usd)
            last_ts[key] = o.timestamp
        try:
            ts = datetime.fromisoformat(o.timestamp) if o.timestamp else None
        except (ValueError, TypeError):
            ts = None
        if ts is None:
            continue
        d = ts.date()
        if d == today:
            today_total += delta
        if d.isocalendar()[:2] == (iso_year, iso_week):
            week_total += delta
        if (d.year, d.month) == cur_month:
            month_total += delta

    return {
        "today": today_total,
        "week":  week_total,
        "month": month_total,
        "all":   all_total,
    }


def _aggregate_costs_by_mode(
    outcomes: list[CycleOutcome],
) -> dict[str, float]:
    """Sum honest per-cycle deltas grouped by mode."""
    if not outcomes:
        return {}
    by_mode: dict[str, float] = {}
    last_cum: dict[tuple[str, str], float] = {}
    last_ts: dict[tuple[str, str], str] = {}
    SESSION_WINDOW_S = 4 * 3600

    for o in outcomes:
        delta = _outcome_cost_delta(o)
        if delta <= 0:
            key = (o.provider or "", o.mode or "")
            prev = last_cum.get(key, 0.0)
            prev_ts = last_ts.get(key, "")
            same_session = bool(prev_ts) and _ts_close(prev_ts, o.timestamp, SESSION_WINDOW_S)
            if same_session and float(o.cost_usd) >= prev:
                delta = float(o.cost_usd) - prev
            else:
                delta = float(o.cost_usd)
            last_cum[key] = float(o.cost_usd)
            last_ts[key] = o.timestamp
        mode = o.mode or "unknown"
        by_mode[mode] = by_mode.get(mode, 0.0) + delta
    return by_mode


def _render_cost_insights(period: dict, by_mode: dict[str, float]) -> str:
    """Render the period cards + per-mode bar list."""
    period_cards = []
    for label, key, color in [
        ("Today",     "today", "#3b82f6"),
        ("This week", "week",  "#10b981"),
        ("This month","month", "#6366f1"),
        ("All time",  "all",   "#8b5cf6"),
    ]:
        value = float(period.get(key, 0.0) or 0.0)
        period_cards.append(
            f'<div style="flex:1;min-width:120px;padding:10px 12px;'
            f'background:#fff;border:1px solid #e5e7eb;border-radius:8px;'
            f'border-left:3px solid {color};">'
            f'<div style="font-size:10px;font-weight:600;color:#6b7280;'
            f'text-transform:uppercase;letter-spacing:0.4px;">{_html.escape(label)}</div>'
            f'<div style="font-size:16px;font-weight:700;color:#111827;'
            f'margin-top:2px;">{_fmt_cost(value)}</div></div>'
        )
    cards_html = (
        '<div style="display:flex;flex-wrap:wrap;gap:8px;margin-bottom:12px;">'
        + "".join(period_cards) + '</div>'
    )

    if not by_mode:
        return cards_html
    max_v = max(by_mode.values()) or 1.0
    rows = []
    for mode, v in sorted(by_mode.items(), key=lambda kv: -kv[1]):
        pct = max(2.0, (v / max_v) * 100.0)
        rows.append(
            f'<div style="display:flex;align-items:center;gap:8px;'
            f'margin:3px 0;font-size:11px;">'
            f'<div style="width:96px;color:#374151;">{_html.escape(mode)}</div>'
            f'<div style="flex:1;background:#f3f4f6;border-radius:4px;'
            f'height:14px;position:relative;overflow:hidden;">'
            f'<div style="position:absolute;left:0;top:0;height:100%;'
            f'width:{pct:.1f}%;background:#6366f1;"></div></div>'
            f'<div style="width:80px;text-align:right;color:#6366f1;'
            f'font-weight:600;">{_fmt_cost(v)}</div></div>'
        )
    bars_html = (
        '<div style="background:#f9fafb;border:1px solid #e5e7eb;'
        'border-radius:6px;padding:8px 12px;">'
        '<div style="font-size:10px;font-weight:600;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;margin-bottom:4px;">'
        'Cost by mode</div>'
        + "".join(rows) + '</div>'
    )
    return cards_html + bars_html


def _aggregate_stats(outcomes: list[CycleOutcome]) -> dict:
    """Compute summary statistics (totals, success rate, avg cost/duration).

    ``total_cost`` is the honest sum of per-cycle Δ (uses cost_usd_delta
    when present, falls back to session-aware diff for legacy entries).
    """
    n = len(outcomes)
    if n == 0:
        return {
            "n": 0, "passes": 0, "fails": 0, "partials": 0,
            "success_rate": 0.0, "total_cost": 0.0,
            "avg_cost": 0.0, "avg_duration": 0.0,
        }
    passes = sum(1 for o in outcomes if o.verdict == "PASS")
    fails = sum(1 for o in outcomes if o.verdict == "FAIL")
    partials = sum(1 for o in outcomes if o.verdict == "PARTIAL")
    rated = passes + fails + partials
    success_rate = (passes / rated) if rated else 0.0
    total_cost = _outcomes_total_delta(outcomes)
    durations = [o.duration_s for o in outcomes if o.duration_s > 0]
    return {
        "n": n,
        "passes": passes,
        "fails": fails,
        "partials": partials,
        "success_rate": success_rate,
        "total_cost": total_cost,
        "avg_cost": total_cost / n,
        "avg_duration": (sum(durations) / len(durations)) if durations else 0.0,
    }


def _render_summary(stats: dict) -> str:
    """Render the aggregate-stats card row as HTML."""
    if stats["n"] == 0:
        return (
            '<div style="padding:12px;background:#f9fafb;border-radius:6px;'
            'color:#6b7280;font-size:12px;">'
            'Keine Outcome-Einträge gefunden. Sobald der Agent Cycles abschließt, '
            'erscheinen sie hier.</div>'
        )
    cards = [
        ("Total runs", str(stats["n"]), "#3b82f6"),
        ("Pass", str(stats["passes"]), "#10b981"),
        ("Fail", str(stats["fails"]), "#ef4444"),
        ("Partial", str(stats["partials"]), "#f59e0b"),
        ("Success rate", f"{stats['success_rate'] * 100:.0f}%", "#3b82f6"),
        ("Total cost", _fmt_cost(stats["total_cost"]), "#6366f1"),
        ("Avg cost/run", _fmt_cost(stats["avg_cost"]), "#8b5cf6"),
        ("Avg duration", _fmt_duration(stats["avg_duration"]), "#14b8a6"),
    ]
    cells = []
    for label, value, color in cards:
        cells.append(
            f'<div style="flex:1;min-width:120px;padding:10px 12px;'
            f'background:#fff;border:1px solid #e5e7eb;border-radius:8px;'
            f'border-left:3px solid {color};">'
            f'<div style="font-size:10px;font-weight:600;color:#6b7280;'
            f'text-transform:uppercase;letter-spacing:0.4px;">{label}</div>'
            f'<div style="font-size:16px;font-weight:700;color:#111827;'
            f'margin-top:2px;">{value}</div></div>'
        )
    return (
        '<div style="display:flex;flex-wrap:wrap;gap:8px;margin-bottom:12px;">'
        + "".join(cells) + '</div>'
    )


def _outcome_has_drilldown(o: CycleOutcome) -> bool:
    """Decide whether an outcome warrants the ▸ expand toggle.

    True when there's content the summary row truncates: long task
    text, denied commands, error_type, or any retries. Otherwise the
    summary line tells the whole story and we don't add the toggle.
    """
    if (o.task or "") and len(o.task) > 80:
        return True
    if o.denied_commands:
        return True
    if o.error_type:
        return True
    if int(o.retries or 0) > 0:
        return True
    return False


def _render_outcome_drilldown(o: CycleOutcome) -> str:
    """Render the expanded detail block under a row's ▸ toggle."""
    parts: list[str] = []
    if o.task:
        parts.append(
            f'<div style="margin:6px 0 2px 0;font-size:10px;color:#6b7280;'
            f'text-transform:uppercase;letter-spacing:0.4px;">Task (full)</div>'
            f'<div style="font-size:11px;color:#1f2937;'
            f'white-space:pre-wrap;word-break:break-word;">'
            f'{_html.escape(o.task)}</div>'
        )
    if o.denied_commands:
        cmd_html = "".join(
            f'<li style="margin:2px 0;">{_html.escape(c)}</li>'
            for c in o.denied_commands
        )
        parts.append(
            f'<div style="margin:8px 0 2px 0;font-size:10px;color:#ef4444;'
            f'text-transform:uppercase;letter-spacing:0.4px;">'
            f'Denied commands ({len(o.denied_commands)})</div>'
            f'<ul style="margin:0;padding-left:18px;font-size:11px;'
            f'color:#7f1d1d;font-family:monospace;">{cmd_html}</ul>'
        )
    if o.error_type:
        parts.append(
            f'<div style="margin:8px 0 2px 0;font-size:10px;color:#dc2626;'
            f'text-transform:uppercase;letter-spacing:0.4px;">Error</div>'
            f'<div style="font-size:11px;color:#7f1d1d;font-family:monospace;">'
            f'{_html.escape(o.error_type)}</div>'
        )
    if int(o.retries or 0) > 0:
        parts.append(
            f'<div style="margin:8px 0 2px 0;font-size:10px;color:#f59e0b;'
            f'text-transform:uppercase;letter-spacing:0.4px;">Retries</div>'
            f'<div style="font-size:11px;color:#92400e;">{int(o.retries)}</div>'
        )
    if o.timestamp:
        parts.append(
            f'<div style="margin:8px 0 2px 0;font-size:10px;color:#6b7280;'
            f'text-transform:uppercase;letter-spacing:0.4px;">Timestamp</div>'
            f'<div style="font-size:11px;color:#374151;font-family:monospace;">'
            f'{_html.escape(o.timestamp)}</div>'
        )
    return (
        '<div style="background:#f9fafb;border-top:1px solid #e5e7eb;'
        'padding:8px 14px;">' + "".join(parts) + '</div>'
    )


def _render_timeline(outcomes: list[CycleOutcome], limit: int = 100) -> str:
    """Render the recent-cycles timeline.

    Layout: a vertical list of <details> elements, one per outcome.
    The summary line is a flex row that mimics the previous table
    columns; clicking it expands the drilldown panel underneath.

    Outcomes with nothing additional to show (short task, no denied
    commands, no error, no retries) skip the toggle entirely so the
    row stays clean.
    """
    if not outcomes:
        return '<div style="padding:12px;color:#6b7280;font-size:12px;">' \
               'Keine Einträge passen zum Filter.</div>'
    rows = []
    for o in reversed(outcomes[-limit:]):
        verdict_color = _VERDICT_COLORS.get(o.verdict, "#9ca3af")
        verdict_html = (
            f'<span style="background:{verdict_color};color:#fff;'
            f'padding:1px 8px;border-radius:10px;font-size:10px;'
            f'font-weight:700;">{_html.escape(o.verdict or "—")}</span>'
        )
        denied = ""
        if o.denied_commands:
            denied = (
                f'<span style="margin-left:6px;color:#ef4444;font-size:10px;">'
                f'⛔ {len(o.denied_commands)}</span>'
            )
        err = (
            f'<span style="margin-left:6px;color:#dc2626;font-size:10px;">'
            f'✕ {_html.escape(o.error_type)}</span>' if o.error_type else ""
        )
        retries = (
            f'<span style="margin-left:6px;color:#f59e0b;font-size:10px;">'
            f'↻ {o.retries}</span>' if o.retries else ""
        )

        # Common cell styling for the flex-row summary
        cell = "padding:6px 8px;font-size:11px;"
        time_cell = (
            f'<span style="{cell}font-family:monospace;color:#6b7280;'
            f'min-width:120px;">{_html.escape(_fmt_ts(o.timestamp))}</span>'
        )
        verdict_cell = (
            f'<span style="{cell}min-width:100px;">'
            f'{verdict_html}{denied}{err}{retries}</span>'
        )
        mode_cell = (
            f'<span style="{cell}color:#374151;min-width:80px;">'
            f'{_html.escape(o.mode or "—")}</span>'
        )
        provider_cell = (
            f'<span style="{cell}color:#374151;min-width:140px;">'
            f'{_html.escape(o.provider or "—")}/{_html.escape(o.model or "")}'
            f'</span>'
        )
        class_cell = (
            f'<span style="{cell}color:#374151;min-width:90px;">'
            f'{_html.escape(o.task_class or "—")}</span>'
        )
        cost_cell = (
            f'<span style="{cell}color:#6366f1;text-align:right;'
            f'min-width:70px;">{_fmt_cost(o.cost_usd)}</span>'
        )
        dur_cell = (
            f'<span style="{cell}color:#14b8a6;text-align:right;'
            f'min-width:60px;">{_fmt_duration(o.duration_s)}</span>'
        )
        task_cell = (
            f'<span style="{cell}color:#1f2937;flex:1;min-width:0;'
            f'overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">'
            f'{_html.escape(_truncate(o.task, 80))}</span>'
        )

        summary_row = (
            f'<summary style="display:flex;align-items:center;gap:4px;'
            f'border-bottom:1px solid #f3f4f6;cursor:pointer;'
            f'list-style:none;">'
            f'{time_cell}{verdict_cell}{mode_cell}{provider_cell}'
            f'{class_cell}{cost_cell}{dur_cell}{task_cell}'
            f'</summary>'
        )

        if _outcome_has_drilldown(o):
            row = (
                f'<details style="border-bottom:1px solid #f3f4f6;">'
                f'{summary_row}'
                f'{_render_outcome_drilldown(o)}'
                f'</details>'
            )
        else:
            # No drilldown content — render a plain row that doesn't
            # advertise an expand toggle (cursor stays default).
            row = (
                f'<div style="display:flex;align-items:center;gap:4px;'
                f'border-bottom:1px solid #f3f4f6;">'
                f'{time_cell}{verdict_cell}{mode_cell}{provider_cell}'
                f'{class_cell}{cost_cell}{dur_cell}{task_cell}'
                f'</div>'
            )
        rows.append(row)

    header = (
        '<div style="display:flex;align-items:center;gap:4px;'
        'background:#f9fafb;border-bottom:2px solid #e5e7eb;'
        'font-size:10px;color:#6b7280;text-transform:uppercase;'
        'letter-spacing:0.4px;font-weight:600;">'
        '<span style="padding:6px 8px;min-width:120px;">Time</span>'
        '<span style="padding:6px 8px;min-width:100px;">Verdict</span>'
        '<span style="padding:6px 8px;min-width:80px;">Mode</span>'
        '<span style="padding:6px 8px;min-width:140px;">Provider</span>'
        '<span style="padding:6px 8px;min-width:90px;">Class</span>'
        '<span style="padding:6px 8px;min-width:70px;text-align:right;">Cost</span>'
        '<span style="padding:6px 8px;min-width:60px;text-align:right;">Duration</span>'
        '<span style="padding:6px 8px;flex:1;">Task</span>'
        '</div>'
    )
    return (
        '<div style="max-height:520px;overflow-y:auto;border:1px solid #e5e7eb;'
        'border-radius:6px;background:#fff;font-family:-apple-system,'
        'BlinkMacSystemFont,sans-serif;">'
        f'{header}{"".join(rows)}'
        '</div>'
    )


def _options_with_blank(values: list[str]) -> list[tuple[str, str]]:
    """Build ipywidgets dropdown options with a leading blank entry."""
    seen = sorted({v for v in values if v})
    return [("(alle)", "")] + [(v, v) for v in seen]


def create_tab(ctx, history_path: Path | None = None):
    """Build the Agent Activity tab.

    Parameters
    ----------
    ctx : DashboardContext
        Standard dashboard context (unused, kept for tab API consistency).
    history_path : Path, optional
        Override the outcome history file (mainly for tests).
    """
    # ---- widgets -------------------------------------------------------
    title = widgets.HTML(
        value='<h3 style="margin:0;color:#111827;">Agent Activity</h3>'
              '<p style="margin:4px 0 12px 0;color:#6b7280;font-size:12px;">'
              'Outcome history aller Agent-Cycles (Lese-Modus).</p>'
    )
    provider_dd = widgets.Dropdown(
        options=[("(alle)", "")], value="", description="Provider:",
        layout=widgets.Layout(width="200px"),
        style={"description_width": "70px"},
    )
    mode_dd = widgets.Dropdown(
        options=[("(alle)", "")], value="", description="Mode:",
        layout=widgets.Layout(width="200px"),
        style={"description_width": "70px"},
    )
    verdict_dd = widgets.Dropdown(
        options=[("(alle)", ""), ("PASS", "PASS"), ("FAIL", "FAIL"),
                 ("PARTIAL", "PARTIAL")],
        value="", description="Verdict:",
        layout=widgets.Layout(width="200px"),
        style={"description_width": "70px"},
    )
    class_dd = widgets.Dropdown(
        options=[("(alle)", "")], value="", description="Class:",
        layout=widgets.Layout(width="200px"),
        style={"description_width": "70px"},
    )
    refresh_btn = widgets.Button(
        description="Refresh", icon="refresh", button_style="info",
        layout=widgets.Layout(width="110px"),
    )
    live_html = widgets.HTML(value=_format_live_state({}))
    summary_html = widgets.HTML(value="")
    timeline_html = widgets.HTML(value="")
    cost_html = widgets.HTML(value="")
    mcp_html = widgets.HTML(value="")
    hooks_html = widgets.HTML(value="")
    settings_html = widgets.HTML(value="")
    worktrees_html = widgets.HTML(value="")
    schedules_html = widgets.HTML(value="")

    # ---- state ---------------------------------------------------------
    state: dict = {
        "all_outcomes": [],
        "_live_timer": None,
        "_live_stop": False,
        # Cached counts so we can update Accordion section titles cheaply
        "_counts": {"mcp": 0, "hooks": 0, "settings": 0, "worktrees": 0, "schedules": 0},
    }

    # ---- A1: Live snapshot collector + auto-refresher ------------------

    def _collect_live_snapshot() -> dict:
        """Pull the current engine + jobs status from ctx (read-only)."""
        snap: dict = {
            "engine_status": "missing",
            "current_mode": "",
            "current_model": "",
            "current_provider": "",
            "session_cost_usd": 0.0,
            "session_turns": 0,
            "active_jobs": {},
            "active_session_id": "",
            "perm_profile": "",
        }
        agent_state = getattr(ctx, "agent_state", None) or {}
        engine = agent_state.get("engine") or getattr(ctx, "agent_engine", None)

        if engine is not None:
            snap["engine_status"] = "streaming" if agent_state.get("streaming") else "idle"
            snap["current_mode"] = getattr(engine, "mode", "") or ""
            snap["current_model"] = getattr(getattr(engine, "client", None), "model", "") or ""
            snap["current_provider"] = getattr(engine, "provider", "") or ""
            try:
                snap["session_cost_usd"] = float(getattr(engine, "cost_usd", 0.0) or 0.0)
            except (TypeError, ValueError):
                pass
            try:
                snap["session_turns"] = len(getattr(engine, "messages", []) or []) // 2
            except Exception:
                pass
            snap["active_session_id"] = getattr(engine, "session_id", "") or ""

        snap["perm_profile"] = agent_state.get("_perm_profile", "") or ""

        backend = getattr(ctx, "backend", None)
        if backend is not None and hasattr(backend, "list_jobs"):
            try:
                from collections import Counter as _C
                jobs = list(backend.list_jobs() or [])
                if jobs:
                    statuses = _C(
                        (str(getattr(j, "status", "") or "?")).upper() for j in jobs
                    )
                    snap["active_jobs"] = dict(statuses.most_common())
            except Exception:
                pass

        return snap

    def _refresh_live() -> None:
        try:
            live_html.value = _format_live_state(_collect_live_snapshot())
        except Exception:
            pass

    def _live_loop() -> None:
        """Periodic refresh — 5 s tick, kooperative stop via _live_stop."""
        import threading as _th
        if state.get("_live_stop"):
            return
        _refresh_live()
        timer = _th.Timer(5.0, _live_loop)
        timer.daemon = True
        state["_live_timer"] = timer
        timer.start()

    def _reload() -> None:
        outcomes = load_outcomes(path=history_path, max_entries=500)
        state["all_outcomes"] = outcomes
        # Refresh dropdown options based on what we have
        provider_dd.options = _options_with_blank([o.provider for o in outcomes])
        mode_dd.options = _options_with_blank([o.mode for o in outcomes])
        class_dd.options = _options_with_blank([o.task_class for o in outcomes])
        # Refresh MCP server inventory
        try:
            servers = discover_mcp_servers()
            mcp_html.value = render_mcp_overview_html(servers)
            state["_counts"]["mcp"] = len(servers or [])
        except Exception:
            mcp_html.value = ""
            state["_counts"]["mcp"] = 0
        # Refresh Claude hook inventory (read-only).
        try:
            project_path = (Path.cwd() / ".claude" / "settings.json")
            hooks = discover_hooks(project_path=project_path)
            hooks_html.value = render_hooks_html(hooks)
            state["_counts"]["hooks"] = len(hooks or [])
        except Exception:
            hooks_html.value = ""
            state["_counts"]["hooks"] = 0
        # Refresh Claude settings inventory (read-only).
        try:
            project_path = (Path.cwd() / ".claude" / "settings.json")
            views = discover_settings(project_path=project_path)
            settings_html.value = render_settings_html(views)
            state["_counts"]["settings"] = len(views or [])
        except Exception:
            settings_html.value = ""
            state["_counts"]["settings"] = 0
        # Refresh git worktree inventory (read-only).
        try:
            wts = list_worktrees(Path.cwd())
            worktrees_html.value = render_worktrees_html(wts)
            state["_counts"]["worktrees"] = len(wts or [])
        except Exception:
            worktrees_html.value = ""
            state["_counts"]["worktrees"] = 0
        # Refresh scheduled-task inventory (read-only).
        try:
            scheds = discover_schedules()
            schedules_html.value = render_schedules_html(scheds)
            state["_counts"]["schedules"] = len(scheds or [])
        except Exception:
            schedules_html.value = ""
            state["_counts"]["schedules"] = 0
        _render()
        _refresh_live()
        _update_section_titles()

    def _render() -> None:
        outcomes = _filter_outcomes(
            state["all_outcomes"],
            provider=provider_dd.value,
            mode=mode_dd.value,
            verdict=verdict_dd.value,
            task_class=class_dd.value,
        )
        summary_html.value = _render_summary(_aggregate_stats(outcomes))
        timeline_html.value = _render_timeline(outcomes)
        # A3 — cost insights mirror the same filtered set
        cost_html.value = _render_cost_insights(
            _aggregate_costs_by_period(outcomes),
            _aggregate_costs_by_mode(outcomes),
        )

    # ---- A2: Accordion layout with live count badges -------------------

    filter_row = widgets.HBox(
        [provider_dd, mode_dd, verdict_dd, class_dd, refresh_btn],
        layout=widgets.Layout(gap="8px", flex_flow="row wrap"),
    )
    recent_runs = widgets.VBox([filter_row, summary_html, timeline_html])
    cost_insights = widgets.VBox([cost_html])
    configuration = widgets.VBox([mcp_html, hooks_html, settings_html])
    local_state = widgets.VBox([worktrees_html, schedules_html])

    accordion = widgets.Accordion(
        children=[recent_runs, cost_insights, configuration, local_state],
        selected_index=0,
    )
    accordion.set_title(0, "Recent runs")
    accordion.set_title(1, "Cost insights")
    accordion.set_title(2, "Configuration")
    accordion.set_title(3, "Local state")

    def _update_section_titles() -> None:
        c = state["_counts"]
        n_outcomes = len(state['all_outcomes'])
        accordion.set_title(0, f"Recent runs ({n_outcomes} cycles)")
        # Cost-insights title gets the honest "all-time" total in the badge
        all_time = _outcomes_total_delta(state["all_outcomes"])
        accordion.set_title(1, f"Cost insights ({_fmt_cost(all_time)} all-time)")
        accordion.set_title(
            2,
            f"Configuration ({c['mcp']} MCP · {c['hooks']} hooks · {c['settings']} settings)"
        )
        accordion.set_title(
            3,
            f"Local state ({c['worktrees']} worktrees · {c['schedules']} schedules)"
        )

    # ---- wiring --------------------------------------------------------
    refresh_btn.on_click(lambda _b: _reload())
    for dd in (provider_dd, mode_dd, verdict_dd, class_dd):
        dd.observe(lambda _change: _render(), names="value")

    _reload()
    _live_loop()  # kicks off the 5-s refresh chain

    return widgets.VBox(
        [title, live_html, accordion],
        layout=widgets.Layout(padding="8px"),
    )


__all__ = [
    "create_tab",
    "_aggregate_costs_by_mode",
    "_aggregate_costs_by_period",
    "_aggregate_stats",
    "_filter_outcomes",
    "_format_live_state",
    "_outcome_cost_delta",
    "_outcomes_total_delta",
    "_render_cost_insights",
    "_render_summary",
    "_render_timeline",
]
