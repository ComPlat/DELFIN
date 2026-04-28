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


def _render_timeline(outcomes: list[CycleOutcome], limit: int = 100) -> str:
    """Render the recent-cycles timeline as an HTML table."""
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
        rows.append(
            '<tr style="border-bottom:1px solid #f3f4f6;">'
            f'<td style="padding:6px 8px;font-family:monospace;font-size:11px;'
            f'color:#6b7280;">{_html.escape(_fmt_ts(o.timestamp))}</td>'
            f'<td style="padding:6px 8px;">{verdict_html}{denied}{err}{retries}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#374151;">'
            f'{_html.escape(o.mode or "—")}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#374151;">'
            f'{_html.escape(o.provider or "—")}/{_html.escape(o.model or "")}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#374151;">'
            f'{_html.escape(o.task_class or "—")}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#6366f1;'
            f'text-align:right;">{_fmt_cost(o.cost_usd)}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#14b8a6;'
            f'text-align:right;">{_fmt_duration(o.duration_s)}</td>'
            f'<td style="padding:6px 8px;font-size:11px;color:#1f2937;">'
            f'{_html.escape(_truncate(o.task, 80))}</td>'
            '</tr>'
        )
    header = (
        '<tr style="background:#f9fafb;border-bottom:2px solid #e5e7eb;">'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Time</th>'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Verdict</th>'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Mode</th>'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Provider</th>'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Class</th>'
        '<th style="padding:6px 8px;text-align:right;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Cost</th>'
        '<th style="padding:6px 8px;text-align:right;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Duration</th>'
        '<th style="padding:6px 8px;text-align:left;font-size:10px;color:#6b7280;'
        'text-transform:uppercase;letter-spacing:0.4px;">Task</th>'
        '</tr>'
    )
    return (
        '<div style="max-height:520px;overflow-y:auto;border:1px solid #e5e7eb;'
        'border-radius:6px;background:#fff;">'
        '<table style="width:100%;border-collapse:collapse;font-family:'
        '-apple-system,BlinkMacSystemFont,sans-serif;">'
        f'<thead>{header}</thead><tbody>{"".join(rows)}</tbody></table></div>'
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
    summary_html = widgets.HTML(value="")
    timeline_html = widgets.HTML(value="")
    mcp_html = widgets.HTML(value="")
    hooks_html = widgets.HTML(value="")
    settings_html = widgets.HTML(value="")
    worktrees_html = widgets.HTML(value="")
    schedules_html = widgets.HTML(value="")

    # ---- state ---------------------------------------------------------
    state: dict = {"all_outcomes": []}

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
        except Exception:
            mcp_html.value = ""
        # Refresh Claude hook inventory (read-only).
        try:
            project_path = (Path.cwd() / ".claude" / "settings.json")
            hooks = discover_hooks(project_path=project_path)
            hooks_html.value = render_hooks_html(hooks)
        except Exception:
            hooks_html.value = ""
        # Refresh Claude settings inventory (read-only).
        try:
            project_path = (Path.cwd() / ".claude" / "settings.json")
            views = discover_settings(project_path=project_path)
            settings_html.value = render_settings_html(views)
        except Exception:
            settings_html.value = ""
        # Refresh git worktree inventory (read-only).
        try:
            wts = list_worktrees(Path.cwd())
            worktrees_html.value = render_worktrees_html(wts)
        except Exception:
            worktrees_html.value = ""
        # Refresh scheduled-task inventory (read-only).
        try:
            scheds = discover_schedules()
            schedules_html.value = render_schedules_html(scheds)
        except Exception:
            schedules_html.value = ""
        _render()

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

    # ---- wiring --------------------------------------------------------
    refresh_btn.on_click(lambda _b: _reload())
    for dd in (provider_dd, mode_dd, verdict_dd, class_dd):
        dd.observe(lambda _change: _render(), names="value")

    _reload()

    filter_row = widgets.HBox(
        [provider_dd, mode_dd, verdict_dd, class_dd, refresh_btn],
        layout=widgets.Layout(gap="8px", flex_flow="row wrap"),
    )
    return widgets.VBox(
        [title, filter_row, summary_html, timeline_html,
         mcp_html, hooks_html, settings_html, worktrees_html, schedules_html],
        layout=widgets.Layout(padding="8px"),
    )


__all__ = ["create_tab", "_aggregate_stats", "_filter_outcomes",
           "_render_summary", "_render_timeline"]
