"""Per-turn agent metrics for iterative behaviour optimisation.

Records one JSONL line per finished agent turn into
``~/.delfin/agent_metrics.jsonl`` so we can answer:

- "Did this model get cheaper at the same task after we changed
  its profile?"  — diff median ``cost_usd_delta`` per task class
  before / after a profile edit.
- "Which model handles tab-open requests cheapest while still
  completing them?" — group by ``model`` × ``task_class``.
- "Is the agent's tool-routing accuracy improving?" — tool error
  rate over time per model.
- "How often is the agent silently producing nothing?" — silent-
  exit rate per model.

The recorder is best-effort + append-only: the file can be deleted
at any time, never blocks the agent path on failure, trimmed at 50k
entries (one entry ≈ 250 B → ~12 MB cap). The schema is forward-
compatible: new fields appended as null-defaulted columns.

Slash command: ``/agents metrics`` prints per-model aggregates.
"""

from __future__ import annotations

import json
import os
import time
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


_LOG_PATH = Path.home() / ".delfin" / "agent_metrics.jsonl"
_MAX_LINES = 50_000


@dataclass
class TurnMetrics:
    """One finished agent turn. Mirrors the per-turn footer plus a
    few aggregate-friendly columns."""

    ts: float = 0.0                          # unix timestamp
    session_id: str = ""
    model: str = ""                          # active model
    profile_name: str = ""                   # registered profile or 'default'
    mode: str = ""                           # solo / dashboard / plan / quick / ...
    effort: str = ""                         # low / medium / high / xhigh
    duration_s: float = 0.0                  # wall-clock
    input_tokens: int = 0                    # this turn only (delta)
    output_tokens: int = 0
    cost_usd: float = 0.0                    # this turn only (delta)
    tool_calls: int = 0                      # number of tool dispatches this turn
    tool_errors: int = 0                     # tools that returned {"error": ...}
    continuation_rounds: int = 0             # auto-exec rounds in the dashboard
    silent_exit: bool = False                # worker finished with zero tokens
    cooperative_stop: bool = False           # stale-kill watchdog fired
    task_class: str = ""                     # simple / moderate / complex (heuristic)
    user_corrected: bool = False             # user's NEXT message was a redirect
    extras: dict = field(default_factory=dict)  # for future fields


def record_turn(metrics: "TurnMetrics", *, path: Path | None = None) -> None:
    """Append one TurnMetrics record to the JSONL log. Best-effort —
    never raises."""
    p = path or _LOG_PATH
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        rec = asdict(metrics)
        if rec.get("ts", 0.0) == 0.0:
            rec["ts"] = time.time()
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        # Lazy trim
        try:
            if p.stat().st_size > 18_000_000:
                lines = p.read_text(encoding="utf-8").splitlines()
                if len(lines) > _MAX_LINES:
                    tail = lines[-_MAX_LINES:]
                    p.write_text("\n".join(tail) + "\n", encoding="utf-8")
                os.chmod(p, 0o600)
        except OSError:
            pass
    except OSError:
        return


def read_turns(
    *, last_n: int | None = None, path: Path | None = None,
) -> list[dict]:
    """Load metric records. Newest entries are at the end."""
    p = path or _LOG_PATH
    if not p.exists():
        return []
    out: list[dict] = []
    try:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    except OSError:
        return []
    if last_n is not None and last_n > 0:
        return out[-last_n:]
    return out


def aggregate_by_model(
    records: list[dict] | None = None,
    *,
    since_s: float | None = None,
    path: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Return per-model aggregates suitable for ``/agents metrics``.

    Fields per model:
      n_turns                  total turns observed
      avg_duration_s           mean wall-clock
      avg_cost_usd             mean cost per turn
      total_cost_usd           sum
      avg_tokens_out           mean output tokens
      tool_error_rate          tools_errors / max(1, tool_calls)
      silent_exit_rate         silent / total
      cooperative_stop_rate    coop-stops / total
      user_correction_rate     turns where user's next msg redirected
      avg_continuation_rounds  mean auto-exec rounds in dashboard
    """
    if records is None:
        records = read_turns(path=path)
    if since_s is not None:
        cut = time.time() - float(since_s)
        records = [r for r in records if float(r.get("ts") or 0) >= cut]
    by_model: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_model[r.get("model") or "?"].append(r)
    out: dict[str, dict[str, Any]] = {}
    for model, rows in by_model.items():
        n = max(1, len(rows))
        total_tool = sum(int(r.get("tool_calls") or 0) for r in rows)
        total_tool_err = sum(int(r.get("tool_errors") or 0) for r in rows)
        silent = sum(1 for r in rows if r.get("silent_exit"))
        coop = sum(1 for r in rows if r.get("cooperative_stop"))
        usercorr = sum(1 for r in rows if r.get("user_corrected"))
        cost_total = sum(float(r.get("cost_usd") or 0) for r in rows)
        dur_total = sum(float(r.get("duration_s") or 0) for r in rows)
        out_tok_total = sum(int(r.get("output_tokens") or 0) for r in rows)
        cont_total = sum(int(r.get("continuation_rounds") or 0) for r in rows)
        out[model] = {
            "n_turns": len(rows),
            "avg_duration_s": dur_total / n,
            "avg_cost_usd": cost_total / n,
            "total_cost_usd": cost_total,
            "avg_tokens_out": out_tok_total / n,
            "tool_error_rate": total_tool_err / max(1, total_tool),
            "silent_exit_rate": silent / n,
            "cooperative_stop_rate": coop / n,
            "user_correction_rate": usercorr / n,
            "avg_continuation_rounds": cont_total / n,
        }
    return out


def compare_windows(
    *,
    model: str,
    older_window_s: float = 7 * 86_400,
    newer_window_s: float = 86_400,
    path: Path | None = None,
) -> dict[str, Any] | None:
    """Compare the LAST ``newer_window_s`` seconds against the
    ``older_window_s`` seconds preceding them for a given model.
    Used to answer: "did this model's behaviour improve after the
    last profile change?"

    Returns ``None`` if either window has fewer than 5 turns
    (statistically too thin to compare).
    """
    records = read_turns(path=path)
    if not records:
        return None
    now = time.time()
    newer_cut = now - newer_window_s
    older_cut = newer_cut - older_window_s
    rows_old = [
        r for r in records
        if r.get("model") == model
        and older_cut <= float(r.get("ts") or 0) < newer_cut
    ]
    rows_new = [
        r for r in records
        if r.get("model") == model
        and float(r.get("ts") or 0) >= newer_cut
    ]
    if len(rows_old) < 5 or len(rows_new) < 5:
        return None
    agg_old = aggregate_by_model(rows_old)[model]
    agg_new = aggregate_by_model(rows_new)[model]

    def _delta(field_name: str, direction: str) -> dict[str, float]:
        """direction='lower_is_better' or 'higher_is_better'."""
        old = float(agg_old.get(field_name) or 0)
        new = float(agg_new.get(field_name) or 0)
        change = new - old
        if direction == "lower_is_better":
            improved = change < 0
        else:
            improved = change > 0
        return {
            "old": old, "new": new, "delta": change,
            "improved": improved,
        }

    return {
        "model": model,
        "n_old": len(rows_old),
        "n_new": len(rows_new),
        "avg_cost_usd":         _delta("avg_cost_usd",         "lower_is_better"),
        "avg_duration_s":       _delta("avg_duration_s",       "lower_is_better"),
        "tool_error_rate":      _delta("tool_error_rate",      "lower_is_better"),
        "silent_exit_rate":     _delta("silent_exit_rate",     "lower_is_better"),
        "user_correction_rate": _delta("user_correction_rate", "lower_is_better"),
        "avg_continuation_rounds": _delta(
            "avg_continuation_rounds", "lower_is_better",
        ),
    }


__all__ = [
    "TurnMetrics",
    "record_turn",
    "read_turns",
    "aggregate_by_model",
    "compare_windows",
]
