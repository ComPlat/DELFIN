"""What did recent agent work cost, per model?

A read-only diagnostic over the two real metric stores:

* ``delfin.agent.agent_metrics`` — the JSONL of finished ``TurnMetrics``
  records (``read_turns`` / ``aggregate_by_model`` / ``compare_windows``).
  This is the ONLY store that carries ``cost_usd``.
* ``delfin.agent.turn_metrics`` — per-session turn logs with
  ``input_tokens`` / ``output_tokens`` / ``cached_tokens`` but no cost.

Findings baked into the design (each surfaced by the report itself, not
hidden):

* ``aggregate_by_model`` does NOT aggregate ``input_tokens``, even though
  ``TurnMetrics.input_tokens`` is recorded. This module therefore sums
  input tokens itself from the raw records instead of trusting the
  aggregate.
* The turn_metrics session logs cannot be joined to a model on their own
  (``model`` is recorded but cost is not), so the warm-vs-cold view is
  built from the agent_metrics records' own input/output token sums;
  per-session cache behaviour stays a turn_metrics question and is
  reported only when sessions exist.

Usage: ``python -m delfin.agent.report_costs``
"""

from __future__ import annotations

import sys
import time
from typing import Any

from .agent_metrics import aggregate_by_model, compare_windows, read_turns
from . import turn_metrics

# The recent window for the headline table: 7 days.
RECENT_WINDOW_S = 7 * 86_400


def _safe_records(records: list[dict]) -> list[dict]:
    """Drop records whose numeric metric fields fail to coerce.

    aggregate_by_model() raises ValueError on a malformed numeric
    field (e.g. cost_usd="x"); the "never raise" contract here means
    such records are excluded rather than crashing the report.
    """
    out = []
    for r in records:
        try:
            for key in ("ts", "cost_usd", "delegated_cost_usd",
                        "duration_s"):
                float(r.get(key) or 0)
            for key in ("input_tokens", "output_tokens", "tool_calls",
                        "tool_errors", "delegate_count"):
                int(r.get(key) or 0)
        except (TypeError, ValueError):
            continue
        out.append(r)
    return out


def collect(
    *, since_s: float = RECENT_WINDOW_S, path: Any = None,
) -> dict:
    """Gather per-model cost/token aggregates over the recent window.

    Returns a dict (never raises on missing data — empty stores yield an
    empty table):

    ``models``: model -> {n_turns, input_tokens, output_tokens,
    total_cost_usd, total_delegated_cost_usd, avg_cost_usd}
    ``window_s``: the window actually used
    ``window_comparison``: model -> compare_windows() result (models with
    enough history only)
    ``turn_metrics_sessions``: number of session logs found by
    turn_metrics (informational)
    """
    records = _safe_records(read_turns(path=path))
    if since_s is not None:
        cut = time.time() - float(since_s)
        records = [r for r in records if float(r.get("ts") or 0) >= cut]

    aggregates = aggregate_by_model(records, path=path)

    models: dict[str, dict[str, Any]] = {}
    for model, agg in aggregates.items():
        rows = [r for r in records if r.get("model") == model]
        # Input tokens are in the records but NOT in the aggregate —
        # sum them here. Missing values count as 0, not as missing data,
        # matching aggregate_by_model's own convention for output tokens.
        in_tok = sum(int(r.get("input_tokens") or 0) for r in rows)
        out_tok = sum(int(r.get("output_tokens") or 0) for r in rows)
        models[model] = {
            "n_turns": int(agg.get("n_turns") or len(rows)),
            "input_tokens": in_tok,
            "output_tokens": out_tok,
            "total_cost_usd": float(agg.get("total_cost_usd") or 0.0),
            "total_delegated_cost_usd": float(
                agg.get("total_delegated_cost_usd") or 0.0),
            "avg_cost_usd": float(agg.get("avg_cost_usd") or 0.0),
        }

    comparison: dict[str, Any] = {}
    for model in models:
        result = compare_windows(model=model, path=path)
        if result is not None:
            comparison[model] = result

    sessions = 0
    try:
        # Informational only: how many session logs exist. turn_metrics
        # has no "list sessions" helper, so count via the metrics dir.
        from .turn_metrics import _DIR  # private but stable in-repo
        sessions = len(list(_DIR.glob("*.jsonl")))
    except Exception:
        sessions = 0

    return {
        "models": models,
        "window_s": since_s,
        "window_comparison": comparison,
        "turn_metrics_sessions": sessions,
    }


def format_text(data: dict) -> str:
    """Render collect()'s output as a compact table."""
    models = data.get("models") or {}
    lines: list[str] = []
    if not models:
        return "No agent turn records in the recent window."

    window_days = (data.get("window_s") or RECENT_WINDOW_S) / 86_400.0
    lines.append(f"Agent costs per model (last {window_days:g} days):")
    header = (
        f"{'model':<28}{'turns':>6}{'in_tok':>12}{'out_tok':>12}"
        f"{'cost$':>9}{'deleg$':>9}{'avg$/turn':>10}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    total_cost = 0.0
    for model in sorted(models, key=lambda m: -models[m]["total_cost_usd"]):
        row = models[model]
        total_cost += row["total_cost_usd"]
        lines.append(
            f"{model[:28]:<28}{row['n_turns']:>6}"
            f"{row['input_tokens']:>12}{row['output_tokens']:>12}"
            f"{row['total_cost_usd']:>9.2f}"
            f"{row['total_delegated_cost_usd']:>9.2f}"
            f"{row['avg_cost_usd']:>10.4f}"
        )
    lines.append("-" * len(header))
    lines.append(f"Total direct cost: ${total_cost:.2f}")

    comparison = data.get("window_comparison") or {}
    if comparison:
        lines.append("")
        lines.append("Recent-vs-earlier window (needs >=5 turns per window):")
        for model, cmp_data in sorted(comparison.items()):
            for field_name in (
                "avg_cost_usd", "avg_delegated_cost_usd",
                "total_cost_usd_incl_delegated",
            ):
                delta = cmp_data.get(field_name) or {}
                change = delta.get("delta")
                improved = delta.get("improved")
                if change is None:
                    continue
                mark = ""
                if improved is True:
                    mark = " (improved)"
                elif improved is False:
                    mark = " (worse)"
                lines.append(
                    f"  {model[:28]:<28}{field_name[:30]:<30}"
                    f"{change:+.4f}{mark}"
                )
    else:
        lines.append("")
        lines.append("No model has >=5 turns in both windows; "
                     "no recent-vs-earlier comparison possible.")

    n_sessions = data.get("turn_metrics_sessions") or 0
    if n_sessions:
        lines.append(f"turn_metrics session logs: {n_sessions}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    data = collect()
    text = format_text(data)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
