"""What did recent agent work cost, per model?

A read-only diagnostic over the two real metric stores:

* ``delfin.agent.agent_metrics`` — the JSONL of finished ``TurnMetrics``
  records (``read_turns`` / ``aggregate_by_model`` / ``compare_windows``).
  This is the ONLY store that carries ``cost_usd``.
* ``delfin.agent.turn_metrics`` — per-session turn logs with
  ``input_tokens`` / ``output_tokens`` / ``cached_tokens`` but no cost.

Findings baked into the design (each surfaced by the report itself, not
hidden):

* ``aggregate_by_model`` did NOT aggregate ``input_tokens`` when this was
  written, and raised ``ValueError`` on a malformed numeric field. Both
  are fixed at the source (it reports ``total_tokens_in`` and drops what
  it cannot read), so the sums here no longer work around anything --
  this module keeps its own per-row pass because its table is per model
  AND per window, which the aggregate does not provide together.
* A cost of 0.00 is not a free run: ``agent_metrics`` records what each
  turn's cost MEANS (measured, non-billing, unpriced), and a total of
  zero over unpriced turns says only that nobody could price them.
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
from .attribution import COLUMNS

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


# --------------------------------------------------------------------------
# Usage by attribution (molecule / workflow / task / agent)
# --------------------------------------------------------------------------
#
# The per-model view above answers "what did a model cost". It cannot answer
# "what did THIS MOLECULE cost, and which agent spent the tokens in THIS
# WORKFLOW", because a model is not an attribution. The writers put the
# columns on every row (``attribution.COLUMNS``, imported below rather than
# respelled here), and a row that names none of its ids is reported as
# unattributed rather than shown as an empty table (2026-11-26).
#
# Two rules the fields below are built around:
#
# * A 0.00 over unpriced turns is "not measured", not "free". The cost is
#   therefore always accompanied by what it MEANS (``price_state`` on the
#   record), and a group where nobody could price the turns reports
#   "unpriced" rather than "$0.00".
# * Absent is not zero. ``cached_tokens`` is not on these records at all and
#   ``reasoning_tokens`` is recorded by nobody; both are reported as None
#   ("not recorded") instead of 0, which would read as "measured zero
#   cache hits". A sum that silently counted absence as zero is the same
#   defect as an unpriced turn counted as free.

#: The attribution a row can carry. Read only if present; absent -> null.
#: Taken from the writer's own tuple rather than respelled: the two lists
#: are the two halves of one contract, and when they drifted apart the
#: writer emitted ``agent_role``/``subagent_name``/``parent_session`` while
#: this read ``role``/``subagent``/``parent_session_id``, so ``--by agent``
#: could never group a real row (measured: 2 of 2 rows unattributed).
USAGE_KEYS = COLUMNS

#: Optional counters whose absence must stay visible as None.
OPTIONAL_TOKEN_KEYS = ("cached_tokens", "reasoning_tokens")

#: Groups smaller than this are shown but marked: a mean over two turns is
#: a coincidence with a decimal point, and hiding it loses the small group
#: that is often the whole question.
MIN_GROUP_ROWS = 3

#: The bucket every grouping uses for a row that names none of its ids.
UNATTRIBUTED = "(unattributed)"

#: Grouping name -> the record fields that make up its key, in order. The
#: first one present and non-empty wins; none present -> UNATTRIBUTED.
GROUP_FIELDS: dict[str, tuple[str, ...]] = {
    "molecule": ("molecule_id",),
    "workflow": ("workflow_id",),
    "task": ("task_id",),
    # "which agent spent the tokens": the delegate that did the work, else
    # the role the turn ran as. model is NOT used -- a model is not an
    # agent, and falling back to it would report every row as attributed.
    "agent": ("subagent_name", "agent_role"),
}

USAGE_BY_CHOICES = tuple(GROUP_FIELDS)


def _present(value: Any) -> bool:
    """A value that is actually recorded (not None, not the empty string)."""
    if value is None:
        return False
    return bool(str(value).strip())


def _opt_int(value: Any) -> int | None:
    """An int, or None when the record does not carry one.

    Distinct from ``_num``, which folds "absent" into 0. Here the
    difference matters: absent means nobody recorded it.

    Deliberately NOT named ``_as_float``-style: ``data_checks._as_float``
    raises ``DataCheckError`` on a bad value, while this module's contract
    is to never raise on missing data. Same input, opposite contract -- so
    a distinct name keeps the two from being "deduped" into one.
    """
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _opt_float(value: Any) -> float | None:
    """A float, or None -- never raises (see _opt_int on the naming)."""
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _group_key(record: dict, fields: tuple[str, ...]) -> str:
    """The row's group id, or UNATTRIBUTED when it names none."""
    for field_name in fields:
        value = record.get(field_name)
        if _present(value):
            return str(value).strip()
    return UNATTRIBUTED


def _cost_meaning(rows: list[dict]) -> str:
    """What this group's cost figure MEANS, from the records themselves.

    ``price_state`` is written per turn: "measured" (the backend said so or
    a published rate applied), "non_billing" (a real zero) or "unknown"
    (no rate, nothing observed). A record written before the field existed
    carries "" and is reported as "unstated" rather than assumed measured.

    Delegated spend counts too. ``delegated_cost_usd`` sums only the
    delegates that could be priced, and ``delegates_unpriced`` is the store
    saying how many it could not (agent_metrics.py:252 -- its own comment
    says counting those as free is the defect the field exists to prevent).
    A group whose turn is priced but whose five delegates are not is
    therefore NOT "measured": the figure it prints covers part of what was
    spent. Ignoring the count here reproduced exactly the defect this
    module's header forbids, one level down.
    """
    if any((_opt_int(r.get("delegates_unpriced")) or 0) > 0 for r in rows):
        priced = any(
            str(r.get("price_state") or "").strip() in
            ("measured", "non_billing")
            or (_opt_float(r.get("cost_usd")) or 0.0) > 0
            or (_opt_float(r.get("delegated_cost_usd")) or 0.0) > 0
            for r in rows)
        return "mixed" if priced else "unpriced"
    states = [str(r.get("price_state") or "").strip() for r in rows]
    n_unknown = sum(1 for s in states if s == "unknown")
    n_non_billing = sum(1 for s in states if s == "non_billing")
    n_measured = sum(1 for s in states if s == "measured")
    n_unstated = sum(1 for s in states if not s)
    if n_unknown and (n_measured or n_non_billing or n_unstated):
        return "mixed"
    if n_unknown:
        return "unpriced"
    if n_non_billing and (n_measured or n_unstated):
        return "mixed"
    if n_non_billing:
        return "non-billing"
    if n_measured and n_unstated:
        return "mixed"
    if n_measured:
        return "measured"
    return "unstated"


def _sum_optional(rows: list[dict], key: str) -> int | None:
    """Sum a counter that may be recorded by nobody.

    None when NO row carries it -- that is "not recorded", and printing 0
    there would claim a measured zero. Otherwise the sum, with absent rows
    contributing nothing (matching the store's own convention).
    """
    values = [_opt_int(r.get(key)) for r in rows]
    if all(v is None for v in values):
        return None
    return sum(v for v in values if v is not None)


def collect_usage(
    *,
    by: str = "agent",
    since_s: float = RECENT_WINDOW_S,
    path: Any = None,
    min_rows: int = MIN_GROUP_ROWS,
) -> dict:
    """Aggregate the same stores by attribution instead of by model.

    Returns (never raises on missing data)::

        by                 the grouping used
        window_s           the window actually applied
        groups             group id -> per-group sums (see _group_totals)
        n_rows             rows inside the window
        n_unattributed     rows whose group is UNATTRIBUTED
        unattributed_share n_unattributed / n_rows (0.0 when there are none)
        columns_present    which attribution fields the store really carries
        min_rows           the mark threshold the groups were tested against

    An empty store yields n_rows 0 and no groups; that is reported as "no
    records", which is a different statement from "no attribution".
    """
    if by not in GROUP_FIELDS:
        raise ValueError(
            f"unknown grouping {by!r}; expected one of {USAGE_BY_CHOICES}")
    fields = GROUP_FIELDS[by]

    records = _safe_records(read_turns(path=path))
    n_stored = len(records)
    if since_s is not None:
        cut = time.time() - float(since_s)
        records = [r for r in records if float(r.get("ts") or 0) >= cut]

    buckets: dict[str, list[dict]] = {}
    for record in records:
        buckets.setdefault(_group_key(record, fields), []).append(record)

    groups: dict[str, dict[str, Any]] = {}
    for key, rows in buckets.items():
        groups[key] = _group_totals(rows, min_rows=min_rows)

    n_rows = len(records)
    n_unattributed = len(buckets.get(UNATTRIBUTED, []))
    return {
        "by": by,
        "window_s": since_s,
        # How many rows the store holds at all, next to how many survived
        # the window. Without this, an empty file and a file whose rows are
        # all older than the cut both read as "0 rows", and the reader
        # cannot tell "nothing was ever recorded" from "nothing in this
        # window". Those are different answers to different questions.
        "n_stored": n_stored,
        "groups": groups,
        "n_rows": n_rows,
        "n_unattributed": n_unattributed,
        "unattributed_share": (
            n_unattributed / n_rows if n_rows else 0.0
        ),
        "columns_present": [
            k for k in USAGE_KEYS if any(_present(r.get(k)) for r in records)
        ],
        "min_rows": min_rows,
    }


def _group_totals(rows: list[dict], *, min_rows: int) -> dict[str, Any]:
    """The per-group sums, each with what it means."""
    cost_usd = sum(_opt_float(r.get("cost_usd")) or 0.0 for r in rows)
    delegated = sum(
        _opt_float(r.get("delegated_cost_usd")) or 0.0 for r in rows)
    wall = sum(_opt_float(r.get("duration_s")) or 0.0 for r in rows)
    tool_calls = sum(_opt_int(r.get("tool_calls")) or 0 for r in rows)
    # One model call per turn, plus each auto-exec continuation round the
    # dashboard recorded. Stated as a basis because it is an inference
    # from what a turn is, not a counter the record carries.
    continuations = sum(
        _opt_int(r.get("continuation_rounds")) or 0 for r in rows)
    n = len(rows)
    return {
        "n_rows": n,
        "below_threshold": n < min_rows,
        "input_tokens": sum(_opt_int(r.get("input_tokens")) or 0 for r in rows),
        "output_tokens": sum(
            _opt_int(r.get("output_tokens")) or 0 for r in rows),
        "cached_tokens": _sum_optional(rows, "cached_tokens"),
        "reasoning_tokens": _sum_optional(rows, "reasoning_tokens"),
        "model_calls": n + continuations,
        "model_calls_basis": "turns + continuation_rounds",
        "tool_calls": tool_calls,
        "wall_time_s": wall,
        "cost_usd": cost_usd,
        "delegated_cost_usd": delegated,
        "cost_total_usd": cost_usd + delegated,
        # How many delegates went unpriced is said next to the figure it
        # qualifies: a 0.00 delegated sum over N unpriced delegates is not
        # a measured zero. None (not 0) when no row carries the field, so
        # an old store is "not recorded" rather than a measured zero.
        "delegates_unpriced": _sum_optional(rows, "delegates_unpriced"),
        "cost_meaning": _cost_meaning(rows),
        "unpriced_rows": sum(
            1 for r in rows
            if str(r.get("price_state") or "").strip() == "unknown"),
        "non_billing_rows": sum(
            1 for r in rows
            if str(r.get("price_state") or "").strip() == "non_billing"),
    }


def _fmt_cost(row: dict) -> str:
    """A cost figure that never reads as free when nobody priced it."""
    meaning = row.get("cost_meaning") or ""
    amount = f"${row.get('cost_total_usd', 0.0):.2f}"
    if meaning == "unpriced":
        return f"{amount} unpriced"
    if meaning == "non-billing":
        return f"{amount} non-billing"
    if meaning == "unstated":
        return f"{amount} unstated"
    if meaning == "mixed":
        return f"{amount} mixed"
    return amount


def _fmt_opt(value: int | None) -> str:
    return "-" if value is None else str(value)


def format_usage(data: dict) -> str:
    """Render collect_usage()'s output as a table, with its caveats."""
    groups = data.get("groups") or {}
    n_rows = int(data.get("n_rows") or 0)
    by = data.get("by") or "agent"
    if not n_rows:
        # An empty store and a store whose rows all fall outside the window
        # are different statements. Printing one sentence for both told a
        # reader with a full store and a narrow --since that nothing had
        # ever been recorded.
        if int(data.get("n_stored") or 0):
            return (f"No agent turn records in the window "
                    f"({int(data['n_stored'])} older records exist; widen "
                    f"--since to include them).")
        return "No agent turn records in the recent window."

    window_days = (data.get("window_s") or RECENT_WINDOW_S) / 86_400.0
    min_rows = int(data.get("min_rows") or MIN_GROUP_ROWS)
    lines = [f"Agent usage by {by} (last {window_days:g} days, "
             f"{n_rows} turns):"]

    header = (
        f"{by:<24}{'rows':>6}{'in_tok':>13}{'out_tok':>13}"
        f"{'cached':>9}{'reason':>8}{'m_calls':>9}{'t_calls':>9}"
        f"{'wall_s':>9}{'cost':>18}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    ordered = sorted(
        groups.items(),
        key=lambda kv: (kv[0] == UNATTRIBUTED, -kv[1]["cost_total_usd"],
                        -kv[1]["n_rows"]),
    )
    for key, row in ordered:
        name = key[:24]
        mark = "*" if row.get("below_threshold") else ""
        lines.append(
            f"{name:<24}{row['n_rows']:>6}"
            f"{row['input_tokens']:>13}{row['output_tokens']:>13}"
            f"{_fmt_opt(row['cached_tokens']):>9}"
            f"{_fmt_opt(row['reasoning_tokens']):>8}"
            f"{row['model_calls']:>9}{row['tool_calls']:>9}"
            f"{row['wall_time_s']:>9.1f}"
            f"{_fmt_cost(row):>18}"
            f"{mark}"
        )
    lines.append("-" * len(header))

    if any(r.get("below_threshold") for r in groups.values()):
        lines.append(f"* fewer than {min_rows} rows: shown, not comparable "
                     f"with the groups above.")

    # The absence of a field is stated, not rendered as a zero.
    missing = [k for k in OPTIONAL_TOKEN_KEYS
               if all(r.get(k) is None for r in groups.values())]
    if missing:
        lines.append(
            "not recorded by any row, so not shown as 0: "
            + ", ".join(missing))
    lines.append("model calls = turns + continuation_rounds (a turn makes "
                 "at least one call; the record carries no per-call count).")

    columns = data.get("columns_present") or []
    share = float(data.get("unattributed_share") or 0.0)
    if not columns:
        lines.append("")
        lines.append(
            f"No attribution columns on any record: all {n_rows} turns "
            f"({share:.0%}) are unattributed. Nothing can be grouped by "
            f"{by} until the writers record "
            f"{', '.join(GROUP_FIELDS.get(by, ()))}.")
    else:
        lines.append("")
        lines.append(
            f"Attribution columns present: {', '.join(columns)}. "
            f"Unattributed: {data.get('n_unattributed', 0)} of {n_rows} "
            f"turns ({share:.0%}).")
    return "\n".join(lines)


def usage_main(
    *, by: str = "agent", since_s: float = RECENT_WINDOW_S,
    path: Any = None, as_json: bool = False,
    min_rows: int = MIN_GROUP_ROWS,
) -> int:
    """The `delfin-agent usage` body. Read-only; 0 on success."""
    if by not in GROUP_FIELDS:
        print(f"unknown grouping {by!r}; expected one of "
              f"{', '.join(USAGE_BY_CHOICES)}", file=sys.stderr)
        return 2
    data = collect_usage(by=by, since_s=since_s, path=path,
                         min_rows=min_rows)
    if as_json:
        import json as _json
        print(_json.dumps(data, indent=2, sort_keys=True))
    else:
        print(format_usage(data))
    return 0


def main(argv: list[str] | None = None) -> int:
    data = collect()
    text = format_text(data)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
