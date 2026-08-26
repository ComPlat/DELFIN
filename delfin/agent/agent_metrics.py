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


# ---------------------------------------------------------------------------
# What a turn delegated
#
# A turn's own tokens are billed to the parent model. The runs it hands to
# sub-agents are billed separately and were recorded in a different file
# (``~/.delfin/subagent_telemetry.jsonl``) that nothing ever joined back to
# a turn. Every number a user reads for a turn -- the live terminal line,
# ``/cost``, ``/usage``, ``get_status`` -- was therefore the cost of the
# part of the turn the parent model produced itself, and a turn that handed
# the work to five sub-agents reported a few cents.
#
# The delegated half is carried BESIDE the direct half, never folded into
# it: "this session cost $6" and "this session cost $6, of which $5.40 was
# delegated" call for different actions, and one figure cannot say which.
# ---------------------------------------------------------------------------


@dataclass
class DelegateSpend:
    """Tokens and USD spent by delegated runs, never mixed with the parent's.

    ``cost_usd`` is MEASURED spend only, exactly like the engine's own
    counter: a delegate on a model with no published rate adds nothing to
    it and increments ``unpriced`` instead. Adding an unpriced run as 0.0
    would report "we cannot price this" as "this was free", which is the
    same defect as an applied bound that is never named.
    """

    count: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    cost_usd: float = 0.0
    unpriced: int = 0

    def add_payload(self, payload: Any, *, provider: str = "") -> bool:
        """Book one finished delegate from the runner's own return value.

        Returns False when there was nothing to book, so a caller cannot
        mistake a malformed payload for a run that cost nothing.
        """
        if not isinstance(payload, dict):
            return False
        try:
            tin = max(0, int(payload.get("input_tokens") or 0))
        except (TypeError, ValueError):
            tin = 0
        try:
            tout = max(0, int(payload.get("output_tokens") or 0))
        except (TypeError, ValueError):
            tout = 0
        self.count += 1
        self.input_tokens += tin
        self.output_tokens += tout
        try:
            from . import pricing
            usd = pricing.resolve(
                str(payload.get("model") or ""), provider).cost(tin, tout)
        except Exception:
            usd = None
        if usd is None:
            self.unpriced += 1
        else:
            self.cost_usd += float(usd)
        return True

    def merge(self, other: "DelegateSpend") -> None:
        """Add another tally into this one, field for field."""
        self.count += int(other.count)
        self.input_tokens += int(other.input_tokens)
        self.output_tokens += int(other.output_tokens)
        self.cost_usd += float(other.cost_usd)
        self.unpriced += int(other.unpriced)

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Any) -> "DelegateSpend":
        """Rebuild from a persisted record. A missing or unreadable record
        is a fresh tally, never an exception: a session saved before this
        existed still has to load."""
        if not isinstance(data, dict):
            return cls()
        out = cls()
        for name in ("count", "input_tokens", "output_tokens", "unpriced"):
            try:
                setattr(out, name, max(0, int(data.get(name) or 0)))
            except (TypeError, ValueError):
                pass
        try:
            out.cost_usd = max(0.0, float(data.get("cost_usd") or 0.0))
        except (TypeError, ValueError):
            pass
        return out


@dataclass
class DelegateLedger:
    """The three questions a delegated run's cost has to answer.

    ``turn``
        Delegates the CURRENT turn started and waited for. Reset at the
        start of every turn.
    ``session``
        Every delegate this session, whenever it finished. ``turn`` and
        ``background`` are both subsets of it.
    ``background``
        The part of ``session`` that finished outside any turn -- a
        backgrounded delegate returns long after the turn that started it,
        so there is no turn it can honestly be charged to.

    A delegate is booked into ``session`` and into exactly ONE of ``turn``
    or ``background``, which is what makes double-counting impossible
    rather than merely unlikely.
    """

    turn: DelegateSpend = field(default_factory=DelegateSpend)
    session: DelegateSpend = field(default_factory=DelegateSpend)
    background: DelegateSpend = field(default_factory=DelegateSpend)

    def as_dict(self) -> dict:
        return {
            "turn": self.turn.as_dict(),
            "session": self.session.as_dict(),
            "background": self.background.as_dict(),
        }

    @classmethod
    def from_dict(cls, data: Any) -> "DelegateLedger":
        """Rebuild a session's delegated spend.

        ``turn`` is deliberately NOT restored: a session that has just
        been resumed has no turn in flight, so carrying a per-turn bucket
        would charge the first turn after the resume for delegates that
        belonged to a turn in another process.
        """
        if not isinstance(data, dict):
            return cls()
        return cls(
            turn=DelegateSpend(),
            session=DelegateSpend.from_dict(data.get("session")),
            background=DelegateSpend.from_dict(data.get("background")),
        )


def format_turn_delegation(direct_cost_usd: float,
                           spend: DelegateSpend | None) -> str:
    """The one line a turn owes when it delegated -- and nothing when it
    did not.

    A "delegated $0.00" line after every turn is noise that trains a
    reader to skip the place the real number appears, so the empty case
    returns the empty string and the caller prints nothing at all.

    Tokens are stated unconditionally because they are measured on every
    provider; the USD split only appears when there is a rate behind it,
    so a KIT/local run says what it spent instead of printing a $0.0000
    that reads as free.
    """
    if spend is None or int(getattr(spend, "count", 0) or 0) <= 0:
        return ""
    parts = [
        f"this turn delegated {spend.count} sub-agent"
        f"{'' if spend.count == 1 else 's'}: "
        f"↑{spend.input_tokens:,} ↓{spend.output_tokens:,} tokens"
    ]
    if spend.cost_usd > 0:
        direct = float(direct_cost_usd or 0.0)
        parts.append(
            f"${direct:.4f} direct + ${spend.cost_usd:.4f} delegated "
            f"= ${direct + spend.cost_usd:.4f}")
    if spend.unpriced:
        parts.append(
            f"{spend.unpriced} on a model with no published rate "
            f"(not in the figure)")
    return " · ".join(parts)


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
    # What the turn DELEGATED, kept apart from what it spent itself. The
    # columns above describe the parent model only; without these a turn
    # that handed the work to five sub-agents is recorded as one of the
    # cheapest turns in the log, and every aggregate built on cost_usd
    # ranks "delegates a lot" as "cheap".
    #
    # Same measured-only rule as cost_usd: delegated_cost_usd counts
    # delegates whose model has a published rate (or a real zero), and
    # delegates_unpriced counts the rest instead of adding them as free.
    delegate_count: int = 0                  # delegates this turn waited for
    delegated_input_tokens: int = 0
    delegated_output_tokens: int = 0
    delegated_cost_usd: float = 0.0          # measured delegate spend
    delegates_unpriced: int = 0              # delegates with no rate on record
    extras: dict = field(default_factory=dict)  # for future fields

    @property
    def total_cost_usd(self) -> float:
        """Direct + delegated.

        Derived rather than stored: a total written beside its parts is a
        third number that can contradict them after any edit, and the
        JSONL record is the parts.
        """
        return float(self.cost_usd) + float(self.delegated_cost_usd)


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
      avg_cost_usd             mean DIRECT cost per turn
      total_cost_usd           sum of the direct cost
      total_delegated_cost_usd sum of what those turns delegated
      avg_delegated_cost_usd   mean delegated cost per turn
      total_cost_usd_incl_delegated  direct + delegated
      delegating_turns         turns that delegated at least one run
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
        # The half the aggregates could not see. Kept as its own sum so
        # "this model got dearer" stays distinguishable from "this model
        # started delegating" -- the two call for opposite changes.
        deleg_total = sum(
            float(r.get("delegated_cost_usd") or 0) for r in rows)
        delegating = sum(
            1 for r in rows if int(r.get("delegate_count") or 0) > 0)
        dur_total = sum(float(r.get("duration_s") or 0) for r in rows)
        out_tok_total = sum(int(r.get("output_tokens") or 0) for r in rows)
        cont_total = sum(int(r.get("continuation_rounds") or 0) for r in rows)
        out[model] = {
            "n_turns": len(rows),
            "avg_duration_s": dur_total / n,
            "avg_cost_usd": cost_total / n,
            "total_cost_usd": cost_total,
            "total_delegated_cost_usd": deleg_total,
            "avg_delegated_cost_usd": deleg_total / n,
            "total_cost_usd_incl_delegated": cost_total + deleg_total,
            "delegating_turns": delegating,
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
        # Direct spend can fall while the total rises, because the work
        # moved to sub-agents. Comparing only the direct figure reports
        # that as an improvement.
        "avg_delegated_cost_usd": _delta(
            "avg_delegated_cost_usd", "lower_is_better"),
        "total_cost_usd_incl_delegated": _delta(
            "total_cost_usd_incl_delegated", "lower_is_better"),
        "avg_duration_s":       _delta("avg_duration_s",       "lower_is_better"),
        "tool_error_rate":      _delta("tool_error_rate",      "lower_is_better"),
        "silent_exit_rate":     _delta("silent_exit_rate",     "lower_is_better"),
        "user_correction_rate": _delta("user_correction_rate", "lower_is_better"),
        "avg_continuation_rounds": _delta(
            "avg_continuation_rounds", "lower_is_better",
        ),
    }


__all__ = [
    "DelegateLedger",
    "DelegateSpend",
    "format_turn_delegation",
    "TurnMetrics",
    "record_turn",
    "read_turns",
    "aggregate_by_model",
    "compare_windows",
]
