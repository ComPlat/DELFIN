"""Long-horizon benchmark metrics for one session, from its records only.

Package F2. A pure reducer in the same shape as longrun_metric (KG4): every
input is a record a REAL writer produced (turn_metrics.record,
tool_trace.record, compaction_log.record_compaction), the module does no
I/O, and anything the records cannot attest is reported as ``None`` — never
guessed, never silently a zero.

F2 MEASURES compaction and long sessions; it does NOT change them (F3 and
operator-B do). These are the yardstick every compaction change must beat
(E2-holdout style): adopt a change only if it measures better here.

The questions it answers (the brief, item 2):
  * late recall of planted facts (exact numbers count only if exact),
  * decisions/refusals still honoured late (a decision must still bind at
    the probe — the same survival definition as fact recall),
  * goal drift and redone work (KG4 already measures them; F2 reuses the
    definitions on its own inputs),
  * tokens, cost, cache hit rate, number and size of compactions, wall time.

Writers of the inputs (one line each, as in longrun_metric's docstring):
  * ``turns``              -- turn_metrics.record / turn_metrics.read.
  * ``tools``              -- tool_trace.record / tool_trace.read.
  * ``compactions``        -- compaction_log.record_compaction /
                              compaction_log.read_compactions.
  * ``planted_facts``      -- the task spec (the task yaml names its planted
                              facts and their late probes). A fact is
                              recalled late when its EXACT value still
                              appears in the session's records at/after the
                              probe position: ``probe`` is a 1-based record
                              position, and a rounded value is not a recall.
Meaning of a record ``position``: the 1-based index of the tool record in
timestamp order. The probe sits at the position the task yaml states (a
fact planted at turn 3 and needed at turn 90 has a probe near the end), so
'a value that survives to the late records' scores whereas 'only early
records carry it' does not.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

__all__ = ["LongHorizonResult", "measure", "to_json", "to_text"]


def _norm(text: Any) -> str:
    """Whitespace-normalised exact text, the same narrow normalisation
    longrun_metric uses for commands: two spellings that differ only in
    whitespace are the same call."""
    return " ".join(str(text or "").split())


def _ts(record: Dict[str, Any]) -> Optional[float]:
    ts = record.get("ts")
    if ts is None:
        return None
    try:
        return float(ts)
    except (TypeError, ValueError):
        return None


def _identity(record: Dict[str, Any]) -> str:
    """The redo-work identity: (tool name, normalised input). Two records
    are 'the same call' only when both name the same tool and the same
    arguments."""
    tool = str(record.get("tool") or "").strip()
    arg = _norm(record.get("input"))
    if not tool and not arg:
        return ""
    return f"{tool}::{arg}"


@dataclass
class LongHorizonResult:
    """One session's long-horizon picture. ``None`` under a name marks a
    value the session's records cannot attest (unmeasurable), never zero."""

    turns: int = 0
    tool_calls: int = 0
    compactions: int = 0

    # --- planted facts / decisions: do they survive to the late records?
    facts_recalled_late: int = 0
    facts_total: int = 0
    decisions_honoured_late: int = 0
    decisions_total: int = 0

    # --- work
    redone_work: int = 0          # same tool call + args repeated post-compaction
    drift_actions: int = 0
    drift_fraction: Optional[float] = None

    # --- cost / throughput from turn_metrics records
    input_tokens: int = 0
    output_tokens: int = 0
    cached_tokens: int = 0
    cache_hit_rate: Optional[float] = None
    total_ms: int = 0
    wall_time_s: Optional[float] = None

    # --- compactions
    compactions_size: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def recall_late(self) -> float:
        """Fraction of planted facts recalled at/after their probe (0.0
        when there are none). Retained under this name for the red test."""
        if not self.facts_total:
            return 0.0
        return self.facts_recalled_late / self.facts_total


def measure(
    turns: List[Dict[str, Any]],
    tools: List[Dict[str, Any]],
    compactions: List[Dict[str, Any]],
    *,
    planted_facts: Optional[List[Dict[str, Any]]] = None,
    decisions: Optional[List[Dict[str, Any]]] = None,
) -> LongHorizonResult:
    """Reduce one session's records to the long-horizon picture. Pure.

    ``planted_facts``: [{"fact_id", "value", "probe"}]; a fact is recalled
    late when its EXACT value appears in a tool record at position >=
    ``probe`` (1-based position in ts order). ``decisions``:
    [{"declaration_id", "subject", "probe"}]; honoured late when the
    subject still appears once at position >= ``probe``.
    """
    tools = list(tools or [])
    res = LongHorizonResult(turns=len(turns or []),
                            tool_calls=len(tools),
                            compactions=len(compactions or []))

    # --- tokens / cost / cache. Sum over turns with a recorded value.
    in_tok = sum(t.get("input_tokens") or 0 for t in (turns or []))
    out_tok = sum(t.get("output_tokens") or 0 for t in (turns or []))
    cach_tok = sum(t.get("cached_tokens") or 0 for t in (turns or []))
    res.input_tokens, res.output_tokens, res.cached_tokens = in_tok, out_tok, cach_tok
    res.total_ms = sum(t.get("total_ms") or 0 for t in (turns or []))
    if in_tok > 0:
        res.cache_hit_rate = cach_tok / in_tok
    else:
        res.cache_hit_rate = None   # no turns attest a cache -> unmeasurable

    # --- wall time: first turn start -> last turn end (both observed).
    tss = [_ts(t) for t in (turns or [])]
    tss = [x for x in tss if x is not None]
    ts_tools = [_ts(t) for t in tools]
    ts_tools = [x for x in ts_tools if x is not None]
    all_ts = sorted(tss + ts_tools)
    if len(all_ts) >= 2:
        res.wall_time_s = all_ts[-1] - all_ts[0]

    # --- planted facts: exact value recalled at/after the probe position.
    facts = planted_facts or []
    res.facts_total = len(facts)
    for fact in facts:
        value = _norm(fact.get("value"))
        probe = fact.get("probe") or fact.get("probe_turn") or 1
        if not value:
            continue
        for idx, rec in enumerate(tools, start=1):
            if idx < probe:
                continue
            # "exact numbers count only if exact": the exact value token must be
            # present. "-37.460" is contained in "E = -37.460" (recalled) but
            # NOT in "E = -37.46" (a rounded emission does not recall it).
            if value in _norm(rec.get("output")) or value in _norm(rec.get("input")):
                res.facts_recalled_late += 1
                break   # one recall per fact

    # --- decisions / refusals: the subject still appears late (a decision
    # that must bind late has the same survival definition as a fact).
    decs = decisions or []
    res.decisions_total = len(decs)
    for dec in decs:
        subject = _norm(dec.get("subject"))
        probe = dec.get("probe") or 1
        if not subject:
            continue
        for idx, rec in enumerate(tools, start=1):
            if idx < probe:
                continue
            if subject in _norm(rec.get("output")) or subject in _norm(rec.get("input")):
                res.decisions_honoured_late += 1
                break

    # --- redone work: a tool call that ran BEFORE a compaction boundary and
    # the same (tool, normalised input) that runs AGAIN after it. The
    # identity is (tool name, normalised input) — the same arguments on a
    # different tool are not the same call. The boundary is the compaction
    # record's ts (the sole marker; reviewer constraint 2). Instance-counting
    # like longrun_metric: each post-boundary occurrence of an identity that
    # also occurred pre-boundary counts.
    comps_ts = sorted(x for x in (_ts(c) for c in (compactions or []))
                      if x is not None)
    if comps_ts:
        boundary = comps_ts[0]
        seen_before: set = set()
        for rec in tools:
            ts = _ts(rec)
            if ts is not None and ts < boundary:
                seen_before.add(_identity(rec))
        for rec in tools:
            ts = _ts(rec)
            if ts is not None and ts > boundary:
                key = _identity(rec)
                if key and key in seen_before:
                    res.redone_work += 1

    # --- compactions size (messages_compacted) is kept raw.
    for c in (compactions or []):
        res.compactions_size.append({
            "kind": c.get("kind"),
            "messages_compacted": c.get("messages_compacted"),
            "tokens_before": c.get("tokens_before"),
            "tokens_after": c.get("tokens_after"),
        })

    return res


def to_json(result: LongHorizonResult) -> str:
    """Machine-readable report (JSON). Pure serialiser — no I/O."""
    import json
    return json.dumps(asdict(result), ensure_ascii=False, sort_keys=True)


def to_text(result: LongHorizonResult) -> str:
    """Short one-line-per-metric text report. Human-readable twin of
    to_json, still pure."""
    return "\n".join(f"{k}: {v}" for k, v in asdict(result).items()
                     if str(v) not in ("[]", ""))
