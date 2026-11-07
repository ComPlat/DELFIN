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


def _exact_hit(text: Any, value: Any) -> bool:
    """Exact-token containment. A value is recalled only when it appears as
    its own whitespace-delimited token (or contiguous token run) in the text
    — so "28.5" is NOT recalled by an emission of "128.5", and the minus
    sign of "-37.460" is not eaten by a larger token "1-37.460". Exact
    numbers count only if EXACT (brief item 2)."""
    tv = _norm(value).split()
    if not tv:
        return False
    toks = _norm(text).split()
    width = len(tv)
    return any(toks[i:i + width] == tv for i in range(len(toks) - width + 1))


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
    actions: Optional[List[Dict[str, Any]]] = None,
) -> LongHorizonResult:
    """Reduce one session's records to the long-horizon picture. Pure.

    ``planted_facts``: [{"fact_id", "value", "probe"}]; a fact is recalled
    late when its EXACT value appears in a tool record at position >=
    ``probe`` (1-based position in ts order). ``decisions``:
    [{"declaration_id", "subject", "probe"}]; honoured late when the
    subject still appears once at position >= ``probe``.

    ``actions``: the caller's per-action goal verdicts [{"ts", "in_goal"}];
    the reducer never judges ``in_goal`` itself (the stated goal is the
    caller's knowledge, exactly as in longrun_metric). Drift is the share of
    judged actions that are out of goal.
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
      # a hit rate is a fraction in [0,1]; cached>input is unphysical data,
      # clamp so the rate can never exceed 1.0 (never a misleading >1).
      res.cache_hit_rate = min(1.0, cach_tok / in_tok)
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
            # "exact numbers count only if exact": the exact value token must
            # appear as its own token. "-37.460" is recalled by "E = -37.460"
            # but NOT by "E = -37.46" (rounded), "128.5" or "1-37.460".
            if _exact_hit(rec.get("output"), value) or _exact_hit(rec.get("input"), value):
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
            if _exact_hit(rec.get("output"), subject) or _exact_hit(rec.get("input"), subject):
                res.decisions_honoured_late += 1
                break

    # --- goal drift: the KG4 definition, reused on this reducer's inputs
    # (longrun_metric.py:171-179). The caller judges each action's in_goal;
    # the reducer only counts the share that drifted. Untested elsewhere:
    # drift is the share of JUDGED actions (in_goal not None) that are out.
    judged = [a for a in (actions or [])
              if a.get("in_goal") is not None]
    if judged:
        drift = sum(1 for a in judged if not a.get("in_goal"))
        res.drift_actions = drift
        res.drift_fraction = drift / len(judged)
    else:
        res.drift_actions = 0
        res.drift_fraction = None   # no judged action -> unmeasurable

    # --- redone work: a tool call that runs in a LATER compaction epoch than
    # (strictly after a compaction boundary) and whose identity was already
    # seen in an EARLIER epoch == work redone. Every boundary counts (a call
    # first seen between compaction1 and compaction2 and repeated after
    # compaction2 is rework across the second boundary, not the first). Each
    # post-boundary occurrence is counted once. Identity = (tool, input).
    comps_ts = sorted(x for x in (_ts(c) for c in (compactions or []))
                      if x is not None)
    if comps_ts:
        seen_all_prior: set = set()   # identities seen in strictly earlier epochs
        current_epoch: set = set()    # identities seen in the epoch we are in
        ci = 0
        ordered = sorted((_ts(r), i)
                         for i, r in enumerate(tools)
                         if _ts(r) is not None)
        for ts, i in ordered:
            while ci < len(comps_ts) and comps_ts[ci] < ts:
                seen_all_prior |= current_epoch
                current_epoch = set()
                ci += 1
            key = _identity(tools[i])
            if not key:
                continue
            if key in seen_all_prior:
                res.redone_work += 1
            current_epoch.add(key)

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
