"""Eval loop — turn real usage outcomes into measured improvement.

Stufe 4 of the program: the agent's outcome history (every PASS/FAIL the
dashboard records) is mined for **recurring failure patterns**, recurring
patterns are scaffolded into **draft benchmark tasks** (same draft dir
the 🐞→task flow uses), and a **report** is written so the team sees the
agent's health at a glance.  Combined with the committed suite +
``optimize_check`` this closes the frontier-style loop:
measure → find weak spots → pin them as tasks → drive them to zero.

Cost & consent rules (same philosophy as the job monitor):
- Everything this module does by default is **LLM-free (0 tokens)**:
  outcome mining, task scaffolding and report writing are pure file work.
- Running the live benchmark suite costs tokens, so it is NOT part of
  the default run — the report just prints the command for it.
- The entrypoint is gated by ``agent.eval_loop.enabled`` (default on —
  the default pass is LLM-free file analytics; set false to disable).

Run once (cron / dashboard scheduler / by hand)::

    python -m delfin.agent.eval_loop
"""

from __future__ import annotations

import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


_REPORTS_DIR = Path.home() / ".delfin" / "eval_reports"
# Same draft location the 🐞→task scaffolder uses — one review inbox.
_TASK_DRAFTS_DIR = Path.home() / ".delfin" / "bug_tasks"

# A pattern is "recurring" when at least this many failures share the
# same fingerprint within the analysis window.
_RECURRENCE_THRESHOLD = 3


@dataclass
class FailurePattern:
    fingerprint: str              # "<task_class>|<error_type>|<mode>"
    count: int
    task_class: str
    error_type: str
    mode: str
    examples: list[str] = field(default_factory=list)   # sample task texts


def eval_settings(settings: dict | None = None) -> dict:
    if settings is None:
        try:
            from delfin.user_settings import load_settings
            settings = load_settings()
        except Exception:
            settings = {}
    cfg = ((settings or {}).get("agent") or {}).get("eval_loop") or {}
    return {
        "enabled": bool(cfg.get("enabled", True)),
        "window": int(cfg.get("window", 200) or 200),
        "threshold": int(cfg.get("threshold", _RECURRENCE_THRESHOLD)
                         or _RECURRENCE_THRESHOLD),
    }


# ---------------------------------------------------------------------------
# Outcome mining (LLM-free)
# ---------------------------------------------------------------------------

def analyze_outcomes(
    outcomes: list[Any] | None = None,
    *,
    window: int = 200,
    threshold: int = _RECURRENCE_THRESHOLD,
) -> list[FailurePattern]:
    """Mine the outcome history for recurring failure fingerprints.

    Fingerprint = (task_class, error_type, mode).  Only FAIL verdicts
    count.  Returns patterns sorted by count, most frequent first.
    """
    if outcomes is None:
        try:
            from delfin.agent.outcome_tracker import load_outcomes
            outcomes = load_outcomes(max_entries=window)
        except Exception:
            outcomes = []

    def _get(o: Any, key: str) -> str:
        if isinstance(o, dict):
            return str(o.get(key) or "")
        return str(getattr(o, key, "") or "")

    fails = [o for o in outcomes if _get(o, "verdict").upper() == "FAIL"]
    buckets: dict[str, list[Any]] = {}
    for o in fails:
        fp = "|".join((
            _get(o, "task_class") or "unclassified",
            _get(o, "error_type") or "unknown",
            _get(o, "mode") or "any",
        ))
        buckets.setdefault(fp, []).append(o)

    patterns: list[FailurePattern] = []
    for fp, items in buckets.items():
        if len(items) < threshold:
            continue
        task_class, error_type, mode = fp.split("|", 2)
        patterns.append(FailurePattern(
            fingerprint=fp,
            count=len(items),
            task_class=task_class,
            error_type=error_type,
            mode=mode,
            examples=[_get(o, "task")[:120] for o in items[:3]],
        ))
    patterns.sort(key=lambda p: p.count, reverse=True)
    return patterns


# ---------------------------------------------------------------------------
# Draft-task scaffolding (LLM-free; human reviews before committing)
# ---------------------------------------------------------------------------

def scaffold_task_from_pattern(pattern: FailurePattern) -> dict:
    """Build a draft benchmark task pinning a recurring failure pattern.

    The maintainer fills in expected/forbidden signals (TODO markers),
    exactly like the 🐞→task flow — drafts are never auto-committed.
    """
    ts = time.strftime("%Y%m%d%H%M%S")
    safe = pattern.fingerprint.replace("|", "_").replace("/", "-")[:48]
    prompt = (pattern.examples[0] if pattern.examples
              else "(no example recorded — fill in a representative prompt)")
    return {
        "id": f"recur_{ts}_{safe}",
        "task_class": "regression",
        "mode": pattern.mode if pattern.mode != "any" else "solo",
        "prompt": prompt,
        "expected_signals": [
            {"pattern": "TODO-what-the-correct-answer-must-contain",
             "against": "text"},
        ],
        "max_duration_s": 90,
        "max_cost_usd": 0.30,
        "max_tool_calls": 8,
    }


def write_pattern_drafts(patterns: list[FailurePattern],
                         dest_dir: Path | None = None) -> list[Path]:
    """Write one draft YAML per pattern (reuses the bug-task renderer)."""
    from delfin.agent.bug_report import task_to_yaml
    d = dest_dir or _TASK_DRAFTS_DIR
    d.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []
    for p in patterns:
        task = scaffold_task_from_pattern(p)
        path = d / f"{task['id']}.yaml"
        path.write_text(
            task_to_yaml(task, source_report=f"outcome-pattern {p.fingerprint} "
                                             f"(x{p.count})"),
            encoding="utf-8",
        )
        out.append(path)
    return out


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

# Telemetry window the report summarises (days).
_HEALTH_WINDOW_DAYS = 7


def _tool_health_lines(window_days: int = _HEALTH_WINDOW_DAYS) -> list[str]:
    """"Tool health" report section: top tools by error rate over the
    window — the consumer the tool-trace telemetry has been missing."""
    lines = [f"## Tool health (last {window_days} days)"]
    try:
        from delfin.agent.tool_trace import aggregate_tool_stats
        stats = aggregate_tool_stats(window_days)
    except Exception:
        stats = {}
    rows = [(tool, s) for tool, s in stats.items() if s.get("calls")]
    if not rows:
        lines.append("- no tool calls recorded in window")
        return lines
    total_calls = sum(s["calls"] for _, s in rows)
    total_errors = sum(s.get("errors", 0) for _, s in rows)
    lines.append(f"- {total_calls} call(s), {total_errors} error(s) "
                 f"across {len(rows)} tool(s)")
    rows.sort(key=lambda kv: (kv[1].get("error_rate", 0.0),
                              kv[1].get("calls", 0)), reverse=True)
    for tool, s in rows[:5]:
        snips = s.get("top_error_snippets") or []
        extra = f" — e.g. \"{snips[0]}\"" if snips else ""
        lines.append(
            f"- **{tool}**: {s.get('error_rate', 0.0) * 100:.0f}% errors "
            f"({s.get('errors', 0)}/{s['calls']} calls), "
            f"avg {s.get('avg_duration_ms', 0)}ms{extra}")
    return lines


def _turn_health_lines(window_days: int = _HEALTH_WINDOW_DAYS) -> list[str]:
    """"Turn health" report section: ttft/stall summary from turn_metrics."""
    lines = [f"## Turn health (last {window_days} days)"]
    try:
        from delfin.agent.turn_metrics import aggregate_turn_stats
        s = aggregate_turn_stats(window_days)
    except Exception:
        s = {}
    if not s.get("turns"):
        lines.append("- no turns recorded in window")
        return lines
    # Crashes and never-started turns were counted and not printed, which
    # is the same defect as not counting them: a run of turns that died,
    # or that waited out a silent endpoint, left this section reading
    # exactly like a healthy week. They are the three kinds of bad turn
    # and they are told apart on purpose -- a crash HAS an error, a
    # never-started turn produced no token at all, a stall produced one
    # late -- so naming them separately is what makes the section
    # actionable rather than atmospheric.
    lines.append(f"- turns: {s['turns']}, stalls: {s.get('stalls', 0)}, "
                 f"crashes: {s.get('crashes', 0)}, "
                 f"never started: {s.get('never_started', 0)}, "
                 f"stopped: {s.get('stopped_count', 0)}")
    # The sample size sits next to the average it belongs to. A mean over
    # three turns out of ninety reads exactly like a mean over ninety.
    sample = s.get("ttft_sample")
    of = (f" over {sample} of {s['turns']} turns" if sample is not None
          else "")
    lines.append(f"- ttft: avg {s.get('avg_ttft_ms', 0) / 1000:.1f}s, "
                 f"p90 {s.get('p90_ttft_ms', 0) / 1000:.1f}s{of}")
    return lines


def _benchmark_drift_lines() -> list[str]:
    """"Benchmark drift" note: latest vs previous persisted run. Empty when
    fewer than two runs exist or anything fails (best-effort)."""
    try:
        from delfin.agent import benchmark as bm
        runs = bm.list_runs()
        if len(runs) < 2:
            return []
        ab = bm.ab_compare(str(runs[-2]), str(runs[-1]))
        note = bm.format_ab_note(ab)
        lines = ["## Benchmark drift (latest vs previous run)"]
        lines.append(f"- runs: `{runs[-2].name}` → `{runs[-1].name}`")
        lines += [f"- {ln}" for ln in note.splitlines() if ln.strip()]
        return lines
    except Exception:
        return []


def build_report(
    *,
    outcomes: list[Any] | None = None,
    window: int = 200,
    threshold: int = _RECURRENCE_THRESHOLD,
    drafts: list[Path] | None = None,
) -> str:
    """Render the eval report (markdown). LLM-free."""
    if outcomes is None:
        try:
            from delfin.agent.outcome_tracker import load_outcomes
            outcomes = load_outcomes(max_entries=window)
        except Exception:
            outcomes = []

    def _get(o: Any, key: str) -> str:
        if isinstance(o, dict):
            return str(o.get(key) or "")
        return str(getattr(o, key, "") or "")

    verdicts = Counter(_get(o, "verdict").upper() or "NONE" for o in outcomes)
    total = len(outcomes)
    fails = verdicts.get("FAIL", 0)
    patterns = analyze_outcomes(outcomes, window=window, threshold=threshold)

    # Suite integrity (free, deterministic).
    try:
        from delfin.agent.optimize_check import run_checks
        issues = run_checks()
        n_err = sum(1 for i in issues if i.severity == "error")
        n_warn = sum(1 for i in issues if i.severity == "warn")
        integrity = (f"OK — 0 errors, {n_warn} warning(s)" if n_err == 0
                     else f"BROKEN — {n_err} error(s), {n_warn} warning(s)")
    except Exception as exc:
        integrity = f"check failed: {exc}"

    lines = [
        f"# DELFIN agent eval report — {time.strftime('%Y-%m-%d %H:%M')}",
        "",
        "## Outcome window",
        f"- outcomes analysed: {total} (last {window})",
        f"- verdicts: " + ", ".join(f"{k}={v}" for k, v in
                                     sorted(verdicts.items())) if total else
        "- no outcomes recorded yet",
        f"- failure rate: {fails}/{total}"
        + (f" ({100.0 * fails / total:.0f}%)" if total else ""),
        "",
        "## Recurring failure patterns "
        f"(>= {threshold} same-fingerprint fails)",
    ]
    if patterns:
        for p in patterns:
            lines.append(f"- **{p.fingerprint}** — {p.count}x; e.g. "
                         f"\"{p.examples[0] if p.examples else ''}\"")
    else:
        lines.append("- none — no recurring failure above threshold")
    lines += ["", "## Benchmark suite integrity", f"- {integrity}", ""]
    # Telemetry consumption: surface the collected tool-trace / turn-metric
    # health plus benchmark drift, so the daily pass reads what is recorded.
    lines += _tool_health_lines() + [""]
    lines += _turn_health_lines() + [""]
    drift = _benchmark_drift_lines()
    if drift:
        lines += drift + [""]
    if drafts:
        lines.append("## Draft tasks scaffolded (review before committing)")
        lines += [f"- `{d}`" for d in drafts]
        lines.append("")
    lines += [
        "## Next steps",
        "- review drafts in `~/.delfin/bug_tasks/`, complete the TODO "
        "signals, move keepers into `delfin/agent/pack/benchmark/`",
        "- live benchmark run (costs tokens, run deliberately): "
        "`python -m delfin.agent.benchmark_runner`",
    ]
    return "\n".join(lines)


def run_eval(
    *,
    settings: dict | None = None,
    reports_dir: Path | None = None,
    drafts_dir: Path | None = None,
) -> Path:
    """One full LLM-free eval pass: mine → scaffold → report. Returns the
    report path."""
    cfg = eval_settings(settings)
    patterns = analyze_outcomes(window=cfg["window"], threshold=cfg["threshold"])
    drafts = write_pattern_drafts(patterns, dest_dir=drafts_dir) if patterns else []
    report = build_report(window=cfg["window"], threshold=cfg["threshold"],
                          drafts=drafts)
    d = reports_dir or _REPORTS_DIR
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"eval_{time.strftime('%Y%m%d')}.md"
    path.write_text(report, encoding="utf-8")
    return path


def maybe_run_scheduled(
    settings: dict | None = None,
    *,
    min_interval_hours: float = 24.0,
    reports_dir: Path | None = None,
) -> Path | None:
    """Opportunistic scheduled pass: run the LLM-free eval at most once per
    ``min_interval_hours``, gated by the same opt-in.

    Called from session-end paths (CLI run, dashboard) so the loop
    actually closes without anyone remembering a cron job — previously
    the module had NO automatic caller at all. Never raises; returns the
    report path when a pass ran, else None."""
    cfg = eval_settings(settings)
    if not cfg["enabled"]:
        return None
    d = reports_dir or _REPORTS_DIR
    stamp = d / ".last_run"
    now = time.time()
    try:
        last = float(stamp.read_text(encoding="utf-8").strip() or 0.0)
    except (OSError, ValueError):
        last = 0.0
    if now - last < min_interval_hours * 3600.0:
        return None
    try:
        path = run_eval(settings=settings, reports_dir=d)
    except Exception:
        return None
    try:
        d.mkdir(parents=True, exist_ok=True)
        stamp.write_text(str(now), encoding="utf-8")
    except OSError:
        pass
    return path


def main() -> int:
    cfg = eval_settings()
    if not cfg["enabled"]:
        print("eval_loop is disabled (agent.eval_loop.enabled=false). "
              "Enable it in settings to run scheduled evals — the default "
              "pass is LLM-free (0 tokens).")
        return 2
    path = run_eval()
    print(f"eval report written: {path}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
