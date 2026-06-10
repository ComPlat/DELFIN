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
- The entrypoint is gated by ``agent.eval_loop.enabled`` (default False)
  so nothing runs on anyone's machine without an explicit opt-in.

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
        "enabled": bool(cfg.get("enabled", False)),
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
