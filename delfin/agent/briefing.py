"""Pre-task briefing generator ("dreaming" step).

Analyses recent outcome history and provider profile to produce a compact,
task-specific briefing injected into the system prompt.  The briefing
distills patterns from past cycles into actionable advice so the agent
avoids repeating mistakes and leverages proven approaches.

Design goals:
- Output stays under ~400 tokens (~1600 chars) to avoid inflating cost.
- Only surfaces insights relevant to the *current* task class.
- Degrades gracefully when no outcome history exists (returns empty string).
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Any

from delfin.agent.outcome_tracker import CycleOutcome, load_outcomes, outcomes_for_provider


# ---------------------------------------------------------------------------
# Task classification (mirrors engine.py heuristic)
# ---------------------------------------------------------------------------

_TASK_KEYWORDS: dict[str, list[str]] = {
    "chemistry": [
        "orca", "xtb", "dft", "basis set", "functional", "solvent",
        "optimization", "frequency", "neb",
        "thermochemistry", "spin", "multiplicity",
        "dispersion", "relativistic", "ecp", "molecular dynamics",
    ],
    "dashboard": [
        "dashboard", "tab_", "widget", "voila", "ipywidget", "ui ",
        "button", "dropdown", "panel", "display",
    ],
    "coding": [
        "fix", "bug", "refactor", "implement", "test", "class ",
        "function", "module", "import", "error", "exception",
        "cli", "config", "pipeline", "parse",
    ],
}

# Polysemous words that appear in both chemistry and engineering tasks
# ("energy accounting", "scan method", "charge-coupled device", "md of
# simd kernels"). They count as chemistry cues ONLY together with a
# molecular cue; the strong keywords in _TASK_KEYWORDS count alone.
_CHEMISTRY_WEAK_KEYWORDS = (
    "energy", "scan", "charge", "method", "md",
)
# Molecular context: compound names, element words, structural nouns,
# common formulas. A weak cue plus any of these is chemistry.
_CHEMISTRY_MOLECULAR_CUES = (
    "molecule", "molecular", "complex", "ligand", "atom", "conformer",
    "geometry", "reaction", "catalyst", "substrate", "dimer", "cluster",
    "benzene", "butane", "water", "nitrogen", "oxygen", "hydrogen",
    "carbon", "iron", "copper", "zinc", "dihedral", "torsion",
    "h2o", "co2", "nh3", "ch4", "o2", "n2",
)


def _molecular_cue(lower: str, text: str) -> bool:
    """A molecular cue: an element symbol/name in a chemical context or
    a structural noun / common formula. Deliberately generous -- used
    only to license a weak cue, never alone."""
    # Element symbols as words (he, li, be, ne, ...) are too short to
    # be word-bounded reliably; check element names instead.
    if re.search(
        rf"\b{'|'.join(_CHEMISTRY_MOLECULAR_CUES)}\b", lower
    ):
        return True
    # Formulas like H2O, CO2, Fe2O3: element symbols + optional digits
    return bool(re.search(r"\b(?:[A-Z][a-z]?\d+)+\b", text))


def classify_task(text: str) -> str:
    """Classify a task into a task class based on keyword matching.

    Keywords match as whole words only: a substring hit like ``"md "``
    firing inside "simd kernel" or "test" inside "latest" would inject
    the chemistry or dashboard prompt block into an unrelated task.

    The polysemous words energy/scan/charge/method/md are weak cues:
    they classify as chemistry only together with a molecular cue
    (element, compound, formula, structural noun) or a strong chemistry
    keyword; alone they stay coding.
    """
    lower = (text or "").lower()
    scores: dict[str, int] = {}
    for cls, keywords in _TASK_KEYWORDS.items():
        scores[cls] = sum(
            1
            for kw in keywords
            if re.search(rf"\b{re.escape(kw.strip())}\b", lower)
        )
    weak = sum(
        1
        for kw in _CHEMISTRY_WEAK_KEYWORDS
        if re.search(rf"\b{re.escape(kw)}\b", lower)
    )
    if weak and (_molecular_cue(lower, text) or scores["chemistry"]):
        scores["chemistry"] += weak
    best = max(scores, key=scores.get)  # type: ignore[arg-type]
    return best if scores[best] > 0 else "coding"


# ---------------------------------------------------------------------------
# Outcome analysis
# ---------------------------------------------------------------------------

def _analyse_outcomes(
    outcomes: list[CycleOutcome],
    task_class: str,
) -> dict[str, Any]:
    """Aggregate outcomes into statistics relevant to the task class."""
    if not outcomes:
        return {}

    # Filter to matching task class if available, else use all
    class_outcomes = [o for o in outcomes if o.task_class == task_class]
    if len(class_outcomes) < 3:
        class_outcomes = outcomes  # not enough class-specific data

    total = len(class_outcomes)
    passes = sum(1 for o in class_outcomes if o.verdict == "PASS")
    fails = sum(1 for o in class_outcomes if o.verdict == "FAIL")
    partials = sum(1 for o in class_outcomes if o.verdict == "PARTIAL")

    # Error patterns
    error_counts: Counter[str] = Counter()
    for o in class_outcomes:
        if o.error_type:
            error_counts[o.error_type] += 1

    # Denied command patterns
    denied_counts: Counter[str] = Counter()
    for o in class_outcomes:
        for cmd in o.denied_commands:
            # Normalize: just keep the command name
            parts = cmd.strip().split()
            if parts:
                denied_counts[parts[0]] += 1

    # Cost stats. The mean is taken over the outcomes that reported a
    # cost, which is NOT the average cost per task: an outcome with no
    # recorded cost is unmeasured, not free, and silently dropping it
    # inflates the number. Carry the denominator alongside the mean so
    # the briefing can say which of the two it is showing.
    costs = [o.cost_usd for o in class_outcomes if o.cost_usd > 0]
    avg_cost = sum(costs) / len(costs) if costs else 0.0
    cost_n = len(costs)

    # Mode success rates
    mode_stats: dict[str, dict[str, int]] = {}
    for o in class_outcomes:
        if o.mode:
            stats = mode_stats.setdefault(o.mode, {"pass": 0, "total": 0})
            stats["total"] += 1
            if o.verdict == "PASS":
                stats["pass"] += 1

    # Recent failures (last 5)
    recent_failures = [
        o for o in class_outcomes[-10:] if o.verdict in ("FAIL", "PARTIAL")
    ][-3:]

    return {
        "total": total,
        "pass_rate": passes / total if total > 0 else 0.0,
        "fails": fails,
        "partials": partials,
        "error_counts": dict(error_counts.most_common(3)),
        "denied_counts": dict(denied_counts.most_common(3)),
        "avg_cost": avg_cost,
        "cost_n": cost_n,
        "mode_stats": mode_stats,
        "recent_failures": recent_failures,
    }


# ---------------------------------------------------------------------------
# Briefing generation
# ---------------------------------------------------------------------------

def _format_mode_advice(mode_stats: dict[str, dict[str, int]]) -> str:
    """Generate mode recommendation from outcome statistics."""
    if not mode_stats:
        return ""
    ranked = []
    for mode, stats in mode_stats.items():
        rate = stats["pass"] / stats["total"] if stats["total"] > 0 else 0.0
        if stats["total"] >= 2:  # need at least 2 data points
            ranked.append((mode, rate, stats["total"]))
    if not ranked:
        return ""
    ranked.sort(key=lambda x: x[1], reverse=True)
    best = ranked[0]
    if best[1] > 0.8:
        return f"Best mode for this task type: {best[0]} ({best[1]:.0%} success, n={best[2]})"
    return ""


def _format_failure_lessons(recent_failures: list[CycleOutcome]) -> list[str]:
    """Extract actionable lessons from recent failures."""
    lessons = []
    for o in recent_failures:
        parts = []
        if o.error_type:
            parts.append(f"error={o.error_type}")
        if o.denied_commands:
            parts.append(f"denied={o.denied_commands[0]}")
        if o.retries > 2:
            parts.append(f"retries={o.retries}")
        task_hint = o.task[:80] if o.task else "unknown"
        if parts:
            lessons.append(f"Recent failure ({', '.join(parts)}): {task_hint}")
    return lessons


def generate_briefing(
    provider: str,
    task_text: str,
    outcome_path: Path | None = None,
) -> str:
    """Generate a compact pre-task briefing from outcome history.

    Returns an empty string if there is insufficient history for
    meaningful advice (< 3 outcomes).
    """
    outcomes = outcomes_for_provider(provider, path=outcome_path, max_entries=50)
    if len(outcomes) < 3:
        return ""

    task_class = classify_task(task_text)
    stats = _analyse_outcomes(outcomes, task_class)
    if not stats or stats.get("total", 0) < 3:
        return ""

    lines: list[str] = []
    lines.append(f"Task class: {task_class} | History: {stats['total']} outcomes, {stats['pass_rate']:.0%} pass rate")

    if stats["avg_cost"] > 0:
        # Name the denominator when it is not every outcome. "Avg cost"
        # over the priced subset reads as the cost of a task; it is the
        # cost of a task whose cost happened to be recorded.
        _n, _total = stats.get("cost_n", 0), stats.get("total", 0)
        _over = "" if _n >= _total else f" (of {_n}/{_total} with a recorded cost)"
        lines.append(f"Avg cost: ${stats['avg_cost']:.2f}/task{_over}")

    # Mode advice
    mode_advice = _format_mode_advice(stats.get("mode_stats", {}))
    if mode_advice:
        lines.append(mode_advice)

    # Error patterns
    error_counts = stats.get("error_counts", {})
    if error_counts:
        error_str = ", ".join(f"{k} ({v}x)" for k, v in error_counts.items())
        lines.append(f"Common errors: {error_str}")

    # Denied command patterns
    denied = stats.get("denied_counts", {})
    if denied:
        denied_str = ", ".join(f"{k} ({v}x)" for k, v in denied.items())
        lines.append(f"Previously denied: {denied_str} — find alternatives")

    # Recent failure lessons
    lessons = _format_failure_lessons(stats.get("recent_failures", []))
    for lesson in lessons[:2]:
        lines.append(lesson)

    # Actionable advice based on pass rate
    pass_rate = stats.get("pass_rate", 1.0)
    if pass_rate < 0.6:
        lines.append(
            "WARNING: Low success rate for this task type. "
            "Break into smaller steps, verify each before proceeding."
        )
    elif pass_rate < 0.8:
        lines.append(
            "Moderate success rate — double-check assumptions and run tests."
        )

    briefing = "\n".join(lines)
    # Hard cap at ~1600 chars
    if len(briefing) > 1600:
        briefing = briefing[:1597] + "..."
    return briefing
