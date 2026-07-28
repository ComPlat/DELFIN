"""Safe-optimisation pre-flight check for the DELFIN agent.

Optimising the agent means editing prompts, benchmark tasks, or
ground-truth — all of which can silently break the optimisation set
itself (a bad regex, a duplicate task id, a curated keyword that drifted
out of the manual).  This module is the cheap, deterministic safety net
to run *before committing* such a change:

    python -m delfin.agent.optimize_check

It needs no live model and no network.  Exit code 0 = safe to ship,
1 = problems found (printed).  Checks:

1. **Benchmark integrity** — every task loads, has a unique id + a prompt,
   a known mode, at least one expected signal, and every signal pattern
   compiles as a regex with a valid ``against`` channel.
2. **Ground-truth integrity** — the curated fact-verify keyword pools are
   still backed by the extracted manual (reuses ``generate_fact_tasks``'s
   build-fail validation), so author memory can't drift from the manual.
3. **Agent prompts present** — the role prompts the modes depend on exist
   and are non-empty.

This makes the regression discipline a one-command gate for anyone
optimising the agent.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass
from typing import Mapping


_KNOWN_MODES = {
    "dashboard", "solo", "research", "quick", "reviewed",
    "cluster", "tdd", "full",
}
_VALID_AGAINST = {"text", "action", "tool_name", "any"}
_REQUIRED_AGENT_PROMPTS = ("solo_agent.md", "dashboard_agent.md")


@dataclass
class Issue:
    severity: str   # "error" | "warn"
    where: str
    message: str

    def __str__(self) -> str:
        icon = "❌" if self.severity == "error" else "⚠️"
        return f"{icon} [{self.where}] {self.message}"


def strict_from_env(env: Mapping[str, str] | None = None) -> bool:
    """True when ``DELFIN_OPTIMIZE_STRICT`` is set to a truthy value
    (1/true/yes/on, case-insensitive)."""
    src = os.environ if env is None else env
    return (str(src.get("DELFIN_OPTIMIZE_STRICT", "")).strip().lower()
            in {"1", "true", "yes", "on"})


def check_benchmark_tasks(*, strict: bool = False) -> list[Issue]:
    issues: list[Issue] = []
    try:
        from delfin.agent.benchmark import load_tasks
        tasks = load_tasks()
    except Exception as exc:
        return [Issue("error", "benchmark", f"load_tasks() failed: {exc}")]

    if not tasks:
        return [Issue("error", "benchmark", "no tasks loaded")]

    seen_ids: set[str] = set()
    for t in tasks:
        tid = getattr(t, "id", "") or "<no-id>"
        if tid in seen_ids:
            issues.append(Issue("error", tid, "duplicate task id"))
        seen_ids.add(tid)
        if not (getattr(t, "prompt", "") or "").strip():
            issues.append(Issue("error", tid, "empty prompt"))
        mode = getattr(t, "mode", "")
        if mode not in _KNOWN_MODES:
            issues.append(Issue("warn", tid, f"unknown mode '{mode}'"))
        exp = list(getattr(t, "expected_signals", []) or [])
        if not exp:
            issues.append(Issue("warn", tid, "no expected_signals (task can't pass)"))
        for sig in exp + list(getattr(t, "forbidden_signals", []) or []):
            against = getattr(sig, "against", "any")
            if against not in _VALID_AGAINST:
                issues.append(Issue("error", tid, f"invalid against='{against}'"))
            pat = getattr(sig, "pattern", "")
            try:
                re.compile(pat)
            except re.error as exc:
                issues.append(Issue("error", tid, f"bad regex {pat!r}: {exc}"))
            if "TODO" in pat:
                # In strict mode (A/B gate) an unfilled ground truth is a
                # FAILURE — the loop must not run against placeholder tasks.
                sev = "error" if strict else "warn"
                msg = "expected signal still contains a TODO placeholder"
                if strict:
                    msg += " (strict mode: fill in the ground truth first)"
                issues.append(Issue(sev, tid, msg))
    return issues


def check_ground_truth() -> list[Issue]:
    issues: list[Issue] = []
    try:
        from delfin.agent import generate_fact_tasks as g
    except Exception as exc:
        return [Issue("error", "ground-truth", f"import failed: {exc}")]
    for program in g._PROGRAM_BLOCK_TESTS:
        try:
            # Raises if a curated keyword drifted out of the manual, or a
            # forbid entry is actually a real manual keyword.
            g.build_tasks_for(program)
        except Exception as exc:
            issues.append(Issue("error", f"ground-truth/{program}", str(exc)))
    return issues


def check_agent_prompts() -> list[Issue]:
    issues: list[Issue] = []
    from pathlib import Path
    base = Path(__file__).resolve().parent / "pack" / "agents"
    for name in _REQUIRED_AGENT_PROMPTS:
        p = base / name
        if not p.is_file():
            issues.append(Issue("error", "prompts", f"missing {name}"))
        elif not p.read_text(encoding="utf-8").strip():
            issues.append(Issue("error", "prompts", f"empty {name}"))
    return issues


def run_checks(*, strict: bool = False) -> list[Issue]:
    """All checks. ``strict=True`` escalates TODO-placeholder benchmark
    tasks from warning to error (the A/B-loop gate)."""
    return (
        check_benchmark_tasks(strict=strict)
        + check_ground_truth()
        + check_agent_prompts()
    )


def main(argv: list[str] | None = None) -> int:
    args = list(argv if argv is not None else sys.argv[1:])
    strict = strict_from_env() or ("--strict" in args)
    issues = run_checks(strict=strict)
    errors = [i for i in issues if i.severity == "error"]
    warns = [i for i in issues if i.severity == "warn"]
    mode_note = " (strict mode)" if strict else ""
    for i in issues:
        print(i)
    if errors:
        print(f"\n❌ {len(errors)} error(s), {len(warns)} warning(s) — "
              f"NOT safe to ship.{mode_note}")
        return 1
    print(f"\n✅ Safe to optimise — 0 errors, {len(warns)} warning(s)."
          f"{mode_note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
