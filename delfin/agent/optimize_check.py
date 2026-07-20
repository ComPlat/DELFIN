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

import re
from dataclasses import dataclass


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


def check_benchmark_tasks() -> list[Issue]:
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
                issues.append(Issue("warn", tid, "expected signal still contains a TODO placeholder"))
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


def run_checks() -> list[Issue]:
    return (
        check_benchmark_tasks()
        + check_ground_truth()
        + check_agent_prompts()
    )


def main() -> int:
    issues = run_checks()
    errors = [i for i in issues if i.severity == "error"]
    warns = [i for i in issues if i.severity == "warn"]
    for i in issues:
        print(i)
    if errors:
        print(f"\n❌ {len(errors)} error(s), {len(warns)} warning(s) — NOT safe to ship.")
        return 1
    print(f"\n✅ Safe to optimise — 0 errors, {len(warns)} warning(s).")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
