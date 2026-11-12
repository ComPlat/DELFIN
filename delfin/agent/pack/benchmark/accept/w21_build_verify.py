#!/usr/bin/env python3
"""Acceptance for the wave-21 workflow-BUILDING tasks (package WB0).

The agent's deliverable is a DELFIN pipeline spec (spec.json) that
answers the chemist's brief (brief.txt), plus answer.txt where it
justifies every deviation from the defaults. This script grades the SPEC
with delfin/agent/workflow_grader.grade() — the same function the
platform exposes — and never runs a calculation.

Checks per task id:
  * <id>/spec.json exists and parses (a dict with "name" and "steps").
  * <id>/answer.txt exists (the justification text).
  * grade(spec, brief, answer) returns status "pass".
  * the spec uses at least one orca step (the request asks for physics
    an xtb-only chain cannot deliver) — enforced by the grader's
    coverage rules, not duplicated here.

usage: w21_build_verify.py <workspace> [<task_id>]
Without a task id (the runner's form, benchmark_runner.py:1349) the
task folder is the one under the workspace that holds a spec.json;
ambiguous when several do.
Exit 0 = accepted, 1 = failed, 2 = usage error.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_TASK_IDS = {
    "w21_build_conformers", "w21_build_redox", "w21_build_ts",
    "w21_build_spin", "w21_build_cycle", "w21_build_spectrum",
}


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def _resolve_task(ws: Path, task_id: str | None) -> Path:
    """The task folder, given the runner's one-argument form.

    benchmark_runner.run_acceptance (benchmark_runner.py:1349) launches
    this script with the workspace as the ONLY argument. The setup
    (setup/w21_build_fixtures.py) pre-creates one folder per task with
    only start.xyz/brief.txt; the attempt writes deliverables
    (spec.json, answer.txt, ...) only where the agent worked. With no
    task id, the folder holding a spec.json IS the task. The
    two-argument form is kept (extend, never lose).
    """
    if task_id is not None:
        if task_id not in _TASK_IDS:
            print(f"usage: unknown task id {task_id}", file=sys.stderr)
            sys.exit(2)
        return ws / task_id
    candidates = sorted(
        p.parent for p in ws.rglob("spec.json")
        if p.parent.name in _TASK_IDS)
    if not candidates:
        _fail("no task folder in the workspace holds a spec.json -- "
              "nothing to grade")
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        _fail(f"several task folders hold deliverables ({names}) -- "
              "ambiguous without a task id")
    return candidates[0]


def main(argv: list[str]) -> None:
    if len(argv) < 2 or len(argv) > 3:
        print("usage: w21_build_verify.py <workspace> [<task_id>]",
              file=sys.stderr)
        sys.exit(2)
    ws = Path(argv[1])
    task_id = argv[2] if len(argv) == 3 else None
    d = _resolve_task(ws, task_id)
    task_id = d.name
    if not d.is_dir():
        _fail(f"{d} missing")
    brief_p = d / "brief.txt"
    if not brief_p.is_file():
        _fail(f"{brief_p} missing — the grader needs the chemist's brief: "
              "request-dependent rules (conformers, relativity, coverage) "
              "grade against it")
    brief = brief_p.read_text(encoding="utf-8")
    spec_p = d / "spec.json"
    if not spec_p.is_file():
        _fail(f"{spec_p} missing — the agent built no workflow spec")
    try:
        spec = json.loads(spec_p.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _fail(f"{spec_p} unreadable: {e}")
    if not isinstance(spec, dict) or "steps" not in spec:
        _fail(f"{spec_p} is not a pipeline spec (needs name + steps)")

    answer_p = d / "answer.txt"
    if not answer_p.is_file():
        _fail(f"{answer_p} missing — unjustified deviations cannot be graded")
    answer = answer_p.read_text(encoding="utf-8")

    from delfin.agent import workflow_grader
    verdict = workflow_grader.grade(spec, brief, answer=answer)
    if verdict["status"] != "pass":
        print("verdict reasons:", file=sys.stderr)
        for r in verdict["reasons"]:
            print(f"  - {r}", file=sys.stderr)
        _fail(f"{task_id}: workflow_grader graded the spec fail")

    print(f"ok: {task_id}: spec graded pass "
          f"({len(verdict['findings'].get('lint', []))} lint findings, "
          f"estimate {verdict['findings']['cost']['core_hours_estimate']:.2f} "
          f"core-hours)")


if __name__ == "__main__":
    main(sys.argv)
