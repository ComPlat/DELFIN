#!/usr/bin/env python3
"""Acceptance for the wave-21 real-run workflow tasks (package WB0).

Three of the build tasks run FOR REAL (tiny, <= 20 atoms, through
DELFIN's submit path, inside the practice budget). The agent's
deliverable has two halves, and BOTH are graded:

1. the SPEC (spec.json + answer.txt) — graded with
   workflow_grader.grade() exactly like the build-only tasks
   (w21_build_verify.py);
2. the RUN EVIDENCE — graded with CB0A's whole-run interface:
   delfin.agent.run_evidence.collect(workspace/<task_id>) +
   accept/w20_common.grade_run(ev, max_cores=2). A run the agent
   claims but whose inputs/outputs do not show it fails; a run whose
   charge/multiplicity diverges from CONTROL fails; a run inside the
   core cap that terminated normally passes.

usage: w21_run_verify.py <workspace> [<task_id>]
Without a task id (the runner's form, benchmark_runner.py:1349) the
task folder is the one under the workspace that holds a spec.json;
ambiguous when several do.
Exit 0 = accepted, 1 = failed, 2 = usage error.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# task id -> (brief, required result.json keys)
_TASKS = {
    "w21_run_methanol": (
        "Optimise methanol in the gas phase and confirm the geometry is "
        "a true minimum; report the thermochemical correction at 298 K.",
        ("energy_hartree", "n_imag", "thermal_corr_298k"),
    ),
    "w21_run_water_ir": (
        "Compute the harmonic frequencies and IR intensities of water "
        "and report the frequency of the bending mode.",
        ("energy_hartree", "n_imag", "bend_freq_cm"),
    ),
    "w21_run_hcl_spectrum": (
        "Compute the infrared spectrum of hydrogen chloride and report "
        "the stretch frequency.",
        ("energy_hartree", "n_imag", "stretch_freq_cm"),
    ),
}


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def _resolve_task(ws: Path, task_id: str | None) -> tuple[str, Path]:
    """The task id and folder, given the runner's one-argument form.

    benchmark_runner.run_acceptance (benchmark_runner.py:1349) launches
    this script with the workspace as the ONLY argument. The setup
    (setup/w21_run_fixtures.py) pre-creates one folder per task with
    only start.xyz/brief.txt; the attempt writes deliverables
    (spec.json, run_evidence.json, ...) only where the agent worked.
    With no task id, the folder holding a spec.json IS the task. The
    two-argument form is kept (extend, never lose).
    """
    if task_id is not None:
        if task_id not in _TASKS:
            print(f"usage: unknown task id {task_id}", file=sys.stderr)
            sys.exit(2)
        return task_id, ws / task_id
    candidates = sorted(
        p.parent for p in ws.rglob("spec.json")
        if p.parent.name in _TASKS)
    if not candidates:
        _fail("no task folder in the workspace holds a spec.json -- "
              "nothing to grade")
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        _fail(f"several task folders hold deliverables ({names}) -- "
              "ambiguous without a task id")
    return candidates[0].name, candidates[0]


def main(argv: list[str]) -> None:
    if len(argv) < 2 or len(argv) > 3:
        print("usage: w21_run_verify.py <workspace> [<task_id>]",
              file=sys.stderr)
        sys.exit(2)
    ws = Path(argv[1])
    task_id = argv[2] if len(argv) == 3 else None
    task_id, d = _resolve_task(ws, task_id)
    brief, required = _TASKS[task_id]

    if not d.is_dir():
        _fail(f"{d} missing")

    # --- half 1: the SPEC, through the real grader -----------------------
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
        print("spec verdict reasons:", file=sys.stderr)
        for r in verdict["reasons"]:
            print(f"  - {r}", file=sys.stderr)
        _fail(f"{task_id}: workflow_grader graded the spec fail")

    # --- half 2: the RUN EVIDENCE, through CB0A's interface --------------
    from delfin.agent import run_evidence
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = run_evidence.collect(d)
    run_verdict = w20_common.grade_run(ev, max_cores=2)
    if run_verdict["verdict"] != "pass":
        print("run evidence problems:", file=sys.stderr)
        for p in run_verdict["problems"]:
            print(f"  - {p}", file=sys.stderr)
        _fail(f"{task_id}: the run's evidence does not show a real, "
              "terminated, in-budget run")

    # --- the reported numbers exist in the run's own results -------------
    rj = d / "result.json"
    if not rj.is_file():
        _fail(f"{rj} missing — the run reported no numbers")
    try:
        data = json.loads(rj.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _fail(f"{rj} unreadable: {e}")
    if not isinstance(data, dict):
        _fail(f"{rj} is not an object")
    for key in required:
        if key not in data or data[key] is None:
            _fail(f"{rj} lacks required field {key}")
        if isinstance(data[key], str):
            _fail(f"{key} must be a number with its unit in the key, "
                  f"got prose: {data[key]!r}")

    # reported numbers must come from the outputs the run produced
    for key in required:
        if "energy" in key or "corr" in key:
            if not w20_common.reported_value_matches_output(
                    ev, key, kind="energy"):
                _fail(f"{key} is not traceable to the run's outputs")

    if "n_imag" in data and float(data["n_imag"]) != 0.0:
        _fail(f"claimed a minimum but n_imag={data['n_imag']}")

    print(f"ok: {task_id}: spec graded pass AND run evidence graded pass "
          f"({len(ev.steps)} steps collected, "
          f"estimate {verdict['findings']['cost']['core_hours_estimate']:.2f} "
          f"core-hours)")


if __name__ == "__main__":
    main(sys.argv)
