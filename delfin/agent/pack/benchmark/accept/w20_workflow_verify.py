#!/usr/bin/env python3
"""Acceptance for the wave-20 workflow tasks (package CB0B).

One script, keyed by task id (passed as the molecule key slot). It
grades the INTERNAL-CONSISTENCY level offline: the workspace must
carry the numbers the task asked for, written by the agent's run, and
-- for a trapped task -- the answer must name the designed pitfall.
It never re-runs ORCA/xtb (level-b reproduction is the runner's job
with the validated defaults; this script only checks the reported
evidence is complete and consistent).

Checks per task id:
  * <id>/result.json exists and parses, with the numeric fields the
    prompt named (energies carry their units as string suffixes).
  * <id>/answer.txt names the trap when the task carries one (the
    trap text's first keyword must appear).
  * geometry claims come with a frequency count (minimum tasks).

usage: w20_workflow_verify.py <workspace> [<task_id>]
Without a task id (the runner's form, benchmark_runner.py:1349) the
task folder is the one under the workspace that holds a result.json;
ambiguous when several do.
Exit 0 = accepted, 1 = failed, 2 = usage error.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# task id -> (required result.json keys, trap keyword or '')
TASKS = {
    "w20_classic_opt_freq_1": (["energy_hartree", "n_imag", "thermal_corr_298k"], "saddle"),
    "w20_redox_1": (["energy_oxidized_hartree", "energy_reduced_hartree", "potential_v"], ""),
    "w20_excited_states_1": (["wavelength_nm", "oscillator_strength"], "visible"),
    "w20_imaginary_mode_1": (["n_imag", "energy_change_hartree"], "saddle"),
    "w20_spin_states_1": (["relative_energies_hartree", "s2_values"], "contamination"),
    "w20_conformer_prestep_1": (["conformer_spread_hartree", "lowest_smiles"], ""),
    "w20_screening_1": (["solvation_free_energy_kj_mol"], ""),
}


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def _resolve_task(ws: Path, task_id: str | None) -> tuple[Path, str]:
    """The task folder and id, given the runner's one-argument form.

    benchmark_runner.run_acceptance (benchmark_runner.py:1349) launches
    this script with the workspace as the ONLY argument. The setup
    (setup/w20_rest.py, setup/w20_classic_opt_freq.py) pre-creates one
    folder per task with only start.xyz/brief.txt; the attempt writes
    deliverables (result.json, answer.txt, ...) only where the agent
    worked. So with no task id, the folder that holds a deliverable IS
    the task. Kept for the two-argument form (extend, never lose).
    """
    if task_id is not None:
        if task_id not in TASKS:
            _fail(f"unknown task id {task_id}")
        return ws / task_id, task_id
    candidates = sorted(
        p.parent for p in ws.rglob("result.json")
        if p.parent.name in TASKS)
    if not candidates:
        _fail("no task folder in the workspace holds a result.json -- "
              "nothing to grade")
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        _fail(f"several task folders hold deliverables ({names}) -- "
              "ambiguous without a task id")
    return candidates[0], candidates[0].name


def main(argv: list[str]) -> None:
    if len(argv) < 2 or len(argv) > 3:
        print("usage: w20_workflow_verify.py <workspace> [<task_id>]",
              file=sys.stderr)
        sys.exit(2)
    ws = Path(argv[1])
    task_id = argv[2] if len(argv) == 3 else None
    d, task_id = _resolve_task(ws, task_id)
    required, trap_key = TASKS[task_id]
    if not d.is_dir():
        _fail(f"{d} missing")
    rj = d / "result.json"
    if not rj.is_file():
        _fail(f"{rj} missing -- the run reported no numbers")
    try:
        data = json.loads(rj.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        _fail(f"{rj} unreadable: {e}")
    if not isinstance(data, dict):
        _fail(f"{rj} is not an object")
    for key in required:
        if key not in data or data[key] is None:
            _fail(f"{rj} lacks required field {key}")
        val = data[key]
        if isinstance(val, str):
            _fail(f"{key} must be a number with its unit in the key, got prose: {val!r}")

    answer_p = d / "answer.txt"
    answer = answer_p.read_text(encoding="utf-8").lower() if answer_p.is_file() else ""
    if trap_key and trap_key not in answer:
        _fail(f"trap not reported: expected '{trap_key}' mentioned in answer.txt")

    if "n_imag" in required and float(data["n_imag"]) != 0.0:
        _fail(f"claimed a minimum but n_imag={data['n_imag']}")

    print(f"ok: {task_id} internally consistent ({len(required)} fields, "
          f"trap {'checked' if trap_key else 'n/a'})")


if __name__ == "__main__":
    main(sys.argv)
