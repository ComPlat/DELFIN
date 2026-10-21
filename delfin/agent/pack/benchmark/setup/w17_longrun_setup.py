#!/usr/bin/env python3
"""Setup for the four wave-17 long-run fixtures (package KG4).

Builds four sibling fixture folders under <workspace>/w17_longrun/:

* lr_wait_for_job/    -- one slow job whose result must be collected
  (slow_report.py, sleeps a FIXED constant -- never a machine-speed
  probe -- then writes result.json). The idle/waiting rubric is read
  from the session's metrics records, not from a wall clock.
* lr_steered_goal/    -- a two-part task where the person steers the
  goal mid-way: the README asks for the FIRST part, and the scripted
  person message changes the goal to the second. The steering message
  lives in the fixture folder as STEERING_MESSAGE.txt for the accept
  script to check it was followed; the benchmark harness injects it.
* lr_forced_compaction/ -- a task with a large notes corpus to digest,
  so the session compacts mid-way; the digest must survive it.
* lr_restart_midway/  -- a task whose checkpoint file (state.json)
  marks where the work stopped; the restart must continue from it.

Exit non-zero on any failure: a task whose precondition was not built
is reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_SAMPLES = "index,enthalpy_kj_per_mol\n1,42.31\n2,42.05\n3,42.18\n" \
           "4,42.37\n5,41.98\n6,42.25\n"
_SLEEP_S = 20.0

_SLOW_SCRIPT = '''"""Compute the mean enthalpy of samples.csv -- slowly, on purpose.

Sleeps 20 s, then writes result.json. Run it in the background and
collect the result when it lands.
"""
import csv
import json
import time

TIME_S = {sleep:.1f}


def main():
    time.sleep(TIME_S)
    with open("samples.csv", newline="", encoding="utf-8") as fh:
        vals = [float(r["enthalpy_kj_per_mol"]) for r in csv.DictReader(fh)]
    mean = sum(vals) / len(vals)
    with open("result.json", "w", encoding="utf-8") as fh:
        json.dump({{"n": len(vals), "mean_kj_per_mol": round(mean, 2)}}, fh)
    print(f"{{len(vals)}} samples, mean {{mean:.2f}} kJ/mol")


if __name__ == "__main__":
    main()
'''

_STEERING_TXT = (
    "PERSON (mid-task steering, new goal): forget the summary part -- the "
    "goal is now only this: read catalog.txt and write ONE line per entry "
    "to answer.txt in the form '<letter>: <value> <unit>'; entries without "
    "a value are written as '<letter>: absent'. Nothing else is wanted.\n"
)

_CATALOG = "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent\n"

# -- lr_forced_compaction: a corpus large enough that a session digesting
#    it compacts mid-way. Built at runtime (n_blocks blocks), never a
#    hard-coded blob; the task is to produce digest.md summarising the
#    decisions recorded across the notes.
_NOTES_TMPL = ("## Note {i}\n\nDecision {i}: the cutoff for batch {i} is "
               "{cut:.2f} eV. Rationale: derived from run {i} of the "
               "screening series (see batch_{i}.dat).\n")

# -- lr_restart_midway: a checkpoint that marks work already done, so a
#    restarted session continues instead of starting over. state.json
#    records the done/remaining split; accept checks the remaining work
#    is finished AND the already-done part is not redone.
_CHECKPOINT = {
    "done": ["alpha", "beta"],
    "remaining": ["gamma", "delta"],
    "outputs_done": {
        "alpha": "2.10 kJ/mol",
        "beta": "5.80 kJ/mol",
    },
}

_README_RESTART = """# Restart task

state.json marks the work already done (alpha, beta with their
energies) and what remains (gamma, delta). Continue the run from the
checkpoint: write energy.txt with one line per molecule
'<name>: <energy>', alpha and beta with the recorded values, gamma and
delta computed from energies.csv by averaging its two columns.

Do NOT redo alpha or beta computations -- their values are recorded.
"""


def _mean_expected() -> float:
    rows = [ln.split(",") for ln in _SAMPLES.strip().splitlines()[1:]]
    return sum(float(r[1]) for r in rows) / len(rows)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w17_longrun_setup.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "w17_longrun"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    # Fixture value is re-checked against the arithmetic before the
    # fixture is accepted -- the graded value is measured, not invented.
    mean = _mean_expected()
    if abs(mean - 42.19) > 5e-3:
        print(f"fixture mean {mean:.4f} kJ/mol, expected 42.19",
              file=sys.stderr)
        return 1

    wait_dir = root / "lr_wait_for_job"
    wait_dir.mkdir(parents=True)
    (wait_dir / "samples.csv").write_text(_SAMPLES, encoding="utf-8")
    (wait_dir / "slow_report.py").write_text(
        _SLOW_SCRIPT.format(sleep=_SLEEP_S), encoding="utf-8")

    steer_dir = root / "lr_steered_goal"
    steer_dir.mkdir()
    (steer_dir / "catalog.txt").write_text(_CATALOG, encoding="utf-8")
    (steer_dir / "STEERING_MESSAGE.txt").write_text(
        _STEERING_TXT, encoding="utf-8")
    (steer_dir / "README.txt").write_text(
        "Task (original goal): summarise catalog.txt into summary.txt "
        "(one paragraph).\n"
        "NOTE: the person may steer the goal mid-way; a steering message "
        "replaces this goal for the rest of the run.\n",
        encoding="utf-8")

    comp_dir = root / "lr_forced_compaction"
    comp_dir.mkdir()
    notes = []
    for i in range(1, 41):
        notes.append(_NOTES_TMPL.format(
            i=i, cut=0.5 + 0.05 * (i % 7)))
    (comp_dir / "notes.md").write_text("\n".join(notes), encoding="utf-8")
    (comp_dir / "README.txt").write_text(
        "Task: read notes.md and write digest.md, one line per decision: "
        "'Decision <n>: cutoff <x> eV'. The notes are long; keep the "
        "digest complete even if your context is compacted mid-way.\n",
        encoding="utf-8")

    restart_dir = root / "lr_restart_midway"
    restart_dir.mkdir()
    (restart_dir / "state.json").write_text(
        json.dumps(_CHECKPOINT, indent=2) + "\n", encoding="utf-8")
    (restart_dir / "energies.csv").write_text(
        "molecule,col1,col2\ngamma,1.20,1.30\ndelta,3.00,3.20\n",
        encoding="utf-8")
    (restart_dir / "README.txt").write_text(
        _README_RESTART, encoding="utf-8")

    print(f"prepared {root} -- lr_wait_for_job (mean {mean:.2f} kJ/mol), "
          f"lr_steered_goal, lr_forced_compaction (40 notes), "
          f"lr_restart_midway")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
