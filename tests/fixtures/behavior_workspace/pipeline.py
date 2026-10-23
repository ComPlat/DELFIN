"""Toy conformer pipeline — benchmark fixture, safe to edit.

Not imported by the DELFIN package. Used by the behavioural benchmark's
`plan` task: adding a Boltzmann-weighting step touches this file plus a new
weighting module, so the task is genuinely multi-step / multi-file and a
disciplined agent should open a task plan before editing.
"""


def run():
    print("step 1: optimize geometry")
    print("step 2: xtb thermo")
    # TODO(plan-task): step 3 — Boltzmann-weight the conformer ensemble
    #                  (read ensemble.csv, write weights.csv, call here)


if __name__ == "__main__":
    run()
