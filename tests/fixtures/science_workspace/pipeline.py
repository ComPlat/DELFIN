"""Conformer pipeline — fixture, safe to edit.

Not imported by the DELFIN package. Two steps exist; the third is the
one a real ensemble needs and is deliberately absent.
"""


def step_optimize():
    print("step 1: optimize geometry")


def step_thermo():
    print("step 2: xtb thermo")


def run():
    step_optimize()
    step_thermo()
    # TODO: step 3 — Boltzmann-weight the ensemble.
    #   read ensemble.csv, write weights.csv (conformer, population),
    #   populations at 298.15 K, and call it from here.


if __name__ == "__main__":
    run()
