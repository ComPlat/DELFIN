# Science-benchmark workspace

Everyday computational-chemistry material for the `science_*` benchmark
tasks: a conformer ensemble and four xtb outputs.

The files carry the things real output carries and demonstrations leave
out:

* `ensemble.csv` is unsorted, and `c6` repeats `c3`'s energy exactly --
  a duplicate geometry, which changes a population sum if it is counted
  twice and not if it is recognised.
* `run_d.out` did NOT converge (`abnormal termination`). It still holds a
  total energy and a gap, so anything that reads the numbers without
  reading the termination line will happily average it in.

Nothing here is imported by the DELFIN package; it is fixture data, reset
between replicates by the benchmark's workspace guard.

## Constants and conditions

Stated here because a result that cannot be reproduced is not a result,
and the integrity rules this fixture is built to exercise say so: methods
belong with results.

* Energies in `ensemble.csv` are electronic energies in **hartree**.
* Boltzmann populations are at **298.15 K**.
* 1 hartree = **2625.4996392852 kJ/mol**; R = **8.314462618e-3
  kJ/(mol·K)**, so RT = 2.478960 kJ/mol at that temperature.
* `run_*.out` carry the total energy in **Eh** and the HOMO-LUMO gap in
  **eV**, as xtb writes them.

Reported by kit.deepseek-v4-flash, asked what it would change after
building the weighting step: the temperature and the conversion factor
lived only inside the code that used them, so a reader could not check a
population without reading the implementation first.
