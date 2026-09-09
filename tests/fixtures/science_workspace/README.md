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
