# Offline biomolecular foundation example

This synthetic I/O example runs ordinary DELFIN adapters through the normal
pipeline: structure summary → ligand preparation → pose RMSD. The ligand
artifact is wired to the RMSD reference input by the pipeline. It requires
DELFIN and RDKit; no predictor, docking executable, network or model download
is involved.

From the repository root, using the Python environment that contains DELFIN:

```sh
python -m examples.biomolecular_foundation.run_example --output /tmp/bio-example
```

Choose an empty output directory. Inputs are never overwritten. The output
contains the step artifacts and `example_report.json`, including input hashes
and relative artifact paths.

The receptor is an incomplete synthetic backbone fragment, not an experimental
structure or a scientifically usable docking receptor. The reference is ethane;
the supplied poses translate it by 0, 1 and 3 Å. Heavy-atom RMSD without alignment
should reproduce those distances, with two poses under the default 2 Å cutoff.
The `DEMO_SCORE` values are arbitrary demonstration numbers, not energies or
affinities. They deliberately rank the farthest pose first: score ranking and
reference-pose agreement are different questions.

This checks file semantics, discovery, wiring and reporting. It does not validate
a live docking engine, receptor preparation or a physical model.
