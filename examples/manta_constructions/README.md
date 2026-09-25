# MANTA construction examples

Multi-frame XYZ files produced by DELFIN's MANTA construction from a SMILES, exactly as the
current `main` builds them with the shipped champion settings (the same construction the
`delfin manta` CLI and the dashboard use). No geometry optimisation of any kind was applied:
every frame is the deterministic construction output, no force field, no xTB, no DFT.

## What is in here

* `index.tsv` — one row per example: `id`, the input `smiles`, `n_atoms` (heavy atoms + H),
  `n_frames` (isomers / conformers in the manifold).
* `structures/<id>.xyz` — the manifold of that example as concatenated XYZ frames. The comment
  line of every frame carries the example id, the frame index and DELFIN's own frame label
  (coordination polyhedron, isomer or conformer tag).

## How the examples were chosen

* Metal complexes and organics with more than 30 atoms.
* Only systems whose every frame passes DELFIN's internal geometry battery without a hard
  finding (no torn or collapsed bond, no clash, no broken hybridisation, no misplaced
  hydrogen). The selection tells nothing about how many systems fail that bar.
* At most 1000 examples; this release holds the number listed in `index.tsv`.

## What is NOT in here

No crystal structures, no crystal coordinates and no numbers derived from crystal data. The
SMILES is the only input; the frames are DELFIN's own construction. The examples are
regenerated with every landing of a new construction on `main`.

## Reproduce one example

```bash
delfin-manta "<smiles from index.tsv>" -o my_example
```

The construction config defaults to `champion` and no optimisation is run. The output is
byte-identical to the shipped file when run on the same `main` commit (the
construction is deterministic).
