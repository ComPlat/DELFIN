# MANTA split tooling

The MANTA constructor used to be one module (`delfin/smiles_converter.py`,
40 389 lines) plus `delfin/manta/assemble_complex.py` (6 833 lines). Both were
split into `delfin/manta/*.py` modules by pure moves: every statement keeps its
text, the new modules only gain the imports they need, and the old paths
re-export every moved name so existing imports keep working. The files here are
the instrument that proves the split changed nothing, and the scripts that did
the moving.

All commands use the interpreter DELFIN runs with (the one that has RDKit and
Open Babel) and are run from the repository root.

## 1. The identity harness

`split_identity.py` builds the SMILES in `split_identity_smiles.tsv` (70 public
examples from `examples/manta_constructions`, 14 from the repository tests) the
way the `delfin manta` CLI and the dashboard build them: the shipped champion
construction environment from `delfin.cli_manta.construction_env`, a fresh
subprocess per SMILES with `PYTHONHASHSEED=0`, the same four lines
`delfin.common.manta_build` runs, the CLI's default build arguments. For every
SMILES it hashes the multi-frame XYZ, the frame labels, the error string and the
FF-free ISO trace.

```bash
# compare the working tree against the recorded manifest (sha256 per artefact)
python tools/split_identity.py

# the quick set: only SMILES whose recorded build took at most 150 s
python tools/split_identity.py --quick 150

# compare the working tree against a reference commit, built on the fly
# (git archive of that commit into a temporary directory, never a checkout)
python tools/split_identity.py --ref-commit 96c2a60a

# the legacy path (no construction switches), tree vs reference commit
python tools/split_identity.py --ref-commit 96c2a60a --construction default

# record a new reference manifest from the working tree
python tools/split_identity.py --record
```

`--tree PATH` tests another checkout (it is put in front of `PYTHONPATH`),
`--only id1,id2` restricts the set, `--workers N` sets the parallelism,
`--keep-ref-outputs DIR` keeps the artefacts of a `--record` run so a later
compare can print the first differing line. The exit status is 0 only when
every SMILES is identical (`identity: N/N identical`).

The manifest `split_identity_manifest.json` was recorded at `origin/main`
96c2a60a. The second instrument is the MANTA test selection:

```bash
pytest -m "not slow" tests/test_cli_manta.py tests/test_cli_dashboard_parity.py \
    tests/test_user_smiles_suite.py tests/test_an_isomer_is_the_molecule_that_was_asked_for.py \
    tests/test_smiles_converter_regressions.py tests/test_shim_integrity.py \
    tests/test_no_regression_undefined_names.py tests/test_every_frame_carries_the_whole_smiles.py \
    tests/test_a_manta_setting_reaches_the_builder.py
```

## 2. How a module was moved

Three scripts, in this order:

1. `split_depgraph.py <module.py> <depgraph.json>` parses the monolith, lists
   every top-level definition with the module-level names it references, the
   strongly connected components of that graph and the author's section
   banners, and writes the graph as JSON.
2. `split_plan_smiles_converter.py <depgraph.json> <out_dir>` and
   `split_plan_assemble_complex.py <depgraph.json> <out_dir>` turn the group
   specification (line ranges of the original file plus explicit moves) into
   `<out_dir>/names/<module>.txt` and `<out_dir>/plan.json`. The groups are
   ordered so that every module only needs names that already left the
   monolith; `split_extract.py --plan plan.json` checks exactly that
   ("closure check") and prints every name a group would still need from the
   shim.
3. `split_extract.py --shim <monolith> --plan-step plan.json <module> --registry
   split_registry.json` performs one move: it extracts the group's units
   verbatim into `delfin/manta/<module>.py`, writes that module's imports
   (standard library and third-party lines taken from the monolith's own
   header, `delfin.manta` lines for names already extracted), replaces the
   units in the monolith by a re-export block, and records name -> module in the
   registry. `--logger-name delfin.smiles_converter` emits a logger with that
   name in every module (same log records as before); `--origin` is the path
   named in the module docstring; `--insert-after NAME` places the re-export
   block after that unit of the shim.

After each move: `ruff check --select F821,F811,F601,F402,F823` on the new
module and the shim, an import of `delfin.smiles_converter`, the identity
harness, the test selection, one commit.

Two things cannot move as text and were followed by hand, both in
`_smiles_to_xyz_isomers_impl` and both behind switches that are off in the
champion configuration: the write of the module global
`_ITER84_SIGMA_CAPS_OVERRIDE` (now set on `delfin.manta.chelate_templates`,
where its reader lives) and the temporary patch of
`_verify_topology_from_graph` through `sys.modules` (now aimed at
`delfin.manta.topo_isomers`, where the function and all its call-time readers
live).
