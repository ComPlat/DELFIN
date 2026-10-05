# Construction batches: MANTA, Architector, molSimplify and epic-MACE on many SMILES

`delfin cluster` builds 3D structures for a list of `ID;SMILES` lines -- from a handful to
hundreds of thousands -- with one of four constructors:

| Tool | What is built | Runs in |
|---|---|---|
| `manta` | the shipped MANTA construction: every isomer and conformer frame, deterministic | DELFIN's own environment |
| `architector` | [Architector](https://github.com/lanl/Architector) (`full`: 10 symmetries x 10 conformers; `default`: its defaults) | its own environment |
| `molsimplify` | [molSimplify](https://github.com/hjkgrp/molSimplify), one structure per geometry of the coordination number | its own environment |
| `mace` | [epic-MACE](https://github.com/EPiCs-group/epic-mace) (`paper`: OH / SP; `extended`: + SPY, TBP, TET, SAN), every stereomer, ten conformers each | its own environment (Python 3.7) |

The list is cut into shards; one shard is one Slurm array task on one node (or one local run on a
workstation). A task that dies keeps every finished system and continues when it is resubmitted.
The result is one archive in DELFIN's multi-frame xyz format.

The external builders are never imported into DELFIN: each runs as a separate process started
with the interpreter of its own environment (`delfin/cluster_bench/tool_workers/`, which import
nothing from DELFIN). Nothing in this repository is input data -- lists, selections and results
live in your run directory.

## 1. Install

DELFIN itself (Python >= 3.10). For byte-for-byte reproducible MANTA builds, use the pinned
versions of `env/reference-environment.lock`:

```bash
git clone https://github.com/ComPlat/DELFIN.git
python -m venv delfin_env            # or micromamba / conda
delfin_env/bin/pip install -e ./DELFIN -c ./DELFIN/env/reference-environment.lock
delfin_env/bin/python -m delfin cluster --help
```

The external builders (only the ones you need):

```bash
# epic-MACE: Python 3.7 + RDKit 2020.09 + epic-MACE at the pinned commit efb5778e (GPL-3.0),
# built with micromamba under ~/.delfin/ai_tools/.mamba_env/epic_mace -- found automatically
delfin_env/bin/python -m delfin.installer --install epic-mace

# Architector and molSimplify: either into DELFIN's environment ...
delfin_env/bin/pip install 'delfin-complat[ai-complex]'
# ... or, recommended for comparisons, into an environment of their own with fixed versions:
micromamba create -y -p ./bench_env -c conda-forge python=3.11 rdkit openbabel numba
./bench_env/bin/pip install architector==0.0.10 molSimplify==2.0.0
```

Which interpreter a tool runs in: `--tool-python` if given, else `DELFIN_ARCHITECTOR_PYTHON` /
`DELFIN_MOLSIMPLIFY_PYTHON` / `DELFIN_MACE_PYTHON`, else the environment DELFIN's installer built
for the tool (epic-MACE), else DELFIN's own interpreter. `prepare` refuses an interpreter that does
not have the tool, before anything is written.

The same tools are available interactively in the dashboard (MANTA / ARCHITECTOR / MOLSIMPLIFY /
MACE buttons) and in CONTROL (`smiles_converter=MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE`).

## 2. Prepare the list

One system per line, `ID;SMILES` (`ID|SMILES` is accepted too):

```text
cisplatin;[Cl][Pt-2]([Cl])([NH3+])[NH3+]
hexaaqua_fe;[OH2+][Fe-3]([OH2+])([OH2+])([OH2+])([OH2+])[OH2+]
```

IDs: letters, digits and `_ . + -`, unique. A single bad line aborts `prepare` and writes nothing.
An optional selection file (`--select`) lists IDs, one per line; only those are built, in that
order.

## 3. Run on a Slurm cluster

```bash
PY=delfin_env/bin/python

# 1. run directory: shards, per-shard specs (external builders) and manifest.json
$PY -m delfin cluster prepare --tool mace --mode paper --input list.txt --run-dir runs/mace_list

# 2. a pilot of two shards in its own output folder, then the whole run
$PY -m delfin cluster slurm runs/mace_list --run pilot --array 0-1 --submit
$PY -m delfin cluster slurm runs/mace_list --throttle 20 --submit

# 3. progress; incomplete shards are printed as an --array value
$PY -m delfin cluster status runs/mace_list -v

# 4. merge into one archive
$PY -m delfin cluster collect runs/mace_list
```

`prepare` options: `--mode` (MANTA `champion`/`builder`, Architector `full`/`default`, MACE
`paper`/`extended`), `--shard-size` (MANTA 500, others 250), `--timeout` (per system, seconds,
default 21600), `--speed-factor`, `--workers` (concurrent builds per node: MANTA 36 x 7 threads,
others 48), `--threads`, `--repeat N`, `--specs`, `--label`. A run directory is never overwritten.

`slurm` options: `--time` (job wall time, default `72:00:00`), `--cpus` (48), `--mem` (MANTA 72G,
others 80G), `--throttle N` (`%N`, array tasks at once), `--array`, `--partition`, `--account`,
`--setup 'module load ...'` (repeatable), `--python` (DELFIN's interpreter on the nodes),
`--speed-factor`, `--submit`. Without `--partition`, `--submit` asks `sbatch --test-only` which of
the configured partitions (`DELFIN_SLURM_PARTITIONS`) accept the job. Without `--submit` only
`RUN/slurm/<tool>_<set>_<run>.sbatch` is written; submit it with `sbatch`.

The defaults fit a 48-core node with a 72 h limit (JUSTUS 2 is one such cluster). On another
cluster set `--cpus`, `--mem`, `--time` and `--workers` to its nodes. Put the run directory on a
file system the compute nodes see (a workspace, not `$HOME` if that is small);
`DELFIN_CLUSTER_WORK_ROOT=$TMPDIR` moves the tools' scratch directories to node-local disk.

**Without Slurm** (a workstation), build the shards one after the other:

```bash
for k in $(seq 0 $((N_SHARDS - 1))); do $PY -m delfin cluster run-shard runs/mace_list --shard $k; done
```

## 4. From the dashboard

Submit Job tab, panel **Construction batch (MANTA / ARCHITECTOR / MOLSIMPLIFY / MACE)** under the
Batch SMILES field. Builder, Mode, Run name and List file are always shown; everything else sits
in a folded **Advanced** section. Every field has a tooltip and a grey one-line explanation, and
shows its real default (no hidden 0 or empty value):

| Field | Default shown | `delfin cluster` argument |
|---|---|---|
| Builder | MANTA | `prepare --tool` |
| Mode | the builder's default, marked "(default)", with one line on what it builds | `--mode` |
| Run name | empty = `<builder>_<YYYYMMDD>_<N>mol` (shown as the placeholder, `_2`, `_3` when taken); a warning when the folder already exists | run directory `<calculation folder>/construction_batch/<run name>` (`--run-dir`) |
| List file | empty = the Batch SMILES field above, `name;SMILES[;...]` lines cut to `name;SMILES` and saved as `<run name>.input.txt` beside the run directory | `--input` |
| *Advanced:* Only these IDs (file), Reuse specs (JSONL file) | empty = all molecules / specs made from the SMILES | `--select`, `--specs` |
| Builder environment (python) | found automatically (shown as the placeholder) | `--tool-python` |
| Molecules per job (shard size) | the builder's default (MANTA 500, others 250) | `--shard-size` |
| Time limit per molecule (s) | 21600, shown as hours and with the slowness factor applied | `--timeout` |
| Cluster slowness factor (x) | 1.0 | `--speed-factor` |
| Max. parallel jobs | 40 (0 = no limit) | `slurm --throttle` |
| Determinism check: rebuild N molecules | 0 = off | `--repeat` |
| Job wall time (hh:mm:ss) | auto = 72 h, with the computed worst case of the biggest shard (a warning when that exceeds 72 h) | `slurm --time` |

Options at their default are left out of the command, as a user typing it would leave them out.
Under the fields an estimate (molecules, shards, at most core-hours and wall time with the
parallel jobs, every molecule at its limit) follows every change.

**Prepare** checks the list first (invalid lines of the Batch SMILES field by their line number,
an empty list, an existing run folder) and then runs `prepare`; it shows a one-line summary
(molecules, shards, at most core-hours, ~wall time with the parallel jobs, output folder) before
the details. A refusal of the CLI is shown as its first line (a missing builder environment as
the install command). **Submit** runs `slurm --submit` for the main set (and the repeat set) on a
Slurm backend, with the partitions DELFIN uses for every job and the site's `DELFIN_MODULES`
loaded on the node; on the local backend it builds the shards one after the other in the
background. **Status** shows one progress line per shard set (done/total, ok/timeout/failed,
shards complete and in progress) followed by `status -v`; **Refresh** updates that line only.
**Collect** runs `collect` (plus `repeat-stats` when there is a repeat set) and says where the
archive is and how many molecules have frames; a set already collected is not collected again.

Each button runs the very same code as the command line and prints the command it ran, so a run
from the dashboard can be repeated without it, and the same settings give the same run directory
either way (`tests/test_a_construction_batch_from_the_dashboard_is_the_cli_run.py` checks this
file by file; only the creation time in `manifest.json` differs).

## 5. Per-system limit, resume, determinism

**Limit per system** = ceil(`--timeout` x speed factor), in exact decimal arithmetic. It is not the
job's wall time. A system at the limit is killed together with its whole process group and recorded
as `timeout`; the limit decides only *which* systems finish, never the bytes of a finished one.
Calibrate the factor on the pilot (cluster core speed vs. your workstation) and use the same factor
for every tool of a comparison.

**Resume:** a system with a final record is skipped. A task that hit the wall time or lost its
node is resubmitted with the same index and the *same* speed factor (`run-shard` refuses another):

```bash
$PY -m delfin cluster slurm runs/mace_list --array 3,17-19 --submit
```

**Provenance:** `prepare` records the DELFIN commit, a content hash of `delfin/**/*.py`, the
interpreters and package versions of DELFIN's and the tool's environment, the sha256 of every
worker and shard, the MANTA construction switches and all settings in `manifest.json`.
`run-shard` checks all of it before it builds and refuses (exit 2, reason in
`refused_<time>.json`) when anything differs. `collect` reports chunks built with different code,
environment or limit as problems; `n_problems` must be 0.

**Determinism:** MANTA builds with `PYTHONHASHSEED=0`, BLAS/OpenMP single-threaded and a fixed
number of build threads (`--threads`, default 7 -- keep it fixed when bytes are compared).
Identical SMILES are built once and served to the other IDs with their own ID in the header.
`prepare --repeat N` adds a second set of N systems (one per distinct SMILES, chosen
deterministically) that is built a second time:

```bash
$PY -m delfin cluster slurm runs/mace_list --set repeat --submit
$PY -m delfin cluster collect runs/mace_list --set repeat
$PY -m delfin cluster repeat-stats runs/mace_list     # identical / different / only one / neither
```

MANTA and epic-MACE reproduce byte for byte in the same environment (epic-MACE exposes no random
seed, so this is measured, not guaranteed); Architector and molSimplify can differ from run to run.

## 6. Output

```
RUN/manifest.json                      settings, provenance, shard table (sha256)
RUN/shards/<set>/shard_NNNN.txt        ID;SMILES of each shard; specs_NNNN.jsonl (external builders)
RUN/out/<set>_<run>/chunk_NNNN/        written by run-shard (archive/, logs/, status, DONE.json)
RUN/collected/<label>/
    archive_<label>/<ID>.xyz           all frames of one system
    archive_<label>/_meta/<ID>.json    external builders: status, class, frames, wall time, RSS
    build_<label>.json                 class per ID
    buildtime_<label>.json             seconds per ID
    buildmem_<label>.json              peak memory per ID, GB
    summary_<label>.json               classes, coverage, missing shards, limit, provenance, problems
    resubmit.txt                       shards still missing, as an --array value
```

`<ID>.xyz` is a multi-frame xyz in Angstrom: per frame the atom count, the comment line
`<ID> frame<k> <label>` (external builders append `tool=<tool> E=<energy or na>`), then one
`Symbol x y z` line per atom. Frames are in the tool's order (MANTA: its isomer/conformer order;
MACE: geometry by geometry, stereomer, conformers by force-field energy). Atom order is the tool's
own; compare structures by graph isomorphism, not by atom index. epic-MACE's hapto centroid dummy
atoms are removed.

Classes:

| Class | Meaning |
|---|---|
| `ok` | at least one frame |
| `timeout` | the per-system limit was reached |
| `empty` | the tool ran and returned no frame (no structure, exception, crash) |
| `fail` | the tool's own input limits refuse the system (no geometry for the CN, metal unsupported, ...) |
| `not_expressible` | the SMILES cannot be written in the tool's input language (not cut into metal + free ligands + oxidation state 0..8; for epic-MACE also a site count without a MACE geometry). Not a failure of the tool and not in its denominator |

MANTA reads the SMILES directly and knows `ok`, `empty`, `timeout` and `fail`.
`coverage_of_expressible` = ok / (n - not_expressible); `coverage_of_all` = ok / n.

How the external builders get their input: the SMILES is cut into metal and free ligands with
DELFIN's `split_complex_smiles` (the same cut the dashboard buttons and `smiles_converter` use):
each metal-donor bond is cut as written, coordination number and oxidation state follow from it,
and the tool receives metal, oxidation state, CN and the ligands with their donor atoms. For
epic-MACE the ligands are written with mapped donors and put on the central atom with its own
`ComplexFromLigands`; a group of mutually bonded donors (a hapto ligand) becomes one centroid site.
Ligand SMILES spelling can depend on the RDKit version; to reuse the exact specs of an earlier run,
pass its `specs_*.jsonl` (concatenated) with `--specs`.

## 7. Licences

epic-MACE is GPL-3.0, Architector is BSD-3-Clause, molSimplify is GPL-3.0. DELFIN neither bundles
nor imports them: they are installed separately by the user and run as external programs in their
own environments. epic-MACE is installed from the pinned upstream commit efb5778e.
