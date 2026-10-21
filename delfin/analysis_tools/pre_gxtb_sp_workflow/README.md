# DELFIN Pre-Screening: `finalensemble.xyz` → g-xTB Opt → ORCA DFT SP

**Standalone, non-invasive workflow module**, suitable for later integration with DELFIN. It does **not** modify the clustering code or an existing DELFIN checkout. Uses Python 3 standard library, Bash, and Slurm.

## Chemistry and selection

1. Read **every** conformer from an ORCA/GOAT multi-XYZ `finalensemble.xyz` (no geometry clustering/deduplication).
2. Optimize all structures independently using **g-xTB** (`xtb --gxtb --opt`, gas phase).
3. Calculate **ωB97X-D4/def2-TZVP RIJCOSX def2/J CPCM(dichloromethane) TightSCF** single-point energies using ORCA 6.1.1 on each successfully g-xTB-optimized geometry. There is deliberately **no** ORCA `Opt` or `Freq` keyword.
4. Validate termination and geometry, and generate `results/single_point_energies.csv`, `per_conformer_status.csv`, and g-xTB end geometries sorted by the DFT-SP energy (`gxtb_geometries_sorted_by_dft_sp.xyz`). All original zero-based GOAT indices are retained.

**Scientific caveat:** This ranks single-point energies evaluated at **g-xTB minima**, not fully relaxed DFT energies or DFT Gibbs energies; selecting a few candidates for DFT Opt/Freq is a subsequent, separate workflow.

## On KIT UC3

The user already has g-xTB v2.0.1 under `~/software/gxtb-v2.0.1/xtb-6.7.1/bin/xtb`, plus ORCA 6.1.1 under `~/software/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`.

```bash
# Once this zip has been unpacked on the HPC:
cd ~/DELFIN_pre_gxtb_sp_workflow
python3 prescreen.py prepare \
  --ensemble /path/to/GOAT_run/molecule.finalensemble.xyz \
  --submit
```

The script prints a **new unique** `pre_gxtb_sp_*` folder next to the given ensemble. `--submit` immediately queues both arrays and an automatic collector. Omit `--submit` for preparation only. The same scripts work for any positive number of conformers sharing the same atom order (not hardcoded to 34 or 38).

Default charge/multiplicity: `0 1`; override with `--charge N --multiplicity M` if necessary. Paths to binaries can be overridden using `--gxtb PATH --orca PATH`.

## Array configuration, consistent with `cluster_specs(3).md`

| Stage | Slurm | CPU/task | RAM/task | Walltime | Parallelism |
|---|---|---:|---:|---|---|
| g-xTB geometry opt | `cpu`, account `ka`, `ntasks=1` | 4 | 4 GB | 2 h | Entire array; **no `%` limit** |
| ORCA DFT SP | `cpu`, account `ka`, `ntasks=8`, `ntasks-per-node=8`, `cpus-per-task=1` | 8 MPI ranks | 40 GB | 24 h | Entire array; **no `%` limit** |
| Postprocessor | `cpu`, account `ka` | 1 | 1 GB | 15 min | 1 |

`%pal nprocs 8` matches `--ntasks=8`; `%maxcore 4000`. ORCA runs by its full binary path without a prepended `srun`/`mpirun`. Compute jobs use node-local scratch and copy results back on exit. A read-only output-log snapshot is attempted before walltime if a job lasts that long. Each job has a lock, separate working directory, exit code and log; no existing run directory is overwritten. A frozen `jobs.<stage>.submission.txt` defines the array bounds.

The SP array has `afterok:<gxtb-array-job-id>`: it starts only if **all** g-xTB array tasks succeeded. Result collection runs `afterany:<sp-array-job-id>`, writes partial status even if an SP task fails, and returns failure if not all are valid. If the entire g-xTB array does not succeed, the dependent SP jobs remain blocked; run the manual reporter to diagnose:

```bash
python3 prescreen.py report --run-dir /absolute/path/to/pre_gxtb_sp_...
```

## Inspect

```bash
squeue -u "$USER"
sacct -j <gxtb-array-id>,<sp-array-id> -X --array --format=JobID,JobName,State,Elapsed,AllocCPUS,AllocTRES
```

The `submission.json` in the new run directory lists all three job IDs. Once everything has ended:

```bash
cd /path/to/pre_gxtb_sp_...
python3 prescreen.py report --run-dir .
head -n 12 results/single_point_energies.csv
```

`results/single_point_energies.csv` is **0-based GOAT index** and sorted by lowest SP energy; ΔE is in kcal/mol. All original GOAT frames and provenance are copied into the run directory. No cleanup of chemical results is automated.

## Test (without HPC)

```bash
python3 -m unittest discover -s tests -v
```

Tests use fake `xtb`, `orca` and `sbatch` binaries to validate staging, controls, submission dependencies, geometry preservation, error detection, both SLURM wrappers, and final result generation. These tests do **not** establish chemistry or actual MPI usability on UC3; that is confirmed by the preceding real 34/34 benchmark.
