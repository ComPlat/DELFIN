| same GOAT, submitted as a user does (`submit_calculation`, no sbatch flags added) | routed by DELFIN: `--gres=gpu:1`, the four configured 30-min GPU partitions; ran on an H100 | terminated normally; 5 conformers within 2.5 kcal/mol; `check_goat_e2e` ok | 518.1 s (job 10 min) |
# MLIP GOAT on a SLURM GPU cluster (orca-mlips)

How to run a GOAT conformer search (and an MLIP optimisation) driven by an
MLIP through ORCA's ExtOpt interface on a GPU SLURM cluster, using DELFIN's
own installer, its own submit paths and its end-to-end checker, and what was
measured when it was run for real.

> Measurement honesty: everything under "Measured" below ran for real on a
> SLURM cluster with ORCA 6.1.1 and the orca-mlips environment built by
> DELFIN's installer -- ethanol, closed-shell, neutral.  Earlier live jobs
> stopped on the blockers D1-D4 ("Known limits"); each is now fixed in
> DELFIN or shown to be someone else's limit, and the status is stated per
> item.  MLIP energies are a pre-screen, never a final result.

## What this runs

| Run                         | Method in the ORCA input                                  |
|-----------------------------|-----------------------------------------------------------|
| GOAT (CONTROL.txt)          | `! GOAT ExtOpt` + `%method ProgExt "<env>/bin/orca-mlips-<backend>"` |
| Opt / Opt+Freq (block)      | `! Opt ExtOpt` / `! Opt NumFreq ExtOpt`, same `%method`   |
| ExtOpt job (extopt_submit)  | any orca-mlips run type, one sbatch job on a GPU node     |

The MLIP is called through ORCA's `ProgExt`/`Ext_Params` mechanism: ORCA
asks the MLIP for energies and gradients, GOAT drives the conformer search,
and the ensemble comes back in `<base>.finalensemble.xyz` with
`<base>.globalminimum.xyz` for the lowest conformer.  orca-mlips starts a
resident model server on the first call (the `[mlip-server] Listening ...`
line in the ORCA output) and stops it with its parent ORCA process.

## Prerequisites

1. **orca-mlips installed by DELFIN's installer**, in an environment of its
   own (Python 3.12, a CUDA build of PyTorch, the aimnet2 and orb plugins):

   ```
   python -m delfin.installer --install orca-mlips
   ```

   It is built under the AI tools root (`~/.delfin/ai_tools/.mamba_env/
   orca_mlips`, or `DELFIN_AI_TOOLS_ROOT`) on the shared file system, so a
   GPU job sees it; `DELFIN_ORCA_MLIPS_PYTHON` points DELFIN at another
   environment.  `--update orca-mlips` and `--repair orca-mlips` work as for
   every other tool; `--status` reads its presence from the distribution
   metadata of that environment (`detect_orca_mlips()` does the same and
   never imports torch).  Every version is pinned in
   `delfin/ai_tools/install_ai_tools.sh` (`ORCA_MLIPS_*`).
   DELFIN writes the **absolute** `<env>/bin/orca-mlips-<backend>` as
   `ProgExt` (`delfin.mlp_tools.orca_mlips.resolve_progext`), so the compute
   node needs nothing on its `PATH`.
2. **A clean weight licence.** The installer brings only the two plugins
   whose weights are licence-clean for any use:
   - `aimnet2` -- MIT
   - `orb` -- Apache-2.0
   `uma` (FAIR Chemistry terms, gated) and `mace` (Academic Software
   License, non-commercial) are not installed: their console scripts exist
   in the environment, their model stacks do not.
3. **A GPU partition.** The CONTROL path runs GOAT on the node the DELFIN
   job runs on (`global_opt_mlip_device=cuda` refuses to start without a
   visible GPU); the ExtOpt submit path needs `runtime.slurm.gpu_partitions`
   (or `DELFIN_MLP_CLUSTER_GPU_PARTITIONS`).  GOAT has no local CPU
   fallback.

## Run one

**GOAT from CONTROL.txt.**  In a CONTROL.txt from `delfin --define`, set

```
global_optimizer=GOAT
global_opt_method=mlip:aimnet2
global_opt_mlip_device=cuda
```

and submit the folder as any DELFIN job (dashboard, agent, or
`delfin.api.submit_calculation`).  The submission sees the MLIP GOAT in
CONTROL.txt and sends the job with `--gres=gpu:1` to the configured GPU
partitions (`runtime.slurm.gpu_partitions`, `DELFIN_SLURM_GPU_PARTITIONS`,
`DELFIN_MLP_CLUSTER_GPU_PARTITIONS`, the site profile, else scontrol):
those whose `MaxTime` cannot hold the time limit are left out (one
scontrol query), the rest asked with `sbatch --test-only`, and the job goes
to all that fit.  PAL, maxcore and the time limit are never changed; no GPU
partition that can take it is a refusal, not a CPU job.  The GOAT runs in `<xTB_method>_GOAT/`
(`XTB_GOAT.inp`, `output_XTB_GOAT.out`, `XTB_GOAT.finalensemble.xyz`,
`XTB_GOAT.globalminimum.xyz`), and the global minimum becomes the start
geometry of the rest of the workflow.

**Opt or Opt+Freq as a building block.**  `run_step("orca_mlip_opt_opi",
elements=..., coordinates=..., run_type="opt_freq", backend="aimnet2",
device="cpu")` writes `! Opt NumFreq ExtOpt` and runs ORCA; the step fails
when frequencies were asked for and the output holds none.  The block runs
the MLIP on the CPU of the node it runs on; it refuses `device=cuda`.

**ExtOpt job on a GPU node.**  `delfin.mlp_tools.extopt_submit`: build an
`ExtOptJobRequest` (`provider="orca_mlips"`, `mlips_backend`, `device=
"cuda"`) and call `submit_extopt_job(..., submit=True)`.  It writes the
input and the sbatch script and submits through DELFIN's scheduler client;
before that it leaves out every configured GPU partition whose `MaxTime`
cannot hold `time_limit` (one `scontrol show partition -o`), and refuses
with each partition's `MaxTime` when none can.  The script carries
`--partition`, `--gres`, `--nodes=1`, `--ntasks=1`, `--cpus-per-task`,
DELFIN's OpenMPI settings (`delfin.dashboard.mpi_env`, oversubscribe on)
and ORCA by its full path.  `submit=False` (the default) is a dry run.

A rendered GOAT input (CONTROL path, PAL=8):

```
! GOAT ExtOpt

%pal
    nprocs 8
end

%method
  ProgExt "<AI tools root>/.mamba_env/orca_mlips/bin/orca-mlips-aimnet2"
  Ext_Params "--device cuda"
end

*xyzfile 0 1 XTB_GOAT_start.xyz
```

## Verify a finished run

`delfin.mlp_tools.goat_e2e_check.check_goat_e2e(run_dir, base_name=...)`
checks, from the run directory alone:

- the ensemble `<base>.finalensemble.xyz` exists and is non-empty;
- the global minimum `<base>.globalminimum.xyz` is present;
- the ORCA output terminated cleanly -- `<base>.out`, or
  `output_<base>.out` as the CONTROL path names it (`base_name="XTB_GOAT"`);
- **the job logs are free of the known failure shapes** (Open MPI "not
  enough slots", an `Error (ORCA_MAIN)/ProgExt not found` abort, an
  `Error (ORCA/EXT): external program failed` MLIP exit, "unrecognized
  arguments"): when one is present the checker names it;
- the ensemble energies are finite and sorted (GOAT writes the lowest
  conformer first);
- a CUDA/GPU line in a log of the job directory (log evidence, a warning
  when absent; the line it finds in an ORCA output is the `--device cuda`
  the MLIP was started with);
- wall time from ORCA's `TOTAL RUN TIME` line (ORCA 6), else `Total Job
  Time (sec)`;
- whether a live-run `.pid` marker was left in the job directory.

`check_goat_e2e` is a **completion and well-formedness** check.  It does
not certify the ensemble as a final result: MLIP energies are a
pre-screen, and the physical verdict is a DFT re-ranking of the lowest
conformers (`energy_source="mlip"`, `final=False`).

## Measured

Ethanol (9 atoms, neutral singlet), ORCA 6.1.1, orca-mlips 1.3.0 with
torch 2.15.0 (CUDA 13.2 build), aimnet 0.2.6, orb-models 0.5.5.

| Run | Where | Outcome | ORCA run time |
|-----|-------|---------|---------------|
| `--install orca-mlips` | CPU node | environment built, `--status`: installed | 162 s (whole installer) |
| GOAT, CONTROL.txt, `mlip:aimnet2`, `cuda`, PAL 8 | 1 x H100, 30 min dev partition | terminated normally; 6 conformers within 2.5 kcal/mol, global minimum written | 568.8 s (job 10 min incl. staging) |
| same GOAT, submitted as a user does (`submit_calculation`, no sbatch flags added) | routed by DELFIN: `--gres=gpu:1` to the four configured 30-min GPU partitions; ran on an H100 | terminated normally; 5 conformers within 2.5 kcal/mol; `check_goat_e2e` ok | 518.1 s (job 10 min) |
| `orca_mlip_opt_opi`, `opt_freq`, aimnet2, `cpu`, 4 cores | CPU node | converged; NumFreq: 21 vibrational modes, none imaginary (lowest 220 cm^-1) | 27.7 s |
| ExtOpt Opt, `orb`, `cuda`, 4 cores, via `submit_extopt_job` | 1 x H100, 30 min dev partition | converged, terminated normally (server `OrbMolEvaluator`); no "unrecognized arguments" | 34.9 s (job 63 s) |

- `check_goat_e2e(..., base_name="XTB_GOAT")` on the GOAT run: `ok=True`,
  `terminated_normally=True`, 6 conformers, energies sorted, relative
  energies 0.00 / 0.08 / 0.90 / 1.01 / 2.39 / 2.49 kcal/mol.
- GPU evidence: the model server was started with `--device cuda` on the
  job's one GPU and answered every call (`[mlip-server] Listening ...
  evaluator=AIMNet2Evaluator`, then `Server is ready`); the evaluator loads
  the model onto `cuda` and fails when that device is absent.  No GPU
  utilisation counter was recorded (the scheduler accounts `gres/gpu=1`,
  not its use).
- orb (D4 check): the same command line without `--backend` was
  accepted by `orca-mlips-orb`.
- D1, live: for the default 04:00:00 the first version of the check dropped
  one 30-minute partition and kept the second unread -- its own scontrol
  query came back throttled by the scheduler client.  The check now asks
  once for all partitions (unit-tested; re-run live in that submit: one query read
  all four MaxTimes; for 04:00:00 it refuses with each of them).
- GOAT ran 8 workers with one process each on `%pal nprocs 8` and finished.
  Whether the workers shared one model server or started one each cannot
  be read back afterwards (worker outputs are deleted unless
  `KEEPWORKERDATA`).

## Known limits

- **D1 -- walltime vs partition: fixed in DELFIN** (ExtOpt submit path).
  A 30-minute GPU partition refused the default 04:00:00.  The submit now
  drops partitions too short for the time limit and refuses when none is
  left; the time limit itself is never changed.  The CONTROL path is
  submitted through `submit_delfin.sh`, whose caller chooses `--time`.
- **D2 -- MPI slots: fixed in DELFIN.**  `--ntasks=1` is one Open MPI slot;
  the ExtOpt job script now sets DELFIN's OpenMPI environment
  (oversubscribe), as every DELFIN job runner does.  The CONTROL path ran
  8 GOAT workers this way.
- **D3 -- ProgExt executable not found: fixed in DELFIN.**  ProgExt is the
  absolute path of the installed environment's command.
- **D4 -- MLIP CLI argument rejected: fixed in DELFIN.**  DELFIN rendered
  `--backend <name>`, an option no orca-mlips plugin has (argparse:
  "unrecognized arguments: --backend", reproduced with orca-mlips 1.3.0 for
  aimnet2 and orb); the command itself is the backend.  `--precision` is
  rendered for orb only.
- **`! Opt Freq ExtOpt` runs no frequencies** in ORCA 6.1.1 (optimisation
  only, normal termination, no warning).  The block writes `NumFreq`.
- **Worker sharing** of one model server is not shown (see Measured).
- **Open-shell / charged systems** are not exercised: ethanol is
  closed-shell neutral.  Whether ORB and AIMNet2 handle a charged,
  open-shell input (e.g. `[Cu(NH3)4]2+`) is not established.
