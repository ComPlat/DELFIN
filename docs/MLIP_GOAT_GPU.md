# MLIP GOAT on a SLURM GPU cluster (orca-external-tools, orca-mlips)

## In short: what DELFIN does when CONTROL says `global_opt_method=mlip:<backend>`

Everything else is decided by DELFIN, in one place
(`delfin/mlp_tools/mlip_goat.py`, `plan_goat`), from measurements, and
written to the job log (`GOAT provider: ...`, `GOAT GPU before/after the
run: ...`) and to `XTB_GOAT.method.json`.  Every choice has a CONTROL key.

| Decision | Default | From | Override |
|---|---|---|---|
| Provider | orca-external-tools (OET) when installed with the backend, else orca-mlips | OET keeps one model server per job; orca-mlips started one model per GOAT worker | `global_opt_mlip_provider` |
| OET server | one `oet_server <backend> -b 127.0.0.1:<free port> -n W`, started before ORCA, stopped after it | OET README "Server"; ORCA 6.1 tutorial "ORCA with External Methods" | -- |
| ProgExt | `<oet env>/bin/delfin-oet-client`: runs the documented `oet_client -b <address>`; the address is read from the run directory, so the input does not change between runs (smart recalc) | -- | `mlip_fast_client` |
| Workers W (OET) | `%pal nprocs W` = GOAT workers = server model processes; W = PAL, at most what the GPU (2 GiB per model planned, 0.8 GiB measured) and the job memory (3 GiB per model) hold, a multiple of 4 | oet_server refuses a request for more cores (NCores) than its `-n`; the initial optimisation asks for all of `%pal` ("Requested 32 cores but only 8 total available", measured) | `goat_max_workers` |
| Workers (orca-mlips) | `%goat nworkers 8`, `%pal nprocs 8` (fewer on a small GPU) | 40 workers ran the GPU out of memory (orca-mlips restarted its server 259 times), 8 ran | `goat_max_workers` |
| GPU partition, `--gres` | the configured GPU partitions that hold the time limit (one scontrol query), asked with `sbatch --test-only` | D1 | `runtime.slurm.gpu_partitions`, `DELFIN_MLP_CLUSTER_GPU_PARTITIONS` |
| Speed switches | off: slim client (curl) and ORCA staged to local disk measured no faster | ethanol, see "Speed" | `mlip_fast_client`; `DELFIN_STAGE_ORCA` |

The ORCA input is the one the ORCA 6.1 tutorial describes for an MLIP GOAT
(`! GOAT ExtOpt`, `%method ProgExt "<client>"`), with the model and device
per request in `Ext_Params` (`-d cuda`; `-m <model>` when
`global_opt_mlip_model` is set).  The GOAT settings are ORCA's own; DELFIN
changes no method, convergence or GOAT keyword to gain speed.

orca-external-tools (FACCTs, github.com/faccts/orca-external-tools) is
GPL-3.0 (with a separate commercial licence from FACCTs).  DELFIN installs it
as a separate program in an environment of its own
(`python -m delfin.installer --install orca-external-tools`, pinned release
v2.1.0 held to its SHA-256, the aimnet2 extra) and only starts its commands;
none of its code is shipped with DELFIN.  Its `mace` (non-commercial
weights) and `uma` (gated weights; its dependencies conflict with aimnet2)
extras are not installed.

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
   License, non-commercial) are opt-in: installed only when named, each
   stating its terms first.

   ```
   python -m delfin.installer --install orca-mlips-mace   # into the same environment
   python -m delfin.installer --install orca-mlips-uma    # an environment of its own
   ```

   MACE (`mace-torch`, `e3nn==0.4.4`) goes into the orca-mlips environment
   beside its torch.  UMA gets `.mamba_env/orca_mlips_uma`: `fairchem-core`
   pins `torch~=2.13.0` and `e3nn>=0.5`, which neither that torch nor MACE's
   e3nn allows (pip: ResolutionImpossible); ExtOpt then names that
   environment's `orca-mlips-uma`.  UMA weights also need your own
   HuggingFace login (`hf auth login`, older huggingface_hub:
   `huggingface-cli login`) and the accepted model terms at
   https://huggingface.co/facebook/UMA.
3. **A GPU partition.** The CONTROL path runs GOAT on the node the DELFIN
   job runs on (`global_opt_mlip_device=cuda` refuses to start without a
   visible GPU); the ExtOpt submit path needs `runtime.slurm.gpu_partitions`
   (or `DELFIN_MLP_CLUSTER_GPU_PARTITIONS`).  GOAT has no local CPU
   fallback.

## Run one

**GOAT from the dashboard.**  The Submit tab has a *GOAT method* selector
above the CONTROL.txt editor: xTB (the `xTB_method`), AIMNet2, ORB, MACE or
UMA, each with whether it is installed here, the exact install command and
the licence or HuggingFace login it needs (read from the environments'
metadata; torch is never imported).  Choosing one writes
`global_opt_method=` (and `global_optimizer=GOAT` when GOAT is not selected
yet); editing the text moves the selector; a value it does not know is shown
as such and left in the text.  The optional *Model* field writes
`global_opt_mlip_model=`.  Only GOAT changes: every later step of the
CONTROL runs as before, but an MLIP GOAT sends the whole DELFIN job to a GPU
partition.

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
| GOAT, `mlip:orb`, chosen in the Submit tab's selector and sent with its SUBMIT JOB handler, PAL 8, 2 h | routed by DELFIN to the long CUDA partitions only; ran on an A100 | terminated normally; 4 conformers within 1.13 kcal/mol; `check_goat_e2e` ok | 1584 s |
| same, `mlip:mace` (opt-in MACE extra, academic use) | same routing; ran on an A100 | terminated normally; 4 conformers within 1.13 kcal/mol; `check_goat_e2e` ok | 1672 s |
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

- The same ORB and MACE GOATs with a 30-minute limit did not finish (the
  job's walltime sync stopped them in global iteration 1): on an A100 an
  ethanol MLIP GOAT took 26-28 min.  A first ORB submit had been routed to an
  AMD MI300 partition and stopped at once (no NVIDIA GPU); AMD-named
  partitions are now left out of an MLIP GOAT's candidates.
- Two lines with failure words appear in finished runs and are no failure:
  ORCA's "The lambda equations have not converged!" (one RFO step) and
  orca-mlips' "[Errno 98] Address already in use" (GOAT workers starting the
  model server at once); `check_goat_e2e` reports them as notes.

## Speed: what was measured, and when MLIP GOAT beats xTB GOAT

All runs through the CONTROL path, ORCA 6.1.1, AIMNet2 (default model), one
GPU, 30-minute GPU partitions unless noted; "step" is ORCA's `Time for
complete geometry iter`, "E+G" its `Time for energy+gradient` (the external
call).  A GOAT's total work depends on the worker count (GOAT gives each
worker at least ~15 optimisations for ethanol, ~34 for the 81-atom
chromophore), so total run times of runs with different W are not
comparable; steps and optimisations per minute are.

Ethanol (9 atoms):

| Setup | GPU | Workers | Step (E+G) | GOAT run time |
|---|---|---|---|---|
| orca-mlips, PAL 8 | H100 | 8 | -- | 9 min 28 s |
| orca-mlips, PAL 40, capped (`%pal nprocs 8`, `nworkers 8`) | H100 | 8 | 0.91 s | 8 min 33 s, one server start |
| orca-mlips, PAL 40, capped, `nprocs 40` (first version) | A100 | 8 | 1.46 s | not finished in 25 min |
| orca-external-tools, documented client | H100 | 8 | 0.85 s (0.22 s) | 7 min 36 s; 8 model processes held 6.4 GiB GPU |
| same, ORCA staged to node-local disk | H100 | 8 | 0.85 s | 7 min 7 s (within run-to-run spread) |
| same, slim client (curl), run directory on the network home | H100 | 8 | 0.93 s | 13 min 10 s |
| orca-external-tools, PAL 40, job memory 78 GiB | A100 | 24 | 1.24 s (0.36 s) | not finished in 25 min (3x the optimisations of 8 workers) |

An 81-atom push-pull chromophore, gas phase (an archived xTB GOAT of the same input: GFN2-xTB,
40 workers on 40 cores, 2 h 27 min, 162 conformers, 6.0 kcal/mol window,
~7000 optimisations, ~48 per minute):

| Setup | GPU | Workers | Step (E+G) | Optimisations per minute | Projected GOAT |
|---|---|---|---|---|---|
| orca-mlips, PAL 40, uncapped | H100 | 40 | -- | -- | GPU out of memory after 25 min |
| orca-mlips, 8 workers | A100 | 8 | 1.0-1.5 s | 3.8 | ~28 h |
| orca-external-tools, PAL 40 (DELFIN's default), only `global_opt_method` set | H100 | 40 | 1.49 s (0.22 s) | ~24 | ~5 h; 40 model processes held 32 GiB GPU, 73 GB job memory |

Where the time goes: the external energy+gradient is 0.22 s per step with
orca-external-tools; the rest of every step -- 0.6 s
for ethanol, 1.3 s with 40 workers on 40 cores for the chromophore -- is
ORCA's own per-step work outside the optimiser (the optimiser itself takes
1 ms).  A faster client (curl instead of a Python start-up) and ORCA on
node-local disk did not shorten the step, so both stay off.

Recommendation, from these numbers:

- **Small and mid-size molecules, many conformers, GPU at hand:** MLIP GOAT
  with orca-external-tools is practical (ethanol: 7-9 min) and gives
  ωB97M-D3-trained energies instead of GFN2-xTB.  It is not faster than xTB.
- **Large chromophores (~80 atoms) on a time-limited GPU partition:** xTB
  GOAT on 40 CPU cores is about twice as fast (2.5 h against ~5 h projected)
  and fits CPU partitions; MLIP GOAT needs a GPU partition with a multi-hour
  limit.  The MLIP route is worth it when the AIMNet2 ranking or geometry is
  what is wanted, not for speed.
- **Geometry (bond-length alternation of the push-pull bridge):** a
  CAM-B3LYP-D3(BJ)/def2-SVP optimisation of the lowest conformer gives
  0.051 A; GFN2-xTB 0.023 A, AIMNet2 0.087 A.  Neither pre-screen gets it;
  re-optimise the selected conformers at the DFT level before computing
  NLO properties.  Heavy-atom RMSD to the DFT minimum: xTB 0.16 A, AIMNet2
  0.43 A; the bridge torsions of all three agree within ~12 degrees.
- **Ranking:** on this chromophore the AIMNet2 ranking of the ten lowest
  xTB conformers kept the same top three as xTB and as CAM-B3LYP-D3(BJ);
  without dispersion CAM-B3LYP prefers a different conformer family by
  5.5 kcal/mol.

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
- **The ORCA 6.1 tutorial's GOAT example** (`! EXTOPT GOAT PAL8` with
  `oet_server aimnet2 --nthreads 4`; faccts.de/docs/orca/6.1/tutorials/
  workflows/extopt.html) is refused by orca-external-tools 2.1.0: its server
  counts each request's NCores against `--nthreads`, and the initial
  optimisation asks for all 8 ("Requested 32 cores but only 8 total
  available" in our PAL-32 run).  DELFIN writes `%pal nprocs` equal to the
  server's `-n`.
- **Worker sharing -- answered.** orca-mlips does not share: with many
  workers each started its own model server (the GPU ran out of memory at
  40); DELFIN caps it at 8. orca-external-tools shares one server per job.
- **Open-shell / charged systems** are not exercised: ethanol is
  closed-shell neutral.  Whether ORB and AIMNet2 handle a charged,
  open-shell input (e.g. `[Cu(NH3)4]2+`) is not established.
