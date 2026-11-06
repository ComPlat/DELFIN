# MLIP GOAT on a SLURM GPU cluster (orca-mlips)

How to run a GOAT conformer search driven by an MLIP through ORCA's
ExtOpt interface on a GPU SLURM cluster, using DELFIN's own submit path
and its end-to-end checker, and what is known to limit it right now.

> Measurement honesty: the exact run here is a **dry run** through
> DELFIN's own path (`plan_extopt_job(submit=False)`) -- the rendered
> input and sbatch script are shown and checked below.  The **live GPU
> run and its wall times are not yet measured**; they are blocked on the
> orca-mlips install action described in "Prerequisites" and will be
> appended to this file once the first real job lands on a GPU node.

## What this runs

| Run        | Method in the ORCA input                                   |
|------------|------------------------------------------------------------|
| GOAT       | `! GOAT ExtOpt` + `%method ProgExt "orca-mlips-<backend>"` |

The MLIP is called through ORCA's `ProgExt`/`Ext_Params` mechanism:
ORCA asks the MLIP for energies and gradients, GOAT drives the conformer
search, and the ensemble comes back in `<base>.finalensemble.xyz` with
`<base>.globalminimum.xyz` for the lowest conformer.  The one-node
`--nodes=1` line in the job script keeps the ORCA process and the
orca-mlips resident model server on the same GPU node (the MLIP side
is orca-mlips' auto-server, stopped with the job by the EXIT trap).

## Prerequisites

1. **orca-mlips installed** in the worker job environment, giving the
   console entry points `orca-mlips-<backend>` for the backend you
   choose.  Detect it with `detect_orca_mlips()` (reads distribution
   metadata only; never imports torch).
2. **A clean weight licence.** Use a backend whose weight licence is
   compatible with your use:
   - `orb`  -- Apache-2.0 (verified)
   - `aimnet2` -- MIT (verified)
   - `uma` -- FAIR Chemistry terms, **gated** on HuggingFace: needs the
     user's own HF login + acceptance of the terms before use.
   - `mace` -- Academic Software License, **non-commercial**.
   A weight whose licence forbids your use must not be used; prefer
   `orb` or `aimnet2` unless the gated/non-commercial terms are accepted.
3. **A GPU partition configured** for DELFIN: `runtime.slurm
   .gpu_partitions` (or `DELFIN_MLP_CLUSTER_GPU_PARTITIONS`).  GOAT has
   **no local CPU fallback**; without a GPU partition the job is refused.

## Run one

DELFIN's own path (`delfin.mlp_tools.extopt_submit`): build an
`ExtOptJobRequest` with `run_type="goat"`, `provider="orca_mlips"`, the
backend and `device="cuda"`, and submit through
`submit_extopt_job(..., submit=True)`.  This writes the input and the
sbatch script, then submits through DELFIN's scheduler client -- never a
hand-written sbatch.  With `submit=False` (the default) it is a dry run:
input + script are written, nothing is submitted.

A rendered GOAT input looks like:

```
! GOAT ExtOpt

%pal
    nprocs 8
end

%method
  ProgExt "orca-mlips-orb"
  Ext_Params "--backend orb --device cuda"
end

%goat
    nworkers 4
end

*xyzfile 0 1 <base>.xyz
```

The job script carries `--partition` (the GPU partition), `--gres=gpu:1`,
`--nodes=1` (same-node guarantee), `--cpus-per-task` shared among the
GOAT workers (`%goat nworkers` <= `%pal nprocs`), a sensible `--time`,
and an EXIT trap that stops the model server.  Note: the current render
does **not** emit `#SBATCH --export=ALL`; if your MLIP environment
depends on variables inherited through SLURM's export, either set them
in `runtime.mlp.env_activation` or add the export line explicitly
(reported as finding F1 in .gate/a4_defects.md).

## Verify a finished run

`delfin.mlp_tools.goat_e2e_check.check_goat_e2e(run_dir, base_name=...)`
checks, from the run directory alone:

- the ensemble `<base>.finalensemble.xyz` exists and is non-empty;
- the global minimum `<base>.globalminimum.xyz` is present;
- `<base>.out` terminated cleanly (the ORCA termination marker, no
  SCF/geometry failure line);
- the ensemble energies are finite and sorted (GOAT writes the lowest
  conformer first);
- whether the MLIP server ran on the GPU (a CUDA/GPU line in any
  `*.log`/`*.out` in the job directory -- log evidence, not assumption);
- wall time recorded from ORCA's `Total Job Time (sec)` line;
- whether a live-run `.pid` marker was left in the job directory.

`check_goat_e2e` is a **completion and well-formedness** check.  It does
not certify the ensemble as a final result: MLIP energies are a
pre-screen, and the physical verdict is a DFT re-ranking of the lowest
conformers.  The report's `energy_source="mlip"` and `final=False` say
exactly that.

## Known limits

- **Worker sharing is unverified.** Orca-mlips starts one resident
  server for the parent ORCA process; whether several parallel GOAT
  workers share one server cleanly has not been exercised on a real GPU
  job.  Cap the workers with `goat_max_workers` until the first real
  job shows it (worker outputs, GPU memory, server processes).
- **Weight licences** restrict the usable backends (see Prerequisites).
- **No local CPU fallback** for GOAT: a GPU partition is required.
- **Live wall times** are pending the first GPU run (the orca-mlips
  install action).  This page will record the measured wall time, GPU
  hours and the CPU-efficiency line of the job alongside the procedure
  once the run has happened.
