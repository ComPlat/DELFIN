# DELFIN ORCA Recovery and Retry Logic

DELFIN can retry failed ORCA jobs automatically when
`enable_auto_recovery=yes`.

Core behavior:

- classify the failure from the ORCA output
- build a recovery strategy for the detected error type
- write `input.retryN.inp`
- continue from the latest usable state when possible

Main settings:

```ini
enable_auto_recovery=yes
max_recovery_attempts=3
```

Common recovery actions:

- add `MOREAD` and `%moinp`
- update coordinates from the latest `.xyz`
- modify `%scf` or `%geom`
- reduce PAL after MPI or memory-related failures
- skip frequency steps after repeated frequency/LEANSCF failures

Main error classes:

- `SCF_NO_CONVERGENCE`
- `LEANSCF_NOT_CONVERGED` -- in ORCA 6 every SCF runs in orca_leanscf, so this
  is the ordinary SCF failure (all 472 archived ones were the main SCF); it
  gets the SCF escalation and keeps FREQ. Only an SCF that failed inside the
  analytic Hessian gets the frequency path
- `TRAH_SEGFAULT`
- `DIIS_ERROR`
- `GEOMETRY_NOT_CONVERGED`
- `MPI_CRASH`
- `FREQUENCY_FAILURE`
- `MEMORY_ERROR`
- `TRANSIENT_SYSTEM_ERROR`
- `INPUT_ERROR` -- ORCA refused the input (e.g. a multiplicity that cannot
  match the electron count); named, not retried
- `ESD_RATE_UNPHYSICAL`, `ESD_WINDOW_TRUNCATED` -- see below

A run that ends with "ORCA TERMINATED NORMALLY" is also retried when its
result is not one: an optimisation that ran out of cycles, a collapsed
TD-DFT root, an ESD rate that is negative or whose correlation function was
cut off before it decayed. When the retries are used up the last result is
kept, as before, and the ESD reports mark such a rate `[NOT A RESULT: ...]`.
An output that recalc kept is not judged again.

Representative standard SCF path:

```text
attempt 1: SlowConv + MaxIter 300
attempt 2: VerySlowConv + KDIIS + DampFac/DampErr
attempt 3+: VerySlowConv + CNVSOSCF true + stronger damping
```

How a retry is written (`delfin/orca_recovery.py`, through
`delfin/common/orca_input.py`, which keeps every line it does not change):

- only the job that failed is changed; in a file with `$new_job` that is the
  job whose `JOB NUMBER` banner ORCA printed last. `%tddft` settings go to a
  job that has a `%tddft` block, never into one without (in an optimisation
  that would make ORCA optimise an excited state)
- a repeated block is folded the way ORCA reads it (later value wins) before
  it is changed, so `%scf maxiter 125 end` + `%scf BrokenSym 1,1 end` become
  one `%scf` with both
- `%scf` switches are ORCA's own values: `CNVSOSCF true` (a bare `SOSCF`
  line opens a sub-block), `Convergence Tight` / `VeryTight`
- a job without `%base` gets `%base "<job>"`, so the retry's `.gbw`/`.xyz`/
  `.hess` replace the failed run's instead of landing beside them
- `MORead` uses the job's own `.gbw` (copied to `<job>_old.gbw`); without one
  the job keeps the guess DELFIN gave it, or ORCA's own if that file is gone
- an optimisation restarts from its own `<job>.xyz` if it holds the same
  atoms; per-atom `NewGTO ... end` stays on its atom
- MPI fixes also pass OpenMPI settings to the rerun
- a retry that would not read as ORCA input is not written

Measured on 543 archived inputs x every strategy and attempt (14 350
retries): the dictionary rebuild this replaced wrote 7 519 wrong ones (nested
`%scf`, lost `BrokenSym`, dropped metal basis, invalid `Convergence` values,
`%tddft` in the wrong job); ORCA 6.1.1 reads every retry written now.

State is tracked in:

```text
.delfin_recovery_state.json
```

Retry generation preserves comments, geometry structure, inline basis-set
directives, `%basis`, `%ecp`, and additional `$new_job` sections.

ESD rates (`delfin/common/esd_numerics.py`): ORCA chooses the time window of
the correlation function from the linewidth (2934 fs for LINEW 50 cm-1 on
6.1.1). The older template's `ESD_MAXTIME=12000` (290 fs) cut it off at 6.5 %
of its amplitude -- formaldehyde ISC 5.70e7 instead of 2.46e7 s-1, IC
negative -- and is now read as `auto`. A window that is still cut short is
handed back to ORCA on the retry; a negative rate is rerun with more points,
as ORCA's own warning asks. Rate jobs run with TightSCF, like the state jobs
they start from.

