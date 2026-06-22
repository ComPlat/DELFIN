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
- `LEANSCF_NOT_CONVERGED`
- `TRAH_SEGFAULT`
- `DIIS_ERROR`
- `GEOMETRY_NOT_CONVERGED`
- `MPI_CRASH`
- `FREQUENCY_FAILURE`
- `MEMORY_ERROR`
- `TRANSIENT_SYSTEM_ERROR`

Representative standard SCF path:

```text
attempt 1: SlowConv + MaxIter 300
attempt 2: VerySlowConv + KDIIS + DampFac/DampErr
attempt 3+: VerySlowConv + SOSCF + stronger damping
```

State is tracked in:

```text
.delfin_recovery_state.json
```

Retry generation preserves comments, geometry structure, inline basis-set
directives, `%basis`, `%ecp`, and additional `$new_job` sections where
possible.
