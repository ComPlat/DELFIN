---
domains: code
---
# Recalc failed jobs
> Find calculations that need recalc and prepare a confirmation list.

The user wants to identify failed or stalled calculations so they can be
recalculated.  Follow this protocol:

1. Run `/recalc check-all` to scan every folder under `calc/`.
2. Run `/analyze status` to cross-reference convergence + completion data.
3. Group the findings into three categories with explicit folder paths:
   - **Definitely needs recalc**: SCF not converged, ORCA error, missing
     output file expected by the workflow.
   - **Possibly needs recalc**: started but no progress for an unusually
     long time, partial output truncated.
   - **Probably fine**: completed normally, just no DELFIN_Data.json yet.
4. Stop.  Do NOT submit any recalc.  Present the categorised list and ask
   the user which group(s) to submit.

4b. Route the failures to DELFIN's own healing, not to hand edits

For every "definitely needs recalc" folder, name the route the user
should take instead of re-submitting the same thing:

- A run DELFIN's recovery can heal (SCF/DIIS/TRAH/MPI/memory, plain
  geometry) → `/recalc` on that folder: the retry carries DELFIN's
  modified parameters automatically.
- `SCAN_CONSTRAINT_ERROR` / `SCAN_STEP_NOT_CONVERGED`
  (`RecoveryStrategy.reason_for_user()`): a recalc of the identical
  input fails identically — present the reason and the user's options
  instead of proposing a recalc.
- A converged optimisation with an imaginary mode → IMAG
  (`delfin --imag`, or the `imag` workflow on the calculation's
  .inp/.out), not a recalc.  The default noise floor
  (`allow_imaginary_freq` = -50 cm-1) skips shallow modes silently;
  healing a named mode needs an explicit `allow_imaginary_freq`
  SHALLOWER than that mode (e.g. 20 for -21.86 cm-1), proposed through
  the `/control` path with its reason and only applied after the user
  confirms.

Safety rules:
- Never call `/recalc auto` automatically — always wait for explicit user
  confirmation per folder or per group.
- Never call `/cancel all` — only suggest it if the user asks.
- For every folder you propose, include the suspected failure mode in
  one sentence so the user can decide.

Output format: English prose, table-like layout, max 30 folders shown
(if more, say "and N more similar cases").
