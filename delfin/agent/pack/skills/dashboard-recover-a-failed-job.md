---
domains: code
---
# Recover a failed job in the dashboard
> Diagnose first, recalc second, never overwrite the old run.

1. **Diagnose before touching anything**: `ACTION: /analyze errors`
   and `ACTION: /analyze status` on the folder; `ACTION: /calc tail`
   the output's end for the real error line. Name the failure mode in
   one sentence (SCF divergence, walltime, basis/OOM, geometry) —
   a recalc of an undiagnosed folder repeats the failure.
2. **Check the failure is not a trap**: an OOM needs more memory or
   fewer cores, not a rerun; a walltime abort needs a time limit, not
   a rerun; an unconverged SCF may need a different solver — not the
   same settings again. `/recalc` only when the same settings can
   actually finish.
3. **Recalc, never overwrite**: `ACTION: /recalc <folder>` prepares a
   fresh run; the old folder and its outputs stay untouched. Never
   propose deleting or moving a failed folder to "clean up" — its
   output is the diagnosis.
4. **One folder at a time** (one destructive action per response).
   `/recalc auto` is a bulk op the dashboard blocks without the
   person's explicit "alle/neuberechnen" request — do not try to
   route around that.
5. **Confirm with the person**: `/recalc` is state-changing (safety
   rule 1). If the request asked for the recovery ("recalc that
   failed job"), the person's words are the yes. Otherwise diagnose,
   propose, `ACTION: /done`, and wait.
6. **Verify after**: `ACTION: /jobs` for the new run, `ACTION:
   /analyze status` on the recalc folder when it lands. Then
   `ACTION: /done`.
