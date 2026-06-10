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
   - **Probably fine**: completed normally, just no DELFIN_data.json yet.
4. Stop.  Do NOT submit any recalc.  Present the categorised list and ask
   the user which group(s) to submit.

Safety rules:
- Never call `/recalc auto` automatically — always wait for explicit user
  confirmation per folder or per group.
- Never call `/cancel all` — only suggest it if the user asks.
- For every folder you propose, include the suspected failure mode in
  one sentence so the user can decide.

Output format: English prose, table-like layout, max 30 folders shown
(if more, say "and N more similar cases").
