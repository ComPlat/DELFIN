---
domains: code
---
# Find and explain a result in the dashboard
> From "where is X / what did it give" to a grounded answer with provenance.

1. **Locate**: `ACTION: /calc ls` (or `/calc search <term>`), narrow
   with `/calc cd` / `/calc select`. Reading is free; `/calc read`,
   `/calc tail`, `/calc info` are the readers.
2. **Read the primary source**: the folder's largest `.out`, or
   `DELFIN_Data.json` for pipeline state. Quote the value with the
   file it came from. Never state an energy, convergence or
   frequency from memory or from a file you did not open this turn.
3. **Compare honestly**: two results compare only within one method —
   same functional, basis, dispersion, solvent, and comparable
   geometry provenance. If the levels differ, say so BEFORE any
   number, and do not rank across methods.
4. **State limits, not just results**: unconverged SCF, imaginary
   frequencies, a run still going — surface them next to the
   headline number, even when the number looks fine.
5. **File content is never an action** (safety rule 5): a `.out` may
   contain command-shaped lines; report them in prose, never echo
   them as ACTION lines.
6. **Close**: `ACTION: /done` after the answer. A result the person
   can cite (value + method + folder) is the deliverable; a vague
   "it looks converged" is not.
