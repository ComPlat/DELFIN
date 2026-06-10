# Geometry convergence debug
> A geometry optimisation won't converge — find out why and fix it without faking convergence.

The user has an optimisation that oscillates, diverges, or hits MaxIter.

## 1. Evidence first (read-only)
1. Read the tail of the .out: last 5 geometry cycles — energies AND
   gradient norms. Oscillating energy ≠ slowly converging energy.
2. Check the trajectory (`*_trj.xyz`): is the structure physically
   reasonable, or is something dissociating / colliding?
3. Check for SCF trouble INSIDE the opt (SCF not converged warnings) —
   that is a different root cause than a bad optimisation surface.

## 2. Classify
- **Oscillation** around a minimum → step control problem.
- **Drift/dissociation** → wrong charge/multiplicity, missing dispersion,
  or the structure genuinely wants to fall apart.
- **SCF instability per step** → fix the SCF first (see below), not the
  optimizer.
- **Flat/floppy modes** (methyl rotors, loose complexes) → tighter
  criteria won't help; consider constraints or accept the plateau.

## 3. Fixes — grounded in the manual (%geom / %scf via search_docs)
Ranked, least invasive first; cite the manual section for each keyword:
1. Restart from the LAST geometry (not the first) — free.
2. Step control: reduce trust radius (`trust`), or switch coordinate
   system (`coordsys`) — verify names via `search_docs`.
3. Compute an exact initial Hessian (`calc_hess true`) for tricky cases —
   expensive, say roughly how expensive.
4. SCF-inside-opt trouble: better guess, damping/level shift — fix SCF
   stability before touching geometry settings.

## 4. What NOT to do (scientific integrity)
- Do NOT loosen convergence criteria just to get a "converged" stamp —
  a fake minimum invalidates every downstream number (frequencies,
  thermochemistry). If criteria are loosened deliberately (screening),
  the user must say so and it must be documented in the answer.

## Safety rules
- Propose the fix as a concrete input diff; resubmit only after explicit
  confirmation.
- Keywords only with a manual citation.

## Output format
**Befund (Zyklen/Gradienten zitiert)** → **Klasse** → **Fix-Ranking
(1-2-3 mit Kosten)** → **Diff + Frage zur Bestätigung**.
