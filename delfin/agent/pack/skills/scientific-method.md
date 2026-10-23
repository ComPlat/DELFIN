# Scientific method — how to report a DELFIN result

Concrete rules for every answer that states a quantum-chemistry result.
Priority of the person: correctness stands above speed, cost and
autonomy. A less correct result is never a gain.

## Use the workflow as it is
- Start from DELFIN's validated workflows (CONTROL `method=`, the
  templates in `delfin/tools/templates/`). They encode what to compute.
- A deviation from a workflow default needs a stated reason in the
  answer: what was changed, from what, and why (for example "TightSCF
  instead of NormalSCF because the preliminary run had SCF
  oscillations"). A deviation without a reason is a defect.

## State the level of theory, always
- An energy or property is reported with its unit, method, basis set
  and solvent (if any). "The energy is -1234.5" is not a result;
  "-1234.567890 Eh at B3LYP/def2-SVP, solvent water" is.
- Do not mix levels of theory in one comparison. Two numbers compare
  only at the same functional, basis set, dispersion correction and
  solvation model.

## Check before you report
- Run `/check` (result_critic) on every finished output before quoting
  its numbers. It flags imaginary frequencies, spin contamination,
  SCF and geometry convergence and abnormal termination.
- An imaginary mode on a minimum, an unconverged SCF, spin
  contamination or a crashed run is REPORTED next to the number, never
  hidden. A caveat you omit is a wrong result.
- Never claim a run that did not happen. A number without a run or a
  cited output is an estimate.

## Mark estimates as estimates
- "Roughly", "about", "estimated" marks a value as an estimate. Say
  what it rests on (an analogy, a smaller basis set, a literature
  value) and how wrong it could plausibly be.
- An unmarked number that no output supports is treated as made up.

## Uncertainty honesty
- Quote only the precision the method supports. Do not write six
  decimals for a value whose method error is 1 kcal/mol.
- If you could not verify something, say so and name what would decide
  it. An honest "unverified" is better than a confident guess.

## Never loosen
- Convergence criteria are never loosened to make a run finish.
  If a run does not converge, report that and the diagnostics, not a
  relaxed re-run.
