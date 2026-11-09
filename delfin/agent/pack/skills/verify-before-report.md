---
domains: research
---

# verify-before-report

Every number you put into a final answer must be checked against a real
source before the answer is sent. A number out of nowhere — or a number whose
evidence is only its own claim — is an estimate, not a result, and must be
marked as such. This applies to every answer, not only the ones that ask for
verification.

## When to use this skill

Use it whenever your answer contains a physical quantity that a calculation
or an experiment produced: an energy, a frequency, a distance, a dipole, a
charge, a HOMO/LUMO level, a time. Use it in particular when the number is
traced to an output file — the check must point at the actual line, not at
the claim that said the number.

## The rule, in one sentence

A number is *reported* only if you can point at a file, a line, a parsed
value and a unit of the run that produced it; otherwise it is dropped, or
said to be unverifiable.

## How to check yourself before the final answer

1. **List every number+unit in your draft answer.** Scan the prose for
   anything that looks like a result. Counts, years, parameter settings and
   job identifiers are not results — do not chase them.
2. **For each number, point at the source.** For a DELFIN run the source is a
   file in the calculation folder. The evidence ledger carries a `ref`; for a
   the run it is the calculation folder, and the value is read from its output
   file. Name `file:line` where the value stands.
3. **Convert and compare within a tolerance.** The claim's units and the
   file's units need not match — convert (Eh, kcal/mol, kJ/mol, eV; cm-1;
   Angstrom/bohr) and compare with a tolerance appropriate to the quantity
   (an energy to about 1 kcal/mol, a frequency to 1 cm-1, a distance to
   0.01 Angstrom). A difference of unit, of sign, or a rounding beyond the
   tolerance is a *mismatch*.
4. **Classify, then act on the verdict.** Use DELFIN's own verification
   (for a number traced to a run, `claim_locator`):
   - *verified* — the value is at a real `file:line` within tolerance:
     report it, citing the location.
   - *mismatch* — the output computes the quantity but not at the claimed
     value: do not report the number as it stands; say it does not match the
     run.
   - *no-source* — the run does not report the quantity: do not report the
     number; say it could not be verified.
   - *estimate* — a number with no run or no citable source: keep it only if
     you label it an estimate, and say why no source is available.
5. **Never trust the evidence text as the source.** The text of an evidence
   record is what was claimed; the *file* is what happened. A number that
   appears in the evidence text but not in the real output file is not
   verified.

## What you must never do

- Never report a number you have not located, or locate only in an answer or
  an evidence text.
- Never re-parse a program's output with a fresh regex when a DELFIN parser
  already reads it; reuse `delfin.energies` and the typed extractors.
- Never quietly round, change units, or flip a sign to make a number fit; a
  mismatch is a finding to report, not a value to repair.
- Never mark a number verified on the strength of a source you did not
  actually read this session.
