## What this changes

<!-- What the change does, and why. If it fixes something, what was wrong. -->

## How it was verified

<!--
Unit tests are necessary and not sufficient here. If this touches how ORCA,
xTB or CREST are driven, say what you ran and what came out — they do things
the code cannot be read for. Two examples from this repository: ORCA looks for
a restart file in its working directory rather than beside the input, and it
words an out-of-memory abort differently than the pattern that was waiting for
it. Both were found by running, neither by reading.

A measurement beats an assurance: "benzene Hessian, 161 s cold, 8 s resumed,
0 displacements recomputed" says more than "the restart works".
-->

- [ ] `pytest -q -m "not slow"` passes
- [ ] Ran against the real programs, where this affects how they are driven
- [ ] Added the test that would have caught this

## Does it change what is computed?

<!--
Scheduling, restarts, logging and recovery may be reorganised freely.
Functionals, basis sets, solvation, convergence criteria and CONTROL keys
change results or break saved input files, and belong in their own change with
their own justification.
-->

- [ ] No — same numbers, different plumbing
- [ ] Yes — and the justification is above
