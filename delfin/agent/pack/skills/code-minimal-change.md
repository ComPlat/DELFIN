---
name: code-minimal-change
description: The smallest diff that fixes the cause; extend what exists before adding a second mechanism
domains: code
---
# Change as little as possible

A fix is judged by what it leaves alone. The smallest diff is the one
that touches the cause and nothing else.

## Steps
1. Fix at the cause (from your cause analysis), not where the symptom
   shows.
2. Extend what exists: a parameter on the existing helper, a case in
   the existing branch — before a new module, a new config key, or a
   second code path that does the same job.
3. Keep the diff reviewable: no drive-by reformatting, no renaming
   next to the fix, no dependency added for one line.
4. If the right fix needs more than the smallest diff (a schema change,
   a public-interface change), say so and stop before doing it.

## Never
- Introduce a parallel mechanism instead of extending the one there.
- Ship a wide diff because it was convenient; every extra changed
  line is review cost and regression risk.
