---
name: code-verify-and-report
description: Red before the change, green after it, and say what was not run
domains: code
---
# Verify, then report

A change is done when the evidence says so, not when the code looks
right.

## Steps
1. Before the change: the failing test or run (the reproduction).
   Keep its output.
2. After the change: the same test or run, now passing. Paste both
   outputs; "red before, green after" needs both.
3. Run the covering checks the project itself uses (its suite, its
   linter). Report unrelated failures, do not fix them unasked.
4. Report what was NOT run: the checks you skipped and why (no
   environment, out of scope). An honest "untested" beats a claimed
   green.

## Shell
- Quote values with shell metacharacters (`$ & > | * ^`): an unquoted
  value is expanded, globbed or redirected — never passed bare.
- A pipeline's exit status is its LAST command's. A failing first
  command hides behind `| tail`; read the failure from the output.
- A command that takes longer than a few seconds runs in the
  background; report progress and the final result, and wait on the
  condition instead of polling it.
- Never discard uncommitted work. `git checkout --`, `git restore`,
  `git clean` and `rm` throw the user's lines away: commit (a WIP
  commit), stash, or ask — never delete.

## Never
- Report a success you did not see on screen or in output.
- Let a skipped check look like a passed one.
- Fix unrelated failing tests silently — they are the user's call.
