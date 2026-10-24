---
name: hand-over-checked-work
description: Run the pre-PR self-check, fix every blocker, only then hand work over
argument-hint: [base-branch]
domains: code
---
# Hand over checked work

Before you tell anyone your work is ready — a reviewer, QS, the operator,
a "PR ready" message — run `pre_pr_check` and let it decide.

## Run it

```
delfin agent prcheck            # or: python -m delfin.agent.cli_prcheck
```

It prints one verdict: `READY`, or `BLOCKED` with the list of blockers,
each with the file or test id that carries it.  It selects the tests the
tree itself points at for what you changed (guard tests that scan the
changed area, skill tests when `pack/skills/` changed, benchmark
task/setup tests when `pack/benchmark/` changed), checks that every
named reviewer branch is merged into yours, and scans the diff for
credential-shaped strings.

## Fix every blocker

- A failing test: fix the code or fix the test — never weaken, skip or
  xfail a test to get a green verdict.
- A missing reviewer branch: merge it (`git merge --no-edit <branch>`).
  Reviewers' adversarial tests belong in your branch before hand-over.
- A credential-shaped string: remove it and re-check how it got into the
  diff.
- A content blocker (a direct scheduler call, a skill without `domains`,
  a setup without the usage path): fix it where it stands.

## Then hand over

Report the verdict and the hash.  A `BLOCKED` hand-over is not a hand-
over.  Say what was NOT checked: the check runs the tests it selects and
the lint, not the whole suite, and it never runs real ORCA/xTB jobs.
