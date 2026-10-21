---
title: Run the failing check on the base commit before blaming the change
topic: git-and-prs
scope: always
provenance: delfin-maintainers | 2026-10-20 | pr:337
review_by: 2027-04-18
---
Situation: A test or check fails after a change, and it is not obvious that the change touched what failed.
Lesson: Run the same check on the base commit (a worktree of main or the commit before the change) and report both results side by side. Only a check that passes on the base and fails on the change is a regression of the change.
Why: Host-dependent and pre-existing failures look exactly like regressions. A control run costs one command and separates the two; skipping it has sent reviews after failures that were already on main.
