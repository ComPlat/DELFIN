---
title: Run the repository's guard tests before PR ready
topic: git-and-prs
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | pr:376
review_by: 2027-12-31
---
Situation: A package passes all of its own tests and still fails the repository's CI, because a guard test outside the package's own set catches what the package touched (a direct scheduler call, a path rule, a network or secret scan).
Lesson: Before declaring a change ready, run through the gate every guard test that matches what the change touches — scheduler, paths, network, secrets — in addition to the package's own tests.
Why: Guard tests exist for the crossings between packages, so a package's own suite cannot see them; on PR 376 the package's 33 tests passed while the direct-scheduler-call guard failed the PR in CI.
