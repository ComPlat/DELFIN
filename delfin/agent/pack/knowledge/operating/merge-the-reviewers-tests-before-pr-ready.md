---
title: Merge the reviewer's tests before PR ready
topic: git-and-prs
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:tests/test_mm1_shipped_files.py
review_by: 2027-12-31
---
Situation: A change is declared ready while the reviewer's adversarial tests still live only on the reviewer's branch, so the pull request carries none of them and CI never runs them.
Lesson: Before PR ready, merge the reviewer's branch into the builder branch, so the pull request carries the reviewer's tests and CI runs them; the reviewer's findings are then enforced by the build, not by memory.
Why: A finding that lives only in a review message is re-litigated every round; a merged test runs on every later change to the same code.
