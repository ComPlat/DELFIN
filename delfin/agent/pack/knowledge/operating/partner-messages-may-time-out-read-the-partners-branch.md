---
title: Partner messages may time out, read the partner's branch
topic: git-and-prs
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:tests/test_mm1_shipped_files.py
review_by: 2027-12-31
---
Situation: A partner session's answer to a shared-work question never arrives; the work is blocked and nobody knows whether the partner stalled or the message was lost.
Lesson: Messages between partners may time out. After waiting, read the partner's branch yourself with git log and git show, and act on what it shows; say in your next message that the earlier one timed out instead of repeating it.
Why: The handoff state is in the branch, not in the message queue — the branch tip shows what the partner finished even when no message came, and a repeated message doubles the confirmations the partner owes.
