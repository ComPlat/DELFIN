---
name: swarm-orchestrate
description: When a DELFIN agent swarm pays off and when one session is better, and how to build one (swarm.yaml, dry-run start, triage, merge)
domains: agent-orchestration, code
---

# Orchestrating an agent swarm

## When a swarm pays off

Use a swarm when the work splits into packages that touch DISJOINT
files, each small enough for one session, and each needs its own
review loop. Good shapes: two independent features in different
directories; a core package plus a UI package over a stable contract;
parallel audits of unrelated modules.

Use ONE session instead when: the work shares files (merge conflicts
burn more time than parallelism saves); the task is smaller than the
swarm's coordination cost (spec + review + QS on both sides); or the
result is one number one workflow produces -- a second session adds
latency, not speed.

## How to build one

1. `delfin-agent swarm init` writes a commented `swarm.yaml`.
2. Edit it: goal, packages with DISJOINT write globs, tasks, roles
   (builder / reviewer / qs), limits (max_sessions, memory gate,
   start spacing).
3. `delfin-agent swarm plan "<goal>"` validates and lists members.
4. `delfin-agent swarm start --dry-run` shows exactly what would start.
   Confirm with `swarm start` (real worktrees + branches).
5. `delfin-agent swarm status [--json]` shows phases and health.
6. Approvals queue through `swarm triage`: only plain partner messages
   and a member's own test output are auto-answered; everything else
   waits for you.
7. `delfin-agent swarm merge <pkg>` runs the pre-PR check and PRINTS
   the merge command -- merging is your action, never the agent's.

## Safety invariants (never relaxed)

- Every member is a full DELFIN session with the same permission gate
  as your own session; no member is more permissive.
- Write scopes are disjoint and refuse protected paths (the safety
  frame and the chemistry core are not deliverables).
- Permission mode is never bypass.
- Auto-answers are logged; escalated items wait for you.
