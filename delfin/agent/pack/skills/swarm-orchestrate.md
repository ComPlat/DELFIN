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

## Waves: many agents on one task list

The whole method, with the reasons behind each step, is in
docs/ORCHESTRATING_WAVES.md; `wave_compose` builds the task list of the
shape every wave has had (arms of different families, an integrator per
package, a lander, four reviewers of another family, a fixer).

A wave is for the orchestrator who runs many agents at once on work it
already split: a task list, one caged agent session per task in its own
worktree and branch, started as soon as what it waits for has ended.

1. Write the tasks: `[{"id", "brief", "base", "provider", "model",
   "after", "after_min", "after_wait_s", "collect", "no_realrun"}]`.
   `after` names `run/tid` (or `run/tid.final`); `collect` copies another
   member's `scratch/` into `scratch/arms/<tid>/`.
2. `delfin-agent swarm wave start RUN tasks.json --repo .` starts what is
   ready; `delfin-agent swarm wave guard install` lets the scheduler
   daemon run the tick every minute, which starts the rest, runs the
   stage (real runs outside the cage, the resume that reads them, the
   final mark), records events and guards the host.
3. `delfin-agent swarm wave watch RUN` waits for the next events (a
   member ended, went silent, crashed, became final); `review` shows a
   member's PR.md, commits and real runs; `resume --note` steers it.
4. `delfin-agent swarm wave emergency-stop "why"` stops every member and
   holds until a person runs `clear-stop`.
5. `delfin-agent swarm wave mcp` serves all of it over MCP (`wave_start`,
   `wave_watch`, `wave_review`, ...), for an orchestrator that is itself
   an agent.

Members commit and cannot push: review the branch, then land it.
