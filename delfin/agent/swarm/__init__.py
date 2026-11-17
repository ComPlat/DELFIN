"""swarm: user-built agent swarms (package SW1).

Every member is a full DELFIN session with its own worktree, write scope
and permission gate; no swarm member is ever more permissive than the
user's own session.

Design notes (ideas cited, not code):
- CrewAI: roles + tasks declared up front (here: swarm.yaml, spec.py).
- LangGraph: explicit orchestration with human-in-the-loop interrupts
  (here: everything state-changing goes through the user; start is a
  dry-run-able plan, merge only prints the command).
- OpenAI Agents SDK: handoffs and guardrails (here: the existing
  permission gate, write_allow_globs; triage auto-answers only a
  whitelist).
- AutoGen: group chat (here: session_messages between members).
- Claude Code sub-agents / agent teams: reviewer/QS roles over a
  builder.

Library-first contract (SW3 lead, SW2 dashboard): every function takes
plain JSON-able inputs and returns JSON-able outputs. Read-only:
validate, plan, status, triage_list. Side-effecting: start, pause,
resume, triage_answer, prcheck. Events (phase change, approval waiting,
member stalled/ended) are emitted through ``status.emit_event``.
"""

from delfin.agent.swarm import spec, status  # noqa: F401

__all__ = ["spec", "status"]
