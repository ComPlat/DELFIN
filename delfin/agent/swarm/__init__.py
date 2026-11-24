"""swarm: user-built agent swarms (package SW1).

Every member is a full DELFIN session with its own worktree, write scope
and permission gate; no swarm member is ever more permissive than the
user's own session.

Design notes:
- roles and tasks declared up front (swarm.yaml, spec.py);
- explicit orchestration with the user in the loop: everything
  state-changing goes through the user; start is a dry-run-able plan,
  merge only prints the command;
- handoffs and guardrails: the existing permission gate,
  write_allow_globs; triage auto-answers only a whitelist;
- members talk through session_messages;
- reviewer and QS roles over a builder.

Waves (wave.py, stages.py, realrun.py, guard.py, friction.py,
mcp_server.py) are the orchestrator's unit beside the spec: a task list
whose members start as their dependencies end, are watched through an
event log, get real runs outside their cage, and are steered through the
same verbs on the command line (``delfin-agent swarm wave``) and over MCP.

Library-first contract (SW3 lead, SW2 dashboard): every function takes
plain JSON-able inputs and returns JSON-able outputs. Read-only:
validate, plan, status, triage_list. Side-effecting: start, pause,
resume, triage_answer, prcheck. Events (phase change, approval waiting,
member stalled/ended) are emitted through ``status.emit_event``.
"""

from delfin.agent.swarm import spec, status  # noqa: F401

__all__ = ["spec", "status"]
