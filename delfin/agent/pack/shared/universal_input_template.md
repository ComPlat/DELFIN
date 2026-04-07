# Universal Input Template

Use this exact input shape for every active agent.

```text
Task:
Goal:
Why it matters for DELFIN:
Task class:
Scope:
Out of scope:
Constraints:
Affected files/modules:
Acceptance criteria:
Known risks:
Dependencies:
Current plan:
What this agent should deliver:
```

## DELFIN guidance

- `Task class` should use one of:
  `bugfix`, `refactor`, `new workflow`, `runtime/local`, `runtime/cluster`,
  `testing/reliability`, `docs/UX`, `agent architecture`,
  `strategic architecture`
- `Affected files/modules` should name real DELFIN files
- `Acceptance criteria` should be testable
- `Known risks` should mention scheduler, runtime, recovery, and API risks when relevant
