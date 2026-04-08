# Universal Input Template

Use this exact input shape for every active agent.

```text
Task:
Goal:
Why it matters for DELFIN:
Task class:
Goal lock:
Scope:
Out of scope:
Constraints:
Affected files/modules:
Acceptance criteria:
Stage gates:
Evidence / oracle:
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
- `Stage gates` should be small enough that success or failure is obvious
- `Evidence / oracle` should state how "better" will be proven
- `Known risks` should mention scheduler, runtime, recovery, and API risks when relevant
