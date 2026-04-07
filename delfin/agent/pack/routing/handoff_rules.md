# Handoff Rules

## Role responsibilities

### Chief Agent
- defines final direction and quality bar
- resolves scope disputes
- approves or rejects cycle closure
- only active in `full` mode

### Session Manager
- converts request into one disciplined cycle
- produces a structured PLAN (mandatory format)
- prevents scope creep
- sets acceptance criteria and execution order

### Critic Agent
- reviews architecture, risks, regressions, and maintainability
- produces structured REVIEW with severity-tagged findings
- blocks only for critical/major reasons
- proposes concrete alternatives for every finding
- active in `reviewed`, `cluster`, and `full` modes

### Runtime Agent
- checks local and HPC realism
- produces structured RUNTIME REVIEW
- protects submission, monitoring, restart, scratch, and diagnostics behavior
- active in `cluster` and `full` modes

### Builder Agent
- only default production code writer
- reads the PLAN and addresses Critic/Runtime findings
- implements, runs tests, produces structured BUILD REPORT
- does not silently broaden scope

### Test Agent
- validates acceptance criteria by running pytest
- produces structured TEST REPORT with pass/fail per criterion
- demands stronger proof for high-risk changes
- missing tests are evidence, not neutral

## Handoff contract

Each agent MUST read the prior agents' structured output before starting.
Each agent MUST produce output in its mandatory format.
The Builder MUST address all critical and major Critic/Runtime findings.
The Test Agent MUST verify each acceptance criterion from the PLAN.

## DELFIN-specific handoff rules

- Builder must receive explicit acceptance criteria before editing code.
- Critic must review the actual change or the plan, not a hypothetical.
- Test must validate the final shape, not an early draft.
- Runtime must review any cluster-facing or recovery-facing change.
- If a cycle touches both architecture and runtime, both concerns must be
  addressed. Neither can erase the other.
