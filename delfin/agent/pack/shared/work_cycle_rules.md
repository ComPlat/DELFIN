# Work Cycle Rules

These rules apply to every DELFIN work cycle.

## Core rules

1. DELFIN works in closed implementation cycles.
2. Prefer the smallest useful agent set.
3. Only one agent may modify production code unless explicitly overruled.
4. Test is a hard gate unless explicitly waived by the requester.
5. Scope must stay bounded.
6. Runtime, scheduler, and recovery changes are high-risk by default.
7. The user goal and success metric are locked once the plan is issued.
8. Downstream agents may challenge the plan, but must not silently redefine it.

## Execution rules

- Split oversized requests into current cycle and next cycle.
- Define acceptance criteria before implementation.
- Break complex work into small stage gates with exit evidence.
- Do not let builder work begin before scope and criteria are clear.
- If the current gate fails, stop and surface the blocker.

## DELFIN priority order

1. correctness
2. reproducibility
3. runtime safety
4. maintainability
5. user-facing convenience
