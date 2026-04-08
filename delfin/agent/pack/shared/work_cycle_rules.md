# Work Cycle Rules

These rules apply to every DELFIN work cycle.

## Core rules

1. DELFIN works in closed implementation cycles.
2. Prefer the smallest useful agent set.
3. Only one agent may modify production code unless explicitly overruled.
4. Test is a hard gate unless explicitly waived by the requester.
5. Scope must stay bounded.
6. Architecture-improving solutions are preferred when they stay practical.
7. Runtime, scheduler, and recovery changes are high-risk by default.
8. The user goal and success metric are locked once the Session Manager plan is issued.
9. Downstream agents may challenge the plan, but must not silently redefine it.

## Execution rules

- Split oversized requests into current cycle and next cycle.
- Define acceptance criteria before implementation.
- Break complex work into small stage gates with exit evidence.
- Make hidden tradeoffs explicit.
- Do not use extra agents unless they materially improve the cycle.
- Do not let research dominate implementation.
- Do not let builder work begin before scope and acceptance criteria are clear.
- Prefer proving one subgoal at a time over broad speculative changes.
- If the current gate fails, stop and surface the blocker instead of optimizing a different proxy.

## DELFIN-specific priority order

1. correctness
2. reproducibility
3. runtime safety
4. maintainability
5. user-facing convenience
6. feature breadth

## Hard safeguards

- Any scheduler, runtime, recovery, or cluster change needs explicit validation.
- Any public API change needs explicit compatibility review.
- Any change that grows `cli.py`, `pipeline.py`, `parallel_classic_manually.py`,
  or `runtime_setup.py` should justify why it does not worsen architecture.
- If a task cannot be closed safely in one cycle, stop at a good boundary and
  define the next cycle.
