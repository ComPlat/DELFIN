# Goal Decomposition Rules

These rules exist to stop agent drift and "solved a different problem" failures.

## Contract-first execution

1. Lock the real goal before proposing implementation.
2. State which proxy metrics are acceptable and which are misleading.
3. Break the cycle into small stage gates with explicit exit evidence.
4. Advance only when the current gate is actually satisfied.
5. If a gate cannot be satisfied, stop, report the blocker, and ask instead of drifting.

## Required shape of a good plan

A good Session Manager plan must contain:

- a locked goal statement
- explicit in-scope and out-of-scope boundaries
- acceptance criteria
- stage gates
- evidence or oracle for each important gate
- the concrete failure mode to avoid

## Anti-patterns

- redefining `works better` into `no longer crashes`
- replacing a quality goal with a throughput goal without approval
- broad refactors before the failing behavior is reproduced
- claiming progress without evidence for the current gate
- using fallback output as proof that the original goal was met

## Role expectations

- Session Manager: define the gates and lock the target.
- Research/Critic/Runtime: challenge weak proxies, missing oracles, or gate order.
- Builder: execute one gate at a time and report gate status explicitly.
- Reviewer/Test: verify the implementation against the locked goal, not a rewritten one.
