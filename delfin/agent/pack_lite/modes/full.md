# Mode: full

Use this for high-confidence closure, not for normal daily tasks.

## Objective

Decide whether a DELFIN change set is good enough to ship, merge, or close as a
milestone.

## Route

Chief -> Session Manager -> Runtime -> Critic -> Builder -> Test

## Use this for

- release hardening
- final milestone validation
- broad user-facing change batches
- cross-cutting architecture plus runtime changes

## Why Chief is included

Release closure is a strategic decision, not just a local patch decision.

## Why Runtime and Critic come before Builder

At release time, all review agents (Runtime for operational realism, Critic for
architecture) must inform the Builder before code is written. This avoids
expensive rework discovered only after implementation.

## Acceptance expectation

- strategic scope is clear
- architecture concerns are addressed
- runtime and cluster concerns are explicitly reviewed
- implementation is complete
- test evidence is strong
- the release decision is based on the locked target, not on a substituted proxy
