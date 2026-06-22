# Mode: quick

This is the recommended daily operating mode for DELFIN.

## Objective

Close normal implementation cycles cheaply and cleanly.

## Route

Session Manager -> Builder -> Test

## Why this is the default

- one planning pass
- one implementation agent
- one validation pass
- low coordination overhead
- enough discipline for most bounded tasks

## Use this for

- small bugfixes
- medium changes in isolated modules
- docs and UX maintenance
- small refactors without runtime or architecture sensitivity
- DELFIN agent-pack maintenance

## Do not use this if

- scheduler correctness is central
- local versus SLURM behavior could diverge
- recovery or submission behavior changes
- the task touches DELFIN monoliths in a risky way
- the cycle is effectively release gating

## Acceptance expectation

- implementation is complete
- tests or validation evidence are present
- no known unreviewed high-risk behavior changes remain
- the original user goal was met, not replaced by an easier proxy
