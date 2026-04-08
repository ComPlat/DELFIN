# Mode: reviewed

Use this when DELFIN architecture or user-facing behavior could degrade if the
cycle is treated as a normal small fix.

## Objective

Add an architectural counterweight before implementation without paying the full
cost of a large multi-agent cycle.

## Route

Session Manager -> Critic -> Builder -> Test

## Use this for

- risky refactors
- orchestration-file changes
- changes in config semantics
- public API changes
- behavior changes that need explicit architectural review

## Typical DELFIN triggers

- `delfin/cli.py`
- `delfin/pipeline.py`
- `delfin/parallel_classic_manually.py`
- `delfin/config.py`
- `delfin/api.py`

## Why Critic comes before Builder

In DELFIN, risky architecture work often becomes expensive when criticism comes
too late. Early critique is cheaper than rework after the patch has spread.

## Acceptance expectation

- critical architectural objections are addressed
- implementation stays bounded
- validation covers the changed behavior
- success is evaluated against the real goal, not a weaker fallback metric
