# Example Cycle Kickoff

Use this when starting a new DELFIN cycle.

```text
Task:
Refactor DELFIN runtime path handling without changing public behavior.

Goal:
Reduce hidden branching in runtime resolution and improve testability.

Why it matters for DELFIN:
Runtime resolution is part of DELFIN's platform value across local and HPC setups.

Task class:
refactor

Scope:
`delfin/qm_runtime.py`, `delfin/runtime_setup.py`, related tests.

Out of scope:
Dashboard redesign, new tool integrations, packaging changes.

Constraints:
No public CLI behavior change.
No regression in local or SLURM-backed paths.
One production code writer only.

Affected files/modules:
`delfin/qm_runtime.py`
`delfin/runtime_setup.py`
`tests/test_qm_runtime.py`
`tests/test_runtime_setup.py`

Acceptance criteria:
- existing behavior preserved
- logic becomes easier to test
- tests cover the changed paths
- runtime risks explicitly reviewed

Known risks:
tool detection regressions, path normalization changes, cluster-specific assumptions

Dependencies:
existing runtime tests, DELFIN runtime conventions

Current plan:
Session Manager -> Runtime -> Builder -> Critic -> Test

What this agent should deliver:
role-specific output using the DELFIN_AGENT schemas
```
