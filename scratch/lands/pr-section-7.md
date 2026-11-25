### The two guard-protected patches, verified line by line

Both patches were re-checked against the bundle at `6d4516389` in this turn,
and neither is applied:

| patch | `git apply --check` | in the tree? | what it changes |
|---|---|---|---|
| `scratch/api_client-lab-run.patch` | **clean** | **no** | `_LAB_TOOL_NAMES` (:12347) gains `lab_run`; the `_DocToolExecutor` dispatch tuple (:18502) gains `lab_run`; `_is_honoured_dry_run` (:8977) learns `lab_run`'s argument-conditional `dry_run` |
| `scratch/eln-reads-in-the-gate-allowlist.patch` | **clean** | **no** | the read-only allow-list gains `eln_find`, `eln_read`, `eln_pending` |

Why they are patches and not commits: the guard's own predicate decides it,
not a guess —

```
$ PYTHONPATH=. python -c "from delfin.agent.kit_confirm import _is_protected_path; \
    print(_is_protected_path('delfin/agent/api_client.py'))"
PROTECTED
```

`delfin/agent/api_client.py` is the first entry of
`_DEFAULT_PATH_PROTECTED_GLOBS` (api_client.py:9860) — the file that decides
what the agent's own gate allows. Changing it is a self-modification by
definition, and its confirmation channel is `perms.confirm_callback`
(api_client.py:10579), which this harness leaves unset. Rule 13 names this
case and its route: ship the diff, ship the test, say it needs applying.
Going round it with `bash` is the one thing that rule forbids.

**Measured consequence, and it is the rule-5 shape exactly:** `lab_run` is
*advertised* to the model in the `lab-specialist` preset and *not reachable*
in the executor. An allow-list name is not a capability:

```
tests/test_a_lab_experiment_runs_through_the_kernel.py   14 passed, 5 failed
  FAILED test_the_tool_reaches_the_planner_and_the_runner
  FAILED test_the_tool_refuses_a_plan_it_judged_bad_and_writes_nothing
  FAILED test_a_hazardous_step_stops_the_run_and_does_not_claim_it_ran
  FAILED test_plan_mode_judges_a_plan_and_does_not_run_it
  FAILED test_every_tool_the_specialist_is_offered_actually_dispatches
  every one: {'error': 'Unknown tool: lab_run. Did you mean: lab_run, lab_read, lab_set?'}
```

The five are l2's own tests and they assert the capability the patch
supplies; the patch was reviewed and corrected by l2's integrator (its
report rejects the arm's original `lab-run-dispatch.patch`, which put
`lab_run` in `_PLAN_READONLY_TOOLS` — a NAME-only read-only set — and the
integrator measured plan mode calling it with `dry_run: false` and the
simulator's stirrer coming back `power: 'on'`). The replacement decides on the
argument instead, which is why it is worth applying rather than reverting.

This is the one thing a person must do that this run cannot:
`git apply scratch/api_client-lab-run.patch scratch/eln-reads-in-the-gate-allowlist.patch`
in a tree that permits a protected write, then
`pytest tests/test_a_lab_experiment_runs_through_the_kernel.py` — expected
19 passed, against the 14 passed / 5 failed measured here.
