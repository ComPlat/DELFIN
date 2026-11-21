# FL1: what is flaky on main, with the denominator

Measured 2026-11-21 over the **150 most recent `main` CI runs**
(2026-11-19 to 2026-11-21): **144 success, 5 failure, 1 with no
conclusion**. One of the five is a dependabot docker update, so four are
test failures.

"Flaky" here means *failed in N of M runs*, not *failed once*. A test that
failed once has one observation, and a quarantine without a named cause is
not a quarantine — it is a test switched off.

## Not flaky: six tests that failed in 2 of 2 slow runs

| file | tests |
|---|---|
| `tests/test_fffree_chelate_backbone.py` | `test_determinism_with_flag`, `test_large_kappa3_native_with_flag` |
| `tests/test_fffree_metalloid_donor.py` | `test_distibine_build_is_deterministic`, `test_distibine_builds_coordinated_frame` |
| `tests/test_generalist_tasks_score_the_end_state.py` | `test_a_scripted_attempt_is_judged_on_the_end_state`, `test_the_correct_solution_grades_above_the_wrong_one` |

The same six failed in the slow run of 2026-11-20 (`44730210721`) **and**
of 2026-11-21 (`44897851637`). Two of two is not a flake.

**Control on a developer host, current `main` (`009497359`): 96 passed, 2
xfailed, 0 failed** — all three files green. So the code is not broken; the
difference is the environment, and the slow workflow does install a real
xtb 6.7.1, so a missing tool is not the explanation either.

**What the log could say, and why it could not.** The failure reduces to

```
>  assert isinstance(res, list) and len(res) >= 1
E  assert (False)
E   +  where False = isinstance(None, list)
```

`converter_backend._fffree_isomers` declines a system by returning `None`,
and it **already names the reason**: `_scope_no` emits
`[FFREE_SCOPE] <REASON>` (`DECOMPOSE_NONE`, `LIGAND_SMILES_FAIL`,
`ENUM_ERROR`, …). It is silent unless `DELFIN_FFFREE_ISO_TRACE` is set, and
the slow job did not set it. So the instrument existed, was correct, and
was switched off in the one place its output was needed.

Driven here to confirm the reason is emitted and names the decline:

```
[FFREE_SCOPE] DECOMPOSE_POLYNUCLEAR
[FFREE_SCOPE] DECOMPOSE_NONE
```

**Change:** the slow job sets `DELFIN_FFFREE_ISO_TRACE=1`. Value `1` (no
path) writes to stderr, which pytest captures per test and prints in the
failure report — so the reason lands in the log of the run that failed and
nowhere else. `tests/test_fl1_a_declined_build_names_its_reason.py` pins
the flag the workflow sets to the flag the code reads, because this
instrument has been lost once already: `_ff_trace_write`'s docstring
records that it wrote to stderr while its caller discarded stderr, so it
"wrote into the void for a week".

The next slow failure of these six names its own cause. Until then the
cause is unknown, and this file does not guess it.

## One observation each, and therefore no verdict

| test | CI | local |
|---|---|---|
| `test_nothing_an_agent_command_starts_outlives_it.py::test_a_background_job_outlives_the_turn_that_started_it` | 1 failure (`44916554135`), the other CI runs of that commit green | **5/5 passed** |
| `test_the_classifier_waits_as_long_as_the_setting_says.py::test_the_same_verdict_without_the_setting_times_out_and_asks` | 1 failure (`44685614736`) | **5/5 passed** |

One failure in roughly seventy-five runs each, no named cause, and five
local passes. That is not enough to call either one flaky and not enough to
quarantine either one. The denominator is recorded here so the next
occurrence is the **second** observation rather than a fresh surprise.

A third, from 2026-11-21 outside this window and already judged the same
way: `test_a_kept_session_comes_back_in_the_full_chat` failed one CI shard,
passed on a rerun with no change, and passes locally on two checkouts —
three observations before it was called flaky.

## Not reproduced: a rewritten benchmark baseline

A claim that a suite run had silently rewritten the tracked
`tests/fixtures/benchmark_baseline.json` with failure rows is **withdrawn**.
The file has exactly one commit in its history (`14705269c`, 2026-08-15);
the tracked copy holds 7 entries and 0 failure rows; five swarm worktrees
from the preceding two days agree and none has it modified; and nothing in
the repository sets `UPDATE_BENCHMARK_BASELINE`, which is the only way
`_save_baseline` is reached.

What remains from that reading is a design weakness rather than an
incident: the update path writes whatever the run just produced, with no
check that the run was sound, so one re-baseline on a host without xtb
would record failures **as** the baseline and every later comparison would
pass. It has never fired.

## What this did not measure

Three limits, each one a boundary of this measurement rather than a
finding:

- **The cause of the six slow failures is not measured.** This work turned
  the builder's reason channel on; it did not learn why the builder
  declines in CI and passes here. The next slow failure will name its own
  cause, and until it runs that cause is unknown.
- **Whether the two single-failure tests are flaky is not verifiable from
  what exists.** One CI failure and five local passes each is one
  observation in each direction. A denominator was recorded precisely
  because a verdict could not be.
- **The environment difference itself could not be verified** from the
  runner side: the slow workflow installs a real xtb 6.7.1, so the obvious
  missing tool is ruled out, but dependency versions under
  `constraints-ci.txt` were not compared against this host's, and no
  claim is made about them.
