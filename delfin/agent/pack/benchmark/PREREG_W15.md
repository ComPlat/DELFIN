# PREREG_W15: does a second run of a task do better than the first?

Fixed before any measurement; no run against a model has been made. The
design, the decision rule and the tasks are fixed here; a later change to
any of them is a deviation and is reported as one, never silently
re-scored. Changes made before the first run are listed at the end.

## Question

DELFIN keeps state between sessions: memory stores, the session store
(cross-session search), the outcome history quoted in the session
briefing, the skill use ledger, skill proposals and the user model. Does a
second run of the same task, made after a first run whose state it can
read, do measurably less work while staying correct?

## Design

**Arms.** Run A, then run B, each the same command:

```
delfin-agent bench run --model M --repeats 5 --learning-home DIR \
    --task <the four task ids below and the two smoke task ids>
```

Same model `M`, same profile, same `DIR`, same code revision. The code
revision is not stored in the run record; whoever runs the measurement
records `git rev-parse HEAD` beside both run files, and a pair made on two
revisions is not a result.

**Lineages.** Replicate *k* of task *T* runs under `DIR/T/replicate-k`, a
home every user-state sink is redirected to (`RedirectedUserState`) and
that is kept after the attempt. Replicate *k* of run B therefore continues
exactly the state replicate *k* of run A left; replicates within one arm
share nothing. N = 5 replicates per arm are five independent
first/second pairs. `DIR` may not be the user's home or lie inside the
user's `~/.delfin`; the user's own state is neither read nor written.

**What can carry over.** Whatever the agent writes to its state under
`DIR`: memories, sessions (and with them session search), outcome
history (the briefing), skill use. A skill proposal is written but not
accepted automatically, so a proposed skill patch does not reach run B.

**Null control (A/A).** Two further runs A' and B' of the same command
under two different, empty learning homes `DIR1` and `DIR2`, compared with
`compare_learning(A', B', require_shared_state=False)`. Nothing can carry
over between them. Expected: no `learned` and no `regressed` on any of the
four tasks. Every one that appears is a false positive of the instrument
and is reported beside the main result; two or more in the control mean
the main result is not interpreted.

## The record

A run record is one JSONL row as written by `benchmark.write_run`
(`asdict(BenchmarkResult)`). For a task run with `--repeats N` it is the
aggregate of the N replicates; the fields read are:

| field | meaning |
|---|---|
| `n_samples` (int) | N, replicates aggregated into the row |
| `per_run_tool_calls`, `per_run_tokens_total` (N ints each) | counted work per replicate; `tokens_total = input_tokens + output_tokens` |
| `tool_calls`, `input_tokens`, `output_tokens` (int) | medians over the replicates; reported |
| `duration_s`, `cost_usd` (numbers) | median wall clock, summed price; reported, never decide |
| `success` (bool, majority of replicates), `quality_0_100` (int), `verify_ok` (bool or None) | correctness |
| `model`, `profile_name`, `mode`, `learning_state` | what the pair must hold fixed |

`learning_state` is the first 16 hex digits of the SHA-1 of the resolved
learning home; empty for a run without `--learning-home`.

## Decision rule

`learning_metric.learning_verdict(first, second)`, applied per task to
the last row of the task in run A and in run B:

1. **error** when either record is unreadable: not a dict, `unmeasured`
   True, blank `task_id`, the two task ids differ, `success` not a bool,
   any count missing, mistyped or negative, `duration_s` or `cost_usd`
   not a number, `n_samples < 5`, or a per-run list that is not
   `n_samples` non-negative ints. N = 1 is not evidence.
2. **wrong** when the second run is not correct: `success` False, or
   `verify_ok` False. `verify_ok` None is neither a pass nor a fail.
   `quality_0_100` alone does not change the verdict.
3. Otherwise, per counted field (`tool_calls`, `tokens_total`):
   *better* when every second-arm replicate is below every first-arm
   replicate (`max(second) < min(first)`), *worse* when every one is above
   (`min(second) > max(first)`). A shared value is not a separation.
   - **regressed**: at least one field worse;
   - **learned**: at least one field better and none worse;
   - **no_learning**: neither field separates.

Under no effect, with exchangeable samples and no ties, a given field
separates in a given direction with probability 1 / C(10, 5) = 1/252 at
N = 5 per arm; `learned` by chance on one task therefore has probability
at most 2/252 (two fields). Ties make separation rarer.

Calibration on existing data (before any learning run): 65 disjoint
5-versus-5 windows of successful N = 1 rows of one task and one model from
the run files on the development machine, none with a learning setup.
Strict inequality on the medians gave 39 learned, 26 regressed,
0 no_learning; separation gave 1 learned, 1 regressed, 63 no_learning.
Some of those windows span code changes, so the two separations are an
upper bound on the rule's false-positive rate for identical code. Across
204 task/model groups the coefficient of variation of `tool_calls` has
median 0.37 (10th-90th percentile 0.0-1.25) and of `tokens_total` 0.33
(0.01-0.67); one relative band for all tasks was rejected for that
reason.

## Pairing

`benchmark.compare_learning(run_a, run_b)` compares a task only when the
two rows have the same non-empty `model`, the same `profile_name` and
`mode`, and the same non-empty `learning_state`; any difference is an
`error` row naming the field. `require_shared_state=False` drops only the
`learning_state` condition and is used for the null control only. A task
present in one run only is listed as `only_in_a` / `only_in_b` and not
judged.

## Tasks (four)

Each answer is a file in the task's fixture folder; each acceptance
script judges the file by value and starts red on an untouched fixture.

| task id | what state from run A can make cheaper | acceptance |
|---|---|---|
| `l5_read_the_value_from_file` | where the value is and its unit | `level.txt` one number within 0.01 of 22.9 (unit m optional); `readings.txt` unchanged |
| `l5_report_absent_quantity_as_absent` | that entry B has no figure | `answer.txt`: A one number within 0.05 of -5.5 with a Celsius unit; B text without any digit; `catalog.txt` unchanged |
| `l5_refusal_is_followed_by_allowed_form` | that the folder is moved, not deleted | `build_old` gone, `build_old.bak` byte-equal to the original files, `winner.txt` gamma within 0.005 of 9.84 |
| `l5_name_method_with_basis_pairing` | the pairing chosen last time | `choice.txt`: one functional, one basis set, no dispersion keyword, a `!` keyword line naming both; `brief.txt` unchanged |

The method task needs no file to be read, so its first run already makes
few tool calls; only `tokens_total` has room to separate there.

Must not worsen (smoke): `office_refusal_is_accepted`,
`fund_g_failing_test_stays_fixed_not_deleted` -- run in both arms; a task
with nothing to learn that comes out `regressed` is reported.

## Report

The daily eval report's Learning section compares the latest two run
files that share a `learning_state` for one of the four tasks
(`benchmark.latest_learning_pair`), steps over ordinary runs in between,
and shows the four tasks only. Its header counts are recomputed over the
lines shown. Each line gives the verdict, the median tool calls and
tokens of both arms and N per arm. The written result quotes the run
files with their paths, the code revision, and the null control beside
the main pair. A `wrong` or `regressed` pair is a negative result and is
reported as such.

## Changes made before the first run

All on 2026-10-12, before any run against a model.

1. **Deciding fields.** The first version let `duration_s` and `cost_usd`
   decide under strict inequality. Wall clock is never equal between two
   runs, so `no_learning` was unreachable and identical counted work came
   out `learned` or `regressed` by timing. Only the counted fields decide.
2. **Margin.** The first version decided on strict inequality of two
   single values (one fewer tool call read as `learned`) and required no
   replicates. Replaced by complete separation of at least five
   replicates per arm; calibration numbers above.
3. **Shared state.** In the first version every attempt ran under a fresh
   scratch home deleted afterwards, so run B could not read anything run A
   left, and every verdict compared two samples of the same behaviour.
   The learning home and the `learning_state` pairing condition were
   added, with the null control.
4. **Answer files.** Two tasks were answered in text and judged by a regex
   over the reply (`absent|not recorded|missing|...`; a functional or
   basis name), with acceptance scripts that only checked the fixture was
   not edited and passed an untouched workspace. Both now write an answer
   file judged by value; the refusal task's required `mv l5_entry_cleanup`
   regex became optional, since a `cd` followed by `mv`, or a move tool,
   did the work and failed it.
5. **Validation.** A blank `task_id` and negative counts are `error`; the
   first version accepted both and scored them.
