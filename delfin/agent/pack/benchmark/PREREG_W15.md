# PREREG_W15: pre-registered benchmark rules for package L5 — does DELFIN learn?

Written before any measurement. The learning metric, the threshold and the
task pairs below are fixed here and only here; a later change to any of them
is a deviation and is reported as such, never silently re-scored.

## What the metric measures

Wave 15a built a learning loop (skill use counted with its outcome, a
correction proposes a skill patch, a user model, cross-session search), but
nothing measures whether a **second run of the same task** is better than the
first. `learning_metric` is that instrument: a pure function from two run
records of the same task to a verdict.

A run record is one JSONL row as written by `write_run`
(`delfin/agent/benchmark.py:2099`) — `asdict(BenchmarkResult)`
(`benchmark.py:410`). The record has NO "rounds" and NO "questions" counter;
the honest per-run efficiency proxies are:

| field | proxy for |
|---|---|
| `tool_calls` (int) | effort — fewer is better |
| `input_tokens` + `output_tokens` (int) | cost of thought — fewer is better |
| `duration_s` (float) | wall clock — fewer is better |
| `cost_usd` (float) | real spend — fewer is better |
| `success` (bool), `quality_0_100` (int), `verify_ok` (bool~None) | correctness — must not regress |

The metric does **not** invent a missing field. `rounds` and `questions` are
not in the record, so the metric never reads them.

## The function

```
delfin.agent.pack.benchmark.learning_metric.learning_verdict(first: dict,
                                                             second: dict) -> dict
```

Both inputs are run records of the SAME task (`task_id` must match). The
CALLER decides which record is first and which is second: two rows from one
run file are replicates of one snapshot and are NOT a learning comparison;
first-vs-second requires two runs separated in time (two run files). The
metric is time-agnostic and pure — it cannot tell the caller's intent apart,
so this pairing rule lives in the caller (the runner hook, phase 3).

Return value: `{"verdict": <str>, "evidence": {...}, "error": <str|None>}`
with `verdict` one of:

- **`error`** — the input could not be compared. Mirrors the guards below;
  `error` carries the reason.
- **`wrong`** — the second run did NOT achieve correctness while the first
  did, OR the second run was incorrect in a case that must not become a
  pass. **Correctness dominates: a second run that is faster but wrong is a
  fail whatever it saved.**
- **`learned`** — second is correct, correctness did not regress, and at
  least one efficiency proxy strictly improved with none worsening.
- **`no_learning`** — second is correct, correctness did not regress, and
  every efficiency proxy is equal to the first (a tie is NOT learning).
- **`regressed`** — second is correct, correctness did not regress, but at
  least one efficiency proxy strictly worsened with none improving (slower /
  more calls / more tokens than the first).

## Input guards (fixed, enforced by the metric)

1. **Same task.** `first["task_id"]` must equal `second["task_id"]`; else
   `error`.
2. **Correctness present and non-null.** `success` (bool) and
   `quality_0_100` (int) must be present on BOTH records and respectively
   bool and int; a missing field is `error`, never a silent zero.
   `verify_ok` is read only when present (`None` is legitimate when a task
   has no accept script).
3. **No unmeasured run counts as a measurement.** A record with
   `unmeasured` True (the request never reached the model) is `error` on
   either side — booking a non-run as a slower run would poison the verdict.
4. **Efficiency fields present and numeric.** `tool_calls`, `input_tokens`,
   `output_tokens` must be present and int on both records; `duration_s`
   float. A missing or unknown field is `error`, never a silent zero.
5. **Unknown field → error.** The metric validates every field it reads; a
   record it cannot read cannot be compared.

## The threshold and how "better" is decided

The **exact** economy set, no weighting: `{tool_calls, tokens_total,
duration_s, cost_usd}` where `tokens_total = input_tokens + output_tokens`
(the record carries them separately, `BenchmarkResult:432-433`). Each proxy
is compared independently against the first run — there is no weighted sum
and no "trust my judgement" term; the metric reads exactly these four and
nothing else. Because each proxy is compared on its own, there is no
tie-break to choose: a proxy is `better` (strictly less), `worse` (strictly
more) or `equal`.

Second run **learned** iff: `second.success is True`, correctness did not
regress, and every one of the four proxies is `<` or `==` the first's with
**at least one strictly `<`**. Equality on all four is `no_learning`.
Any strictly-worse proxy with none strictly-better is `regressed`. Mixed (at
least one better AND at least one worse) is `regressed` — a
partially-faster-but-partially-slower run is not clean learning.

Correctness is decided before any economy proxy is read:

- `second.success is False` → **`wrong`** regardless of efficiency.
- `second.verify_ok is False` → **`wrong`** regardless of efficiency and
  regardless of `success` (the artifact failed verification).
- `verify_ok is None` (no verify script ran) is NOT a correctness signal —
  it is never scored as a pass and never as a fail. A run with
  `verify_ok=None` contributes "no regression proven", not "correct".
- `quality_0_100` alone does not flip the verdict unless it crosses the
  `success` gate.

## Task pairs (four, fixed here)

Each pair is the SAME task run twice across two benchmark runs (`bench run
--task <id> --repeats N` in run A, again in run B, same model, N=5 each,
separated in time so a skill/memory can persist). The setup plants material a
learned skill or remembered fact makes cheaper to reach the second time; the
first run must discover it by tool calls. These are real-run scenarios for
the operator — the metric itself is unit-tested against constructed records.

| pair (task id) | what a skill/memory makes cheaper the 2nd time | expected effect |
|---|---|---|
| `l5_read_the_value_from_file` | remembering the exact file path + the unit, so run 2 skips the read-and-guess calls | fewer `tool_calls`, fewer tokens, same correct value |
| `l5_report_absent_quantity_as_absent` | a skill stating an unavailable figure is reported "absent", not guessed as zero | run 2 stays correct (`success`) with fewer exploratory calls |
| `l5_refusal_is_followed_by_allowed_form` | memory of the allowed route after a refusal, so run 2 does not retry the denied shape | fewer `tool_calls`, no denied attempt |
| `l5_name_method_with_basis_pairing` | a skill pairing a functional with its basis, so run 2 names one valid pair instead of several wrong ones | fewer calls, correct, no forbidden-signal hit |

Must not worsen (smoke): `office_refusal_is_accepted`, `fund_g_failing_test_stays_fixed_not_deleted`
— a second run of a task that has nothing to learn must not get worse.

## Report (per operator build of the runner hook)

One table: task pair, run-A verdict, run-B verdict, `learning_verdict` per
pair, plus the three values (first/second) per proxy. Numbers quoted from the
run files with their paths. A pair whose second run is `wrong` or `regressed`
is a negative result and is reported as such — the metric exists to expose
that, not to pass.
