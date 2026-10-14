# Generalist benchmark

14 tasks (`tasks_generalist.yaml`, task class `generalist`) that measure
DELFIN on everyday software work in repositories that are not DELFIN.
Each task is scored on the end state alone: the repository's tests, hidden
tests run on a scratch copy, its files and its git history. No pattern reads
the conversation.

## Tasks

| id | repo | toolchain | measures |
|---|---|---|---|
| `generalist_py_feature` | `textkit/` | Python | a specified feature plus its tests; hidden behaviour tests; whether the agent's own tests catch a wrong implementation |
| `generalist_ts_feature` | `invoice-ts/` | Node >= 22.6 | the same in TypeScript (`npm test` = `node --test`); no dependency added, test script untouched |
| `generalist_red_to_green` | `datewin/` | Python | a failing test given; fix the code, not the test; hidden dates; size of the fix |
| `generalist_refactor` | `shoplib/` | Python | rename and move across 8 files; tests green; no alias and no leftover of the old name anywhere |
| `generalist_conventions` | `notes-app/` | Python | house rules found only in AGENTS.md: test location, changelog line, lint, docstrings |
| `generalist_merge_conflict` | `configlib/` | Python, git | finish a stopped merge keeping both sides; merge commit with the right parents; no markers |
| `generalist_contribution` | `durfmt/` | Python, git | fix plus PR preparation per CONTRIBUTING.md: branch name, Conventional-Commit subjects, main untouched, PR description file not committed |
| `generalist_large_repo` | `inventory/` (153 files) | Python | locate a root cause behind two entry points; minimal edit; regression test |
| `generalist_investigate` | `ratelimit/` | Python | answer a question with file:line evidence, scored against a key; change nothing |
| `generalist_scope` | `ledger/` | Python | what must not happen: vendored code untouched, the user's uncommitted work kept byte for byte, no commit |
| `generalist_rust_feature` | `unitconv/` | cargo | a parser feature in a crate; hidden integration test; no dependency added |
| `generalist_rust_red_to_green` | `stats/` | cargo | a failing integration test; fix the code, not the test |
| `generalist_go_feature` | `coll/` | go | a generic function plus tests; hidden external test; gofmt and vet graded |
| `generalist_go_red_to_green` | `netcfg/` | go | a failing test; fix the bound, not the test |

Each task has a German `prompt` (the default) and an English `prompt_en`.

## Scoring

The acceptance script (`accept/generalist_*_verify.py`) prints one line per
check, `QUALITY: n` (share of check weight that held, 0-100) and `OK` or
`FAIL`. Required checks decide pass/fail; graded checks only lower the
quality. In the score card, `success_pts` (40) follows the exit status and
the routing share (20) is `20 * QUALITY / 100`, because these tasks declare
no signals. `verify_quality` is kept on every record and as the median over
replicates.

Not run is not failed: a host without node, cargo or go makes the setup exit
3 (`toolchain unavailable`), and the attempt is recorded as unmeasured.

## Toolchains

Tools are looked up on PATH with `$DELFIN_TOOLCHAIN_BIN` first; unset, it is
`~/micromamba/envs/delfin-toolchains/bin` when that directory exists. The
runner prepends it for the whole attempt (setup, the agent's own commands,
acceptance) and restores PATH afterwards. Everything runs offline
(`GOTOOLCHAIN=local`, `GOPROXY=off`, `CARGO_NET_OFFLINE=true`); no fixture
has dependencies.

## Running

From the DELFIN checkout (the fixture lives under
`tests/fixtures/generalist_workspace/`):

```sh
delfin-agent bench run --model kit.glm-5.3 --task 'generalist_*' --repeats 5
delfin-agent bench run --model kit.glm-5.3 --task 'generalist_*' --repeats 5 --lang en
delfin-agent bench run --model kit.deepseek-v4-flash --task generalist_scope,generalist_refactor
```

Quote the glob so the shell does not expand it. `--lang en` runs the English
prompts under the ids `<id>.en`, so German and English results never share a
baseline. `bench audit <run>` shows each failed check's line from the
acceptance output.

## A/B procedure

1. Calibrate on main first: `--repeats 5` per task on each model (GLM and
   DeepSeek). A task whose pass rate is 0/5 or 5/5 on both models cannot
   show a change in that direction; note it before comparing anything.
2. Arms: the branch under test against `origin/main` at a recorded commit,
   both in their own checkout (`git rev-parse HEAD` into the output of each
   arm). Same model, same `--lang`, history isolated, N >= 5 per arm per
   task. Alternate the arm order per batch, so a backend slow window does
   not land on one arm.
3. Decide on the pre-registered rule, per model: keep the change when the
   generalist pass count rises by >= 3 of 70 (14 tasks x 5) on both models,
   no single task loses >= 2 of 5, and median wall time per task rises by no
   more than 1.5x. Report `verify_quality` medians beside the pass counts:
   a change that moves quality without moving passes is evidence of
   direction, not of a result.
4. Unmeasured attempts (toolchain missing, endpoint down) are excluded from
   both arms and counted in the report; rerun them rather than imputing.

## Offline validation

`tests/test_generalist_tasks_score_the_end_state.py` builds every fixture
with the real setup script and judges three scripted attempts with the real
acceptance script: a correct one (must pass, quality >= 80), a plausible
wrong one (must fail, lower quality), and an untouched repository (must
fail); three tasks have a second wrong attempt. The CI job `toolchains`
runs it with node, cargo and go installed and a missing toolchain counted
as a failure.
