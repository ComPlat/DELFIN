# PACKAGE g1: the land gate refuses a branch without its proof, and fixes the signers — SHOULD LAND

One line: land arm a1's proof-and-signers gate on a1's branch, with three
measured repairs on top (probe reach, squash body, the wave's own table
format). The alternative arm a2 was measured on the same instrument and is
not taken: it refuses a **new** test file — every member's normal case — and
accepts a member-written probe log with no stage record.

Branch: `wave/w54/g1-int`, commits `c96cd3d38` (a1, cherry-picked),
`8510b9035`, `41dc4c5c6`, `384eee5ab`, and `09944bb26` (current
`origin/main` merged in) on the base `fe59cd49b`.

## Control, on the untouched tree (before any edit)

Worktree `scratch/ctl` at `origin/main` = `fe59cd49b`, module provenance
printed on every run (`delfin.__file__` inside that worktree).

| suite | on `origin/main` |
|---|---|
| tests/test_a_wave_branch_lands_through_the_review_gate.py | 48 passed |
| tests/test_a_wave_is_run_through_delfin.py | 42 passed |
| tests/test_delfins_own_git_reads_go_through_safe_git.py | 8 passed |
| **the three together** | **98 passed** |
| a1's new test file, on base | **13 failed, 2 passed** (the 2 are its controls) |
| a2's new test file, on base | **15 failed, 2 passed** (the 2 are its controls) |

## Decision table

| item | decision | commit | reason |
|---|---|---|---|
| a1: PROOF — the gate RUNS the decision-table test on a `git archive` of the base and in the worktree; reads the stage's real-run rows from `wave.events` (outside the cage) | take | `c96cd3d38` | red on base 13/2, green on a1 15/0; the only arm that reads the run record the member cannot write |
| a1: SIGNERS — committer checked like the author; Co-Authored-By must be the agent account at the agent's NUMBER | take | `c96cd3d38` | refuses account `999999+` for an account `424242+`; a2 does the same and I kept a1's (one mechanism, not two) |
| a1: SQUASH MESSAGE composed before the review gate, trailers dropped, one agent line | take | `c96cd3d38` | refuses a model name in the subject; site name refused with the host pinned |
| a1 fix: PROBE REACH — the stage row says the run happened, the log must show a `tool_use` on a file the branch changed | combine (repaired) | `8510b9035` | measured: a1 ACCEPTED "stage row present, log a chat" (case B). Neither half is trusted alone: the row is outside the cage, the log inside |
| a1 fix: the squash message is the PR **body** too | combine (repaired) | `8510b9035` | measured: a1 ACCEPTED a model name in the body (case D) — the body is most of what a squash publishes |
| a1 fix: the decision table word `item` (rule 17's own word) accepted like `package` | combine (repaired) | `41dc4c5c6` | measured: the gate REFUSED the table format `DELIVERY_RULES` 17 tells every member to write, calling it "no proof table" |
| compose: rule 19/20 say what the gate now does (reach, body, a NEW test counts) | combine | `8510b9035` | a rule the gate enforces but the brief does not state is a refusal members cannot predict |
| a2: whole arm | reject | — | measured below: refuses a new test file, and fails open on a member-written log |
| wave 48 h1's `LANDING_RULES` (A–F, lander-side) | not taken | — | a different object: h1 carries landing lessons into the lander's brief. This package's rule belongs in `DELIVERY_RULES` (member-side); rules 19/20 sit next to 17, which is where a member looks. No conflict, no merge needed |

## The measurement that decided between the arms

One instrument, one case set, both arms, each measured in its own tree with
the module path printed (`scratch/measure_arms.py`), cases driven through
`land()`/`check()` as the product calls them:

| case | wanted | a1 | a2 | this branch |
|---|---|---|---|---|
| A existing red-on-base test, stage row + log touching it | accept | accept | accept | accept |
| B same, but the log touched nothing | refuse | **accept** | refuse | refuse |
| C no stage row, only a member-written `scratch/realrun/*.jsonl` | refuse | refuse | **accept** | refuse |
| D model name in the PR body | refuse | **accept** | refuse | refuse |
| E a NEW test file (every member's case) | accept | accept | **refuse** | accept |
| F no decision table at all | refuse | refuse | refuse | refuse |
| **wrong verdicts** | 0 | **2/6** | **2/6** | **0/6** |

a2's "0 Co-Authored-By lines" refusal on case A in an earlier run of the same
instrument was a probe artefact: a2 requires the trailer in the body, so the
case body now carries it for every case. Stated because the first table I
read would otherwise overstate a2's defect.

a2's own gate, on its own clean tree, **passes with no `scratch/PR.md` at
all** (`scratch/selfcheck.py`: `PASSES -- 65 passed`), because
`_decision_table_refusal` returns "" when the file is missing. That is why
case C is a fail-open and not a gap in my probe.

## Prediction vs measurement (written before each measurement)

| # | quantity | predicted | measured |
|---|---|---|---|
| 1 | a1's new tests red on base / green on branch | all but the controls red, for their own rule | 13 failed + 2 passed / 15 passed |
| 2 | a2's new tests red on base / green on branch | same | 15 failed + 2 passed / 17 passed |
| 3 | the two arms disagree on some package rule | no — same spec, same red/green | **yes**: 2/6 wrong each, and on *different* cases |
| 4 | a1 repairs needed after cherry-pick | none; the arm reported its rules landed | **three** (B, D, and the house table format) |
| 5 | consumer suites unchanged by the additions | still 98 | 98 on base, 98 on a1, 98 on a2, **118** with the 20-test new file |
| 6 | the gate accepts the wave's own table format | yes — the brief and the gate are one product | **no**: refused as "no proof table". Fixed; both formats now accepted |

The gap between 4 and 6 and what the arms reported is the finding: both arms
proved their own fixtures and neither ran the format the wave ships.

## Suites run, with numbers

| suite | base | this branch |
|---|---|---|
| tests/test_a_wave_branch_lands_through_the_review_gate.py + test_a_wave_is_run_through_delfin.py + test_delfins_own_git_reads_go_through_safe_git.py | 98 passed | 98 passed |
| tests/test_the_land_gate_refuses_a_branch_without_its_proof.py (new) | 13 failed, 2 passed | **22 passed** |
| the four together | — | **120 passed** |
| + test_a_composed_wave_carries_the_landing_lessons + test_the_composer_grows_as_one_piece + test_gh1_skills_carry_wave_rules (h1's suites, after the merge) | — | **147 passed** |
| test_the_composer_grows_as_one_piece, test_the_swarm_subcommand_takes_its_own_arguments, test_gh1_skills_carry_wave_rules, test_a_dead_helper_is_not_kept_defined, test_swarm_start_job, test_swarm_start_tab, test_swarm_watch_lists_every_running_session, test_a_swarm_job_member_works_unattended, test_a_wave_directory_is_ignored_by_the_tracked_file | not measured | **109 passed** |
| test_no_secrets_are_committed, test_a_test_is_bound_to_delfin_not_to_a_machine, test_the_security_wave_holds_together, test_knowledge_private_and_swarm, test_a_registry_forgets_the_gone, test_a_dead_helper_is_not_kept_defined | not measured | **95 passed** |
| ruff `--select F821,F811,F601,F402,F823` on land.py, compose.py, the new test | — | **all checks passed** |

### The full suite: what was actually run

`tests/` is 57,071 tests (614 deselected by `-m "not slow"`) and takes about
2.8 h serially on this host; there is no `pytest-xdist` in the CI interpreter
and the background-job cap allowed one slot at a time, so an 8-way shard was
not possible. I ran **shard 0 (407 of 3,260 test files, the `test_a*` range)**
to completion on both trees:

| shard 0 of 8 | result |
|---|---|
| this branch | **56 failed, 6409 passed, 177 skipped, 6 errors** (1158 s) |
| `origin/main` control (worktree at `fe59cd49b`) | **56 failed, 6409 passed, 177 skipped, 6 errors** (1124 s) |

Identical counts, and the failure markers sit at the same progress positions
(3 %, 5 %, 7–8 %, 40–41 %, 47–48 %, 61 %, 72 %). None of the 56 is a
regression of this change. The names are environment-bound suites (dashboard
login, service sessions, file-shown-once, git-configuration) that fail on the
base for the same reasons.

**Not run: shards 1–7.** Said plainly per rule 14 rather than claimed. The
suites that consume the change (120 → 147) and the repository-wide scanners
(95) were run completely, and that is the coverage this claim rests on.

## Coordination with wave 48's h1 (the package asked for it)

While I worked, `origin/main` advanced from `fe59cd49b` to `80a43c7ff`, and
`c1e511f8c` ("A composed wave carries the landing lessons … wave 48") changed
**the same file**: `compose.py`, +120 lines, adding a `LANDING_RULES` dict
(A–F) and a `LANDER_RULES`/`.format(land_ref=…, base=…)` mechanism.

My change adds member-side rules **19/20 inside `DELIVERY_RULES`**, which is a
**string** (5,428 chars) — a different object from h1's `LANDING_RULES` dict.
I merged `origin/main` into this branch (`09944bb26`) and verified both
survive, by importing the merged module rather than reading the diff:

```
DELIVERY_RULES  type str, 5428 chars, contains "THE LAND GATE READS YOUR PROOF": True
LANDING_RULES   type dict, keys ['A','B','C','D','E','F']
```

No conflict, no merge needed at the source level: h1 briefs the **lander**, my
rules brief the **member**, and rules 19/20 sit next to 17 where a member
looks. h1's own suites (`test_a_composed_wave_carries_the_landing_lessons`,
`test_the_composer_grows_as_one_piece`, `test_gh1_skills_carry_wave_rules`) are
green on the merged tree (part of the 147).

## Prediction for the real runs (written before they ran)

| probe | prediction: base | prediction: branch |
|---|---|---|
| `proof-refused` (table names a test that also passes on `origin/main`) | the gate says **ready** — nothing on main reads the table, so nothing can refuse it | **refused**, naming the row, with a message like "passes on origin/main too, so it proves nothing" |
| `signers-refused` (committed as another address; agent account with another number) | **ready** — only `%ae` was compared, so the committer and every trailer passed | **refused**, naming the committer or the account number |
| `control-passes` (a branch meeting every rule) | ready | **ready** — the control: the new refusals are a gate, not a blanket |

The risk I accept in predicting these: the probes are model-driven and the
gate needs a `waves_root`, a `state.json` and a test interpreter from the
environment. If the model cannot build the run directory the probe measures
its own obstacle, not the gate — that would show as BOTH sides failing to
reach a verdict, which is distinguishable from the before/after above.

## The real runs: what they measured

One of the three probes produced a before/after. The other two measured the
cage, not the change — stated plainly because the package's step 6 asks
exactly that question and a refused run is a result, not a pass.

| probe | base (`origin/main`) | this branch | did it measure the change? |
|---|---|---|---|
| `proof-refused` | **no verdict**: `max_turns` after 137 calls; the sandbox refused every route to the gate (`is not on the auto-allow list` ×15) | **two verdicts, both from the gate**, on two fixture branches it built | **yes, one side only** — see below |
| `signers-refused` | **no verdict**: `done`, but the run ends `I could not confirm that by running it`; `git init/add/commit` denied (`not_allowlisted`) | **no verdict**: same wall, `I created one empty directory and stopped there` | **no** — both sides measured the sandbox |
| `control-passes` | **no verdict**: `incomplete`, 118 calls, no gate line | **no verdict**: `max_turns`; it reached the gate and the refusal came from the **pre-existing review stage**, not from my rules | **no** — both sides measured the sandbox |

**The one measured pair, quoted from `proof-refused-branch`:**

```
t1/m1: refused -- report-header: tests/test_unchanged.py passes on origin/main too (1 passed in 0.01s), so it proves nothing about this branch
t2/m1: refused -- report-header: no real run reached the change -- no real-run row for inline-xyz
```

Both strings are the branch's own code and were verified against it after the
run: `land.py:357-358` is the f-string that prints the first (`"passes on
{base} too ... so it proves nothing about this branch"`), `land.py:365` prints
the second, and `git show origin/main:delfin/agent/swarm/land.py | grep -c "so
it proves"` returns **0** — the sentences do not exist on the base, so they
cannot have been printed there. Their absence on base is not the same as an
observed base verdict, and I do not claim one.

`t2` is the stronger half of the pair: it got **past** the base-pass check
(its test is absent from `origin/main`, so it is red there by construction) and
was then refused for the real-run requirement. Two different rules, on two
fixture branches differing only in the `test` cell of one table row — so the
proof rule is not a blanket refusal.

**Prediction against measurement:**

| probe | predicted | measured |
|---|---|---|
| `proof-refused` | base ready, branch refused naming the row | branch **refused** naming the row — half right; base never reached a verdict |
| `signers-refused` | base ready, branch refused | **neither side reached the gate** |
| `control-passes` | both ready | **neither side reached ready**; the branch's refusal is the pre-existing review gate |

I predicted the obstacle risk in the paragraph above before the runs, and it is
what happened on four of the six sides: the probe's own cage, not the gate.
`signers-refused-branch` is the clearest evidence for that reading — it says
itself that it could not run the gate and stops, and its probe prompt is
unchanged between the two sides.

**What this does NOT establish.** The signers rule and the squash-message rule
have **no real-run evidence** from this wave. They are carried by tests that
are red on base and green on the branch (13 failed + 2 passed on base → 22
passed on the branch), not by a probe. Anyone reading this should treat those
two rules as test-proven only. The control probe did not reach `ready` on
either side, so the "a branch meeting every rule passes" half of the package
is proved by the test suite's control cases, not by a run.

**No fix was made on the strength of these runs**, because no run showed a
defect in the change: the only verdicts produced are the two the branch is
supposed to produce, and the refusals on the base side are sandbox denials. A
fix here would have been a guess dressed as a measurement.

**One defect the runs did show, in the wave's own instrument (not fixed
here).** `delfin-agent` is not on the PATH these probes run with, and
`python -m delfin.agent.cli` is refused by the script-content scan because
`cli.py` mentions `~/.delfin/credentials.json`. The product's own documented
entry point is therefore unreachable inside the cage without an allow-rule,
which is why four sides never reached the gate. That is outside this package's
files (`delfin/agent/swarm/land.py`, `delfin/agent/swarm/compose.py`) per rule
15, so it is recorded here rather than patched.

## What I did NOT do, and why

- **Did not run shards 1–7 of `tests/`**: no xdist, one background slot at a
  time, ~2.8 h serially. Shard 0 was run on both trees and matched exactly
  (56/6409/177/6 each); the consumer and scanner suites were run completely.
- **Did not make a fix on the strength of the real runs**: no run showed a
  defect in the change. Four of six sides measured the probe's own cage, and
  the only two verdicts produced are the ones the branch is built to produce.
- **Did not run the real probes myself**: the cage removes every key; they are
  requested in `scratch/realrun.json`.
- **Did not re-run the real probes after the merge**: the merge only added
  wave 48's `LANDING_RULES`/main history; `land.py` was untouched by it (the
  `fe59cd49b..origin/main` diff for `land.py` is empty). The signers and
  squash-message rules therefore still have **test-only** evidence, as the
  real-run section says.
- **Did not take a2's signer or squash code** even where it is equivalent: two
  mechanisms for one job is the defect this package exists to remove.
- **Did not change `check()`'s signature meaning for its only caller**:
  `land()` now passes `run`/`tid`; with neither, the proof check refuses
  (fail-closed) rather than skipping — a gate that passes when it cannot see
  the run is the fail-open a2 measured.
- **Did not touch** any file outside `delfin/agent/swarm/land.py`,
  `delfin/agent/swarm/compose.py`, `delfin/agent/swarm/`'s tests, and
  `scratch/`.
- No personal account names, home paths or host names in this report or in
  the commit messages; the agent address quoted in tests is the shipped
  default from `delfin/user_settings.py`.
