# PACKAGE g1: the land gate refuses a branch without its proof, and fixes the signers — SHOULD LAND

One line: land arm a1's proof-and-signers gate on a1's branch, with three
measured repairs on top (probe reach, squash body, the wave's own table
format). The alternative arm a2 was measured on the same instrument and is
not taken: it refuses a **new** test file — every member's normal case — and
accepts a member-written probe log with no stage record.

Branch: `wave/w54/g1-int`, commits `c96cd3d38` (a1, cherry-picked),
`8510b9035`, `41dc4c5c6` on `origin/main` `fe59cd49b`.

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
| the four together | — | **118 passed** |
| test_the_composer_grows_as_one_piece, test_the_swarm_subcommand_takes_its_own_arguments, test_gh1_skills_carry_wave_rules, test_a_dead_helper_is_not_kept_defined, test_swarm_start_job, test_swarm_start_tab, test_swarm_watch_lists_every_running_session, test_a_swarm_job_member_works_unattended, test_a_wave_directory_is_ignored_by_the_tracked_file | not measured | **109 passed** |
| test_no_secrets_are_committed, test_a_test_is_bound_to_delfin_not_to_a_machine, test_the_security_wave_holds_together, test_knowledge_private_and_swarm, test_a_registry_forgets_the_gone | not measured | **94 passed** |

Repository-wide suites: run. **The full `tests/` suite was NOT run to
completion** — the background-job cap was full (8 jobs held by other work in
this session) and the foreground cap is 600 s. Stated plainly rather than
claimed: rule 14.

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

## What I did NOT do, and why

- **Did not run the full `tests/` suite**: background-job cap full, see above.
- **Did not run the real probes myself**: the cage removes every key; they are
  requested in `scratch/realrun.json`.
- **Did not take a2's signer or squash code** even where it is equivalent: two
  mechanisms for one job is the defect this package exists to remove.
- **Did not change `check()`'s signature meaning for its only caller**:
  `land()` now passes `run`/`tid`; with neither, the proof check refuses
  (fail-closed) rather than skipping — a gate that passes when it cannot see
  the run is the fail-open a2 measured.
- **Did not touch** any file outside `delfin/agent/swarm/land.py`,
  `delfin/agent/swarm/compose.py` and the one test file.
- No personal account names, home paths or host names in this report or in
  the commit messages; the agent address quoted in tests is the shipped
  default from `delfin/user_settings.py`.
