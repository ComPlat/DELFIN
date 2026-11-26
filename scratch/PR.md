# PACKAGE p3 — a review finds real defects with evidence, and reports nothing it did not verify

## Locked goal

A review of a diff or branch — reached through `/review` and through the `code-reviewer` subagent tool — must produce **ranked findings that each carry a `file:line`, a concrete failure scenario and a verification verdict**, and must **report none it could not verify**. Success is measured by recall, precision and scenario-carrying rate on seeded defects.

## Items

| item | decision | commit | reason (one line) |
|---|---|---|---|
| Prediction recorded before measuring | done | `scratch/PR.md` (first write) | Rule 3: the gap between prediction and data is the finding. |
| Measurement harness `experiments/parity/p3-a1/measure.py` | added | `f29261cff` | Drives `cli.main` with the model stream scripted; 12 variants. |
| Baseline `before.json` | measured on a pristine `origin/main` worktree | `b0e510576` | 0/12 reached, used, correct. A control run, so the number is about main. |
| Findings contract | added to `subagents.py` beside `verify_subagent_report` | `b0e510576` | A review's findings are a sub-agent contract, exactly like the report verdict beside them. |
| `code-reviewer` preset | now returns `REVIEW_FINDINGS_SCHEMA` | `b0e510576` | It was the one shipped reviewer returning free prose. |
| `/review` | names the answer form it owes | `b0e510576` | It expanded to the same generic briefing as `/explore`. |
| Test | `tests/test_a_review_finding_carries_its_location_and_its_failure.py` | `b0e510576` | Red on base (22 failed, 1 passed), green on branch (23 passed). |

## Predicted vs. measured

| endpoint (12 variants) | predicted before | measured before | predicted after | measured after |
|---|---|---|---|---|
| REACHED | 0/12 | **0/12** | 12/12 | **12/12** |
| USED | 0/12 | **0/12** | 12/12 | **12/12** |
| CORRECT | 0/12 | **0/12** | 10/12 | **12/12** |

**Where I was wrong:** I predicted CORRECT would be 10/12 after, expecting the two adversarial citations (a line past EOF, a path outside the workspace) to need a second pass. They did not — all 12 were correct first time. The prediction was too pessimistic; the check for both is one read of the cited file plus a `resolve()` containment test, which is one pass by construction.

**A defect in my own instrument, found by the data:** the first after-run reported REACHED 0 while USED was 12. That is impossible, and it was the trace hook rather than the product: `sys.settrace` returns a *local* trace function, and returning `None` turns tracing off for that frame, so no nested call was ever seen. Two further corrections followed — matching the contract's *function names* rather than `subagents.py` (a file that exists on main, so a file match would have reported REACHED on the baseline), and running the baseline through a pristine `origin/main` worktree so `before` and `after` use the same corrected script. Every number above comes from that corrected script on both trees.

## What the numbers mean, and what they do not

**Measured here:** whether an unverifiable finding is withheld rather than reported (10 distinct ways of being unverifiable, each named in the withheld list); whether each survivor carries a `file:line` and a scenario; whether survivors are severity-ranked; that an empty review and an unreadable contract report different statuses; that a non-review preset's payload is untouched; that the path runs through `run_subagent` and reaches the parent.

**NOT measured here — recall and precision of a review against real defects.** Those are properties of the model reading the code, and this cage has no provider key. Scripting the model would make the number a measurement of my script. Probes are in `scratch/realrun.json` (prompt, workspace, expected difference, and for each the result that would falsify the claim).

## What I did NOT do

- **Did not complete the whole `tests/` suite.** My check script (`scratch/p3-a1/check.sh`) passes **255 tests**: the review test plus every preset-, slash-command- and schema-consuming suite. The whole-tree run was attempted five ways and did not finish inside the run budget; the last attempt (`scratch/p3-a1/suite.sh`, three name-ordered passes) was still in pass 1 when the budget ran out, and it was stopped rather than left running. **So `tests/` is NOT verified as a whole — say so rather than implying otherwise.** To close it: `sh scratch/p3-a1/suite.sh`, then read `scratch/p3-a1/suite-third{1,2,3}.log`.
- **Seven test modules cannot be collected on this host** (`cluster_tunnel` is not installed). Confirmed **pre-existing**: the same seven fail to collect on an untouched `origin/main` worktree (`--co` there reports the same 7 errors). Excluded by name in `suite.sh`. A host property, not a change.
- **Did not run a real model.** No key; the probes are handed to the integrator instead.
- **Did not touch** `report_verdict` (its `criteria` schema has callers I would have had to reclassify), any permission check, or any file outside the review surface.

## Reachability

`delfin-agent metrics record` on the branch: `unmounted: 11`, `within_ceiling: true` — the same ceiling as before. None of `verify_review_findings`, `attach_review_findings`, `REVIEW_FINDINGS_SCHEMA`, `REVIEW_PRESETS`, `REVIEW_SEVERITIES` is in the unmounted list. Reachability is also shown by a call, not by a name: the 12/12 USED column is the `subagent` tool's payload arriving at the parent with the normalised review on it, and the end-to-end test drives `run_subagent` itself.

## Suites run, with numbers

| command | result |
|---|---|
| `sh scratch/p3-a1/check.sh` (review test + 11 preset/slash/schema suites) | **255 passed** |
| the review test file on a pristine `origin/main` worktree | **22 failed, 1 passed** — the 1 is the control case |
| `ruff check delfin/ tests/ --select F821,F811,F601,F402,F823` | **All checks passed** |
| `pytest tests/ --co` on the `origin/main` worktree | 57529 collected, **7 errors** (missing `cluster_tunnel`) |
| `pytest tests/ --co` on the branch | the same 7 errors — identical, pre-existing |
| `pytest tests/` whole tree | **NOT COMPLETED** — see above |

The red-on-base proof, verbatim from that run: `module 'delfin.agent.subagents' has no attribute 'verify_review_findings'` / `'attach_review_findings'` / `'REVIEW_FINDINGS_SCHEMA'`, `KeyError: 'review'`, and `assert 'FINDINGS' in <the /review directive>`. The control case (`test_the_review_preset_stays_read_only_and_undelegating`) passes on both trees, so the file is not blanket-red.

## Three places I did not take the easy route

- **The pins that said `output_schema is None`.** Two shipped tests asserted every builtin preset had no schema. I did not delete the assertions; I kept the invariant and named the one exception, so a *second* preset acquiring a schema silently still fails them.
- **`REVIEW_FINDINGS_SCHEMA` was used above its definition.** The linter's `F821` caught it; the contract block was moved above `_BUILTIN_PRESETS` rather than the preset moved down.
- **A finding with an unknown severity** ("blocker") is kept and sorted last, not dropped. Dropping it would lose a real finding to punish a vocabulary; only its rank is unknown.

## For the integrator

1. `sh scratch/p3-a1/suite.sh` closes the whole-suite gap.
2. The probes in `scratch/realrun.json` answer recall and precision — the one number this package cannot produce.
3. `scratch/PR.md` and `scratch/realrun.json` are tracked with `git add -f` (`scratch/` is gitignored for untracked files only; earlier waves did the same).
