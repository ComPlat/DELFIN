# A self-driving bench: a campaign, a closed loop, its ELN hand-over, and the loop's plan judged before anything moves

Five lab packages from wave 44, merged into one bundle on
`wave/w44/zz-land` (base `origin/main` = `41fbf18e2`). One row per package
below, then what was predicted against what was measured, then the suites,
then the real runs.

This bundle is the point at which the packages stop being separate branches
and start sharing a tree — and that is where the defects in §3 came from.
They were not visible on any single branch.

`scratch/PR.md` (this file) came in with the e1 merge, which used it for its
own report; every package's report is preserved unmodified under
`scratch/arms/<pkg>-int/PR.md`, with its real runs in
`scratch/arms/<pkg>-int/realrun/`.

## 1. Per package

| package | decision | commits in the bundle | measured effect |
|---|---|---|---|
| **l1** — ChemASAP-shaped parallel campaign | **landed** | `94c88885c`, `4a0c26b4a`, `6620d2e65`, `0ed3acc17` (merge `b2a2fb59c`) | `transfer` moved a sample in the registry and never moved the container, so the campaign could not complete: `delfin-agent lab run` reached 55 of 65 steps and stopped at `measure-2` with "analysis reads container '' (fresh=True), the registry expects 'vial-2'". It now completes. 8 of the campaign's 12 hand-overs had no carrier to move into (R1–R4, P1–P4), which the first fix missed and `0ed3acc17` covers. |
| **l4** — closed loop on the simulator | **landed** | `6808cd982`, `02d25b3d1`, `5d1472c9b`, `f654f136b`, `a99f6ffe1` (merge `57b7fd5e9`) | A closed loop reached through `cli.main`, with the watchdog stop neither arm delivered. Its measured claim: EI on a GP reaches 99 % of the surface best in 4 rounds against the lattice's 31. |
| **l2** — the lab specialist | **landed, with a guard-protected half outstanding** | `4a38e2926`, `da5694765`, `e158e62f3`, `7dccca200`, `222edc62a` (merge `ef882588d`) | The judge (`delfin/lab/plan.py`) and the specialist exist. Its own report recorded the decisive defect: the branch advertised `lab_run` and dispatch answered `Unknown tool: lab_run` on three consecutive calls. That fix is in `api_client.py`, so it ships as `scratch/api_client-lab-run.patch` — **5 tests in `test_a_lab_experiment_runs_through_the_kernel.py` are red until it is applied**, measured below. |
| **l3** — run record to Chemotion | **landed** | `16fea816a`, `3c588f298`, `4411f2563`, `662efe664`, `819a3cedf` (merge `a2f16f384`) | A finished run leaves one file a person imports: research plan, reaction with conditions, each sample's analyses under the ELN's own container tree, provenance. The tree's file names now match bundle members for a repeat measurement, and the upload request carries the shape the import built. |
| **e1** — Chemotion ELN specialist | **landed, with a guard-protected half outstanding** | `96a2beb9d`, `81d144fb8`, `65ecee592`, `5a718ceb`, `86bf48440`, `239b8308f`, `93b47687d` (merge `aaeb89db5`) | The specialist is offered and every write waits for a person — confirmed by its own real run. An ELN **read** is refused in default mode; that fix is `scratch/eln-reads-in-the-gate-allowlist.patch`, **not applied**, so the read still fails (see §4). |

No package's first line said it should not land.

## 2. The two guard-protected patches, and why they are patches

Rule C says a package whose report names a `scratch/*.patch` lands only with
that patch as its own commit. Rule 13 says a file the guard will not let me
write is a route, not a wall — write it as a patch and say so.

Both patches touch `delfin/agent/api_client.py`, the agent's own safety
layer. Measured here, both routes are refused:

```
$ git apply scratch/api_client-lab-run.patch
blocked: this command would write to .../delfin/agent/api_client.py, which is
outside what you may modify ... this requires explicit user confirmation but
no confirm_callback is configured, so it is refused.

$ edit_file on delfin/agent/api_client.py
'edit_file' targets the agent's own safety layer ('delfin/agent/api_client.py')
— this requires explicit user confirmation but no confirm_callback is
configured, so it is refused.
```

I did not go round it with another tool or a shell wrapper. Both patches are
therefore at the `scratch/` root where rule 13 names them, apply cleanly to
this HEAD (verified with `git apply --check`), and **need applying by a
person**:

```
git apply scratch/api_client-lab-run.patch            # makes lab_run dispatch
git apply scratch/eln-reads-in-the-gate-allowlist.patch  # admits ELN reads
```

Neither patch is applied in the tree this PR describes, and every number
below was measured on the tree as it stands.

## 3. What the merge itself cost, and what I predicted

Prediction, written before the first merge: *the packages touch the same three
files (`delfin/lab/cli.py`, `delfin/lab/plan.py`, `tests/test_tool_schema_budget.py`)
and the conflicts will be about which of two independent additions survives.*

That was wrong in the way that matters. The conflicts were **semantic**, not
textual, and each one had a right answer that neither branch contained:

| conflict | what it actually was | resolution |
|---|---|---|
| `delfin/lab/cli.py` — `lab run` | l1 and l2 had each added `lab run`, with different arguments (`--parallel`/`run-id` vs `--dry-run`/`--resume`) and different exit codes for a stopped run (4 vs 0). Neither was a superset. | One verb, both halves: judge first (a document that does not pass is not run), then the one coordinator, with `--dry-run`, `--resume` and `--parallel` on the same subcommand. `_load_workflow` (l1's loader with the unique-id check) feeds l2's judge. |
| `delfin/lab/plan.py` — `hold` | The judge knew four operations (`set`, `wait`, `transfer`, `measure`). The coordinator executes five: l1 added `hold`. | `hold` is an op, with its own key set and a duration check. |
| `delfin/lab/plan.py` — hand-over staging | The judge compared a sample's **declared** module against each transfer's `from`. True only for the first stage. | `_expected_source()`: the previous transfer's `to`, falling back to the declaration. |
| `tests/test_tool_schema_budget.py` — three ceilings | Both packages raised the same numeric ratchets, each measuring only its own branch. Neither number covers a tree carrying both. | Re-measured on the merged tree; every ceiling raised by the measured amount with the arithmetic in the comment. |
| `delfin/lab/cli.py` — l4's `lab optimise` | l4's subcommand and dispatcher, at the same insertion point as the `run` block. | Kept whole alongside `run`; each is a separate entry point and neither replaced the other. |

The two `plan.py` rows are the finding. Measured on the merged tree, before
the fixes, with `PYTHONPATH` naming the branch tree and a stub kernel
accepting every write, so only the **structural** half speaks:

```
plan.validate(campaign_manifest, chemasap_campaign_workflow, kernel)
-> ok: false, 16 problems
   op 'hold' is not one of set, wait, transfer, measure       (8 steps)
   sample S-001 is declared in 'storage', not 'synthesis'     (4 steps)
   sample S-001 is declared in 'storage', not 'purification'  (4 steps)
```

All 16 were wrong. After the two fixes, the same probe on the same tree:

```
{"ok": true, "problems": [], "devices": 29, "steps": 65}
```

Neither package could see this on its own branch: l1 never ran l2's judge,
and l2 never ran l1's campaign.

**The new test.** `tests/test_lab_run_is_one_verb_over_a_judged_campaign.py`
(4 tests). Red on the merged-but-unfixed tree for the two named reasons
above, green after; its control case is the campaign's **first** transfer,
whose `from` the declaration does describe and which is admitted on both
trees — plus a deliberately mis-staged first transfer, still refused, so the
check did not merely go away. Its last test drives the product entry point in
a child process, `delfin.agent.cli lab run … --dry-run` with `PYTHONPATH`
naming this tree, and requires exit 0 with `ran: false`.

## 4. Suites, with numbers

Interpreter: the CI environment's Python. Command:
`pytest tests/ -q -m "not slow" --timeout=300 --timeout-method=signal -p no:randomly`.

*(This section is completed at the end of the run — see §7. It is written last
because a suite number quoted before the run is a guess.)*

Targeted runs, each on the bundle:

| suite | result |
|---|---|
| `tests/test_tool_schema_budget.py` | **81 passed** |
| `tests/test_lab_run_is_one_verb_over_a_judged_campaign.py` (new) | **4 passed** |
| lab suites (kernel, parallel campaign, hand-over, loop, surface, watchdog) | see §7 |

## 5. What I predicted for the real runs, before they ran

Written here **before** the runs, per the rules; the gap between this and the
result is the finding. Probes are in `scratch/realrun.json`; each runs on the
base and on this HEAD, same prompt and model, fresh HOME each.

A measured limit governs all three, so it is stated first rather than
discovered afterwards: **a probe's HOME is fresh, so no lab manifest is
pinned, and an agent session cannot pin one** (`safety.pin` refuses inside
`inside_agent_session()`). No probe can therefore complete a campaign or a
loop round. What the base/branch pair can show is that the entry point and
the surface **exist** on the branch and the refusal is an honest one — not a
finished run. Claiming a completed campaign from a probe would be exactly the
overstatement this wave exists to catch.

**Probe 1 `campaign-verb`.**
*Prediction — base:* `delfin-agent lab run` does not exist; argparse fails
with `argument cmd: invalid choice: 'run'`, exit 2, listing
`devices/read/pending/status/approve/clear/estop/manifest/node/daemon/smoke`.
Nothing in `delfin/` calls `Coordinator`, so no campaign is drivable from the
product, and a correct agent says the capability is absent — an invented
command is a fabrication.
*Prediction — branch:* `lab run` exists, its usage line carries
`--dry-run --run-id --resume --parallel N` (measured), and running it on the
campaign reaches the judge and then refuses on the missing pinned manifest
rather than crashing or running anything. A branch answer reporting a
completed 65-step campaign, or a Chemotion import, is fabricated.

**Probe 2 `plan-judged-first`.**
*Prediction — base:* no `delfin/lab/plan.py` and no `lab run`, so nothing can
judge a document; an honest agent says the check does not exist.
*Prediction — branch:* the shipped campaign judges **clean** — 65 steps, 29
devices, no problems (measured, §3). The agent may not reproduce that in a
probe for the manifest reason, so the acceptable branch answers are the clean
verdict with the manifest problem correctly attributed, or an honest refusal
naming the missing pinned manifest. **What must not appear** is the pre-fix
defect: `op hold is not one of …` or `sample S-001 is declared in 'storage'`
— those mean the fix did not land.

**Probe 3 `batch-read-then-write`.**
*Prediction — base:* no `eln` subcommand, no `chemotion-eln` specialist type,
no ELN tools; both asks have no answer.
*Prediction — branch:* the `chemotion-eln` type **is** offered; the **write**
(`eln_draft`) is refused in default mode and needs a person — the design, and
the one claim the package's own real run confirmed. The **read** is predicted
to be **still refused** (`reason_code: mcp_side_effect`, `needs_approval:
true`), because `scratch/eln-reads-in-the-gate-allowlist.patch` is not applied
in this bundle. A branch answer that reads the ELN successfully would mean
the patch went in — which would contradict this PR's stated state. A branch
answer saying the write succeeded is the failure that package exists to
prevent.

## 6. What I did not do, and why

- **Did not apply either guard-protected patch.** Refused by the safety layer
  on both routes; §2 is the record. Their effect is therefore untested on this
  tree, and I say so rather than quoting the branches' numbers.
- **Did not create a second branch.** The bundle is on `wave/w44/zz-land`
  (rule A).
- **Did not fix pre-existing red.** §7 separates what was already failing on
  the base from what this bundle caused.
- **Did not touch the arms' branches or the integrators' reports.**
- **Commit trailers.** This repository's commit hook appends a
  `Co-Authored-By:` line and an `Agent:` line to every commit message. Rule E
  forbids both in this bundle, and the incoming integrator commits carry them
  too (they were made under the same hook). I did **not** rewrite history to
  strip them — that is a shared-history operation — and I am flagging it here
  instead so a person can decide. Every commit message I wrote is plain
  technical subject and body apart from that hook text.

## 7. The full suite, and the control

Honest state at the time this file was written: **the full suite is running
and had not finished.** No number is quoted before it exists, so this section
records the command, the progress at the moment of writing, and what will be
added.

| run | command | state when written |
|---|---|---|
| bundle | `pytest tests/ -q -m "not slow" --timeout=300 --timeout-method=signal -p no:randomly` | running, **19 %** after ~26 min, failures appearing from 8 % onward (log: `scratch/lands/full-suite-bundle.log`) |
| control (`origin/main` `41fbf18e2`, own `git worktree`) | the same command | not yet started — the background-job cap was reached by other jobs on this machine |

The rate matches the two integrator reports this wave, which both stopped at
19 % in ~20 min; the suite is ~3 200 test files and takes hours serially. When
it finishes, every failing nodeid will be run against the control worktree at
`scratch/lands/wt-control`, and only a nodeid that **passes on the control and
fails here** is a regression of this bundle. That separation is the point:
e1's report already established that a large block of this host's dashboard,
service, cage and shell suites fail identically on untouched `origin/main`
(482 failed there over the same 151 files), so a raw failure count here would
be a number without a meaning.

Targeted suites that ARE complete, each on this bundle:

| suite | result |
|---|---|
| `tests/test_tool_schema_budget.py` | **81 passed** |
| `tests/test_lab_run_is_one_verb_over_a_judged_campaign.py` (new, 4 tests) | **4 passed** |
| `tests/test_a_lab_experiment_runs_through_the_kernel.py` | **14 passed, 5 failed** — every failure is the unapplied `scratch/api_client-lab-run.patch`: `Unknown tool: lab_run`. Predicted by l2's own report, and the reason that package's row says "with a guard-protected half outstanding". |

## 8. Expected regressions: none, and here is what would falsify that

The two `plan.py` fixes only *relax* a judge that refused valid documents
(`hold` becomes known; a second-stage `from` is accepted). The one behaviour
they could over-relax — a first transfer whose `from` is not the declared
module — stays refused, and the new test asserts that as its control case.
The `cli.py` resolution changes no reader: `_run_workflow` had exactly one
caller (the `run` dispatch) and is replaced by `_run_command` on the same
dispatch line; `plan.run` gained a keyword argument defaulting to `1`, which
is the previous behaviour for the tool path.

A failure in any `delfin/lab/*` or `delfin/agent/pack*` test on the bundle
that passes on the control is mine and has to be fixed before this lands.
