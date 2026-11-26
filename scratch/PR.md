# A self-driving bench: a campaign, a closed loop, its ELN hand-over, and the loop's plan judged before anything moves

**DRAFT: this branch does not land as it stands.** 11 tests are red at the
branch HEAD, all of them for one mechanical reason — two product edits live in
`delfin/agent/api_client.py`, which the agent layer refuses to let me write —
and the §2 patches fix exactly those eleven. A person runs the two `git apply`
lines in §2 and they go green; nothing else is red on the branch that is not
also red on `origin/main` (§7).

Five lab packages from wave 44, merged into one bundle, then reviewed four
ways and fixed here. Branch `wave/w44/zz-fix`, based on the lander's
`wave/w44/zz-land`, itself branched from `origin/main` at
`41fbf18e2` and brought up to date with main (`6ce2b6e79`). Every measurement
below is of `bc56dc79b`, **the last commit that changes code**; the commits
after it touch this file only (`git diff --name-only bc56dc79b HEAD` reports
`scratch/PR.md` alone), so a sha or a commit count quoted here would be a
figure this very commit invalidates.

Section 2 is what a person still has to run — three patch files, all
documented with their red and green numbers. Section 3 is the review: every
finding from the four reviewers, one row each, with what was done about it.
Section 7 is the suite, with a control.

## 1. Per package

| package | decision | commits in the bundle | measured effect |
|---|---|---|---|
| **l1** — ChemASAP-shaped parallel campaign | **landed** | `94c88885c`, `4a0c26b4a`, `6620d2e65`, `0ed3acc17` (merge `b2a2fb59c`) | `transfer` moved a sample in the registry and never moved the container, so the campaign could not complete. It now completes: the retained after-logs show `pass 8: rc=0 status=completed done=65`, 4 measurements, ~336 s. **The "before" side is weaker than the package's report claimed** — see review row claims-1: the transcript cited for the 55/65 stop was not retained, and the log named as its source ends on a duplicate-grant `ApprovalError` at `done=41`. The defect itself is real code (`delfin/lab/workflow.py:666` emits `sample mix-up: …`), and 8 of the campaign's 12 hand-overs had no carrier to move into (R1–R4, P1–P4), which the first fix missed and `0ed3acc17` covers. |
| **l4** — closed loop on the simulator | **landed** | `6808cd982`, `02d25b3d1`, `5d1472c9b`, `f654f136b`, `a99f6ffe1` (merge `57b7fd5e9`) | A closed loop reached through `cli.main`, with the watchdog stop neither arm delivered. Its measured claim: EI on a GP reaches 99 % of the surface best in 4 rounds against the lattice's 31. |
| **l2** — the lab specialist | **landed, with a guard-protected half outstanding** | `4a38e2926`, `da5694765`, `e158e62f3`, `7dccca200`, `222edc62a` (merge `ef882588d`) | The judge (`delfin/lab/plan.py`) and the specialist exist. The branch advertises `lab_run` and the dispatcher answers `Unknown tool: lab_run` — measured again on this head. That fix lives in `delfin/agent/api_client.py`, so it ships as `scratch/lab-run-dispatch.patch`, **not applied**; 6 tests stay red until it is (§2). |
| **l3** — run record to Chemotion | **landed** | `16fea816a`, `3c588f298`, `4411f2563`, `662efe664`, `819a3cedf` (merge `a2f16f384`) | A finished run leaves one file a person imports: research plan, reaction with conditions, each sample's analyses under the ELN's own container tree, provenance. The upload request carries the shape the import built. |
| **e1** — Chemotion ELN specialist | **landed, with a guard-protected half outstanding** | `96a2beb9d`, `81d144fb8`, `65ecee592`, `5aa718ceb`, `86bf48440`, `239b8308f`, `93b47687d` (merge `aaeb89db5`) | The specialist is offered and every write waits for a person. An ELN **read** is refused in default mode; the fix is `scratch/eln-reads-in-the-gate-allowlist.patch`, **not applied**, so 5 tests stay red until it is. |

No package's first line said it should not land.

Three defects belonged to the bundle rather than to any one package, and each
was found only here — no package's own branch could see the others' work:

| fix | commit | what it was |
|---|---|---|
| e1's control asserted the lab tool set was exactly five names; l2's merge made it six (`lab_run`) | `6ec548134` | The control's intent — "the ELN work did not disturb the shared dispatcher" — survives with the sixth name recorded; the assertion is still exact equality. |
| the tests named `scratch/lab-run-dispatch.patch`; the file was `scratch/api_client-lab-run.patch` and had never been tracked | `d80dd952e` | The instruction the tests carry could not be followed, and the route did not ship. |
| a wave test started a DELFIN child on the raw environment | `fd603b15f` | The child resolved its `~/.delfin` sinks against the real home; caught by the repository-wide scanner, fixed with `conftest.child_env`. |

Four more came from this review turn: the two token pins (`48f386d0d`),
the text channel's unenforced `max_step` (`25a4ed7fb`), the unguarded read
in `eln_attach` (`fbd812f05`), and two assertions that could not fail
(`a609ae6eb`).

A fifth was found by the widened control pair and by no reviewer: a
key-shaped fixture in `tests/test_eln_tools_are_wired.py` wrote a genuine
PEM private-key header into a tracked file, and
`tests/test_no_secrets_are_committed.py` was **red on this branch and green
on main** because of it (measured: 1 failed 33 passed vs 34 passed). The
refusal that test exercises never reads the bytes — `eln_attach` matches the
path against `*.key` (`delfin/agent/api_client.py:9586`) — so the header
proved nothing. Fixed in `bc56dc79b`; guard clean, 52 passed over the two
files.

## 2. The three guard-protected patches, and what a person must run

Rule C says a package whose report names a `scratch/*.patch` lands only with
that patch as its own commit. Every one of these three touches the agent's
own safety layer, which is in `_DEFAULT_PATH_PROTECTED_GLOBS`
(`delfin/agent/api_client.py:9860`, `:9872`). Measured on this tree, both
routes are refused:

```
$ git apply scratch/lab-run-dispatch.patch
blocked: this command would write to .../delfin/agent/api_client.py, which is
outside what you may modify ... this requires explicit user confirmation but
no confirm_callback is configured, so it is refused.

$ <edit_file on delfin/agent/api_client.py>
'edit_file' targets the agent's own safety layer ('delfin/agent/api_client.py')
— this requires explicit user confirmation but no confirm_callback is
configured, so it is refused.
```

I did not go round it with another tool, a shell wrapper or a second
interpreter. The three patches are at the `scratch/` root, apply cleanly to
this HEAD (`git apply --check`, all three rc=0), and **need applying by a
person**:

```
git apply scratch/lab-run-dispatch.patch               # makes lab_run dispatch
git apply scratch/eln-reads-in-the-gate-allowlist.patch # admits ELN reads
git apply scratch/eln-approve-as-command-text.patch     # refuses `eln approve` as text
```

| patch | what it changes | measured before → after |
|---|---|---|
| *(the two rows below marked † are the red side as measured this turn; their green side was measured by the earlier session that copied the module, not re-run here — the protected file cannot be written to re-measure it)* | | |
| `lab-run-dispatch.patch` † | `_LAB_TOOL_NAMES` (:12347) gains `lab_run`; the `_DocToolExecutor` dispatch tuple (:18516) gains it; `_is_honoured_dry_run` (:8977) learns `lab_run`'s argument-conditional `dry_run` | red side re-measured this turn: 11 failed over the three files, of which the `lab_run` share is 6, with `advertised with no executor and no dispatch: lab_run` → green (earlier session) |
| `eln-reads-in-the-gate-allowlist.patch` † | the read-only allow-list gains `eln_find`, `eln_read`, `eln_pending` | red side re-measured this turn: the 5 `test_an_eln_read_is_reachable_and_its_writes_are_not.py` cases → green (earlier session) |
| `eln-approve-as-command-text.patch` | `answer_command` matches `eln approve\|deny\|answer` as well as `approvals …`, and a Python `grant(` beside its module | 6 of 11 new cases red → 11 green (measured on a copy of the module; the tree is untouched) |

**What is measured on this head, per file** (`-m "not slow"`, CI
interpreter, `-p no:randomly`):

| file | result |
|---|---|
| `test_a_lab_experiment_runs_through_the_kernel.py` | **5 failed, 22 passed** (27 collected) |
| `test_an_eln_read_is_reachable_and_its_writes_are_not.py` | **5 failed, 7 passed** |
| `test_a_tool_reaches_the_function_behind_it.py` | **1 failed, 6 passed** |

So the branch is **not green as it stands**, and it is not green for one
reason: three of its own guards assert a capability whose enabling line is
in a patch a person has to apply. That is stated here rather than hidden.
The 5 + 5 + 1 split matters (review row claims-4): six of the eleven need
`lab-run-dispatch.patch` and five need `eln-reads-in-the-gate-allowlist.patch`,
so a reader who applies only the first still sees five red and should not
read that as the patch having failed.

**The lander's §2b numbers were wrong** and are not repeated here. It
claimed "46 passed, 10 failed" unpatched and "56 passed, 0 failed" patched
over these three files. Those files collect 46 tests in total, so 56 passed
is impossible, and the measured count at this head is **11 failed, 35
passed**. The patches add no test file, so the correct end state is 46
passed, 0 failed.

## 3. The review, one row per finding

Four reviewers read the bundle (behaviour, claims, safety, tests) and
returned 35 findings. Each row below was re-checked here before anything was
done about it — the command or the caller is named in the "confirmed?" cell.
"blocks" findings are first inside each reviewer's block.

### 3a. Behaviour review (6)

| # | sev | what | confirmed here | what was done |
|---|---|---|---|---|
| b1 | blocks | `tool_catalog.py:454` advertises `lab_run` with no dispatch | **yes**: `test_a_tool_reaches_the_function_behind_it.py` fails with `advertised with no executor and no dispatch: lab_run` | patch shipped, not applied (§2) — needs a person |
| b2 | should fix | `server.py:2381`: 3 ELN reads registered but not classified read-only | **yes**: `test_the_reads_are_not_a_back_door_to_the_writes` and `test_every_eln_tool_the_server_registers_is_classified` both fail | patch shipped, not applied (§2) |
| b3 | blocks | prompt-budget row 11757 below the merged head | **yes**: the head measured 11776 | fixed `48f386d0d`; re-measured green (19 passed over the two pin files) |
| b4 | blocks | catalogue pins re-taken from one branch | **yes**: three pins disagreed with the sibling budget file | fixed `48f386d0d` (12170 / 12080); green |
| b5 | should fix | a writable `text` channel declares `max_step`, the kernel silently drops it | **yes**: `manifest.py` kept it; `safety.py` skips text channels | fixed `25a4ed7fb` — the test raised `DID NOT RAISE ManifestError` before, passes now; no shipped manifest affected |
| b6 | should fix | the lander's §7 numbers are a 93-file subset | **yes**: the logs' own headers say `FILES=93`; the "full" logs hold no result | §7 of this text says what it is and what it is not |

### 3b. Claims review (9)

| # | sev | what | confirmed here | what was done |
|---|---|---|---|---|
| c1 | blocks | l1's "before" row (55/65 steps, `measure-2`) has no artifact behind it | **yes**: `grep -rn 'done=55' scratch/` matches only the review; `scratch/arms/l1-int/campaign_final.log` ends on `ApprovalError: lab-f0ba876b is granted` at `done=41` | §1 row now says the observation is reported without a retained transcript, and quotes what the log really ends on |
| c2 | blocks | "46 passed, 10 failed" / "56 passed, 0 failed" are impossible | **yes**: the three files collect 46 tests; measured 11 failed, 35 passed | §2 states the measured counts and the correct end state |
| c3 | should fix | §7 presents a 93-file rerun as the full suite | **yes**, same evidence as b6 | §7 relabelled; the "fixed by the wave: 0" claim scoped to the pair actually run |
| c4 | should fix | "the remaining 11 all fail with one message" | **yes**: 6 are `lab_run`/`Unknown tool`, 5 are ELN-gate failures with different text | §2 splits them by patch |
| c5 | should fix | the tracked `scratch/PR.md` was a stale draft naming a non-existent patch | **yes**: `git show HEAD:scratch/PR.md` differed from the lander's text and named `scratch/api_client-lab-run.patch` | this file replaces it; §7 is measured, not "not yet started" |
| c6 | should fix | `scratch/lands/pr-section-7.md` names a file that does not exist and counts a 27-test file as 19 | **yes**: `git cat-file -e HEAD:scratch/api_client-lab-run.patch` fails; measured `22 passed, 5 failed` | untracked in `9fd4f72ce` rather than kept as a second copy of §2 |
| c7 | should fix | the hash `5a718ceb` resolves to no object | **yes**: `git cat-file -t 5a718ceb` → `Not a valid object name`; `5aa718ceb` is a commit | corrected in §1 |
| c8 | note | e1's report quotes 502 failures / ~8500 passed with no retained run record | **yes**: those logs are bare progress streams and their per-chunk table went to a `mktemp -d` | not rewritten — it is the arm's own report, preserved unmodified; recorded here instead |
| c9 | note | §6 blames a repository commit hook for the trailers | **yes**: `core.hooksPath` and `commit.template` are unset; the hooks dir holds only `*.sample` | §6/§9 reworded: the trailer is injected by the agent layer, not by a hook |

### 3c. Safety review (9)

| # | sev | what | confirmed here | what was done |
|---|---|---|---|---|
| s1 | blocks | the ELN write gate rests on the `DELFIN_AGENT_SESSION` variable alone; `eln approve` is not refused as command text | **yes**: `answer_command("delfin-agent eln approve eln-abc123")` is `None`, while `approvals approve` is refused | `AnswerRefused` closure shipped as `scratch/eln-approve-as-command-text.patch` with a test that is red without it (11 cases: 6 red before, 11 green with the patch applied to a copy). Not applied: protected file. The durable fix belongs with `delfin/lab/approvals.py`, which has the same single-layer shape and is identical on `origin/main` |
| s2 | should fix | `eln_cli.py`'s docstring claims a command-text gate that does not exist | **yes**, same probe | docstring corrected in `65daae37f` to the measured single-layer truth plus the residual risk |
| s3 | should fix | `eln_attach` reads an arbitrary absolute path before approval, with no deny-list | **yes**: no `_is_secret_path`/deny-list call existed on that path | fixed `fbd812f05` — deny-listed paths refused before the read, with a test |
| s4 | blocks | four fast-suite tests fail on the branch and pass on main (the pins) | **yes**: reproduced, then green after the fix | fixed `48f386d0d`; the remaining patch-dependent red is b1/b2/t3 |
| s5 | should fix | a live proxy password is printed in a captured traceback in the package reports | **not confirmed as a credential**: it is the loopback egress proxy's own generated session token (`delfin/agent/egress_proxy.py:1023` compares it against `self.token`), and the log is git-ignored (`.gitignore:178` → `scratch/`) | nothing changed; not a secret, and not in the tree. Worth a redaction in `conftest`/pytest's repr if reports are ever shipped |
| s6 | note | `eln_draft`/`eln_attach` are not `_host_gated`, deliberately | confirmed sound: one approval per write is the gate | none |
| s7 | note | `lab_run` advertised while its dispatch is a patch | confirmed; identical to b1 | recorded open (§2) |
| s8 | note | every commit carries `Co-Authored-By: DELFIN-Agent (<model>-at-<institution>)` and `Agent: w44/…` | **yes**: 56 `Co-Authored-By` lines and 45 `Agent:` lines over 46 commits (`git log origin/main..HEAD --format='%b'`); the same trailer form is on `origin/main` | not rewritten — see §6. Every commit I wrote had no trailer in its message file (`grep -c` = 0) and the trailer appeared in the commit anyway, which is the mechanism |
| s9 | note | new outbound network path, examined and found sound | confirmed: one transport seam, redacted URLs, no TLS verification disabled, no host in the new example manifests | none |

### 3d. Tests review (11)

| # | sev | what | confirmed here | what was done |
|---|---|---|---|---|
| t1 | blocks | three characterisation pins fail on the branch | **yes**: 7 → 4 failed before the fix | fixed `48f386d0d`; 19 passed now |
| t2 | blocks | `test_prompt_token_budget.py:297` row 11757 fails (head 11776) | **yes** | fixed `48f386d0d` |
| t3 | blocks | 10 tests red because two product edits ship as patches | **yes**, and it is 11, not 10 (the lander's 12th file is `test_a_tool_reaches_the_function_behind_it.py`) | §2 — needs a person to apply two patches |
| t4 | should fix | `assert k.release(me) or True` cannot fail | **yes**: true for every value | replaced with the released-locks assertion in `a609ae6eb` |
| t5 | should fix | the six ELN write cases assert only `out is not None`, which a deny-by-default refusal satisfies for an unknown name | **yes** | strengthened in `a609ae6eb`: classification plus the refusal naming the tool |
| t6 | should fix | a docstring says `lab`'s dashboard half "arrives as `scratch/lab-run-dispatch.patch`" — neither patch carries it | **yes**: both patches touch only `api_client.py`; `tab_agent.py` is untouched on this branch | corrected in `9058e5e43`: the picker lists still read four modes, nothing here ships `lab` there, nothing tests its absence |
| t7 | should fix | the lander's §7 is a 93-file rerun presented as the suite | **yes** | see §7 |
| t8 | note | a second `... or True` in the catalogue test, false condition | **yes**: `api_client` is imported at module level, so its condition was always false | removed in `48f386d0d`, with the real AST check above it named as the check that stays |
| t9 | note | four tests in `test_the_tool_surface_is_exercised.py` cannot measure the search surface on this host | **yes**: three skip with `gaierror`, one fails on **both** trees | none required; recorded |
| t10 | note | the parallel-campaign test asserts on wall-clock | confirmed genuine, with its own sequential control | left as is |
| t11 | note | `tests/_lab_support.py` builds the executor with `__new__`, skipping `__init__` | **yes**: the `__new__` instance has none of the six attributes; the product builds `_DocToolExecutor()` | fixed `0c7c22271`; every file using the helper re-run (201 passed, the same 5 `lab_run` reds) |

**Totals:** 35 findings — **10 blocks, 15 should-fix, 10 notes** (behaviour
3/3/0, claims 2/5/2, safety 2/3/4, tests 3/4/4). Confirmed by re-running the
check: **34 of 35** (the exception is s5, which is the harness's own loopback
proxy token, not a credential). Several findings are the same defect seen by
two reviewers and are fixed once: b3/b4/t1/t2/s4 are one commit
(`48f386d0d`), and b1/s7 are one patch. In summary — the pin regression, the
text channel, the `eln_attach` path, the two dead assertions, the docstring
routes and the test helper are fixed in commits on this branch; the
`lab_run` dispatch, the ELN read allow-list and the `eln approve` command
text stay in shipped-but-unapplied patches (§2); the rest are corrected in
this text and in the files it cites, or left as they are with the reason
given in their row.

## 4. Suites, with numbers

Interpreter: the CI environment's Python (CPython 3.11.16, pytest 9.1.1).
Command shape: `PYTHONPATH=. pytest <files> -q -m "not slow" --timeout=300
--timeout-method=signal -p no:randomly`.

| suite | result |
|---|---|
| `test_the_tool_catalogue_is_the_same_after_the_move.py` + `test_prompt_token_budget.py` + `test_tool_schema_budget.py` | **100 passed** |
| the four lab suites touching the files this turn edited | **55 passed** |
| every file using `tests/_lab_support.py` | **201 passed, 5 failed** (all five are the `lab_run` cases of §2) |
| the three patch-dependent files | **11 failed, 35 passed** (§2, per file) |
| `test_no_secrets_are_committed.py` + `test_eln_tools_are_wired.py`, after `bc56dc79b` | **52 passed** (before it: the guard was 1 failed 33 passed) |
| `ruff check delfin/ tests/ --select F821,F811,F601,F402,F823` | **clean** (`All checks passed!`; one pre-existing `# noqa` format warning in `tests/test_branch_dropdown_options.py:63`) |

## 5. What I predicted for the real runs, and what they measured

Written before the runs; the probes are in `scratch/realrun.json`. Same
three probes as the lander's, so the pair is comparable; each runs on the
base and on this HEAD, same prompt and model, fresh HOME each.

**What actually happened** (`scratch/realrun/summary.md`): *this branch has
not had its own real-run pair yet.* The lander's pair ran on
`wave/w44/zz-land`; a branch's own pair is what the wave's rule asks for, so
`scratch/realrun.json` asks for the same three probes again on this HEAD, and
this section is the prediction, written before them. What follows is what the
change should make visible; the measured answer replaces it when the runs
return.

**The measured limit that governs all three.** A probe's HOME is fresh, so no
lab manifest is pinned, and an agent session cannot pin one. Measured on this
tree, not assumed: `lab run <campaign> --dry-run` answers
`manifest refused: cannot read manifest.yaml: [Errno 2] No such file or
directory: <home>/.delfin/lab/manifest.yaml` and exits 1. So no probe can
complete a campaign, and the predictions below are written to that fact.
What separates base from branch is that the verb and the judge **exist** and
refuse honestly by name — not that a run finishes.

**Probe 1 `campaign-verb`.**
*Prediction — base:* `lab` has no `run` verb. Measured by building the
base's parser in a control worktree: its subcommands are approve, clear,
daemon, devices, estop, manifest, node, pending, read, smoke, status. A base
answer that quotes a usage line containing `run`, or reports anything `run`
printed, is fabricated.
*Prediction — branch:* the verb exists. Measured here: `run`'s options are
workflow, `--dry-run`, `--run-id`, `--resume`, `--parallel`. Run on the
shipped campaign it prints the manifest refusal above and exits 1. The
acceptable branch answers are that refusal correctly attributed, preceded by
the real usage line. **What must not appear** is a fabricated 65-step
completion or a Chemotion import.

**Probe 2 `plan-judged-first`.**
*Prediction — base:* no `delfin/lab/plan.py` (checked: absent in the control
worktree) and no `lab run`, so nothing can judge a document; an honest agent
says so. A base answer quoting a problem count has invented it.
*Prediction — branch:* `plan.py` is there and `lab run` judges before it
runs, but the judge is reached behind `Kernel.load()`, which refuses the
missing pinned manifest first — so the acceptable answers are the manifest
refusal correctly attributed, or the judged verdict if the agent reaches the
judge another way. **What must not appear** is the pre-fix defect the
bundle's two plan fixes removed: `op hold is not one of set, wait, transfer,
measure` or `sample S-001 is declared in 'storage', not synthesis`. Those
strings mean the judge is answering and still wrong.

**Probe 3 `batch-read-then-write`.**
*Prediction — base:* no `eln` subcommand (`delfin/lab/eln_cli.py` does not
exist on the base; `delfin/agent/cli.py` has no `eln` parser there) and no
`chemotion-eln` type, so both asks have no answer.
*Prediction — branch:* the subcommand and the type exist, and the two
outcomes differ **on purpose**. The write (`eln_draft`) is refused in
default mode and needs one person's approval. The read is the known gap, and
its state here is measured, not guessed: `eln_find`, `eln_read` and
`eln_pending` are still absent from `_MCP_READONLY_TOOL_BASES` (a grep of
`delfin/agent/api_client.py` finds no ELN name there), because the edit that
adds them is `scratch/eln-reads-in-the-gate-allowlist.patch`, shipped
**unapplied**. So the read is predicted to be **still refused**, with
`reason_code: mcp_side_effect` and `needs_approval: true`. A branch answer
that reads the ELN successfully would mean the patch went in, which would
contradict this text; an answer that says the write succeeded is the failure
that package exists to prevent.

## 6. What I did not do, and why

- **Did not apply any of the three patches.** Both write routes are refused
  by the harness for `delfin/agent/api_client.py` and
  `delfin/agent/approval_guard.py`, and a missing confirmation is not a wall
  to go round with another tool, a shell wrapper or a second interpreter.
  §2 is the record, and the branch does not land until they are commits.
- **Did not rewrite the bundle's history to strip the trailers**, though this
  turn's brief asks for it, and this is the one item of the brief left undone.
  **Measured this turn, not assumed:**
  - A commit message I wrote and verified trailer-free on disk
    (`tail -c` showed the file ends at `…are 52 passed.\n`) still came out of
    `git commit -F` with both attribution lines appended.
  - `git commit --amend --no-verify -F <same clean file>` appended them again,
    so no hook is involved and `--no-verify` cannot reach it.
  - A **throwaway repo** (`scratch/fix-probes/probe-repo`, `git init`) gave the
    same result for `-m` and for `-F`, so the injection is not specific to this
    checkout: it is applied by the execution wrapper around every `git commit`.
  - `core.hooksPath`, the repo's `hooks/` directory (only `pre-push`), the
    global and local git configs and `init.templateDir` carry nothing that
    adds these lines.

  The other route — rewriting the bundle's messages — is refused by the sandbox:
  both `filter-branch` and `filter-repo` are on its deny list, including when a
  helper script merely names them, and the git rules reserve a shared-history
  rewrite for the user in any case. Trying to reach it another way (a second
  spelling, a wrapper, a hand-built commit object) is exactly the route around
  a refusal that must not be taken, so it was not taken. **What this needs:**
  the setting that appends the lines switched off in the agent layer
  (`agent.commit_coauthor` in the user's settings, which
  `_commit_trailer_enabled()` reads), then a maintainer-decided message rewrite
  of this branch — the situation the earlier waves already left backup branches
  for. Recorded here so the decision is visible and made by a person, not
  quietly taken here.
- **Did not commit the widened pair's logs.** They are
  `scratch/fix-probes/{all146.txt,set-branch.txt,set-control.txt,*-all146.log}`
  and are git-ignored; §7 quotes every number from them. The control worktree
  is `/tmp/delfin-wt-5ab4647e` and is outside the branch.
- **Did not delete the reviewers' untracked logs or the untracked artifact
  folders.** They are git-ignored and are other sessions' evidence; the one
  credential-shaped string in them is the harness's own proxy token (row s5).
- **Did not fix pre-existing red.** §7 separates it from what this bundle
  caused, per node.
- **Did not touch the arms' branches or the integrators' reports.**

## 7. The suite, and the control

Two runs, serially, one pytest process at a time.

1. **The branch, the whole fast suite, chunked into 16 file-list slices** —
   one pytest holds several GB of resident modules on this host, so a single
   whole-suite process was killed at 51 % by the memory budget on the
   lander's side as well. The 16 slices are in `scratch/fix-probes/suite/`
   with their per-slice results in `index.txt`. **A chunked run is not a
   whole-suite total**: the slices sum to more tests than a single session
   collects, because some files parametrise differently under a different
   selection, so no whole-suite count is quoted from them.
2. **A pair over the 146 files that held a failing test in the chunked run.**
   This supersedes the 85-file pair below, which was the first, narrower
   attempt: the 85 were the files whose tests failed *under the 85-file
   selection*, and re-running the whole fast suite showed failures in more
   files, so the pair was widened to every file that held a failing test in
   the 16 chunk slices. Both runs used the identical file list
   (`scratch/fix-probes/all146.txt`, verified identical in both trees), the
   same flags, interpreter and host. **The control is a `git worktree` of main
   at `62fdab2d8`** — the commit `enter_worktree(base_ref="origin/main")`
   resolved to when it was created. Today's `origin/main` is `80a43c7ff`;
   the control is therefore **7 commits behind it** (waves 45–50), not one as
   the earlier text said. Stated because it matters for reading the pair: the
   control carries neither the bundle nor those 7 main commits. The comparison
   is a failure **set**, not a count.

   | run | total |
   |---|---|
   | control, 146 files, `62fdab2d8` | 476 failed, 2087 passed, 23 skipped, 37 errors in 3931 s |
   | branch, same 146 files, HEAD `bc56dc79b` | 476 failed, 2087 passed, 23 skipped, 37 errors in 3816 s |

   | relation (unique nodeids) | count |
   |---|---|
   | failing on both — pre-existing on main, not this wave | **509** |
   | failing on the branch only | **2** |
   | failing on the control only | **2** |

   The four asymmetric nodeids, each re-run alone on both trees:

   | nodeid | branch, isolated | control, isolated |
   |---|---|---|
   | `test_a_tool_reaches_the_function_behind_it.py::test_every_advertised_tool_can_actually_be_dispatched` | **red** — `advertised with no executor and no dispatch: lab_run`, i.e. §2's unapplied patch, not a new defect | green |
   | `test_no_secrets_are_committed.py::test_no_tracked_file_carries_a_generic_secret` | **red before `bc56dc79b`**, green after; the key-shaped fixture | **green** (34 passed) |
   | `test_agent_sandbox.py::test_allowlist_accepts_safe_commands[grep -r foo .]` | **red** | **red** |
   | `test_agent_sandbox.py::test_allowlist_accepts_safe_commands[find . -name '*.py']` | **red** | **red** |

   So the widened pair leaves exactly one branch-only red: the `lab_run`
   dispatch that waits for a person to apply `lab-run-dispatch.patch` (§2).
   Nothing else on this branch fails that does not also fail on main. The two
   sandbox parameters fail on both trees and merely happened not to be
   exercised in the branch's file list during the paired run; they are part of
   main's 509. This host's dashboard, service, sandbox and cage suites are
   order- and load-sensitive, which is why the pair is compared as a set.

   **The branch's own contributing-file checks, all green:** the secret guard
   is clean (`scripts/secret_guard.py --staged`, run with `PYTHONPATH=.` so it
   reports this tree and not an installed DELFIN), and ruff with CI's gated
   selection is clean.

**For the record — the earlier, narrower 85-file pair.** The 85 files were the
files whose tests failed under the 85-file selection. Its control-only and
branch-only rows are superseded by the 146-file pair above; kept because the
wording of the lander's §7 referred to them.

| relation | count |
|---|---|
| failing in the 85 files on the control (`origin/main`) | **275** |
| failing in the same 85 files on the branch (HEAD `9fd4f72ce`) | **274** |
| failing on both — pre-existing on main, not this wave | **272** |
| failing on the branch only | **2** |
| failing on the control only | **3** |

Neither of the two branch-only nodeids is a defect this wave introduced —
both were re-run alone, and one of them is the known §2 item:

| nodeid | isolated re-run on the branch |
|---|---|
| `test_a_tool_reaches_the_function_behind_it.py::test_every_advertised_tool_can_actually_be_dispatched` | **red, 1 failed** — this is b1: `lab_run` advertised, its dispatch in the unapplied patch (§2). It is the same failure the patch fixes, not a new defect. |
| `test_a_typed_result_fits_the_window.py::test_the_archive_tables_fit_the_legacy_cap[<lambda>1]` | **green, 3 passed** (whole parametrisation) — order-dependent in the shared file list |

The three control-only nodeids are the mirror image, and every one of them
was re-run on **both** trees afterwards:

| nodeid | branch, isolated | control, isolated |
|---|---|---|
| `test_agent_sandbox.py::test_allowlist_accepts_safe_commands[grep -r foo .]` | **red** | **red** |
| `test_agent_sandbox.py::test_allowlist_accepts_safe_commands[find . -name '*.py']` | **red** | **red** |
| `test_g1_degenerate_loop_fires_in_a_run.py::test_a_changing_failure_set_runs_on` | **green** | **green** |

So the two sandbox parameters fail on both trees (they are part of main's
272, and the branch's file list simply happened not to exercise that
parameterisation during the paired run), and the `g1` test passes on both
when it has the machine to itself. That is the whole of the difference
between the two trees — this host's dashboard, service, sandbox and cage
suites are order- and load-sensitive, which is why the pair is compared as a
set and not as a count.

**Scope of the pair, stated rather than implied:** the 85 files are the
files that failed on the branch. A test that fails on the control and passes
on the branch, inside a file that never failed on the branch, was not re-run
on the control. So the "control only" column is what this pair can see, not
a claim about the whole suite.
