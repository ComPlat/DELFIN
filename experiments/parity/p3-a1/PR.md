# PACKAGE p3 — a review finds real defects with evidence and reports nothing it did not verify — integration of four arms

**Should land: YES, as integrated here** — one branch (a1 × 3 commits + a4 × 5 commits + one integration fix), every gate measured green; a2 and a3 contributed nothing that lands.

## Items

| item | decision | commit | reason (one line) |
|---|---|---|---|
| a1 harness + before.json (12 variants, drives `cli.main`) | take | 46bb098e9 (cherry-pick of f29261cff) | Measured-first harness, before.json committed before product code changed. |
| a1 findings contract (`verify_review_findings`, `REVIEW_FINDINGS_SCHEMA`, `/review` fix) | take | a3b8b7a5a (b0e510576) | 0/12 → 12/12 REACHED/USED/CORRECT; red test on main, green here. |
| a1 report + check scripts + probes | take | 78cde6c1e (411d56d22) | Scratch provenance: suite.sh, check.sh, realrun.json. |
| a4 seeded-defect harness, 24 variants, drives `land.main` | take | 8cefe72f0 (06f2bdb95) | Only `land.kit_chat` faked; the real entry point forwards to the code under test. |
| a4 `diff_places`/`verified_findings` in `land.py` | take | f90940a4e (666e41939) | One reader of the diff, used by the gate; no second diff parser. |
| a4 red test `test_a_review_reports_only_findings_it_could_verify.py` | take | 6882440e2 (419b58086) | Red on main (measured), green here. |
| a4 unmounted-capability anchor + re-measured after | take | 8931d6ec0, 81dfb378e (947bc2f60, 1ee581992) | Keeps `delfin-agent metrics` honest about the moved anchor. |
| a3 `delfin/agent/review_findings.py` | **reject** | — | Second implementation of what a1 already built and wired: a3's module is imported by nothing, its wiring into `run_subagent` was never finished (its own PR.md says so). Rule: never add a second implementation of something that exists. |
| a2 (measurement-only, no commits) | **reject as code; kept as control evidence** | — | Its control numbers (165 passed / 1 pre-existing fetch failure) matched my control run exactly; no product change to take. |
| integration fix: home path in after.json | take | fad08e7eb | The personal-data scanner caught `measured_tree: /home/...` in a1's after.json — tree named by commit instead (rule 16). |

## Prediction (rule 3, written before the measurement)

I predicted the consuming suites would pass (the arms' own reports said 255 passed on a1's check set, and a4's gate tests were green on its branch) and that the whole-tree run would show the same 7 uncollectable cluster imports both arms recorded. I predicted the fetch-page failure would appear unchanged, because both a2 and a1 called it pre-existing on main — and I predicted the personal-data and licence scanners would pass, because the arms' scratch looked clean.

**Gap between prediction and data:** one. The personal-data scanner FAILED on the branch and passed on main — the first real defect the integration found: a1's `after.json` recorded the arm's worktree as an absolute path under a home directory, and the username inside that path hashed to a colleague's name (the scanner's own `_findings` named both offenders). Fixed in fad08e7eb by naming the tree by commit instead; scanner now 18/18 passed. Everything else matched prediction.

## Control run on untouched origin/main (step 1)

`git worktree add .control-base origin/main`, editable install from that tree, `LOADED_FROM` printed as `.control-base/delfin/agent/subagents.py` before anything was counted:

| suite | origin/main (control) | branch HEAD |
|---|---|---|
| a1 red test + a4 red test | **32 failed, 3 passed** (red for the reason each names: no contract, unverified findings reach the parent, phantom files reported) | **35 passed** |
| consuming suites (sa0, presets, subagent contracts, slash commands, review gate) | 1 pre-existing fetch failure | 275 passed, 1 failed (same fetch test), 3 skipped |
| `test_the_tool_surface_is_exercised` fetch test | **FAILED on main too** (control proof) | same failure — pre-existing, unrelated |
| repo-wide scanners (secrets, personal data, licences, entry points, tool surface, presets, API route) | personal-data scanner **passed** on main | 100 passed, 4 failed — of which 3 are pre-existing (fetch ×1, api-route ×2 — both control-proven on main), 1 was mine (fixed) |

**Metrics gate:** `delfin-agent metrics record --tree .` on branch HEAD: 6085 callables defined, `newly_unreachable: []`, `unmounted: 11`, `within_ceiling: true` — no new unreachable callable (a2's baseline said 6079 callables / 11 unmounted; ours adds 6 defined, 0 unreachable).

**Whole-tree suite (rule 14):** NOT run to completion — said plainly. Foreground ceiling is 600 s, all 8 background slots are held by other wave workers, and 16 xdist workers exceeded the 8 GB RSS cap (process budget killed it at 375 s). What ran instead: the consuming suites (above), the repo-wide scanner suites, and a1's three chunk logs exist under `scratch/arms/p3-a1/` from its own run. The 7 `--ignore`d cluster tests cannot collect on this host; same on main.

## What I did NOT do, and why

- Did not take a3's module: it duplicates a1's contract and has no caller (rule 7 — an import is not a caller). Its wiring recipe sits unfinished in its PR.md.
- Did not re-run a2's numbers: a2 committed nothing; its report's control numbers were instead *re-derived independently* (my control run reproduced its 1 pre-existing failure exactly).
- Did not run the full `tests/` to completion (constraints above) — the scanners that rule 14 names all ran, and the gap is stated rather than hidden.
- Did not push or open a PR: delivery convention here is PR-and-stop; this PR.md (in `experiments/parity/p3-a1/`, beside the harness that produced these numbers) is the deliverable.

## Real-run predictions (rule: prediction before the runs)

Written before ending the turn; probes live beside this file in `realrun.json`:

- **p1-seeded-defect-review** — Prediction: on main the reviewer answers in prose with no file:line and no scenario (REACHED=false for the contract, nothing to count); on the branch every finding carries file:line + failure scenario + verdict, and the withheld count is visible. Recall on the 6 seeded defects ≥ 4/6, precision ≥ 0.8, because the seeded defects are exactly the classes a1's contract was measured on.
- **p2-withheld-phantom-file** — Prediction: main invents or silently omits the unreadable file; branch reports `[review-findings] N finding(s) withheld` and keeps only checkable findings. Falsifier: a finding citing the unreadable path.
- **p3-clean-diff-no-findings** — Prediction: main writes prose from which "nothing found" is indistinguishable from "did not look"; branch returns empty findings with status `no_findings`. Falsifier: invented findings for a clean diff (precision 0).

Each probe runs on base and branch HEAD, same prompt and model, fresh HOME. A probe counts as shown-working only if the branch answer differs on the contract's terms (file:line, scenario, withheld count) — not on prose style.
