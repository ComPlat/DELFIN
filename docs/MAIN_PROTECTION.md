# How main is protected

`main` is the branch every DELFIN user installs from, and on shared clusters
a change to it reaches calculations other people are running. It is
protected by four GitHub rulesets. Their configuration lives on GitHub, not in
this repository; this document is the reference for what they hold, why, and
the one sanctioned way around them.

Current state, as GitHub enforces it:

    gh api repos/ComPlat/DELFIN/rules/branches/main

## The rulesets

| Ruleset | Holds | Bypass |
|---|---|---|
| `main-hard-limits` | `main` cannot be deleted or force-pushed | none |
| `main-required-checks` | the required checks below, on a branch that contains the current `main`; CodeQL code scanning | none |
| `protect-main` | a pull request, approved by a code owner (`.github/CODEOWNERS`) after its last push, with every review thread resolved | repository admins |
| `protect-release-tags` | `refs/tags/v*` cannot be moved or deleted | repository admins |

The split between `main-required-checks` and `protect-main` is deliberate.
A review requirement that an admin may bypass still leaves an audited
decision by a person who read the change. A required check that an admin may
bypass leaves nothing: the merge is recorded as passing a gate it never ran.
So the checks hold for admins too, and stepping around them is the explicit,
logged procedure below rather than a flag on a merge command.

## The required checks

`strict_required_status_checks_policy` is on: a pull request is mergeable
only when its branch contains the current `main`, so the combination that
lands is the combination that was tested. Two pull requests that pass
separately can fail together; this happened on 2026-10-15 (#273 and #274,
each green against an older `main`, together broke a test on `main`).

| Check | Workflow | What it establishes |
|---|---|---|
| `tests (py3.11)` | ci.yml | the fast suite on the interpreter production runs |
| `lint (ruff)` | ci.yml | no undefined name, shadowed definition or unreachable branch |
| `container isolation (docker, podman)` | ci.yml | the agent's container backends hold their boundary |
| `generalist toolchains (node, cargo, go)` | ci.yml | the agent's test runners for other languages run inside the cage |
| `reference environment (py3.12, numpy 2.2)` | ci.yml | the code runs on the newer stack users install |
| `analyze (python)`, `analyze (actions)` | codeql.yml | static analysis of the code and of the workflows |
| `dependency review` | dependency-review.yml | no dependency with a high-severity advisory is added |
| `window and service` | tui.yml | the terminal window and the agent service it talks to |

The `code_scanning` rule additionally refuses a merge while CodeQL reports a
new alert of severity *error* or a security alert of *high* or *critical*.

Not required, on purpose:

- `dependency audit` (audit.yml, nightly and on demand) audits the closure
  a user installs (`.[agent,docs]`) and is red while an advisory is open
  and not recorded in `.github/ci/audit-ignore.txt`. A vulnerability
  published upstream turns it red on every branch at once, unrelated to the
  change under review; as a required check, or as a job in ci.yml, whose
  conclusion is the merge verdict, it would stop all merges until upstream
  ships a fix. `dependency review` is the change-scoped half and is
  required.
- `binary (...)` (tui.yml) builds the window for four platforms, two of them
  on macOS runners that are scarce; it reports, it does not gate.
- `slow-tests.yml` and `sandbox-macos.yml` run on a schedule or on demand.

A required check must be produced on every pull request. A workflow with a
`paths:` filter on `pull_request` never starts for a pull request outside
those paths, and GitHub then waits for the check indefinitely. tui.yml
therefore starts on every pull request and decides inside, in its `changes`
job, whether the window tests apply; a job skipped by its `if:` condition
counts as passed.

`tests/test_the_gate_can_still_say_it_passed.py` holds the list above against
the workflows: renaming a job renames its check and silently drops the
requirement, and that test fails first. Changing the list means changing the
ruleset and that test together.

## Emergency procedure (break-glass)

For situations in which the gate itself cannot answer, and waiting would do
more harm than merging untested:

- GitHub Actions is down or a required workflow cannot start;
- a required check is broken for reasons unrelated to any change, and the fix
  for it cannot pass the gate it repairs;
- a security fix must reach users before a full run could finish.

Not for: a slow run, a flaky test, or a deadline. Those wait, or the flaky
test is fixed first.

Only a repository admin can do this. The steps, in order:

1. Write down why: an issue, or the pull request's description, naming which
   of the cases above applies.
2. Disable the checks ruleset:

       gh api -X PUT repos/ComPlat/DELFIN/rulesets/25188885 -f enforcement=disabled

3. Merge the one change that needed it. Nothing else.
4. Re-enable it immediately:

       gh api -X PUT repos/ComPlat/DELFIN/rulesets/25188885 -f enforcement=active

5. Confirm that it is active again: `gh api repos/ComPlat/DELFIN/rules/branches/main`
   lists `required_status_checks` and `code_scanning`.
6. Run CI on the merge commit (`gh workflow run ci.yml --ref main`) and
   record its result next to the reason from step 1.

Every change to a ruleset is recorded in the repository's ruleset history and
the organisation audit log, so the procedure leaves a trace whether or not
the steps are followed. `main-hard-limits` is not part of it: a force-push or
deletion of `main` is never the remedy.
