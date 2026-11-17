---
name: pr
description: Push this branch and open a pull request (into main, or the branch given)
argument-hint: [base-branch]
domains: code
---
# Open a pull request

The user pressed "PR": push this branch and open a pull request into the
base branch from the arguments, or the default branch without one. This is
their request for that push and that pull request.

## Steps
0. Branch hygiene before anything is pushed (measured failures of
   wave 26): the branch carries ONLY the files its commits name --
   probe/scratch files (`test_zz_*`, notes) are never committed; a
   review of another branch happens on a scratch branch, never by
   merging the reviewed branch into this one; `main` (the base) is
   left where it was. `git log --oneline origin/<base>..HEAD` and
   `git diff --name-only origin/<base>..HEAD` show both: fix the
   branch BEFORE the PR, not in review.
1. GitHub CLI: `gh --version` and `gh auth status`. Missing or not logged
   in: stop. Say what is missing and how to fix it -- `gh` can be installed
   through DELFIN's installer when it offers it, and the login is the
   user's to do: they type `! gh auth login`. Do not work around it.
2. Do the `push` steps first (checks, commit, `git push -u origin <branch>`).
   Never open a pull request from the default branch itself.
3. Without write access to the repository (its git role is not
   `maintainer` and the push is refused for permissions), use a fork:
   `gh repo fork --remote --remote-name fork`, push the branch to `fork`,
   and open the pull request from there.
4. Build the body with `delfin.agent.pr_helper.body(<workspace>,
   "origin/<base>", session_id=<your session id from task_list>)`: goal and
   open points from this session's tasks, what/why/how tested from the
   branch's commits, account names and secrets redacted. Write the title
   yourself. Then `gh pr create --base <base> --title ... --body ...` (base
   without the `origin/` prefix).
5. An open pull request for this branch already exists
   (`gh pr view --json url`): do not open a second one; report its link.
6. Merging is a maintainer's step, and only on the user's request: wait
   until every check on the PR's current head is green (`gh pr view <n>
   --json headRefOid,statusCheckRollup`), then `gh pr merge <n>
   --match-head-commit <that sha>`. Never `--admin` or `--auto`. A
   contributor, or a user who asked for the PR route, does not merge from
   here: give the link.

## Report
The pull request link, what it contains, and the check results. When the
user asked for it, say "watching until green" and continue with the
`pr-until-green` skill: it waits for the checks (no tight polling),
repairs reds on this branch within its limits, and stops and asks the
person at its limits. It never merges.
