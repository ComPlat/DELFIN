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
1. GitHub CLI: `gh --version` and `gh auth status`. Missing or not logged
   in: stop. Say what is missing and how to fix it -- `gh` can be installed
   through DELFIN's installer when it offers it, and the login is the
   user's to do: they type `! gh auth login`. Do not work around it.
2. Do the `push` steps first (checks, commit, `git push -u origin <branch>`).
   Never open a pull request from the default branch itself.
3. Without write access to the repository (agent.git_role is not
   `maintainer` and the push is refused for permissions), use a fork:
   `gh repo fork --remote --remote-name fork`, push the branch to `fork`,
   and open the pull request from there.
4. Write the title and body yourself from the commits
   (`git log --oneline origin/main..HEAD`): what changed, why, how it was
   tested. Then `gh pr create --base <base> --title ... --body ...` (base
   without the `origin/` prefix).
5. An open pull request for this branch already exists
   (`gh pr view --json url`): do not open a second one; report its link.

## Report
The pull request link, what it contains, and the check results.
