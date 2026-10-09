---
name: push
description: Check, commit and push this branch to the remote
domains: code
---
# Push this branch

The user pressed "Push": please push the current branch to origin. This is
their request for one push of this branch, now.

## Steps
1. `git status` and `git branch --show-current`.
2. On the default branch (main/master) as a contributor (agent.git_role is
   not `maintainer`): do not push it. Create a branch first, as in the
   `new-branch` skill, and push that.
3. Uncommitted changes of THIS session: commit them with a message that
   says what changed and why. Changes you did not make are not yours to
   commit: list them and ask.
4. Before pushing, run the checks a reviewer would: the tests that cover
   the changed files, and `ruff check` when the project uses it. A failing
   check stops the push: report it and offer the fix.
5. Push: `git push -u origin <branch>`. A rejected push (the remote moved
   on) is not forced: bring the remote in with the `merge` steps, then
   push again.
6. After the push, watch its CI where that is available and call it
   pending until the result arrives.

## Report
Branch, commits pushed, the check results, and the link for a pull
request when the remote prints one.
