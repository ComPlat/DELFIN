---
domains: code
---
# Git: branch, merge, resolve conflicts
> Get work onto the right branch, merge or rebase, and resolve conflicts without losing anyone's changes.

Use this when a merge, rebase, pull or cherry-pick is needed, or when git
reports a conflict.

## Before anything
1. `git status` and `git log --oneline -5`: know the branch, what is staged,
   and which changes are not yours.
2. Changes you did not make stay untouched: no `git checkout -- <file>`, no
   `git reset --hard`, and no `git stash` in a checkout another session works
   in (the prompt lists them). Work that has to be set aside goes into a WIP
   commit on your own branch.
3. `git fetch origin` before comparing with the remote.

## Getting onto the default branch
- Contributor (agent.git_role, the default): never push to main or master.
  Commit on a branch of your own (`git switch -c <user>/<topic>`), push it
  with `git push -u origin <branch>`, and give the user the pull-request link
  the push note names.
- Maintainer: rebase onto the default branch first
  (`git fetch origin && git rebase origin/main`), run the tests that cover the
  change, then push -- once per request of the user.
- After a push its CI is watched; call it pending until the result arrives.

## Resolving a conflict
1. `git diff --name-only --diff-filter=U` lists the conflicted files.
2. For each file read all three versions before choosing: `git show :1:<path>`
   (the common base), `:2:<path>` (ours), `:3:<path>` (theirs). Know what each
   change was for.
3. Edit the file so it keeps the intent of both sides. Take one side whole
   only when the other change is truly superseded, and say so in the report.
4. No markers may remain: `git diff --check`, and no `<<<<<<<` in the file.
5. Run the tests that cover the resolved files, `git add <path>` for each,
   then `git rebase --continue` or `git commit` for a merge.
6. When the intent of a side is unclear, stop and ask the user which behaviour
   they want, showing both versions.
7. To give up cleanly: `git merge --abort` or `git rebase --abort`.

## Worktrees
A session in a worktree (`.delfin/worktrees/...`, branch `session/...`) brings
its work back with `worktree_merge` (applied uncommitted, for review), or by
committing on its branch and merging that branch; the conflict steps above
apply to both.

Report what was merged, each conflict and how it was resolved, and the test
result.
