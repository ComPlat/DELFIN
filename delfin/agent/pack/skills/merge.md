---
name: merge
description: Bring another branch (default: main) into this one, resolve conflicts, run the tests
argument-hint: [branch, e.g. origin/main]
domains: code
---
# Bring another branch into this one

The user pressed "Merge": they want the newest state of a branch in the
branch they are working on. The branch they chose is in the arguments;
without one, it is the default branch. When the arguments name something
that is not a branch, or it is unclear which one is meant, ask. Do not
publish anything; this request does not include sending changes to the
remote.

## Steps
1. `git status` and `git branch --show-current`. Changes you did not make
   stay untouched. If uncommitted changes of THIS session are in the way,
   commit them on this branch as `WIP: <what>` first and say so. If changes
   of someone else are in the way, stop and ask the user what to do.
2. `git fetch origin`. The source is the branch from the arguments
   (prefer its `origin/` version when it exists there), else the default
   branch: `git symbolic-ref --short refs/remotes/origin/HEAD`.
3. When the current branch IS the source's local counterpart (on main,
   bringing in origin/main): `git merge --ff-only <source>`. If that
   refuses, the local branch has commits of its own: stop and show them.
4. Otherwise: `git merge <source>`. Merge, do not rebase: the branch may
   already be on the remote, and rewriting it would break it for everyone
   else.
5. Conflicts: follow the conflict steps of the `git-merge-conflicts` skill.
   Resolve a conflict yourself only when both sides' intent is clear (an
   import added on both sides, the same fix twice, independent additions
   next to each other). When the intent of a side is unclear, stop and ask
   the user, showing both versions side by side. To give up cleanly:
   `git merge --abort`.
6. Run the tests that cover the files the merge changed.

## Report
What came in (`git log --oneline HEAD@{1}..HEAD`, shortened), each conflict
and how it was resolved, and the test result.
