---
name: git-errors
description: Git error text to cause and safe fix; DELFIN checkouts; PR route
domains: code
---
# Git errors and their safe fixes
> Read the error, name its cause, run the standard fix. Never force, never `reset --hard`, never delete a lock, worktree or branch to make an error go away.

Each section: the text git prints, the cause, the fix, what not to do, and
when to ask. A fix that touches another checkout, another session's work
or the remote is the user's decision: say what you would run and why.

## Branch held by another worktree
Text: `fatal: '<branch>' is already used by worktree at '<path>'` (older git: `is already checked out at`).
Cause: one branch name can be checked out in only one worktree; `<path>` holds it.
Fix: `git worktree list` shows who holds what. To be at that branch's commit here, detach: `git fetch origin`, then `git switch --detach origin/<branch>` -- same commit, no name taken, the holder keeps its label.
To have the name itself here, the holder must release it first (`git switch --detach` run in `<path>`): ask the user, that checkout may be another session's or the dashboard's.
Never: `--ignore-other-worktrees`, `--force`, `git worktree remove`, deleting `<path>`.

## Detached HEAD
Text: `You are in 'detached HEAD' state`, `HEAD detached at <rev>`.
Cause: HEAD names a commit, not a branch -- normal after checking out a tag, a remote branch or a commit, and on purpose after a detach.
Fix: reading and building need nothing. Back to a branch: `git switch <branch>`. Commits made here: name them first with `git branch <name>`, then switch. Lost ones are in `git reflog`.
Never: switch away from new commits without naming them.

## index.lock exists
Text: `Unable to create '<path>/index.lock': File exists`, `Another git process seems to be running`.
Cause: a git process is running in this repository, or one died and left its lock.
Fix: `pgrep -a git`. One running: wait, then retry once. None running: tell the user the lock is stale and ask them to remove it (`rm <path>/index.lock`); the gate does not let the agent touch files under `.git`.
Never: remove a lock while a git process runs; retry in a loop.

## Push rejected or branches diverged
Text: `! [rejected] ... (fetch first)`, `(non-fast-forward)`, `Updates were rejected`, `have diverged`, `Need to specify how to reconcile divergent branches`.
Cause: the remote has commits this branch does not.
Fix: `git fetch origin`; your own unpublished commits: `git pull --rebase`; a shared or already published branch: `git merge origin/<branch>`. Resolve, run the tests, then push again -- only a push the user asked for, as in "Contributing a change".
Never: `git push --force` (refused); `--force-with-lease` only on your own branch and with the user's yes.

## Untracked files would be overwritten
Text: `The following untracked working tree files would be overwritten by checkout` (or `by merge`).
Cause: files untracked here are tracked on the target, and git will not replace them.
Fix: look at each listed file. Yours: commit it on your branch, or `mv` it aside. Not yours: ask.
Never: `git clean -f`, `git checkout -f`.

## Local changes would be overwritten
Text: `Your local changes to the following files would be overwritten by checkout` (or `merge`, `pull`), `Please commit your changes or stash them`.
Cause: uncommitted edits in files the operation changes.
Fix: your own changes: commit them (a WIP commit on your branch). Changes you did not make: list them and ask. The agent never uses `git stash` here: the stash is shared by every checkout of the repository (refused).
Never: `git checkout -- <file>`, `git reset --hard`.

## Merge conflict
Text: `CONFLICT (content): Merge conflict in <file>`, `Automatic merge failed`, `<<<<<<<` in a file.
Fix: the `git-merge-conflicts` skill: `git diff --name-only --diff-filter=U`, read both sides, edit, `git diff --check`, `git add <file>`, continue; back out with `git merge --abort` or `git rebase --abort`.
Never: commit markers; take one side whole without saying so. Ask when a side's intent is unclear.

## Not a git repository
Text: `fatal: not a git repository (or any of the parent directories)`, `fatal: not a git repository: <path>/.git/worktrees/<name>`.
Cause: the wrong directory, or a worktree whose `.git` FILE points at a gitdir that was moved or pruned.
Fix: `pwd`; `cat .git` shows the gitdir it points at; `git worktree list` in the main checkout shows whether it is still registered. Re-adding or repairing it is the user's call.
Never: `git init` there (a new, unrelated repository).

## Shallow clone
Text: `shallow update not allowed`, `(grafted)` in `git log`, `unknown revision` for a commit that should exist.
Cause: the clone holds only recent history.
Fix: `git rev-parse --is-shallow-repository`; then `git fetch --unshallow` (or `--deepen=<n>`).
Never: conclude a commit or tag does not exist before checking this.

## LFS pointer files
Text: a file that reads `version https://git-lfs.github.com/spec/v1`, `git-lfs: command not found`, `should have been pointers`.
Cause: the file lives in Git LFS; the checkout has the pointer, not the content.
Fix: `git lfs pull` where git-lfs is installed; otherwise tell the user it is missing.
Never: analyse, edit or commit the pointer text as if it were the data.

## Permission denied (publickey)
Text: `Permission denied (publickey)`, `Could not read from remote repository`.
Cause: the remote does not accept this host's key -- the host or the account, not the repository.
Fix: `git remote -v` names the remote; tell the user the key must be added to the account, or that they push themselves.
Never: vary transports or proxies; generate, copy or read keys.

## DELFIN's own checkouts and servers
Use when a DELFIN path (`.delfin/worktrees/`), a writer lock, a dashboard port or a question about checkouts comes up.
Fix: first run `delfin-agent ops status` (read-only) instead of guessing, then follow the safe fix it names.
The layout: the main checkout; the dashboard's own worktree, which often holds a branch (so `main` may be taken -- see the first section); session worktrees under `.delfin/worktrees/<name>` on `session/...` branches, brought back with `worktree_merge`. Each dashboard (voila) server runs from one checkout on its own port; a kernel whose last window closed ends after a 90 s grace, so a reload within that keeps it. A session's writer lock means another live process has that conversation open: open it there.
Never: kill a server, kernel or process you did not start, or any other user's; remove a worktree, lock or branch because it looks left over -- report it and ask.

## Contributing a change
The git role is set per repository by the user (`/git role`; `delfin-agent ops status` shows it). A contributor, or anyone after the user asked for a pull request, never pushes the default branch.
1. `git switch -c <user>/<topic>`. In DELFIN's own checkout (it stays on main): `enter_worktree(base_ref="origin/main")`, work there, push with `git push -u origin HEAD:<user>/<topic>`, then `gh pr create --head <user>/<topic>`.
2. Commit each step: `git add <paths>`, `git commit -m "<what and why>"`.
3. Push the branch, when asked: `git push -u origin <branch>`.
4. `gh pr create --base main --fill`, and give the user the link.
5. CI: `gh pr view <n> --json headRefOid,statusCheckRollup`; pending until every check on that head is green.
6. Merge only when the user asks and the role is maintainer: `gh pr merge <n> --merge --match-head-commit <headRefOid>`. Otherwise the user merges.
Never: `--admin`, `--auto`, a merge pinned to an older head, a push to main as a contributor.
