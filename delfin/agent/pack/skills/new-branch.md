---
name: new-branch
description: Start a new branch (from main, or from the branch given)
argument-hint: [new-name] [from <branch>]
domains: code
---
# Start a new branch

The user pressed "New branch". The arguments may hold the new name and,
after `from`, the branch to start from; without one it starts from the
default branch.

## Steps
1. `git status`. Uncommitted changes come along to the new branch; say so,
   and say which files.
2. No name given: derive one from what this session is working on, in the
   form `<topic-in-a-few-words>` (lowercase, hyphens), and use it. Only if
   there is nothing to derive it from, ask the user for a name.
3. `git fetch origin`, then `git switch -c <name> <start>` -- the start is
   the `from` branch (its `origin/` version when it exists there), else the
   default branch from `git symbolic-ref --short refs/remotes/origin/HEAD`.
   If the name exists already, say so and ask whether to switch to it.
4. Nothing is sent to the remote by this step. The branch stays local
   until the user asks for it to be published.

## Report
The branch name, what it starts from (`git log --oneline -1`), and which
uncommitted changes came along.
