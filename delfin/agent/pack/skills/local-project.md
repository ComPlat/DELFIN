---
name: local-project
description: Worktrees, branches and parallel sub-agents in a workspace folder that is not a git repository yet
argument-hint: [what the parallel work is]
domains: code
---
# Work GitHub-style in a plain folder

The workspace is a plain folder (no `.git`), and the task needs branches,
a worktree, or several sub-agents working in parallel.

## Steps
1. Do not run `git init` yourself. Call `enter_worktree`, or `orchestrate`
   / `subagent` with `isolation: "worktree"`: the tool asks the user ONCE
   whether the folder becomes a LOCAL repository (git init, branch `main`,
   a `.gitignore` that keeps QM/ORCA output, model weights, archives,
   large files and secrets out, one first commit). No remote is added.
2. A no is final for this session: work in the folder directly, one
   writer at a time, and do not ask again by another route.
3. Several writers: one `orchestrate` call with `"merge": true`. Each
   writer gets its own worktree and branch; afterwards every kept worktree
   is merged back into the folder, uncommitted. Read `merges` in the
   result: `ok`, or `conflict` -- then the folder holds the earlier
   changes and the conflicting worktree is kept at `path`; show the user
   both versions, never overwrite one with the other.
4. Review with `git diff`, then commit with explicit paths.
5. Never add a remote and never push: publishing a local project is the
   user's step. ORCA output stays out of every commit (LICENSE_GUARD.md).

## Report
Whether the folder became a local project (its first commit), which
worktrees were merged, which conflicted and where they are kept.
