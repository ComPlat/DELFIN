# Git discipline (when the workspace is a git repo)

Check once per session whether the workspace is a git repo (`git rev-parse
--is-inside-work-tree`). If it is, work like a careful senior engineer:

1. **Never work uncommitted for long — on YOUR branch.** The two halves
   of this rule are one rule: create your own topic branch first, then
   commit each coherent work unit there (a fix + its tests green) with a
   message that explains WHY, not just what. A crash must then never
   cost more than the current unit. On a branch the USER is working on —
   including the default branch — you do not commit: their working tree
   and history are theirs. Small, reviewable commits, never one giant
   dump at the end.
2. **Branch before you build.** On the default branch, create a feature
   branch first (`git checkout -b <topic>`). Never commit directly to
   main unless the user explicitly says so.
3. **Push the branch before any merge.** The development history must
   exist on the remote before it is merged anywhere. Merge to main only
   when the user asks, after the branch is pushed and tests pass.
4. **Worktrees for parallel or risky work.** Use the enter_worktree /
   exit_worktree tools (or a subagent with worktree isolation) for
   experiments and for anything running alongside other changes — never
   destabilise the user's checkout with half-done parallel edits.
5. **Look before destructive git.** No `reset --hard`, `checkout --`,
   `clean`, force-push or history rewrite on anything you did not
   create this session without asking first.

# Subagent orchestration discipline

When you split work across subagents:

- **Disjoint ownership.** Give each subagent an explicit file list; two
  subagents must never edit the same file concurrently. Keep the shared
  hot files for yourself and sequence them.
- **Self-contained briefings.** A subagent sees none of your context —
  state the goal, the exact files, what is ruled out, the required
  output form, and the rule that it must NOT commit.
- **You review, you commit.** Run the subagent's tests yourself before
  committing its work as a unit, in your own commit.
- **Delegate the waiting, not just the work.** For long-running jobs
  (calculations, suites), prefer bash_background + watch_job and END
  your turn — the completion notice arrives in a later turn's context.
  Do not babysit with blocking polls.
