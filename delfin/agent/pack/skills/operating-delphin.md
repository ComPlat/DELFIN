---
name: operating
description: Work in DELFIN's repo and runtime the measured way — tests, gates, knowledge, memory, partners
domains: code
---
# Operating DELFIN the measured way

> Small steps with evidence. Read what the check actually says before you
> change anything; a claim you cannot cite is not a claim.

## Tests and the gate

- Run only the test files your change touches, through the project's gate,
  never the whole suite: a suite on a shared login node costs other people
  their turn, and one slow test hides the one that matters.
- A red control comes BEFORE the fix, and it runs against the COMMITTED
  tree, not only your working copy — a fix and a test staged separately
  can leave the commit green in your editor and red everywhere else.
- When a test fails, read the failure verbatim before changing anything.
  Fixing the assertion instead of the cause makes the next failure
  quieter, not rarer.

## Knowledge and memory

- Shipped lessons live in the knowledge pack; every entry carries an
  evidence pointer (a test node, a commit, a PR, a DOI). A lesson without
  evidence is a rumour — do not write one, do not rely on one.
- Memory levels, narrowest first: this session, this project, the shipped
  pack, the user's store. When a wider level repeats a narrower one, keep
  the narrower copy and drop the wider — a lesson recalled twice is not
  twice as true.

## Working with partners and reviews

- Messages to partners are asynchronous: one message per event, keep
  working while one waits. If a partner is silent, read their branch —
  `git log`/`git show` on their branch is reading, not intruding.
- A review finding is answered with a fix or with evidence, never by
  weakening the test that found it. A reviewer test that fails is a
  defect with an address.

## Cluster and runtime

- Long commands run in the background and report back; poll nothing in a
  loop. Calculations go through DELFIN's own submit path, never a raw
  scheduler call from the shell.
- Every number you report carries its unit, its method and its source. An
  estimate is labelled an estimate.

## When you must stop

- A refused action is not retried in another spelling. Read the reason;
  the reason names the route that is allowed, or the person to ask.
- If you cannot verify something, say so plainly and name what would
  decide it. "Unverified" written down is progress; a confident guess
  committed is a defect someone finds at 2 a.m.
