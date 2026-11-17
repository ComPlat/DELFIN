# Swarm practice

How a swarm of agent workers runs a shared project well. This is the
practice the operators learned from running waves of sessions; it applies
to any team of cooperating agents.

## Goals first, evidence before claims

- A task states its goal and its success criterion up front. Unfinished
  work is reported as unfinished; a number is only reported with its
  source; a cause is only claimed if it was tested, otherwise it is
  labelled a hypothesis.
- A failing test or run is reproduced and understood before anything is
  changed. A change lands only when a red test (written first) turns
  green, or a measured number improves while nothing else gets worse.

## ONE focus per session

- A session works one task at a time. It does not drift to a different
  problem or silently redefine an easier one. When one step is done, it
  names what it waits for and from whom, and notes the wait where the
  team can see it.

## Scope is enforced, not promised

- Each worker writes only the files it owns; the scope is enforced by the
  gate, not by intention. A refusal names the route, not a spelling: the
  next step is the allowed route or a question, never the same act
  through another tool. Nothing is deleted that another worker may have
  made; a file that looks wrong is read before it is changed.

## Never push or merge as a worker

- A worker never pushes or merges by itself. When a body of work is
  approved, the worker names the head for the operator (or person) to
  push. Getting a change to main is a human action, not an agent action.

## Red first, then the smallest change

- Every change starts as a test that fails on the unchanged code; the
  control is written into the commit message. The fix is the smallest
  change that makes it green. A test that pins old behaviour is not bent:
  it is rewritten along its intent and the old decision is stated.
- A failing test is run on the base state before a change is blamed:
  a passing base and a failing change is a regression of the change; a
  failing base is a pre-existing failure and is reported, not fixed
  silently.

## Small, reviewable commits

- Commit as you go, one finished thought per commit, each with a message
  that says WHY, not just what. Never amend or rewrite a commit someone
  may have merged; make a new commit.

## Reviewers write adversarial tests

- A reviewer reads every commit as it lands, questions the diff, and
  writes adversarial tests against the new code. A finding names a file
  and line and the evidence; the builder fixes it or refutes it with
  evidence, never by weakening a test. A reviewer writes adversarial
  tests; it does not commit product code.

## Verification runs the real path

- Tests run through the project's sanctioned test entry, through the real
  call site, on real-shaped records. A mock that delivers nothing is not
  a test. A metric decided by wall clock, a matcher that never matched a
  real command, a path filter that blocks a real request: all passed
  their own tests once and were wrong.
- Nobody approves on word alone. The quality lead runs the tests
  itself before marking a package green.

## Public text stays clean

- PR text, commit messages and docs carry no account names, user names,
  host names, or machine paths. What a reader cannot reproduce from the
  text is left out.

## Standards stay short

- A rule that is only a list of forbidden phrases is a symptom: ask what
  makes the bad input impossible instead (who may call, what shape a
  value may have, where it may come from). Templates and guides stay
  short so they are followed, not skimmed.

## Keep going without waiting for a nudge

- A waiting session names the message it waits for and from whom, and does
  not fabricate work to look busy. Only a new event (a commit, a ruling, a
  decision) justifies a fresh round of tools. A partner that is genuinely
  waiting on you is told so at once.
