# test-evidence

Measurements that CI cannot reproduce, kept beside the code they are about.

## What belongs here

A measurement whose apparatus is not available to CI, and which is therefore
the only record that something was actually driven rather than reasoned
about:

- a run against a real service, cluster or device — the login node, a SLURM
  batch, the lab Pi, a live model endpoint;
- a run in a real browser, a real pty, a real cage;
- a before/after number taken on a host CI does not have (a machine with a
  QM program installed, a loaded workstation, a different glibc);
- a control run at a baseline, where the point is the *difference* between
  two trees rather than a green suite.

## What does NOT belong here

**A CI result.** The suite is its own evidence, it runs on every commit, and
a copy of its output in a file is a copy that goes stale. Writing one here
would turn the directory into a changelog, which is how this kind of
directory dies.

So: if CI could have produced the number, it does not go in a file. If CI
could not, the file is the only place the number will ever exist.

## The form

One file per change, named after its branch or its pull request:
`test-evidence/<slug>.md`. It states, in this order:

1. **the commit it was measured at** — a measurement without a commit cannot
   be attributed, and `origin/main` moves;
2. **what was driven**, concretely enough to repeat: the command, the host
   kind, the apparatus;
3. **the numbers**, before and after, with units;
4. **what was NOT verifiable from there**, named. This is the part that
   makes the file worth reading: a measurement that claims no limits is
   claiming the wrong thing.

Prose without a number is not evidence. `tests/test_the_evidence_is_a_measurement.py`
holds that: a file here has to carry a commit and at least one figure.

## Why in the repository

Until 2026-11-20 these measurements lived in two places, both outside the
code they prove: the body of a pull request, and a coordination file in
another repository whose log line had grown past two hundred entries. Both
answer "was this driven?" only for as long as someone remembers where to
look. A file beside the change answers it in a `git log`, and travels with
a clone.

## Files

Kept as an index because the directory now holds more than one file,
and the one-line summaries say which measurement a reader wants
without opening three files. Added by the branch that brought the
first three.

- [`agent-session-variables.md`](agent-session-variables.md) — the three
  ambient variables a DELFIN agent session exports, and the tests that
  measured them instead of the module.
- [`sandbox-and-mount-limits.md`](sandbox-and-mount-limits.md) — the seccomp
  listener refusal, the bwrap network-namespace failure, the git grep thread
  ceiling, and `/var/tmp` being read-only.
- [`peer-id-and-socket-path.md`](peer-id-and-socket-path.md) —
  `SO_PEERCRED` answering with the sandbox supervisor's pid instead of the
  peer's, and a `TMPDIR` deep enough to overflow `sockaddr_un`, with the
  before/after over the nine pairing-test files.