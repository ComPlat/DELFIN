# test-evidence

Evidence for measurements in this change that CI cannot reproduce.

A file here is for a number that came from **this host** — a read-only
mount, a process already inside a seccomp filter, a thread budget spent by
an enclosing session — and that a GitHub runner therefore cannot produce.
A CI result does not belong here: the suite is its own evidence.

Files in this directory state, for each measurement:

- the commit it was measured at;
- what was driven (the exact command);
- the numbers before and after;
- what was **not** verifiable from where the measurement was taken.

## Files

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
