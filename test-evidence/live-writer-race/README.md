# Evidence measured on this host

This file records numbers that CI cannot produce, because they come from
thread scheduling on one machine rather than from a test. Everything central
to the change is proven by the committed suite; this is the one measurement
that is not.

## What was driven

`append` writes a record's line and then its head file. A reader that arrives
between the two finds the log a few bytes longer than its head says. Whether
that gap is tampering or a write in flight is the question this change
answers, and the answer only shows up under a live writer.

`probe.py` in this directory runs one writer thread appending records in a
loop and one reader calling `audit_log.verify` in a loop for 3 s, then counts
the reads and the reports that came back not intact, split by problem kind.

## Commit measured

* **before**: `origin/main` at `5f1975075` (`14a08f5f9` is an ancestor; the
  branch was fetched immediately before the run, and
  `git rev-list --count origin/main..HEAD` = 0 in its own worktree).
* **after**: `63ea8ae5a` on `swarm/night5/n5i`, in the task worktree.

Both arms printed the loaded module's path and it was the worktree's own
copy, not an installed one:

```
before: /tmp/delfin-wt-1cba2275/delfin/agent/audit_log.py
after : <the agent's own checkout>/delfin/agent/audit_log.py
```

The second path is written as a placeholder: an absolute path here
encodes a home directory and the site's group-account name, which no
tracked file carries (tests/test_no_personal_data_is_committed.py
refused this file for exactly that). What the measurement needs is
that the two paths DIFFER -- the import came from the worktree, not
from the installed copy -- and that survives the redaction.

The probe is run as `PYTHONPATH=. python3 scratch/probe.py` **from the
worktree root**. That is not a detail: run as a bare script from inside a
subdirectory, Python puts the script's own directory first on `sys.path` and
`import delfin` then resolves to an installed copy instead of the worktree —
the first attempt at this measurement did exactly that and measured code that
was not the worktree's.

## Numbers

Six runs per arm, 3 s per run:

| Arm | Reads | Reported not intact | Rate | Kinds |
|---|---|---|---|---|
| before (`origin/main`) | 10 311 | 1 707 | **16.6 %** | `head` 1705, `unparseable` 6 |
| after (this branch) | 9 311 | 9 | **0.1 %** | `head` 9, `unparseable` 9 |

Per run, so the spread is visible rather than hidden in the total:

```
before: 355/2139  228/1375  258/1500  218/1544  249/1485  399/2268
after :   1/1487    2/1442    1/1720    1/1519    2/1564    2/1579
```

The `head` findings on the before arm are the defect: two readers in three
seconds hit a live log sixteen times in a hundred and were told their log had
been tampered with. On the after arm the finding is no longer the head at all
— it is the other half of the same race, a reader catching a **partial line**,
which is honestly `unparseable` and is not what `_just_grew` claims to fix.

The earlier figure in the branch history ("204 of 1666", then "0 of 1378")
came from single runs with a different harness. These six-run numbers replace
it; the rate is what is stable, not any one count.

## What was NOT verifiable from here

* **The residual `unparseable` rate is not fixed, only reduced and relabelled.**
  A writer's half-written line is still a finding a reader can see. Closing it
  needs the reader to distinguish "torn last line" from "edited line", which is
  a different change from the one made here, and I did not make it.
* **Two processes, not two threads.** The writer here is a thread in the same
  process, so it shares `_LOG_LOCK`. A log appended by a *different process*
  was not exercised; the lock is per-process, so that is a genuinely different
  race and this measurement says nothing about it.
* **One host, one filesystem.** The rate depends on how the filesystem orders
  a write and a metadata update. It is a property of this machine's storage,
  not a constant; another machine may show a different before-rate. The
  direction (16.6 % → 0.1 %) is what the change is responsible for.
* **Not a CI number, and not a substitute for the test.** The committed test
  `test_a_log_longer_than_its_head_because_a_write_was_in_flight` constructs
  the same state deterministically by rewinding the head file after the line
  lands, so CI checks the fix on every run without depending on a race.

The probe source is `probe.py` in this directory, byte-identical to
`scratch/race_measure.py` in the task worktree.
