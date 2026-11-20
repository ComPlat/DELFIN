# Kernel and mount limits of the host the suite ran on

Measured at commits `28a2002b0` (git grep thread ceiling), `78d5fbb06`
(sandbox gate) and `c41ae4a5f` (`/var/tmp`). Baseline for every number:
`origin/main` at `7abca4483` in its own worktree, `delfin.__file__` printed
per run.

All four limits below are properties of **this workstation**. A GitHub
runner has none of them, which is why the tests that met them passed on CI
and failed here. None of them is a bug in the product; each is a test (or a
test gate) reading the host where it meant to read the code.

## 1. The process is already inside a seccomp filter

```
$ grep -E '^(Seccomp|Seccomp_filters|NoNewPrivs):' /proc/self/status
Seccomp:         2
Seccomp_filters: 1
NoNewPrivs:      1
```

Consequence, measured in a forked child so the parent was unchanged:

```
install with SECCOMP_FILTER_FLAG_NEW_LISTENER ->
    OSError: [Errno 16] seccomp filter: Device or resource busy
install the SAME program with flags 0 -> succeeded
```

So a nested *listener* filter is refused while a plain one is accepted. The
product builds a listener, so every caged command was refused:

```
refused: sandbox could not be applied ([Errno 16] seccomp filter:
Device or resource busy); the command was not run.
refused: the cage supervisor could not take over (the child reported
no listener number ('')); the command was not run.
```

Reproduce with `scratch/probe_seccomp.py` (prints `Seccomp`,
`Seccomp_filters`, `filter_supported`, and the install result).

## 2. bwrap cannot unshare the network

```
$ bwrap --dev-bind / / --proc /proc --unshare-pid --new-session \
        --die-with-parent true
rc: 0
$ bwrap --dev-bind / / --unshare-net true
rc: 1
stderr: bwrap: loopback: Failed RTM_NEWADDR: Permission denied
```

So `_bwrap_functional()`, which probes the first form, answered **True**
while the wrap the tests actually build could not run. That is the whole of
the third cause: the gate asked a narrower question than the tests needed.

## 3. git grep cannot create one thread per CPU

384 cores, load average 192:

```
$ python3 -c "import os; print(os.cpu_count(), os.getloadavg())"
384 (192.56, 191.31, 191.66)

$ git -C <repo> grep -q -E '(def|class) present\b'      # tracked file holds it
rc=128  fatal: grep: failed to create thread: Resource temporarily unavailable
$ git -C <repo> -c grep.threads=1 grep -q -E ...        rc=0
$ git -C <repo> -c grep.threads=4 grep -q -E ...        rc=0
$ git -C <repo> -c grep.threads=16 grep -q -E ...       rc=0
```

Every `git grep` behind `safe_git` failed here. `knowledge_curation._git_grep`
maps a non-0/1 exit to `None` ("cannot tell"), so the knowledge statistics
silently left broken function references out instead of failing —
`test_knowledge_curation` was the only test that noticed.

## 4. /var/tmp is mounted read-only

```
$ python3 -c "import os; print(os.path.isdir('/var/tmp'), os.access('/var/tmp', os.W_OK))"
True False
$ python3 -c "import tempfile; tempfile.mkdtemp(dir='/var/tmp')"
OSError: [Errno 30] Read-only file system: '/var/tmp/probe-...'
```

The hide-path fixture hard-coded `/var/tmp` for a lockbox that a cage must
not already mask (a lockbox under `/tmp` is invisible because every cage
mounts a fresh tmpfs there). Five tests came back as setup ERRORS.

## Numbers, before and after

| Check | Baseline `7abca4483` | After |
|---|---|---|
| `test_a_git_read_does_not_ask_for_a_thread_per_cpu.py` (new) | 2 failed, 1 passed | 3 passed |
| `test_knowledge_curation.py` | 1 failed | 24 passed |
| two caged-command files | 11 failed, 879 passed | 733 passed, 158 skipped |
| `test_a_run_can_hide_paths_from_its_cage.py` + `test_no_cage_reaches_anyones_pairing_door.py` | 2 failed, 6 passed, 7 errors | 6 passed, 8 skipped, 1 failed |
| every file naming the sandbox gate (11 files) | 2 failed, 6 passed, 7 errors (subset) | 817 passed, 179 skipped, 1 failed, 5 errors |
| `DELFIN_EXPECT_ISOLATION=1` on the caged files | — | 4 failed, 2 passed (the promise is still enforced) |

The skip lines name the refusing layer, measured:

```
a caged command is refused here (([Errno 16] seccomp filter: Device or
resource busy); the command was not run.); set DELFIN_EXPECT_ISOLATION=1 to
make this a failure
no writable directory outside /tmp on this host
```

## The product change that was tried and reverted

Making `socket_guard.available()` answer "can a listener filter be
installed here" — the honest question — was implemented, measured, and
reverted. On this host the answer is False, which routes a caged command
into the `--unshare-net` path instead; that layer fails here too (limit 2
above), so the failures rose from **11 to 116**. The product's per-layer
fallbacks are correct; only the TEST gate was asking the wrong question.
That measurement is why the fix landed in `tests/conftest.py`.

## Not verifiable from where I stood

- **CI was not observed.** `git fetch` cannot reach the remote here. That a
  runner is neither inside a seccomp filter, nor at load 192, nor on a
  read-only `/var/tmp`, is the premise of this whole task and is not
  something I could confirm by running CI.
- **The skips are skips.** On this host the confinement tests make no
  assertion; they report that the host cannot host the measurement. The
  ordinary-host arm is what CI runs, and I did not run it.
- **The remaining 1 failure and 5 errors** in
  `test_no_cage_reaches_anyones_pairing_door.py` (`KeyError: 'link'`) and in
  the gate-consuming set reproduce identically at the baseline. I did not
  find their cause; they are recorded in `scratch/PR.md` as open, not
  claimed as fixed.
