# A peer id and a socket path the kernel decides, not the code

Measured at commits `7a8e4df5d`, `93649c16d` and `dd8cd07b0` (the changes
this file is evidence for). Baseline for every number below: the branch
point `7abca4483`, in its own worktree at `scratch/baseline-wt` with
`PYTHONPATH` inside it; `delfin.__file__` was printed in the baseline run
and resolved to the worktree, proving which tree was loaded.

Both mechanisms below are properties of **this workstation**. A GitHub
runner has neither, which is why the tests met them here and passed on CI.
Neither is a bug in the product: the guards that read them are correct, and
what was wrong was that the tests took the host's answer for the code's.

## 1. `SO_PEERCRED` does not name the peer

```
$ python3 -c "import os; print(os.getpid())"
7
$ # a socket pair inside THIS one process; the accepted side asks the kernel
client side creds (pid,uid,gid): (7, 1001, 1001)
server side creds (pid,uid,gid): (2, 1001, 1001)
expected peer pid: 7
```

The accepted socket answered pid **2**, not 7. `/proc/2` is an ancestor
(the sandbox supervisor) and its `environ` cannot be read:

```
('from_an_agent', 2, '->', True, ['pid differs',
  "environ unreadable: [Errno 13] Permission denied: '/proc/2/environ'"])
```

`pairing._peer_is_agent` delegates to `Supervisor._from_an_agent`, which
answers True for a pid whose environment it cannot read — deliberately, so
that a process which hides itself is not treated as a window. A truthful
kernel names the peer, the peer is this process, and the guard answers
False.

Consequence, before the change: the door answered `{"ok": false}` and the
invitation read `refused / "the asking window is gone"`.

What proves it: with the kernel's answer supplied, the same test that
failed returns `{'state': 'offered', ...}` and the link reaches the door;
with the host's answer it still returns `refused`. Both measured on
`7a8e4df5d`.

## 2. `TMPDIR` is deep enough to overflow `sockaddr_un`

```
$ python3 -c "import os; print(len(os.environ['TMPDIR']), os.environ['TMPDIR'])"
64 /tmp/delfin-cmd-tmp-lujwhdkr/scratch-hxbyd5bw
```

pytest appends `pytest-<n>/<test name>/b.sock`, which reached **110 bytes**
against a cap near 107:

```
OSError: AF_UNIX path too long
```

With `TMPDIR=/tmp`, two of the three parametrised cases passed and one
still failed — for the peer-id reason above, not this one. That split is
what separated the two causes.

## Numbers, before and after

| Check | At `7abca4483` (baseline worktree) | After |
|---|---|---|
| nine pairing-test files, together | **11 failed**, 64 passed | **78 passed, 1 skipped** |
| `test_no_cage_reaches_anyones_pairing_door.py` | 2 failed, 5 passed | 7 passed, 1 skipped |
| `test_the_messenger_runs_without_root.py` | 3 failed, 11 passed | 14 passed |
| `test_a_phone_name_is_shown_clean.py` | 2 failed, 1 passed | 3 passed |
| `test_the_invitation_link_is_never_persisted.py` + `test_pairing_of_a_second_contact_needs_todays_evidence.py` | 4 failed | 15 passed |

The baseline figure is from the worktree, not from this checkout: running
both and comparing is the whole point, and 11 → 0 is the change.

## Not verifiable from where I stood

- **The CI arm.** `git fetch` cannot reach the remote from this workspace
  (`ssh: Could not resolve hostname github.com`). That a GitHub runner's
  kernel names the real peer, and that its `TMPDIR` is short enough, is the
  premise of this task and of the numbers above — it is not something I
  measured.
- **`origin/main` as it stands now.** It moved to `5f1975075` during this
  work and that commit is not an ancestor of this branch. The control run
  is against the branch point `7abca4483`. Whether `5f1975075` behaves the
  same here was not measured.
- **A GitHub runner's `SO_PEERCRED`.** I inferred it from the fact that
  these tests pass on CI. I did not run a runner.
