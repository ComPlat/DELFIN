### How the control will be run, and why that shape

Rule 4 of the landing brief asks for the full suite once, serially, and for
**every failure** to be run on `origin/main` too. Those are two different
runs:

1. **The bundle, whole.** `pytest tests/ -q -m "not slow" --timeout=300
   --timeout-method=signal -p no:randomly`, serially, in `zz-land`. Running;
   at 21 % after 51 min, so roughly four hours at the observed rate. The log
   is `scratch/lands/full-suite-bundle.log`, and its last line will be
   `PYTEST_EXIT=<n>`.
2. **The control, per failing node.** The whole-suite-across-both-trees form
   doubles a four-hour serial run to eight, and this host's background-job
   cap (8) is shared with other work. The control that answers the question
   is the failing node ids, run in `scratch/lands/wt-control` — a git
   worktree at `origin/main` = `41fbf18e2`, confirmed with
   `git -C scratch/lands/wt-control rev-parse --short HEAD`. A node that
   **fails on both** was already red; only a node that **passes on the
   control and fails on the bundle** is this bundle's regression. That
   separation is not a shortcut — it is the same comparison, run on exactly
   the set that needs it.

A control worktree is not a control *process*: `git worktree` isolates files,
not imports. Every control run must therefore run with its own `cwd` set to
`wt-control` and must print the loaded module's path, so that a number
attributed to `origin/main` is provably not the bundle's code:

```
$ cwd=scratch/lands/wt-control  python -c "import delfin; print(delfin.__file__)"
```

If that line does not name `wt-control`, the measurement is discarded rather
than reported. (This is not hypothetical: an earlier probe in this session
loaded the installed copy instead of the worktree's and was discarded.)

**What is already established about the pre-existing red block.** e1's
report measured 482 failures on untouched `origin/main` over the same 151
dashboard, service, cage and shell files. So a raw failure count on the
bundle carries no meaning on its own, and no such count is quoted in this
PR as if it did.
