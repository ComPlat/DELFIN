# The environment an agent session hands to the suite

Measured at commit `76c3d84fd` (the change this file is evidence for).
Baseline for every number below: `origin/main` at `7abca4483`, in its own
worktree, with `delfin.__file__` printed in each run to prove which tree was
loaded.

## What was driven

The suite was started from a shell that a DELFIN agent session had created,
so it inherited the variables that session exports to every child:

```
$ python3 -c "import os; print({k:v for k,v in os.environ.items() if k.startswith('DELFIN')})"
{'DELFIN_CAGE_HIDE': '/home/localuser/.delfin',   # the operator's home, written as a stand-in
 'DELFIN_ENDPOINT_PER_MODEL_LIMIT': '1',
 'DELFIN_NO_BROWSER': '1',
 'DELFIN_SWARM_RUN': 'night5',
 'DELFIN_AGENT_SESSION': '985466',
 'DELFIN_SESSION_SCRATCH': '/tmp/delfin-cmd-tmp-.../scratch-...',
 'DELFIN_AGENT_LABEL': 'night5/n5l via blablador'}
```

`DELFIN_AGENT_SESSION` is `approval_guard.ENV_MARKER`. It is read as "a
person is not the one acting" by `knowledge_share` (`_refuse_in_agent_session`),
`trust_store`, `file_confirm`, `messenger.pairing`, the agent service
supervisor, and `api_client._classifier_settings` (which ignores
`DELFIN_CLASSIFIER=off` when it is set).

## Numbers

Each command run serially, one file at a time:

```
python -m pytest <file> -q -m "not slow" --timeout=300 --timeout-method=signal
```

| File / test | With the ambient variables | With them removed |
|---|---|---|
| `test_zz_u12_second_opinion.py::test_the_switch_reads_the_setting_and_one_process_can_override` | FAILED — `assert {...} is None` | passed |
| `test_knowledge_agent_surface.py` | 1 failed, 19 passed | 20 passed |
| `test_sx2_session_scratch_export.py` | 1 failed | 2 passed |

The removal was done two ways, and both are in this evidence because they
are not the same measurement:

1. **By hand**, `env -u DELFIN_AGENT_SESSION python3 -m pytest <file>` — the
   one-test-at-a-time proof that the variable is the cause.
2. **By the change**, which pops both variables in `pytest_configure` and
   restores them in `pytest_unconfigure`.

After the change, the same three files: 81 passed across the seven files
this branch touches.

A first attempt used `monkeypatch.delenv` in an autouse fixture and was
rejected by measurement: the classifier test calls `monkeypatch.undo()`
(20 call sites in the suite do), which reverts the shared instance the
autouse fixture had used, putting the variable back mid-test. Measured:
the test still failed. It was also reported by the suite's own `DELFIN_*`
leak detector, which snapshots before autouse fixtures run.

## Not verifiable from where I stood

- **The CI arm of each pair was not run by me.** `git fetch` cannot reach
  the remote from this workspace (`ssh: Could not resolve hostname
  github.com`), so "CI is green on the same commit" is the repository's own
  claim about `origin/main`, not a run I watched. The claim rests on the
  mechanism being unambiguous: on a runner where the variable is absent the
  tests passed before this change and pass after it.
- **Which variable a given shell carries is a property of that shell.**
  Another host may export neither, one, or others. The tests added here
  assert the property that matters — a test body sees no ambient copy, and
  a test that means to be inside an agent session can still set one — and
  both halves were run here.
- **`DELFIN_SWARM_RUN`, `DELFIN_CAGE_HIDE` and the other DELFIN variables in
  the list above were left alone.** They were not shown to change a test's
  answer here, and removing variables without a failing test behind them is
  how a suite stops measuring the thing it was written for.
