# x1 — every module compiles without a warning (integration of w46/x1-a1 and w46/x1-a2)

**Should this land? Yes — the guard and the one-character fix, together, and only
together. The guard is red on this branch until `scratch/api_client-escapes.patch` is
applied; that is the intended state, not a defect of the branch, because rule 13 puts
the fix for the protected `api_client.py` in a patch.**

Base: `origin/main` = `41fbf18e280eda03cb531500d0eb41207e8ebadc`. My workspace was exactly
that commit and clean, so it served as the control tree; the pristine control for the
long suite is a separate `enter_worktree` at the same commit.

## The finding, measured on the untouched tree

```
$ PYTHONPATH=. python3 -c "import py_compile; py_compile.compile('delfin/agent/api_client.py', cfile='/dev/null', doraise=True)"
  File "delfin/agent/api_client.py", line 5330
    """(start, end) of each shell word in cmd[start:end], offsets kept.
    ^^^
SyntaxError: invalid escape sequence '\ '
```

Every `.py` under `delfin/` compiled with `compile(src, path, "exec")` inside
`warnings.simplefilter("error")`: **1288 files, exactly 1 offender** —
`delfin/agent/api_client.py:5330`, message `invalid escape sequence '\ '`. The arms
both reported this and both are right; I re-derived it here rather than taking it.

The brief says line 5333. The compiler blames **5330**, the line the docstring *opens*
on, not the line holding the backslash (5333). A test hard-coding 5333 would have been
red for the wrong line.

## Per-item decisions

| item | decision | commit | reason (one line) |
|---|---|---|---|
| a2's guard `tests/test_no_module_carries_a_syntax_warning.py` (5 commits) | **take** | `e26e0d58b` + `41b66e65e` `14e2995e1` `07231aa1d` `15484ab4d` | reports **every** warning in a file; a1's reports only the first (measured below) |
| a1's guard `tests/test_the_package_compiles_with_no_syntax_warning.py` | **reject** | — | on a file with two bad escapes it names only the first; the package asks for *each* offender |
| the fix `scratch/api_client-escapes.patch` | **take** | `f47d672ed` | both arms produced a **byte-identical** patch; one character on line 5330 |
| a1's `scratch/make_api_client_escapes_patch.py` | **reject** | — | a generator for a one-character fix that is committed once; nothing re-derives it |
| a1's and a2's `scratch/PR.md` | **reject** | — | their reports; this file supersedes them |
| "fix every other offender" | **nothing to fix — measured** | — | 1 offender in 1288 package files; the repo-wide scan the arms ran agrees |

### Why a2's guard over a1's — the measurement that decided it

The package asks the test to "list **each** offender file:line". Given a file with two
undefined escapes, the two instruments disagree:

```
two-escape file -> a1 reports: ["two_bad.py:2: invalid escape sequence '\ '"]
two-escape file -> a2 reports: [(2, "invalid escape sequence '\ '"), (5, "invalid escape sequence '\ '")]
```

a1 runs one `simplefilter("error")` compile per file, which stops at the first warning;
a2 runs an `"always"` pass as well. Everything else a1 has, a2 also has: red on base
naming the offender, two controls that pass on both trees (an invalid escape is
reported; a raw docstring is not), no host lookup, a file-count floor so a moved
`delfin/` cannot pass vacuously.

a2 also carries a test a1 only argues in prose —
`test_the_shell_word_docstring_keeps_the_backslash_it_shows` asserts the docstring
renders **one** backslash and not two, which is the package's "the intended character
must stay" turned into a check rather than a comment.

### Cost of a2's completeness, measured rather than assumed

| pass | time | what it adds |
|---|---|---|
| `read_text` all 1288 files | 0.09 s | — |
| the `"always"` pass | 4.37 s | **every** warning with its line (this is a2's reporting) |
| a2's second `"error"` pass per file | 4.72 s | nothing measurable — see below |
| `_package_ast.source()` per file | 7.43 s | shares the parse with the other tree guards |
| a1's whole scan | 4.61 s | first warning only |
| a2's whole scan | 25.71 s | every warning, 5.6× a1 |

The redundant `"error"` pass never disagreed with the `"always"` pass: **0 of 1288
files**. So a2 buys its completeness with the `"always"` pass (which a1 does not have
for reporting) and pays ~5 s on top for a second compile it does not need. I kept it
anyway: the task specifies `simplefilter("error")`, the pass is what makes an offender
unmissable whatever filter the rest of the run installed, and 25.7 s once per suite is
not a cost worth a second rewrite of another agent's tested work. Recorded here so the
next change has the number.

## Red → green, on the CI interpreter (Python 3.11.16, pytest 9.1.1)

| tree | command | result |
|---|---|---|
| pristine `origin/main` | guard test staged in, `pytest … -q` | **1 failed, 8 passed** — the failure names `delfin/agent/api_client.py:5330: invalid escape sequence '\ '` |
| pristine `origin/main` + `git apply scratch/api_client-escapes.patch` | same | **9 passed in 22.02s** |
| my branch HEAD (patch **not** applied) | same | **1 failed, 8 passed** — same offender |
| `tests/test_syntax_integrity.py` on base | `pytest -q` | **1289 passed** — and its warning summary carries the same `api_client.py:5330` DeprecationWarning twice: the existing parse-based guard is green over the defect |

The control case that passes on both trees, as a test, is
`test_a_raw_string_with_the_same_backslash_is_no_offender` — it would fail if the
instrument flagged the fix the offenders need.

The `-W error` boundary, measured directly:

| tree | command | result |
|---|---|---|
| base | `python3 -W error -c "import delfin.agent.api_client"` | **`SyntaxError: invalid escape sequence '\ '`** |
| base + patch | same | **`IMPORT OK`** |
| base + patch | `py_compile` with `doraise=True` under `-W error` | **`COMPILED CLEAN`** |
| base + patch | `_word_spans.__doc__` | one backslash and one space present, **no** double backslash — the raw prefix changed nothing a reader or caller sees |

## Predictions vs measurements (predictions were written before the runs)

| # | prediction | measured | verdict |
|---|---|---|---|
| P1 | exactly one offender, `api_client.py:5330`, 1288 files | exactly that | **right** |
| P2 | both arms red on base; neither branch carries the fix | both guard tests 1 failed / 8–3 passed on base; both patches unapplied in their branches | **right** |
| P3 | `test_syntax_integrity.py` green on base while printing the warnings | 1289 passed, 2 DeprecationWarnings for the same line | **right** |
| P4 | pre-existing failures in the full suite | confirmed: failures and errors appear from the first percent and are unrelated to this package (my branch changes no module) | **right** |
| P5 | a2's guard is the stronger one | measured: names both offenders where a1 names one; also pins the docstring text | **right** |
| P6 | patch makes `-W error` compile clean and changes no code object; no real run can show a behaviour difference | `SyntaxError` → clean; docstring byte-identical | **right** |

## Real runs: requested, with predictions

`scratch/realrun.json` holds two probes. **Prediction, written before the runs:**

1. **`guard-exists-and-names-offender`** — prompt asks whether a test compiles every
   module under `delfin/` with warnings as errors, and to run it.
   *Base:* the model finds no such test and says so.
   *Branch HEAD:* it finds `tests/test_no_module_carries_a_syntax_warning.py`, runs it,
   and reports exactly `delfin/agent/api_client.py:5330`.
   **The pair differs**: the guard is the delivered capability.
2. **`compile-offenders`** — prompt asks to compile every `delfin/` module under
   `warnings.simplefilter("error")` and list offenders.
   **Prediction: identical on both** — `api_client.py:5330` twice. This probe is
   deliberately a *control*: it documents that the underlying defect is still present on
   branch HEAD, because the fix ships as a patch that a session may not apply. If this
   probe's answer differs between base and branch, something other than my change moved
   the file and the run must be re-read, not believed.

I do **not** predict any probe showing a difference in DELFIN's *behaviour*: the patch
is one character inside a docstring, the code object is unchanged (both arms measured
this; I re-measured the docstring), and a behavioural difference would mean the fix
changed something it must not.

## Rule 13 — why the fix is a patch

`delfin/agent/api_client.py` is the first entry of `_DEFAULT_PATH_PROTECTED_GLOBS`
(`api_client.py:9859-9860`). I attempted the direct edit; the tool refused:

> `'edit_file' targets the agent's own safety layer ('delfin/agent/api_client.py') —
> this requires explicit user confirmation but no confirm_callback is configured, so it
> is refused.`

I did not route around it with bash. The change is `scratch/api_client-escapes.patch`,
committed, one hunk, and **needs applying**:

```
git apply scratch/api_client-escapes.patch
```

## What I did NOT do, and why

- **Did not apply the patch to the branch.** That would be going round the self-mod
  guard with bash (rule 13). Consequence, stated plainly: **a landing of this branch
  alone turns the fast suite red** by exactly one test, which is the guard doing its
  job. Land the patch in the same landing.
- **Did not fix the full suite's pre-existing failures.** They are present on the
  pristine base at the same commit; this package changes no module, so they cannot be
  mine (rule 15).
- **Did not edit any other module.** The scan says there is nothing else to fix.
- **Did not push.** The remote is configured `file:///nonexistent/…` for pushes; the
  branch is local.
- **Did not read `~/.delfin`** (rule 18).

## Suites run, with numbers

| suite | tree | result |
|---|---|---|
| every-`delfin/`-module compile scan | base | 1288 files, **1 offender**: `api_client.py:5330` |
| a1 guard, staged on base | base | **1 failed, 3 passed** |
| a2 guard, staged on base | base | **1 failed, 8 passed** |
| a2 guard (final, 9 tests) | base + patch, py3.11.16 | **9 passed in 22.02 s** |
| a2 guard (final, 9 tests) | branch HEAD, py3.11.16 | **1 failed, 8 passed** (the patch is not applied) |
| `tests/test_syntax_integrity.py` | base | **1289 passed** |
| guard + syntax integrity together | **branch HEAD + patch** (the artifact as it lands) | **1298 passed in 33.75 s** |
| repo-wide scanners: shard split, public API, undefined names, licence guard | branch HEAD + patch | **93 passed in 25.61 s** |
| source-reading population, 368 files in 6 chunks (see below) | branch HEAD + patch **vs** pristine base | **137 vs 137 unique failures — set difference empty both ways** |
| full `pytest tests/ -m "not slow"` | branch HEAD | **did not finish — the 3600 s job timeout killed it at 38 %.** No number. The population above covers every file that reads source; the remainder is unrelated to a change that touches no module |

### The population comparison, in full

The population is every test file that names `api_client` **and** reads source text
(`read_text|getsource|_package_ast|open(`) — 368 files, the ones that could observe a
source-text change. Both trees ran the same 368 files, chunked six ways, on Python
3.11.16.

| chunk | base failures | artifact failures | same set? |
|---|---|---|---|
| aa | 22 | 22 | yes |
| ab | 25 | 25 | yes |
| ac | 27 | 27 | yes |
| ad | 21 (after dropping the one file base does not have) | 21 | yes |
| ae | 30 | 30 | yes |
| af | 12 | 12 | yes |
| **unique union** | **137** | **137** | **`comm` empty in both directions** |

Chunk `ad` first measured 21 vs 1; that was **my** artifact, not the change: `ad`
contained `tests/test_no_module_carries_a_syntax_warning.py`, which does not exist on
base, and pytest aborted the base run with `ERROR: file or directory not found` rather
than reporting failures. Re-run without that file, base reports the same 21. This is the
one place a naive before/after would have read as "the branch broke 20 tests".

The 137 are pre-existing on `origin/main` and unrelated: sandbox/git-config, keyring
routes, grep/cat predicates, managed-policy gates, landlock isolation. This package
changes no module, so it cannot have caused them, and the identical sets say so
directly rather than by argument.
