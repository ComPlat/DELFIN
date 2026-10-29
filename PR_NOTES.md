# PR notes — persistent per-session scratch directory across the cage layers

Task: inside every cage layer (bwrap in `delfin/agent/sandbox.py`, the
api_client bwrap bash cage, Landlock, Seatbelt, and no-isolation), /tmp is
fresh per command and invisible to the file tools, so scratch files written
by one command are gone for the next and agents lose work.

## Design

One helper, one root, four consumers.

1. **`_session_scratch_dir()`** (new, `delfin/agent/api_client.py`, beside
   `_private_tmp_dir` at ~line 28509): a `scratch/` directory under the
   process's own `_private_tmp_dir()` (`tempfile.mkdtemp(prefix=
   "delfin-cmd-tmp-")`, mode 0700, `atexit` rmtree). That single choice
   gives, without new bookkeeping:
   - **0700** (inherited from the parent; `mkdir(mode=0o700)` on the leaf);
   - **never shared between sessions** — `_private_tmp_dir` is per *process*,
     and one session process = one session; a second session process gets its
     own mkdtemp tree (test spawns a fresh interpreter and asserts the dirs
     differ);
   - **pruned when the session ends** — the existing atexit rmtree of the
     parent takes the scratch subtree with it.
2. **File tools** (`api_client.py`): `KitToolPermissions.all_workspace_roots`
   and `all_readable_roots` append the scratch path (also when
   `scope_locked` — it is the session's own private state, not a user grant,
   so this does not widen the lock's promise to anything outside the
   session's own tree). `find_root_for` / `find_readable_root_for` then make
   `read_file` / `write_file` / `list_files` accept absolute paths under the
   scratch dir. Path deny globs still apply (relative under the root), which
   is correct: an agent may not name a file `secrets.pem` there either.
3. **Cages — bound read-write, TMPDIR/TMP/TEMP exported, /tmp untouched:**
   - **Landlock** (`_landlock_argv` ~28887 and `sandbox._unsandboxed_argv`
     ~959): add `--write <scratch>` and pass `--tmpdir <scratch>`
     (`landlock_exec` already sets TMPDIR/TMP/TEMP from `--tmpdir`,
     landlock_exec.py:600-603, so nothing more is needed there). Note
     `_landlock_argv` currently passes `--tmpdir _private_tmp_dir()` — the
     scratch dir is *under* that, so today Landlock commands already share
     the private tmp, but files there are not tool-readable and the helper
     prunes nothing per session; the explicit scratch root makes it visible
     to the tools.
   - **Seatbelt** (`_seatbelt_argv` ~28851): `prof` write_roots and
     allow_sockets get the scratch dir, and `env_prefix` exports
     `TMPDIR=<scratch>` (currently `TMPDIR={tmp}` with tmp = private tmp).
   - **bwrap, api_client bash cage** (~30146): replace the three
     `--setenv TMPDIR/TMP/TEMP /tmp` with the scratch path and add
     `--bind-try <scratch> <scratch>` to `_write_binds` (bind-try, because
     the cage must not refuse to start for a missing dir; the helper creates
     it, but a cage can be built before the first call).
   - **bwrap, `delfin/agent/sandbox.py` `_bwrap_argv`** (~1053): after the
     repo bind, `--bind <scratch> <scratch>` (source exists — the helper
     created it on the host — so a plain bind is safe here; it already
     `--tmpfs /tmp` stays). The env: `run_agent_command` passes
     `_scrubbed_bash_env()`; set `TMPDIR/TMP/TEMP` in that dict to the
     scratch path (only when the cage is actually applied — same condition
     as the argv builds; with `mode="off"` the process env of the user's
     shell is not ours to rewrite, so there set them too? Decision: set them
     unconditionally in the env handed to the *command*, in every mode: the
     env dict is a fresh copy per run, it does not touch the user's shell,
     and pointing scratch at a writable 0700 dir is safe even uncaged).
   - **firejail** (`_firejail_argv` ~1110): `--private-tmp` makes /tmp
     fresh, but the scratch dir is NOT bound into the cage — there is no
     naive way to make the two compose (whitelist/bind attempt: scratch
     lives under a mkdtemp path that can contain glob chars from TMPDIR,
     and --whitelist is per-path and repo only). **Decision (implemented):
     firejail keeps `TMPDIR/TMP/TEMP = /tmp`** — firejail's own private,
     writable tmp — so the behaviour is exactly what predates the scratch
     dir and mktemp keeps working on a firejail host. The per-layer env
     change is done centrally in `sandbox._command_env(cfg.mode)`, which
     applies `_scrubbed_bash_env()` (TMPDIR=scratch) and resets the three
     variables to /tmp for mode "firejail" only. Pinned by
     `test_firejail_command_env_keeps_tmpdir_at_firejail_tmp`.
   - **no-isolation / allowlist fallback** (`run_agent_command` env +
     `_network_only` path): the env copy gets TMPDIR/TMP/TEMP = scratch; no
     filesystem change needed (nothing is confined anyway).
   - **socket-guard Landlock inside bwrap** (~30294-30327): add the scratch
     dir to `fs_part` `--write` list and to the `--allow-socket` list next
     to `/tmp`.
4. **Package-cache env**: `_package_cache_dirs()` builds
   `{k}={tmp}/{sub}` from the TMPDIR value; point the same subdirs under
   the scratch dir (Seatbelt env_prefix and bwrap `--setenv`), so caches
   also survive the session — same rationale, same consumers.

## Security reasoning

- **Least widening.** The only new writable path is a 0700 directory the
  session process created itself, inside its own private tmp tree. Another
  session cannot reach it (different mkdtemp, 0700 of the same uid); another
  *user* never could. Nothing under the workspace, home, archive or calc
  becomes writable through this.
- **No gate bypass.** The file-tool change adds the scratch root to the
  same root lists every other writable root uses, so the deny-glob screen,
  the change journal, the write gate's confirm logic and the bash read scan
  all apply to it exactly as to the workspace. It is deliberately *not*
  added to `add_extra_dir`/grants: no user approval, no scope widening —
  `scope_locked` sessions include it because the lock's promise is "cannot
  leave the workspace", not "cannot leave its own process".
- **Secret hygiene unchanged.** The scratch dir is not in `_HOME_SECRET_DIRS`,
  not hidden, and not on the path deny list. It holds no credentials —
  `_scrubbed_bash_env` keeps provider keys out of commands that write there.
  Landlock `_landlock_hide_paths` must NOT hide it: it is under
  `_private_tmp_dir()`, which the hide-path computation already treats as a
  cage-private path, not a secret location.
- **Lifecycle.** Pruning = process exit (existing atexit rmtree). No
  cross-session persistence is intended; a session that ends takes its
  scratch with it, which is also why this is not a place to stash anything
  that must outlive the session (the workspace remains the durable store).
- **Known limit: the scratch dir is per PROCESS.** `_SESSION_SCRATCH` is a
  module global, so its scope is the *process*, not the logical session.
  Consequence, stated rather than hidden: sub-agents (which share the
  parent process) and the scheduler daemon (several runs in one process)
  all see the SAME scratch directory. They are the same principal doing
  stateful work and cannot harm each other any more than they already can
  — each already has full write access to the whole workspace — so this is
  a convenience boundary, not a security one. Two separate *processes*
  still get two different 0700 dirs (pinned by the fresh-interpreter
  test), which is the actual isolation guarantee. The same-dir test pins
  the sharing so it is a deliberate, reviewed property.
- **Residual risk, stated.** (a) Inside the session, scratch is freely
  readable/writable by commands and tools — that is its purpose; it is not a
  place for secret material and nothing enforces that beyond the existing
  env scrubbing. (b) On hosts where the kernel offers neither bwrap, nor
  firejail, nor Landlock, the "cage" was already only the allow-list; the
  scratch dir changes nothing there. (c) `TMPDIR` now points at a path
  outside /tmp: tools that hard-code /tmp assumptions are unaffected
  (they still get a fresh /tmp), but tools that *derive* state from TMPDIR
  now persist — intended. (d) **firejail is the one layer that does NOT
  bind the scratch dir**: it runs `--private-tmp` and does not whitelist
  the scratch path, so its command env keeps `TMPDIR/TMP/TEMP = /tmp`
  (firejail's private writable tmp), restoring the pre-scratch behaviour.
  Every other layer gets `TMPDIR = scratch`. The argv/env test pins this
  so it cannot silently regress.

## Test (red first, then green)

`tests/test_the_session_scratch_dir_survives_the_command.py` — written
against the UNFIXED code and confirmed failing (7/7 AttributeError). It
covers, per layer, on this host (bwrap and Landlock both functional here,
probed):

- scratch exists, 0700, under the process private tmp; a fresh child
  process gets a different one;
- TMPDIR/TMP/TEMP inside a caged command point at the scratch dir;
- a file written by one caged command is visible to the next command and
  to `read_file`;
- the write gate admits a `write_file` into the scratch dir;
- `_landlock_argv` carries `--tmpdir <scratch>` and `--write <scratch>`
  (exactly one `--tmpdir`, and it is the scratch dir);
- `_landlock_argv` falls back to `--tmpdir <private tmp>` when
  `_session_scratch_dir()` raises (still exactly one `--tmpdir`);
- two `KitToolPermissions` in the same process resolve to the same scratch
  dir (documents the per-process sharing);
- (`test_agent_sandbox.py`) `_command_env("firejail")` keeps
  TMPDIR/TMP/TEMP at /tmp while `bwrap`/`landlock`/`off` get the scratch;
- /tmp stays fresh: a file written to /tmp by one command is gone for the
  next.

## Status

Implemented and committed on `swarm/wave4/t04`, each step test-first:

- `51497b765` tests: session scratch dir survives the command (red).
- `1e050f098` agent: per-session scratch dir bound into every cage, TMPDIR
  points at it.
- `8783bfc3d` tests: fixed three test bugs (setenv parse, child interpreter
  + PYTHONPATH, real-cage gate); Landlock one-tmpdir assert left red.
- `8efe1e4fc` fix(landlock): pass exactly one `--tmpdir` so TMPDIR is the
  scratch dir (keeps private tmp as a separate `--write` root).
- `4d021a013` + `ca6d62b74` firejail: env keeps TMPDIR at /tmp (test, then
  `sandbox._command_env` fix).
- `5e0525fe4` + `98568a559` landlock fallback: `--tmpdir _private_tmp_dir()`
  restored on the scratch-unavailable path (test, then fix).
- `cf3539ba4` tests: two KitToolPermissions in one process share the scratch
  dir (pins the per-process sharing).

Remaining edits were rejected or deferred by review and are captured in the
sections above: there was no approval-window stall on this run — the
api_client and sandbox changes were edited directly in the worktree.

Not done / known from review:
- firejail is not exercised with a real `firejail` binary here (not
  installed); the argv/`_command_env` test is the oracle.
- Three real-cage achieved-run smoke tests fail on this host for host
  reasons (nested-container seccomp refusal and a bwrap `lo` lookup), also
  on `origin/main` — not regressions from this work.
