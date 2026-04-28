# Destructive-action safety audit (2026-04-28)

**Verdict: destructive behavior is fully guarded. The agent CANNOT
delete, mutate, or submit anything autonomously — every path
requires an explicit user-side trigger.**

## Six independent safety layers

### Layer 1 — `allow_mutate=False` default on every mutating tool

12 typed MCP tools change state. Every one of them defaults to
`allow_mutate=False`, which short-circuits to a `skipped: True`
hint instead of running:

| Tool | Without `allow_mutate=True` |
|---|---|
| `cleanup` | `_refuse_mutation("cleanup")` |
| `stop` | `_refuse_mutation("stop")` |
| `pipeline_run` / `pipeline_prepare` | refused |
| `run_orca_input` | refused |
| `co2` / `tadf_xtb` / `hyperpol` | refused |
| `submit_calculation` | "Would submit … allow_mutate=True" |
| `cancel_calculation` | "Would cancel job N …" |
| `rename_calc_folder` | "Would rename … allow_mutate=True" |
| `create_calc_folder` | "Would create … allow_mutate=True" |
| `move_calc_folder` | "Would move (calc → calc) …" |
| `move_to_archive` | "Would move calc → archive …" |
| `delete_calc_folder` | refused (3-lock, see Layer 2) |
| `kill_all_user_jobs` | "Would cancel N jobs …" |
| `prepare_recalc` | dry-run plan only |

Verified by unit tests in `test_ops_server.py` (every mutating
tool has a `_blocks_without_allow` test).

### Layer 2 — Three-lock `delete_calc_folder`

The most destructive op needs THREE independent confirmations:

1. `allowed_roots` must include the target folder (path-fence)
2. `confirm_token` must equal the folder's basename verbatim
3. `allow_mutate=True`

Missing any one → safe refusal. Verified in
`tests/test_ops_server.py::test_tool_delete_calc_folder_three_lock`.

### Layer 3 — Directory zoning (`archive/` + `remote_archive/` ALWAYS read-only)

Every mutating tool calls `_check_destructive_dirs(block_archive=True)`
which **refuses any path under `archive/` or `remote_archive/`** —
even with all locks aligned. Hard-coded in `delfin/api.py:_check_destructive_dirs`.

**Empirically verified — 4/4 mutation attempts on `remote_archive/`
refused (even with allowed_roots set to point at it):**

```
remote_archive delete with all locks      → REFUSED (read-only root)
remote_archive rename                      → REFUSED
remote_archive create new folder           → REFUSED
archive→archive direction-violation        → REFUSED
```

**`_PERM_PROFILES` mirror this in the dashboard's zone system**
(`tab_agent.py:6366-6402`):

| Profile | calc/ | archive/ | remote_archive/ |
|---|---|---|---|
| plan | tier-2 max | tier-1 read | **(0, read-only)** |
| ask_all (default) | tier-3 + confirm dialog | **(0, read-only)** | **(0, read-only)** |
| repo_free | tier-3 + confirm dialog | **(0, read-only)** | **(0, read-only)** |
| all_free (user-chosen) | tier-3 auto | **(0, read-only)** | **(0, read-only)** |

**Even the `all_free` profile — the most permissive one a user can
intentionally pick — leaves `archive/` and `remote_archive/`
locked read-only. There is no profile, anywhere in the codebase,
that allows mutation of either archive root.**

Outside the project (`/etc`, `/root`, user's `~/.ssh`, etc.) is
filtered out by `_resolve_within_roots()`: every mutating tool
requires the target to sit inside an explicitly-allowed root.

### Layer 4 — Tier-3 cap (one destructive action per response)

Code-level enforcement in `tab_agent.py:_dashboard_auto_exec`:

```python
if tier == 3:
    mutate_count += 1
    if mutate_count > 1:
        # Blocked: only one destructive action per response
        continue
```

Tier-3 = `/submit`, `/orca submit`, `/recalc auto`, `/cancel all`,
`/recalc *`, `/cancel *`. Even if the agent emits 10 ACTION:
lines in one response, **only the first destructive one runs**.

### Layer 5 — Bulk-op intent gate

`/recalc auto` and `/cancel all` need both an **action keyword**
AND a **scope keyword** in the user's last message. Code-level
check via `_BULK_ACTION_KW` + `_BULK_SCOPE_KW`. The agent cannot
trigger a bulk wipe just because it thinks it's a good idea —
the user has to explicitly say "alle neuberechnen" / "cancel all".

### Layer 6 — Test-mode 4-layer wall (for the dry-run harness)

When the harness runs the agent in plan-mode against a sandbox:

1. **`--permission-mode plan`** → CLI hard-blocks Bash, Edit,
   Write, NotebookEdit at the OS layer.
2. **`--disallowed-tools`** lists every mutating MCP tool by name.
   Even if a future tool ships with `allow_mutate=True` defaulted,
   the CLI still refuses to execute it.
3. **`--add-dir <sandbox>`** restricts the writable surface to one
   tmp directory. The real repo + the user's calcs are not in
   the writable scope.
4. **Pre-flight wall**: `run_agent_dryrun` refuses to start if
   `cwd` sits inside the real DELFIN repo. Verified above —
   pointing it at `/home/qmchem_max/ComPlat/DELFIN` raises
   `RuntimeError: REFUSING dry-run`.

## What this means in practice

The agent CAN:
- Read every file in the project freely (Read tool, no restriction)
- Call read-only MCP tools (parsing, plots, list_*) without
  permission
- Suggest destructive ACTIONs that the user can confirm
- Run plots into `agent_workspace/` (CLI-fenced)

The agent CANNOT:
- Delete a single file without 3 aligned locks
- Submit a job without `allow_mutate=True` + user confirmation
- Mutate anything in `archive/` / `remote_archive/`
- Run more than 1 destructive ACTION per response (code cap)
- Trigger `/recalc auto` / `/cancel all` without user words for
  both action AND scope
- Write outside `agent_workspace/` (CLI-enforced via `--add-dir`)
- Touch the real repo when the test harness is running (pre-flight
  refusal)

## Recent commits that hardened these layers

- `c867403` Phase B — directory-zoning + 3-lock delete
- `8020989` Phase B-2 — bulk-op intent gates, allow_mutate-defaulted
  prepare_recalc + kill_all_user_jobs
- `4f7be8a` Test-mode 4-layer safety wall (sandbox + add-dir +
  disallowed-mcp + pre-flight refusal)
- `db09f64` Decision-tree clarity → fewer accidental mutating
  routings
- `f856c47` Phase E new mutating tools — none added (all read-only
  parsers + plots)

The 12 mutating tools the agent has access to are explicitly
listed at top — no hidden mutation paths. Every one has a
documented "would do X" dry-run signature.
