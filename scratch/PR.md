# Package e1 (integrator): a Chemotion ELN specialist reads and writes through the API, every write approved

**Verdict: LAND, WITH ONE APPLIED PATCH PENDING AND ONE FIX NOT YET REAL-RUN.**
The specialist is offered, reads and writes are wired, and every write waits
for a person. One real run found a defect four rounds of tests had missed (an
ELN **read** was refused in default mode); its fix ships as a patch a person
must apply, and that fix has had **no real run of its own**.

## 1. Control, before any edit

Untouched tree at `origin/main`, own worktree, `delfin` imported from that
worktree (the loaded `__file__` was printed in every real run):

| suite | result |
|---|---|
| lab/ELN suites (`test_eln_route_gate`, `test_lab_chemotion_exporter`, `test_lab_edge_and_workflow`, `test_the_agpl_boundary_holds`, `test_lab_tools_registered`, `test_lab_safety_kernel`) | **71 passed** |
| repo-wide scanners | **266 passed** |

## 2. The two arms, verified — reports not trusted

| arm | what its branch actually is | verdict |
|---|---|---|
| **w44/e1-a1** (5 commits) | Real client (`eln_api.py`), route gate (`eln_routes.py`), write ledger (`eln_write.py`), tool module, ops-server registration, CLI (`delfin-agent eln`), specialist | **taken as base**, three defects fixed (§3) |
| **w44/e1-a2** | Has a chunked attachment route and a socket-based fake server, but **no reachable entry point**, and its preset names six `eln_*` tools that are not known tools — it fails `test_every_shipped_preset_is_valid` **on its own branch** | **rejected as base**; its attachment route shape informed §3 |

a1's report claimed "15 passed with the patch". Unpatched, its branch is
**40 passed, 2 failed** — the report described a patch it never shipped.

## 3. One row per item

| item | decision | commit | reason |
|---|---|---|---|
| API client: timeouts, pagination, named errors, credential from the store, `display_url` | **take a1** | `96a2beb9d` | Meets the package; credential never in a URL, log or error |
| Route gate, GET/POST/chunked shapes | **take a1** + a2's attachment shape | `81d144fb8`, `65ecee592` | Two route-gate properties survive the added routes |
| Tools: search, read, draft | **take a1** | `96a2beb9d` | Real ops-server tools |
| **`eln_attach`** | **combine** — a1's seam was dead (`send_with_approval` existed for "an attachment upload", nothing behind it) | `5aa718ceb` | The package names attach. Built in a1's idiom so it waits for a person like the draft |
| Specialist preset | **take a1, fix** | `86bf48440` | a1's `tools:` list **narrowed all five ELN tools away** (§4) |
| ELN schemas in the builtin catalogue | **reject** | `86bf48440` | **Inert** (no dispatch branch) and over the budget in four tests. The tools reach the model over MCP, as `parse_orca_output` does |
| Token ratchets | **raise by the measured amount** | `86bf48440` | The `subagent` growth is derived: the new specialist is correctly listed as a writer |
| **ELN reads refused in default mode** | **patch only** (§5) | `scratch/eln-reads-in-the-gate-allowlist.patch` | The real run's finding |
| a2's preset, a2's fake server | **reject** | — | Preset fails a CI invariant; a1's fake server already answers in the API's shapes |

## 4. Predicted vs measured (tests vs reality)

Predicted: the specialist "can reach the tools it names". **Measured: it could
not.** Its preset validated, sat in the registry, and the child's derived
permissions then refused every ELN tool:

```
eln_find  -> "not on this session's --allowed-tools list"
```

Cause: a preset's `tools:` **narrows** the child's surface, and a1 listed only
generic read tools. Fixed in `86bf48440`. Control case passing on **both**
trees — `xyz-editor` is still offered and still reaches its tools, so the fix
changed the ELN preset, not every preset.

## 5. The fix that cannot be committed

`delfin/agent/api_client.py` is the agent's own safety layer. The guard refused
`edit_file` **and** `git apply`. Routing around it is forbidden, so it ships as:

```
scratch/eln-reads-in-the-gate-allowlist.patch   (applies cleanly against this HEAD)
```

**A person must apply it.** Until then, on this branch an ELN **read** is still
refused in `mode=default`.

Measured with the patch's content rebound in memory, guarded file untouched
(`scratch/probe_gate_effect.py`):

| tool | before | after |
|---|---|---|
| `eln_find` / `eln_read` / `eln_pending` | refused | **ADMITTED** |
| `eln_draft` / `eln_attach` / `eln_approve` | refused | refused |

## 6. The real runs, judged against the prediction written before they ran

`scratch/realrun/summary.md`, base and branch HEAD, same prompt and model,
fresh HOME each.

| probe | base | branch | prediction |
|---|---|---|---|
| `specialist-offered` | no `chemotion-eln` type offered | offered; 44 tool calls vs 17 | **held** |
| `draft-write-needs-approval` | invents a draft flow with no tools | names the real contract (`eln_approve`: *"approving an ELN write is a person's command"*), refuses without a person | **held** — write gated |
| `read-question` | — | **`needs_approval`, 1 tool call, answer NOT delivered** | **FAILED** |

The read probe's own log:

```
{"type": "permission_denied", "tool": "eln_find",
 "reason_code": "mcp_side_effect", "needs_approval": true}
```

The agent answered honestly and did not route around it:

> So: the ELN read was blocked by the permission gate, not by the ELN. … I am
> deliberately **not** trying a different ELN entry point to reach the same data.

**Root cause.** `eln_find` sends a GET. The gate decides by
`_MCP_READONLY_TOOL_BASES`, and the ELN reads were never added. A pin test for
exactly this exists — `test_a_read_only_tool_is_callable_in_default_mode` — but
it slices `server.py` between the `Read-only` marker and the `P1 — statistical
plots` marker, and the ELN block sits **below** that slice. The one test that
would have caught it never looked at the ELN tools.

**This fix has had no real run of its own.** The patch is unapplied, so the
real-run evidence above does not cover it; its evidence is the gate-level probe
in §5, not a product run.

## 7. What I did NOT do

- No second `scratch/realrun.json` (instructed not to).
- Did not apply the guarded patch, and did not go round the guard.
- Did not merge, push, or touch the arms' branches.
- Did not adopt a2's preset — it fails `test_every_shipped_preset_is_valid`.
- **The 502 full-suite failures are NOT attributed to this branch**, and are
  not claimed as pre-existing either: they carry zero file-level overlap with
  the 13 files this branch touches and none is in `eln`/`lab`/`chemotion`, but
  name-level overlap is not causation. They are *unattributed*.

## 8. Every suite, with its numbers

| suite | result |
|---|---|
| lab/ELN suites (6 files) | **71 passed** — identical to control |
| repo-wide scanners | **269 passed** (control 266; +3 from the corrected wiring tests) |
| consuming + scanners after `eln_attach` and the patch test | **203 passed** |
| scanner set incl. personal-data + budget + route gate + exporter (7 files) | **155 passed** |
| `test_an_eln_read_is_reachable_and_its_writes_are_not.py` (new) | **red before, green after** — red on the unpatched tree for the named reason |
| **full suite: 3225 files, 12 chunks** | **502 failed, ~8500 passed**; **all 12 chunks exit 1** |

The runner script's exit 0 was its own `echo`, not the suite's — every chunk
exited 1. Rule 10 applied: the number quoted is the chunks' summary lines, not
a truncated tail.

## 9. The full-suite attribution

151 distinct files failed. **None** is among the 13 files this branch touches
and none is in `eln`/`lab`/`chemotion`; the failures concentrate in dashboard,
session, cage and shell suites (`test_a_file_is_shown_once_in_a_turn.py` ×13,
`test_a_dashboard_tab_hosts_its_session_in_the_service.py` ×7), the profile
this host is already known to produce.

The control over exactly those 151 files on `origin/main` is in
`scratch/suite/control_out.log`. The comparison is what makes a failure
attributable; nothing here is claimed as caused by this branch.

## 10. Next step for whoever lands this

Apply `scratch/eln-reads-in-the-gate-allowlist.patch`, then re-run the
`read-question` probe — that is the one claim still only gate-level, not
product-level.

## 9. Attributing the full-suite failures (rule 14)

The full suite ran chunked over `tests/test_*.py` (3225 files, 12 chunks):
**502 failed**, ~8600 passed, every one of the twelve chunks exited 1. The
runner's own exit code was 0 because its last statement was an `echo` — judged
from the chunks' summary lines, not from `exit_code`.

Those 502 were **not attributed** at first: the failing files (151 of them)
share no *name* with the files this branch touches, which is not the same as
being pre-existing.

Control: the same 151 files, on an untouched `origin/main` worktree, same
interpreter, same `-m "not slow"`, same `--timeout=300`:

| side | result |
|---|---|
| origin/main (control) | **482 failed, 2073 passed**, 19 skipped, 39 errors, 3175 s |

So at file granularity nearly all of the branch's 502 already fail on
`origin/main` — dashboard, service, cage, shell and peer-text suites this host
has a recorded history of failing in a sandbox. **Prediction (written before
the measurement):** re-running the branch over exactly those 151 files will
show **0 failures unique to the branch**, and specifically none in any
`eln`/`lab`/`chemotion`/`attach`/`draft` file. If any nodeid fails on the branch
and passes on the control, that is a regression I own.

A first attempt to diff nodeids printed "0 branch failures" — a **false zero**:
the chunk logs lived in a per-command temp directory that had been cleaned up,
so the grep read nothing. An absent input is not a measurement; re-run with the
logs kept inside the workspace.
