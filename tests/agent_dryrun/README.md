# Agent dry-run capability harness

Test environment that asks the **actual dashboard / solo agent** realistic
German queries via the local OAuth-authenticated `claude` CLI, captures
which tool the agent picks, and asserts whether the right tool got
chosen — without ever executing a mutation.

Built to find the bigger class of bug that unit tests can't see: the
agent picks the wrong tool / wrong syntax even though the right one
exists and is well-tested. Concrete example: the user typed `batch
alle initial.xyz Dateien aus calc/`. The agent emitted
`/batch from-calc initial.xyz`, but the optional argument is a folder-
name glob (matching zero folders) — not a file pattern. The agent had
the right tool, used it wrong, and the unit tests of `/batch from-calc`
all passed. Only a behavioral test catches that.

## Safety — the agent CANNOT damage the host

Four layers, defense-in-depth:

1. **`--permission-mode plan`** — the CLI hard-blocks `Bash`, `Edit`,
   `Write`, `NotebookEdit`. The model emits tool-use blocks but the CLI
   denies their execution before they hit the OS.
2. **Mutating MCP tools blocked by name** — every tool that could write
   files / submit jobs / scancel / mkdir / rmtree is added to
   `--disallowed-tools`. Even if a future tool ships with
   `allow_mutate=True` defaulted by accident, plan-mode + the
   disallowed list still block it. List in
   [`runner.py`](runner.py)'s `_hard_blocked_mcp`.
3. **`--add-dir` restricted to the temp sandbox** — the agent's
   writable surface is exactly one tmp tree, never the real DELFIN
   repo or the user's calc folders.
4. **Pre-flight wall** — `run_agent_dryrun` REFUSES to start if `cwd`
   sits inside the real DELFIN repo. A typo can't redirect the agent
   at the user's actual code/data.

## Quickstart

```bash
# Run all 34 capability tests against haiku, capture markdown table
DELFIN_AGENT_DRYRUN=1 python -m tests.agent_dryrun.report \
  --models haiku --out /tmp/report.md

# Multi-model report (compare which queries haiku misses but sonnet hits)
DELFIN_AGENT_DRYRUN=1 python -m tests.agent_dryrun.report \
  --models haiku,sonnet --out /tmp/report.md

# Filter to a single test by substring
DELFIN_AGENT_DRYRUN=1 python -m tests.agent_dryrun.report \
  --models sonnet --include "imag-freq" --out /tmp/quick.md
```

Cost: ~$0.05–0.20 per cell on haiku, ~$0.20–1 on sonnet. A full
34-case haiku sweep is ~$2-3.

## What's tested (34 capability cases)

| Phase | Coverage |
|---|---|
| Phase A | imag-freq, functional comparison, two-folder diff |
| Phase A extras | thermochem, dipole, energy table, ORCA errors, lowest-Gibbs |
| Phase B | rename, create, move (calc/calc), archive (calc/archive), delete (3-lock) |
| Phase B-2 | cancel-all bulk, smart-recalc preparation |
| Phase C | (Options) dropdown lookup |
| Phase D | HOMO/LUMO, UV/Vis, opt convergence + plots |
| ORCA Builder + lifecycle | INP validation, submit, list active jobs, SSH transfer |
| Discovery + recipes | dashboard pattern, widget catalog, DELFIN concepts, describe_tool, list_tools |
| Literature | ORCA-manual lookup, ad-hoc PDF read |
| /batch UX regression | the exact bug from the user's screenshot |

Each case is parametrized over models — set
`DELFIN_AGENT_DRYRUN_MODELS=haiku,sonnet,opus` to run every case
against every model.

## Findings shipped (2026-04-28)

The first sweep on haiku found systematic gaps in tool-selection.
Each was fixed and the fix shipped in a separate commit:

- **`tool_X` prefix in MCP tool names** — wrappers were registered
  with `mcp.tool()(tool_X)`, exposing tools as
  `mcp__delfin-ops__tool_extract_imaginary_frequencies` (duplicate
  prefix). Renamed to plain `mcp__delfin-ops__extract_imaginary_frequencies`.
  Commit: `63ea627`.
- **Haiku fell back to Glob/Grep on `.out` files** instead of calling
  the typed parsers. The dashboard agent's "Speed first" rule said
  "Typed MCP tool > grep > Python script" which read as a soft
  preference; haiku interpreted it as "either is fine". Promoted to
  a HARD rule with explicit negative examples + a per-question
  decision tree mapping every common query class to the FIRST tool
  call. Commit: `0938d37`.
- **`/batch from-calc initial.xyz` filter-vs-content gotcha** — the
  optional arg filters folder names, not file patterns. Strengthened
  the error message + `_DASHBOARD_PATTERNS["batch"]` recipe + made
  the pattern lookup mandatory before any /batch /recalc /cancel
  ACTION. Commit: `7de129c`.
- **Effort dropdown not wired to CLI** — solo and dashboard had a
  level dropdown (low/medium/high/xhigh) used locally for token
  caps but never passed to the actual `claude` CLI subprocess.
  Plumbed through `tab_agent → AgentEngine → create_client →
  CLIClient → cmd['--effort', value]`. Commit: `d1fb7e1`.
- **Solo had no session-boot primer** — dashboard had `S4` (5-line
  context block on first user-send: outcomes, jobs, commits, calc
  folder, branch). Solo started cold. Added `_build_solo_session_boot`
  with a code-work-appropriate payload (outcomes + branch + commits;
  jobs + calc-folder skipped — those are workflow concerns).
  Commit: `f15b401`.

## Architecture

```
test_capabilities.py    pytest tests, parametrized over models
report.py               markdown-table runner with per-cell incremental
                        flush to .partial (sweep-crash safe)
runner.py               run_agent_dryrun() — spawns claude CLI, parses
                        stream-json, returns AgentDryRunResult
sandbox.py              build_default_sandbox() — tmp DELFIN tree with
                        synthetic ORCA outputs (BP86/PBE0/B3LYP, varying
                        imag-freq counts)
```

`has_tool_call(result, *substrings)` is the substring-matching
assertion used in tests — it returns True if the agent attempted any
tool whose name contains one of the substrings, even if plan-mode
later blocked the call. We measure **intent**, not execution.

## When to add a new test

When you add a new typed MCP tool that an end-user would invoke via
natural language, add a CASE entry to `report.py` with a German query
that expresses the intent (not the tool name) and an `_make_picks`
assertion listing the tools that would correctly answer it. Run the
sweep, see which models pick it up, and tune the prompt if needed.
