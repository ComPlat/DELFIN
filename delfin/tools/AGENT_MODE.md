# DELFIN Pipeline-Builder Mode — binding guide

The platform side (this `delfin/tools/` package) is complete and tested. What
remains is wiring an **agent mode** that drives it. The mode itself lives in the
agent framework (`delfin/agent/`) — this file is the contract it binds to.

## 1. Register the MCP server

The agent "sees" the platform tools once `delfin-tools` is in its MCP config
(stdio). Entry point: `delfin-tools-server` (= `python -m delfin.tools.mcp_server`).

```json
{
  "mcpServers": {
    "delfin-tools": { "command": "python", "args": ["-m", "delfin.tools.mcp_server"] }
  }
}
```

## 2. Tool inventory (over MCP)

- **discover** — `get_guide`, `get_manifest`, `list_capabilities`,
  `describe_capability`, `catalog`, `compatible_successors`, `list_keys`,
  `describe_key`
- **build** — `resolve_spec` (fill defaults + show auto-wiring), `validate_spec`
  (structural diagnostics with concrete fixes), `scientific_lint` (physics
  plausibility)
- **extend** — `new_capability_template`, `register_module` (build + integrate a
  missing building block)
- **persist** — `save_application` (→ appears in the Pipelines tab)
- **run** — `run_application`, `submit_application`, `run_status`, `run_metrics`
- **diagnose** — `run_diagnostics` (status, error, events, `~/calc` log files)
- **environment** — `probe`, `install_plan`

## 3. Mode system prompt (paste into the mode)

```text
You are DELFIN's Pipeline-Builder. You work exclusively through the MCP server
"delfin-tools". Call get_guide first, then follow this loop:

1. DISCOVER  – describe_capability / catalog / describe_key (use only allowed
               enum values); compatible_successors for step order.
2. BUILD     – write a pipeline spec to schemas.pipeline_spec; set only the
               required params (defaults + auto-wiring fill the rest). Inspect
               with resolve_spec.
3. VALIDATE  – validate_spec (structural) AND scientific_lint (physics); apply
               every suggested fix until clean.
4. NO TOOL?  – if no building block fits: new_capability_template → fill the
               execute() body → register_module (writes, loads, verifies). Back to 2.
5. PERSIST   – save_application (it then shows in the Pipelines tab).
6. RUN       – submit_application (local or SLURM); results land in ~/calc.
7. DIAGNOSE  – run_diagnostics(run_id): study status/error/logs, then iterate.

Rules: prefer allowed enum/key values; always validate before running;
ORCA and Turbomole are license-restricted and are NEVER installed (guidance only).
```

## 4. What lands where

- Built workflows → **Pipelines tab** (via `save_application`).
- Run results → **`~/calc/<app>_<id>/`** (the Tools tab's run-detail view +
  `run_diagnostics` both read this).
- User-built blocks → **`~/.delfin/adapters/*.py`**; user apps →
  **`~/.delfin/applications/*.json|*.py`** (auto-discovered).

## 5. Boundaries (do not cross)

- The CONTROL.txt-driven workflow, the engines (`run_OCCUPIER` / ESD / IMAG), and
  the agent framework are owned elsewhere — the platform only *derives* building
  blocks from them, never modifies them.
- License-restricted engines (ORCA, Turbomole) are never auto-installed.
