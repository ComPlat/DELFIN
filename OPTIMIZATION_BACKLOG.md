# DELFIN Agent — Optimization Backlog

Durable work-queue for the self-paced optimization loop. The loop draws tasks
in this priority order and appends new findings back here.

## Loop protocol (each iteration)
1. **Sync & baseline** — `git pull`/fetch; run the full suite. Establish that
   the only expected failure is the pre-existing `test_runtime_setup.py::…xtb4stda`
   (missing binary bundle). Anything else red → fix that first.
2. **Pick ONE task** by priority: open bug reports → a regression caught by the
   smoke → the highest-value backlog item below → a weak spot from
   telemetry/traces.
3. **Implement** small + atomic; **test** (full suite + a live check where it
   matters). **Never commit on red.**
4. **Ship** — commit to `main` (no Co-Authored-By trailer), push; move a fixed
   bug report to `…/AGENT_BUGS/Solved`; tick the item here + note any new
   finding.
5. **Repeat.** Keep each iteration independently green and reversible; route
   anything risky/large to review instead of `main`.

## Backlog (curated; newest insight wins)

### Model adaptivity
- [ ] **Weak-model detection misses small multi-modal/reasoning tags** — e.g.
  `qwen3-vl:4b` (4b) isn't matched by `prompt_loader._WEAK_MODEL_PATTERNS`
  (qwen pattern lacks `4b`), so a 4B model gets the full 45-tool surface +
  full prompt. Extend the patterns / prefer the capability layer's window+size
  signal for the weak/strong split.
- [ ] **KIT live context-window detection needs auth** — `model_capabilities.
  _fetch_openai_models` sends no `Authorization` header, so KIT `/v1/models`
  (which carries `max_model_len`) 401s and we fall back to the static 128k.
  Thread the provider key into the probe so KIT windows are detected live.
- [ ] **Live `available_models` into `route_model`** — the mechanism exists
  (`route_model(available_models=…)`) but nothing feeds the live model list, so
  routing can't avoid an unlisted/removed model. Wire the dashboard's fetched
  list in.

### Multimodal / providers
- [ ] **Claude-format vision** — `engine._build_user_message` only emits
  OpenAI `image_url` content (Ollama/KIT/OpenAI). Claude uses Anthropic image
  blocks → currently text-fallback. Add the Anthropic shape for the
  CLI/API backends.
- [ ] **`standard-extern` KIT alias 400s** — surfaced cleanly (no crash) but
  unusable. Detect non-chat aliases and skip/flag them in the model picker.

### Observability / self-improvement
- [ ] **Tool-trace → metrics/UI** — feed `tool_trace` into `agent_metrics`/
  `outcome_tracker`; add a dashboard panel for the live trace (not only
  `/trace`).
- [ ] **MCP resources/prompts UX** — beyond the agent meta-tools, surface MCP
  prompts as slash commands + resources as @-mentions in the dashboard.

### Done (recent, for context)
- [x] read_file/grep/list no longer require the doc index (off-repo subagents).
- [x] reasoning-model token floor; reasoning channel surfaced.
- [x] vision wired (OpenAI-compatible backends).
- [x] memory consolidated on the typed store under `~/.delfin`.
- [x] honesty addendum; parallel-writer auto-isolation; tool tracing.
