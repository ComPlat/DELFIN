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
3. **Implement** small + atomic; **test** — full suite AND a **live check with a
   real model** (qwen via Ollama localhost or KIT `kit.qwen3.5-397b-A17b`):
   don't just build, prove the change actually works end-to-end. **Never commit
   on red.**
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

### Worktree orchestration (focus area)
- [x] Auto-isolate parallel writers into separate worktrees (done earlier).
- [x] **Review**: subagent result now carries a `diff_summary` (changed files +
  --stat + untracked) so the parent SEES what an isolated writer did, not just
  `had_changes`. `worktree.diff_summary`. Live-verified with KIT-qwen.
- [x] **Integrate**: `worktree_merge` tool brings a worktree's changes back
  into the main tree, **clean-apply-or-nothing** (stage full state → patch vs
  branch-point → `git apply --check` → apply only if clean; on conflict the
  target is left untouched, worktree kept). Lands uncommitted for review.
  `worktree.merge_worktree`. Live-verified: KIT-qwen chained
  `enter_worktree → write_file → worktree_merge`, file arrived in the main tree.
  Fan-out-writers → review → merge flow is now complete.

### Done (recent, for context)
- [x] **Busy-poll guard for `bash_status`** (bug 20260615-152119) — the
  `wait_seconds` blocking wait existed but only helped if the model used it;
  models still tight-polled a ~10-min job every 3-4s and exhausted the
  tool-round budget. Now the FIRST status check is an instant snapshot, but a
  re-check of the SAME still-running job *without* `wait_seconds` is throttled
  server-side (early-returns the instant the job ends). Model-independent.
  Live-verified two ways: real subprocess (instant snapshot → throttled
  re-poll returns at job end), and KIT-qwen drove `bash_background → 1×
  bash_status(wait_seconds) → bash_output` in 3 calls (was ~170).
- [x] `/loop` command — recurring agent loop (interval-based, scheduler-backed),
  like Claude Code's /loop. `/loop <5m|2h|1d> <prompt>` · `/loop list` ·
  `/loop stop <id|all>`. Live-verified: scheduler fires a real KIT-qwen turn.
  Also fixed: `/trace` was never routed (missing from `_SLASH_PREFIXES`).
- [x] read_file/grep/list no longer require the doc index (off-repo subagents).
- [x] reasoning-model token floor; reasoning channel surfaced.
- [x] vision wired (OpenAI-compatible backends).
- [x] memory consolidated on the typed store under `~/.delfin`.
- [x] honesty addendum; parallel-writer auto-isolation; tool tracing.
