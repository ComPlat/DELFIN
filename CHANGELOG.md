# Changelog

All notable changes to DELFIN will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-07-23

### Added — Agent (long-session continuity + production parity)

A multi-sprint pass on the dashboard agent. Solo mode is now suitable for
multi-day sessions: nothing the agent learns or does is silently lost,
the user can fork / branch / search across history, and the canonical
slash-command, hook, and MCP surfaces are all live.

**Context management**
- Live context-window bar above the chat input (blue → orange → red
  ramp at 60 % / 80 %; thin marker at the 95 % auto-compact trigger).
- LLM-quality compaction (API backends) with structured Goal / Decisions
  / Open-items summarisation; falls back to extractive on the CLI
  backend to keep the persistent session untouched.
- Tool-result auto-truncation via `_smart_truncate(cap=5000)` so 200 kB
  MCP outputs no longer detonate the next request's input budget; head
  + tail kept so Python tracebacks always survive.
- Live per-turn footer reports duration / token Δ / tool calls / cost Δ
  after every turn.
- Pre-compaction transcript archive (append-only JSONL under
  `~/.delfin/transcript_archive/<sid>.jsonl`). Long sessions never truly
  lose history — only summarise it in-place.
- New `/context`, `/cost`, `/usage` slash commands; status block injected
  into the solo system prompt every turn (msgs / tokens / % of window /
  last compaction).

**Session**
- Full state save+restore: chat, perm profile, provider, model, effort,
  active gate, last-compaction info, subagent panel, pending plan body,
  todo payload (10 fields total — legacy sessions still load).
- `/session ls | restore <id> | search <query> | fork [<id>] | tree
  [<id>] | archive | archive <id>` slash commands.
- `delfin-voila --resume <sid>` (or `--resume latest`) auto-loads a
  saved session at dashboard boot; also honours
  `DELFIN_RESUME_SESSION` env var for wrapper scripts.
- Headless CLI: `python -m delfin.agent.cli run/init/session ls/session
  search`, suitable for CI hooks, nightly summaries, scripted refactors.

**Memory (typed, Claude-Code-compatible layout)**
- Four memory types — `user` / `feedback` / `project` / `reference` —
  with auto-classifier when `/remember` is called without an explicit
  prefix (e.g. *"don't mock the DB"* → feedback).
- `/remember [type:] <text>` writes both the legacy
  `~/.delfin/agent_memory.json` and a typed sidecar file under
  `~/.claude/projects/<slug>/memory/<type>_<slug>.md` with frontmatter,
  plus an indexed line in `MEMORY.md`.
- `[[name]]` wiki-style cross-link resolution at prompt-build time:
  resolved targets become inline markdown links, dangling targets get a
  ``(not yet written)`` annotation.
- `/memories verify` greps every memory file for path-shaped references
  that no longer exist on disk — flags rotten recommendations.

**Plan Mode**
- New `plan` mode in the dropdown locks the permission profile to
  read-only and injects `plan_mode_addendum.md` into the system prompt.
- `ExitPlanMode` tool round-trip persists the approved plan to
  `~/.claude/projects/<slug>/plans/<slug>.md` with frontmatter.
- `/plans` slash command (list / view / delete saved plans).

**Subagents**
- Built-in presets `explore` / `plan` / `code-reviewer` /
  `general-purpose` extended by markdown files in
  `delfin/agent/pack/agents/*_subagent.md` (cannot override built-ins)
  and `~/.delfin/subagents/*_subagent.md` (can override). First
  user-facing custom preset shipped: `chemistry-reviewer`.
- Optional `isolation="worktree"` parameter creates a fresh git
  worktree before deriving sub-perms, so subagent edits never touch
  the user's working tree until reviewed.
- Per-run telemetry JSONL log + `/agents stats` for lifetime
  aggregates (count, mean duration, tokens, errors, truncations).

**Skills + custom slash commands**
- `/skills` lists discovered skills; `/skills <name>` prints the body
  + source path + usage hint. Palette autocomplete bug fixed
  (referenced non-existent fields, crashed on first skill).
- User-defined slash commands as markdown templates: drop
  `~/.delfin/commands/<name>.md` or
  `<workspace>/.delfin/commands/<name>.md` with optional frontmatter
  (`description`, `argument-hint`). `$ARGUMENTS` / `$@` / `$1` / `$2`
  / `${ARGUMENTS}` are substituted; unknown placeholders preserved.

**Project bootstrap**
- `/init` (and `delfin-agent init`) scans the repo for Python /
  Node / Rust / Go indicators + Makefile targets, writes a populated
  `AGENTS.md` with detected build/test/lint commands, and scaffolds
  `.delfin/settings.json`. `AGENTS.md` is auto-loaded by the prompt
  loader on the next turn.

**Hooks** (PreToolUse / PostToolUse / UserPromptSubmit / Stop)
- All four events now fire end-to-end. `UserPromptSubmit` can block
  the send via `decision: block`; non-blocking hook output (exit code,
  stderr, duration) surfaces in chat so the user sees activity.
- `/hooks list / add / remove / dry-run` slash commands edit
  `~/.delfin/settings.json` without hand-editing JSON.

**MCP servers**
- `/mcp list / add / remove / enable / disable / reload` slash
  commands manage `~/.delfin/mcp_servers.json`. Add / remove / toggle
  auto-reload the process-wide registry; new MCP tool surfaces
  register without a dashboard restart.

**Tooling**
- 45 native function-calling tools (`/tools` for the categorised
  browser): filesystem, bash, web, notebook, subagent, planning,
  task, skill, user interaction, scheduling, code navigation, docs,
  permissions. Plus `task_get` for runtime task lookup.
- Diff preview for `Edit` / `Write` / `multi_edit` rendered in the
  approval row (`difflib.unified_diff`, head + tail capped at 120
  lines).

**UX / observability**
- Spinner with four modes — streaming (blue) / gated (orange,
  paused on approval) / queued (grey) / stale (red, no output for 10
  min). Header chip mirrors the activity-tab spinner exactly.
- Silent-exit detection: when the worker finishes with no output
  and no exception, surfaces an informative system message so the
  user isn't left wondering why the spinner just disappeared.
- AskUserQuestion options accept an optional `preview` markdown
  field; single-select questions with previews render side-by-side
  (40 % option list / 60 % preview pane) for visual comparison.
- `/bash watch <id>` streams `bash_background` job output to chat in
  real time; auto-stops on job exit. `/bash jobs` / `/bash unwatch` /
  `/bash kill` for fleet management.
- `solo_agent.md` extended with "Strategies for approaching tasks",
  "Subagents", "Parallel tool calls", "Context management", "Memory",
  "Skills", "Sandbox security boundary" sections.

### Added — Agent (multi-day session hardening)

Follow-up to the parity sprint focused specifically on long-running
multi-day usage: making each turn cheaper, keeping context useful as
it grows, refusing to forget user goals at compaction, and learning
from repeated failures.

**Prompt processing**
- *Progressive disclosure stripper*. Heavy lazy sections of the solo
  role prompt (chemistry decision tree, web research, notebook
  handling, KIT sandbox, project-dev workflow, background jobs) are
  marked with `<!-- module:NAME -->` comments. The prompt loader's
  new `_strip_lazy_modules` runs a task-keyword heuristic and removes
  inactive modules before injection. **Measured: 10 894 → 7 558
  tokens (−31 %)** on a typical no-keyword turn; chemistry / web /
  notebook tasks get their module back. Pipeline roles (quick /
  reviewed / cluster / full) keep the full prompt — they're more
  sensitive to subtle context changes. Opt-out via user setting
  `agent.slim_prompt: false`.
- *Cache-aware section order*. The per-turn-variable session-
  environment block (cwd, git status, recent commits) moved from
  position 3 in the solo prompt to the bottom (just before live
  state + anchor), so the stable role-prompt + project-context
  prefix can stay in the provider prompt cache between turns.

**Context management**
- *Sliding-window proactive compaction*. New
  `AgentEngine._should_slide` / `_slide_window_trim` fire at 70 % of
  the window (configurable via `_SLIDING_WINDOW_PCT`). User messages
  survive verbatim; assistant payloads >2 kB get head + tail
  truncation with a middle marker. The cliff-edge full compaction
  at 95 % becomes much rarer because the sliding window keeps
  trimming proactively.
- *Selective full compaction*. When full compaction does fire,
  user-role messages now contribute up to 400 chars to the summary
  (not just the first line); tool-role messages are dropped except
  for JSON-error markers; assistant text continues to be extractive-
  summarised. Goal: the post-compaction agent remembers WHY it was
  doing what it was doing.

**Task workflow**
- *Task dependency graph*. `TaskStore.create` gains `blocked_by`;
  `TaskStore.update` gains `add_blocked_by` / `remove_blocked_by`
  with automatic reverse-index maintenance on `blocks`.
  Transitioning to `in_progress` is refused while any blocker is
  unfinished. The `task_create` tool schema accepts the new field.
- *Stalled-task detector*. After every worker turn, the dashboard
  checks `TaskStore.find_stalled(max_age_s=600)` against the current
  session; if any in-progress task has been untouched for >10 min the
  user gets a one-shot suggestion to switch into `/mode plan` or
  `/forget` the task.

**Learning loops**
- *Retrospective failure log* (`delfin/agent/failure_log.py`). Every
  tool error gets written to `~/.delfin/failure_log.jsonl` (best-
  effort, append-only, trimmed at 2 000 entries). Helpers
  `read_failures`, `top_recurring` (≥N occurrences in 7 days), and
  `detect_repeat_for_current_task` for retry-loop short-circuits.
  `api_client._dispatch` writes the log automatically when a tool
  returns `{"error": …}`. New `/failures` slash command (default
  `top`; also `recent`, `clear`).
- *BM25-style memory rerank*. `memory_store.format_memory_context`
  now scores every fact against the query via the standard BM25
  saturation curve (k1=1.2, b=0.75) with recency as a tiebreaker.
  Replaces the previous token-set Jaccard score, which surfaced
  unrelated facts whenever they happened to share stopword-like
  tokens. Materially more relevant memories reach the prompt on
  long sessions where facts have accumulated.

### Added
- **Unified Calculator Factory** (`delfin/calculators.py`): 34 computational backends (8 MLP, 20 QM, 3 semi-empirical, 3 MD) accessible through a single `create_calculator()` interface returning standard ASE Calculator objects.
- **ML Potential backends**: CHGNet, M3GNet/MatGL, SchNetPack, NequIP/Allegro, and ALIGNN added alongside existing ANI-2x, AIMNet2, and MACE-OFF, with CUDA device validation and automatic CPU fallback.
- **AI/ML tools integration** (`delfin/ai_tools/`): 21 tools across Foundation Models (MoLFormer, Uni-Mol, ChemBERTa), Generative (REINVENT4, SyntheMol), Conformers (GeoMol, torsional-diffusion), Crystal Generation (MatterGen, CDVAE), Retrosynthesis (AiZynthFinder, RXNMapper, LocalRetro), Screening (DeepChem, ADMETlab), Metal Complex ML (molSimplify, architector), Visualization (plotly), and Wrapper Libraries (pymatgen, QCEngine, MDAnalysis, pymolpro).
- **Analysis tools**: cclib (output parsing), nglview (3D visualization), and Packmol (solvation box builder) with wrapper modules.
- **Per-tool Install/Update buttons** in Dashboard Settings tab for all pip-installable tools, with outdated-package detection via `pip list --outdated`.
- **Auto-detection of 88 external programs**: Turbomole, Gaussian, NWChem, Q-Chem, GAMESS, Molpro, VASP, Quantum ESPRESSO, CP2K, FHI-aims, GROMACS, LAMMPS, AMBER, NAMD, OpenMM, and many more via PATH search (compatible with cluster module systems).
- **ARCHITECTOR button** in Submit Job tab: converts metal-complex SMILES to 3D XYZ structures using the architector library, with isomer navigation and 3D viewer display.
- Shell-based installers for all tool categories (`install_mlp_tools.sh`, `install_ai_tools.sh`, `install_analysis_tools.sh`) with per-tool environment variable control.

### Fixed
- Fixed xTB calculator missing OrcaProfile (caused ASE BadConfiguration error).
- Fixed MOPAC calculator unreachable UHFSINGLET branch.
- Fixed OpenMM calculator dead code and confusing error flow.
- Fixed OpenMolcas factory importing wrong module (openmx instead of openmolcas).
- Fixed AIMNet2 and M3GNet silently ignoring `device` parameter (GPU requests had no effect).
- Fixed ALIGNN ignoring `model_name` parameter.
- Fixed Dashboard Update button using `install_cmd` instead of `pip_pkg_name` (wrong package updated).
- Fixed Dashboard pip command parsing using `str.split()` instead of `shlex.split()` (broke on complex package specs).
- Fixed operator precedence in wrapper library version lambdas (`or` vs `if/else`).
- Fixed missing ADMETlab and LocalRetro in AI tools installer script.
- Added `logging.exception()` to Dashboard refresh error handlers (errors were silently swallowed).
- Added installer return code check in `setup_runtime()` for QM tools.

## [1.1.9] - 2026-03-16

### Added
- New runtime and setup controls in the Settings tab, including configurable `Calculations`/`Archive` paths, backend selection, ORCA path overrides, ORCA scan/selection, and system-aware local resource detection.
- New packaged `qm_tools` workflow controls in the Settings tab: `Prepare qm_tools`, `Install qm_tools`, and `Update qm_tools`.
- New bwUniCluster actions in the Settings tab: `Setup bwUniCluster`, `Verify bwUniCluster`, and `Full bwUni install`.
- New packaged runtime resources for wheel/PyPI installs, including generic submit templates and the packaged bwUniCluster installer.
- New English setup documentation in `docs/SETTINGS_AND_SETUP.md`.

### Changed
- DELFIN runtime configuration is now stored in `~/.delfin_settings.json`, keeping user-specific paths and runtime choices outside the git repository.
- Local execution now uses a bundled Python local runner fallback when shell-based local runner scripts are not present.
- `delfin-voila` now stages packaged notebooks safely, avoids fragile server-side browser assumptions, and behaves better in remote/server sessions.
- Local and cluster setup flows are now exposed directly in the GUI while preserving compatibility with editable installs, existing repo-based installations, and established `software/` layouts.
- Packaged installs now include the resources needed for runtime setup instead of assuming a neighboring source checkout.

### Fixed
- Fixed packaged `delfin-voila` launch failures caused by notebook root-directory issues, hidden staging paths, and server/browser launch behavior.
- Fixed missing packaged runtime resources for PyPI/wheel installs, including local runner and setup assets that were previously only available from a source checkout.
- Fixed Settings-tab issues around saving configurable workspace/runtime paths and improved masking/handling of local sensitive transfer settings.

### Documentation
- Added detailed documentation for Settings, runtime controls, setup buttons, local installs, PyPI installs, and bwUniCluster workflows.

## [1.1.0] - 2026-03-03

### Added
- Dashboard/GUI expansion across Submit/Recalc/ORCA Builder/TURBOMOLE Builder/Calculations/Archive tabs, including improved archive-to-calculations workflows and MO plotting controls.
- Explicit SMILES conversion modes in the Submit tab (`CONVERT SMILES`, `QUICK CONVERT SMILES`, `CONVERT SMILES + UFF`) with improved handling of metal complexes.
- `delfin-guppy` workflow improvements: quick conversion as additional start geometry and post-XTB topology validation.
- Extended coordination-chemistry support in SMILES conversion (including additional coordination numbers and improved topology/isomer enumeration robustness).
- `ESD_T1_opt` toggle support for controlling UKS vs TDDFT T1 optimization behavior.
- Smart recalc controls and fingerprint-based skip logic for avoiding unnecessary ORCA reruns.

### Changed
- Isomer deduplication logic was tightened and then made workflow-aware: dashboard and GUPPY now preserve broader labeled variant diversity for SMILES isomer sets.
- Dashboard molecule viewer and trajectory UX refined (layout stability, playback reliability, print mode controls, cube/MO interaction refinements).
- Cleanup/copy-back behavior around `.orca_iso*` and GOAT/XTB side products hardened for better restart safety.
- HPC runtime I/O overhead reduced via caching and scratch/runtime optimizations.

### Fixed
- Multiple SMILES conversion regressions affecting metal complexes (fragment topology, ring-count checks, hydrogen handling, and fallback strategy behavior).
- Archive/statistics browser issues (path extraction, folder navigation toggles, and clipboard/export actions).
- Several ORCA input/output propagation issues in submit/recalc flows (including missing copied `.inp` and timeout/abort recovery paths).
- Additional CO2 coordinator input-generation edge cases (basis assignment and convergence keyword handling).

## [1.0.4] - 2025-02-XX

### IMAG improvements
- CLI: add `delfin --imag` entry point to rerun the IMAG workflow on existing results and regenerate DELFIN.txt.
- IMAG inputs now reuse the original ORCA template (incl. `%QMMM`, per-atom `NewGTO`) while stripping `MORead` / `%moinp` hints to force fresh guesses.
- `%pal` and `%maxcore` are inherited from the source frequency job so IMAG iterations use the same resource allocation.
- Classic, manually, and OCCUPIER pipelines call the refactored IMAG routine for initial and redox steps, copying original inputs/outputs for traceability.
- README updated with `--imag` documentation and CONTROL option `IMAG_scope` (initial/all).

### Added
- Optional excited-state dynamics (ESD) module with standalone or post-redox execution, including S0/S1/T1/T2 optimisations plus ISC/IC scheduling in a dedicated `ESD/` directory.
- ORCA input builders for ESD workflows (`esd_module.py`, `esd_input_generator.py`) and CLI/pipeline switches (`ESD_modul`, `states`, `ISCs`, `ICs`) with documentation.
- Automatic per-run logging: the CLI attaches a global `delfin_run.log` and OCCUPIER subprocesses emit an `occupier.log` alongside ORCA outputs.
- CLI help now documents QM/XTB splitting via `$` markers and clarifies how `parallel_workflows` toggles between parallel and sequential scheduling modes.
- `delfin --purge` command to wipe all intermediates (keeps CONTROL.txt + main input) with interactive confirmation.

### Changed
- Persist QMMM split detection (`$` separator) via shared cache so OCCUPIER/Classic/Manually runs keep the `QM/XTB` flag even if later geometries omit the marker.
- All geometry writers now pass the geometry path to the splitter so cache lookups work across workflow steps.
- OCCUPIER FoB scheduling always uses the global job manager; sequential runs simply cap `max_jobs` to one, ensuring consistent PAL enforcement.
- `delfin --purge` now deletes only recognized DELFIN artifacts (OCCUPIER folders, ORCA inputs/outputs, logs) and refuses to touch unknown files.
- Global job manager initialization is now idempotent and refuses to tear down an active pool when jobs are still running; configuration differences are detected via resource signatures.
- Improved banner and logging output for the global scheduler, including explicit parallel-mode reporting.
- Added `orca_parallel_strategy` CONTROL option allowing OCCUPIER runs to force serial ORCA execution (`threads`/`serial`) when MPI backends are unstable.

### Fixed
- Ensured newly spawned OCCUPIER steps reuse the cached QM range, restoring the `%QMMM` block for oxidation/reduction inputs.
- Sanitized PAL/maxcore/pal_jobs parsing so subprocesses inherit consistent limits from CONTROL.txt and `parallel_workflows`.

## [1.0.3] - 2025-01-XX

### Added
- **Global job manager singleton** for centralized resource coordination across all workflows
- **Automatic PAL splitting** for parallel oxidation/reduction workflows (e.g., PAL=12 → 6 cores per workflow)
- **Thread-safe workflow execution** with proper locking and coordination mechanisms
- **Subprocess PAL coordination** via `DELFIN_CHILD_GLOBAL_MANAGER` environment variable
- **Bootstrap mechanism** for OCCUPIER subprocesses to respect global resource limits
- New module: `global_manager.py` - Singleton GlobalJobManager class
- New module: `thread_safe_helpers.py` - Thread-safe workflow preparation functions
- New module: `cluster_utils.py` - Automatic cluster resource detection (SLURM/PBS/LSF)
- **Verification script** `verify_global_manager.py` for testing singleton behavior

### Changed
- **Refactored parallel execution** to use shared `DynamicCorePool` across all workflows
- **OCCUPIER subprocesses** now receive allocated PAL via environment variables instead of reading original CONTROL.txt
- **PAL management** centralized - read once at startup, propagated consistently to all jobs
- **CONTROL.txt updates** in subprocess folders now include reduced PAL values
- Improved logging with visual banners for parallel workflow execution
- Updated `parallel_classic_manually.py` to use global pool with intelligent fallback
- Updated `parallel_occupier.py` for global pool integration

### Fixed
- **Fixed critical issue**: Double core allocation when ox/red workflows run in parallel
- **Resolved race conditions** in parallel workflow execution through global coordination
- **Prevented PAL over-allocation** - total CPU usage never exceeds configured PAL
- Improved thread safety in file operations (CONTROL.txt, geometry files)
- Fixed subprocess resource management to respect global limits

### Documentation
- Added "Global Resource Management" section to README.md
- Expanded "Parallel Processing" section in docs/methodology.md with implementation details
- Updated project layout in README.md with new modules
- Added architectural details for global job manager pattern

## [1.0.2] - 2025-01-26

### Added
- Initial public release
- OCCUPIER, classic, and manually calculation modes
- Automated spin-state prediction and redox potential calculations
- Integration with ORCA 6.1.1, xTB, and CREST
- Support for up to 3 sequential oxidation/reduction steps
- Parallel workflow execution for oxidation/reduction steps
- Cluster environment detection (SLURM/PBS/LSF)
- CPCM/SMD implicit solvation models
- TD-DFT excited state calculations
- Broken-symmetry DFT for transition metal complexes

### Documentation
- Comprehensive README.md with installation and usage instructions
- Detailed methodology.md documentation
- Example job submission scripts for SLURM/PBS/LSF
- CITATION.cff for academic citation

---

## Release Notes

### Version 1.0.3 - Global Resource Management

This release introduces a major architectural improvement: **global job management**. Previously, when oxidation and reduction workflows ran in parallel, each could potentially allocate the full PAL (number of CPU cores), leading to double allocation (e.g., 2×12=24 cores when only 12 were available).

**Key Improvements:**
- Single source of truth for PAL allocation
- Automatic core splitting when workflows run in parallel
- Guaranteed compliance with cluster resource limits
- Thread-safe coordination between all workflows

**Impact:**
- More efficient cluster resource utilization
- No more job failures due to over-subscription
- Predictable performance in parallel execution mode

**Migration:**
No changes required to existing CONTROL.txt files or workflows. The global job manager is automatically initialized and works transparently.
