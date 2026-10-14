# TUI rewrite: implementation plan

Companion to [TUI_REWRITE.md](TUI_REWRITE.md), which holds the agreed scope.
The phases build on each other, and each one ends with something you can run.

Sizes are relative: **S** small, **M** medium, **L** large.

## Principles

- **Upstream stays easy to merge.** New Python code lives in new files under
  `delfin/agent/service/`. Edits to existing upstream files are limited to
  registering the new subcommands, so merging `upstream/main` stays easy.
- **No new Python dependencies.** `asyncio`, `json` and `socket` from the
  standard library are enough.
- **TypeScript dependencies are pinned exactly:** `@opentui/core` 0.6.4,
  `@opentui/solid` 0.6.4, `solid-js` 1.9.12, `web-tree-sitter` 0.25.10.
- **A scripted fake engine exists from phase 1 on.** The TUI can then be built
  and tested without an API key and without token cost.
- **New-TUI sessions are ordinary DELFIN sessions.** They use the same session
  store, reports and safety gates. So `delfin-agent chat --resume`,
  `delfin-agent report`, `delfin-agent sessions`, `delfin-agent stop-all` and
  the dashboard keep working with them.
- **Approvals are never answered automatically.** A session whose question has
  no window attached waits until someone attaches.
- **Smoothness is measured, not hoped for.** Phase 7 turns the performance
  budgets into tests.

## Phase 0: Spike. Is OpenTUI's inline mode good enough? (S)

**Goal:** a go/no-go on OpenTUI's split-footer mode before anything is built
on it.

**Build:** `tui/spike/inline.tsx`, driven by fake turns:

- streaming markdown, committed to the scrollback block by block, with the
  unfinished tail live in the footer;
- tool rows and a syntax-highlighted diff;
- an approval panel and a full-height agent view, so the footer grows and
  shrinks;
- divider + replay, a 1,000-line replay, resize, and a clean exit.

**Test:**

- scripted in tmux, capturing and checking both the screen and the history;
- by you in GNOME Terminal, in the VS Code terminal, over SSH from the host,
  and on macOS.

**Done when:**

- growing the footer loses or duplicates no scrollback rows (OpenTUI bug
  #1630);
- streaming shows no visible flicker;
- a resize leaves a clean footer;
- a 1,000-line replay takes under a second;
- ←, shift+tab and esc work inside tmux.

The decision goes into TUI_REWRITE.md. We fall back to Ink + React only on a
blocker that cannot be worked around.

**You can try:** `bun spike/inline.tsx` in `tui/`.

## Phase 1: Protocol and session worker, in Python (M)

**Goal:** one agent session that can be driven over a pipe, with an exact
event stream.

**Build** (in `delfin/agent/service/`):

- `protocol.py`:
  - JSON-RPC 2.0, one JSON message per line; very long lines must not break
    the reader;
  - the protocol version, and every request and event as a typed definition;
  - a JSON Schema export for the TypeScript side;
  - `docs/TUI_PROTOCOL.md` documents it.
- `worker.py`, one session around `AgentEngine`:
  - the engine is built the way `delfin-agent chat` builds it: same launch
    checks, banner facts, permissions and isolation;
  - `stream_response` callbacks become events, and tool-call IDs are assigned
    first-in-first-out;
  - it handles prompt, queue, steer and cancel, plus changes of approval mode,
    model and effort;
  - text slash commands go through `repl_commands.dispatch`;
  - the session is saved to the normal session store after every turn.
- `approvals.py`:
  - permission, plan and question callbacks become `*.asked` events;
  - they offer the same options `terminal_confirm.options_for` gives today;
  - each blocks until its answer arrives.
- `eventlog.py`: an append-only JSONL log per session, and a replay of the
  last N turns.
- `fake_engine.py`: scripted turns covering text, thinking, tools, diffs,
  approvals, errors and long output.

**Tests** (pytest, named like the repository's other tests):

- the event order of a scripted turn;
- the approval round trip;
- cancel and steer;
- a replay from the event log equals the live stream;
- a worker that exits saves its session.

**Done when:** a scripted client runs a full turn against both the fake engine
and a real KIT model.

## Phase 2: Supervisor and launcher, in Python (M)

**Goal:** sessions live in a background service and outlive the window.

**Build:**

- `rundir.py`:
  - run directory `${TMPDIR:-/tmp}/delfin-$UID/` with mode 0700;
  - a single instance per host via `flock`;
  - the socket at mode 0600, with a peer-UID check;
  - a token in `~/.delfin/run/<host>.json`.
- `supervisor.py`:
  - accepts clients and starts one worker per session;
  - routes requests and fans events out to every attached client;
  - holds pending approvals until they are answered; the first answer wins,
    and the other windows are told;
  - keeps a roster for reconnects;
  - registers in `agent_processes`, so `where` and `stop-all` see it.
- `launcher.py` and the new command `delfin-agent tui`: connect to the
  supervisor or start it detached, wait until it is ready, then hand the
  terminal to the TUI.
- `delfin-agent service status|stop|logs`.

**Tests:**

- two clients on one session;
- a client exits and its session keeps working;
- a reattach replays the session;
- recovery from a stale lock and a stale socket;
- `stop-all` ends everything.

**Done when:** a scripted client starts three sessions, disconnects,
reconnects and finds them still running.

**Needs from you:** a quick check on bwUniCluster that a detached process
survives logout, the way tmux does.

## Phase 3: TUI foundation, in TypeScript (M)

**Goal:** the skeleton every later feature plugs into, structured like
OpenCode's `packages/tui`.

**Build** (in `tui/src/`):

- `transport/`: a JSON-RPC client over the Unix socket, with reconnect and the
  version handshake.
- `protocol/`: types generated from the Python schema; a test fails when the
  two drift apart.
- `store/`: a Solid store of sessions, messages, parts, pending approvals and
  usage. Events are applied in batches every ~16 ms.
- `scrollback/`: an ordered commit queue, streaming markdown surfaces, and the
  row renderers (user prompt, tool, diff, notice, divider, banner).
- `footer/`, `keymap/`, `theme/` and `config/`: the footer shell, a keymap
  registry whose defaults are today's keys, the theme tokens, and the user
  config file `~/.delfin/tui.json`.
- A test harness: OpenTUI's headless renderer with text-frame snapshots, plus
  a fake event stream.

**Done when:** a scripted fake conversation renders the same way in the
snapshot tests and in a real terminal.

This phase needs only the protocol from phase 1, because it runs on the fake
event stream. It can start before phase 2 is finished.

## Phase 4: One session, done well (L)

**Goal:** parity with `delfin-agent chat` for a single session, good enough
for daily use.

**Build:**

- **Composer:** multi-line editing, history, ctrl+r search, `@path` and
  `/command` completion, and placeholders for large pastes.
- **Answer:** markdown committed per finished block, the live tail shown in the
  footer, and thinking shown dimmed or hidden.
- **Tool rows**, one renderer per tool type:
  - shell, search and web;
  - read, write, edit and patch, with diffs;
  - subagents, DELFIN/ORCA tools, and refusals.
- **Approvals** in the footer:
  - today's keys `y n A k R e a ?`;
  - plan approval with `d e n`;
  - numbered answers for questions.
- **Status line and banner:** a status line with spinner, elapsed time, model,
  approval mode, tokens and cost, and the startup banner.
- **Keys:**
  - esc interrupts, enter queues a message, ctrl+g steers it into the running
    turn;
  - shift+tab, ctrl+o, ctrl+t and ctrl+l as today;
  - the ctrl+c ladder, and ctrl+d to leave.
- **Commands:**
  - a ctrl+p palette;
  - native pickers for model, permissions, effort, mode and resume;
  - every other command runs in Python.
- **Resize:** after a resize settles, clear the screen and the terminal's
  scrollback, then reprint the session at the new width, capped at a
  configurable number of turns.

**Done when:** a parity checklist built from today's `/help` and `/keys` is
complete, and a real day of work on the KIT model goes through the new TUI.

## Phase 5: Agent view and many sessions (M)

**Goal:** the headline feature.

**Build:**

- **Opening it:** ← on an empty prompt expands the footer into the list,
  grouped into Needs input / Working / Done, plus Pinned. Each row shows a
  state glyph, the name, a fixed summary, the age and the cost.
- **Keys in the view:**
  - ↑↓ to select;
  - → or Enter to open a session: divider + replay of its last few turns;
  - Space to peek and reply;
  - ctrl+x to stop (a second ctrl+x deletes), ctrl+r to rename, ctrl+t to pin;
  - alt+1…9 to jump;
  - typing + Enter starts a new session in the window's folder.
- **An "Elsewhere" group** for other DELFIN sessions: the dashboard,
  `delfin-agent chat`, other login nodes. These rows are read-only; you can
  message them through the session inbox.
- **Notifications:** a needs-input badge in every session's footer, and a
  terminal notification while the terminal is unfocused.

**Done when:**

- three real sessions run side by side;
- approvals reach you from any of them;
- closing and reopening the window finds them all.

## Phase 6: Lifecycle and robustness (M)

**Goal:** it behaves well on HPC and survives failures.

**Build:**

- **Idle handling:**
  - an unattached session waiting for you sleeps after about an hour and wakes
    when you reply;
  - the supervisor exits once nothing is left.
- **Crash recovery:** workers outlive a dead supervisor briefly, then stop.
  Sessions resume from their event log, from any host.
- **Upgrades:** a version handshake, and a restart onto the new version that
  moves idle sessions over.
- **A warm spare worker,** so a new session opens instantly. Importing the
  engine takes about 0.2 s, but building one with its MCP servers takes about
  a second.
- **Resize under stress:** resize storms, and replay of very long sessions.
- **HPC details:**
  - a guard for a `noexec` `/tmp`;
  - pointers between login nodes;
  - `stop-all` integration.

**Tests:**

- `kill -9` of the supervisor and of workers;
- dropped SSH connections;
- resize storms;
- two hosts sharing `~/.delfin`.

**Done when:** the chaos tests pass, and a week of use on bwUniCluster loses no
session.

## Phase 7: Polish, docs and shipping (M)

**Goal:** it looks good, feels smooth, and installs everywhere we support.

**Build:**

- **A design pass:** the theme (dark or light, taken from the terminal),
  spacing, glyphs, `NO_COLOR`.
- **Performance budgets checked in tests:**
  - keystroke to screen in under 16 ms;
  - no full redraw while streaming;
  - footer visible within 300 ms of launch while the service is running;
  - a 1,000-line replay in under 1 s.
- **Compiled binaries** (`bun build --compile`) for linux-x64, linux-arm64,
  macOS x64 and macOS arm64 (glibc ≥ 2.17, macOS ≥ 13). They ship as platform
  wheels or a checksum-verified download.
- **CI:** `bun test`, pytest for the service, and a binary build matrix.
- **User docs:** a section in `USER_MANUAL.md` and a reference for
  `~/.delfin/tui.json`.

**Done when:** a fresh Linux install and a fresh Mac run `delfin-agent tui`
without Bun installed.

## Later

- an ACP adapter for editors (`delfin-agent acp`);
- an opt-in git worktree per session;
- session summaries written by a model;
- making the TUI the default for `delfin-agent`;
- offering the work upstream.

## Working rules

- **Independent work:** the phases are carried out without stopping. Each
  gate decision (OpenTUI go/no-go, HPC design, parity) is made and written
  down here with its reasoning.
- **Every phase is checked three ways:**
  - unit tests;
  - end-to-end tests on this VM: tmux, an emulated xterm.js, GNOME Terminal
    (with screenshots), SSH to localhost;
  - code-review agents, whose findings are fixed before the phase closes.
- **Real-model tests** use free KIT-hosted models only (`kit.deepseek-v4-flash`,
  `kit.glm-5.3-flash`, `kit.mistral-small-4-119b-a8b`). No `azure.*`,
  `google.*` or other priced model.
- **Checks needing other hardware or a human eye** (macOS, bwUniCluster,
  look-and-feel) go into a single checklist handed over at the end.
- **Upstream:** merge `upstream/main` at the start of every phase. Run the
  affected part of the DELFIN suite before every push, and the full suite at
  the end of each phase.
- **Commits:** straight to `tui-rewrite`, one per finished step, with
  sentence-style messages as in the rest of the repository, pushed after every
  step.

## Status

| Phase | State |
|---|---|
| 0 Spike | done: GO (c35d85e1; exact-row matrix and three startup fixes in phase 0b, see TUI_REWRITE.md) |
| 1 Protocol and worker | done (b314a988, reviews applied in cc1e12e0) |
| 2 Supervisor and launcher | done (96789513; two reviews, 12 bugs fixed in 0d592018); bwUniCluster logout check is on the checklist |
| 3 TUI foundation | done: transport, generated protocol types, store, rows, footer, keymap, theme, config; the recorded conversation renders the same headless and live in tmux; a review's 12 bugs fixed in 2d37e333 |
| 4 One session | done: composer (6291c296), tool rows (702c9cc0), approvals and status line (2d37e333, 39e18945), palette, pickers and resize policy (4aef3b0a), /editor /keys /exit and the parity checklist (75627d84, docs/TUI_PARITY.md); a day with a real KIT model; a review's 15 bugs fixed in 33f58e19, a service re-review's in f57e6498 |
| 5 Agent view | built: grouped list on the terminal's alternate screen (esc gives the screen back exactly; GNU screen grows the footer instead), peek and reply, open, new session, stop, rename, pin, elsewhere, badge and title; live check `check_agents_tmux.py`; with the free KIT model, three real sessions side by side, an approval answered from the agent view and a new window finding them all ($0.00); a review's 13 bugs fixed in 7220059a. Decision: stopping takes ctrl+x twice (one key must not end a working session) and a third ctrl+x on the ended row removes it |
| 6 Lifecycle | done: dormant sessions -- a session no window shows sleeps after an hour unused (never while a turn, a question or a background shell of it runs) and wakes when opened or sent a message; what a stopped or crashed service ran comes back asleep in the next one; a window whose service vanished starts the next one and goes on, one told the service stops on purpose leaves; a hangup ends the window; an idle service of an older release is moved over (39a704da, b5a06be4). A review of phases 6 and 7 (fixed in 0e802f70 and 3475a86c): an update asks each session what a restart would cut and never cuts a turn, a question, a background shell or a subagent (a service from before dry runs is not asked at all); the roster follows every state, posture and name change at once, so a crash keeps them; nothing is imported from a project folder or the home; a restart no longer tells windows their sessions ended. Decision: no warm spare worker -- a worker imports in 0.07 s, and the slow part, the engine with its MCP servers, depends on the workspace and cannot be built ahead. Live checks `check_lifecycle_tmux.py`, `check_resize_tmux.py` |
| 7 Polish and shipping | done: `NO_COLOR` (mono theme), the first frame no longer waits for the terminal's theme (remembered, corrected in the background), budgets as tests (keystroke 0.5 ms, 1,020 rows written again in 0.6 s, launch to footer about 255 ms with the binary), `bun build --compile` binary that ignores a project's bunfig/.env, CI workflow, the manual's section and docs/TUI_CHECKLIST.md. A build counts as stale when a source, a patch, package.json or bun.lock is newer, or a source was deleted (0e802f70, 3475a86c). Not run here: the macOS and linux-arm64 binaries (CI only) |
