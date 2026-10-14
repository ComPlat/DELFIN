# TUI rewrite: `delfin-agent tui`

Status: agreed scope, 2026-10-13. Branch `tui-rewrite` of
`the16thpythonist/DELFIN`; fork-only for now.

## Goal

A new terminal UI for the DELFIN agent that looks good and feels smooth.
It combines OpenCode's way of building a TUI (TypeScript, OpenTUI, SolidJS)
with Claude Code's way of working with sessions: one window, an agent view on
`←`, and sessions that keep running after the window closes.

## Decisions

| Topic | Decision |
|---|---|
| Stack | TypeScript on Bun, OpenTUI + SolidJS (OpenCode's stack). |
| Screen | Full screen (since 2026-10-14, after a first try): the window takes the terminal's alternate screen, the conversation scrolls inside it (wheel, PgUp/PgDn), the composer stays at the bottom, a drag selects and copies (OSC 52), and the shell's lines come back as they were. Inline stays as an option (`--inline`, tui.json `"screen"`, the default under GNU screen): finished output goes into the terminal's own scrollback; only a footer (composer, status, running tools, approvals, menus) is live. OpenTUI `screenMode: "alternate-screen"` and `"split-footer"`; blueprints: OpenCode's TUI and `opencode mini` (`anomalyco/opencode`, `packages/tui/src/mini/`). |
| Look and feel | Must look good and feel smooth: no flicker, no torn lines, responsive while the agent works. "DELFIN blue, restrained": neutral text, with the logo's steel blue (`#6290C0`) as the single accent; green, amber and red only for state; truecolor with a fallback; the light or dark variant follows the terminal. |
| Resize | Full screen: the rows re-wrap in place; nothing is written again. Inline: re-wrap by replay, like `opencode mini` -- after a resize settles, the screen and the terminal's scrollback are cleared and the session is reprinted at the new width (capped at a configurable number of turns). |
| Sessions | Every session is part of the agent view from the moment it starts; there is no "backgrounding" step. Sessions keep running when the window closes, and the next window re-attaches. A session left idle (waiting for your next message) and unattached sleeps after about an hour and wakes when you reply. |
| Agent view | `←` on an empty prompt. A Claude Code-style list grouped into Needs input / Working / Done. It covers every session on this machine across projects, plus other open DELFIN sessions (dashboard, `delfin-agent chat`) as read-only rows. It is drawn on the terminal's alternate screen (decided in phase 5): closing it gives the screen back exactly as it was, where growing the footer pushed the session's rows into history. Under GNU screen, whose alternate screen is off by default, it grows the footer. |
| Switching | A labelled divider, then a replay of the selected session's last few turns into the scrollback. Nothing is erased. |
| Keys | Today's keys (esc, shift+tab, ctrl+g, ctrl+o, ctrl+t) plus a ctrl+p command palette. Every binding lives in one configurable registry. |
| Edits | Sessions edit in place (no git worktrees). |
| Protocol | Our own JSON-RPC 2.0 (newline-delimited JSON) over a Unix socket, with events modelled on OpenCode v2. An ACP adapter for editors comes later. |
| Rollout | A new command, `delfin-agent tui`. `delfin-agent chat`, `-p`, pipes and JSON output stay as they are. |
| Platforms | Linux workstations, HPC login nodes (SSH, tmux), macOS. Not Windows. |

## Architecture

```
 terminal ── delfin-agent tui        Bun binary: OpenTUI + SolidJS, a thin client
                 │  JSON-RPC over /tmp/delfin-$UID/agent.sock
                 ▼
            agent service            Python, one per host, started on demand, detached
                 │  pipes
         ┌───────┼────────┐
       worker  worker   worker       one Python process per session, wraps AgentEngine
                 │
       event log per session         replay, resume after a crash, resume from another host
```

- **Launcher** (`delfin-agent tui`, Python): connects to the service, or
  starts it (detached; a single instance per host via `flock` on the run
  directory). It then hands the terminal to the TUI binary along with the
  socket path.
- **Service**:
  - owns the sessions, passes requests to their workers, and fans events out
    to every attached window;
  - holds permission requests until someone answers them; the first answer
    wins, and the other windows are told;
  - exits after N idle minutes;
  - writes `~/.delfin/run/<host>.json`, so a window on another login node can
    say where the sessions run.

  This copies Claude Code's on-demand supervisor (`claude daemon`), which has
  a lock file, a control socket with an auth token, one process per session,
  and a roster file used to reconnect. Like it, we accept a client only
  through a socket in a 0700 directory, after checking the peer's UID and a
  token kept in the run file.

  We also copy two of its behaviours:
  - after an upgrade, the service restarts onto the new version and moves
    idle sessions over; sessions that are working, waiting on you or attached
    are not interrupted;
  - if the service dies, its workers wait briefly for a new one, then stop;
    the next attach or reply resumes them from the event log.
- **Worker**, one per session:
  - runs `AgentEngine.stream_response` and turns its callbacks into events;
  - blocks on a permission callback until the answer arrives;
  - is a separate process because DELFIN changes process-wide state
    (`os.chdir`, `os.environ`), and so that a stuck session can be killed on
    its own.
- **Event log**: append-only, one per session. It feeds the replay on attach
  (the last N turns), and lets a session resume from any host after the
  service died.
- **TUI**, structured like OpenCode's: transport → events batched every
  10–16 ms → normalised Solid store → views. Finished blocks are committed to
  the scrollback; everything still changing stays in the footer.

## Protocol sketch (to be firmed up in milestone 1)

**Requests (TUI → service)**

- `initialize` (version handshake)
- `session.list`, `session.create {cwd, model?, mode?}`
- `session.attach {id, replay: {turns}}`, `session.detach`
- `session.prompt {id, text, queue?}`, `session.steer {id, text}`,
  `session.cancel {id}`
- `session.setMode`, `session.setModel`, `session.setEffort`
- `session.command {id, line}`: runs one of the existing Python slash
  commands and returns text
- `session.stop`, `session.delete`
- `permission.reply {requestId, optionId}`

**Events (service → TUI)**

- `session.created`, `session.updated`
- `session.state {running | idle | needs_input | sleeping | error}`
- `message.started`, `message.ended`
- `text.delta`, `reasoning.delta`
- `tool.started`, `tool.ended {callId, name, input, result, meta}`
- `permission.asked`, `permission.resolved`
- `question.asked`, `question.resolved`
- `notice`
- `usage {tokens, cost}`
- `presence.updated` (other DELFIN sessions)

The engine reports tool results by tool name, without a call ID. Until the
engine gets real IDs, the worker assigns them first-in-first-out per tool
name, as `repl.run_turn` does today.

## What is reused from the current Python TUI

- **Engine callbacks:** `repl.run_turn` / `RenderItem` become worker events.
- **Approval options:** `terminal_confirm.options_for` provides them
  (y/n/A/k/R/e/a/?) and they travel with `permission.asked`, so both UIs offer
  the same choices.
- **Slash commands:** `repl_commands.dispatch` runs the text-only ones in
  Python (`session.command`). The interactive ones (model, permissions,
  sessions) get native pickers in the TUI.
- **Launch checks:** the `launch_guard` checks and the banner facts are shown
  when a session is created.
- **Session data:** `session_store` and `session_presence` provide saved
  sessions and other DELFIN sessions.

## Milestones

The phased implementation plan is in [TUI_PLAN.md](TUI_PLAN.md).

## Defaults chosen (change any of them)

- **Session rows** show a fixed summary (title or last prompt, plus current
  activity), not a model-written one, so they cost no tokens.
- **Approvals** are answered inside the session (→), as in Claude Code, not
  from the agent view's peek panel.
- **New sessions** started from the agent view use the window's folder and
  the default approval mode (plan) unless one is configured.
- **Code layout**: the TypeScript app goes in `tui/` (a Bun workspace) and the
  Python service in `delfin/agent/service/`.
- **Tests**: OpenTUI's headless test renderer with text-frame snapshots, a
  fake event stream for UI tests, and pytest for the service.

## Phase 0 findings (2026-10-13)

**Decision: GO.** We build on OpenTUI 0.6.4 + SolidJS in split-footer mode.
The native library is built from source with upstream PR #1665 until
upstream merges it, and three startup fixes of our own are carried as a bun
patch.

### How the scrollback is checked

The spike (`tui/spike/inline.tsx --auto`) plays a scripted session through
the real scrollback core (`tui/src/scrollback/`):
- four agent turns with tool rows, a highlighted diff and a streamed
  markdown answer;
- approval panels of random height, which make the footer grow and shrink
  at random;
- the agent view, which grows the footer close to full height;
- a session switch, a realistic 30-turn replay and a 1,000-line replay.

Every row it commits carries a unique tag (`{t00001}`, …), and its report
lists every tag the scrollback must hold, in order.

Each harness gets the terminal's real history out of the terminal:
- tmux: `capture-pane` on a private server;
- xterm.js: its buffer, live;
- VTE: its own text API;
- XTerm: `print-everything`;
- GNU screen: `hardcopy -h`.

`tui/spike/scrollback_check.py` then requires:
- exactly the expected tags: nothing missing, nothing twice, nothing out of
  order;
- the shell lines printed before the spike still intact above it;
- after a resize, exactly one replay;
- no spike errors, nothing on stderr, and a clean exit.

`tui/spike/run_checks.sh` runs the whole matrix.

Two things the first, marker-based harness could not see:
- Brackets are wrong for tags: markdown reads `[t00001]` as a link and the
  renderer conceals the brackets.
- tmux applies `history-limit` only to panes created after it is set. The
  limit must be set in the same command that creates the session.

### Results (patched build, final)

| Terminal | Cases | Result |
|---|---|---|
| xterm.js 6.0, live (VS Code, JupyterLab, browser shells) | replies delayed 0, 150, 300, 1000, 2000 ms × start rows 1, 9, 21, 32, 35 × 2 fuzz seeds (50 runs) | all rows exactly once |
| tmux 3.4 | start rows 1, 21 and the bottom row; resized 110→80 and 110→140 mid-run; over SSH from the bottom row; keys | all pass; resize gives exactly one replay |
| VTE 0.76 (GNOME Terminal) | start rows 1, 21, 35; 80×24 from the bottom | all pass |
| XTerm 390 | start rows 1, 21, 35; 80×24 from the bottom | all pass |
| GNU screen 4.9.1 | start rows 1 and 35 | all pass |
| GNU screen 4.9.1, resized | 110→80 mid-run | 52 rows twice: screen has no ED3, so a replay cannot erase its scrollback (see phase 4) |

Negative control: with the prebuilt native library, xterm.js loses 49 and 54
rows (two seeds) when the footer grows. That is opentui#1630, and the checker
names every missing row.

### What we carry on top of OpenTUI 0.6.4

- **Native: upstream PR #1665** (`tui/patches/opentui-v0.6.4-pr1665-footer-growth-scrollback.patch`).
  - It grows the footer by saving the cursor, moving to the footer's last row,
    writing line feeds and restoring the cursor (DECSC, CUP, n LF, DECRC),
    instead of `CSI n S`. xterm.js deletes the rows that `CSI n S` scrolls
    off, while line feeds push them into history.
  - Its JS half (`scrollIntoHistory`) is in the bun patch.
  - `tui/scripts/build-opentui-native.sh` builds it at a pinned commit with a
    pinned Zig, runs its tests and installs it. `--check` verifies the
    installed library. Re-run it after every `bun install`.
- **JS: the bun patch** (`tui/patches/@opentui%2Fcore@0.6.4.patch`, applied
  by `bun install`), which carries three startup fixes:
  1. **Starting near the bottom row overwrote the line above.**
     - Setup reserves the footer's rows with line feeds after the
       cursor-position query was sent. On the bottom rows those line feeds
       scroll the screen.
     - The startup seed then took the pre-scroll row, and the first commit
       wrote over the last line of output. The usual case for this is
       `delfin-agent tui` typed into a terminal that has already scrolled.
     - Prebuilt and patched native libraries both had the bug, in every
       terminal.
     - Fix: seed from the row as it is after the reservation. If that row is
       under the footer, push it into history first.
  2. **A fast terminal's reply arrived before setup finished.**
     - That happens with XTerm, which answers every startup query. Fix 1's
       scroll then has to wait until setup completes. It is written at once
       after setup, before the first frame.
  3. **A slow reply (over 1.5 s) overwrote the screen from row 1.**
     - The fallback seed was row 1. Now a missing reply counts as "cursor on
       the bottom row".
     - At worst that leaves blank rows above the footer; it never overwrites.
     - The reply timeout itself is raised from 120 ms to 1.5 s
       (`OTUI_SPLIT_CPR_TIMEOUT_MS`): an earlier review showed that 300 ms of
       latency overwrote shell lines.
- These three fixes should go upstream as a PR once the user has seen them
  (left for the user: it means posting on a third-party repository).

### Measurements

- A 1,000-line replay takes 400–470 ms, and a realistic 30-turn replay with
  markdown, diffs and highlighting about 0.85–0.9 s. Both were measured in
  tmux on a loaded machine (load average 15–20).
- Idle CPU use is zero.
- The exit is clean: the footer disappears and the transcript stays.

### Lessons the implementation must follow

1. **Render on demand, and run live while scrollback output is queued.**
   Demand-only rendering stalls on a burst of commits. Always-on rendering
   costs CPU and never lets `renderer.idle()` resolve.
   (`tui/src/scrollback/queue.ts`)
2. **`resize` also fires when the footer changes height.** Only a change of
   terminal width may trigger the clear-and-replay; anything else would wipe
   scrollback for nothing. (`resize.ts`)
3. **Commit whole items, and replay as one queue task.**
   - Record each item when it is written, not when it is queued.
   - Then a replay can neither miss an item nor write one twice, and
     nothing lands between the wipe and the reprint. (`transcript.ts`)
4. **Solid's JSX compiler collapses runs of spaces in literal text.**
   Significant spacing must be written as string expressions.
5. **Never exit the process from `onDestroy`.** OpenTUI flushes queued
   scrollback after destroy. Exiting early drops it.
6. **The diff renderer rejects hunks whose line counts are wrong.** Validate
   diffs with jsdiff's `parsePatch`, and fall back to plain text.
7. **Markdown heading styles need one key per level**
   (`markup.heading.1` … `.6`).
8. **The live tail in the footer must be rendered as markdown.** Raw source
   such as `| --- |` looks broken.
9. **OpenTUI bundles no Python tree-sitter grammar.** Phase 4 registers
   Python, plus ORCA/shell where available.
10. **Colours over SSH:** truecolour is only detected when `COLORTERM`
    reaches the remote host. Otherwise output falls back to 256 colours.
11. **`OPENTUI_GRAPHICS=0`.** The kitty-graphics probe retitles tmux panes,
    and this TUI draws no images.
12. **GNU screen cannot erase its scrollback (no ED3).**
    - On a resize under screen, print a divider ("resized to N columns —
      the rows above keep their old width") instead of replaying.
    - Detect screen by `STY`. TERM doesn't tell it apart, because tmux sets
      `screen*` too.
13. **A block on its way from the live tail into the rows is shown by one
    of them in every frame** (found 2026-10-14; `handoff.test.tsx`).
    - A new markdown element draws nothing (or, `streaming`, its raw text)
      until its highlighting lands a few ms later. Keep the finished block
      in the footer (`held`), drawn by the element that drew it while it
      streamed, until its rows are drawn; on the full screen lay the new
      element out at no height until it is highlighted.
    - Inline, the rows and the footer reach the terminal on two channels:
      commit a held block's rows at a frame's start (a frame callback), so
      the footer has let go of it before that frame is painted.
    - Don't draw the live block `streaming`: OpenTUI's preview of each
      delta is the raw text ("## ") until it is highlighted again.
    - OpenTUI highlights in a frame's pass over what is visible: a block
      culled (below the view) or hidden (behind the agent view) is never
      highlighted, and blank when it comes into view. Start it by hand
      (`whenHighlighted`), after passing the element's style on.

## Risks

- **OpenTUI's inline mode is young.**
  - Its native path dates from 2026-03.
  - We carry upstream PR #1665 for bug #1630, plus three startup fixes (see
    "Phase 0 findings"). The native part has to be rebuilt for each OpenTUI
    upgrade until upstream has both.
  - `tui/spike/run_checks.sh` is the test that tells whether an upgrade
    still works.
- **HPC login nodes** have process reapers, several login nodes and a shared
  NFS `$HOME`. The event log and resumable sessions mitigate this. Sessions
  survive logout exactly as well as tmux does on that cluster.
- **tmux has no kitty keyboard protocol.** Keys stay simple: letters, Enter,
  Esc and Ctrl+letter.
- **Upstream churn.** Upstream changes the engine and the REPL daily, so keep
  the Python changes small and merge `upstream/main` often.
- **Shipping.**
  - Binaries are about 90 MB per platform.
  - Bun extracts its native library to a temporary directory, which fails on
    a `noexec` `/tmp`; point `TMPDIR` at a cache directory instead.
