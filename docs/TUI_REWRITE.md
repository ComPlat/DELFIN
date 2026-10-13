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
| Screen | Inline. Finished output goes into the terminal's own scrollback; only a footer (composer, status, running tools, approvals, menus) is live. OpenTUI `screenMode: "split-footer"`; blueprint: `opencode mini` (`anomalyco/opencode`, `packages/tui/src/mini/`). |
| Look and feel | Must look good and feel smooth: no flicker, no torn lines, responsive while the agent works. "DELFIN blue, restrained": neutral text, with the logo's steel blue (`#6290C0`) as the single accent; green, amber and red only for state; truecolor with a fallback; the light or dark variant follows the terminal. |
| Resize | Re-wrap by replay, like `opencode mini`. After a resize settles, the screen and the terminal's scrollback are cleared and the session is reprinted at the new width (capped at a configurable number of turns). |
| Sessions | Every session is part of the agent view from the moment it starts; there is no "backgrounding" step. Sessions keep running when the window closes, and the next window re-attaches. A session left idle (waiting for your next message) and unattached sleeps after about an hour and wakes when you reply. |
| Agent view | `←` on an empty prompt. A Claude Code-style list grouped into Needs input / Working / Done. It covers every session on this machine across projects, plus other open DELFIN sessions (dashboard, `delfin-agent chat`) as read-only rows. |
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

**Decision: GO.** We build on OpenTUI 0.6.4 + SolidJS in split-footer mode,
using a native library built from source with one patch until upstream fixes
issue #1630.

The spike (`tui/spike/inline.tsx`) runs a scripted session: two agent turns
with tool rows, a highlighted diff, an approval panel, a streamed markdown
answer, the agent view (the footer grows close to full height), a session
switch, and a 1,000-line replay. Afterwards the scrollback was captured from
each terminal and checked: every expected row present exactly once and in
order.

| Terminal | Prebuilt OpenTUI | Patched OpenTUI |
|---|---|---|
| tmux 3.4 (also over SSH) | all rows kept | all rows kept |
| xterm.js 6.0 (VS Code, JupyterLab, browser shells) | 14 rows lost when the footer grows | all rows kept |
| VTE 0.76 (GNOME Terminal) | all rows kept | all rows kept |
| XTerm 390 | all rows kept | all rows kept |

- **The patch:** `tui/patches/opentui-v0.6.4-footer-growth-linefeeds.patch`
  grows the footer with line feeds instead of `CSI n S`.
  `tui/scripts/build-opentui-native.sh` builds and installs it. Re-run the
  script after every `bun install`.
- **Speed:** a 1,000-line replay takes 400–480 ms; idle CPU use is zero.
- **Resize:** shrinking the terminal from 110 to 80 or 90 columns mid-run
  re-wrapped the whole session correctly, including over SSH.
- **Exit:** the exit is clean. The footer disappears and the transcript stays.
- **Keys** behave correctly in tmux: ←, shift+tab, esc, and ← with text typed.

Lessons the implementation must follow:

1. **Render on demand, and run live while scrollback output is queued.**
   Demand-only rendering stalls on a burst of commits. Always-on rendering
   costs CPU and never lets `renderer.idle()` resolve.
2. **`resize` also fires when the footer changes height.** Only a change of
   terminal width may trigger the clear-and-replay. Anything else would wipe
   scrollback for nothing.
3. **Solid's JSX compiler collapses runs of spaces in literal text.**
   Significant spacing must be written as string expressions.
4. **Never exit the process from `onDestroy`.** OpenTUI flushes queued
   scrollback after destroy. Exiting early drops it.
5. **The diff renderer rejects hunks whose line counts are wrong.** Validate
   diffs, and fall back to plain text.
6. **Markdown heading styles need one key per level**
   (`markup.heading.1` … `.6`).
7. **The live tail in the footer must be rendered as markdown.** Raw source
   such as `| --- |` looks broken.
8. **OpenTUI bundles no Python tree-sitter grammar.** Phase 4 registers
   Python, plus ORCA/shell where available.
9. **Colours over SSH:** truecolour is only detected when `COLORTERM` reaches
   the remote host. Otherwise output falls back to 256 colours.

## Risks

- **OpenTUI's inline mode is young.** Its native path dates from 2026-03, and
  bug #1630 is patched locally (see "Phase 0 findings"). The patch has to be
  rebuilt for each OpenTUI upgrade until upstream fixes the bug.
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
