# TUI rewrite: what is left for you

The rewrite ran its phases back to back (docs/TUI_PLAN.md). This is the one
list of what needs your hardware, your eyes or your decision. Each item says
why it could not be closed here.

## Needs your hardware or your eyes

1. **macOS** (Terminal.app, iTerm2; Apple silicon and Intel).
   - Run `delfin-agent tui` there: the window, the agent view (alternate
     screen), approvals, a resize, ctrl+d and a second window.
   - The patched native library is built for Linux only
     (`tui/scripts/build-opentui-native.sh` has no Darwin target yet); the
     macOS binaries in CI use OpenTUI's prebuilt library, which lacks the
     footer-growth fix (opentui#1630). That matters for xterm.js terminals
     (VS Code) on a Mac.
   - The watch-only rule for an agent's own processes uses `ps -wwE` on
     macOS: check that a window can act and that `delfin-agent service
     stop` typed inside an agent's shell is refused.
   - Full screen: Terminal.app takes no OSC 52 (the window says so after a
     drag); a selection of the terminal's own needs fn (Terminal.app) or
     option (iTerm2) held while dragging.
2. **bwUniCluster** (or another cluster).
   - Does the service survive a logout (systemd's `KillUserProcesses`)?
     The docs say "as long as the node allows".
   - Two real login nodes sharing the home: the window's note about the
     other node's sessions, and each node keeping its own service.
   - The cluster's own tmux and GNU screen versions.
   - The binary there: it needs glibc 2.17 at most (checked: Bun's runtime
     and the native library), but whether the node's kernel suits Bun is
     only known by running it.
3. **xterm.js terminals**: VS Code Remote-SSH and JupyterLab's terminal.
   The scrollback was checked in xterm.js headless (phase 0); a look in the
   real ones is still due -- full screen too: the wheel, a drag that copies
   (OSC 52; VS Code may ask to allow the clipboard), shift+drag.
4. **A look at the design** in your usual terminal: colours ("DELFIN blue,
   restrained"), the agent view, approval panels, light and dark themes,
   `NO_COLOR=1`.
5. **GNU screen**: the window is inline there by default (screen's
   alternate screen is off). With `altscreen on`, `"screen": "full"` in
   tui.json should work like elsewhere; inline, set
   `DELFIN_TUI_AGENTS_SCREEN=alternate` and check the agent view gives the
   screen back; without it the view grows the footer.

## Decisions for you

6. **CI**: `.github/workflows/tui.yml` runs the window's and the service's
   tests, the live tmux checks and builds the binaries; it passes on the
   fork, where each binary's smoke test starts it (macOS and linux-arm64
   included). The macOS jobs run on `macos-14` (Apple silicon) and
   `macos-15-intel`; GitHub retires runner images now and then, so a label
   may need updating.
7. **OpenTUI upstream**: we carry opentui PR #1665 and three startup fixes
   of our own (`tui/patches/`). Offering the three fixes upstream means
   posting on their repository: your call.
8. **DELFIN upstream issues found on the way** (not fixed here, they are
   upstream's code):
   - `repl.py` `_answer_question` shows a dict choice with `str(o)`;
   - the REPL's path completer ignores a trailing `/`;
   - under `DELFIN_SCRATCH_STATE`, session reports and session-index rows
     still go to the real `~/.delfin`;
   - the engine starts its own MCP servers as `python -m
     delfin.tools.mcp_server` and `python -m delfin.ops_server`
     (`agent/mcp_client.py`) in the project folder, and `-m` imports from
     the working directory first: a `json.py` in the project folder (one
     the agent wrote earlier, say) is imported by them when they start --
     by default outside any containment (`agent.mcp_isolation` is off).
     The REPL does the same. `-P` in their args (Python 3.11+) closes it;
     the service and its workers now start in the runtime directory for
     the same reason (3475a86c);
   - tests that build an agent engine leave its MCP servers running when
     they end: one full run left 20 of them (from `test_mcp_http`,
     `test_the_rules_can_change_while_the_agent_works`,
     `test_two_hundred_skills_through_the_public_path`,
     `test_u11_whole_request_is_counted_once`), plus multiprocessing
     helpers of `test_a_completion_notice_is_not_destroyed_before_it_arrives`;
   - `delfin-agent`'s start spends ~66 ms importing `api_client` for the
     security notices, and the scratch-state path imports ipywidgets
     (~90 ms): the window's launch would be faster without them.
9. **When to make the window the default** for `delfin-agent` (today it is
   `delfin-agent tui`, beside the REPL).

## Known limits

10. Same user, no cage: an agent's own process can connect to the service;
    it is told apart by process ancestry and may only watch.
11. Notifications are the terminal title and the hints line (◉): OpenTUI
    offers no bell or OSC 9 here.
12. `NO_COLOR`: everything visible is in the terminal's own colours; OpenTUI
    still writes a colour code for blank cells (invisible).
13. Parity with the REPL (docs/TUI_PARITY.md): marking jobs with arrows and
    cycling suggestions with Tab work differently.
14. No warm spare worker: a worker imports in 0.07 s; the slow part, the
    engine with its MCP servers, depends on the workspace.
15. The binary is ~108 MB (it carries the Bun runtime); the linux-arm64
    and macOS binaries are built and smoke-tested in CI only, not used here.
16. An update that changes the protocol between window and service: an
    older service that still runs sessions is never stopped by the new
    window -- it cannot be asked what a stop would cut (a turn, a
    background shell). The window names it; `delfin-agent service stop`
    ends it when you choose (its sessions save and come back asleep). An
    older service of the same protocol is asked, session by session, and
    moved over by itself once nothing would be cut; one from before this
    branch's last review cannot be asked (it would end the session it is
    asked about) and is moved over only when no session of it is live.
17. Code in an answer is coloured only for the grammars OpenTUI ships
    (JavaScript, TypeScript, Markdown, Zig): Python, shell and ORCA input
    are drawn plain. Adding them means shipping their tree-sitter grammars
    (a `.wasm` and its highlight queries each) with the window.
18. ctrl+o differs from the REPL on purpose: it shows or hides every tool
    call's details (the window draws a call as its headline and what came
    back); the last result in full is in the palette.
19. What needs you (docs/TUI_ATTENTION.md) comes from the attention inbox
    and the job monitor's findings, read every 2 s. A dashboard session
    that waits for an answer right now shows only once its question timed
    out and was parked (phase 3: presence that says it waits). A prepared
    job fix is applied in the dashboard (`/fix`), not from the window.
20. A resumed conversation is shown from the agent service's event log
    when that has its turns; when it has none (the session ran in the
    dashboard or the REPL, or it is a copy), from the store: its last 20
    turns. A session that ran in a window and then went on in the
    dashboard shows the window's turns only -- the model has all of them.
