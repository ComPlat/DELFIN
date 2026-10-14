# TUI rewrite: what is left for you

The rewrite ran its phases back to back (docs/TUI_PLAN.md). This is the one
list of what needs your hardware, your eyes or your decision, and what to
clean up on the VM. Each item says why it could not be closed here.

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
2. **bwUniCluster** (or another cluster).
   - Does the service survive a logout (systemd's `KillUserProcesses`)?
     The docs say "as long as the node allows".
   - Two real login nodes sharing the home: the window's note about the
     other node's sessions, and each node keeping its own service.
   - The cluster's own tmux and GNU screen versions.
3. **xterm.js terminals**: VS Code Remote-SSH and JupyterLab's terminal.
   The scrollback was checked in xterm.js headless (phase 0); a look in the
   real ones is still due.
4. **A look at the design** in your usual terminal: colours ("DELFIN blue,
   restrained"), the agent view, approval panels, light and dark themes,
   `NO_COLOR=1`.
5. **GNU screen with `altscreen on`**: set `DELFIN_TUI_AGENTS_SCREEN=alternate`
   and check the agent view gives the screen back; without it (screen's
   default) the view grows the footer and the prompt sits at the top after
   it closes.

## Decisions for you

6. **CI**: `.github/workflows/tui.yml` runs the window's and the service's
   tests, the live tmux checks and builds the binaries. GitHub keeps Actions
   off on forks until you enable them (Settings → Actions). The macOS
   runner labels may need updating.
7. **OpenTUI upstream**: we carry opentui PR #1665 and three startup fixes
   of our own (`tui/patches/`). Offering the three fixes upstream means
   posting on their repository: your call.
8. **DELFIN upstream issues found on the way** (not fixed here, they are
   upstream's code):
   - `repl.py` `_answer_question` shows a dict choice with `str(o)`;
   - the REPL's path completer ignores a trailing `/`;
   - under `DELFIN_SCRATCH_STATE`, session reports and session-index rows
     still go to the real `~/.delfin`;
   - `delfin-agent`'s start spends ~66 ms importing `api_client` for the
     security notices, and the scratch-state path imports ipywidgets
     (~90 ms): the window's launch would be faster without them.
9. **When to make the window the default** for `delfin-agent` (today it is
   `delfin-agent tui`, beside the REPL).

## Clean up on this VM

10. Bun in `~/.bun` (and its PATH lines in `~/.bashrc`), and the OpenTUI
    build cache `~/.cache/delfin-tui` (2.9 GB: Zig and three builds).
11. Packages installed for the terminal matrix: xterm, xvfb, imagemagick,
    xclip, screen (see `/var/log/apt/history.log`).
12. The SSH test key `~/.ssh/delfin_tui_localhost_test` and its
    `from="127.0.0.1,::1"` line in `~/.ssh/authorized_keys`.
13. In the real `~/.delfin`: two old test session records (0cc8db05…,
    32cfddf4…), and audit-log entries of the real-model test sessions
    (four of e974faa0, two of the agent-view check's) -- kept on purpose:
    an audit log is not edited. The tests' session reports and
    session-index rows were removed after each run.
14. `/home/jonas/project/.scratch-tui` (scripts, review repros, logs) and
    the test folders under `/tmp` (gone at the next reboot).
15. Processes: during the work, three orphaned test windows of an earlier
    review run, two stale `bun test` runs (one had spun at 100% CPU for
    hours: OpenTUI swallows the SIGTERM `timeout` sends) and one forgotten
    private test service were found and stopped; a faulty check helper
    (since fixed) also sent SIGHUP to up to three `bun src/main.tsx`
    processes it had not started -- by all signs the same kind of orphans.
    Nothing of yours was touched.

## Known limits

16. Same user, no cage: an agent's own process can connect to the service;
    it is told apart by process ancestry and may only watch.
17. Notifications are the terminal title and the hints line (◉): OpenTUI
    offers no bell or OSC 9 here.
18. `NO_COLOR`: everything visible is in the terminal's own colours; OpenTUI
    still writes a colour code for blank cells (invisible).
19. Parity with the REPL (docs/TUI_PARITY.md): marking jobs with arrows and
    cycling suggestions with Tab work differently.
20. No warm spare worker: a worker imports in 0.07 s; the slow part, the
    engine with its MCP servers, depends on the workspace.
21. The binary is ~108 MB (it carries the Bun runtime); linux-arm64 and
    macOS binaries are built in CI only, not run here.
