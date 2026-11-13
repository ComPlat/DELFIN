# Real-page check of the agent tab

`scripts/page_check.py` starts the real dashboard against a scripted model,
opens it in headless Chromium and records what the agent tab shows. It is
the check to run when a change touches the dashboard's agent tab, its
sidebar or the agent service the dashboard talks to: the unit tests do not
render the page.

## Running it

The script runs on an interpreter with `playwright` and its Chromium
(`pip install playwright && playwright install chromium`). The dashboard
runs on the interpreter named by `--server-python` or
`DELFIN_PAGE_CHECK_SERVER_PYTHON`, which needs voila and DELFIN's dashboard
dependencies. The two may be different environments.

```bash
export DELFIN_PAGE_CHECK_SERVER_PYTHON=$MAMBA_ROOT_PREFIX/envs/delfin-ci/bin/python
nice -n 10 ionice -c3 $MAMBA_ROOT_PREFIX/envs/delfin-voila/bin/python \
    scripts/page_check.py --out /tmp/page-check --scenario composer
```

| Option | Meaning |
|---|---|
| `--out DIR` | output directory (required) |
| `--tree PATH` | DELFIN checkout the dashboard runs from; default: the checkout holding the script |
| `--scenario` | `composer` (default), `service-session`, `sidebar` or `keep` |
| `--timeout S` | hard deadline for the whole check; default 1800 |
| `--server-python P` | dashboard interpreter; default `$DELFIN_PAGE_CHECK_SERVER_PYTHON`, else the script's own |

Exit code 0: the scenario passed. 1: it ran and failed (the record says
where). 2: usage error.

## Scenarios

- **composer**: opens the tab whose title contains "Agent", waits for the
  composer (`.delfin-agent-input textarea`), selects the `ollama` provider
  and the scripted model, sends one prompt and waits until the reply has
  stopped changing for 5 s. Passes when the composer is visible and the
  scripted reply ("received ...") is on the page.
- **service-session**: starts a private agent service, creates a session
  there and runs one turn, then opens the dashboard's agent tab and
  records the sidebar's running-sessions list (`.delfin-running-item`).
  Passes when at least one running item is listed.
- **sidebar**: a private service with a front session, a background run (a
  folder under `.delfin/worktrees`) and an ended, saved session. The agent
  tab is photographed at 1366x768, 1920x1080, 1280x720, 2560x1440, 900x1400
  and 390x844 (`page_*.png`, `sidebar_*.png`, `controls_*.png`); at each size
  the record says whether the page itself scrolls, how far the inner area
  scrolls, and which controls of the sidebar or the session's control row
  overlap. Then the resume search: opened, typed into, Up, Esc, and Enter on
  a typed filter (`resume_*.png`): whether its list shows, which row is
  highlighted, what Enter opened and whether the field emptied.
- **keep**: the keep switches measured by process. Page A runs two
  sessions in their own folders and switches "Keep session" on for one;
  the page is left (pagehide, the close beacon) and closed. Recorded: when
  page A's kernel process ended (its pid, read from `/proc` by the kernel's
  connection file -- never the id), which sessions the agent service runs
  then (the kept one, keep true; the other not), the kept-session records
  and hand-over log, whether the server stopped with its last window (it is
  started again, as a person would), and what a new page B lists as running.
  B then switches the kept session off and is closed: the session must end
  with it. Run on a tree without the switches, the same steps are the
  control and the record says which switch was missing. The window grace is
  8 s here (`DELFIN_SESSION_GRACE_SECONDS`).

## Output

`DIR/record.json` (also printed): scenario, `tree_head`, `ok`, `error`,
`timings_s` (server ready, page loaded, composer visible, reply, wall),
the tab titles, the composer count, the reply text or the running items,
uncaught page errors, and a `cleanup` report. Screenshots (`opened.png`,
`answered.png` or `running.png`, `stuck.png` on a timeout) and
`logs/{voila,service,model}.log` sit beside it.

## Isolation and cleanup

- The dashboard's environment is built from nothing: no `VSCODE_*`,
  `BROWSER`, `TERM_PROGRAM` (they make the dashboard open a desktop
  browser tab), no `LD_LIBRARY_PATH`, no provider keys of the caller.
- It gets its own `HOME`, `DELFIN_SCRATCH_STATE` and workspace under
  `DIR`, which are removed afterwards; nothing is written to the caller's
  `~/.delfin`.
- `DELFIN_RUNTIME_DIR` is a 0700 directory directly under `/tmp`. A Unix
  socket path is limited to 108 bytes; a longer runtime path breaks the
  agent service's socket with no error on the page. The script refuses a
  path over 100 bytes.
- The agent tab is only built when a backend is configured: a placeholder
  `ANTHROPIC_API_KEY` makes it visible, and the page is switched to the
  scripted `ollama` endpoint (`tests/agent_terminal_probes/fake_model.py`,
  via `OLLAMA_HOST`).
- Waits are long (up to 420 s for the composer, 600 s for the first
  render) because the first render takes minutes on a loaded host.
- Cleanup runs in `finally` and on the deadline (`SIGALRM`, `SIGTERM`,
  `SIGINT`): the dashboard's process group is terminated, an agent service
  on the private socket is asked for `service.shutdown`, any remaining
  process whose environment carries this run's `DELFIN_RUNTIME_DIR`
  (kernels leave the dashboard's process group) is terminated and then
  killed, and the runtime directory is removed. No other process is
  touched.
- Do not attach a CDP debugger (`Debugger.enable`/`pause`) to this page;
  it hangs.

`tests/test_the_page_check_ends_what_it_started.py` checks the argument
parsing, the environment and the cleanup without a browser.
