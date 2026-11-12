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
| `--scenario` | `composer` (default) or `service-session` |
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
