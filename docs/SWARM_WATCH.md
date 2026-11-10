# SWARM_WATCH — watching a swarm and one member session

The person asked (2026-11-06): *"in the dashboard or TUI, more intuitive, to
see what they do -- and to step into each one and see what it does."* Today a
swarm is visible as status rows only (the SW2a tab, `delfin-agent /swarm`),
and watching one member means `delfin-agent watch --session <name>` in a
terminal. This package (SWV) adds a read-only watch: the same table in the
dashboard and the TUI, with a click/Enter into a member's live feed.

Everything in this package is **read-only**: observing a swarm or a member
never starts, pauses, answers, prompts or types into any session. Steering
(answering an approval, pausing, resuming) stays in the hands of an SW3 lead
with the person's confirmation.

## The data path — one API, three surfaces

Dashboard, TUI and an SW3 lead all call the **same** read-only functions in
`delfin/agent/swarm/observe.py`:

- `swarm_overview(source)` -> one row per member: package, role, model,
  state, commits on its branch beyond the base, last tool call (name + short
  argument), last activity age, pending-approvals count.
- `member_feed(session, since)` -> the member's recent steps (tool calls with
  result status, refusals, commits), newest last.

`source` is either an SW1 swarm spec (a `{"packages": …}` dict, the shape
`delfin/agent/swarm/status.py` reads) or a **plain roster file** mapping
member name -> `{package, role}`, so an operator-run wave with no spec shows
too.

The step stream comes from the **existing** stores, never a second trace:
`observe.member_feed` reuses the audit-log step stream that
`delfin-agent watch` already reads (`audit_log.steps_since`), and the state
comes from the session-health record (`session_health.collect()` shape).
Liveness is derived from state files and records — **never** a process grep.

## In the dashboard

The **Swarm Watch** tab shows the overview table with a bounded live refresh
(one refresh per interval, no `squeue` per refresh, cooperative stop). Click
a member to open its detail pane: the live feed (follows like `watch`), its
commits and its pending approvals. Refusals are highlighted. Everything is
read-only: there is no script box into a member session and no approve
button of its own (the approval UI is reached by *link* only, and answering
an approval is an SW3-lead action guarded by the person's confirmation).

## In the TUI (`/swarm`)

The same overview table; **Enter** on a member opens its live feed, **Esc**
returns to the table. Keys are documented in the footer. Attaching is like
`tmux -r`: you can look, you cannot type.

## The read-only guarantee

Observing never changes a session. On the data side this is structural: the
watch functions only call existing read-only stores and render them; they
return JSON-shaped dicts and render HTML. On the UI side it is tested: the
rendered watch HTML contains no approve/answer button and no input box (see
`tests/test_swv_watch.py`). A member's feed is shown only to the user who
owns the session — the same ownership check the service uses; the watch
module only renders the data it is handed, and the service layer enforces
ownership server-side.

## Honest boundaries

- **"Waiting for a message."** The task's state vocabulary includes a member
  waiting on a person's message. The read-only health record distinguishes
  approval waits, model waits, pauses and ended turns, but it does not
  separately distinguish a plain "waiting for a person's reply" from
  "working". Rather than invent a liveness the store does not prove, the
  state falls back to **working** for that case (fail-closed). This is a
  documented limit, not a hidden assumption.
- **Secrets.** Nothing is read beyond what the audit log already shows; the
  redaction layer keeps masking. The watch never lifts redaction.

## Files

- `delfin/agent/swarm/observe.py` — the read-only data API.
- `delfin/dashboard/tab_swarm.py` — the Swarm Watch tab (`watched_rows_html`,
  `watched_feed_html`, `create_watch_tab`) and the existing SW2a manage tab.
- `tui/src/swarm/**`, `tui/test/swarm/**` — the TUI watch view.
- `tests/test_swv_*.py` — red controls (observe contract, read-only HTML).

See also: `docs/TUI_PARITY.md`, `docs/MAIN_PROTECTION.md`.
