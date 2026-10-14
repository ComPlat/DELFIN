# Everything that needs you, in one place (plan)

**Status:** a plan, not built yet. The decisions at the end are open.

The agent view shows what this window's sessions need. DELFIN knows a lot
more that waits for you, and today it is spread over four places:
- the dashboard's `/attention` command;
- the job monitor's banners;
- desktop notifications that are gone once seen;
- sessions in other terminals.

The goal is one list in the window, fed by everything DELFIN already
records, that you can act on where it makes sense.

## What DELFIN records today

| Source | Where | Who writes it | What it holds | How an entry ends |
|---|---|---|---|---|
| This service's sessions | the agent service, in memory | every session the window runs | open questions, approvals, plans; turns that failed; sessions that ended | answered in the window |
| The attention inbox | `~/.delfin/attention_inbox.jsonl` (`delfin.agent.attention`) | see below | blocking: `confirm_pending`, `question_pending`, `plan_pending`; notices: `run_failed`, `run_finished`, `budget_warning`, `security_alert` | `answer_item`, `dismiss_item`, `clear_all`; notices expire after 2-7 days by kind, at most 200 kept |
| Job monitor findings | `~/.delfin/monitor_findings.jsonl` (`job_monitor`) | the monitor daemon | a failed job: state, error signatures, the diagnosis session, a summary, a prepared fix | the dashboard's banner and its fix button |
| Session presence | `~/.delfin/session_presence/` | every open session (dashboard, REPL, window) | where a session runs and what it is called | gone 15 min after its last heartbeat |

Who writes into the attention inbox today:
- **The self-modification guard** (dashboard, REPL): `confirm_pending`.
- **A question that timed out** with nobody there (dashboard, REPL): `question_pending`. The question is parked. An answer given later goes to the next session in that workspace (`drain_resolved` matches by workspace).
- **The job monitor**: `run_failed` for a failed watched job, with the diagnosis; `run_finished` for a finished one; also watch-list problems.
- **The scheduler**: a scheduled run failed or finished; a schedule switched itself off.
- **The engine**: `budget_warning`, when a run's budget is nearly spent or
  reached, a turn is stopped by the cost circuit-breaker, or a USD ceiling
  cannot be enforced (a model with no published price).
- **Security events**: `security_alert`, the first time the sandbox blocks
  a kind of action of a tool.

The file is shared between processes, guarded by a cross-process lock, and
rewritten atomically. On a cluster it is shared by every login node,
because the home is.

**A bug found on the way:** the scheduler daemon emits `watchpost_alert`,
which is not one of the inbox's kinds. `emit_attention` raises, the caller
swallows the error, and watchpost alerts never reach the inbox. Fix: add
the kind (urgent, a notice) or send them as `security_alert`.

TUI sessions themselves add nothing to the inbox except budget and security
notices. Their questions go through the service, and the window already
shows them.

## What the window will show

### In the agent view: NEEDS YOU and NOTICES

```
 ◆ AGENTS  4 sessions · 3 need you
  NEEDS YOU
  › ◉ scan of benzene      question: which functional for the excited states?   dashboard · uc3n991    2h
    ✗ job 4711 failed      OUT_OF_MEMORY · diagnosis ready · a fix is prepared    ~/calc/benzene        12m
    ◉ tui-p5-C             waiting for you: bash -> created_by_c.txt              here                   0s
  WORKING
    ⠙ opt of Fe(CO)5       running crest …                                                              3m
  DONE
    ✓ tui-p5-A             Item one                                                                    17s
  NOTICES
    · job 4702 finished    COMPLETED · ~/calc/water                                                     1h
    · budget nearly spent  long refactor · 80% of $2.00                                                 3h
  ELSEWHERE · watch and message
    ◇ ORCA builder         dashboard · uc3n991
```

- **NEEDS YOU** replaces NEEDS INPUT. It holds what blocks a session or
  needs a decision:
  - this service's sessions that wait for an answer;
  - parked questions and confirmations from any DELFIN session, dashboard and REPL included;
  - failed jobs and failed scheduled runs;
  - security alerts.

  Blocking items come first, then failures, newest first in each.
- **NOTICES**, just above ELSEWHERE, holds what happened and needs nothing:
  finished jobs and runs, budget warnings. It is quiet: faint, no badge,
  and it empties itself as the inbox lets notices expire.
- **A row** shows a mark (◉ waiting, ✗ failed, · notice), what it is about
  (the session's name, `job 4711`), one line on what it needs, where
  (`here`, `dashboard · node`, a folder), and its age.

### Keys

They are the agent view's, and work the same on every row:

| Key | On a row of NEEDS YOU or NOTICES |
|---|---|
| space | peek: the whole question and its choices, the job's diagnosis summary and the prepared fix, the alert's detail |
| → or enter | open it: a session of this service in this window; a job's diagnosis session (saved by the monitor), resumed in this window; a parked question's workspace as a new session that gets the answer on its first turn |
| typing + enter | on a parked question: the answer (`attention.answer_item`), delivered to the next session in that workspace. On a session of this service: a reply, as today |
| ctrl+x | dismiss. A notice goes at once. A blocking item asks twice, because whoever waits on it gets no answer |

Applying a prepared job fix (more memory per core, a longer time limit)
stays in the dashboard in the first version: the peek shows the fix and
says where to apply it. Doing it from the window needs the dashboard's
job backend in the service, which is phase 3.

### Seen from anywhere in the window

- **The hints line** counts what blocks across all sources:
  `◉ 3 waiting · ✗ 1 failed`. The badge is shown today only for this
  service's sessions.
- **The terminal title** carries the count: `(3) DELFIN · name`.
- **A new blocking item or failure** while you are in a conversation leaves
  one line in the footer for six seconds, for example: `a question from
  'scan of benzene' (dashboard) waits -- ← to answer`. Desktop notifications
  and webhooks already come from `emit_attention` itself; the window adds
  nothing there.

## How it works

### In the agent service

- **A watcher thread** checks the two files' size and modification time
  every 2 s. Polling, not inotify, because the home on a cluster is a
  network file system. When a file changed, it reads it with the module's
  own tolerant reader (`attention._load_events`, `job_monitor.load_findings`).
  It then builds the items and sends them to every window as one event,
  `attention.updated {items}`. The list is small: blocking items plus at
  most 200 notices.
- **Merging:**
  - an inbox record whose `session_id` is a session this service runs is
    folded into that session's row (no row twice);
  - a finding and its `run_failed` record (the same job id) become one item;
  - expired notices are left out, as the inbox's own read does.
- **Requests:**
  - `attention.list` returns the items, for a window's start and after a reconnect.
  - `attention.answer {id, text}` calls `attention.answer_item`.
  - `attention.dismiss {id, blocking?}` calls `attention.dismiss_item`.

  A watch-only client (an agent's own process) gets the list but cannot
  answer or dismiss, as with sessions today. `serviceVersion` becomes 3; a
  window shows the new groups only with a service that says 3.
- **Hosts:** an inbox record says nothing about where it came from.
  `emit_attention` gains a `host` field (`socket.gethostname()`). This is
  additive, and older records show no host.

### In the window

- The store gains the items. The agent view gains the two groups, the peek
  texts and the keys above. The hints line, the title and the footer note
  read the same items.
- Nothing is polled by the window. It gets `attention.updated` while it is
  connected and `attention.list` after a reconnect.

## Edge cases

- **A parked question whose session ended:** the answer waits in the inbox
  and goes to the next session started in that workspace. The row says
  "answered: given to the next session in ~/calc/benzene".
- **Old blocking items** (a dashboard closed hours ago) show their age and
  can be dismissed. Nothing blocking is dismissed by age, which is the
  inbox's own rule.
- **Several windows** see the same list. A dismiss in one removes the row
  from all of them within 2 s.
- **Login nodes sharing a home** see each other's items, with the host.
  Answering one from another node is safe: the answer goes into the
  shared file.
- **A torn or locked file:** the window keeps the last list it had. The
  reader skips torn lines, as the module's own does.

## Phases

1. **Service:**
   - the watcher, the items, the merging, the three requests and the event;
   - the `host` field and the watchpost fix;
   - unit tests with a fake inbox and fake findings files, including several processes writing at once;
   - docs/TUI_PROTOCOL.md updated.
2. **Window:**
   - the groups, the rows, the peek texts, the keys, the hints-line count, the title and the footer note;
   - render and key tests;
   - a live tmux check that writes a question and a finding into a private home's files and answers and dismisses them from the window;
   - the manual.
3. **Later, your call:**
   - presence records that say a dashboard session waits right now, instead of only after its question timed out;
   - applying a prepared job fix from the window;
   - the jobs list in the agent view (suggestion 1).

## Decisions for you

1. **Where:** NEEDS YOU and NOTICES inside the agent view (recommended,
   one place for oversight), or a screen of its own.
2. **What counts:** everything, with notices kept quiet (recommended), or
   only what blocks plus failures.
3. **Acting:** answer and dismiss from the window (recommended), or
   show only and point to the dashboard.
4. **New items while you work:** the footer note plus the count
   (recommended), or the count only.
