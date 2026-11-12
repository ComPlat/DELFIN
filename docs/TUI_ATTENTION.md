# Everything that needs you, in one place (plan)

**Status:** phases 1 and 2 built (2026-10-14); phase 3 is for later.
Decided with the user on 2026-10-14: → opens every item (a copy when the
session is still open elsewhere, a new session with the item as its
context when there is no session); new sessions wait for the user; the
note in a conversation is quiet, at the bottom right.

What building it added to the plan:
- **A resumed conversation is shown.** The window's replay comes from the
  agent service's event log, which has nothing of a session that ran in
  the dashboard or the REPL, nor of a copy. The worker now shows such a
  conversation as it opens: a notice and its last 20 turns, without the
  rows the harness wrote. Before, such a session opened on its banner alone.
- **`--fork-session` copies again.** `cli._open_session` read the copy's
  id from `session_store.fork_session` as a record; it is a string, the
  read raised, and every `--fork-session` (the window's copies included)
  resumed the original with a warning. Its test had a stand-in returning
  a record; a test against the store itself now pins it.
- **A session open elsewhere is looked up in the store too**, so → opens a
  copy of it rather than a new session; a session the store does not have
  yet is asked again at the next look.
- **A prepared fix** is the job monitor's assessment's `fix` (the
  assessment is always there, the fix only when one could be prepared).
- **Names:** a question with no session is named for its folder ("a
  question in benzene"); a session started about an item is named after
  it ("job 4711 failed").
- A `#` note before the first prompt leaves the brief for the prompt.
- Live check: `tui/test/e2e/check_attention_tmux.py`.

What the two reviews of it found, and what changed (each with a test):
- **A late answer went to the next session anywhere.** The engine drained
  answered items without its workspace, so an answer to a question about
  one project was handed to a session working in another, mid-turn, as
  "act on them now". It now drains with the folder the question was
  parked in (paths compared as paths); an item with no folder still goes
  to the next session.
- **The monitor's diagnosis session opened empty.** It was saved for the
  dashboard's display only, so a resumed engine had nothing of it. It is
  now saved as engine history too (its prompt marked as the monitor's),
  with the job's folder.
- **The brief goes with the first turn that reaches the model**, attached
  where a prompt becomes a turn: Esc right after enter, or a first turn
  the endpoint refuses before any output, leaves it for the next prompt.
  It is kept with the session's settings, so a session about an item that
  sleeps (or moves to a new service) before its first prompt still hands
  it over, and shows its context row once.
- **Items are masked** as a session's events are (`redact.Redactor`,
  moved out of the worker): a confirmation's command line can hold a
  token, and every window and watching agent process receives the items.
- **The inbox's other writers take its cross-process lock.** Noting a
  delivery, pruning and clearing rewrote the file under a thread lock
  only: a desktop notification noted by the job monitor could undo a
  dismissal made in the window, or drop another process's new event.
- **Smaller:** a session asleep in this service is no longer left out of
  the list (it may be open in a dashboard by now), and a session starting
  or ending here is noticed at the next poll; a saved title is read again
  when its record changes; the list carries a version (from the service's
  start on), so an update that crosses an `attention.list` answer on the
  socket is not overwritten; the history is shown before the session can
  run anything; a copy gets an id of its own (two in one second shared
  one); a brief ends with a line no brief holds by chance.
- **In the window:** a refused answer comes back with its peek, so enter
  tries it again instead of starting a session with it; ctrl+x pressed
  twice by habit dismisses one notice, not the one that moves up under
  it; `/rename` on an item renames nothing; a peek whose item was answered
  elsewhere closes, or says so while something is typed for it; of items
  that arrive together the one that waits names itself ("(+1 more)"); a
  newer open is never taken over by an older one's fallback; the notes on
  how an item was opened come after the session is in view.

The rest of this document is the plan as decided.

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
| The attention inbox | `~/.delfin/attention_inbox.jsonl` (`delfin.agent.attention`) | see below | blocking: `confirm_pending`, `question_pending`, `plan_pending`; notices: `run_failed`, `run_finished`, `budget_warning`, `security_alert`; `cluster_login` (a cluster tool found no login: the window opens its masked login field for the item's `cluster`; never answered with text, closed by the login) | `answer_item`, `dismiss_item`, `clear_all`; notices expire after 2-7 days by kind, at most 200 kept |
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
  › ◉ scan of benzene      question: which functional for the excited states?   dashboard · login02    2h
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
    ◇ ORCA builder         dashboard · login02
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
| → or enter | jump into it, wherever it came from (next section) |
| typing + enter | on a parked question: the answer (`attention.answer_item`), delivered to the next session in that workspace. On a session of this service: a reply, as today |
| ctrl+x | dismiss. A notice goes at once. A blocking item asks twice, because whoever waits on it gets no answer |

### Jumping into an item (→)

Every row opens in this window. The dashboard, the REPL and the window
all save sessions in one store (`~/.delfin/agent_sessions/`), so the
window can open a session that started anywhere else.

| The row | → does |
|---|---|
| A session of this service | opens it, as today |
| A session elsewhere that is closed (no longer open in the dashboard or terminal) | resumes that same session here: its whole history, the same id |
| A session still open elsewhere | opens a **copy** here: its history up to its last saved turn, under a new id. Opening the original twice would make two programs write one conversation; the session lock would then silently skip one side's saves. The row and the banner say "a copy of … (still open in the dashboard on login02)" |
| A parked question | opens its session, as above (resumed or a copy). You answer there as an ordinary message, and the inbox entry is closed, so the answer is never delivered a second time |
| A failed job the monitor diagnosed | resumes the monitor's diagnosis session, which holds the diagnosis and the prepared fix |
| Anything with no session behind it: a failed or finished job, a scheduled run, an alert, a budget notice whose session is gone | starts a **new session** in the item's folder, with the item as its opening context, shown as one row you can expand with ctrl+o. For a failed job the context is the job id, the SLURM state, the folder, the error lines found, the diagnosis and the prepared fix. A suggested prompt stands in the prompt box ("Job 4711 failed with OUT_OF_MEMORY. What went wrong, and what should we change?"). Nothing is sent, and nothing costs tokens, until you press enter. The dashboard's `/session handoff new` already works this way |

How the window does it, with what exists:
- **Resume:** `session.create {resume: id}`.
- **Copy:** when the session is still open elsewhere, the worker refuses the
  resume ("being written by another process; open a copy"), and the window
  asks for `session.create {resume: id, fork: true}`. A session on another
  login node is open when its presence record is fresh: the lock's pid
  cannot be checked from here.
- **New with context:** `session.create` gains `brief` (new, phase 2). The
  worker shows the brief as the session's first row and puts it in front of
  the first prompt. The suggested prompt is the window's own: it fills the
  prompt box.

Applying a prepared job fix (more memory per core, a longer time limit)
stays in the dashboard in the first version: the peek shows the fix and
says where to apply it. Doing it from the window needs the dashboard's
job backend in the service, which is phase 3.

### Seen from anywhere in the window

- **Quiet, at the bottom right:** the hints line ends with the count, right
  aligned and faint, with only the marks in colour: `◉ 3 · ✗ 1`. Today the
  hints line counts only this service's sessions, at its left.
- **A new item** that waits on you or failed names itself there for a few
  seconds, in the same quiet place, for example `◉ relax Fe complex asks`,
  then the count takes its place again. Nothing appears above the prompt,
  and nothing moves in the conversation.
- **The terminal title** carries the count: `(3) DELFIN · name`.
- Desktop notifications and webhooks already come from `emit_attention`
  itself; the window adds nothing there.

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
  texts and the keys above. The bottom-right count and note, and the title,
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
   - the groups, the rows, the peek texts, the keys, the bottom-right count and note, the title;
   - → on every row: resume, a copy, the diagnosis session, or a new session with the item as its context (`session.create` gains `brief`) and a suggested prompt;
   - render and key tests;
   - a live tmux check that writes a question and a finding into a private home's files and answers and dismisses them from the window;
   - the manual.
3. **Later, your call:**
   - presence records that say a dashboard session waits right now, instead of only after its question timed out;
   - applying a prepared job fix from the window;
   - the jobs list in the agent view (suggestion 1).

## Decisions

Taken with the user (2026-10-14):
- → jumps into every item. A session still open elsewhere opens as a
  **copy**: a take-over would need the dashboard to honour a hand-over
  request, which can come later.
- A new session for an item **waits for you**: its context is in, a
  suggested prompt is typed, nothing is sent.
- The note in a conversation is **quiet, at the bottom right**.

Standing as recommended, unless changed:
- NEEDS YOU and NOTICES live inside the agent view.
- Everything counts, with notices kept quiet.
- Answer and dismiss from the window.
