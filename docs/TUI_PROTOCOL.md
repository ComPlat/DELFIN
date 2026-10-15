# TUI protocol

How the terminal UI (`delfin-agent tui`, TypeScript) talks to an agent
session (Python). The single definition of every message is
`delfin/agent/service/protocol.py`. `python -m delfin.agent.service.protocol`
prints it as JSON Schema (requests, results and events), and the TUI's
types are generated from that. This page explains what the messages mean.

## Transport

- **Wire format.** JSON-RPC 2.0, one JSON object per line, UTF-8. JSON
  escapes every control character, so a newline never occurs inside a
  message. A line that is not a JSON object is reported with a parse-error
  response and skipped; it never ends the session.
- **Text that is not UTF-8.** A file name that is not valid UTF-8 reaches
  Python with lone surrogates in it. It travels as the JSON escape for the
  same code unit (`\udce9`), so the line stays valid JSON and valid UTF-8.
- **Directions.** The client sends requests (`id`, `method`, `params`). Each
  one gets a `result` or an `error`.
- **Events.** The session reports what happens as notifications with
  `method: "event"`. Their `params` carry:
  - `type`, the event name;
  - `session`, the session id;
  - `ts`, unix time;
  - `seq`, on durable events only (see Events): increasing per session and
    never reused, also across a restart of the worker;
  - the event's own fields.
- **Order.** Events arrive in the order things happened. An approval asked
  from a tool thread, or a message queued from the window, is placed behind
  the output the turn had already produced, so a question never arrives
  ahead of the tool call it is about.
- **Two hops.** A window talks to the agent service over a Unix socket
  (see "The service"); the service talks to one worker process per session
  over the worker's stdin and stdout. Session requests and events are the
  same on both hops: the service adds `sessionId` to route a request and
  passes events through unchanged.

## Lifecycle

1. The worker starts and emits `worker.ready`, `session.info` and
   `session.state {idle}`.
2. The client sends `initialize {protocolVersion}`. Nothing else is accepted
   before it but `session.shutdown`, so a session still starting can be ended. A version mismatch is refused with code -32001.
   - The result carries the session info, the current state, and every
     approval still waiting.
   - A window that attaches later sends `initialize` again, and then
     `session.replay` to draw what it missed.
3. Prompts, answers and mode changes follow, in any order.
4. The session ends in one of these ways. Either way no further turn
   starts, lines still queued are dropped (never run: a queued `!command`
   does not run on the way out), and approvals still open are refused,
   never approved. The session is saved, and `session.ended` is the last
   event.
   - `session.shutdown`: answered at once, then the session ends the way
     ctrl+d at the prompt ends a terminal session. A running turn is told
     to stop first.
   - The client's end of the pipe closing, or SIGTERM: the way a closed
     terminal ends a session -- a hangup, exit code 129 (143 for SIGTERM).
     It does not wait for a running tool call or `!command`. Background
     shells and MCP servers end too, and the process is gone within the
     REPL's leave deadline (10 s).
   - `/exit` sent as a prompt.

   `session.ended.reason` says which: `shutdown`, `eof`, `exit`, `signal`
   or `closed`.

## Requests

| Method | Params | Result | Meaning |
|---|---|---|---|
| `initialize` | `protocolVersion`, `client?` | `protocolVersion`, `pid`, `session`, `state`, `pending[]` | Handshake. Repeatable: a window attaching later gets the state and the open approvals. |
| `session.prompt` | `text`, `how?` | `accepted`, `queued?`, `steered?` | A line exactly as the user typed it. See "Prompts". |
| `session.cancel` | | `stopping`, `refused` | Esc. See "Cancelling". |
| `session.setMode` | `mode` | `approval` | Approval posture: `plan`, `default` or `acceptEdits`. Never `bypassPermissions`; that needs a new session started with `--unattended`. |
| `session.status` | | engine counters, `state`, `approval`, `queued` | The engine's token and cost counters, plus the session's state. |
| `session.replay` | `turns?` | `events[]` | The durable events from the start of the `turns`-th last turn (0 means all). See "Replay". |
| `session.shutdown` | `ifIdle?`, `dryRun?` | `closing`, `busy?` | End the session as described in Lifecycle. `ifIdle`: only when nothing runs (a turn, a question, a background shell, a subagent); otherwise `closing: false` and `busy` says what runs. `dryRun`: only say so, never close (an update deciding whether a restart would cut anything). |
| `permission.reply` | `requestId`, `choice?`, `answers?`, `reason?`, `persist?` | `applied`, `message` | Answer an open approval. See "Approvals". |
| `session.commands` | | `commands[]`: `name`, `summary`, `category`, `takesArgs`, `source` | Every `/command` a line may start with. See "Completion and pickers". |
| `session.completePath` | `fragment` | `paths[]` | What `@fragment` completes to, under the workspace only. |
| `session.choices` | `what`: `model`, `effort` or `mode` | `current`, `choices[]`: `value`, `detail?` | What a picker offers. |
| `session.toggleTasks` | | `showTasks` | ctrl+t: the open task list at the end of each turn, on or off, as in the REPL (a notice says which). |
| `session.setTitle` | `title` | `title` | What the session is called: what it is saved under and every window shows (a `session.updated` follows); one plain line of at most 80 characters (escapes and control characters become spaces). Its address -- the `name` it started with, or its id, by which mail, approvals and its presence record reach it -- does not change. Answered for a sleeping session too, without waking it. |

### Completion and pickers

The window completes `/commands` and `@paths` while the user types, and
opens pickers for the model, the effort and the agent mode. The lists are
the session's:

- `session.commands` lists the commands in the order a typed line tries
  them: the REPL's builtins, the user's and the workspace's own commands
  (`.delfin/commands/*.md`), skills, then the subagent shortcuts
  (`/explore`, `/review`, `/plan`, `/delegate`). A name that an earlier
  kind already takes is left out, because typing it runs the earlier one.
  `source` says which kind: `builtin`, `user`, `skill` or `subagent`.
- `session.completePath` follows the REPL completer's rules: paths under
  the workspace that start with the fragment (a leading `@` is ignored),
  directories with a trailing `/`, dot files only when the fragment asks
  for them, at most 50. Nothing outside the workspace is offered. A
  fragment ending in `/` lists that folder, so completion can go into
  one; what was typed before the last `/` is kept as typed.
- `session.choices` lists the models the session's provider serves (asked
  from its `/v1/models`, and reused for ten minutes; `detail` says where a
  model runs and, when it has a price, what it costs), the effort levels,
  or the agent modes. It is answered on a thread of its own, so a slow
  model list never holds up the next request. When a KIT session's list
  cannot be read, it offers the curated roster
  (`model_capabilities._KIT_PREFERENCE`) with a detail saying so; the
  current model is always among the choices.

A pick is applied by sending what the user would type: `/model <name>`,
`/effort <level>`, `/mode <name>` as a `session.prompt`, so it is queued
behind a running turn like any typed command. `/model <name>` is the
REPL's own command (`repl_commands._model`, then the client's
`switch_model`): on an OpenAI-compatible endpoint a name its
`/v1/models` does not list as a chat model is refused, with the nearest
names, and the model stays; when the list cannot be read the switch is
made and its notice says the name was not checked. A switch changes the
model only -- a `session.updated` with `model` follows; the approval
posture, the permission rules and the managed policy stay as they are. The approval posture is
`session.setMode`.

### Prompts

A prompt is handled as the REPL handles a line: a message, a `/command`,
`!shell` or a `# note`. A message is one message whatever it contains; the
REPL's ways of typing several lines into one prompt (a `"""` block, a
trailing backslash) do not apply. Where it goes depends on what the session
is doing when it arrives:

- **At the prompt:** it is handled next (`{accepted}`); several sent at
  once are handled in the order sent.
- **During a turn:** it is queued and goes out when the turn ends
  (`{accepted, queued}`), or, with `how: "steer"`, it goes into the running
  turn (`{accepted, steered}`). A steer is taken only while the model's turn
  is live; after that, and on a backend that takes nothing mid-turn, it is
  queued instead. A steer the turn never read (it ended first) goes out as
  the next message, with a notice saying so.
- **Between a turn and an automatic follow-up** (a continuation after the
  token limit, waiting out an endpoint outage, announced work): it is
  queued and goes before the follow-up, which is dropped. A message
  therefore also ends the wait during an outage at once.
- **While the session closes:** refused with -32003.

The text is stripped first, as a typed line is (`"  !ls"` is a shell
escape). Every message is shown once, when it arrives: `user.message` with
`source: "user"` at the prompt, `"queued"` behind a turn, `"steer"` into
one. It is not shown again when it runs.

### Cancelling

`session.cancel` ends the turn that is running, or keeps the turn of the
line in progress (one just sent and not started yet) from starting. A tool
call already running finishes first. It also stops the automatic follow-ups
of that turn, including the wait during an outage. Every open approval is
refused, and so is every one the ending turn still asks.

A cancel names what is in progress when it arrives, nothing later: a
message sent after it runs normally, and what it put in place (the refusal
of questions, the engine's stop request) ends when the session is back at
the prompt. A running `!command` or `/command` cannot be cancelled, as in
the REPL. `stopping` is false when nothing was running.

### Replay

A turn starts at its `turn.started`, together with the messages sent right
before it, which are what started it. A message queued or steered while an
earlier turn ran stays inside that turn, where it happened, and does not
count as a turn of its own. The first event of a replay is a `session.info`
as the session is now (current approval and model), whatever the log holds;
it is made for the replay and carries `seq: 0`.

## Events

Durable events are written to the session's event log and come back in a
replay; only they carry `seq`. Deltas, progress lines, state changes and
usage samples are live only; the matching `*.ended` and `turn.ended` events
carry their final values.

| Event | Fields | Durable |
|---|---|---|
| `worker.ready` | `protocolVersion`, `pid` | |
| `session.info` | `sessionId`, `workspace`, `provider`, `model`, `mode`, `approval`, `title?`, `banner[]`, `notices[]`, `tasks[]?` | ✓ |
| `session.updated` | the fields that changed: `approval?`, `model?`, `effort?`, `title?`, `tasks[]?` | ✓ |
| `session.state` | `state`: `idle`, `running`, `needs_input` or `ended` | |
| `session.ended` | `reason`, `exitCode` | ✓ |
| `user.message` | `text`, `source`: `user`, `queued`, `steer` or `harness` | ✓ |
| `turn.started` | `turnId`, `source`: `user` or `harness` | ✓ |
| `turn.ended` | `turnId`, `durationMs`, `inputTokens`, `outputTokens`, `costUsd`, `toolCalls`, `error?`, `stopped?` (ended by a cancel) | ✓ |
| `text.started` / `text.delta` / `text.ended` | `partId`; `delta`; full `text` | ended ✓ |
| `reasoning.started` / `.delta` / `.ended` | as for text | ended ✓ |
| `tool.started` | `turnId`, `callId`, `name`, `input` | ✓ |
| `tool.ended` | `turnId`, `callId`, `name`, `output`, `meta` | ✓ |
| `tool.denied` | `turnId`, `callId?`, `name` | ✓ |
| `notice` | `text` | ✓ |
| `context` | `title`, `text`: what the session was started to look into (`session.create`'s brief), shown as one row | ✓ |
| `wait` | `text`: a progress line the next one replaces | |
| `error` | `text` | ✓ |
| `usage` | `inputTokens`, `outputTokens`, `costUsd`: session totals, sampled during a turn | |
| `permission.asked` | `requestId`, `kind`, `tool`, `preview`, `options[]`, `flags`, `suggestion?`, `readDir?`; for a question also `question`, `header?`, `choices[]`, `choiceDetails[]`, `multiSelect` | ✓ |
| `permission.resolved` | `requestId`, `decision`, `by`: `client`, `outside`, `expired` or `abort` | ✓ |

**Field details:**

- `session.state` is derived from what the session is doing: `needs_input`
  while an approval is open; `running` whenever the session is busy on its
  own (a turn, a command, a follow-up, waiting out an outage); `idle` at the
  prompt.
- `session.updated` follows every change of posture, model or effort,
  however it happened: a plan approval, "stop asking about writes",
  `/permissions` or `/model` sent as a prompt, `session.setMode`.
- `tasks` is the session's task list as it is now, whole: the rows the
  dashboard and the REPL draw (`task_ticker.rows`), in their order (in
  progress, blocked, pending, completed), each `num` (the number the user
  reads), `status` (`pending`, `in_progress`, `blocked`, `completed`) and
  `label` (the active form while in progress, the reason appended while
  blocked). A label is the model's text: it arrives as one plain line,
  escapes and control characters made spaces, at most 160 characters. The
  worker looks at the task store on every tick and sends a
  `session.updated` with `tasks` when the rows changed; a replay's header
  carries the current list.
- `harness` marks something nobody typed: a finished background job, a
  message from another session, a continuation.
- `turnId` is `t1`, `t2`, … per session; a resumed session numbers on.
- `partId` and `callId` are unique within the session.
- `tool.ended` carries the head of the output that the engine forwards.
  `meta.chars` is the true size, `meta.truncated` says whether the output
  was cut, and `meta.ok`/`meta.error` say whether the call succeeded.
- `notice` is harness speech: command output and reports. The output of one
  command (`/help`, `!ls`) is one notice whose lines keep their structure,
  blank lines included.

**Grouping rules:**

- A run of answer deltas is one text part. It is closed by the next tool
  call, notice, approval, or the end of the turn, so text and tool rows
  interleave in the order they happened.
- The engine reports tool results by tool name, so results are matched to
  calls first-in-first-out per name.

**Credentials.** Every event except the live stream of the model's own words
(`text.delta`, `reasoning.delta`) and plain numbers goes through the secret
redactor before it is sent or kept: the exact values of the credentials
DELFIN holds, and anything shaped like a credential. So does every result
and error message (an open question's preview in `initialize`, a persisted
pattern, an error naming the command). A token in a command or a tool's
output therefore never reaches the window or the event log. The
identifiers (`session`, `callId`, `requestId`, …) are left alone.

## Approvals

The worker asks through the broker the terminal REPL uses. So a question
offers exactly the options the REPL offers for it (`terminal_confirm.options_for`),
and each key means what it means there.

| kind | answer |
|---|---|
| `confirm` | `choice`: `y` yes, `n` no, `A` always this exact command, `k` always commands like this (see `suggestion`), `R` always read this directory (see `readDir`), `e` stop asking about writes this session, `a` abort (refuse and end the turn), `?` explain |
| `plan` | `choice`: `d` approve, then default; `e` approve, then acceptEdits; `n` reject |
| `ask` | `answers`: the chosen labels out of `choices`; none dismisses the question |

**Rules for answers:**

- A key the request did not offer is refused with code -32602, and nothing
  is applied. So is a question answered with a `choice`, an answer that is
  not one of its labels or is given twice, and more than one answer to a
  question that is not `multiSelect`.
- An unknown or already-resolved `requestId` is refused with code -32002.
- `A` and `k` save a rule that is global and permanent. The REPL asks a
  second time before saving it; here the second yes is `persist: true`.
  Without it the command is approved just this once, and `message` says so.
- `reason` goes with a refusal to the model, the same as in the REPL. With a
  rejected plan, where the engine has no room for a reason, it goes into the
  running turn as a steered message ("Plan rejected: …").
- `applied` is false when the request ended another way in the meantime
  (answered from outside, aborted); nothing of this answer was applied then.

Nothing is answered on the user's behalf. A request with no window attached
waits, and the session's state is `needs_input` until the request is
answered.

## Errors

| Code | Meaning |
|---|---|
| -32700 | parse error (a malformed line) |
| -32601 | unknown method |
| -32602 | invalid params, or an answer the request did not offer |
| -32603 | internal error |
| -32001 | protocol version mismatch |
| -32002 | unknown or already-resolved approval |
| -32003 | not allowed: before `initialize`, or while the session closes |

## The service

`delfin-agent tui` connects to the agent service on this host, starting it
when none runs (`delfin.agent.service.launcher`). The service
(`delfin.agent.service.supervisor`) is one process per user and host,
detached from the terminal, so sessions keep working after their window
closes; the next window attaches to them again.

- **Where:** the socket `/tmp/delfin-<uid>/agent.sock` (mode 0600, in a
  0700 directory that must be the user's own; `DELFIN_RUNTIME_DIR` moves
  it), one service per host by `flock`. On a cluster, `$HOME` is shared but
  a socket is not: `~/.delfin/agent_service/run/<host>.json` says where a
  service runs, so a window on another login node can say so.
- **Who:** only the user's own processes (the peer's uid is checked).
  Through this socket a window answers approvals and sends prompts, so an
  agent's own commands must not reach it. Three layers keep them away:
  1. **The cage.** Every command an agent runs is in the process cage
     (bubblewrap; the Landlock socket guard where bubblewrap cannot run;
     Seatbelt on macOS). The cage masks the socket's directory as a
     "session door", like the tmux and screen sockets, so the command does
     not see the socket at all (`api_client._session_doors`).
  2. **Watch-only.** A connection from a process that looks started by an
     agent session may watch but not act. That is a descendant of a worker
     (on Linux each worker is a child subreaper, so a process a command
     detaches stays its descendant while the worker runs), one that carries
     `DELFIN_AGENT_SESSION`, or, on Linux, one whose environment cannot be
     read. `hello` says `mayAct: false`, and the
     connection cannot create, prompt, cancel, change modes, answer
     approvals, send messages or stop the service (-32003).
  3. **The gate.** The socket and run files are answer locations for the
     permission gate (`approval_guard`): a command that names them is
     refused, like one that names the approvals directory.

  The limit, as for upstream's approval files: on a host where no cage can
  run, a command that deliberately hides where it came from (detaching and
  scrubbing its environment) can still connect. No secret would help: there
  is none the user's window can hold that the user's other processes cannot
  read.
- **When it ends:** `service.shutdown` (`delfin-agent service stop`),
  SIGTERM, an emergency stop (`delfin-agent stop-all`), or `--idle-exit`
  seconds (default 900) with no session and no window. Every session gets a
  `session.shutdown` first and saves; what does not end in time gets
  SIGTERM, then SIGKILL. The service lives under no dashboard lifeline: a
  TUI started from a dashboard terminal keeps its sessions when the
  dashboard ends.

| Method | Params | Result |
|---|---|---|
| `hello` | `protocolVersion`, `client?` | `protocolVersion`, `serviceVersion`, `release`, `pid`, `host`, `mayAct`, `sessions[]`. Must come first. `serviceVersion` says what the service does beyond the protocol (2: a dry-run `session.shutdown`; a service before it ignores `dryRun` and would end the session, so a launcher never asks it. 3: the `attention.*` requests and event, and `session.create`'s `brief`). |
| `sessions.list` | | `sessions[]` (summaries), `elsewhere[]`: other open DELFIN sessions (dashboards, terminals) with `key`, `sessionId?`, `title?`, `workspace`, `host` |
| `session.create` | `cwd`, `name?`, `title?`, `model?`, `provider?`, `effort?`, `mode?`, `permissionMode?`, `resume?`, `continue?`, `fork?`, `addDirs?`, `readDirs?`, `isolate?`, `noIsolate?`, `fake?`, `fakeDelay?`, `fakeSave?`, `attach?`, `brief?`, `briefTitle?` | `sessionId`, `session`, `state`, `pending[]`. Starts a worker with chat's own flags; the creating window is attached unless `attach: false`. A resumed or copied session whose conversation this service's event log has nothing of (it ran in the dashboard or the REPL; a copy has a new id) shows it as it opens: a notice, then its last 20 turns as `user.message`, `tool.started`/`tool.ended` and `text.ended` events (durable: a later attach replays them), without the rows the harness wrote. `fake` sessions (tests) are scripted and not saved, unless `fakeSave`; resumed, they restore what the store kept, and `fork` copies them under a new id. `name` is the session's address (chat's `-n`), `title` what it is called (the name when not given; a session resumed without either comes back under the title it was saved with, never under the store's stand-in, its first message). A session that cannot start is refused with its reason (-32005); one another process writes is refused with a hint to open a copy (`resume` with `fork`). `brief`: what the session is started to look into (an item of the agent view's list): the worker gets it in `initialize`, not on its command line (every user of a login node can read that), shows it as a durable `context` event (`title` is `briefTitle`) and puts it in front of the first ordinary prompt (not a `/command` or `!` line); the `user.message` event keeps what was typed. |
| `session.attach` | `sessionId`, `turns?` | `session`, `state`, `pending[]`, `summary`, `events[]`: the replay of the last `turns` turns (all for 0), then live events. |
| `session.detach` | `sessionId` | `detached` |
| `session.message` | `key`, `text` | `sent`: the text goes into that session's inbox, as a message from another session (`session_messages`) |
| `sessions.past` | `workspace?`, `query?`, `limit?` | `sessions[]`: saved sessions, newest first, with `sessionId`, `title`, `workspace`, `model`, `provider`, `updatedAt`, `messages`, `costUsd` and `running` (whether it runs in this service). `workspace` scopes the list, `query` searches titles and first messages (the session store's own listing and search). |
| `service.status` | | `pid`, `host`, `startedAt`, `sessions[]`, `windows`, `socket` |
| `service.shutdown` | | `stopping`, `sessions` |
| `attention.list` | | `items[]`: what needs the user, from every DELFIN process (below); `version`: grows with every change, from the service's start on (seconds × 10⁶): a window keeps the newest list or event |
| `attention.answer` | `id`, `text` | `ok`, `error?`: a parked question's answer, written into the attention inbox; the next session that starts a turn in its workspace gets it, once (`attention.drain_resolved`, scoped by the engine to its own workspace; an item with no workspace goes to the next session anywhere). Not for a watch-only connection |
| `attention.dismiss` | `id` | `ok`, `error?`: cleared without an answer, never handed to a model. Not for a watch-only connection |

**What needs the user** (`attention_feed.py`). The service reads the
attention inbox (`~/.delfin/attention_inbox.jsonl`: questions a session
parked, confirmations, failed and finished runs, budget and security
notices, from dashboards, the REPL, the job monitor and the scheduler) and
the job monitor's findings (`~/.delfin/monitor_findings.jsonl`), through the
modules' own readers. It looks at their size and age every 2 s while a
window is connected (and reads them again at least every minute), and sends
`attention.updated` to every window when the items differ. An item:
`id` (the inbox's), `kind` (question, confirm, plan, failed, alert,
finished, budget, notice), `group` (needs: what blocks someone or failed;
notices), `blocking`, `name` (what it is about), `text` (one line), `detail`,
`workspace`, `host` (which node wrote it), `sessionId`, `sessionTitle?`,
`sessionOpen?` (another process has the session open: a window opens a
copy), `sessionSaved?` (the store has it: a window resumes it),
`createdAt`, `options[]`, and for a job `job`: `id`, `state`, `folder`,
`signatures`, `diagnosisSession`, `summary`, `fix` (the job monitor's
assessment; the prepared fix, if any, under its own `fix`). The inbox
decides what is open: a finding shows only while its job's `run_failed`
entry is pending. An entry of a session this service runs right now is
left out (the session says it itself; looked at again when one starts or
ends), as are expired notices. Every text is masked as a session's events
are (`redact.Redactor`): the records hold what other processes wrote, and
every window and watching agent process receives them.

Every session request of the worker protocol (`session.prompt`,
`session.cancel`, `session.setMode`, `session.status`, `session.replay`,
`session.shutdown`, `permission.reply`, `session.commands`,
`session.completePath`, `session.choices`, `session.toggleTasks`,
`session.setTitle`) is sent with `sessionId` and
forwarded to that session's worker; its answer comes back unchanged.

Events from the service itself, to every window:

| Event | Fields |
|---|---|
| `session.summary` | `summary`: `sessionId`, `title?`, `workspace`, `provider`, `model`, `approval`, `state`, `pending` (open questions), `costUsd`, `createdAt`, `lastActivity`, `pid`, and what a row says about it, taken from its events (no model writes them): `lastPrompt?` (the last message sent to it), `activity?` (what it does now: "bash  pytest -q", "writing…"), `waiting?` (the question it waits on), `answer?` (the first line of its last answer). Sent whenever one of them changes: what the agent view draws. |
| `session.removed` | `reason`, `exitCode`: the session's worker ended (shut down, exited, or failed). |
| `attention.updated` | `items[]`: all of them, as `attention.list` has them; `version`, as there. |

**Attaching while a session runs.** The replay is read from the session's
event log in the same step of the service's loop in which the window is
subscribed. An event the worker sends after that is either in the replay --
then the window does not get it again (its `seq` is not above the
replay's last) -- or arrives live. Live-only events (deltas) may continue a
part the replay already ended; a window ignores deltas for a part it has
seen end.

**Answers.** An approval may be answered from any window attached to the
session. The worker takes the first answer; a later one is refused with
-32002, and every window gets `permission.resolved`.

| Code | Meaning |
|---|---|
| -32004 | no such session in this service |
| -32005 | the session's worker has ended, or did not start (the message says why) |

## Debugging

- `DELFIN_TUI_TRACE=1` writes every request the worker receives to its
  stderr (the worker's log).
- `DELFIN_TUI_STACKS=1` makes `kill -USR1 <worker pid>` write every
  thread's stack to the same place, for a session that hangs.

## Versioning

`PROTOCOL_VERSION` is an integer and must match exactly. Adding an optional
field or a new event type does not change it; a client ignores fields and
event types it does not know. Renaming or removing a field, or changing the
meaning of one, does change it.
