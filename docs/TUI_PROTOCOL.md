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
- **Phase 1 setup.** One worker process runs one session over its stdin and
  stdout. Phase 2 puts the supervisor between window and worker, and
  sessions are then reached over a Unix socket. The messages stay the same.

## Lifecycle

1. The worker starts and emits `worker.ready`, `session.info` and
   `session.state {idle}`.
2. The client sends `initialize {protocolVersion}`. Nothing else is accepted
   before it. A version mismatch is refused with code -32001.
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
     terminal ends a session. Background shells and MCP servers end too,
     and the process is gone within the REPL's leave deadline (10 s).
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
| `session.shutdown` | | `closing` | End the session as described in Lifecycle. |
| `permission.reply` | `requestId`, `choice?`, `answers?`, `reason?`, `persist?` | `applied`, `message` | Answer an open approval. See "Approvals". |

### Prompts

A prompt is handled as the REPL handles a line: a message, a `/command`,
`!shell` or a `# note`. A message is one message whatever it contains; the
REPL's ways of typing several lines into one prompt (a `"""` block, a
trailing backslash) do not apply. Where it goes depends on what the session
is doing when it arrives:

- **At the prompt:** it is handled next (`{accepted}`).
- **During a turn:** it is queued and goes out when the turn ends
  (`{accepted, queued}`), or, with `how: "steer"`, it goes into the running
  turn (`{accepted, steered}`). A backend that takes nothing mid-turn gets
  it queued instead.
- **Between a turn and an automatic follow-up** (a continuation after the
  token limit, waiting out an endpoint outage, announced work): it is
  queued and goes before the follow-up, which is dropped. A message
  therefore also ends the wait during an outage at once.
- **While the session closes:** refused with -32003.

A queued message is shown once, as `user.message {source: "queued"}` when
it is queued, and not again when it runs.

### Cancelling

`session.cancel` ends the turn that is running, or the one about to start
from a prompt just sent. A tool call already running finishes first. It
also stops the automatic follow-ups of that turn, including the wait during
an outage. Every open approval is refused, and so is every one the ending
turn still asks. `stopping` is false when nothing was running.

### Replay

A turn starts at its `turn.started`, together with the messages sent right
before it, which are what started it. A message queued or steered while an
earlier turn ran stays inside that turn, where it happened, and does not
count as a turn of its own. The first event of a replay is a `session.info`
as the session is now (current approval and model), whatever the log holds.

## Events

Durable events are written to the session's event log and come back in a
replay; only they carry `seq`. Deltas, progress lines, state changes and
usage samples are live only; the matching `*.ended` and `turn.ended` events
carry their final values.

| Event | Fields | Durable |
|---|---|---|
| `worker.ready` | `protocolVersion`, `pid` | |
| `session.info` | `sessionId`, `workspace`, `provider`, `model`, `mode`, `approval`, `title?`, `banner[]`, `notices[]` | ✓ |
| `session.updated` | the fields that changed: `approval?`, `model?`, `effort?`, `title?` | ✓ |
| `session.state` | `state`: `idle`, `running`, `needs_input` or `ended` | |
| `session.ended` | `reason`, `exitCode` | ✓ |
| `user.message` | `text`, `source`: `user`, `queued`, `steer` or `harness` | ✓ |
| `turn.started` | `turnId`, `source`: `user` or `harness` | ✓ |
| `turn.ended` | `turnId`, `durationMs`, `inputTokens`, `outputTokens`, `costUsd`, `toolCalls`, `error?` | ✓ |
| `text.started` / `text.delta` / `text.ended` | `partId`; `delta`; full `text` | ended ✓ |
| `reasoning.started` / `.delta` / `.ended` | as for text | ended ✓ |
| `tool.started` | `turnId`, `callId`, `name`, `input` | ✓ |
| `tool.ended` | `turnId`, `callId`, `name`, `output`, `meta` | ✓ |
| `tool.denied` | `turnId`, `callId?`, `name` | ✓ |
| `notice` | `text` | ✓ |
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
DELFIN holds, and anything shaped like a credential. A token in a command or
a tool's output therefore never reaches the window or the event log. The
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
