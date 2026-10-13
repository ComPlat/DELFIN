# TUI protocol

How the terminal UI (`delfin-agent tui`, TypeScript) talks to an agent
session (Python). The single definition of every message is
`delfin/agent/service/protocol.py`. `python -m delfin.agent.service.protocol`
prints it as JSON Schema, and the TUI's types are generated from that.
This page explains what the messages mean.

## Transport

- **Wire format.** JSON-RPC 2.0, one JSON object per line, UTF-8. JSON
  escapes every control character, so a newline never occurs inside a
  message. A line that is not a JSON object is reported with a parse-error
  response and skipped; it never ends the session.
- **Directions.** The client sends requests (`id`, `method`, `params`). Each
  one gets a `result` or an `error`.
- **Events.** The session reports what happens as notifications with
  `method: "event"`. Their `params` carry:
  - `type`, the event name;
  - `session`, the session id;
  - `seq`, increasing per session and never reused, also across a restart
    of the worker;
  - `ts`, unix time;
  - the event's own fields.
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
4. `session.shutdown`, or the client's end of the pipe closing, ends the
   session. Approvals still open are refused, never approved. The session is
   saved, and `session.ended` is the last event.

## Requests

| Method | Params | Meaning |
|---|---|---|
| `initialize` | `protocolVersion`, `client?` | Handshake. Repeatable: a window attaching later gets the state and the open approvals. |
| `session.prompt` | `text`, `how?` | A line exactly as the user typed it: a message, `/command`, `!shell` or `# note`, handled as the REPL handles it. When idle it starts a turn. During a turn it is queued (`how: "queue"`, the default; it runs when the turn ends) or sent into the running turn (`how: "steer"`). |
| `session.cancel` | | Esc: end the running turn. A tool call already running finishes first, and open approvals are refused. |
| `session.setMode` | `mode` | Approval posture: `plan`, `default` or `acceptEdits`. Never `bypassPermissions`; that needs a new session started with `--unattended`. |
| `session.status` | | The engine's token and cost counters, plus `state`, `approval` and `queued`. |
| `session.replay` | `turns?` | The durable events from the start of the `turns`-th last exchange (0 means all), with `session.info` first. |
| `session.shutdown` | | End the session as described in Lifecycle. |
| `permission.reply` | `requestId`, `choice?`, `answers?`, `reason?` | Answer an open approval. See below. |

## Events

The type marks a durable event: one that is written to the session's
event log and comes back in a replay. Deltas, progress ticks and usage
samples are live only; the matching `*.ended` and `turn.ended` events carry
their final values.

| Event | Fields | Durable |
|---|---|---|
| `worker.ready` | `protocolVersion`, `pid` | |
| `session.info` | `sessionId`, `workspace`, `provider`, `model`, `mode`, `approval`, `title?`, `banner[]`, `notices[]` | ✓ |
| `session.updated` | `approval?`, `model?`, `effort?`, `title?` | ✓ |
| `session.state` | `state`: `idle`, `running`, `needs_input` or `ended` | |
| `session.ended` | `reason`, `exitCode` | ✓ |
| `user.message` | `text`, `source`: `user`, `queued`, `steer` or `harness` | ✓ |
| `turn.started` | `turnId`, `source` | ✓ |
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
| `permission.asked` | `requestId`, `kind`, `tool`, `preview`, `options[]`, `flags`, `question?`, `choices[]?`, `suggestion?`, `readDir?` | ✓ |
| `permission.resolved` | `requestId`, `decision`, `by`: `client`, `outside`, `expired` or `abort` | ✓ |

**Field details:**

- `harness` marks a turn nobody typed: a finished background job, a message
  from another session, or the continuation of a turn cut off at the token
  limit.
- `turnId` is `t1`, `t2`, … per session.
- `partId` and `callId` are unique within the session.
- `tool.ended` carries the head of the output that the engine forwards.
  `meta.chars` is the true size, `meta.truncated` says whether the output
  was cut, and `meta.ok`/`meta.error` say whether the call succeeded.
- `notice` is harness speech: command output and reports. Lines keep their
  structure.

**Grouping rules:**

- A run of answer deltas is one text part. It is closed by the next tool
  call, notice, or the end of the turn, so text and tool rows interleave in
  the order they happened.
- The engine reports tool results by tool name, so results are matched to
  calls first-in-first-out per name.

## Approvals

The worker asks through the broker the terminal REPL uses. So a question
offers exactly the options the REPL offers for it (`terminal_confirm.options_for`),
and each key means what it means there.

| kind | keys |
|---|---|
| `confirm` | `y` yes, `n` no, `A` always this exact command, `k` always commands like this (see `suggestion`), `R` always read this directory (see `readDir`), `e` stop asking about writes this session, `a` abort (refuse and end the turn), `?` explain |
| `plan` | `d` approve, then default; `e` approve, then acceptEdits; `n` reject |
| `ask` | `answers`: the chosen strings out of `choices` |

**Rules for answers:**

- A key the request did not offer is refused with code -32602, and nothing
  is applied.
- An unknown or already-resolved `requestId` is refused with code -32002.
- `reason` goes with a refusal to the model, the same as in the REPL.

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
| -32003 | not allowed yet (before `initialize`) |

## Versioning

`PROTOCOL_VERSION` is an integer and must match exactly. Adding an optional
field or a new event type does not change it; a client ignores fields and
event types it does not know. Renaming or removing a field, or changing the
meaning of one, does change it.
