"""The wire protocol between the terminal UI and the agent service.

JSON-RPC 2.0, one JSON object per line, UTF-8 (newline-delimited JSON).
The client sends requests; the service answers each one and sends what
happens in a session as notifications whose method is ``event`` and whose
params carry a ``type``.

This module is the single definition of every message. The TypeScript side
does not keep a second copy: ``json_schema()`` is exported and the TUI's
types are generated from it, so the two cannot drift apart silently.

Field specs below are deliberately small: a type name, ``?`` for an
optional field, ``[type]`` for a list. ``object`` is any JSON object.
"""

from __future__ import annotations

import json
import threading
from typing import IO, Any, Callable, Iterator

PROTOCOL_VERSION = 1

#: A line longer than this is refused rather than read into memory. Tool
#: output travels inside events, so the limit is generous; a single event
#: that size is already a bug upstream of here.
MAX_LINE_BYTES = 32 * 1024 * 1024

# JSON-RPC error codes: the standard ones, then ours (-32000..-32099 is the
# range the specification leaves to the application).
PARSE_ERROR = -32700
INVALID_REQUEST = -32600
METHOD_NOT_FOUND = -32601
INVALID_PARAMS = -32602
INTERNAL_ERROR = -32603
VERSION_MISMATCH = -32001
UNKNOWN_REQUEST = -32002
NOT_ALLOWED = -32003

# ---------------------------------------------------------------------------
# Requests: client -> worker
# ---------------------------------------------------------------------------

REQUESTS: dict[str, dict[str, str]] = {
    # Must be the first request. The worker answers with its version and
    # the session it runs; a client that speaks another version is refused.
    # brief: what the session was started to look into (session.create's).
    "initialize": {"protocolVersion": "integer", "client": "string?",
                   "brief": "string?", "briefTitle": "string?"},
    # A line exactly as the user typed it: a message, a /command, !shell or
    # a # note. While a turn runs, how="queue" (the default) waits for the
    # turn to end and how="steer" sends it into the running turn.
    "session.prompt": {"text": "string", "how": "string?"},
    # Esc: end the running turn. A tool call already running finishes first.
    "session.cancel": {},
    "session.setMode": {"mode": "string", "confirmed": "boolean?"},
    "session.status": {},
    # The durable events of the last ``turns`` turns, for a window that
    # attaches to a session already in progress.
    "session.replay": {"turns": "integer?"},
    # ifIdle: only if nothing would be cut -- no turn, nothing queued, no
    # question open, no background shell or subagent running (the service
    # putting an unused session to sleep); else answered closing: false,
    # with why. dryRun: only say what it would cut (closing: false, busy).
    "session.shutdown": {"ifIdle": "boolean?", "dryRun": "boolean?"},
    # The answer to a permission.asked event: ``choice`` is one of the
    # option keys it offered (confirm, plan); ``answers`` the chosen labels
    # of a question (ask). ``text`` is a line typed for a question (ask
    # only): an option when it names one, else the answer as typed
    # (approval_answers.typed_answer); refused for a confirm or a plan.
    # ``persist`` is the second "yes" the keys that save a permanent rule
    # (A, k) need; without it they approve once.
    "permission.reply": {"requestId": "string", "choice": "string?",
                         "answers": "[string]?", "reason": "string?",
                         "persist": "boolean?", "text": "string?",
                         "command": "string?"},
    # What the window completes and picks from (choices.py): the session's
    # /commands, the paths under its workspace for an @fragment, and the
    # models, efforts or agent modes a picker offers (what: model | effort
    # | mode), or what may follow a /command (what: args, command: "/git").
    # A pick is applied by sending the command, as typed.
    "session.commands": {},
    "session.completePath": {"fragment": "string"},
    "session.choices": {"what": "string", "command": "string?"},
    # ctrl+t: show or hide the open task list at the end of each turn, as
    # in the REPL (its notice says which).
    "session.toggleTasks": {},
    # The session's name (the agent view's ctrl+r): what it is saved under.
    "session.setTitle": {"title": "string"},
}

# ---------------------------------------------------------------------------
# Results: what each request answers with
# ---------------------------------------------------------------------------

RESULTS: dict[str, dict[str, str]] = {
    "initialize": {"protocolVersion": "integer", "pid": "integer",
                   "session": "object", "state": "string",
                   "pending": "[object]"},
    # Exactly one of: started next (accepted only), queued behind the
    # running turn, or steered into it.
    "session.prompt": {"accepted": "boolean", "queued": "boolean?",
                       "steered": "boolean?"},
    "session.cancel": {"stopping": "boolean", "refused": "integer"},
    "session.setMode": {"approval": "string"},
    "session.status": {"state": "string", "approval": "string",
                       "queued": "integer"},
    "session.replay": {"events": "[object]"},
    "session.shutdown": {"closing": "boolean", "busy": "string?"},
    "permission.reply": {"applied": "boolean", "message": "string"},
    # name, summary, category, takesArgs, source: builtin | user | skill |
    # subagent -- in the order a typed line tries them.
    "session.commands": {"commands": "[object]"},
    "session.completePath": {"paths": "[string]"},
    # choices: value, detail? (where a model runs, what it costs); usage:
    # what a /command takes when its values cannot all be listed (args)
    "session.choices": {"current": "string", "choices": "[object]", "usage": "string?"},
    "session.toggleTasks": {"showTasks": "boolean"},
    "session.setTitle": {"title": "string"},
}

# ---------------------------------------------------------------------------
# Events: worker -> client (method "event"; every event also carries
# "session" and "ts", and every durable one a "seq")
# ---------------------------------------------------------------------------

EVENTS: dict[str, dict[str, str]] = {
    "worker.ready": {"protocolVersion": "integer", "pid": "integer"},
    # What a window shows at the top of a session.
    "session.info": {"sessionId": "string", "workspace": "string",
                     "provider": "string", "model": "string",
                     "mode": "string", "approval": "string",
                     "title": "string?", "banner": "[string]",
                     "notices": "[string]", "tasks": "[object]?",
                     "bypassBlocked": "string?", "bypassConfirmed": "boolean?",
                     "auto": "string?", "background": "[object]?"},
    # bypassBlocked: why shift+tab steps over bypassPermissions ("" when it
    # may stop there: the managed policy allows it and isolation holds
    # each command; posture_cycle.bypass_refusal); bypassConfirmed: its
    # first step was confirmed in this session. A window that is told
    # neither steps over it.
    # tasks: the session's task list (task_ticker.rows) as it is now, in
    # display order (in progress, blocked, pending, completed), each row
    # num, status, label -- sent whole whenever any of it changed.
    # auto: auto mode's one-line state (loop, iteration n/max, last
    # result, how many items need the user); "" while it is off.
    # background: what this session still has out -- shells, watched
    # jobs, background agents, wake-ups (background_view.snapshot) -- each
    # row group, id, kind, label, url, state, and since / at in epoch
    # seconds; no ages, so it is sent when the set changes, not each
    # second. The window draws the ages.
    "session.updated": {"approval": "string?", "model": "string?",
                        "effort": "string?", "title": "string?",
                        "tasks": "[object]?", "bypassBlocked": "string?",
                        "bypassConfirmed": "boolean?", "auto": "string?",
                        "background": "[object]?"},
    # idle | running | needs_input | ended. "running" covers everything the
    # session does on its own: a turn, a follow-up, waiting out an outage.
    "session.state": {"state": "string"},
    "session.ended": {"reason": "string", "exitCode": "integer"},
    # source: user | queued (sent while a turn ran; goes out after it) |
    # steer (sent into the running turn) | harness (nobody typed it: a
    # finished job, a message from another session, a continuation)
    "user.message": {"text": "string", "source": "string"},
    "turn.started": {"turnId": "string", "source": "string"},
    # stopped: ended by a cancel (Esc), not on its own.
    "turn.ended": {"turnId": "string", "durationMs": "integer",
                   "inputTokens": "integer", "outputTokens": "integer",
                   "costUsd": "number", "toolCalls": "integer",
                   "error": "string?", "stopped": "boolean?"},
    "text.started": {"turnId": "string", "partId": "string"},
    "text.delta": {"partId": "string", "delta": "string"},
    "text.ended": {"turnId": "string", "partId": "string", "text": "string"},
    "reasoning.started": {"turnId": "string", "partId": "string"},
    "reasoning.delta": {"partId": "string", "delta": "string"},
    "reasoning.ended": {"turnId": "string", "partId": "string",
                        "text": "string"},
    "tool.started": {"turnId": "string", "callId": "string",
                     "name": "string", "input": "object"},
    # output is the head the engine forwards; meta carries ok/error, the
    # true size (chars) and whether it was truncated.
    "tool.ended": {"turnId": "string", "callId": "string", "name": "string",
                   "output": "string", "meta": "object"},
    "tool.denied": {"turnId": "string", "callId": "string?",
                    "name": "string"},
    # Harness speech: command output, status reports, anything that is not
    # the model's answer. Lines keep their structure.
    "notice": {"text": "string"},
    # A progress line the next one replaces (waiting for the endpoint).
    "wait": {"text": "string"},
    "error": {"text": "string"},
    "usage": {"inputTokens": "integer", "outputTokens": "integer",
              "costUsd": "number"},
    # kind: confirm | plan | ask. options are the keys the answer may use,
    # exactly the ones the terminal REPL offers for the same request. A
    # question (ask) carries its choices as labels -- the answer names
    # them -- with a description each, and whether several may be picked.
    "permission.asked": {"requestId": "string", "kind": "string",
                         "tool": "string", "preview": "string",
                         "options": "[object]", "flags": "object",
                         "question": "string?", "header": "string?",
                         "choices": "[string]?",
                         "choiceDetails": "[string]?",
                         "multiSelect": "boolean?",
                         "suggestion": "string?", "readDir": "string?",
                         # The tool call's own one-line description (the
                         # bash tool requires one), masked; the model's
                         # statement, not a check of the command.
                         "purpose": "string?",
                         # A shell request's command, for a window to put
                         # in the prompt to change (key c); absent when it
                         # holds a secret the preview masks.
                         "command": "string?"},
    # by: client | outside | expired | abort
    "permission.resolved": {"requestId": "string", "decision": "string",
                            "by": "string"},
    # What a session was started with to look into (session.create's
    # brief: a failed job, a finished run, an alert from the agent view):
    # shown as one row, and put in front of the session's first prompt.
    "context": {"title": "string", "text": "string"},
}

# ---------------------------------------------------------------------------
# The service: window -> agent service (supervisor.py)
# ---------------------------------------------------------------------------

#: Requests a window sends to the service itself. Session requests (above)
#: are sent to the service too, with "sessionId", and forwarded.
SERVICE_REQUESTS: dict[str, dict[str, str]] = {
    "hello": {"protocolVersion": "integer", "client": "string?"},
    "sessions.list": {},
    # Start a session: a worker with these settings (chat's own flags).
    # name: the address the session is reached by (`-n`: mail, approvals,
    # presence) for its whole life; title: what it is called (`/rename`),
    # the name when not given.
    "session.create": {"cwd": "string", "name": "string?", "title": "string?",
                       "model": "string?",
                       "provider": "string?", "effort": "string?",
                       "mode": "string?", "permissionMode": "string?",
                       "resume": "string?", "continue": "boolean?",
                       "fork": "boolean?", "addDirs": "[string]?",
                       "readDirs": "[string]?", "isolate": "boolean?",
                       "noIsolate": "boolean?", "fake": "boolean?",
                       "fakeDelay": "number?", "fakeSave": "boolean?", "attach": "boolean?",
                       # What to look into (the context event); its title
                       # names the row.
                       "brief": "string?", "briefTitle": "string?"},
    "session.attach": {"sessionId": "string", "turns": "integer?"},
    "session.pending": {"sessionId": "string"},
    "session.detach": {"sessionId": "string"},
    # A message for a session elsewhere (dashboard, another terminal).
    "session.message": {"key": "string", "text": "string"},
    # Saved sessions to resume, newest first: the session store's, scoped
    # to a workspace or matched against a query, each with whether it runs
    # in this service now.
    "sessions.past": {"workspace": "string?", "query": "string?",
                      "limit": "integer?"},
    "service.status": {},
    # restart: a new service follows at once (an update); windows reconnect
    # to it instead of leaving.
    "service.shutdown": {"restart": "boolean?"},
    # What needs the user, from every DELFIN process (attention_feed.py):
    # the items of ATTENTION_ITEM. answer: a parked question's answer, for
    # the next session in its workspace; dismiss: cleared, never answered.
    "attention.list": {},
    "attention.answer": {"id": "string", "text": "string"},
    "attention.dismiss": {"id": "string"},
    # Swarm management (package SW2b; the logic is SW1's
    # delfin/agent/swarm/, the wrappers delfin/agent/service/swarm_rpc.py).
    # spec: the spec object, or specYaml: swarm.yaml text. swarm.start is
    # a dry run unless confirm=true -- the TUI sends confirm only after
    # the user approved the plan it showed.
    "swarm.status": {"spec": "object?", "specYaml": "string?"},
    "swarm.plan": {"spec": "object?", "specYaml": "string?"},
    "swarm.validate": {"spec": "object?", "specYaml": "string?"},
    "swarm.start": {"spec": "object?", "specYaml": "string?",
                    "dryRun": "boolean?", "confirm": "boolean?"},
    "swarm.pause": {},
    "swarm.resume": {},
    "swarm.triage.list": {"queue": "[object]"},
    "swarm.triage.answer": {"item": "object", "answer": "string"},
}

SERVICE_RESULTS: dict[str, dict[str, str]] = {
    # release: the DELFIN release the service runs (a window of a newer one
    # moves an idle service over to it).
    "hello": {"protocolVersion": "integer", "serviceVersion": "integer",
              "release": "string?", "pid": "integer", "host": "string",
              "mayAct": "boolean", "sessions": "[object]"},
    "sessions.list": {"sessions": "[object]", "elsewhere": "[object]"},
    # events: what the session said since its worker started, up to the
    # point the window is subscribed at (later events arrive live).
    "session.create": {"sessionId": "string", "session": "object",
                       "state": "string", "pending": "[object]",
                       "events": "[object]"},
    "session.attach": {"session": "object", "state": "string",
                       "pending": "[object]", "events": "[object]",
                       "summary": "object"},
    "session.pending": {"pending": "[object]"},
    "session.detach": {"detached": "boolean"},
    "session.message": {"sent": "boolean"},
    # sessionId, title, workspace, model, provider, updatedAt, running
    "sessions.past": {"sessions": "[object]"},
    "service.status": {"pid": "integer", "host": "string",
                       "startedAt": "number", "sessions": "[object]",
                       "dormant": "[object]?",
                       "windows": "integer", "socket": "string"},
    "service.shutdown": {"stopping": "boolean", "sessions": "integer"},
    "attention.list": {"items": "[object]", "version": "integer?"},
    "attention.answer": {"ok": "boolean", "error": "string?"},
    "attention.dismiss": {"ok": "boolean", "error": "string?"},
    # The status schema SW1 documents (status.py STATUS_SCHEMA); plan
    # carries the launch plan plus dry_run_text; start is a dry run
    # unless confirm=true. pause/resume answer available=false until SW1
    # ships the verbs. triage.answer: answered, item_id, logged.
    "swarm.status": {"schema_version": "integer", "goal": "string",
                     "packages": "[object]",
                     "waiting_for_user": "[object]"},
    "swarm.plan": {"goal": "string", "members": "[object]",
                   "start_spacing_s": "integer?", "memory_gate": "boolean?",
                   "dry_run_text": "string"},
    "swarm.validate": {"valid": "boolean", "errors": "[string]"},
    "swarm.start": {"dry_run": "boolean", "goal": "string",
                    "members": "[object]?", "dry_run_text": "string?",
                    "started": "[object]", "refused": "string?"},
    "swarm.pause": {"available": "boolean", "reason": "string"},
    "swarm.resume": {"available": "boolean", "reason": "string"},
    "swarm.triage.list": {"items": "[object]"},
    "swarm.triage.answer": {"answered": "boolean", "item_id": "string?",
                            "logged": "string?"},
}

#: Events from the service itself, to every window.
SERVICE_EVENTS: dict[str, dict[str, str]] = {
    # sessionId, title?, workspace, provider, model, approval, state,
    # pending (open questions), costUsd, createdAt, lastActivity, pid. A
    # summary's state may also be "sleeping": the session has no worker (it
    # slept, or a service before this one stopped or died) and wakes when
    # used -- attached, sent a prompt, resumed; a replay of it reads its log.
    "session.summary": {"summary": "object"},
    # The service stops on purpose -- `service stop`, an emergency stop, its
    # idle end, a signal -- and says so before it closes: a window then
    # leaves (its sessions are saved) instead of starting a new service, as
    # it does when a service vanishes without a word. restart: a new one
    # follows at once; reconnect to it.
    "service.stopping": {"reason": "string", "restart": "boolean"},
    "session.removed": {"reason": "string", "exitCode": "integer?"},
    # The attention items changed (a file the service watches did, or one
    # was answered or dismissed): all of them, as attention.list has them.
    "attention.updated": {"items": "[object]", "version": "integer?"},
}

#: An attention item (attention.list, attention.updated). kind: question |
#: confirm | plan | failed | alert | finished | budget | notice; group:
#: needs | notices; blocking: someone waits for an answer. name: what it is
#: about (a session's title, "job 4711 failed"); text: one line on what it
#: needs; detail: what a peek shows. sessionId: the session it belongs to,
#: if any -- sessionOpen: another process has it open (open a copy),
#: sessionSaved: the store has it (resume it). job: a job's id, state,
#: folder, error signatures, the diagnosis session, summary and fix.
ATTENTION_ITEM: dict[str, str] = {
    "id": "string", "kind": "string", "group": "string", "blocking": "boolean",
    "name": "string", "text": "string", "detail": "string", "workspace": "string",
    "host": "string", "sessionId": "string", "sessionTitle": "string?",
    "sessionOpen": "boolean?", "sessionSaved": "boolean?", "createdAt": "number",
    "options": "[string]", "job": "object?"}

#: The fields of a session summary (session.summary, sessions.list). The
#: last five say what the session is about in a few words, for the agent
#: view's rows (taken from its events; no model writes them): the first
#: message sent to it (an untitled session's name), the last one, what it
#: does now, the question it waits on, and the first line of its last
#: answer.
SESSION_SUMMARY: dict[str, str] = {
    "sessionId": "string", "title": "string?", "workspace": "string",
    "provider": "string", "model": "string", "approval": "string",
    "state": "string", "pending": "integer", "costUsd": "number",
    "createdAt": "number", "lastActivity": "number", "pid": "integer",
    "firstPrompt": "string?", "lastPrompt": "string?", "activity": "string?",
    "waiting": "string?", "answer": "string?"}

UNKNOWN_SESSION = -32004
WORKER_GONE = -32005

#: Events kept in the session's durable log (eventlog.py) and replayed to
#: a window that attaches later. Only these carry a sequence number. Deltas,
#: progress ticks, state changes and usage samples are live-only: the
#: ended/turn records carry their final state.
DURABLE_EVENTS = frozenset({
    "session.info", "session.updated", "session.ended",
    "user.message", "turn.started", "turn.ended",
    "text.ended", "reasoning.ended",
    "tool.started", "tool.ended", "tool.denied",
    "notice", "error", "permission.asked", "permission.resolved", "context",
})


# ---------------------------------------------------------------------------
# Building messages
# ---------------------------------------------------------------------------

#: The longest session title: what the session store keeps of one.
TITLE_MAX = 80

#: The longest task label a window is sent (a task row is one line).
TASK_LABEL_MAX = 160


def clean_title(text: Any) -> str:
    """A session title as it is kept: one line, control and invisible
    characters (escapes, zero-width marks) made spaces, single spaces, at
    most TITLE_MAX characters. "" when nothing printable is left."""
    return clean_line(text, TITLE_MAX)


def clean_line(text: Any, limit: int) -> str:
    """``text`` as one plain line of at most ``limit`` characters: whole
    escape sequences, then every control and invisible character, made
    spaces; runs of spaces made one."""
    import re
    import unicodedata

    # Whole escape sequences first (a colour, a title, a cursor move), then
    # any control character left.
    plain = re.sub(r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07\x1b]*(?:\x07|\x1b\\)?|.?)", " ", str(text or ""))
    chars = [" " if unicodedata.category(ch) in ("Cc", "Cf", "Zl", "Zp") else ch for ch in plain]
    return " ".join("".join(chars).split())[:limit].rstrip()


def request(msg_id: int | str, method: str, params: dict | None = None) -> dict:
    return {"jsonrpc": "2.0", "id": msg_id, "method": method,
            "params": dict(params or {})}


def result(msg_id: int | str | None, value: Any = None) -> dict:
    return {"jsonrpc": "2.0", "id": msg_id,
            "result": {} if value is None else value}


def error(msg_id: int | str | None, code: int, message: str,
          data: Any = None) -> dict:
    err: dict[str, Any] = {"code": int(code), "message": str(message)}
    if data is not None:
        err["data"] = data
    return {"jsonrpc": "2.0", "id": msg_id, "error": err}


def event(params: dict) -> dict:
    return {"jsonrpc": "2.0", "method": "event", "params": params}


def encode(message: dict) -> bytes:
    """One line: compact JSON, no raw newlines inside, newline at the end.

    ``ensure_ascii`` stays off so text arrives as text; JSON escapes every
    control character, so a newline in a value can never end the line.

    Never raises on content. A file name that is not UTF-8 reaches Python
    as a lone surrogate (``os.fsdecode``); strict encoding would raise on
    it and take the session down with one odd file name in a tool result.
    ``backslashreplace`` writes it as ``\\udcXX`` -- inside a JSON string
    that is the JSON escape for the same code unit, so the line stays valid
    JSON and decodes back to what was sent. A value JSON has no type for
    becomes its text.
    """
    return (json.dumps(message, ensure_ascii=False, separators=(",", ":"),
                       default=str)
            + "\n").encode("utf-8", errors="backslashreplace")


# ---------------------------------------------------------------------------
# Reading and writing lines
# ---------------------------------------------------------------------------

class LineWriter:
    """Thread-safe: one message is one write, so lines never interleave."""

    def __init__(self, stream: IO[bytes]):
        self._stream = stream
        self._lock = threading.Lock()
        self.closed = False

    def send(self, message: dict) -> bool:
        data = encode(message)
        with self._lock:
            if self.closed:
                return False
            try:
                self._write_all(data)
                self._stream.flush()
            except (BrokenPipeError, OSError, ValueError):
                # The reader went away. Later sends are dropped rather than
                # raised into whatever thread happened to be emitting.
                self.closed = True
                return False
        return True

    def _write_all(self, data: bytes) -> None:
        """All of ``data``, however many writes that takes.

        An unbuffered stream may take part of a line (a signal arriving
        mid-write) or, when the far end made the pipe non-blocking, none
        of it for now; both are waited out instead of cutting the line.
        """
        view = memoryview(data)
        while view:
            try:
                n = self._stream.write(view)
            except BlockingIOError as exc:
                n = int(getattr(exc, "characters_written", 0) or 0)
                _wait_writable(self._stream)
            if n is None:
                _wait_writable(self._stream)
                continue
            view = view[n:]


def _wait_writable(stream: IO[bytes], timeout: float = 1.0) -> None:
    import select
    import time

    try:
        select.select([], [stream.fileno()], [], timeout)
    except (AttributeError, OSError, ValueError):
        time.sleep(0.01)


def read_messages(stream: IO[bytes],
                  on_bad_line: Callable[[str], None] | None = None,
                  max_bytes: int = MAX_LINE_BYTES) -> Iterator[dict]:
    """Yield each JSON object read from ``stream`` until EOF.

    A line that is not a JSON object is reported through ``on_bad_line``
    and skipped: one malformed message must not end the session. A line
    over ``max_bytes`` is skipped the same way, without being held whole.
    """
    while True:
        line = stream.readline(max_bytes + 1)
        if not line:
            return
        if len(line) > max_bytes and not line.endswith(b"\n"):
            # Drain the rest of the oversized line.
            while True:
                more = stream.readline(max_bytes + 1)
                if not more or more.endswith(b"\n"):
                    break
            if on_bad_line:
                on_bad_line(f"line longer than {max_bytes} bytes skipped")
            continue
        text = line.decode("utf-8", errors="replace").strip()
        if not text:
            continue
        try:
            message = json.loads(text)
        except ValueError as exc:
            if on_bad_line:
                on_bad_line(f"not JSON: {exc}")
            continue
        if not isinstance(message, dict):
            if on_bad_line:
                on_bad_line("not a JSON object")
            continue
        yield message


# ---------------------------------------------------------------------------
# Validation and the schema the TypeScript side is generated from
# ---------------------------------------------------------------------------

_JSON_TYPES = {"string": str, "integer": int, "number": (int, float),
               "boolean": bool, "object": dict}


def _check_value(spec: str, value: Any) -> bool:
    if spec.startswith("[") and spec.endswith("]"):
        return isinstance(value, list) and all(
            _check_value(spec[1:-1], v) for v in value)
    expected = _JSON_TYPES[spec]
    if spec in ("integer", "number") and isinstance(value, bool):
        return False
    return isinstance(value, expected)


def invalid_params(method: str, params: Any) -> str:
    """Why ``params`` do not fit ``method``, or "" when they do."""
    spec = REQUESTS.get(method)
    if spec is None:
        return f"unknown method {method!r}"
    return params_problem(spec, params)


def params_problem(spec: dict[str, str], params: Any) -> str:
    """Why ``params`` do not fit ``spec`` (this module's field notation),
    or "" when they do."""
    if not isinstance(params, dict):
        return "params must be an object"
    for name, kind in spec.items():
        optional = kind.endswith("?")
        kind = kind.rstrip("?")
        if name not in params or params[name] is None:
            if not optional:
                return f"missing {name!r}"
            continue
        if not _check_value(kind, params[name]):
            return f"{name!r} must be {kind}"
    return ""


def _field_schema(kind: str) -> dict:
    if kind.startswith("["):
        return {"type": "array", "items": _field_schema(kind[1:-1])}
    return {"type": kind}


def _object_schema(spec: dict[str, str], extra: dict | None = None) -> dict:
    props: dict[str, Any] = {}
    required: list[str] = []
    for name, kind in spec.items():
        optional = kind.endswith("?")
        props[name] = _field_schema(kind.rstrip("?"))
        if not optional:
            required.append(name)
    for name, schema in (extra or {}).items():
        props[name] = schema
        required.append(name)
    return {"type": "object", "properties": props, "required": required,
            "additionalProperties": True}


def json_schema() -> dict:
    """Every request, result and event as JSON Schema; the TUI's types
    come from it."""
    common = {"session": {"type": "string"}, "ts": {"type": "number"}}

    def event_schema(type_: str, spec: dict[str, str]) -> dict:
        extra = {"type": {"const": type_}, **common}
        if type_ in DURABLE_EVENTS:
            extra["seq"] = {"type": "integer"}
        return _object_schema(spec, extra)

    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "delfin-agent TUI protocol",
        "protocolVersion": PROTOCOL_VERSION,
        "requests": {m: _object_schema(s) for m, s in REQUESTS.items()},
        "results": {m: _object_schema(s) for m, s in RESULTS.items()},
        "events": {t: event_schema(t, s) for t, s in EVENTS.items()},
        "durableEvents": sorted(DURABLE_EVENTS),
        "service": {
            "requests": {m: _object_schema(s) for m, s in SERVICE_REQUESTS.items()},
            "results": {m: _object_schema(s) for m, s in SERVICE_RESULTS.items()},
            "events": {t: event_schema(t, s) for t, s in SERVICE_EVENTS.items()},
            "sessionSummary": _object_schema(SESSION_SUMMARY),
            "attentionItem": _object_schema(ATTENTION_ITEM),
        },
        "errors": {name: value for name, value in globals().items()
                   if name.isupper() and isinstance(value, int)
                   and -32768 <= value <= -32000},
    }


if __name__ == "__main__":  # python -m delfin.agent.service.protocol
    print(json.dumps(json_schema(), indent=2))
