"""The agent service: the process the sessions live in, one per user and host.

    python -m delfin.agent.service.supervisor [--idle-exit SECONDS]

Started by ``delfin-agent tui`` when it is not running (``launcher.py``),
detached from the terminal, so a session keeps working after its window
closes and the next window attaches to it again -- the way tmux keeps a
shell, and the way agent CLIs keep their background agents.

    window ──socket── supervisor ──pipe── worker (one per session)
    window ──socket──/            \\──pipe── worker

* Every session is one worker process (``worker.py``) speaking the TUI
  protocol on its pipes. The supervisor starts it, hands it requests from
  any window, and fans its events out to every window attached to it.
* A window that attaches gets the session's last turns from its event log
  first (``eventlog.py``), then the live events; the sequence numbers drop
  what it would otherwise see twice.
* Every window hears a summary of every session whenever one changes
  (state, posture, cost, open questions) -- what the agent view draws.
* An approval is answered by whichever window answers first; the worker
  refuses a second answer, and every window hears ``permission.resolved``.
* A process an agent session started (its shell, a background job, a
  child of those) is answered watch-only: it cannot answer approvals,
  send prompts or end sessions through this socket (``_from_an_agent``).
  It is recognised by ancestry (each worker is a subreaper on Linux, so
  a detached process stays below its worker while the worker runs) or by
  the agent marker in its environment. A process that has left both --
  its worker ended, or it started outside the worker's tree and dropped
  the marker -- is not recognised; on hosts where the cage masks the
  socket (bubblewrap, the socket guard, Seatbelt) it cannot connect.

The service ends when ``delfin-agent service stop`` says so, on SIGTERM, on
an emergency stop (``stop_all``), and on its own after ``--idle-exit``
seconds with no session and no window. Sessions end with it: each worker
saves its session, so a later window can resume it.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import re
import signal
import socket
import sys
import time
from pathlib import Path
from typing import Any

from . import attention_feed, eventlog, protocol, rundir

#: What this service does beyond its protocol, for a launcher deciding
#: what it may ask. 2: a dry-run session.shutdown only says what it would
#: cut -- a service before it ignores dryRun and would end the session.
#: 3: attention.list/answer/dismiss and attention.updated (what needs the
#: user, from every DELFIN process), and session.create's brief.
SERVICE_VERSION = 3

#: Requests a window sends to the service itself (the one definition is
#: protocol.SERVICE_REQUESTS, which the TUI's types are generated from).
SERVICE_REQUESTS = protocol.SERVICE_REQUESTS

#: Session requests, forwarded to the session's worker (with "sessionId").
FORWARDED = frozenset({"session.prompt", "session.cancel", "session.setMode",
                       "session.status", "session.replay", "session.shutdown",
                       "permission.reply", "session.commands",
                       "session.completePath", "session.choices",
                       "session.toggleTasks", "session.setTitle"})

#: What a watch-only connection may not do.
ACTING = frozenset({"session.create", "session.prompt", "session.cancel",
                    "session.setMode", "session.shutdown", "permission.reply",
                    "session.message", "service.shutdown", "session.toggleTasks",
                    "session.setTitle", "attention.answer", "attention.dismiss",
                    # A swarm is created, started and answered by the
                    # user, never by an agent session: the same rule that
                    # keeps an agent from answering its own approvals
                    # (attention.answer) keeps it from launching a swarm
                    # or answering the swarm's queue for itself.
                    "swarm.start", "swarm.pause", "swarm.resume",
                    "swarm.triage.answer"})

UNKNOWN_SESSION = protocol.UNKNOWN_SESSION
WORKER_GONE = protocol.WORKER_GONE

#: A window that stays this far behind for CLIENT_STUCK_S is dropped (it
#: can attach again and replay): a stuck window must not make the service
#: hold its output. A burst -- two replays asked at once -- is not stuck:
#: a window that reads catches up within a second.
CLIENT_BACKLOG_BYTES = 32 * 1024 * 1024
CLIENT_STUCK_S = 15.0
#: Behind by this much, a window is dropped at once.
CLIENT_BACKLOG_HARD_BYTES = 256 * 1024 * 1024

#: Windows at once; more is a runaway process, not a person.
MAX_CLIENTS = 128

#: The longest line a window may send before its hello.
HELLO_LINE_BYTES = 64 * 1024

#: How long a session gets to start (MCP servers, a resumed session).
CREATE_TIMEOUT_S = 180.0

#: How long the sessions get to end and save when the service stops.
SHUTDOWN_GRACE_S = 20.0

#: How long a forwarded request may wait for its session's answer. Every
#: one is answered at once except a model list (session.choices: the
#: provider is asked) and a replay of a long log.
FORWARD_TIMEOUT_S = 120.0

#: The service keeps this much of each session in memory (the replay
#: budget), so an attach never reads the log in the event loop.
REPLAY_BUDGET_BYTES = eventlog.REPLAY_BUDGET_BYTES

#: A session no window shows, with nothing to do and nothing to ask, sleeps
#: after this long: its worker ends (the session is saved) and it wakes when
#: it is used again -- opened, answered, sent a message. 0: never.
#: DELFIN_AGENT_SLEEP_AFTER_S sets it.
SLEEP_AFTER_S = 3600.0

#: How long a dormant session (asleep, or left by a service that stopped or
#: died) is kept to be woken; a later service forgets older ones.
DORMANT_KEEP_S = 7 * 86400.0

#: The create settings a session is resumed with -- its brief too: a
#: session about an item waits for the user, and may sleep (or move to a
#: new service) before its first prompt; the worker's log says whether
#: the brief went out already.
RESUME_SETTINGS = ("cwd", "name", "title", "model", "provider", "effort", "mode", "permissionMode",
                   "addDirs", "readDirs", "isolate", "noIsolate", "fake", "fakeDelay", "fakeSave",
                   "brief", "briefTitle")

#: How long a wake waits for the worker of a service that died to finish
#: saving the session (it ends on its own when its pipe closes).
OLD_WORKER_WAIT_S = 30.0

#: How often the attention records are looked at (their size and age: a
#: change there means reading them); on a cluster they are on a network
#: file system, where a timed look is surer than a change notification.
ATTENTION_POLL_S = 2.0
#: Read again at least this often though nothing changed: a notice
#: expires, a session elsewhere closes.
ATTENTION_FULL_S = 60.0

#: Summary changes the roster may take with the next write (half a second
#: on): a turn's activity, many a second. Every other change -- state,
#: posture, name, a question -- is written at once.
ACTIVITY_EVENTS = frozenset({"turn.started", "tool.started", "text.started"})


LINUX = sys.platform.startswith("linux")


def _ppid(pid: int) -> int:
    if not LINUX:
        import subprocess
        try:
            out = subprocess.run(["ps", "-o", "ppid=", "-p", str(pid)],
                                 capture_output=True, text=True, timeout=2)
            return int(out.stdout.strip() or 0)
        except (OSError, ValueError, subprocess.SubprocessError):
            return 0
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return 0
    # The command name is in parentheses and may contain spaces.
    try:
        return int(stat.rsplit(")", 1)[1].split()[1])
    except (IndexError, ValueError):
        return 0


def _peer_pid(sock: socket.socket) -> int:
    try:
        if LINUX and hasattr(socket, "SO_PEERCRED"):
            import struct
            creds = sock.getsockopt(socket.SOL_SOCKET, socket.SO_PEERCRED,
                                    struct.calcsize("3i"))
            return struct.unpack("3i", creds)[0]
        if sys.platform == "darwin":
            # SOL_LOCAL (0), LOCAL_PEERPID (2): <sys/un.h>
            return int(sock.getsockopt(0, 2))
    except OSError:
        pass
    return 0


def _gist(text: str, limit: int = 120) -> str:
    """The first line that says something, cut to ``limit``."""
    for line in str(text or "").splitlines():
        line = " ".join(line.split()).lstrip("#*->` ").strip()
        if line:
            return line if len(line) <= limit else line[:limit - 1] + "…"
    return ""


def _answer_gist(text: str, limit: int = 120) -> str:
    """An answer in a line: its first sentence rather than its heading (an
    answer opens with "## Summary" as often as not), without the markdown
    that marks emphasis or code."""
    lines = [" ".join(line.split()) for line in str(text or "").splitlines()]
    lines = [line for line in lines if line.strip("#*->` ")]
    if not lines:
        return ""
    prose = [line for line in lines if not line.startswith("#")]
    line = re.sub(r"\*\*|__|`", "", (prose or lines)[0].lstrip("#*->` ")).strip()
    return line if len(line) <= limit else line[:limit - 1] + "…"


def _tool_gist(event: dict) -> str:
    """A tool call in a few words: its name and the argument that matters."""
    name = str(event.get("name") or "tool")
    if name.startswith("mcp__"):
        parts = name.split("__", 2)
        if len(parts) == 3:
            name = f"{parts[1]}:{parts[2]}"
    args = event.get("input") or {}
    if isinstance(args, dict):
        for key in ("command", "path", "file_path", "pattern", "query", "url", "description", "name"):
            value = args.get(key)
            if isinstance(value, str) and value.strip():
                return f"{name}  {_gist(value, 80)}"
    return name


def _asked_gist(asked: dict) -> str:
    """The question a session waits on, in a few words."""
    if asked.get("kind") == "ask":
        return _gist(str(asked.get("question") or "a question"))
    if asked.get("kind") == "plan":
        return "the plan"
    preview = str(asked.get("preview") or "")
    target = ""
    for line in preview.splitlines():
        if line.startswith("+++ "):
            target = line[4:].removeprefix("b/").strip()
            break
    return f"{asked.get('tool') or 'a tool'} {target or _gist(preview.lstrip('$ '), 70)}".strip()


def _environ_of(pid: int) -> bytes | None:
    """The environment of ``pid`` (one of the user's own processes), or
    None when it cannot be read. Linux reads /proc; macOS asks ps, which
    shows the environment of the user's own processes (-E)."""
    if LINUX:
        try:
            return Path(f"/proc/{pid}/environ").read_bytes()
        except OSError:
            return None
    import subprocess
    try:
        out = subprocess.run(["ps", "-wwE", "-o", "command=", "-p", str(pid)],
                             capture_output=True, timeout=2)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0 or not out.stdout.strip():
        return None
    return out.stdout


class WorkerFailed(Exception):
    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


class Session:
    """One worker process and what the service knows about it."""

    def __init__(self, supervisor: "Supervisor", proc: asyncio.subprocess.Process,
                 settings: dict):
        self.supervisor = supervisor
        self.proc = proc
        self.settings = settings
        self.id = ""
        self.info: dict = {}
        self.state = "starting"
        self.pending: dict[str, dict] = {}
        self.cost = 0.0
        self.created_at = time.time()
        self.last_activity = self.created_at
        self.ended_reason = ""
        #: Being put to sleep: its end makes it dormant, not gone.
        self.sleeping = False
        #: Not put to sleep before this (its worker said it was busy).
        self.sleep_not_before = 0.0
        # What the agent view's row says (SESSION_SUMMARY's last five).
        self.first_prompt = ""
        self.last_prompt = ""
        self.activity = ""
        self.answer = ""
        self.subscribers: dict["Client", int] = {}    # client -> last seq replayed
        self._next_id = 0
        self._waiting: dict[int, asyncio.Future] = {}
        self._early: list[dict] = []                  # events before it was registered
        #: The highest sequence number read from the worker. What the
        #: service knows (state, open questions) is exactly the events up
        #: to here: a replay is cut here and a window subscribed here, so
        #: what it is told and what it is shown always agree.
        self.seen_seq = 0
        #: The newest durable events, whole turns, within the replay
        #: budget: what an attach replays (see ``replay``).
        self.recent = eventlog.RecentEvents(REPLAY_BUDGET_BYTES)
        self.done = asyncio.get_running_loop().create_future()
        #: Called with every event, in the service's loop (a schedule's run
        #: follows its session this way: schedule_runner).
        self.listeners: list = []

    # -- talking to the worker ------------------------------------------------------
    async def request(self, method: str, params: dict | None = None,
                      timeout: float | None = None) -> Any:
        if self.proc.stdin is None or self.proc.returncode is not None:
            raise WorkerFailed(WORKER_GONE, "the session has ended")
        self._next_id += 1
        msg_id = self._next_id
        future = asyncio.get_running_loop().create_future()
        self._waiting[msg_id] = future
        loop = asyncio.get_running_loop()
        deadline = None if timeout is None else loop.time() + timeout
        try:
            self.proc.stdin.write(protocol.encode(
                protocol.request(msg_id, method, params or {})))
            # Under the same time limit as the answer: a worker that stops
            # reading its stdin (a full pipe) must not hold this request --
            # or the service's stop -- for ever.
            if deadline is None:
                await self.proc.stdin.drain()
            else:
                await asyncio.wait_for(self.proc.stdin.drain(), max(0.0, deadline - loop.time()))
        except (BrokenPipeError, ConnectionError, RuntimeError):
            self._waiting.pop(msg_id, None)
            raise WorkerFailed(WORKER_GONE, "the session has ended") from None
        except asyncio.TimeoutError:
            self._waiting.pop(msg_id, None)
            raise
        try:
            if deadline is None:
                return await future
            return await asyncio.wait_for(future, max(0.0, deadline - loop.time()))
        finally:
            self._waiting.pop(msg_id, None)

    async def read(self) -> None:
        """Every line the worker writes: answers and events, in order."""
        assert self.proc.stdout is not None
        try:
            while True:
                try:
                    line = await self.proc.stdout.readline()
                except ValueError:
                    # A line over the reader's limit: the worker bounds its
                    # own lines, so this is damage; skip to the next one.
                    continue
                if not line:
                    break
                try:
                    message = json.loads(line)
                except ValueError:
                    continue
                if not isinstance(message, dict):
                    continue
                if message.get("method") == "event":
                    self.on_event(message, message.get("params") or {})
                    continue
                future = self._waiting.get(message.get("id"))
                if future is None or future.done():
                    continue
                if "error" in message:
                    err = message.get("error") or {}
                    future.set_exception(WorkerFailed(
                        int(err.get("code", protocol.INTERNAL_ERROR)),
                        str(err.get("message", "error"))))
                else:
                    future.set_result(message.get("result"))
        finally:
            await self.proc.wait()
            reason = self.ended_reason or f"the worker exited ({self.proc.returncode})"
            for future in list(self._waiting.values()):
                if not future.done():
                    future.set_exception(WorkerFailed(WORKER_GONE, reason))
            self.state = "ended"
            if not self.done.done():
                self.done.set_result(self.proc.returncode)
            self.supervisor.session_ended(self)

    def apply(self, event: dict, *, count_cost: bool = True) -> bool:
        """What ``event`` changes in the service's view of the session;
        True when the summary changed."""
        kind = str(event.get("type") or "")
        if kind == "session.state":
            self.state = str(event.get("state") or self.state)
            return True
        if kind in ("session.info", "session.updated"):
            self.info.update({k: v for k, v in event.items()
                              if k not in ("type", "session", "seq", "ts")})
            return True
        if kind == "permission.asked":
            self.pending[str(event.get("requestId"))] = event
            return True
        if kind == "permission.resolved":
            return self.pending.pop(str(event.get("requestId")), None) is not None
        if kind == "user.message":
            self.last_prompt = _gist(str(event.get("text") or ""))
            self.first_prompt = self.first_prompt or self.last_prompt
            return True
        if kind == "turn.started":
            self.activity = "thinking…"
            return True
        if kind == "tool.started":
            self.activity = _tool_gist(event)
            return True
        if kind == "text.started":
            self.activity = "writing…"
            return True
        if kind == "text.ended":
            gist = _answer_gist(str(event.get("text") or ""))
            if gist:
                self.answer = gist
            return False
        if kind == "turn.ended":
            if count_cost:
                self.cost += float(event.get("costUsd") or 0.0)
            self.activity = ""
            return True
        if kind == "session.ended":
            self.ended_reason = str(event.get("reason") or "")
        return False

    def on_event(self, message: dict, event: dict) -> None:
        self.last_activity = time.time()
        changed = self.apply(event)
        for listener in list(self.listeners):
            try:
                listener(event)
            except Exception:
                pass
        seq = event.get("seq")
        if isinstance(seq, int) and seq > self.seen_seq:
            self.seen_seq = seq
            self.recent.add(event)
        if not self.id:
            # Not registered yet (still starting): the create answer
            # carries these (Supervisor._m_session_create).
            self._early.append(event)
            return
        data = protocol.encode(message)
        for client, replayed in list(self.subscribers.items()):
            if isinstance(seq, int) and seq <= replayed:
                continue                     # already in that window's replay
            client.send_bytes(data)
        if changed:
            self.supervisor.broadcast_summary(
                self, roster_now=str(event.get("type") or "") not in ACTIVITY_EVENTS)

    # -- what windows see ---------------------------------------------------------------
    def summary(self) -> dict:
        return {"sessionId": self.id,
                "title": (self.info.get("title") or self.settings.get("title")
                          or self.settings.get("name") or None),
                "workspace": self.info.get("workspace") or self.settings.get("cwd", ""),
                "provider": self.info.get("provider", ""),
                "model": self.info.get("model", ""),
                "approval": self.info.get("approval", ""),
                "state": self.state,
                "pending": len(self.pending),
                "costUsd": round(self.cost, 6),
                "createdAt": round(self.created_at, 3),
                "lastActivity": round(self.last_activity, 3),
                "pid": self.proc.pid,
                "firstPrompt": self.first_prompt or None,
                "lastPrompt": self.last_prompt or None,
                "activity": self.activity or None,
                "waiting": _asked_gist(next(iter(self.pending.values()))) if self.pending else None,
                "answer": self.answer or None}

    def resume_settings(self) -> dict:
        """What resumes this session as it is now: its create settings, with
        the model and the posture it has now -- a posture only if it is one
        shift+tab reaches (plan, default, acceptEdits): a wake must not
        bypass what was asked when it started."""
        out = {k: v for k, v in self.settings.items()
               if k in RESUME_SETTINGS and v not in (None, "", [])}
        out["cwd"] = str(self.settings.get("cwd") or self.info.get("workspace") or "")
        # Renamed, or switched to another effort or agent mode since it
        # started: it wakes as it is -- under the address it started with
        # (``name``), called what it is called now (``title``).
        for key, info_key in (("title", "title"), ("effort", "effort"), ("mode", "mode")):
            if self.info.get(info_key):
                out[key] = str(self.info[info_key])
        approval = str(self.info.get("approval") or "")
        if approval in ("plan", "default", "acceptEdits"):
            out["permissionMode"] = approval
        if not out.get("fake"):
            for key in ("model", "provider"):
                if self.info.get(key):
                    out[key] = str(self.info[key])
        return out

    def header(self) -> dict:
        event = {"type": "session.info", "session": self.id, "seq": 0,
                 "ts": round(time.time(), 3)}
        event.update({k: v for k, v in self.info.items() if v is not None})
        return event

    def replay(self, turns: int) -> list[dict]:
        """The last ``turns`` turns (0: all the service holds), up to the
        last event it has read -- the same point ``pending`` and ``state``
        stand at -- from memory, never from the log in the event loop."""
        events = eventlog.tail_turns(self.recent.events(), turns, header=self.header())
        if self.recent.dropped and (turns <= 0 or turns > self.recent.turns()):
            # Older turns did not fit the replay budget: say so, where they
            # would have been (right after the header).
            events.insert(1, {**eventlog.OLDER_LEFT_OUT, "session": self.id,
                              "ts": round(time.time(), 3)})
        return events

    async def load_history(self) -> None:
        """What the log already held when this service took the session
        on (a resumed session's earlier turns), merged under the events
        read live since. Read on a thread: a long log must not stall the
        service."""
        try:
            # Read from its end, cut to the budget and sized off the event
            # loop: the loop only splices what fits (a 200 MB log must not
            # stall every window, nor be held in memory whole).
            older, cut = await asyncio.to_thread(eventlog.read_newest,
                                                 eventlog.log_path(self.id),
                                                 REPLAY_BUDGET_BYTES)
        except Exception:
            return
        self.recent.merge_older_sized(older, cut)
        # A resumed session is named by what it was first asked (as far
        # back as the replay window reaches), not by the first message of
        # this run.
        first = next((e for e, _ in older if e.get("type") == "user.message"), None)
        if first is not None:
            self.first_prompt = _gist(str(first.get("text") or "")) or self.first_prompt


def _elsewhere(ours: set[str]) -> list[dict]:
    """The other open DELFIN sessions on this machine and the cluster's
    shared home -- dashboards and terminal sessions, to watch and message --
    leaving out ``ours`` (the sessions this service runs)."""
    out = []
    try:
        from .. import session_presence
        for record in session_presence.open_sessions():
            sid = str(record.get("session_id") or "")
            if sid and sid in ours:
                continue
            out.append({"key": record.get("key", ""),
                        "sessionId": sid or None,
                        "title": record.get("title") or None,
                        "workspace": record.get("workspace", ""),
                        "host": record.get("host", ""),
                        "lastActivity": record.get("updated_at")})
    except Exception:
        pass
    return out


class Client:
    """One window's connection."""

    def __init__(self, supervisor: "Supervisor", reader: asyncio.StreamReader,
                 writer: asyncio.StreamWriter, may_act: bool, peer_pid: int):
        self.supervisor = supervisor
        self.reader = reader
        self.writer = writer
        self.may_act = may_act
        self.peer_pid = peer_pid
        self.greeted = False
        self.closed = False
        self.name = ""                    # what it said it is, in its hello
        self._stuck_check: asyncio.TimerHandle | None = None

    def send(self, message: dict) -> None:
        self.send_bytes(protocol.encode(message))

    def send_bytes(self, data: bytes) -> None:
        if self.closed:
            return
        try:
            self.writer.write(data)
        except (ConnectionError, RuntimeError):
            self.close()
            return
        transport = self.writer.transport
        if transport is None:
            return
        behind = transport.get_write_buffer_size()
        if behind > CLIENT_BACKLOG_HARD_BYTES:
            self.close(abort=True)
        elif behind > CLIENT_BACKLOG_BYTES and self._stuck_check is None:
            self._stuck_check = asyncio.get_running_loop().call_later(
                CLIENT_STUCK_S, self._still_behind)

    def _still_behind(self) -> None:
        self._stuck_check = None
        transport = self.writer.transport
        if transport is not None and \
                transport.get_write_buffer_size() > CLIENT_BACKLOG_BYTES:
            # Dropped, backlog and all: a plain close would keep the whole
            # backlog in memory until the window read it.
            self.close(abort=True)

    def close(self, *, abort: bool = False) -> None:
        if self.closed:
            return
        self.closed = True
        if self._stuck_check is not None:
            self._stuck_check.cancel()
            self._stuck_check = None
        try:
            transport = self.writer.transport
            if abort and transport is not None:
                transport.abort()
            else:
                self.writer.close()
        except Exception:
            pass
        self.supervisor.client_gone(self)


class Supervisor:
    def __init__(self, *, idle_exit_s: float = 900.0,
                 worker_cmd: list[str] | None = None,
                 sleep_after_s: float = SLEEP_AFTER_S):
        self.sessions: dict[str, Session] = {}
        self.starting: set[Session] = set()
        self.clients: set[Client] = set()
        #: Sessions without a worker that wake when used: asleep, or left by
        #: a service that stopped or died. id -> {"summary", "settings",
        #: "since", "why", "oldPid"}; kept in the roster for a later service.
        self.dormant: dict[str, dict] = {}
        self._waking: dict[str, asyncio.Future] = {}
        self._roster_timer: asyncio.TimerHandle | None = None
        self.idle_exit_s = float(idle_exit_s)
        self.sleep_after_s = float(sleep_after_s)
        self.worker_cmd = worker_cmd or [sys.executable, "-m",
                                         "delfin.agent.service.worker"]
        self.started_at = time.time()
        self.quiet_since = time.monotonic()
        self._spawned = 0
        self.stopping: asyncio.Event | None = None
        #: Set once the service begins to stop, by any path, even before
        #: serve() runs (an emergency stop noticed while starting up): no
        #: session starts from here on and nothing acts.
        self.closing = False
        #: Why it stops (said to the windows), and whether a new one follows.
        self.stop_reason = "stop"
        self.restart = False
        self._creating = 0                # session.create calls in flight
        #: The tasks the service runs (requests, readers, connections):
        #: what has to end before it is gone.
        self._tasks: set[asyncio.Task] = set()
        self.server: asyncio.base_events.Server | None = None
        #: What needs the user (attention_feed), as last sent to the windows,
        #: and its version: a window keeps the newest of a list and an
        #: event that cross on the socket.
        self.attention_items: list[dict] = []
        #: From the service's start (seconds, times a million) on: a later
        #: service's lists are newer than any of an earlier one's.
        self.attention_version = int(time.time()) * 1_000_000
        self._attention_fp: tuple | None = None
        self._attention_at = 0.0
        self._attention_titles = attention_feed.Titles()
        self._attention_lock = asyncio.Lock()
        #: The user's schedules of this host, run as sessions here.
        from .schedule_runner import ScheduleRunner
        self.schedules = ScheduleRunner(self)

    def spawn(self, coro) -> asyncio.Task:
        task = asyncio.create_task(coro)
        self._tasks.add(task)
        task.add_done_callback(self._tasks.discard)
        return task

    # -- the socket -----------------------------------------------------------------------
    async def serve(self, path: Path) -> None:
        self.stopping = asyncio.Event()
        if self.closing:
            self.stopping.set()
        # Before the first window can ask: what a service before this one
        # left to be woken. The roster says from the start that this one
        # runs, and what it holds (a launcher of another release reads it).
        if self._load_dormant():
            self._write_roster()
        try:
            path.unlink()                       # stale: the lock is ours
        except FileNotFoundError:
            pass
        old = os.umask(0o177)
        try:
            self.server = await asyncio.start_unix_server(
                self._accept, path=str(path), limit=protocol.MAX_LINE_BYTES + 1)
        finally:
            os.umask(old)
        os.chmod(path, 0o600)
        idle = self.spawn(self._watch_idle())
        attention = self.spawn(self._watch_attention())
        timed = self.spawn(self.schedules.run())
        try:
            await self.stopping.wait()
        finally:
            self.closing = True
            idle.cancel()
            attention.cancel()
            timed.cancel()
            self.server.close()
            # Said before anything closes: a window leaves (or, for a
            # restart, reconnects) instead of starting a new service.
            data = protocol.encode(protocol.event({
                "type": "service.stopping", "session": "", "ts": round(time.time(), 3),
                "reason": self.stop_reason, "restart": self.restart}))
            for client in list(self.clients):
                if client.greeted:
                    client.send_bytes(data)
            # What ends with the service is saved; a later service offers
            # it again, to be woken (an update, a `service stop`).
            # Kept before the sessions end: a service killed while they
            # save loses none of them. A restart's windows are let go
            # first -- they reconnect to the next service, and the end of
            # a session that only sleeps across it is not theirs to show.
            for record in [self._dormant_record(s, "the agent service stopped")
                           for s in self.sessions.values() if s.id]:
                self.dormant[record["summary"]["sessionId"]] = record
            if self.restart:
                for session in self.sessions.values():
                    session.subscribers.clear()
            self._write_roster()
            await self.end_all_sessions()
            for client in list(self.clients):
                client.close()
            try:
                path.unlink()
            except FileNotFoundError:
                pass
            if self.dormant:
                self._write_roster(live=False)
            else:
                self._remove_roster()
            await self._settle_tasks()

    async def _settle_tasks(self) -> None:
        """Let what the service still runs end before the loop closes:
        requests whose window is gone, readers of ended workers."""
        others = [t for t in self._tasks if t is not asyncio.current_task()]
        if others:
            _done, pending = await asyncio.wait(others, timeout=2.0)
            for task in pending:
                task.cancel()
            if pending:
                await asyncio.wait(pending, timeout=1.0)

    async def _accept(self, reader: asyncio.StreamReader,
                      writer: asyncio.StreamWriter) -> None:
        sock = writer.get_extra_info("socket")
        uid = rundir.peer_uid(sock) if sock is not None else None
        if uid is not None and uid != os.getuid():
            writer.close()
            return
        if len(self.clients) >= MAX_CLIENTS:
            # A runaway process opening connections must not take the
            # service -- and the sessions in it -- down with it.
            writer.close()
            return
        pid = _peer_pid(sock) if sock is not None else 0
        client = Client(self, reader, writer, not self._from_an_agent(pid), pid)
        self.clients.add(client)
        self.quiet_since = time.monotonic()
        connection = asyncio.current_task()
        if connection is not None:
            self._tasks.add(connection)
            connection.add_done_callback(self._tasks.discard)
        try:
            while not client.closed:
                try:
                    line = await reader.readline()
                except (ValueError, ConnectionError):
                    break
                if not line:
                    break
                if not client.greeted and len(line) > HELLO_LINE_BYTES:
                    break                     # nobody says hello in megabytes
                try:
                    message = json.loads(line)
                except ValueError:
                    client.send(protocol.error(None, protocol.PARSE_ERROR, "not JSON"))
                    continue
                if not isinstance(message, dict):
                    client.send(protocol.error(None, protocol.PARSE_ERROR,
                                               "not a JSON object"))
                    continue
                # Each request runs on its own, so a slow session never holds
                # a window's other requests; requests to one session still
                # reach it in the order they came (request() writes first).
                self.spawn(self._dispatch(client, message))
        finally:
            client.close()

    def _from_an_agent(self, pid: int) -> bool:
        """Whether ``pid`` may have been started by an agent session: a
        descendant of one of the workers, a process carrying the marker
        every agent child process gets (approval_guard.ENV_MARKER), or one
        whose environment cannot be read (a process that made itself
        unreadable is not a window).

        Best effort by nature: a process of the same user can hide where it
        came from. Each worker is a child subreaper on Linux
        (``subreaper.keep_descendants``), so a process an agent command
        detaches (``setsid -f``, a double fork) is re-parented to the
        worker and the walk below still finds it while that worker runs;
        on other platforms, or once the worker has ended, such a process
        with the marker removed is answered as a window. What keeps an
        agent's commands away from this socket is the cage they run in:
        the bubblewrap process cage, the socket guard or Seatbelt mask the
        socket's directory as a "session door" (api_client._session_doors),
        like the tmux and screen sockets."""
        if pid <= 0:
            return True                       # who it is cannot be told: watch only
        workers = {s.proc.pid for s in [*self.sessions.values(), *self.starting]}
        p, hops = pid, 0
        while p > 1 and hops < 64:
            if p in workers:
                return True
            p, hops = _ppid(p), hops + 1
        try:
            from ..approval_guard import ENV_MARKER
        except ImportError:
            return False
        environ = _environ_of(pid)
        if environ is None:
            return True                       # unreadable: fail closed
        return f"{ENV_MARKER}=".encode() in environ

    def client_gone(self, client: Client) -> None:
        self.clients.discard(client)
        for session in self.sessions.values():
            session.subscribers.pop(client, None)
        if not self.clients:
            self.quiet_since = time.monotonic()

    # -- requests ------------------------------------------------------------------------
    async def _dispatch(self, client: Client, message: dict) -> None:
        msg_id = message.get("id")
        method = str(message.get("method") or "")
        params = message.get("params")
        params = {} if params is None else params

        def answer(value: Any) -> None:
            if msg_id is not None:
                client.send(protocol.result(msg_id, value))

        def fail(code: int, text: str) -> None:
            if msg_id is not None:
                client.send(protocol.error(msg_id, code, text))

        if not isinstance(params, dict):
            return fail(protocol.INVALID_PARAMS, "params must be an object")
        if method != "hello" and not client.greeted:
            return fail(protocol.NOT_ALLOWED, "hello first")
        if method in ACTING and not client.may_act:
            return fail(protocol.NOT_ALLOWED,
                        "this connection comes from an agent session's own "
                        "process: it may watch, not act (only the user "
                        "answers, sends and ends)")
        if self.closing and method in ACTING and method != "service.shutdown":
            # Its connections stay open while the sessions end and save; a
            # session started now would outlive the service.
            return fail(protocol.NOT_ALLOWED, "the agent service is stopping")
        try:
            if method in FORWARDED:
                return answer(await self._forward(client, method, params))
            if method not in SERVICE_REQUESTS:
                return fail(protocol.METHOD_NOT_FOUND, f"no method {method!r}")
            problem = protocol.params_problem(SERVICE_REQUESTS[method], params)
            if problem:
                return fail(protocol.INVALID_PARAMS, problem)
            # Package SW2b: swarm.* methods live in swarm_rpc.py, not
            # on this class -- dispatch them by module. Their handlers
            # take only params (no client, no await) and are sync, so
            # they run in a thread; every state-changing one confirms
            # inside swarm_rpc itself (swarm.start double gate).
            if method.startswith("swarm."):
                import asyncio as _asyncio
                from . import swarm_rpc as _swarm_rpc
                fn = getattr(_swarm_rpc, method.replace(".", "_"))
                return answer(await _asyncio.to_thread(fn, params))
            handler = getattr(self, "_m_" + method.replace(".", "_"))
            return answer(await handler(client, params))
        except WorkerFailed as exc:
            return fail(exc.code, exc.message)
        except _Refusal as exc:
            return fail(exc.code, str(exc))
        except Exception as exc:                    # noqa: BLE001
            return fail(protocol.INTERNAL_ERROR, f"{type(exc).__name__}: {exc}")

    async def _forward(self, client: Client, method: str, params: dict) -> Any:
        sid = str(params.get("sessionId") or "")
        rest = {k: v for k, v in params.items() if k != "sessionId"}
        if sid in self._waking and method == "session.shutdown":
            # Being woken: stop what the wake brings up, not a record.
            try:
                await asyncio.shield(self._waking[sid])
            except Exception:
                pass                     # the wake failed: whatever is left is stopped below
        if sid not in self.sessions and sid in self.dormant:
            # What a dormant session can answer without a worker; anything
            # else wakes it.
            record = self.dormant[sid]
            if method == "session.replay":
                problem = protocol.params_problem(protocol.REQUESTS[method], rest)
                if problem:
                    raise _Refusal(protocol.INVALID_PARAMS, problem)
                return {"events": await asyncio.to_thread(
                    _dormant_replay, record, int(rest.get("turns") or 0))}
            if method == "session.shutdown":
                if rest.get("dryRun"):
                    return {"closing": False, "busy": ""}     # asleep: nothing to cut
                self._forget_dormant(sid, "stopped while asleep")
                return {"closing": True}
            if method == "session.status":
                return {"state": "sleeping", "approval": record["summary"].get("approval", ""),
                        "queued": 0}
            if method == "session.setTitle":
                # Kept as the worker keeps it (its _m_session_setTitle): the
                # title changes, the address does not.
                title = protocol.clean_title(rest.get("title"))
                if not title:
                    raise _Refusal(protocol.INVALID_PARAMS, "a name needs some text")
                record["summary"]["title"] = title
                record["settings"]["title"] = title
                self._make_dormant(record)
                # Saved under it too: stopped before it wakes, the session
                # is still resumed by that name later.
                try:
                    from .. import session_store
                    await asyncio.to_thread(session_store.rename_session, sid, title)
                except Exception:
                    pass
                return {"title": title}
            if method == "session.setMode" and str(rest.get("mode") or "") in ("plan", "default", "acceptEdits"):
                record["summary"]["approval"] = str(rest["mode"])
                record["settings"]["permissionMode"] = str(rest["mode"])
                self._make_dormant(record)
                return {"approval": str(rest["mode"])}
            if method == "session.cancel":
                return {"stopping": False, "refused": 0}
        session = await self._awake(client, sid)
        if method == "session.replay":
            # From memory, as an attach: the service holds the replay
            # budget of every session -- what the worker's log would give --
            # cut where its state stands, and nothing reads a log (an agent
            # view peeks at a session on every key).
            problem = protocol.params_problem(protocol.REQUESTS[method], rest)
            if problem:
                raise _Refusal(protocol.INVALID_PARAMS, problem)
            return {"events": session.replay(int(rest.get("turns") or 0))}
        try:
            return await session.request(method, rest, timeout=FORWARD_TIMEOUT_S)
        except asyncio.TimeoutError:
            raise _Refusal(protocol.INTERNAL_ERROR,
                           f"the session did not answer {method} within "
                           f"{FORWARD_TIMEOUT_S:.0f} s") from None

    async def _m_hello(self, client: Client, params: dict) -> dict:
        version = int(params.get("protocolVersion", 0) or 0)
        if version != protocol.PROTOCOL_VERSION:
            raise _Refusal(protocol.VERSION_MISMATCH,
                           f"the service speaks protocol {protocol.PROTOCOL_VERSION},"
                           f" this window {version}")
        client.greeted = True
        client.name = str(params.get("client") or "")[:64]
        return {"protocolVersion": protocol.PROTOCOL_VERSION,
                "serviceVersion": SERVICE_VERSION, "release": _release(), "pid": os.getpid(),
                "host": rundir.hostname(), "mayAct": client.may_act,
                "sessions": self._summaries()}

    async def _m_sessions_list(self, _client: Client, _params: dict) -> dict:
        # The presence records are files, on a cluster's shared home: read
        # on a thread (an agent view asks every few seconds).
        elsewhere = await asyncio.to_thread(_elsewhere, {s.id for s in self.sessions.values()})
        return {"sessions": self._summaries(), "elsewhere": elsewhere}

    async def _m_sessions_past(self, _client: Client, params: dict) -> dict:
        limit = max(1, min(200, int(params.get("limit") or 30)))
        rows = await asyncio.to_thread(past_sessions, str(params.get("workspace") or ""),
                                       str(params.get("query") or ""), limit)
        for row in rows:
            row["running"] = row["sessionId"] in self.sessions or row["sessionId"] in self.dormant
        return {"sessions": rows}

    async def _m_session_create(self, client: Client, params: dict) -> dict:
        cwd = Path(str(params["cwd"])).expanduser()
        if not str(params["cwd"]).strip() or not cwd.is_absolute():
            # A relative folder would be taken relative to the service's
            # own working directory ($HOME), not the window's.
            raise _Refusal(protocol.INVALID_PARAMS, "cwd must be an absolute path")
        if not cwd.is_dir():
            raise _Refusal(protocol.INVALID_PARAMS, f"no such folder: {cwd}")
        for key in ("addDirs", "readDirs"):
            for path in params.get(key) or []:
                # As for cwd: a relative grant would be resolved against the
                # service's own folder ($HOME) and grant that.
                if not Path(str(path)).expanduser().is_absolute():
                    raise _Refusal(protocol.INVALID_PARAMS, f"{key} must be absolute paths: {path!r}")
        if self.closing:
            raise _Refusal(protocol.NOT_ALLOWED, "the agent service is stopping")
        resume = str(params.get("resume") or "")
        if params.get("continue") and not resume and not params.get("fork"):
            # The folder's latest session may be one this service runs or
            # keeps asleep (its window closed). One session, one worker: it
            # is resumed by id, which attaches below instead of starting a
            # second worker that would find the session's lock taken.
            latest = await asyncio.to_thread(past_sessions, str(cwd), "", 1)
            known = str(latest[0]["sessionId"]) if latest else ""
            if known and (known in self.sessions or known in self.dormant
                          or known in self._waking):
                resume = known
                if known in self.sessions:
                    session = await self._awake(client, known)
                    if params.get("attach", True) and not client.closed:
                        session.subscribers[client] = session.seen_seq
                    return {"sessionId": session.id, "session": session.info,
                            "state": session.state, "pending": list(session.pending.values()),
                            "events": session.replay(0)[1:]}
        if resume and (resume in self.dormant or resume in self._waking):
            # Woken by the resume picker, as an attach wakes it: with the
            # settings it had. What the answer carries is its replay.
            session = await self._awake(client, resume)
            if params.get("attach", True) and not client.closed:
                session.subscribers[client] = session.seen_seq
            return {"sessionId": session.id, "session": session.info, "state": session.state,
                    "pending": list(session.pending.values()), "events": session.replay(0)[1:]}
        if resume and (resume in self.sessions or any(
                str(s.settings.get("resume") or "") == resume for s in self.starting)):
            # One session, one worker: two would append to the same log.
            raise _Refusal(protocol.NOT_ALLOWED,
                           f"session {resume} runs here already; attach to it")
        self._creating += 1
        try:
            return await self._create(client, params, cwd)
        finally:
            self._creating -= 1

    async def _create(self, client: Client, params: dict, cwd: Path) -> dict:
        args = [*self.worker_cmd, "--cwd", str(cwd)]
        for key, flag in (("name", "-n"), ("title", "--title"), ("model", "--model"),
                          ("provider", "--provider"), ("effort", "--effort"),
                          ("mode", "--mode"), ("permissionMode", "--permission-mode"),
                          ("resume", "--resume")):
            if params.get(key):
                args += [flag, str(params[key])]
        for key, flag in (("continue", "--continue"), ("fork", "--fork-session"),
                          ("isolate", "--isolate"), ("noIsolate", "--no-isolate"),
                          ("fake", "--fake")):
            if params.get(key):
                args.append(flag)
        if params.get("unattended"):
            # A user schedule's run (schedule_runner): chat's ceilings.
            args.append("--unattended")
        for key, flag in (("maxBudgetUsd", "--max-budget-usd"),
                          ("maxRunSeconds", "--max-run-seconds")):
            if params.get(key):
                args += [flag, str(float(params[key]))]
        if params.get("fake") and params.get("fakeDelay") is not None:
            args += ["--fake-delay", str(float(params["fakeDelay"]))]
        if params.get("fake") and params.get("fakeSave"):
            args.append("--fake-save")
        for path in params.get("addDirs") or []:
            args += ["--add-dir", str(path)]
        for path in params.get("readDirs") or []:
            args += ["--read-dir", str(path)]
        logs = eventlog.service_dir() / "logs"
        from .. import state_paths
        state_paths.ensure_dir(logs)
        self._spawned += 1
        log_name = logs / (f"worker-{time.strftime('%Y%m%d-%H%M%S')}-"
                           f"{os.getpid()}-{self._spawned}.log")
        log_fd = os.open(log_name, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            # Started in the runtime directory (ours alone, 0700): `-m` puts
            # the working directory first on sys.path, and a json.py of the
            # user's home would be imported instead of the standard one.
            proc = await asyncio.create_subprocess_exec(
                *args, stdin=asyncio.subprocess.PIPE,
                stdout=asyncio.subprocess.PIPE, stderr=log_fd,
                start_new_session=True, limit=protocol.MAX_LINE_BYTES + 1,
                cwd=str(rundir.ensure_runtime_dir()))
        finally:
            os.close(log_fd)
        session = Session(self, proc, dict(params))
        self.starting.add(session)
        self.spawn(session.read())
        try:
            if self.closing:
                # The service began to stop while the worker started; it
                # ends with the rest (end_all_sessions takes `starting`).
                raise _Refusal(protocol.NOT_ALLOWED, "the agent service is stopping")
            hello = {"protocolVersion": protocol.PROTOCOL_VERSION,
                     "client": "delfin-agent service"}
            if str(params.get("brief") or "").strip():
                # Handed over in the initialize request -- not in argv,
                # which every user of a login node can read.
                hello["brief"] = str(params["brief"])
                hello["briefTitle"] = str(params.get("briefTitle") or "")
            init = await session.request("initialize", hello, timeout=CREATE_TIMEOUT_S)
        except (WorkerFailed, asyncio.TimeoutError) as exc:
            self.starting.discard(session)
            if proc.returncode is None:
                proc.kill()
            reason = session.ended_reason or str(exc) or "the session did not start"
            raise _Refusal(WORKER_GONE, reason) from None
        # The worker's state at its answer, then what it said after it (read
        # since, in order): the answer is older than those events.
        session.id = str(init["session"]["sessionId"])
        session.info.update(init.get("session") or {})
        session.state = str(init.get("state") or "idle")
        session.pending = {str(a.get("requestId")): a for a in init.get("pending") or []}
        for event in session._early:
            session.apply(event, count_cost=False)
        first_live = min((int(e.get("seq", 0) or 0) for e in session._early
                          if int(e.get("seq", 0) or 0) > 0), default=0)
        session._early = []
        self.starting.discard(session)
        self.sessions[session.id] = session
        self._write_roster()
        if params.get("resume") or params.get("continue") or params.get("fork"):
            # The turns of the session before this worker: for the replay.
            await session.load_history()
        # What the worker has said since it started -- in the answer, since
        # the window cannot place an event of a session it has not heard of
        # yet -- and from there on, live.
        events = [e for e in session.recent.events()
                  if first_live and int(e.get("seq", 0) or 0) >= first_live]
        if params.get("attach", True) and not client.closed:
            session.subscribers[client] = session.seen_seq
        if self.closing:
            # Registered: end_all_sessions ends it with the others.
            raise _Refusal(protocol.NOT_ALLOWED, "the agent service is stopping")
        # Every window hears of it -- after the creator has its answer.
        asyncio.get_running_loop().call_soon(self.broadcast_summary, session)
        return {"sessionId": session.id, "session": session.info,
                "state": session.state, "pending": list(session.pending.values()),
                "events": events}

    async def _m_session_attach(self, client: Client, params: dict) -> dict:
        session = await self._awake(client, str(params["sessionId"]))
        # The replay ends at the last event the service has read, where its
        # state and open questions stand, and the window is subscribed right
        # there, in the same step of the loop: every later event arrives
        # live, none twice, and what the answer says agrees with what it
        # shows.
        events = session.replay(int(params.get("turns", 0) or 0))
        session.subscribers[client] = session.seen_seq
        return {"session": session.info, "state": session.state,
                "pending": list(session.pending.values()), "events": events,
                "summary": session.summary()}

    async def _m_session_pending(self, _client: Client, params: dict) -> dict:
        # Read the current requests without attaching, replaying or waking a worker.
        session = self.sessions.get(str(params["sessionId"]))
        return {"pending": list(session.pending.values()) if session else []}

    async def _m_session_detach(self, client: Client, params: dict) -> dict:
        session = self.sessions.get(str(params["sessionId"]))
        if session is not None:
            session.subscribers.pop(client, None)
        return {"detached": True}

    async def _m_session_message(self, _client: Client, params: dict) -> dict:
        from .. import session_messages
        text = str(params["text"]).strip()
        if not text:
            raise _Refusal(protocol.INVALID_PARAMS, "empty message")
        session_messages.send(str(params["key"]), text, from_key="tui",
                              from_title="the TUI")
        return {"sent": True}

    async def _m_service_status(self, _client: Client, _params: dict) -> dict:
        def windows(clients) -> int:
            # TUI windows, not this request or other probes.
            return sum(1 for c in clients if c.greeted and c.name == "delfin-tui")

        return {"pid": os.getpid(), "host": rundir.hostname(),
                "startedAt": round(self.started_at, 3),
                # Each with the windows attached to it, and whether its
                # worker saves it (`service status`, `service end`).
                "sessions": [{**s.summary(), "windows": windows(s.subscribers),
                              "saves": _saves(s.settings)}
                             for s in self.sessions.values()],
                "dormant": [{**r["summary"], "windows": 0, "saves": _saves(r.get("settings") or {})}
                            for sid, r in self.dormant.items() if sid not in self.sessions],
                "windows": windows(self.clients),
                "socket": str(rundir.socket_path())}

    async def _m_service_shutdown(self, _client: Client, params: dict) -> dict:
        count = len(self.sessions) + len(self.starting)
        # Closing now, so a request already on its way (a create sent right
        # behind this) is refused; stopping on the next step, after this
        # answer is out.
        self.closing = True
        self.restart = bool(params.get("restart"))
        asyncio.get_running_loop().call_soon(self.stop, "restart" if self.restart else "stop")
        return {"stopping": True, "sessions": count}

    # -- what needs the user ---------------------------------------------------------------
    async def _m_attention_list(self, _client: Client, _params: dict) -> dict:
        stale = time.monotonic() - self._attention_at >= ATTENTION_POLL_S
        await self._refresh_attention(force=stale)
        return {"items": self.attention_items, "version": self.attention_version}

    async def _m_attention_answer(self, _client: Client, params: dict) -> dict:
        text = str(params["text"]).strip()
        if not text:
            raise _Refusal(protocol.INVALID_PARAMS, "an answer needs some text")
        done = await asyncio.to_thread(attention_feed.answer, str(params["id"]), text)
        await self._refresh_attention(force=True)
        return {"ok": bool(done.get("ok")), "error": str(done.get("error") or "") or None}

    async def _m_attention_dismiss(self, _client: Client, params: dict) -> dict:
        done = await asyncio.to_thread(attention_feed.dismiss, str(params["id"]))
        await self._refresh_attention(force=True)
        return {"ok": bool(done.get("ok")), "error": str(done.get("error") or "") or None}

    async def _watch_attention(self) -> None:
        """While a window is there to tell: the attention records looked at
        every ATTENTION_POLL_S, and what changed sent to every window."""
        while True:
            await asyncio.sleep(ATTENTION_POLL_S)
            if self.closing or not any(c.greeted for c in self.clients):
                continue
            try:
                await self._refresh_attention()
            except Exception:
                pass                  # a record that cannot be read now is read at the next look

    async def _refresh_attention(self, *, force: bool = False) -> None:
        """The items again when a record changed (or force, or a minute
        went by); every window hears of them when they differ."""
        async with self._attention_lock:
            # The sessions this service runs say what they need themselves;
            # one that sleeps here may be open in a dashboard by now, whose
            # questions are the list's. A session starting or ending here
            # changes what is left out: looked at again at the next poll.
            own = frozenset(self.sessions)
            fp = (await asyncio.to_thread(attention_feed.fingerprint), own)
            now = time.monotonic()
            if not force and fp == self._attention_fp and now - self._attention_at < ATTENTION_FULL_S:
                return
            items = await asyncio.to_thread(attention_feed.snapshot, own,
                                            titles=self._attention_titles)
            self._attention_fp = fp
            self._attention_at = now
            if items == self.attention_items:
                return
            self.attention_items = items
            self.attention_version += 1
            version = self.attention_version
        data = protocol.encode(protocol.event({"type": "attention.updated", "session": "",
                                               "ts": round(time.time(), 3), "items": items,
                                               "version": version}))
        for client in list(self.clients):
            if client.greeted:
                client.send_bytes(data)

    # -- sessions ---------------------------------------------------------------------------
    def broadcast_summary(self, session: Session, *, roster_now: bool = False) -> None:
        """A session's summary changed: every window hears it, and the
        roster follows -- at once when what changed is what a later
        service or launcher acts on (``roster_now``: its state, posture,
        name, a question), before any window can see it."""
        if not session.id:
            return
        if roster_now and not self.closing:
            if self._roster_timer is not None:
                self._roster_timer.cancel()
                self._roster_timer = None
            self._write_roster()
        else:
            self._roster_soon()
        message = protocol.event({"type": "session.summary", "session": session.id,
                                  "ts": round(time.time(), 3),
                                  "summary": session.summary()})
        data = protocol.encode(message)
        for client in list(self.clients):
            if client.greeted:
                client.send_bytes(data)

    def session_ended(self, session: Session) -> None:
        self.starting.discard(session)
        if session.id and session.sleeping and not self.closing \
                and self.sessions.get(session.id) is session:
            # Put to sleep: saved, without a worker, until it is used.
            self.sessions.pop(session.id, None)
            self._make_dormant(self._dormant_record(session, "asleep: unused for a while"))
            if not self.sessions:
                self.quiet_since = time.monotonic()
            return
        if session.id and self.sessions.pop(session.id, None) is not None:
            # Ended for a restart, it sleeps across it: the next service
            # offers it, and its windows reconnect to it -- for them it
            # goes on, and was never removed.
            if not (self.restart and session.id in self.dormant):
                message = protocol.event({"type": "session.removed",
                                          "session": session.id,
                                          "ts": round(time.time(), 3),
                                          "reason": session.ended_reason or _how_it_died(session.proc.returncode),
                                          "exitCode": session.proc.returncode})
                data = protocol.encode(message)
                for client in list(self.clients):
                    if client.greeted:
                        client.send_bytes(data)
            self._write_roster()
        if not self.sessions:
            self.quiet_since = time.monotonic()

    async def end_all_sessions(self) -> None:
        """Every session ends the way a window's shutdown ends it, and saves;
        what does not end in time gets SIGTERM, then SIGKILL. A session
        still starting when the service began to stop is waited for and
        ended too, so no worker outlives the service."""
        self.closing = True
        for _round in range(5):
            sessions = [*self.sessions.values(), *self.starting]
            if not sessions:
                if not self._creating:
                    return
                # A create between its checks and its worker: let it get to
                # one or the other (it sees `closing`), then look again.
                await asyncio.sleep(0.05)
                continue
            await self._end(sessions)

    async def _end(self, sessions: list[Session]) -> None:
        async def ask(session: Session) -> None:
            if not session.id:
                # Still starting (its initialize not answered): nothing of
                # it to save, and it may not read its stdin before its
                # engine is built -- it ends at once.
                try:
                    session.proc.send_signal(signal.SIGTERM)
                except ProcessLookupError:
                    pass
                return
            try:
                await session.request("session.shutdown", timeout=5)
            except Exception:
                pass
        # All at once: each may take its time, and the grace is shared.
        await asyncio.gather(*(ask(s) for s in sessions))
        await self._wait_for(sessions, SHUTDOWN_GRACE_S)
        for session in sessions:
            if session.proc.returncode is None:
                try:
                    session.proc.send_signal(signal.SIGTERM)
                except ProcessLookupError:
                    pass
        await self._wait_for(sessions, 12.0)
        for session in sessions:
            if session.proc.returncode is None:
                try:
                    session.proc.kill()
                except ProcessLookupError:
                    pass
        await self._wait_for(sessions, 5.0)

    async def _wait_for(self, sessions: list[Session], timeout: float) -> None:
        waiting = [s.done for s in sessions if not s.done.done()]
        if waiting:
            await asyncio.wait(waiting, timeout=timeout)

    def _roster_soon(self) -> None:
        """The roster again half a second on, for what changes in bursts (a
        turn's activity). What a later service -- after a crash -- or a
        launcher of another release deciding whether a restart would cut
        anything acts on is written at once (broadcast_summary)."""
        if self._roster_timer is not None or self.closing:
            return
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            return

        def write() -> None:
            self._roster_timer = None
            self._write_roster()

        self._roster_timer = loop.call_later(0.5, write)

    def _write_roster(self, *, live: bool = True) -> None:
        """Which sessions run here (with what resumes them) and which are
        dormant: for `service status`, and for a later service -- after a
        crash too -- to offer them again. ``live``: this service still runs
        (a later one must not take its sessions)."""
        try:
            path = rundir.roster_file()
            from .. import state_paths
            state_paths.ensure_dir(path.parent)
            state_paths.write_text_atomic(path, json.dumps(
                {"pid": os.getpid() if live else 0, "socket": str(rundir.socket_path()),
                 "updatedAt": round(time.time(), 3),
                 "sessions": [{**s.summary(), "settings": s.resume_settings()}
                              for s in self.sessions.values() if s.id],
                 "dormant": list(self.dormant.values())},
                indent=2))
            state_paths.secure_file(path)
        except Exception:
            pass

    def _remove_roster(self) -> None:
        """Every session ended with the service and saved: nothing for a
        later service to report (a crash leaves the roster behind)."""
        try:
            rundir.roster_file().unlink()
        except OSError:
            pass

    async def _watch_idle(self) -> None:
        assert self.stopping is not None
        limits = [self.idle_exit_s] + ([self.sleep_after_s] if self.sleep_after_s > 0 else [])
        tick = min(30.0, max(0.5, min(limits) / 4))
        while True:
            await asyncio.sleep(tick)
            self._sleep_unused()
            if self.sessions or self.starting or self.clients or self.schedules.keeps_up():
                continue
            # Dormant sessions keep nothing running: they wait in the roster.
            if time.monotonic() - self.quiet_since >= self.idle_exit_s:
                self.stop("idle")
                return

    def _sleep_unused(self) -> None:
        """Put to sleep what no window shows and has nothing to do: idle,
        no question open, nothing for SLEEP_AFTER_S."""
        if self.sleep_after_s <= 0 or self.closing:
            return
        now = time.time()
        for session in list(self.sessions.values()):
            if (session.state == "idle" and not session.pending and not session.subscribers
                    and not session.sleeping and now >= session.sleep_not_before
                    and now - session.last_activity >= self.sleep_after_s):
                session.sleeping = True
                self.spawn(self._put_to_sleep(session))

    async def _put_to_sleep(self, session: Session) -> None:
        try:
            answer = await session.request("session.shutdown", {"ifIdle": True}, timeout=30)
        except Exception:
            answer = None
        if isinstance(answer, dict) and answer.get("closing") is False:
            # Busy after all (a background shell runs): asked again after
            # another while, never cut.
            session.sleeping = False
            session.sleep_not_before = time.time() + self.sleep_after_s
            return
        try:
            await asyncio.wait_for(asyncio.shield(session.done), 60)
        except asyncio.TimeoutError:
            if session.proc.returncode is None:
                session.proc.send_signal(signal.SIGTERM)

    # -- dormant sessions -----------------------------------------------------------------
    def _summaries(self) -> list[dict]:
        return ([s.summary() for s in self.sessions.values() if s.id]
                + [dict(r["summary"]) for r in self.dormant.values()
                   if r["summary"]["sessionId"] not in self.sessions])

    def _dormant_record(self, session: Session, why: str) -> dict:
        summary = {**session.summary(), "state": "sleeping", "pending": 0, "pid": 0,
                   "activity": None, "waiting": None}
        return {"summary": summary, "settings": session.resume_settings(),
                "since": round(time.time(), 3), "why": why, "oldPid": 0}

    def _make_dormant(self, record: dict) -> None:
        sid = record["summary"]["sessionId"]
        self.dormant[sid] = record
        self._write_roster()
        message = protocol.event({"type": "session.summary", "session": sid,
                                  "ts": round(time.time(), 3), "summary": dict(record["summary"])})
        data = protocol.encode(message)
        for client in list(self.clients):
            if client.greeted:
                client.send_bytes(data)

    def _forget_dormant(self, sid: str, reason: str) -> None:
        if self.dormant.pop(sid, None) is None:
            return
        self._write_roster()
        data = protocol.encode(protocol.event({"type": "session.removed", "session": sid,
                                               "ts": round(time.time(), 3), "reason": reason}))
        for client in list(self.clients):
            if client.greeted:
                client.send_bytes(data)

    def _load_dormant(self) -> bool:
        """What a service before this one left to be woken: its dormant
        sessions, and -- when it died without saying so -- the ones it ran
        (their workers saved them when its pipes closed). False when the
        roster is another live service's (not to be taken, nor written)."""
        try:
            data = json.loads(rundir.roster_file().read_text())
        except (OSError, ValueError):
            return True
        if not isinstance(data, dict):
            return True
        try:
            old = int(data.get("pid") or 0)
        except (TypeError, ValueError):
            old = 0
        if old and old != os.getpid() and _alive_service(old):
            return False                 # another service's, alive: not ours to take
        if data.get("socket") not in (None, str(rundir.socket_path())):
            return False
        now = time.time()
        for record in data.get("dormant") or []:
            try:
                if isinstance(record, dict) and _usable_record(record) \
                        and now - float(record.get("since") or 0) < DORMANT_KEEP_S:
                    self.dormant[str(record["summary"]["sessionId"])] = record
            except (TypeError, ValueError, KeyError):
                continue                 # a damaged record costs itself only
        for entry in data.get("sessions") or []:
            try:
                if not isinstance(entry, dict) or not isinstance(entry.get("settings"), dict):
                    continue
                entry = dict(entry)
                settings = entry.pop("settings")
                record = {"summary": {**entry, "state": "sleeping", "pending": 0, "pid": 0,
                                      "activity": None, "waiting": None},
                          "settings": settings, "since": round(now, 3),
                          "why": "the agent service ended unexpectedly",
                          "oldPid": int(entry.get("pid") or 0)}
                if _usable_record(record):
                    self.dormant[str(record["summary"]["sessionId"])] = record
            except (TypeError, ValueError, KeyError):
                continue
        return True

    async def _old_worker_gone(self, record: dict) -> None:
        """A worker of a service that died saves its session as its pipe
        closes: a wake before it is done would resume what is being saved."""
        pid = int(record.get("oldPid") or 0)
        deadline = time.monotonic() + OLD_WORKER_WAIT_S
        while pid and _alive_worker(pid) and time.monotonic() < deadline:
            await asyncio.sleep(0.2)

    async def _awake(self, client: Client, sid: str) -> Session:
        """The session ``sid`` with a worker: as it runs, or woken."""
        session = self.sessions.get(sid)
        if session is not None and session.sleeping:
            # Being asked to sleep right now: it either goes (then it is
            # woken) or says it is busy and stays -- never wait for an end
            # that a refusal means will not come.
            while session.sleeping and not session.done.done():
                await asyncio.sleep(0.05)
            if session.done.done():
                await asyncio.sleep(0)       # let session_ended file it as dormant
            session = self.sessions.get(sid)
        if session is not None:
            return session
        waking = self._waking.get(sid)
        if waking is not None:
            await asyncio.shield(waking)
        elif sid in self.dormant:
            waking = asyncio.get_running_loop().create_future()
            self._waking[sid] = waking
            try:
                await self._wake(client, sid)
                waking.set_result(None)
            except BaseException as exc:
                waking.set_exception(exc)
                waking.exception()       # retrieved: no "never retrieved" warning
                raise
            finally:
                self._waking.pop(sid, None)
        session = self.sessions.get(sid)
        if session is None:
            raise _Refusal(UNKNOWN_SESSION, f"no session {sid!r} here")
        return session

    async def _wake(self, client: Client, sid: str) -> None:
        if self.closing:
            raise _Refusal(protocol.NOT_ALLOWED, "the agent service is stopping")
        if not client.may_act:
            raise _Refusal(protocol.NOT_ALLOWED, "an agent session's own process may not wake a session")
        record = self.dormant[sid]
        await self._old_worker_gone(record)
        params = {**record["settings"], "resume": sid, "attach": False}
        cwd = Path(str(params.get("cwd") or "")).expanduser()
        if not cwd.is_absolute() or not cwd.is_dir():
            raise _Refusal(protocol.INVALID_PARAMS, f"its folder is gone: {cwd}")
        self._creating += 1
        try:
            created = await self._create(client, params, cwd)
        finally:
            self._creating -= 1
        if created["sessionId"] != sid:
            # Resumed as something else (it is not in the session store any
            # more): that is no wake.
            other = self.sessions.get(created["sessionId"])
            if other is not None:
                try:
                    await other.request("session.shutdown", timeout=10)
                except Exception:
                    pass
            self._forget_dormant(sid, "it could not be resumed")
            raise _Refusal(WORKER_GONE, "the session could not be resumed: it is not saved any more")
        self.dormant.pop(sid, None)
        self._write_roster()

    def stop(self, reason: str = "") -> None:
        """Stop, from any path: a request, a signal, the emergency stop, the
        idle timer. Before serve() runs too: it stops at once then. The
        first reason given is the one windows are told."""
        if reason and not (self.stopping is not None and self.stopping.is_set()):
            self.stop_reason = reason
        self.closing = True
        if self.stopping is not None:
            self.stopping.set()


def _saves(settings: dict) -> bool:
    """Whether a session's worker saves it when it ends: every one but a
    scripted (fake) session started without fakeSave."""
    return not (settings.get("fake") and not settings.get("fakeSave"))


def _release() -> str:
    """The DELFIN release this service runs: a launcher of another one
    moves an idle service over (launcher.ensure_service)."""
    try:
        from ... import __version__
        return str(__version__)
    except Exception:
        return ""


def _usable_record(record: dict) -> bool:
    summary = record.get("summary")
    settings = record.get("settings")
    return (isinstance(summary, dict) and isinstance(settings, dict)
            and bool(summary.get("sessionId")) and bool(settings.get("cwd")))


def _dormant_replay(record: dict, turns: int) -> list[dict]:
    """A dormant session's last turns, from its log (no worker to ask)."""
    summary = record["summary"]
    sid = str(summary["sessionId"])
    header = {"type": "session.info", "session": sid, "seq": 0, "ts": round(time.time(), 3),
              **{k: summary[k] for k in ("workspace", "provider", "model", "approval", "title")
                 if summary.get(k) is not None}}
    try:
        older, cut = eventlog.read_newest(eventlog.log_path(sid), REPLAY_BUDGET_BYTES)
    except ValueError:
        return [header]
    events = eventlog.tail_turns([e for e, _ in older], turns, header=header)
    if cut and (turns <= 0 or turns > sum(1 for e, _ in older if e.get("type") == "turn.started")):
        events.insert(1, {**eventlog.OLDER_LEFT_OUT, "session": sid, "ts": round(time.time(), 3)})
    return events


def _cmdline_has(pid: int, needle: bytes) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    if not LINUX:
        return True                          # alive, and probably ours: wait
    try:
        return needle in Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return False


def _alive_service(pid: int) -> bool:
    return _cmdline_has(pid, b"delfin.agent.service.supervisor")


def _alive_worker(pid: int) -> bool:
    return _cmdline_has(pid, b"delfin.agent.service.worker")


def _how_it_died(code: int | None) -> str:
    """What ended a worker that said nothing about it."""
    if code is None:
        return "the worker ended"
    if code < 0:
        try:
            name = signal.Signals(-code).name
        except ValueError:
            name = f"signal {-code}"
        return f"the worker was killed ({name})"
    return f"the worker exited ({code})"


def past_sessions(workspace: str, query: str, limit: int) -> list[dict]:
    """Saved sessions to resume, newest first: the session store's own
    listing (scoped to ``workspace``) or its search (``query``)."""
    from .. import session_store

    if query:
        rows = session_store.find_sessions(query, limit=limit)
    else:
        rows = session_store.list_sessions(limit=limit, workspace=workspace or None)
    out = []
    for row in rows[:limit]:
        out.append({"sessionId": str(row.get("session_id") or ""),
                    "title": str(row.get("title") or ""),
                    "workspace": str(row.get("workspace") or ""),
                    "model": str(row.get("model") or ""),
                    "provider": str(row.get("provider") or ""),
                    "updatedAt": float(row.get("updated_at") or 0),
                    "messages": int(row.get("message_count") or 0),
                    "costUsd": float(row.get("cost_usd") or 0.0)})
    return out


class _Refusal(Exception):
    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="python -m delfin.agent.service.supervisor",
                                description="The agent service (one per user and host)")
    p.add_argument("--idle-exit", type=float, default=900.0,
                   help="end after this many seconds with no session and no window")
    p.add_argument("--sleep-after", type=float, default=None,
                   help="put a session no window shows to sleep after this many idle "
                        "seconds (0: never; default DELFIN_AGENT_SLEEP_AFTER_S or 3600)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    from .. import state_paths

    state_paths.scratch_state_from_environment()
    args = parse_args(argv)
    from .. import cli as _cli
    from .. import lifeline as _lifeline
    from .. import process_guard as _process_guard

    _cli.guard_principles_or_exit()
    lock = rundir.ServiceLock()
    try:
        lock.acquire()
    except rundir.ServiceBusy:
        print("the agent service already runs on this host", file=sys.stderr)
        return 0
    # This process starts the sessions that hold the model's key, and
    # answers for them: unreadable by the commands they run, and on the
    # register the emergency stop reads.
    _process_guard.protect("agent service")
    sleep_after = args.sleep_after
    if sleep_after is None:
        try:
            sleep_after = float(os.environ.get("DELFIN_AGENT_SLEEP_AFTER_S") or SLEEP_AFTER_S)
        except ValueError:
            sleep_after = SLEEP_AFTER_S
    supervisor = Supervisor(idle_exit_s=args.idle_exit, sleep_after_s=sleep_after)
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        try:
            loop.add_signal_handler(sig, supervisor.stop, f"signal {sig.name}")
        except (NotImplementedError, RuntimeError, ValueError):
            pass
    # An emergency stop (`delfin-agent stop-all`) ends the sessions -- each
    # worker watches it too -- and then this service.
    _lifeline.watch(lambda: loop.call_soon_threadsafe(supervisor.stop, "emergency stop"))
    path = rundir.socket_path()
    rundir.write_run_file({"pid": os.getpid(), "host": rundir.hostname(),
                           "socket": str(path), "startedAt": rundir.started_at(),
                           "protocolVersion": protocol.PROTOCOL_VERSION,
                           "serviceVersion": SERVICE_VERSION})
    try:
        loop.run_until_complete(supervisor.serve(path))
    finally:
        rundir.remove_run_file(os.getpid())
        lock.release()
        loop.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
