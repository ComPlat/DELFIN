"""The supervisor's side of a set of running sessions, as data.

A supervisor of many sessions (a person with a script, or an operator agent
running a wave of them) used to read their terminals: screen text matched for
"esc to interrupt" or "[y] yes", keys typed into a pane to answer, notes typed
into a prompt that only took them while it was idle. Every change to the
terminal broke it, and answering an approval by keystroke from another agent
went around the rule that a person answers (#202, #207).

This module reads what DELFIN already records and hands it over as events and
state:

* sessions that are open (``session_presence``), their branch and activity;
* turns that ended (``turn_metrics``), steps that were refused and commits
  that were made (the audit log);
* approvals that are waiting and that closed (the terminal and file
  approval rooms);
* messages left in session inboxes (``session_messages``).

And it lets a supervisor act in the ways that cannot widen anything: leave a
message (shown as from the operator, never as the user), put a waiting
approval in front of the person (attention inbox), pause a session. It never
answers an approval: there is no function here that writes an answer.

Events are derived from those records by a reader that keeps a cursor on
disk, so a reader that restarts with the same cursor neither repeats nor
drops an event.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Iterator, Optional

SCHEMA = "delfin.ops_event/1"
OPERATOR_KEY = "operator"
OPERATOR_TITLE = "operator (supervisor of this wave, not the user)"
TEXT_MAX = 2000
DEFAULT_IDLE_S = 900
_CURSOR_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
# ``git commit`` as one segment of a command line: alone, or after && ; ||.
_COMMIT_RE = re.compile(r"(?:^|&&|;|\|\|)\s*git\s+(?:-C\s+\S+\s+)?commit\b")


# --------------------------------------------------------------------------
# where things are (resolved per call, so a redirected home redirects them)
# --------------------------------------------------------------------------

def _home() -> Path:
    return Path.home() / ".delfin"


def _cursor_dir() -> Path:
    return _home() / "ops_cursors"


def _sessions_dir() -> Path:
    try:
        from . import session_store
        return Path(session_store._SESSIONS_DIR)
    except Exception:
        return _home() / "agent_sessions"


def _turn_metrics_dir() -> Path:
    try:
        from . import turn_metrics
        return Path(turn_metrics._DIR)
    except Exception:
        return _home() / "turn_metrics"


def _inbox_dir() -> Path:
    try:
        from . import session_messages
        return Path(session_messages._DIR)
    except Exception:
        return _home() / "session_inbox"


def _bounded(text: Any, limit: int = TEXT_MAX) -> str:
    """A field value for output: a string, secrets scrubbed, length capped."""
    s = str(text if text is not None else "")
    try:
        from .output_guard import scrub_secrets
        s = scrub_secrets(s)
    except Exception:
        pass
    return s[:limit]


# --------------------------------------------------------------------------
# reads
# --------------------------------------------------------------------------

_OTHER_HOST_WINDOW_S = 6 * 3600


def _open_sessions() -> list[dict]:
    """Open sessions by their process, not by their heartbeat.

    The presence heartbeat goes quiet while a session waits at an approval
    or runs one long tool call; read by heartbeat, such a session was
    reported as exited and opened again. On this host a session is open
    while its process exists; on another host while its record is recent.
    Reads only: the presence module's own listing also reaps old records.
    """
    import socket
    try:
        from . import session_presence
        files = sorted(Path(session_presence._DIR).glob("*.json"))
    except Exception:
        return []
    host = socket.gethostname()
    now = time.time()
    out = []
    for f in files:
        try:
            rec = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(rec, dict) or not rec.get("key"):
            continue
        if rec.get("host") == host:
            try:
                pid = int(rec.get("pid") or 0)
            except (TypeError, ValueError):
                continue
            if pid <= 0:
                continue
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                continue
            except OSError:
                pass                     # another user's process: alive
        elif now - float(rec.get("updated_at") or 0) > _OTHER_HOST_WINDOW_S:
            continue
        out.append(rec)
    return out


def pending_approvals() -> list[dict]:
    """Every approval any session is waiting on, one row per id."""
    rows: dict[str, dict] = {}
    try:
        from . import terminal_confirm
        for r in terminal_confirm.pending_for_supervisors():
            if r.get("id"):
                rows[str(r["id"])] = r
    except Exception:
        pass
    try:
        from . import file_confirm
        for r in file_confirm.pending():
            if r.get("id"):
                rows.setdefault(str(r["id"]), r)
    except Exception:
        pass
    # A request names its session by id; report the name `ops` takes.
    names = {str(s.get("session_id") or ""): str(s.get("key") or "")
             for s in _open_sessions()} if rows else {}
    out = []
    for rid, r in rows.items():
        sid = str(r.get("session_id") or "")
        out.append({
            "id": rid,
            "session": str(r.get("session_key") or names.get(sid) or sid),
            "session_id": str(r.get("session_id") or ""),
            "tool": str(r.get("tool") or ""),
            "summary": _bounded(r.get("command") or r.get("path")
                                or r.get("preview") or "", 300),
            "opened_at": float(r.get("asked_at") or 0.0),
        })
    return sorted(out, key=lambda x: x["opened_at"])


def _last_turn(session_id: str) -> Optional[dict]:
    if not session_id:
        return None
    path = _turn_metrics_dir() / f"{session_id}.jsonl"
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - 8192))
            tail = fh.read().decode("utf-8", "replace").splitlines()
        for line in reversed(tail):
            try:
                rec = json.loads(line)
            except ValueError:
                continue
            if isinstance(rec, dict):
                return rec
    except OSError:
        return None
    return None


def _git(workspace: str, *args: str) -> str:
    if not workspace or not Path(workspace).is_dir():
        return ""
    try:
        out = subprocess.run(["git", "-C", workspace, *args],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def _activity(session_id: str) -> float:
    """The newest trace of work: a turn ended, a mid-turn checkpoint, a save."""
    times = [0.0]
    last = _last_turn(session_id)
    if last:
        try:
            times.append(float(last.get("ts") or 0.0))
        except (TypeError, ValueError):
            pass
    for name in (f"{session_id}.turn.json", f"{session_id}.json"):
        try:
            times.append((_sessions_dir() / name).stat().st_mtime)
        except OSError:
            pass
    return max(times)


def _inbox_count(key: str) -> int:
    try:
        from . import session_messages
        path = session_messages._inbox(key)
        with open(path, "rb") as fh:
            return sum(1 for line in fh if line.strip())
    except Exception:
        return 0


def sessions(*, with_git: bool = True) -> list[dict]:
    """One row per open session on this host, newest activity first."""
    approvals = pending_approvals()
    now = time.time()
    out = []
    for rec in _open_sessions():
        sid = str(rec.get("session_id") or "")
        key = str(rec.get("key") or "")
        ws = str(rec.get("workspace") or "")
        mine = [a for a in approvals
                if a["session"] in (key, sid) or a["session_id"] == sid]
        last = _last_turn(sid) or {}
        activity = _activity(sid)
        row = {
            "name": key,
            "session_id": sid,
            "title": _bounded(rec.get("title") or "", 80),
            "workspace": ws,
            "branch": str(rec.get("branch") or ""),
            "pid": rec.get("pid"),
            "host": rec.get("host"),
            "open_approvals": [a["id"] for a in mine],
            "last_turn_end": last.get("ts"),
            "last_turn": {k: last.get(k) for k in
                          ("stopped", "error", "tool_calls", "model")
                          if k in last},
            "idle_s": int(now - activity) if activity else None,
            "inbox": _inbox_count(key),
        }
        if with_git:
            row["head"] = _git(ws, "rev-parse", "--short", "HEAD")
            dirty = _git(ws, "status", "--porcelain", "--untracked-files=no")
            row["dirty_files"] = len(dirty.splitlines()) if dirty else 0
        out.append(row)
    out.sort(key=lambda r: r["idle_s"] if r["idle_s"] is not None else 10**12)
    return out


def session(name: str) -> Optional[dict]:
    name = str(name or "").strip()
    for row in sessions():
        if name and name in (row["name"], row["session_id"],
                             row["session_id"][:8]):
            return row
    return None


# --------------------------------------------------------------------------
# events
# --------------------------------------------------------------------------

def _cursor_path(name: str) -> Path:
    if not _CURSOR_NAME_RE.match(name or ""):
        raise ValueError(f"invalid cursor name {name!r}")
    return _cursor_dir() / f"{name}.json"


def load_cursor(name: str = "default") -> dict:
    try:
        data = json.loads(_cursor_path(name).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def save_cursor(state: dict, name: str = "default") -> None:
    path = _cursor_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(state), encoding="utf-8")
    os.chmod(tmp, 0o600)
    os.replace(tmp, path)


def _read_new_lines(path: Path, pos: dict) -> list[str]:
    """Lines appended to *path* since the remembered (inode, offset)."""
    try:
        st = path.stat()
    except OSError:
        return []
    key = str(path)
    ino, off = pos.get(key, [None, 0])
    if ino != st.st_ino or st.st_size < off:
        off = 0                          # replaced or truncated: start over
    if st.st_size == off:
        pos[key] = [st.st_ino, off]
        return []
    with open(path, "rb") as fh:
        fh.seek(off)
        data = fh.read()
    end = data.rfind(b"\n")
    if end < 0:
        pos[key] = [st.st_ino, off]
        return []
    pos[key] = [st.st_ino, off + end + 1]
    return data[:end].decode("utf-8", "replace").splitlines()


def _receipts() -> dict:
    try:
        from . import session_messages
        data = session_messages._receipts()
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _inbox_messages() -> dict:
    """Messages still waiting in any inbox, by id (inboxes are capped)."""
    out: dict = {}
    try:
        files = sorted(_inbox_dir().glob("*.jsonl"))
    except OSError:
        return out
    for f in files:
        try:
            with open(f, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    try:
                        m = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(m, dict) and m.get("id"):
                        m.setdefault("to", f.stem)
                        out[str(m["id"])] = m
        except OSError:
            continue
    return out


def _first_run_positions(state: dict) -> None:
    """A new cursor starts at 'now': history is not replayed as events."""
    pos: dict = {}
    for directory, pattern in ((_turn_metrics_dir(), "*.jsonl"),):
        try:
            for p in directory.glob(pattern):
                if p.name.endswith(".requests.jsonl"):
                    continue
                st = p.stat()
                pos[str(p)] = [st.st_ino, st.st_size]
        except OSError:
            pass
    state["pos"] = pos
    try:
        from . import audit_log
        _rows, state["audit_off"] = audit_log.steps_since(0)
    except Exception:
        state["audit_off"] = 0
    state["sessions"] = sorted(r.get("key", "") for r in _open_sessions())
    state["approvals"] = {a["id"]: {"session": a["session"], "tool": a["tool"]}
                          for a in pending_approvals()}
    state["messages"] = sorted(set(_receipts()) | set(_inbox_messages()))
    state["since"] = time.time()
    state["idle"] = []
    state["started"] = True


def poll(state: dict, *, idle_s: int = DEFAULT_IDLE_S) -> list[dict]:
    """Events since *state*; updates *state* in place. Never answers anything."""
    if not state.get("started"):
        state.setdefault("seq", 0)
        _first_run_positions(state)
        return []
    events: list[dict] = []
    pos = state.setdefault("pos", {})
    now = time.time()

    open_rows = _open_sessions()
    by_id = {str(r.get("session_id") or ""): str(r.get("key") or "")
             for r in open_rows}
    keys = sorted(str(r.get("key") or "") for r in open_rows)
    before = set(state.get("sessions") or [])
    for k in keys:
        if k not in before:
            events.append({"kind": "session_opened", "session": k})
    for k in sorted(before - set(keys)):
        events.append({"kind": "session_exited", "session": k})
    state["sessions"] = keys

    approvals = {a["id"]: a for a in pending_approvals()}
    seen = state.get("approvals") or {}
    if isinstance(seen, list):
        seen = {rid: {} for rid in seen}
    for rid, a in approvals.items():
        if rid not in seen:
            events.append({"kind": "approval_open", **a})
    for rid in sorted(set(seen) - set(approvals)):
        events.append({"kind": "approval_closed", "id": rid,
                       "session": (seen.get(rid) or {}).get("session", ""),
                       "tool": (seen.get(rid) or {}).get("tool", "")})
    state["approvals"] = {rid: {"session": a["session"], "tool": a["tool"]}
                          for rid, a in approvals.items()}

    try:
        files = sorted(_turn_metrics_dir().glob("*.jsonl"))
    except OSError:
        files = []
    for p in files:
        if p.name.endswith(".requests.jsonl"):
            continue
        sid = p.name[:-len(".jsonl")]
        for line in _read_new_lines(p, pos):
            try:
                rec = json.loads(line)
            except ValueError:
                continue
            end = ("stopped" if rec.get("stopped") else
                   "error" if rec.get("error") else "answer")
            events.append({
                "kind": "turn_ended", "session": by_id.get(sid, ""),
                "session_id": sid, "end": end,
                "tool_calls": rec.get("tool_calls"),
                "error": _bounded(rec.get("error") or "", 300),
                "turn_ts": rec.get("ts"),
            })

    # Messages: those still waiting in an inbox, and those already taken
    # (a receiver takes a message within seconds and its inbox line is gone;
    # the receipt, written on taking, keeps from/to/sent_at but no text).
    found = _inbox_messages()
    for mid, r in _receipts().items():
        if isinstance(r, dict) and mid not in found:
            found[mid] = {"id": mid, **r}
    seen_msgs = set(state.get("messages") or [])
    since = float(state.get("since") or 0.0)
    new = [m for mid, m in found.items()
           if mid not in seen_msgs and float(m.get("sent_at") or 0) >= since]
    for m in sorted(new, key=lambda m: float(m.get("sent_at") or 0)):
        events.append({
            "kind": "message", "id": m["id"],
            "from": _bounded(m.get("from") or "", 80),
            "to": _bounded(m.get("to") or "", 80),
            "sent_at": m.get("sent_at"), "taken_at": m.get("taken_at"),
            "text": _bounded(m.get("text") or ""),
        })
        seen_msgs.add(m["id"])
    state["messages"] = sorted(m for m in seen_msgs if m in found)

    try:
        from . import audit_log
        rows, state["audit_off"] = audit_log.steps_since(
            int(state.get("audit_off") or 0))
    except Exception:
        rows = []
    # Work the turn records do not show yet: a tool step, a message taken.
    busy = dict(state.get("busy") or {})
    for row in rows:
        if row.get("session_id"):
            busy[str(row["session_id"])] = now
    for m in found.values():
        to, taken = str(m.get("to") or ""), m.get("taken_at")
        if to and taken:
            busy[to] = max(float(busy.get(to) or 0), float(taken))
    for row in rows:
        sid = str(row.get("session_id") or "")
        decision = str(row.get("decision") or "")
        command = str(row.get("command") or "")
        if decision and decision != "ok":
            events.append({
                "kind": "denied", "session": by_id.get(sid, ""),
                "session_id": sid, "tool": row.get("tool"),
                "decision": decision, "detail": _bounded(command, 300),
            })
        elif (_COMMIT_RE.search(command)
              and row.get("exit_code") == 0):
            cwd = str(row.get("cwd") or row.get("workspace") or "")
            events.append({
                "kind": "commit", "session": by_id.get(sid, ""),
                "session_id": sid,
                "branch": _git(cwd, "rev-parse", "--abbrev-ref", "HEAD"),
                "head": _git(cwd, "log", "-1", "--format=%h %s")[:200],
            })

    idle_before = set(state.get("idle") or [])
    idle_now = set()
    waiting = {a["session"] for a in approvals.values()} | {
        a["session_id"] for a in approvals.values()}
    for r in open_rows:
        sid = str(r.get("session_id") or "")
        key = str(r.get("key") or "")
        activity = max(_activity(sid), float(busy.get(sid) or 0),
                       float(busy.get(key) or 0))
        if (activity and now - activity >= idle_s
                and key not in waiting and sid not in waiting):
            idle_now.add(key)
            if key not in idle_before:
                events.append({"kind": "session_idle", "session": key,
                               "session_id": sid,
                               "idle_s": int(now - activity)})
    state["idle"] = sorted(idle_now)
    live = {str(r.get("session_id") or "") for r in open_rows} | {
        str(r.get("key") or "") for r in open_rows}
    state["busy"] = {k: v for k, v in busy.items() if k in live}

    for ev in events:
        state["seq"] = int(state.get("seq") or 0) + 1
        ev.update({"schema": SCHEMA, "seq": state["seq"], "ts": now})
    return events


def matches(event: dict, spec: str) -> bool:
    """``kind=approval_open,session=s11,idle_s>=1800``: every term must hold."""
    for term in [t.strip() for t in str(spec or "").split(",") if t.strip()]:
        m = re.match(r"^([A-Za-z_]+)\s*(>=|<=|=)\s*(.*)$", term)
        if not m:
            return False
        key, op, want = m.groups()
        have = event.get(key)
        if op == "=":
            if str(have if have is not None else "") != want:
                return False
        else:
            try:
                a, b = float(have), float(want)
            except (TypeError, ValueError):
                return False
            if (op == ">=" and a < b) or (op == "<=" and a > b):
                return False
    return True


def follow(cursor: str = "default", *, idle_s: int = DEFAULT_IDLE_S,
           interval: float = 5.0, seconds: float = 0.0) -> Iterator[dict]:
    """Events as they happen; the cursor is saved after every batch."""
    state = load_cursor(cursor)
    deadline = time.monotonic() + (seconds if seconds > 0 else 1e12)
    while True:
        batch = poll(state, idle_s=idle_s)
        save_cursor(state, cursor)
        for ev in batch:
            yield ev
        if time.monotonic() >= deadline:
            return
        time.sleep(interval)


def wait(spec: str, *, cursor: str = "default", timeout: float = 0.0,
         idle_s: int = DEFAULT_IDLE_S, interval: float = 5.0
         ) -> Optional[dict]:
    """The first event that matches *spec*, or None when *timeout* passes.

    Events that do not match are consumed: a waiter that wants every event
    reads ``follow`` instead.
    """
    for ev in follow(cursor, idle_s=idle_s, interval=interval,
                     seconds=timeout):
        if matches(ev, spec):
            return ev
    return None


# --------------------------------------------------------------------------
# actions that cannot widen anything
# --------------------------------------------------------------------------

class OpsRefused(RuntimeError):
    """An ops action that may not run from here."""


def _refuse_inside_agent() -> None:
    """A session must not speak as the operator or pause another session."""
    try:
        from . import approval_guard
        if approval_guard.inside_agent_session():
            raise OpsRefused("ops actions are run by the supervisor, not "
                             "from inside an agent session")
    except ImportError:
        pass


def send(name: str, text: str) -> dict:
    """Leave *text* in session *name*'s inbox, from the operator.

    The receiving session reads it as a message from another session, never
    from the user; it does not count as a typed message and lifts nothing.
    """
    _refuse_inside_agent()
    name = str(name or "").strip()
    if not name or name == OPERATOR_KEY:
        raise OpsRefused("send needs a session name")
    from . import session_messages
    return session_messages.send(name, str(text or ""),
                                 from_key=OPERATOR_KEY,
                                 from_title=OPERATOR_TITLE)


def escalate(approval_id: str, note: str = "") -> str:
    """Put a waiting approval in front of the person. It stays open."""
    _refuse_inside_agent()
    rows = {a["id"]: a for a in pending_approvals()}
    a = rows.get(str(approval_id))
    if a is None:
        raise OpsRefused(f"nothing waiting with id {approval_id!r}")
    from . import attention
    return attention.emit_attention(
        "confirm_pending",
        session_id=a["session_id"],
        title=_bounded(f"{a['session'] or a['session_id']} waits: "
                       f"{a['tool']} {a['summary']}", 200),
        detail=_bounded(("Operator note: " + note) if note else
                        "The operator passed this approval on to you.", 800),
    )


def pause(name: str, reason: str) -> dict:
    """Pause session *name* (it stops at its next safe point)."""
    _refuse_inside_agent()
    name = str(name or "").strip()
    if not name:
        raise OpsRefused("pause needs a session name")
    from . import session_pause
    return session_pause.pause(name, _bounded(reason or "operator", 200))
