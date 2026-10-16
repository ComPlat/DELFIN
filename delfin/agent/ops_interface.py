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
drops an event it was handed (``follow`` saves the unread rest of a batch
with the cursor).

Limits. Events are what one poll sees against the last: an approval that
opened and closed between two polls, and a message its recipient took
between two polls (``take`` unlinks the inbox), never appear. The
actions refuse when the process carries the agent-session marker
(``approval_guard.ENV_MARKER``); an agent shell that strips it can still
send as the operator, which is the limit ``approval_guard`` names for
its own refusals -- such a message still reaches the session under the
session-message header, below the user's authority. Every file read
here is opened without following a link, non-blocking, and only when it
is a regular file of this account; git runs only in a workspace of this
account and never compares its working tree (``_GIT_SAFE``).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import stat
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
#: A session id as it names files (turn_metrics/<id>.jsonl): no separator,
#: no "..". A presence record is written by a session; its id is not
#: trusted to be a file name.
_SESSION_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,128}$")
#: Bytes one poll reads from one record file. A record file appended by
#: megabytes between two polls is read over several polls, and one poll's
#: memory does not depend on what a session wrote.
_READ_MAX = 1 << 20
#: Bytes of a file's head whose hash tells it from a later file that the
#: filesystem gave the same inode number (take() unlinks an inbox and the
#: next send creates it again).
_HEAD_BYTES = 4096
#: A filter is at most this long and has at most this many terms.
_SPEC_MAX = 1000
_SPEC_TERMS = 16
_TERM_RE = re.compile(r"^([A-Za-z_]{1,32})\s*(>=|<=|=)\s*(.*)$")
#: Unread events a cursor keeps (the rest of a batch a waiter stopped in).
_QUEUE_MAX = 10000
#: Message files read by `send --file`, at most (the inbox keeps 4000
#: characters of a message).
FILE_MAX = 64 * 1024


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
    """A field value for output: a string, secrets scrubbed, terminal
    control removed, length capped."""
    s = str(text if text is not None else "")[:limit * 4]
    try:
        from .output_guard import scrub_secrets
        s = scrub_secrets(s)
    except Exception:
        pass
    from .session_messages import _plain
    return _plain(s)[:limit]


def dumps(obj: Any) -> str:
    """One JSON line for output, ASCII only.

    ASCII because a record may hold a lone surrogate (JSON allows
    ``\\ud800``) that a UTF-8 stdout refuses: ``print`` raised
    UnicodeEncodeError in the middle of a batch, and a reader whose cursor
    is saved after printing met the same event on every run. Escaped,
    every line is valid JSON whatever the text, and no C1 control reaches
    a terminal raw.
    """
    return json.dumps(obj, ensure_ascii=True, default=str)


def _my_uid() -> int:
    getuid = getattr(os, "getuid", None)
    return getuid() if getuid else -1


def _open_record(path: Path) -> Optional[int]:
    """A read-only descriptor for *path*, or None.

    None unless it is a regular file of this account opened without
    following a final symlink: a FIFO named like a session's record
    blocked ``open`` (and the supervisor with it), a link sent the reader
    to any file the account can read, and a file of another account in
    a directory a mount left group-writable is no record of this one.
    """
    flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
             | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_CLOEXEC", 0))
    try:
        fd = os.open(str(path), flags)
    except OSError:
        return None
    try:
        st = os.fstat(fd)
        uid = _my_uid()
        if stat.S_ISREG(st.st_mode) and (uid < 0 or st.st_uid == uid):
            return fd
    except OSError:
        pass
    os.close(fd)
    return None


def _read_at(fd: int, offset: int, size: int) -> bytes:
    chunks, done = [], 0
    while done < size:
        chunk = os.pread(fd, size - done, offset + done)
        if not chunk:
            break
        chunks.append(chunk)
        done += len(chunk)
    return b"".join(chunks)


def _session_file(directory: Path, session_id: str, suffix: str
                  ) -> Optional[Path]:
    if not _SESSION_ID_RE.match(session_id or ""):
        return None
    return directory / f"{session_id}{suffix}"


# --------------------------------------------------------------------------
# reads
# --------------------------------------------------------------------------

def _open_sessions() -> list[dict]:
    try:
        from . import session_presence
        return list(session_presence.open_sessions())
    except Exception:
        return []


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
    out = []
    for rid, r in rows.items():
        out.append({
            "id": rid,
            "session": _bounded(r.get("session_key") or r.get("session_id")
                                or "", 128),
            "session_id": _bounded(r.get("session_id") or "", 128),
            "tool": _bounded(r.get("tool") or "", 80),
            "summary": _bounded(r.get("command") or r.get("path")
                                or r.get("preview") or "", 300),
            "opened_at": _float_or_none(r.get("asked_at")) or 0.0,
        })
    return sorted(out, key=lambda x: x["opened_at"])


def _last_turn(session_id: str) -> Optional[dict]:
    path = _session_file(_turn_metrics_dir(), session_id, ".jsonl")
    if path is None:
        return None
    fd = _open_record(path)
    if fd is None:
        return None
    try:
        size = os.fstat(fd).st_size
        start = max(0, size - 65536)
        tail = _read_at(fd, start, size - start)
    except OSError:
        return None
    finally:
        os.close(fd)
    lines = tail.decode("utf-8", "replace").splitlines()
    if start:
        lines = lines[1:]                # the first one is cut
    for line in reversed(lines):
        try:
            rec = json.loads(line)
        except ValueError:
            continue
        if isinstance(rec, dict):
            return rec
    return None


#: Read-only git in a directory a session works in. Its ``.git/config``
#: is writable by the session's commands, and this process runs outside
#: every sandbox: ``git status`` ran the configured core.fsmonitor and
#: clean filter (measured), and so did ``diff-index``/``diff-files``
#: for an entry the index holds as racily clean. Only commands that
#: never compare the working tree run here (rev-parse, log, ls-files
#: without -m); the -c settings outrank the repository's own; no
#: optional locks, so the supervisor never holds the session's
#: index.lock; no GIT_* variable of this process redirects it.
_GIT_SAFE = ("--no-pager", "--no-optional-locks",
             "-c", "core.fsmonitor=false", "-c", "core.hooksPath=/dev/null",
             "-c", "log.showSignature=false", "-c", "core.pager=cat")


def _workspace_ok(workspace: str) -> bool:
    """A directory of this account: git runs in nothing else."""
    try:
        if not workspace or not os.path.isdir(workspace):
            return False
        uid = _my_uid()
        return uid < 0 or os.stat(workspace).st_uid == uid
    except OSError:
        return False


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env["GIT_TERMINAL_PROMPT"] = "0"
    return env


def _git(workspace: str, *args: str) -> str:
    if not _workspace_ok(workspace):
        return ""
    try:
        out = subprocess.run(["git", *_GIT_SAFE, "-C", workspace, *args],
                             capture_output=True, text=True, timeout=10,
                             env=_git_env(), stdin=subprocess.DEVNULL)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


_DEBUG_LINE_RE = re.compile(rb"^  (mtime|size): (\d+)(?::(\d+))?")
#: Index entries compared at most; a larger checkout reports None.
_DIRTY_MAX_ENTRIES = 50000


def _dirty_count(workspace: str) -> Optional[int]:
    """Tracked files whose size or mtime differ from what the index
    recorded, or None when that cannot be told.

    Compared here from ``git ls-files --debug -z`` (which prints the
    index's stat data and compares nothing) rather than by git: every git
    command that compares the working tree may run the repository's clean
    filter (see _GIT_SAFE). A file touched without a change counts, a
    change that kept size and mtime does not (git's racy case).
    """
    if not _workspace_ok(workspace):
        return None
    try:
        out = subprocess.run(
            ["git", *_GIT_SAFE, "-C", workspace, "ls-files", "--debug", "-z"],
            capture_output=True, timeout=10, stdin=subprocess.DEVNULL,
            env=_git_env())
    except Exception:
        return None
    if out.returncode != 0:
        return None
    # ls-files lists what is under the directory it runs in, with paths
    # relative to it.
    segments = out.stdout.split(b"\0")
    if len(segments) - 1 > _DIRTY_MAX_ENTRIES:
        return None
    entries: list[tuple[bytes, int, int, int]] = []
    path = segments[0] if segments else b""
    for seg in segments[1:]:
        lines = seg.split(b"\n")
        mtime_s = mtime_ns = size = -1
        rest_at = 0
        for i, line in enumerate(lines):
            if not line.startswith(b"  "):
                rest_at = i
                break
            m = _DEBUG_LINE_RE.match(line)
            if m and m.group(1) == b"mtime":
                mtime_s, mtime_ns = int(m.group(2)), int(m.group(3) or 0)
            elif m and m.group(1) == b"size":
                size = int(m.group(2))
            rest_at = i + 1
        if path and size >= 0:
            entries.append((path, mtime_s, mtime_ns, size))
        path = b"\n".join(lines[rest_at:])
    dirty = 0
    base = os.fsencode(workspace)
    for rel, mtime_s, mtime_ns, size in entries:
        try:
            st = os.lstat(os.path.join(base, rel))
        except OSError:
            dirty += 1
            continue
        if stat.S_ISDIR(st.st_mode):
            continue                      # a submodule's gitlink
        same_time = (int(st.st_mtime) == mtime_s
                     and (mtime_ns == 0 or st.st_mtime_ns % 10**9 == mtime_ns))
        # The index keeps 32 bits of each (st_size and st_mtime mod 2**32).
        same_time = same_time or (
            (int(st.st_mtime) & 0xFFFFFFFF) == mtime_s and (
                mtime_ns == 0 or st.st_mtime_ns % 10**9 == mtime_ns))
        if (st.st_size & 0xFFFFFFFF) != size or not same_time:
            dirty += 1
    return dirty


def _activity(session_id: str) -> float:
    """The newest trace of work: a turn ended, a mid-turn checkpoint, a save."""
    times = [0.0]
    last = _last_turn(session_id)
    if last:
        try:
            times.append(float(last.get("ts") or 0.0))
        except (TypeError, ValueError):
            pass
    for suffix in (".turn.json", ".json"):
        path = _session_file(_sessions_dir(), session_id, suffix)
        try:
            if path is not None:
                times.append(path.lstat().st_mtime)
        except OSError:
            pass
    return max(times)


def _inbox_count(key: str) -> int:
    try:
        from . import session_messages
        fd = _open_record(session_messages._inbox(key))
    except Exception:
        return 0
    if fd is None:
        return 0
    try:
        data = _read_at(fd, 0, _READ_MAX)
    except OSError:
        return 0
    finally:
        os.close(fd)
    return sum(1 for line in data.splitlines() if line.strip())


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
            "name": _bounded(key, 80),
            "session_id": _bounded(sid, 128),
            "title": _bounded(rec.get("title") or "", 80),
            "workspace": _bounded(ws, 1000),
            "branch": _bounded(rec.get("branch") or "", 200),
            "pid": _int_or_none(rec.get("pid")),
            "host": _bounded(rec.get("host") or "", 200),
            "open_approvals": [a["id"] for a in mine],
            "last_turn_end": _float_or_none(last.get("ts")),
            "last_turn": _turn_summary(last),
            "idle_s": int(now - activity) if activity else None,
            "inbox": _inbox_count(key),
        }
        if with_git:
            row["head"] = _git(ws, "rev-parse", "--short", "HEAD")[:40]
            row["dirty_files"] = _dirty_count(ws)
        out.append(row)
    out.sort(key=lambda r: r["idle_s"] if r["idle_s"] is not None else 10**12)
    return out


def _int_or_none(value: Any) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float_or_none(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _turn_summary(last: dict) -> dict:
    """The fields of a turn record a supervisor reads, bounded."""
    out: dict = {}
    if "stopped" in last:
        out["stopped"] = bool(last.get("stopped"))
    if "error" in last:
        out["error"] = _bounded(last.get("error") or "", 300)
    if "tool_calls" in last:
        out["tool_calls"] = _int_or_none(last.get("tool_calls"))
    if "model" in last:
        out["model"] = _bounded(last.get("model") or "", 120)
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
    path = _cursor_path(name)
    fd = _open_record(path)
    if fd is None:
        return {}
    try:
        data = json.loads(_read_at(fd, 0, 64 * _READ_MAX).decode("utf-8"))
    except (OSError, ValueError):
        return {}
    finally:
        os.close(fd)
    return data if isinstance(data, dict) else {}


def save_cursor(state: dict, name: str = "default") -> None:
    """Write the cursor whole and owner-only.

    ``state_paths.write_text_atomic``: a fresh mkstemp file (0600 from
    creation, never an existing path, so no planted link is followed),
    then ``os.replace``; the directory is made 0700.
    """
    from .state_paths import ensure_dir, write_text_atomic
    path = _cursor_path(name)
    ensure_dir(path.parent)
    write_text_atomic(path, json.dumps(state))


def _head_digest(fd: int, length: int) -> str:
    return hashlib.sha256(_read_at(fd, 0, length)).hexdigest()[:32]


def _mark(fd: int, st: os.stat_result, offset: int) -> list:
    """The remembered position: [inode, offset, head length, head hash]."""
    head = min(offset, _HEAD_BYTES)
    return [st.st_ino, offset, head, _head_digest(fd, head) if head else ""]


def _resume_offset(fd: int, st: os.stat_result, entry: Any) -> int:
    """Where to go on reading: the remembered offset if this is still the
    file it was taken in, else 0.

    The same file means the same inode, not shorter than the offset, and
    the same head bytes. The head check is what a reused inode number
    needs: take() unlinks an inbox, the next send creates a new file, and
    the filesystem may give it the old number -- its first message was
    skipped as already read when it was as long as the old one.
    """
    if not isinstance(entry, list) or len(entry) < 2:
        return 0
    try:
        ino, off = entry[0], int(entry[1])
    except (TypeError, ValueError):
        return 0
    if ino != st.st_ino or off < 0 or st.st_size < off:
        return 0
    if len(entry) >= 4 and entry[2]:
        try:
            head = int(entry[2])
        except (TypeError, ValueError):
            return 0
        if head > off or _head_digest(fd, head) != entry[3]:
            return 0
    return off


def _read_new_lines(path: Path, pos: dict) -> list[str]:
    """Complete lines appended to *path* since the remembered position.

    At most _READ_MAX bytes per call; the rest is read by the next one.
    A line longer than that is skipped whole rather than held: a reader
    that waited for its newline would never read past it.
    """
    key = str(path)
    fd = _open_record(path)
    if fd is None:
        pos.pop(key, None)
        return []
    try:
        st = os.fstat(fd)
        off = _resume_offset(fd, st, pos.get(key))
        if st.st_size == off:
            pos[key] = _mark(fd, st, off)
            return []
        data = _read_at(fd, off, min(_READ_MAX, st.st_size - off))
        end = data.rfind(b"\n")
        if end < 0:
            if len(data) < _READ_MAX:
                pos[key] = _mark(fd, st, off)       # a line being written
                return []
            # One line longer than a whole read: find its end and skip it.
            skip = off + len(data)
            while skip < st.st_size:
                more = _read_at(fd, skip, min(_READ_MAX, st.st_size - skip))
                nl = more.find(b"\n")
                if nl >= 0:
                    pos[key] = _mark(fd, st, skip + nl + 1)
                    return []
                skip += len(more)
                if not more:
                    break
            pos[key] = _mark(fd, st, off)
            return []
        pos[key] = _mark(fd, st, off + end + 1)
    except OSError:
        return []
    finally:
        os.close(fd)
    return data[:end].decode("utf-8", "replace").splitlines()


def _record_files() -> tuple[list[Path], list[Path]]:
    """(turn record files, inbox files), sorted."""
    def _glob(directory: Path) -> list[Path]:
        try:
            return sorted(p for p in directory.glob("*.jsonl")
                          if not p.name.endswith(".requests.jsonl"))
        except OSError:
            return []
    return _glob(_turn_metrics_dir()), _glob(_inbox_dir())


def _first_run_positions(state: dict) -> None:
    """A new cursor starts at 'now': history is not replayed as events."""
    pos: dict = {}
    turns, inboxes = _record_files()
    for p in turns + inboxes:
        fd = _open_record(p)
        if fd is None:
            continue
        try:
            st = os.fstat(fd)
            pos[str(p)] = _mark(fd, st, st.st_size)
        except OSError:
            pass
        finally:
            os.close(fd)
    state["pos"] = pos
    try:
        from . import audit_log
        # The end of the log, without parsing it: steps_since(0) read the
        # whole log into memory to learn where it ends.
        state["audit_off"] = audit_log._default_log_path().stat().st_size
    except Exception:
        state["audit_off"] = 0
    state["sessions"] = sorted(str(r.get("key") or "")
                               for r in _open_sessions())
    state["approvals"] = sorted(a["id"] for a in pending_approvals())
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
            events.append({"kind": "session_opened",
                           "session": _bounded(k, 80)})
    for k in sorted(before - set(keys)):
        events.append({"kind": "session_exited", "session": _bounded(k, 80)})
    state["sessions"] = keys

    approvals = {a["id"]: a for a in pending_approvals()}
    seen = set(state.get("approvals") or [])
    for rid, a in approvals.items():
        if rid not in seen:
            events.append({"kind": "approval_open", **a})
    for rid in sorted(seen - set(approvals)):
        events.append({"kind": "approval_closed", "id": rid})
    state["approvals"] = sorted(approvals)

    files, inboxes = _record_files()
    # Positions of files that are gone are dropped: a cursor kept for
    # months otherwise names every inbox and turn record that ever was.
    live = {str(p) for p in files + inboxes}
    for key in [k for k in pos if k not in live]:
        pos.pop(key, None)
    for p in files:
        sid = p.name[:-len(".jsonl")]
        for line in _read_new_lines(p, pos):
            try:
                rec = json.loads(line)
            except ValueError:
                continue
            if not isinstance(rec, dict):
                continue
            end = ("stopped" if rec.get("stopped") else
                   "error" if rec.get("error") else "answer")
            events.append({
                "kind": "turn_ended", "session": by_id.get(sid, ""),
                "session_id": _bounded(sid, 128), "end": end,
                "tool_calls": _int_or_none(rec.get("tool_calls")),
                "error": _bounded(rec.get("error") or "", 300),
                "turn_ts": _float_or_none(rec.get("ts")),
            })

    # session_messages shows a session only the mail it sent (status, ls).
    # From inside an agent session, a message event says who wrote to
    # whom and keeps the text back.
    withhold = _inside_agent()
    for p in inboxes:
        for line in _read_new_lines(p, pos):
            try:
                msg = json.loads(line)
            except ValueError:
                continue
            if not isinstance(msg, dict):
                continue
            events.append({
                "kind": "message", "id": _bounded(msg.get("id") or "", 64),
                "from": _bounded(msg.get("from") or "", 80),
                "to": _bounded(msg.get("to") or p.stem, 80),
                "text": "" if withhold else _bounded(msg.get("text") or ""),
                "sent_at": _float_or_none(msg.get("sent_at")),
            })
            if withhold:
                events[-1]["text_withheld"] = True

    try:
        from . import audit_log
        rows, state["audit_off"] = audit_log.steps_since(
            int(state.get("audit_off") or 0))
    except Exception:
        rows = []
    for row in rows:
        sid = str(row.get("session_id") or "")
        decision = str(row.get("decision") or "")
        command = str(row.get("command") or "")
        if decision and decision != "ok":
            events.append({
                "kind": "denied", "session": by_id.get(sid, ""),
                "session_id": _bounded(sid, 128),
                "tool": _bounded(row.get("tool") or "", 80),
                "decision": _bounded(decision, 80),
                "detail": _bounded(command, 300),
            })
        elif (re.match(r"^\s*git\s+commit\b", command)
              and row.get("exit_code") == 0):
            cwd = str(row.get("cwd") or row.get("workspace") or "")
            events.append({
                "kind": "commit", "session": by_id.get(sid, ""),
                "session_id": _bounded(sid, 128),
                "branch": _bounded(_git(cwd, "rev-parse", "--abbrev-ref",
                                        "HEAD"), 200),
                "head": _bounded(_git(cwd, "log", "-1", "--format=%h %s"),
                                 200),
            })

    idle_before = set(state.get("idle") or [])
    idle_now = set()
    waiting = {a["session"] for a in approvals.values()} | {
        a["session_id"] for a in approvals.values()}
    for r in open_rows:
        sid = str(r.get("session_id") or "")
        key = str(r.get("key") or "")
        activity = _activity(sid)
        if (activity and now - activity >= idle_s
                and key not in waiting and sid not in waiting):
            idle_now.add(key)
            if key not in idle_before:
                events.append({"kind": "session_idle", "session": key,
                               "session_id": sid,
                               "idle_s": int(now - activity)})
    state["idle"] = sorted(idle_now)

    for ev in events:
        state["seq"] = int(state.get("seq") or 0) + 1
        ev.update({"schema": SCHEMA, "seq": state["seq"], "ts": now})
    return events


def parse_spec(spec: str) -> list[tuple[str, str, str]]:
    """``kind=approval_open,session=s11,idle_s>=1800`` as (key, op, value)
    terms. Raises ValueError for a filter that cannot match as written:
    a term without an operator, an empty key, a non-number after ``>=``
    or ``<=``, more than _SPEC_TERMS terms, more than _SPEC_MAX
    characters. No evaluation: three operators, string or float compare.
    """
    text = str(spec or "")
    if len(text) > _SPEC_MAX:
        raise ValueError(f"the filter is longer than {_SPEC_MAX} characters")
    raw = text.split(",")
    terms = [t.strip() for t in raw]
    if not text.strip() or any(not t for t in terms):
        raise ValueError("the filter has an empty term")
    if len(terms) > _SPEC_TERMS:
        raise ValueError(f"the filter has more than {_SPEC_TERMS} terms")
    out = []
    for term in terms:
        m = _TERM_RE.match(term)
        if not m:
            raise ValueError(f"{term!r} is not key=value, key>=n or key<=n")
        key, op, want = m.groups()
        if op == "=" and want.startswith(("=", "<", ">")):
            raise ValueError(f"{term!r}: the operators are =, >= and <=")
        if op != "=":
            try:
                float(want)
            except ValueError:
                raise ValueError(f"{term!r}: {op} needs a number") from None
        out.append((key, op, want))
    return out


def matches(event: dict, spec: Any) -> bool:
    """Whether *event* holds every term of *spec* (a string or the terms
    from :func:`parse_spec`). An unparsable string matches nothing."""
    if isinstance(spec, str) or spec is None:
        try:
            terms = parse_spec(spec or "")
        except ValueError:
            return False
    else:
        terms = spec
    for key, op, want in terms:
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


def next_events(state: dict, *, idle_s: int = DEFAULT_IDLE_S) -> list[dict]:
    """The events a reader gets next: those a previous reader of this
    cursor left unread (a waiter that stopped in the middle of a batch),
    then the new ones. Updates *state*; the caller saves it."""
    queued = state.pop("queue", None)
    batch = [e for e in queued if isinstance(e, dict)] if isinstance(
        queued, list) else []
    return batch + poll(state, idle_s=idle_s)


def follow(cursor: str = "default", *, idle_s: int = DEFAULT_IDLE_S,
           interval: float = 5.0, seconds: float = 0.0) -> Iterator[dict]:
    """Events as they happen.

    The cursor is saved before each event is handed out, holding the
    events of the batch that come after it. A reader that stops after
    any event -- ``wait`` on a match, a consumer that ends -- leaves the
    rest for the next reader of the cursor instead of consuming them.
    """
    _cursor_path(cursor)                     # a bad name fails before work
    interval = max(0.05, float(interval or 0.0))
    state = load_cursor(cursor)
    deadline = time.monotonic() + (seconds if seconds > 0 else 1e12)
    while True:
        batch = next_events(state, idle_s=idle_s)
        if not batch:
            save_cursor(state, cursor)
        for i, ev in enumerate(batch):
            state["queue"] = batch[i + 1:i + 1 + _QUEUE_MAX]
            save_cursor(state, cursor)
            yield ev
        state.pop("queue", None)
        if time.monotonic() >= deadline:
            if batch:
                save_cursor(state, cursor)
            return
        time.sleep(interval)


def wait(spec: str, *, cursor: str = "default", timeout: float = 0.0,
         idle_s: int = DEFAULT_IDLE_S, interval: float = 5.0
         ) -> Optional[dict]:
    """The first event that matches *spec*, or None when *timeout* passes.

    Raises ValueError for a filter :func:`parse_spec` refuses, before
    reading anything: a typo waited for ever (or until the timeout) and
    looked like a quiet wave. Events before the match that do not match
    are consumed; those after it stay for the next reader.
    """
    terms = parse_spec(spec)
    gen = follow(cursor, idle_s=idle_s, interval=interval, seconds=timeout)
    try:
        for ev in gen:
            if matches(ev, terms):
                return ev
    finally:
        gen.close()
    return None


# --------------------------------------------------------------------------
# actions that cannot widen anything
# --------------------------------------------------------------------------

class OpsRefused(RuntimeError):
    """An ops action that may not run from here."""


def _inside_agent() -> bool:
    try:
        from . import approval_guard
        return approval_guard.inside_agent_session()
    except ImportError:
        return False


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

    The receiving session reads it under the session-message header
    ("[Message from the session ... not from the user"), which every
    reader of user intent skips: it does not count as a typed message,
    lifts nothing and is never a slash command. ``session_messages.send``
    scrubs credentials and terminal control and caps the length.

    Refused: an address no session has (a typo queued a message nobody
    would take, which is what the session_message tool refuses too), and
    the operator's own mailbox in any spelling.
    """
    _refuse_inside_agent()
    name = str(name or "").strip()
    if not name or name.lower() == OPERATOR_KEY:
        raise OpsRefused("send needs the name of a session")
    from . import session_messages
    if not session_messages.deliverable(name):
        raise OpsRefused(f"no session is known as {name!r} "
                         "(delfin-agent ops sessions lists them)")
    return session_messages.send(name, str(text or ""),
                                 from_key=OPERATOR_KEY,
                                 from_title=OPERATOR_TITLE)


def _workspace_of(approval: dict) -> str:
    for rec in _open_sessions():
        if approval["session_id"] and str(
                rec.get("session_id") or "") == approval["session_id"]:
            return str(rec.get("workspace") or "")
        if approval["session"] and str(rec.get("key") or "") == \
                approval["session"]:
            return str(rec.get("workspace") or "")
    return ""


def escalate(approval_id: str, note: str = "") -> str:
    """Put a waiting approval in front of the person. It stays open.

    The notice carries the workspace of the waiting session and no
    session id. With the session id, the engine reads the entry as a
    confirmation window that expired while the person was away and tells
    the model to ask again -- while the request is still open. Without a
    workspace, a reply the person types to the notice is drained by
    whichever session looks first. The detail names the command that
    answers; replying to the notice answers nothing.
    """
    _refuse_inside_agent()
    rows = {a["id"]: a for a in pending_approvals()}
    a = rows.get(str(approval_id))
    if a is None:
        raise OpsRefused(f"nothing waiting with id {approval_id!r}")
    rid = a["id"]
    from . import attention
    how = (f"Answer it with delfin-agent approvals show {rid}, then "
           f"delfin-agent approvals approve {rid} or deny {rid} (or at "
           "the session). A reply to this notice does not answer it.")
    detail = (f"Operator note: {_bounded(note, 400)}\n" if note else
              "The operator passed this approval on to you.\n") + how
    return attention.emit_attention(
        "confirm_pending",
        title=_bounded(f"{a['session'] or a['session_id'][:8]} waits on "
                       f"{rid}: {a['tool']} {a['summary']}", 200),
        detail=_bounded(detail, 800),
        workspace=_workspace_of(a) or None,
    )


def pause(name: str, reason: str) -> dict:
    """Pause session *name* (it stops at its next safe point)."""
    _refuse_inside_agent()
    name = str(name or "").strip()
    if not name:
        raise OpsRefused("pause needs a session name")
    from . import session_pause
    return session_pause.pause(name, _bounded(reason or "operator", 200))
