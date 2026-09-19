"""Approvals for a session with nobody at its terminal.

The terminal broker is bound only when ``stdin`` is a tty. A session
driven headlessly -- ``delfin-agent chat -p ...``, a supervisor running
several at once -- therefore has no callback at all, and the one thing
that always needs a human answer, the Self-Modification Guard, comes back
as "no confirm_callback is configured, so it is refused". That is honest,
and it means such a session cannot touch the agent's own safety layer:
the most valuable work there is exactly what it cannot do.

This broker gives that session somebody to ask. The request is written to
a directory; whoever is watching writes the answer beside it; the waiting
thread reads it and carries on. Same contract as the other two brokers --
a bound ``callback`` the gate can reach ``last_timed_out`` through, and
that flag per THREAD, because requests overlap.

It is a security boundary, so it is built to refuse rather than to guess:

* Nothing is approved that was not asked. The answer must name the
  request's id, and an answer whose file is older than the request is
  ignored -- otherwise a pre-written file would approve the next thing
  that happens to come along.
* The answer must come from this user. A file owned by somebody else, or
  writable by group or world, is not read. Nor is a symlink: the target
  of one is not what the directory listing showed.
* Anything unreadable, malformed, or unrecognised is NOT an approval. The
  only value that approves is the exact string.
* Silence is not a refusal. On timeout the flag goes up and the gate
  says the user is away, so that a path nobody looked at does not close
  for the rest of the session.

The directory is the whole interface. A supervisor is any process that
can read it -- the CLI in ``delfin-agent approvals``, or a person with an
editor.
"""

from __future__ import annotations

import json
import os
import stat as _stat
import threading
import time
import uuid
from pathlib import Path
from typing import Optional

#: The one value that approves. Anything else -- including True, 1, "yes"
#: and a missing field -- does not.
APPROVE = "approve"
DENY = "deny"

#: How long a request waits before it is called absent rather than answered.
#: The same 300 s the dashboard has always used.
DEFAULT_TIMEOUT_S = 300.0

_POLL_S = 0.25


def requests_dir() -> Path:
    """Where requests and answers live.

    Resolved per call rather than at import, so the state-path table can
    redirect it: a test that approved something must not leave an answer
    in the user's real directory, where a later session would read it.
    """
    return Path.home() / ".delfin" / "approvals"


class FileConfirmBroker:
    """Ask over the filesystem; block until answered, refused, or absent."""

    def __init__(self, *, session_id: str = "", timeout_s: float = 0.0,
                 poll_s: float = _POLL_S) -> None:
        self.session_id = str(session_id or "") or uuid.uuid4().hex[:8]
        self.timeout_s = float(timeout_s or DEFAULT_TIMEOUT_S)
        self.poll_s = float(poll_s or _POLL_S)
        self._timed_out = threading.local()

    # -- the contract the gate reads -------------------------------------
    @property
    def last_timed_out(self) -> bool:
        """Whether THIS thread's last request expired rather than being
        refused. Per thread because requests arrive on whichever thread
        is running tools and they overlap; one flag for all of them
        recorded an expiry as a refusal and closed a path for good."""
        return bool(getattr(self._timed_out, "value", False))

    @last_timed_out.setter
    def last_timed_out(self, value: bool) -> None:
        self._timed_out.value = bool(value)

    def callback(self, tool_name: str, args: dict, preview: str) -> bool:
        """Bound on purpose: the gate reads last_timed_out off __self__."""
        self.last_timed_out = False
        try:
            request_id, path = self._write_request(tool_name, args, preview)
        except OSError:
            # No way to ask is not permission to act.
            return False
        answer = self._wait_for(request_id, path)
        if answer is None:
            self.last_timed_out = True
            return False
        return answer

    # -- writing the question --------------------------------------------
    def _write_request(self, tool_name, args, preview) -> tuple[str, Path]:
        room = _own_dir(requests_dir())
        request_id = f"{int(time.time())}-{uuid.uuid4().hex[:8]}"
        record = {
            "id": request_id,
            "session_id": self.session_id,
            "tool": str(tool_name or ""),
            "path": str((args or {}).get("path") or ""),
            "preview": str(preview or "")[:20000],
            "pid": os.getpid(),
            "host": _hostname(),
            "asked_at": time.time(),
            "timeout_s": self.timeout_s,
        }
        path = room / f"{request_id}.request.json"
        # Written whole and then moved into place, so a watcher never
        # reads half a question and answers it.
        tmp = room / f".{request_id}.partial"
        tmp.write_text(json.dumps(record, indent=2), encoding="utf-8")
        os.chmod(tmp, 0o600)
        tmp.replace(path)
        return request_id, path

    # -- waiting for the answer ------------------------------------------
    def _wait_for(self, request_id: str, request_path: Path) -> Optional[bool]:
        """True/False once answered, None when the time runs out."""
        answer_path = request_path.with_name(f"{request_id}.answer.json")
        deadline = time.monotonic() + self.timeout_s
        asked_at = _mtime(request_path)
        while True:
            decision = _read_answer(answer_path, request_id, asked_at)
            if decision is not None:
                _retire(request_path, answer_path, decision)
                return decision
            if time.monotonic() >= deadline:
                _retire(request_path, answer_path, None)
                return None
            time.sleep(self.poll_s)


# -- answering, from the other side -----------------------------------------

def pending(room: Optional[Path] = None) -> list:
    """Every question waiting for an answer, oldest first."""
    room = Path(room) if room else requests_dir()
    out = []
    try:
        entries = sorted(room.glob("*.request.json"))
    except OSError:
        return out
    for path in entries:
        record = _read_json(path)
        if not record:
            continue
        if path.with_name(f"{record.get('id')}.answer.json").exists():
            continue
        record["_file"] = str(path)
        out.append(record)
    return out


def answer(request_id: str, decision: bool, *, room: Optional[Path] = None,
           by: str = "") -> bool:
    """Write the answer to *request_id*. True when it was written."""
    room = Path(room) if room else requests_dir()
    request_path = room / f"{request_id}.request.json"
    if not request_path.is_file():
        return False
    record = {"id": str(request_id),
              "decision": APPROVE if decision else DENY,
              "by": str(by or ""), "at": time.time()}
    path = room / f"{request_id}.answer.json"
    tmp = room / f".{request_id}.answer.partial"
    try:
        _own_dir(room)
        tmp.write_text(json.dumps(record), encoding="utf-8")
        os.chmod(tmp, 0o600)
        tmp.replace(path)
    except OSError:
        return False
    return True


# -- the careful parts ------------------------------------------------------

def _read_answer(path: Path, request_id: str,
                 asked_at: float) -> Optional[bool]:
    """The decision in *path*, or None when there is no answer to read.

    Every reason to doubt the file returns None: doubt is not approval,
    and it is not refusal either -- the question simply stays open.
    """
    try:
        info = path.lstat()
    except OSError:
        return None
    if _stat.S_ISLNK(info.st_mode):
        return None                      # not what the listing showed
    if info.st_uid != os.getuid():
        return None                      # somebody else wrote it
    if info.st_mode & (_stat.S_IWGRP | _stat.S_IWOTH):
        return None                      # anyone could have written it
    if asked_at and info.st_mtime < asked_at - 1.0:
        return None                      # written before the question
    record = _read_json(path)
    if not record or str(record.get("id") or "") != str(request_id):
        return None                      # an answer to something else
    decision = record.get("decision")
    if decision == APPROVE:
        return True
    if decision == DENY:
        return False
    return None                          # unrecognised is not an answer


def _own_dir(room: Path) -> Path:
    room.mkdir(parents=True, exist_ok=True)
    try:
        os.chmod(room, 0o700)
    except OSError:
        pass
    return room


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _mtime(path: Path) -> float:
    try:
        return path.stat().st_mtime
    except OSError:
        return 0.0


def _retire(request_path: Path, answer_path: Path,
            decision: Optional[bool]) -> None:
    """Move the finished exchange out of the pending listing.

    Kept, not deleted: what was approved and by whom is the record of how
    the agent's own safety layer came to be changed.
    """
    room = request_path.parent
    done = room / "answered"
    try:
        done.mkdir(parents=True, exist_ok=True)
        os.chmod(done, 0o700)
        stamp = "timeout" if decision is None else (
            APPROVE if decision else DENY)
        if request_path.exists():
            request_path.replace(done / f"{request_path.stem}.{stamp}.json")
        if answer_path.exists():
            answer_path.replace(done / answer_path.name)
    except OSError:
        pass


def _hostname() -> str:
    import socket
    try:
        return socket.gethostname()
    except OSError:
        return ""


__all__ = ["APPROVE", "DENY", "DEFAULT_TIMEOUT_S", "FileConfirmBroker",
           "answer", "pending", "requests_dir"]
