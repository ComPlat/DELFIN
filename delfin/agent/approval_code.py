"""The one-time code that makes an out-of-process approval answer the person's.

Input: a request id, issued by the process that waits on it. Output: a
six-character code held only in that process's memory and shown only on
its human channel (the terminal frame, the dashboard row, the headless
session's controlling terminal), and the check an answer file must pass
to be read as an answer.

Threat model. An answer file is accepted when it is this account's, not
group/world-writable, not a symlink, not older than the question and
naming its id (``file_confirm._checked_record``). Every process of the
account can write such a file, including an agent's shell that got past
the command-text gates of ``approval_guard``, and a process the agent
started outside DELFIN. The code is what such a process does not have:

* it is generated with ``secrets`` in the waiting process and never
  written to the request file, an answer file, a log, an audit or
  security event, a tool result, a transcript or stream-json output;
* an answer carries ``proof = HMAC-SHA256(code, id NUL decision NUL
  detail)``, not the code, where *detail* is the refusal reason or the
  chosen labels -- so a signed refusal cannot be rewritten into an
  approval, and an answer file read by the agent after the fact gives
  away nothing reusable (the request is retired by the first valid one);
* an answer file without a valid proof is moved aside, recorded as an
  ``approval_code_refused`` security event and counted; the third
  expires the request (fail closed, the gate reports "user away").

Guessing: 32**6 = 1.07e9 codes, three tries per request, so a blind
forger succeeds with p = 2.8e-9 per request and every try is recorded.

A run with no controlling terminal (a batch or cluster job) has no
screen to show the code on. Its code goes to ``codes_dir()/<id>.code``
(0700/0600), which every agent read route refuses and every cage masks;
``delfin-agent approvals`` run by the person (no ``DELFIN_AGENT_SESSION``)
reads it from there, so such a run stays answerable from another shell.
The file is deleted when the request is answered or expires. Without any
isolation a same-account program can read that file: for this case the
protection is the read gates, the same trust level as an unsigned answer
file had before the code existed.

``codes_dir()`` is the per-user runtime directory where one exists
(``$XDG_RUNTIME_DIR`` or ``/run/user/<uid>``, 0700, outside every
workspace), else ``~/.delfin/approval_codes`` (and always that inside a
batch job, whose node is not the person's). Readers look in both. Limit,
measured on the Landlock path: Landlock cannot hide a path that lies
inside a root it grants for writing without refusing every new file at
that root (it only grants, beneath a path), so where the code directory
lies inside a writable root -- a granted ``~/.delfin``, a workspace that
holds ``$HOME`` -- on a host with Landlock and no usable bubblewrap, the
cage does not hide it and a command that builds the path at run time can
read a waiting code. There the read gates are the protection (the path
deny globs, the bash read scan, ``approval_guard``), and the file is
0600 and lives only while its request waits. bubblewrap (a tmpfs mask
after the binds) and Seatbelt (a deny rule) hide it in every layout.

What the code does not cover (see ``docs/USER_MANUAL.md``): a process
that can read the waiting process's memory (``ptrace`` / ``/proc/<pid>/
mem`` of the same uid: Linux Yama ``ptrace_scope`` >= 1 limits that to
ancestors, so an agent child cannot read its DELFIN parent; at 0 any
same-uid process can); a process that sees the human channel itself (an
outside agent driving the DELFIN terminal through a pty sees the frame);
and the CLI backend's own hook requests (``cli_gate``), which run inside
the external agent's process tree and are answered by DELFIN's relay
without a code (``exempt``).

Rejected alternatives: the code in the request file with mode 0600 (the
agent runs as the same uid and reads it); a code typed only through a
tty check in ``approvals approve`` (a pty is a tty, and the check runs in
the forger's process); a long-lived per-session key shared with the CLI
(it would have to live in a file or the environment, both of which the
agent's shell reaches); the code itself in the answer (a refusal could be
turned into an approval in the poll window).

In-process callers (the dashboard and terminal brokers' own helpers, the
tests) may answer with ``code=None``: the trust boundary is the process,
and a caller inside it can already resolve the request directly.
"""

from __future__ import annotations

import hashlib
import hmac
import os
import re
import secrets
import stat as _stat
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

#: No 0/O, 1/I: the person reads it off one screen and types it on another.
ALPHABET = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"
LENGTH = 6
#: The answer files without a valid proof a request survives.
MAX_WRONG = 3

OK = "ok"
WRONG = "wrong"


@dataclass
class _Entry:
    code: str
    wrong: int = 0
    spent: bool = False


_LOCK = threading.Lock()
_ISSUED: "dict[str, _Entry]" = {}
_EXEMPT: "set[str]" = set()


def issue(request_id: str) -> str:
    """A fresh code for *request_id*, kept in this process. Output: the code."""
    code = "".join(secrets.choice(ALPHABET) for _ in range(LENGTH))
    with _LOCK:
        _ISSUED[str(request_id)] = _Entry(code)
    return code


def exempt(request_id: str) -> None:
    """Mark *request_id* as answered without a code (``cli_gate`` hook)."""
    with _LOCK:
        _EXEMPT.add(str(request_id))


def code_for(request_id: str) -> str:
    """The code this process issued for *request_id*, or ""."""
    with _LOCK:
        entry = _ISSUED.get(str(request_id))
    return entry.code if entry else ""


def forget(request_id: str) -> None:
    """Drop a finished request's code, and its code file. Never raises."""
    with _LOCK:
        _ISSUED.pop(str(request_id), None)
        _EXEMPT.discard(str(request_id))
    remove_code_file(request_id)


# -- the code file: the human channel of a run with no terminal -------------

_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,80}$")


#: Environment variables that mark a batch-scheduler job. Its node's
#: runtime directory is not where the person runs ``approvals``.
_BATCH_JOB_VARS = ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "JOB_ID",
                   "FLUX_JOB_ID")
#: The per-user runtime directory when ``XDG_RUNTIME_DIR`` is not set.
_RUN_USER_FMT = "/run/user/{uid}"


def _private_dir(path: str) -> Optional[Path]:
    """*path* when it is a real directory (no symlink) of this uid with no
    group/world bits, else None."""
    if not path or not os.path.isabs(path):
        return None
    try:
        info = os.lstat(path)
    except OSError:
        return None
    if (not _stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid()
            or info.st_mode & 0o077):
        return None
    return Path(path)


def runtime_root() -> Optional[Path]:
    """The per-user runtime directory (a tmpfs, 0700, never a workspace),
    or None: ``$XDG_RUNTIME_DIR``, else ``/run/user/<uid>``, each only
    when :func:`_private_dir` accepts it. None inside a batch job
    (:data:`_BATCH_JOB_VARS`): there the directory is the compute node's,
    and the person answers from another node."""
    if any(os.environ.get(v) for v in _BATCH_JOB_VARS):
        return None
    xdg = os.environ.get("XDG_RUNTIME_DIR", "").strip()
    for candidate in (xdg, _RUN_USER_FMT.format(uid=os.getuid())):
        found = _private_dir(candidate)
        if found is not None:
            return found
    return None


def home_codes_dir() -> Path:
    """``~/.delfin/approval_codes``: the location without a runtime
    directory, and the one every reader also looks in."""
    return Path.home() / ".delfin" / "approval_codes"


def _default_codes_dir() -> Path:
    root = runtime_root()
    if root is not None:
        return root / "delfin" / "approval_codes"
    return home_codes_dir()


def codes_dir() -> Path:
    """Where this process writes code files. Resolved per call (tests move
    HOME); redirected as a user-state sink (``state_paths``).

    ``<runtime_root()>/delfin/approval_codes`` where a private runtime
    directory exists, else ``~/.delfin/approval_codes``. The runtime
    directory first because it lies under no plausible write root: a
    Landlock cage cannot hide a path inside a root it grants for writing
    without taking that root's new-file rights away
    (``api_client._landlock_hide_paths``), so a code directory inside a
    granted ``~/.delfin`` or a workspace holding ``$HOME`` is masked by
    bubblewrap and Seatbelt but not by Landlock. Rejected: the runtime
    directory everywhere (a batch job's node is not the person's), and a
    fixed ``~/.delfin`` location (granted roots reach it).

    A deny-read location for every agent route: the path deny globs
    (read_file, grep/glob, the bash read scan, cli_gate), the command-text
    gate (``approval_guard``), and the cages (``sandbox._HOME_SECRET_DIRS``
    and ``sandbox.runtime_masked_dirs``: a bubblewrap tmpfs, a
    Landlock/Seatbelt hide)."""
    return _default_codes_dir()


def _read_dirs() -> list[Path]:
    """Where a code file may be: :func:`codes_dir`, and ``~/.delfin``'s
    when the writer may have had no runtime directory (another session,
    a batch job). Only for the unredirected location: a redirected one
    (tests, a benchmark attempt) is the only place."""
    first = codes_dir()
    out = [first]
    root = runtime_root()
    if root is not None and first == root / "delfin" / "approval_codes":
        fallback = home_codes_dir()
        if fallback != first:
            out.append(fallback)
    return out


def ensure_codes_dir() -> Optional[Path]:
    """Create the directory (0700) so a cage built now masks it, rather
    than skipping a path that appears later. None when it cannot be made."""
    d = codes_dir()
    try:
        d.mkdir(parents=True, exist_ok=True)
        os.chmod(d, 0o700)
        return d
    except OSError:
        return None


def _code_path(request_id: str) -> Optional[Path]:
    rid = str(request_id or "")
    return codes_dir() / f"{rid}.code" if _ID_RE.match(rid) else None


def write_code_file(request_id: str, code: str) -> Optional[Path]:
    """Write *code* for *request_id*, mode 0600, whole. Output: the path,
    or None where it cannot be written."""
    path = _code_path(request_id)
    if path is None or ensure_codes_dir() is None:
        return None
    tmp = path.with_name(f".{path.name}.partial")
    try:
        fd = os.open(str(tmp), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, "w", encoding="ascii") as fh:
            fh.write(normalize(code) + "\n")
        os.chmod(tmp, 0o600)
        tmp.replace(path)
        return path
    except OSError:
        return None


def read_code_file(request_id: str) -> str:
    """The code in *request_id*'s file, normalised, or "". Only a regular
    file of this uid, not group/world-accessible, is read."""
    rid = str(request_id or "")
    if not _ID_RE.match(rid):
        return ""
    for d in _read_dirs():
        path = d / f"{rid}.code"
        try:
            info = path.lstat()
            if (not _stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid()
                    or info.st_mode & 0o077):
                continue
            code = normalize(path.read_text(encoding="ascii", errors="ignore"))
        except OSError:
            continue
        if code:
            return code
    return ""


def remove_code_file(request_id: str) -> None:
    """Delete *request_id*'s code file. Never raises."""
    path = _code_path(request_id)
    if path is None:
        return
    try:
        path.unlink()
    except OSError:
        pass


def normalize(text) -> str:
    """What the person typed, compared case- and separator-free."""
    return "".join(c for c in str(text or "").upper() if c.isalnum())


def display(code: str) -> str:
    """``ABCDEF`` as ``ABC-DEF``: easier to read off and to type."""
    c = normalize(code)
    return f"{c[:3]}-{c[3:]}" if len(c) == LENGTH else c


def detail_of(record: dict) -> str:
    """The part of an answer besides its decision that the proof binds."""
    decision = str(record.get("decision") or "")
    if decision == "choose":
        answers = record.get("answers")
        if not isinstance(answers, list):
            return ""
        return "\n".join(str(a) for a in answers)
    if decision == "deny":
        return str(record.get("reason") or "")
    return ""


def proof(code: str, request_id: str, decision: str, detail: str = "") -> str:
    """HMAC-SHA256 keyed with the code over id, decision and detail (hex)."""
    msg = "\0".join((str(request_id), str(decision), str(detail or "")))
    return hmac.new(normalize(code).encode("ascii", "ignore"),
                    msg.encode("utf-8"), hashlib.sha256).hexdigest()


def sign(record: dict, code: str) -> dict:
    """*record* with its ``proof`` field set for *code*. Returns the record."""
    record["proof"] = proof(code, str(record.get("id") or ""),
                            str(record.get("decision") or ""),
                            detail_of(record))
    return record


def required(request_id: str) -> bool:
    """False only for a request this process exempted."""
    with _LOCK:
        return str(request_id) not in _EXEMPT


def check(request_id: str, record: dict) -> str:
    """``OK`` when *record* is signed with this process's code for the
    request (or the request is exempt), else ``WRONG``.

    A request this process neither issued nor exempted is ``WRONG``: the
    reader of a request is the process that asked it, so an unknown id
    is never an answer.
    """
    rid = str(request_id)
    with _LOCK:
        if rid in _EXEMPT:
            return OK
        entry = _ISSUED.get(rid)
    if entry is None or entry.spent:
        return WRONG
    given = str(record.get("proof") or "")
    want = proof(entry.code, rid, str(record.get("decision") or ""),
                 detail_of(record))
    return OK if hmac.compare_digest(given, want) else WRONG


def note_wrong(request_id: str) -> int:
    """Count one answer without a valid proof. Output: the count so far;
    at ``MAX_WRONG`` the request is spent."""
    with _LOCK:
        entry = _ISSUED.get(str(request_id))
        if entry is None:
            return MAX_WRONG
        entry.wrong += 1
        if entry.wrong >= MAX_WRONG:
            entry.spent = True
        return entry.wrong


def spent(request_id: str) -> bool:
    """True once a request has seen ``MAX_WRONG`` wrong answers."""
    with _LOCK:
        entry = _ISSUED.get(str(request_id))
    return bool(entry and entry.spent)


__all__ = ["ALPHABET", "LENGTH", "MAX_WRONG", "OK", "WRONG", "issue",
           "exempt", "code_for", "forget", "normalize", "display",
           "detail_of", "proof", "sign", "required", "check", "note_wrong",
           "spent", "runtime_root", "home_codes_dir", "codes_dir",
           "ensure_codes_dir", "write_code_file",
           "read_code_file", "remove_code_file"]
