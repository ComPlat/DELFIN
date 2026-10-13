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
import secrets
import threading
from dataclasses import dataclass

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
    """Drop a finished request's code. Never raises."""
    with _LOCK:
        _ISSUED.pop(str(request_id), None)
        _EXEMPT.discard(str(request_id))


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
           "spent"]
