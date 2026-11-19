"""Append-only audit log for every code-modifying KIT-agent action.

Every successful or denied call to write_file / edit_file / multi_edit /
notebook_edit / bash / bash_background / remember_permission(_bundle)
appends one JSON line to ``~/.delfin/audit.log``. Records carry:

* ``ts`` — ISO-8601 UTC timestamp
* ``session_id`` — engine session UUID, if any
* ``tool`` — tool name (lowercase, no MCP prefix)
* ``path`` / ``command`` — the target (file path / shell text)
* ``mode`` — KIT mode active at the time
* ``decision`` — ``auto`` (gate passed without a confirm) /
  ``confirmed`` (user clicked Erlauben) / ``denied`` (gate blocked or
  user clicked Ablehnen) / ``error`` (tool raised)
* ``persistence`` — for remember_* calls, what got written and where
* ``diff_lines_added`` / ``diff_lines_removed`` — for write/edit tools
* ``pid`` — the calling process pid (helps trace bash_background jobs)

Lines are written atomically (single ``write`` call with a final
``\\n``) so concurrent writes from multiple worker threads don't
interleave inside one record. Reads should treat each line
independently — corrupt lines (e.g. from a kill -9 mid-write) are
rare but possible; tooling should ``json.loads`` line-by-line and
skip what doesn't parse.

The log is a HASH CHAIN. Every record carries ``prev``: the SHA-256 of
the previous line's bytes as written (``GENESIS`` for the first line of
a new log; after a rotation, the hash of the rotated log's last line).
A head file beside the log (``audit.log.head``, 0600) holds the hash and
the byte size of the last line written, so a tail cut off the log shows
too. ``verify`` (``delfin-agent audit verify``) reports an edited,
removed or inserted line, a cut tail and a line that does not parse;
lines written before the chain existed are counted as legacy. Appends
are serialised across processes by ``flock`` on ``audit.log.lock``.
Where the head file did not match the log at the time of a write (the
log was changed outside DELFIN), that record carries ``chain_break`` and
a blocked ``audit_chain_broken`` event is recorded.

The chain is KEYED. A record written by a process that holds the account
key carries ``mac`` as its last field: HMAC-SHA256, under
``trust_store.derived_key("audit_log")``, of the line's bytes without
that field (the line up to its final ``}``). And the trust store keeps an
anchor per log (``audit-anchor-<hash of the path>.json``: the last line's
hash and the log's size, signed), written under the same lock as the
line. The anchor pins the active log's bytes; through ``prev`` it pins
every line before the last, and through the first record's ``prev`` the
rotated log it continues. So a same-account process without the key
that rewrites the whole chain and the head file consistently is found
(``mac`` problems, or the log not ending where the anchor says), and so
is one that cuts the tail and puts an earlier head file back.

Migration: records of a log written before the key existed carry
``prev`` and no ``mac``; ``verify`` counts them as unkeyed. From the
first keyed record on, a record without ``mac`` is a problem (written by
a process without the key: one started inside an agent session, or with
the trust store refused). Where verify runs without the key (inside an
agent session) the hash chain is checked and the report says the
signatures were not. A process with the key on the same account that
reads the trust store directly is outside this boundary
(docs/SECURITY_MODEL.md); the administrator's copy
(``managed_policy.append_audit_line``) receives the same lines.

The log is created 0600 in a 0700 directory and opened with
``O_APPEND | O_NOFOLLOW``; a link, a non-regular file, a file with more
than one link or one owned by another account in its place is not
written. A write that fails records one blocked ``audit_write_failed``
event per process (it reaches the attention inbox) instead of vanishing.

The log is ROTATED weekly: when the first line of a calendar week is
appended and the existing log is older than the week boundary, the
old log is renamed to ``audit-YYYY-Www.log`` and a fresh log starts.
This keeps the active file small enough to ``tail`` without pain
while preserving full history under predictable filenames. The changes
report falls back to the immediately preceding rotated file when the
active one has no answer, and names the file it consulted — on the
first day of a new week the active log is empty, and reporting that as
"no recorded changes" hid a whole week of work.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import re
import stat
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

try:
    import fcntl as _fcntl
except ImportError:          # pragma: no cover - non-POSIX
    _fcntl = None


_LOG_LOCK = threading.Lock()

#: ``prev`` of the first line of a log that continues nothing.
GENESIS = "0" * 64
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)
# Set after the first failed write of this process: one event, not one
# per record (a full disk would otherwise flood the panel).
_WRITE_FAILURE_NOTED = False


def _default_log_path() -> Path:
    return Path.home() / ".delfin" / "audit.log"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(ts: datetime) -> str:
    return ts.strftime("%Y-%m-%dT%H:%M:%SZ")


def _isocalendar_key(ts: datetime) -> str:
    """Return 'YYYY-Www' (ISO week) for rotation labelling."""
    iso = ts.isocalendar()
    return f"{iso.year}-W{iso.week:02d}"


def head_path(path: Path) -> Path:
    """The head file of ``path``: hash and size of its last line."""
    path = Path(path)
    return path.with_name(path.name + ".head")


def _line_hash(line: bytes) -> str:
    """SHA-256 hex of one line's bytes, without its newline."""
    return hashlib.sha256(line).hexdigest()


def _open_private(path: Path, flags: int) -> int:
    """``os.open`` that never follows a link at the last component, and
    accepts only a regular file of this account with one link."""
    fd = os.open(str(path), flags | _NOFOLLOW, 0o600)
    try:
        st = os.fstat(fd)
        if (not stat.S_ISREG(st.st_mode) or st.st_nlink > 1
                or (hasattr(os, "getuid") and st.st_uid != os.getuid())):
            raise OSError(f"{path} is not a private regular file")
        if stat.S_IMODE(st.st_mode) != 0o600:
            os.fchmod(fd, 0o600)
    except BaseException:
        os.close(fd)
        raise
    return fd


def _read_private(path: Path) -> bytes:
    fd = _open_private(path, os.O_RDONLY)
    try:
        chunks = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk:
                return b"".join(chunks)
            chunks.append(chunk)
    finally:
        os.close(fd)


def _last_line(fd: int, size: int) -> Optional[bytes]:
    """The last line of the open file (without its newline), or None for
    an empty file. Reads backwards from the end, doubling the window."""
    if size <= 0:
        return None
    window = 64 * 1024
    while True:
        start = max(0, size - window)
        data = os.pread(fd, size - start, start)
        body = data[:-1] if data.endswith(b"\n") else data
        cut = body.rfind(b"\n")
        if cut >= 0 or start == 0:
            return body[cut + 1:]
        window *= 2


def _read_head(path: Path) -> Optional[dict]:
    try:
        head = json.loads(_read_private(head_path(path)).decode("utf-8"))
    except FileNotFoundError:
        return None
    except Exception:
        return {"last": None, "size": None, "unreadable": True}
    return head if isinstance(head, dict) else {"last": None, "size": None}


def _write_head(path: Path, last: str, size: int) -> None:
    from .state_paths import write_text_atomic
    write_text_atomic(head_path(path),
                      json.dumps({"last": last, "size": int(size)}) + "\n")


#: The trust_store purpose of the record and anchor signatures.
_PURPOSE = "audit_log"
_MAC_TAIL = re.compile(rb', "mac": "([0-9a-f]{64})"\}$')


def _anchor_name(path: Path) -> str:
    try:
        real = os.path.realpath(str(path))
    except Exception:
        real = str(path)
    digest = hashlib.sha256(real.encode("utf-8", "surrogateescape"))
    return f"audit-anchor-{digest.hexdigest()[:24]}.json"


def _anchor_bytes(path: Path, last: str, size: int) -> bytes:
    return "\0".join((os.path.realpath(str(path)), str(last),
                       str(int(size)))).encode("utf-8", "surrogateescape")


def _write_anchor(path: Path, last: str, size: int) -> bool:
    """Record in the trust store where ``path`` ends now. Needs the key."""
    from . import trust_store
    tag = trust_store.mac(_PURPOSE, _anchor_bytes(path, last, size))
    if not tag:
        return False
    return trust_store.write_record(_anchor_name(path), {
        "path": os.path.realpath(str(path)), "last": last,
        "size": int(size), "mac": tag})


def _read_anchor(path: Path) -> Optional[dict]:
    try:
        from . import trust_store
        return trust_store.read_record(_anchor_name(path))
    except Exception:
        return None


def _sign_line(line: str, key: Optional[bytes]) -> str:
    """``line`` (one JSON object) with ``mac`` appended as its last field,
    or unchanged without a key."""
    if key is None or not line.endswith("}"):
        return line
    tag = hmac.new(key, line.encode("utf-8"), hashlib.sha256).hexdigest()
    return line[:-1] + f', "mac": "{tag}"}}'


def _audit_key() -> Optional[bytes]:
    try:
        from . import trust_store
        return trust_store.derived_key(_PURPOSE)
    except Exception:
        return None


def _rotate_if_needed(path: Path, now: datetime) -> Optional[str]:
    """Move ``path`` aside when the active log spans across an ISO week.

    Called inside the lock just before each append. If the existing
    log's first-line timestamp belongs to an earlier ISO week than
    ``now``, we rename it to ``audit-YYYY-Www.log`` (with its head file;
    ``-2``, ``-3`` … when that name is taken, so no earlier history is
    replaced) and start fresh. Returns the hash of the rotated log's
    last line, which the new log's first record carries as ``prev``;
    None when nothing was rotated. Cheap: the first line and the tail.
    """
    try:
        fd = _open_private(path, os.O_RDONLY)
    except OSError:
        return None
    try:
        size = os.fstat(fd).st_size
        if size == 0:
            return None
        first_line = os.pread(fd, min(size, 64 * 1024), 0).split(b"\n", 1)[0]
        first = json.loads(first_line.decode("utf-8", "replace"))
        first_ts = datetime.fromisoformat(first["ts"].rstrip("Z"))
        first_ts = first_ts.replace(tzinfo=timezone.utc)
        last = _last_line(fd, size)
    except Exception:
        return None  # Corrupt header — leave it alone, append still works.
    finally:
        os.close(fd)
    if _isocalendar_key(first_ts) == _isocalendar_key(now):
        return None
    base = f"audit-{_isocalendar_key(first_ts)}"
    rotated = path.with_name(f"{base}.log")
    n = 2
    while rotated.exists() or os.path.islink(rotated):
        rotated = path.with_name(f"{base}-{n}.log")
        n += 1
    try:
        path.rename(rotated)
    except OSError:
        return None
    try:
        if head_path(path).exists():
            head_path(path).rename(head_path(rotated))
    except OSError:
        pass
    return _line_hash(last) if last is not None else None


def _append_chained(path: Path, enriched: dict, now: datetime) -> str:
    """Write ``enriched`` as the next link of the chain; return the line.

    Under ``flock`` on the lock file beside the log, so two processes
    that append to one log extend one chain."""
    lock_fd = _open_private(path.with_name(path.name + ".lock"),
                            os.O_RDWR | os.O_CREAT)
    try:
        if _fcntl is not None:
            try:
                _fcntl.flock(lock_fd, _fcntl.LOCK_EX)
            except OSError:
                pass        # a file system without locks: still written
        link = _rotate_if_needed(path, now)
        fd = _open_private(path, os.O_RDWR | os.O_APPEND | os.O_CREAT)
        try:
            size = os.fstat(fd).st_size
            last = _last_line(fd, size)
            prev = _line_hash(last) if last is not None else (link or GENESIS)
            head = _read_head(path)
            rec = dict(enriched)
            rec.pop("prev", None)
            rec.pop("chain_break", None)
            rec.pop("mac", None)
            if head is not None:
                expected = prev if last is not None else None
                if head.get("last") != expected or head.get("size") != size:
                    rec["chain_break"] = {
                        "head_last": head.get("last"),
                        "head_size": head.get("size"), "found_size": size}
            rec["prev"] = prev
            line = _sign_line(json.dumps(rec, ensure_ascii=False),
                              _audit_key())
            data = line.encode("utf-8") + b"\n"
            if size and os.pread(fd, 1, size - 1) != b"\n":
                data = b"\n" + data     # a torn line stays on its own line
            view = memoryview(data)
            while view:
                view = view[os.write(fd, view):]
            _write_head(path, _line_hash(line.encode("utf-8")),
                        size + len(data))
            try:
                _write_anchor(path, _line_hash(line.encode("utf-8")),
                              size + len(data))
            except Exception:
                pass
        finally:
            os.close(fd)
    finally:
        os.close(lock_fd)
    if "chain_break" in rec:
        _note_event("audit_chain_broken",
                    f"{path}: the log was changed outside DELFIN "
                    f"(head {rec['chain_break']['head_size']} bytes, "
                    f"log {size} bytes)")
    return line + "\n"


def _note_event(kind: str, detail: str) -> None:
    """A blocked security event (panel + attention inbox). Never raises."""
    try:
        from . import security_events
        security_events.record(kind, "audit", detail[:300], blocked=True)
    except Exception:
        pass


def _note_write_failure(path: Path, exc: BaseException) -> None:
    global _WRITE_FAILURE_NOTED
    if _WRITE_FAILURE_NOTED:
        return
    _WRITE_FAILURE_NOTED = True
    _note_event("audit_write_failed",
                f"{path}: {type(exc).__name__}: {exc}. Records of this "
                "process are not in the audit log until it can be written.")


def append(
    record: dict,
    *,
    log_path: Optional[Path] = None,
    session_id: str = "",
) -> None:
    """Append one JSON record to the audit log.

    ``record`` is shallow-copied and decorated with the timestamp; the
    caller owns the rest of the schema. ``session_id`` (optional) stamps
    the record with the engine session when the record itself carries
    none — callers that already put ``session_id`` into the record keep
    their value. The record becomes the next link of the hash chain
    (``prev``; see the module docstring). A failure never blocks or
    breaks the tool path it observes, and is not silent either: the
    first one of a process records a blocked ``audit_write_failed``
    security event.
    """
    path = Path(log_path or _default_log_path())
    now = _now()
    enriched = {"ts": _iso(now), **record}
    if session_id and not enriched.get("session_id"):
        enriched["session_id"] = session_id
    line = json.dumps(enriched, ensure_ascii=False) + "\n"
    try:
        from .state_paths import ensure_dir
        with _LOG_LOCK:
            ensure_dir(path.parent)
            # Created 0600 (``_open_private``): a chmod after the first
            # write leaves the file group-readable for the length of that
            # write. A rotated log keeps the mode it had when it was
            # renamed -- the start-up repair sweep is what fixes those.
            line = _append_chained(path, enriched, now)
    except Exception as exc:
        # The action proceeds without logging (a removed home, a full
        # disk, a link in place of the log) -- and the panel says so.
        _note_write_failure(path, exc)
    # The administrator's destination (managed_policy ``audit_log``), in
    # addition to the user's own log; never rotated here.
    try:
        from .managed_policy import append_audit_line
        append_audit_line(line)
    except Exception:
        pass


def _predecessor(path: Path) -> Optional[Path]:
    """The log ``path`` continues: for the active log the newest rotated
    one, for a rotated log the one rotated before it."""
    try:
        rotated = sorted(p for p in path.parent.glob("audit-*.log")
                         if not p.is_symlink())
    except OSError:
        return None
    if path not in rotated:
        return rotated[-1] if rotated else None
    i = rotated.index(path)
    return rotated[i - 1] if i > 0 else None


def verify(log_path: Optional[Path] = None) -> dict:
    """Check one log's hash chain and head. Never raises.

    Output: ``{"path", "ok", "records", "chained", "legacy", "problems",
    "notes"}``; each problem is ``{"line", "kind", "detail"}`` (1-based
    line). Kinds: ``unparseable`` (a torn or edited line), ``chain`` (a
    line's ``prev`` does not match the line before it: that line was
    edited, or lines were removed or inserted), ``unchained`` (a line
    without ``prev`` after the chain started), ``link`` (the first line
    does not continue the log it names), ``head`` (the log does not end
    where the last write left it: lines were cut off or appended outside
    DELFIN), ``head_at_write`` (a write found the head out of step and
    said so in ``chain_break``). ``ok`` is True when there is no
    problem; notes say what could not be checked.
    """
    path = Path(log_path or _default_log_path())
    rep: dict = {"path": str(path), "ok": True, "records": 0, "chained": 0,
                 "legacy": 0, "keyed": 0, "unkeyed": 0, "problems": [],
                 "notes": []}
    try:
        from . import trust_store
        key = trust_store.derived_key(_PURPOSE, create=False)
    except Exception:
        key = None

    def problem(line: int, kind: str, detail: str) -> None:
        rep["problems"].append({"line": line, "kind": kind, "detail": detail})

    try:
        data = _read_private(path)
    except FileNotFoundError:
        if head_path(path).exists():
            problem(0, "head", "the log is missing but its head file exists")
        else:
            rep["notes"].append("no log at this path")
        rep["ok"] = not rep["problems"]
        return rep
    except OSError as exc:
        problem(0, "unreadable", f"{type(exc).__name__}: {exc}")
        rep["ok"] = False
        return rep
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    prev_hash: Optional[str] = None
    started = False
    keyed_seen = False
    for i, raw in enumerate(lines, 1):
        try:
            rec = json.loads(raw.decode("utf-8"))
            if not isinstance(rec, dict):
                raise ValueError("not an object")
        except Exception:
            problem(i, "unparseable", "the line is not a JSON record")
            prev_hash = _line_hash(raw)
            continue
        rep["records"] += 1
        if "prev" not in rec:
            if started:
                problem(i, "unchained", "a record without prev after the "
                                        "chain started")
            else:
                rep["legacy"] += 1
            prev_hash = _line_hash(raw)
            continue
        started = True
        rep["chained"] += 1
        if "mac" in rec:
            keyed_seen = True
            rep["keyed"] += 1
            m = _MAC_TAIL.search(raw)
            if m is None:
                problem(i, "mac", "the signature is not the record's last "
                        "field: the line was rewritten")
            elif key is not None and not hmac.compare_digest(
                    hmac.new(key, raw[:m.start()] + b"}",
                             hashlib.sha256).hexdigest(),
                    m.group(1).decode()):
                problem(i, "mac", "the signature does not match: the record "
                        "was changed by a process without the account key")
        elif keyed_seen:
            rep["unkeyed"] += 1
            problem(i, "unkeyed", "a record without a signature after the "
                    "signed chain started: written by a process without the "
                    "account key (one started inside an agent session) or "
                    "outside DELFIN")
        else:
            rep["unkeyed"] += 1
        got = str(rec.get("prev") or "")
        if prev_hash is not None:
            if got != prev_hash:
                problem(i, "chain", f"prev does not match line {i - 1}: that "
                        "line was edited, or lines were removed or inserted")
        elif got != GENESIS:
            pred = _predecessor(path)
            if pred is None:
                rep["notes"].append("the first record continues a log that "
                                    "is not beside this one")
            else:
                try:
                    pdata = _read_private(pred).rstrip(b"\n")
                    plast = pdata.rsplit(b"\n", 1)[-1]
                except OSError:
                    plast = None
                if plast is None or _line_hash(plast) != got:
                    problem(i, "link", f"the first record does not continue "
                            f"{pred.name}")
        if isinstance(rec.get("chain_break"), dict):
            problem(i, "head_at_write", "when this record was written the "
                    "log did not end where the last write left it: "
                    + json.dumps(rec["chain_break"]))
        prev_hash = _line_hash(raw)
    head = _read_head(path)
    if head is None:
        if started:
            rep["notes"].append("no head file: a cut tail cannot be shown")
    elif head.get("last") != prev_hash or head.get("size") != len(data):
        problem(len(lines), "head", "the log does not end where the last "
                f"write left it ({head.get('size')} bytes then, {len(data)} "
                "now): lines were cut off, or written outside DELFIN")
    anchor = _read_anchor(path)
    if anchor is not None:
        a_last, a_size = anchor.get("last"), anchor.get("size")
        signed = (key is None or hmac.compare_digest(
            hmac.new(key, _anchor_bytes(path, a_last, a_size or 0),
                     hashlib.sha256).hexdigest(), str(anchor.get("mac") or "")))
        if not signed:
            problem(len(lines), "anchor", "the trust store's anchor for this "
                    "log does not carry a valid signature")
        elif a_last != prev_hash or a_size != len(data):
            problem(len(lines), "anchor", "the log does not end where the "
                    f"trust store's anchor says ({a_size} bytes then, "
                    f"{len(data)} now): it was rewritten, cut or extended "
                    "by a process without the account key")
    elif rep["keyed"] and _is_active_log(path):
        rep["notes"].append("the trust store holds no anchor for this log: "
                            "a cut tail with an earlier head file put back "
                            "is not shown")
    if key is None and rep["keyed"]:
        rep["notes"].append("the account key is not available here (inside "
                            "an agent session, or the trust store is "
                            "refused): the signatures were not checked")
    rep["ok"] = not rep["problems"]
    return rep


def _is_active_log(path: Path) -> bool:
    """Whether ``path`` is an active log (not a rotated week), the kind
    that has an anchor."""
    return not Path(path).name.startswith("audit-")


def verify_all(directory: Optional[Path] = None) -> list[dict]:
    """``verify`` for the active log and every rotated one beside it."""
    base = Path(directory) if directory else _default_log_path().parent
    paths = sorted(p for p in base.glob("audit-*.log") if not p.is_symlink())
    active = base / "audit.log"
    if active.exists() or head_path(active).exists():
        paths.append(active)
    return [verify(p) for p in paths]


def format_verify(reports: list[dict]) -> str:
    """Plain text for ``delfin-agent audit verify``."""
    if not reports:
        return "no audit log found"
    out = []
    for rep in reports:
        state = "intact" if rep["ok"] else "NOT intact"
        out.append(f"{rep['path']}: {state} -- {rep['records']} records "
                   f"({rep['chained']} chained, {rep.get('keyed', 0)} signed, "
                   f"{rep['legacy']} legacy)")
        for p in rep["problems"][:20]:
            out.append(f"  line {p['line']}: {p['kind']}: {p['detail']}")
        if len(rep["problems"]) > 20:
            out.append(f"  ... {len(rep['problems']) - 20} more")
        for n in rep["notes"]:
            out.append(f"  note: {n}")
    return "\n".join(out)


def read_last_n(n: int = 100, *, log_path: Optional[Path] = None) -> list[dict]:
    """Return the last ``n`` parseable records (utility for the dashboard)."""
    path = log_path or _default_log_path()
    if not path.exists():
        return []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return []
    out: list[dict] = []
    for line in lines[-(n * 2):]:        # over-read a little for parse errors
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except Exception:
            continue
    return out[-n:]


def _diff_line_counts(old: str, new: str) -> tuple[int, int]:
    """Quick added / removed line counts for the audit record.

    Not a real diff — just compares line sets to avoid bringing
    difflib into a hot path. Good enough for "approximately how big
    was this change" telemetry.
    """
    old_lines = old.splitlines() if old else []
    new_lines = new.splitlines() if new else []
    added = max(0, len(new_lines) - len(old_lines))
    removed = max(0, len(old_lines) - len(new_lines))
    # If the line counts match, count the actually-different lines.
    if added == 0 and removed == 0 and old_lines != new_lines:
        diffs = sum(1 for a, b in zip(old_lines, new_lines) if a != b)
        return diffs, diffs
    return added, removed


def make_record(
    *,
    tool: str,
    decision: str,
    mode: str = "",
    path: str = "",
    command: str = "",
    diff_lines_added: int = 0,
    diff_lines_removed: int = 0,
    persistence: Optional[dict] = None,
    session_id: str = "",
    reason: str = "",
    extra: Optional[dict[str, Any]] = None,
) -> dict:
    """Construct a record dict with the standard fields.

    ``reason`` carries WHY a call was refused. Without it a denial read
    "bash: <command> [denied]" and nothing more -- so after an unattended
    run that stopped having produced nothing, the user could not tell
    whether it was the deny-list, the secret scan, the egress scan, the
    workspace boundary, an auto-allow miss or a hook, and therefore could
    not know what to allow. The gate's own message was handed to the model
    and then discarded; the only place the rule name survived was a
    process-local deque that dies with the process.
    """
    rec: dict[str, Any] = {
        "session_id": session_id,
        "tool": tool,
        "decision": decision,
        "mode": mode,
        "pid": os.getpid(),
    }
    if reason:
        rec["reason"] = str(reason)[:300]
    if path:
        rec["path"] = path
    if command:
        rec["command"] = command[:500]
    if diff_lines_added or diff_lines_removed:
        rec["diff_lines_added"] = diff_lines_added
        rec["diff_lines_removed"] = diff_lines_removed
    if persistence:
        rec["persistence"] = persistence
    if extra:
        rec.update(extra)
    return rec


# ---------------------------------------------------------------------------
# Changes report — the audit log's read side.
#
# "What did you change?" must be answered from the record, not from model
# memory. build_changes_report() aggregates the audit records of one
# session (or a time window) into a small dict; format_changes_report()
# renders that dict as compact markdown for chat display. Both are
# defensive end to end: torn/malformed lines are skipped and neither
# function ever raises.
# ---------------------------------------------------------------------------

_WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "notebook_edit",
    "apply_patch", "publish_report",
    # Document writes. The fillers record the blank form/template as
    # their ``path``; the file they produced is named in the result, so
    # the report shows what was worked on either way.
    "edit_sheet", "fill_pdf_form", "fill_docx_template", "create_docx",
    "fill_series", "merge_pdfs", "split_pdf", "create_pdf",
    # An undo overwrites user files too, and is itself undoable.
    "undo_changes",
    # Files copied in from a cluster (cluster_tools.fetch).
    "cluster_fetch",
})

# Write routes that reach the user's files WITHOUT a pre-image in the
# undo journal. They are listed by name so the gap is a declaration
# somebody has to update rather than something the report discovers by
# being wrong: a file written through one of these can be shown, but it
# cannot be offered as undoable.
#
# * ``bash_background`` — arbitrary shell, and unlike the foreground
#   one it returns before the command has written anything, so there is
#   no moment at which a pre-image and a post-image can be compared.
# * ``fill_series`` — writes one document per table row, none captured.
# * ``publish_report`` — writes its own output tree.
# * MCP file mutators arrive under the BARE native tool name
#   (``write_file`` …) and dispatch to a remote server, so they are
#   indistinguishable in the record from a native write; the journal
#   lookup below is what separates them.
#
# ``bash`` came OFF this list on 2026-09-10. The stated reason for it
# being here — "the write targets are not known before the command runs"
# — was not true: ``_bash_write_targets`` is what the write GATE decides
# on, so the targets are known, and the foreground executor now reads
# each one before running and journals the ones whose bytes moved.
#
# What the target scanner cannot see -- a program the shell RUNS that
# writes files of its own (``python3 build.py``) -- is journalled from a
# snapshot of the writable roots taken around the command
# (``shell_changes``), for ``bash`` and, when the job ends, for
# ``bash_background``, which came off this list for that reason. Its
# limits are per file, not per tool: a file whose pre-image was not held
# is a record without undo, and a root over the file cap is not tracked
# (the command's result says so); the per-file lookup below reports
# both as unrecoverable.
#
# * ``cluster_fetch`` -- copies files in from a cluster, never over an
#   existing file, so there is no pre-image to keep; delete what it wrote.
# * ``cluster_run``, ``cluster_submit``, ``cluster_cancel``,
#   ``cluster_upload`` -- change files and jobs on a cluster, not here: the
#   journal holds nothing of them and undoes none of it.
_UNJOURNALLED_WRITE_TOOLS = frozenset({
    "fill_series", "publish_report", "cluster_fetch",
    "cluster_run", "cluster_submit", "cluster_cancel", "cluster_upload",
})
# The cluster tools act on another machine: their rows carry a command line
# that says what was done there (cluster_tools._audit).
_COMMAND_TOOLS = frozenset({"bash", "bash_background", "cluster_run", "cluster_submit",
                            "cluster_cancel", "cluster_upload"})
_PERSIST_TOOLS = frozenset({"remember_permission", "remember_permission_bundle"})
# A hook is a shell command. It runs outside the permission gate, outside
# the deny-list and outside the sandbox, and it fires on every prompt and
# every tool call -- so it belongs in "what did you change" more than most
# of what is already here. Only the BLOCKING ones appeared, via
# _DENIED_DECISIONS below: a hook that successfully executed matched none
# of the write/command/persist families and fell out of the report
# entirely. The protective refusal was reported and the arbitrary
# execution was not.
_HOOK_TOOLS = frozenset({"hook"})
# "block" is what a PreToolUse hook writes. It was missing here, so a
# hook that stopped an edit produced a record that matched none of the
# write/command/persist families either -- and dropped out of the changes
# report entirely. A refusal nobody can see is a refusal nobody can undo.
_DENIED_DECISIONS = frozenset({"denied", "error", "block"})


def _parse_log(path: Path) -> list[dict]:
    """Every parseable record of one log file, in file order."""
    if not path.exists():
        return []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return []
    out: list[dict] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except Exception:
            continue
        if isinstance(rec, dict):
            out.append(rec)
    return out


def latest_rotated_log(path: Path) -> Optional[Path]:
    """The most recent ``audit-YYYY-Www.log`` beside ``path``, if any."""
    try:
        rotated = sorted(path.parent.glob("audit-*.log"))
    except Exception:
        return None
    return rotated[-1] if rotated else None


def _read_all_records(
    log_path: Optional[Path] = None, *, with_rotated: bool = False,
) -> tuple[list[dict], str]:
    """Records oldest first, plus the name of the rotated file consulted.

    The weekly rotation keeps the active file small, so a full parse is
    cheap. It also means that on the first day of a new week the active
    log is EMPTY: reporting "no recorded changes" from it alone told the
    user their session had changed nothing, hours after it had. When
    ``with_rotated`` is set the immediately preceding rotated file is
    read as well, and its name is returned so the report can say which
    history it looked at.
    """
    path = log_path or _default_log_path()
    out: list[dict] = []
    consulted = ""
    if with_rotated:
        prev = latest_rotated_log(path)
        if prev is not None:
            out.extend(_parse_log(prev))
            consulted = prev.name
    out.extend(_parse_log(path))
    return out, consulted


def _same_dir(a: str, b: str) -> bool:
    """The same directory, whatever spelling each side used."""
    try:
        return os.path.normcase(os.path.realpath(os.path.expanduser(a))) == \
            os.path.normcase(os.path.realpath(os.path.expanduser(b)))
    except Exception:
        return False


def _under_workspace(target: str, workspace: str) -> bool:
    """True unless ``target`` is an absolute path outside ``workspace``.

    Relative paths and empty targets pass — audit records store
    workspace-relative paths for file tools, so only absolute strays
    (e.g. bash cwd escapes into a granted extra dir) are filtered.
    """
    if not target or not workspace:
        return True
    try:
        t = Path(target)
        if not t.is_absolute():
            return True
        t.relative_to(Path(workspace).resolve())
        return True
    except ValueError:
        return False
    except Exception:
        return True


def _journal_paths(
    session_id: Optional[str], journal_paths: Optional[set] = None,
) -> Optional[set]:
    """Absolute paths this session has a restorable pre-image for.

    ``None`` means "not knowable" (no session given / the journal could
    not be read), which the report renders as an unknown rather than as
    a promise it cannot keep.
    """
    if journal_paths is not None:
        return {_norm_path(p) for p in journal_paths}
    if not session_id:
        return None
    try:
        from . import change_journal as _cj
        out = set()
        for rec in _cj.list_changes(session_id) or []:
            if (rec.get("dropped") or rec.get("truncated")
                    or rec.get("lossy")):
                continue        # recorded, but no pre-image to restore
            p = str(rec.get("path", "") or "")
            if p:
                out.add(_norm_path(p))
        return out
    except Exception:
        return None


def _command_entry(rec: dict, command: str) -> dict:
    """One "Commands run" entry, carrying the exit code when known.

    Older records have no ``exit_code`` (the audit only ever described
    the gate's verdict on the attempt); the key is left out entirely
    then, rather than reported as a zero the log never observed.
    """
    entry = {"command": command, "cwd": str(rec.get("cwd", ""))}
    code = rec.get("exit_code")
    if isinstance(code, int):
        entry["exit_code"] = code
    return entry


def _norm_path(path: str) -> str:
    try:
        return os.path.normcase(os.path.abspath(os.path.expanduser(str(path))))
    except Exception:
        return str(path)


def build_changes_report(
    session_id: Optional[str] = None,
    *,
    since_ts: Optional[str] = None,
    last_n: int = 200,
    workspace: Optional[str] = None,
    log_path: Optional[Path] = None,
    journal_paths: Optional[set] = None,
) -> dict:
    """Aggregate audit records into a "what changed" report dict.

    Filtering: when ``session_id`` is given, only that session's records
    count (plus ``since_ts`` if also given); otherwise the time window
    ``since_ts`` applies (ISO strings compare lexicographically). At most
    the newest ``last_n`` matching records are aggregated. ``workspace``
    (optional) drops records whose absolute path/cwd lies outside it.

    Each written file carries ``undoable``, decided by looking the path
    up in the undo JOURNAL rather than by trusting the tool's name. A
    shell command, a series fill or an MCP write appears in the record
    exactly like a native write, so a name-based answer promised an undo
    the journal cannot perform; ``undoable: False`` says so on the line
    itself, and ``unjournalled`` collects those paths.

    Returns::

        {"files_written":        [{"path", "tool", "count", "undoable"}],
         "commands":             [{"command", "cwd", "exit_code"}],
         "hooks_fired":          [{"command", "source", "exit_code",
                                   "event", "decision"}],
         "denied":               [{"tool", "target", "decision", "reason"}],
         "permissions_persisted":[{"tool", "persistence"}],
         "unjournalled":         [paths written with no pre-image],
         "window": {"from_ts", "to_ts", "records", "rotated"}}

    Never raises — on any failure the empty skeleton comes back.
    """
    empty: dict[str, Any] = {
        "files_written": [], "commands": [], "hooks_fired": [], "denied": [],
        "permissions_persisted": [], "unjournalled": [],
        "window": {"from_ts": "", "to_ts": "", "records": 0, "rotated": ""},
    }
    try:
        undoable_paths = _journal_paths(session_id, journal_paths)

        def _matching(records: list[dict]) -> list[dict]:
            out: list[dict] = []
            for rec in records:
                if session_id is not None and session_id != "":
                    if rec.get("session_id", "") != session_id:
                        continue
                if since_ts:
                    if str(rec.get("ts", "")) < str(since_ts):
                        continue
                if workspace:
                    # A record that names its workspace is judged by it,
                    # exactly. Only records from before that field existed
                    # fall back to the path test below -- which cannot
                    # tell a relative path of ours from another
                    # workspace's, and let a benchmark's commands into a
                    # probe's report as its own.
                    tagged = str(rec.get("workspace", "") or "")
                    if tagged:
                        if not _same_dir(tagged, str(workspace)):
                            continue
                    else:
                        target = str(rec.get("path", "") or rec.get("cwd", ""))
                        if not _under_workspace(target, str(workspace)):
                            continue
                out.append(rec)
            return out

        records, rotated = _read_all_records(log_path)
        picked = _matching(records)
        if not picked:
            # An empty answer from the active log alone is the normal
            # state on the first day of a new ISO week, and reporting it
            # as "no recorded changes" is how a whole week of history
            # became invisible with no notice that rotation happened.
            records, rotated = _read_all_records(log_path, with_rotated=True)
            picked = _matching(records)
        try:
            cap = max(1, int(last_n))
        except (TypeError, ValueError):
            cap = 200
        picked = picked[-cap:]

        files: dict[str, dict[str, Any]] = {}
        file_tools: dict[str, set] = {}
        commands: list[dict[str, str]] = []
        hooks_fired: list[dict[str, Any]] = []
        denied: list[dict[str, str]] = []
        persisted: list[dict[str, Any]] = []
        for rec in picked:
            tool = str(rec.get("tool", ""))
            decision = str(rec.get("decision", ""))
            path = str(rec.get("path", ""))
            command = str(rec.get("command", ""))
            # Before the denied branch, and without consuming the record:
            # a blocking hook both RAN and refused, and belongs in both
            # places. Its exit code is what the user needs to see, and
            # its source is what tells them whether it was theirs.
            if tool in _HOOK_TOOLS:
                hooks_fired.append({
                    "command": command,
                    "source": str(rec.get("source", "")),
                    "event": str(rec.get("hook_event", "")),
                    "exit_code": rec.get("exit_code", 0),
                    "decision": decision,
                })
            if decision in _DENIED_DECISIONS:
                denied.append({
                    "tool": tool,
                    "target": path or command,
                    "decision": decision,
                    "reason": str(rec.get("reason", "") or ""),
                })
                # A command that ran and FAILED is still a command that
                # ran; dropping it from "Commands run" would hide it.
                if tool in _COMMAND_TOOLS and command:
                    commands.append(_command_entry(rec, command))
                continue
            if tool in _WRITE_TOOLS and path:
                entry = files.setdefault(
                    path, {"path": path, "tool": "", "count": 0})
                entry["count"] += 1
                file_tools.setdefault(path, set()).add(tool)
            elif tool in _COMMAND_TOOLS and command:
                commands.append(_command_entry(rec, command))
            if tool in _PERSIST_TOOLS or rec.get("persistence"):
                persisted.append({
                    "tool": tool,
                    "persistence": rec.get("persistence") or {},
                })
        unjournalled: list[str] = []
        for path, tools in file_tools.items():
            files[path]["tool"] = "+".join(sorted(tools))
            # Undoability comes from the journal, never from the tool
            # name: the same name covers a native write (journalled) and
            # an MCP write to a remote server (not journalled).
            if undoable_paths is None:
                files[path]["undoable"] = None
            else:
                ok = _norm_path(path) in undoable_paths
                files[path]["undoable"] = ok
                if not ok:
                    unjournalled.append(path)
        return {
            "files_written": list(files.values()),
            "commands": commands,
            "hooks_fired": hooks_fired,
            "denied": denied,
            "permissions_persisted": persisted,
            "unjournalled": unjournalled,
            "window": {
                "from_ts": str(picked[0].get("ts", "")) if picked else "",
                "to_ts": str(picked[-1].get("ts", "")) if picked else "",
                "records": len(picked),
                "rotated": rotated,
            },
        }
    except Exception:
        return empty


def format_changes_report(report: dict) -> str:
    """Render a build_changes_report() dict as compact markdown.

    Chat-display companion of the report dict; tolerates missing keys
    and never raises.
    """
    try:
        files = report.get("files_written") or []
        commands = report.get("commands") or []
        hooks_fired = report.get("hooks_fired") or []
        denied = report.get("denied") or []
        persisted = report.get("permissions_persisted") or []
        window = report.get("window") or {}
        rotated = str(window.get("rotated", "") or "")
        if not (files or commands or hooks_fired or denied or persisted):
            tail = (f" Rotated history ({rotated}) was consulted too."
                    if rotated else
                    " Rotated weekly history was NOT consulted.")
            return ("No recorded changes — the audit log has no matching "
                    "records for this session/window." + tail)
        lines: list[str] = ["### Changes made"]
        if files:
            lines.append("Files written:")
            for f in files:
                undoable = f.get("undoable")
                mark = ("" if undoable is None else
                        "" if undoable else "  — NOT undoable (no pre-image)")
                lines.append(
                    f"- `{f.get('path', '?')}` — {f.get('tool', '?')}"
                    f" ×{f.get('count', 0)}{mark}")
        if commands:
            lines.append("Commands run:")
            for c in commands:
                cwd = c.get("cwd", "")
                suffix = f" (cwd: {cwd})" if cwd else ""
                code = c.get("exit_code")
                if isinstance(code, int) and code != 0:
                    suffix += f" [exit {code}]"
                lines.append(f"- `{c.get('command', '?')}`{suffix}")
        if hooks_fired:
            lines.append("Hooks fired:")
            for h in hooks_fired:
                src = str(h.get("source", "") or "").strip()
                lines.append(
                    f"- `{h.get('command', '?')}` — exit="
                    f"{h.get('exit_code', 0)}"
                    + (f", from {src}" if src else ", source unrecorded"))
        if denied:
            lines.append("Denied / failed:")
            for d in denied:
                # The reason is the whole point of the line. Without it
                # the report says something was refused and leaves the
                # user unable to find out which rule to relax.
                _why = str(d.get("reason", "") or "").strip()
                lines.append(
                    f"- {d.get('tool', '?')}: {d.get('target', '?')}"
                    f" [{d.get('decision', 'denied')}]"
                    + (f" — {_why}" if _why else ""))
        if persisted:
            lines.append("Permissions persisted:")
            for p in persisted:
                lines.append(f"- {p.get('tool', '?')}: "
                             f"{json.dumps(p.get('persistence') or {}, ensure_ascii=False)}")
        unjournalled = report.get("unjournalled") or []
        if unjournalled:
            # The list above is what CHANGED; this is the part an undo
            # cannot take back. Saying it here is the difference between
            # "undo that" restoring everything and restoring a subset
            # the user was never told about.
            lines.append(
                f"Not undoable ({len(unjournalled)}) — written without a "
                "pre-image in the undo journal (shell command, series "
                "fill, MCP server or an earlier session):")
            for p in unjournalled[:20]:
                lines.append(f"- `{p}`")
            if len(unjournalled) > 20:
                lines.append(f"- … {len(unjournalled) - 20} more")
        n = window.get("records", 0)
        frm, to = window.get("from_ts", ""), window.get("to_ts", "")
        if n:
            span = f" ({frm} → {to})" if frm and to and frm != to else (
                f" ({frm})" if frm else "")
            rot = f", incl. rotated {rotated}" if rotated else ""
            lines.append(f"_{n} audit record(s){span}{rot}_")
        return "\n".join(lines)
    except Exception:
        return "No recorded changes — report unavailable."


def steps_since(offset: int = 0, *, log_path: Optional[Path] = None,
                session: str = "", denied_only: bool = False
                ) -> "tuple[list[dict], int]":
    """Records written after *offset*, and the offset to pass next time.

    An offset rather than a count, because a follower that asks for "the
    last 50" prints some records twice and misses others when a burst
    arrives between two looks. The offset is a byte position in the log:
    what has been read is read.

    ``denied_only`` keeps what the gate REFUSED, which is the line a
    supervisor stops on. ``error`` is deliberately not that -- a grep
    that finds nothing exits 1, and reading those as refusals is how a
    first draft of the watcher reported 46 false alarms in one session.

    Never raises: a supervisor's view of the world must not be the thing
    that ends the supervisor.
    """
    path = log_path or _default_log_path()
    out: list[dict] = []
    try:
        size = path.stat().st_size
    except OSError:
        return out, 0
    start = int(offset or 0)
    if start > size:
        start = 0            # rotated out from under us: begin again
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            fh.seek(start)
            text = fh.read()
            end = fh.tell()
    except OSError:
        return out, start
    # A line still being written has no newline yet; leave it for next time.
    if text and not text.endswith("\n"):
        keep = text.rfind("\n")
        if keep < 0:
            return out, start
        end = start + len(text[:keep + 1].encode("utf-8", "replace"))
        text = text[:keep + 1]
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except Exception:
            continue
        if not isinstance(record, dict):
            continue
        if session and not str(record.get("session_id") or "").startswith(session):
            continue
        if denied_only and str(record.get("decision") or "") != "denied":
            continue
        out.append(record)
    return out, end


def render_step(record: dict) -> str:
    """One step, as a supervisor reads it. Never raises."""
    try:
        record = record or {}
        ts = str(record.get("ts") or "")[11:19] or "--:--:--"
        session = str(record.get("session_id") or "")[:8] or "--------"
        tool = str(record.get("tool") or "?")[:22]
        what = str(record.get("command") or record.get("path") or "")
        what = " ".join(what.split())[:96]
        decision = str(record.get("decision") or "")
        if decision == "denied":
            reason = " ".join(str(record.get("reason") or "").split())[:120]
            return f"{ts} {session} DENIED {tool}: {what}\n         └─ {reason}"
        mark = "!" if decision == "error" else " "
        return f"{ts} {session} {mark}      {tool}: {what}"
    except Exception:
        return ""
