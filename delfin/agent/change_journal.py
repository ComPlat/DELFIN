"""Undo journal for agent file edits — pre-image capture + safe revert.

The agent's write path (write_file / edit_file / apply_patch /
notebook_edit) reads each file's old content before overwriting it (for
the chat diff) and then throws it away. Calc workspaces frequently are
NOT git repos, so a single agent edit clobbers uncommitted user work
irrecoverably. This module keeps that pre-image:

Storage layout (per session, under ``~/.delfin/undo/<sid>/``)::

    journal.jsonl               one JSON record per captured change
    <seq:06d>-<prehash8>.txt    full pre-image, text   (content-addressed)
    <seq:06d>-<prehash8>.bin    full pre-image, binary (content-addressed)

Journal record schema::

    {seq, ts, tool, path, pre_hash, post_hash, pre_file,
     created: bool,      # True → the file did not exist before the write
     deleted: bool,      # True → the write REMOVED the file (post = absence)
     truncated: bool,    # True → pre-image exceeded MAX_PRE_IMAGE_BYTES
     lossy: bool,        # True → pre-image did not fit the file's encoding
     binary: bool,       # True → pre-image and hashes are raw bytes
     raw: bool,          # True → byte-exact pre-image, raw-byte hashes
     dropped: bool,      # True → tombstone: pre-image dropped at the cap
     undone: {...},      # set once this record has been reverted
     undo_of: int}       # this record was written BY a revert

Text and binary changes are recorded by separate entry points
(``record_change`` / ``record_binary_change``) because the *inputs*
differ: the text write path hands over the normalised text it edited,
the document writers hand over raw bytes. What is STORED is byte-exact
in both cases. ``record_change`` re-encodes the pre-image in the file's
own convention (encoding / line ending / BOM, read back from the file
the write just produced) and hashes the file's RAW bytes, because the
guard in :func:`revert` compares against raw bytes: hashing the
normalised text instead made undo structurally unavailable for every
CRLF, BOM or non-UTF-8 file and blamed the user for a change nobody
made. Records carrying neither ``raw`` nor ``binary`` are legacy text
records and keep reverting exactly as they did.

Guarantees:

* ``record_change`` never raises — a broken undo store must not block the
  tool path it observes (returns ``None`` on any failure).
* Pre-images and journal rewrites are atomic (temp file + ``os.replace``)
  and 0600; the per-session directory is 0700.
* Pre-images larger than ``MAX_PRE_IMAGE_BYTES`` (2 MB) are stored
  TRUNCATED and flagged; :func:`revert` refuses those entries (a partial
  restore would silently corrupt the file — losing undo for a huge file
  is better than "undoing" it into garbage).
* Per-session cap of ``MAX_ENTRIES_PER_SESSION`` (500) entries; the
  oldest pre-images are unlinked, but a TOMBSTONE stays in the journal
  so a later undo can say the pre-image was dropped at the cap instead
  of returning the same three empty lists as "there was nothing to
  undo". The drop is also reported to the caller that caused it.
* ``revert`` only ever touches a file whose CURRENT sha256 equals the
  recorded ``post_hash``. Anything else — the user or another process
  changed the file since the agent's edit — is reported as a conflict
  and the file is left untouched.
* ``revert`` JOURNALS ITSELF: every file it restores, deletes or
  recreates is written back as a new record (``tool="undo"``), so an
  undo is itself undoable and shows up in the change listing. The
  record it acted on is marked ``undone`` so a second run reports
  "already undone" instead of blaming the user for the agent's own
  undo. Undo records are reachable through ``scope="turn"``; the
  ``last`` and ``session`` scopes walk back through the agent's own
  changes and skip them.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Optional

MAX_PRE_IMAGE_BYTES = 2 * 1024 * 1024   # 2 MB
MAX_ENTRIES_PER_SESSION = 500

_JOURNAL_NAME = "journal.jsonl"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha256_text(text: str) -> str:
    """sha256 hex digest of ``text`` (utf-8)."""
    return hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest()


def sha256_bytes(data: bytes) -> str:
    """sha256 hex digest of raw bytes."""
    return hashlib.sha256(data).hexdigest()


def sha256_file_bytes(path: Path) -> Optional[str]:
    """sha256 of a file's raw bytes, or None if unreadable.

    Distinct from :func:`sha256_file`, which hashes the lossy utf-8
    decoding used by the text write path. A spreadsheet or PDF round-trips
    only as bytes: decoding it with ``errors="replace"`` maps whole ranges
    of byte values onto U+FFFD, so a text hash cannot tell two different
    binaries apart, and a text pre-image cannot restore either of them.
    """
    try:
        return sha256_bytes(path.read_bytes())
    except OSError:
        return None


def sha256_file(path: Path) -> Optional[str]:
    """sha256 hex digest of a file's content, or None if unreadable.

    Decodes the raw bytes (no universal-newline translation — text-mode
    reads would collapse ``\\r\\n`` to ``\\n`` and silently change the
    hash of CRLF files) with the same ``errors="replace"`` policy the
    write path uses when it captures ``old_text``/``new_text``."""
    try:
        text = path.read_bytes().decode("utf-8", errors="replace")
    except OSError:
        return None
    return sha256_text(text)


def _safe_session_id(session_id: Any) -> str:
    """Sanitize a session id for use as a directory name (no traversal).

    Same pattern as ``session_store.save_handoff_brief``: every character
    outside ``[a-zA-Z0-9_-]`` becomes ``_`` (this includes ``.`` and
    ``/``, so ``../../etc`` cannot escape the undo root)."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", str(session_id or "") or "session")[:40]


def _undo_root() -> Path:
    return Path.home() / ".delfin" / "undo"


def _session_dir(session_id: Any) -> Path:
    return _undo_root() / _safe_session_id(session_id)


def _journal_path(session_id: Any) -> Path:
    return _session_dir(session_id) / _JOURNAL_NAME


def _set_file_perms(path: Path) -> None:
    """Best-effort 0600 on a journal/pre-image file (per-user data)."""
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def _atomic_write(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` atomically (temp file + ``os.replace``).

    Same pattern as ``memory_store._atomic_write`` — a reader racing the
    writer must never observe an empty or torn file. ``newline=""``
    disables newline translation so pre-images round-trip byte-exact."""
    import tempfile
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        os.replace(tmp, path)
        _set_file_perms(path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _atomic_restore(path: Path, text: str) -> None:
    """Atomically replace a WORKSPACE file with ``text``, preserving its
    current permission bits (a restored user file must not silently become
    0600 like our internal store files)."""
    import tempfile
    mode: Optional[int] = None
    try:
        mode = os.stat(path).st_mode & 0o7777
    except OSError:
        pass
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".undo.tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        if mode is not None:
            try:
                os.chmod(tmp, mode)
            except OSError:
                pass
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    """Byte-exact counterpart of :func:`_atomic_write`."""
    import tempfile
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
        _set_file_perms(path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _atomic_restore_bytes(path: Path, data: bytes) -> None:
    """Byte-exact counterpart of :func:`_atomic_restore`."""
    import tempfile
    mode: Optional[int] = None
    try:
        mode = os.stat(path).st_mode & 0o7777
    except OSError:
        pass
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".undo.tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        if mode is not None:
            try:
                os.chmod(tmp, mode)
            except OSError:
                pass
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _read_journal(session_id: Any) -> list[dict]:
    """All parseable journal records, sorted by seq ascending.

    Torn/garbage lines (crash mid-append) are skipped, never fatal."""
    p = _journal_path(session_id)
    try:
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    records: list[dict] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except (json.JSONDecodeError, ValueError):
            continue
        if not isinstance(rec, dict) or not isinstance(rec.get("seq"), int):
            continue
        records.append(rec)
    records.sort(key=lambda r: r["seq"])
    return records


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _convention_of(raw: Optional[bytes]) -> tuple[str, str, bool]:
    """The (encoding, line ending, BOM) of the bytes on disk NOW.

    The write path normalises the text it hands to the journal (LF, no
    BOM, decoded) and re-applies the file's convention when it writes.
    The journal therefore has to put that convention back before it
    stores a pre-image, or the restore would rewrite every line ending
    in a Windows-authored file — and the hash of what it stored would
    never match the file it was taken from.

    Detection is delegated to ``text_files``, the module the write path
    itself uses, so the two cannot drift apart.
    """
    if raw is None:
        return "utf-8", "\n", False
    try:
        from . import text_files as _tf
        text = ""
        encoding = _tf._ENCODINGS[-1]
        for candidate in _tf._ENCODINGS:
            try:
                text = raw.decode(candidate)
                encoding = candidate
                break
            except UnicodeDecodeError:
                continue
        bom = text.startswith(_tf._BOM)
        if bom:
            text = text[len(_tf._BOM):]
        return encoding, _tf._detect_newline(text), bom
    except Exception:
        return "utf-8", "\n", False


def _encode_like(
    text: str, encoding: str, newline: str, bom: bool,
) -> tuple[bytes, bool]:
    """``text`` in the given convention, plus whether it fit exactly.

    ``False`` means the pre-image cannot be reproduced: the text handed
    over holds characters the file's own encoding has no room for, which
    is what a caller that read the file with ``errors="replace"`` hands
    over. Storing that anyway and calling it a pre-image is how an undo
    writes U+FFFD into a user's file and reports success — so the flag
    travels with the record and :func:`revert` refuses it.
    """
    out = text.replace("\n", newline) if newline != "\n" else text
    if bom:
        out = "﻿" + out
    try:
        return out.encode(encoding), True
    except (UnicodeEncodeError, LookupError):
        return out.encode("utf-8", errors="replace"), False


def record_change(
    session_id: Any,
    *,
    tool: str,
    path: Any,
    old_text: Optional[str],
    new_text: Optional[str],
) -> Optional[dict]:
    """Capture one file change: store the full pre-image + a journal line.

    Called AFTER the write, so the post state can be hashed from the
    file itself rather than from the text the caller believes it wrote.
    ``old_text is None`` means the file did not exist before the write
    (``created=True``; revert deletes it). Returns
    ``{"seq", "pre_hash", "post_hash", "dropped": [seq, ...]}`` on
    success and ``None`` on ANY failure — the journal must never break
    the write path it observes. ``dropped`` names the entries whose
    pre-image this record pushed out of the per-session cap; it is empty
    in the ordinary case.
    """
    try:
        target = Path(str(path))
        sdir = _session_dir(session_id)
        sdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(sdir, 0o700)
        except OSError:
            pass

        created = old_text is None
        pre_text = "" if created else str(old_text)

        # One read of the file the write just produced: it carries both
        # the post-image hash the revert guard compares against and the
        # convention the pre-image has to be put back in (the writer
        # re-applied that convention, so it is the same one).
        try:
            post_bytes: Optional[bytes] = target.read_bytes()
        except OSError:
            post_bytes = None
        encoding, newline, bom = _convention_of(post_bytes)
        pre_bytes, exact = _encode_like(pre_text, encoding, newline, bom)

        if post_bytes is not None:
            post_hash = sha256_bytes(post_bytes)
        else:
            # Unreadable right after the write (removed, permissions):
            # fall back to hashing what the caller says it wrote, in the
            # same convention, so the guard still has something exact.
            post_text = "" if new_text is None else str(new_text)
            post_hash = sha256_bytes(
                _encode_like(post_text, encoding, newline, bom)[0])

        extra: dict = {"raw": True, "encoding": encoding,
                       "newline": newline, "bom": bom}
        if not exact and not created:
            extra["lossy"] = True
        return _write_record(
            session_id, sdir, target,
            tool=tool, pre_bytes=pre_bytes, post_hash=post_hash,
            created=created, extra=extra,
        )
    except Exception:
        return None


def record_deletion(
    session_id: Any,
    *,
    tool: str,
    path: Any,
    pre_bytes: Optional[bytes],
) -> Optional[dict]:
    """Capture a write that REMOVED a file. Revert recreates it.

    The post state of a deletion is the file's absence, so there is no
    post-image to hash — ``post_hash`` is empty and :func:`revert`
    verifies that the path is still absent before it writes the
    pre-image back. Without this the capture loop held the only copy of
    the file in memory and skipped the journal, which is the one case
    where the journal IS the last copy. Never raises.
    """
    try:
        target = Path(str(path))
        sdir = _session_dir(session_id)
        sdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(sdir, 0o700)
        except OSError:
            pass
        return _write_record(
            session_id, sdir, target,
            tool=tool, pre_bytes=bytes(pre_bytes or b""), post_hash="",
            created=False,
            extra={"raw": True, "deleted": True},
        )
    except Exception:
        return None


def _write_record(
    session_id: Any,
    sdir: Path,
    target: Path,
    *,
    tool: str,
    pre_bytes: bytes,
    post_hash: str,
    created: bool,
    extra: Optional[dict] = None,
) -> dict:
    """Store one byte-exact pre-image + journal line. Shared by the
    text, binary and deletion entry points."""
    pre_hash = sha256_bytes(pre_bytes)
    truncated = len(pre_bytes) > MAX_PRE_IMAGE_BYTES
    stored = b"" if truncated else pre_bytes

    records = _read_journal(session_id)
    seq = (records[-1]["seq"] + 1) if records else 1

    pre_file = f"{seq:06d}-{pre_hash[:8]}.bin"
    _atomic_write_bytes(sdir / pre_file, stored)

    rec = {
        "seq": seq,
        "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "tool": str(tool),
        "path": str(target),
        "pre_hash": pre_hash,
        "post_hash": post_hash,
        "pre_file": pre_file,
        "created": created,
        "truncated": truncated,
    }
    rec.update(extra or {})
    dropped = _append_and_prune(session_id, sdir, records, rec)
    return {"seq": seq, "pre_hash": pre_hash, "post_hash": post_hash,
            "dropped": dropped}


def _append_and_prune(
    session_id: Any, sdir: Path, records: list[dict], rec: dict
) -> list[int]:
    """Append one record and enforce the per-session entry cap.

    The journal is rewritten atomically when pruning, the pre-image
    files of dropped entries are unlinked so the store cannot grow
    without bound, and each dropped entry is replaced by a TOMBSTONE —
    a record that keeps seq/path/tool but has no pre-image. Forgetting
    the entry outright made a later undo of it indistinguishable from
    "nothing to undo": the user asked for a file back and was told, in
    the same three empty lists, that nothing had happened.

    Returns the seqs whose pre-image this append pushed out (empty in
    the ordinary case) so the caller can say so.
    """
    journal = _journal_path(session_id)
    with journal.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    _set_file_perms(journal)

    records.append(rec)
    live = [r for r in records if not r.get("dropped")]
    if len(live) <= MAX_ENTRIES_PER_SESSION:
        return []

    drop = live[: len(live) - MAX_ENTRIES_PER_SESSION]
    dropped_seqs = {int(r["seq"]) for r in drop}
    rewritten: list[dict] = []
    for r in records:
        if int(r.get("seq", -1)) in dropped_seqs:
            rewritten.append({
                "seq": r.get("seq"),
                "ts": r.get("ts", ""),
                "tool": r.get("tool", ""),
                "path": r.get("path", ""),
                "created": bool(r.get("created")),
                "deleted": bool(r.get("deleted")),
                "dropped": True,
                "reason": (
                    f"pre-image dropped at the {MAX_ENTRIES_PER_SESSION}-entry "
                    "session cap"),
            })
        else:
            rewritten.append(r)
    # Tombstones are tiny but must not grow without bound either.
    tomb = [r for r in rewritten if r.get("dropped")]
    if len(tomb) > MAX_ENTRIES_PER_SESSION:
        stale = {id(r) for r in tomb[: len(tomb) - MAX_ENTRIES_PER_SESSION]}
        rewritten = [r for r in rewritten if id(r) not in stale]
    _atomic_write(journal, "".join(
        json.dumps(r, ensure_ascii=False) + "\n" for r in rewritten))
    for old in drop:
        try:
            (sdir / str(old.get("pre_file", ""))).unlink()
        except OSError:
            pass
    records[:] = rewritten
    return sorted(dropped_seqs)


def record_binary_change(
    session_id: Any,
    *,
    tool: str,
    path: Any,
    pre_bytes: Optional[bytes],
) -> Optional[dict]:
    """Capture a change to a file that is not text — spreadsheets, PDFs.

    The text entry point cannot serve these: it takes DECODED text and
    puts it back in the file's text convention, which a spreadsheet or
    PDF does not have. Records written here carry ``binary: True``, hold
    the raw pre-image, and are restored byte-for-byte by :func:`revert`.

    ``pre_bytes is None`` means the file did not exist before the write.
    The post-image hash is taken from the file as it is on disk now, so
    this must be called after the write it records. Like
    :func:`record_change` it never raises.
    """
    try:
        target = Path(str(path))
        sdir = _session_dir(session_id)
        sdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(sdir, 0o700)
        except OSError:
            pass

        created = pre_bytes is None
        return _write_record(
            session_id, sdir, target,
            tool=tool,
            pre_bytes=b"" if created else bytes(pre_bytes or b""),
            post_hash=sha256_file_bytes(target) or "",
            created=created,
            extra={"binary": True},
        )
    except Exception:
        return None


def list_changes(session_id: Any, *, last_n: Optional[int] = None) -> list[dict]:
    """Journal records for a session, oldest first / newest last.

    Torn lines are skipped. ``last_n`` limits to the most recent N."""
    try:
        records = _read_journal(session_id)
    except Exception:
        return []
    if last_n is not None and last_n >= 0:
        records = records[-last_n:] if last_n else []
    return records


def revert(
    session_id: Any,
    *,
    scope: str = "last",
    turn_seqs: Optional[list[int]] = None,
    workspace: Optional[Path] = None,
) -> dict:
    """Undo recorded agent changes, newest first, with conflict safety.

    scope:
        ``"last"``    — the newest entry the agent itself wrote and that
                        has not been undone yet, so repeated calls walk
                        back through the session's changes.
        ``"turn"``    — the caller-supplied ``turn_seqs`` list (the caller
                        tracks turn boundaries; seqs outside the journal
                        are ignored). This is the scope that can also
                        select an ``undo`` record, i.e. undo an undo.
        ``"session"`` — every not-yet-undone agent entry, newest first,
                        so chained edits to one file unwind step by step.

    Safety contract, per entry:

    * The file's CURRENT sha256 must equal the recorded ``post_hash``;
      any mismatch (user/other-process change since the agent's edit)
      is a conflict and the file is NEVER overwritten.
    * ``created`` entries are deleted on revert, ``deleted`` entries are
      recreated — again only when the current state still matches.
    * ``truncated`` (oversize) and ``dropped`` (pre-image pruned at the
      session cap) entries are refused (skipped), with the reason.
    * An entry already undone is skipped saying so, rather than reported
      as a conflict caused by the user.
    * When ``workspace`` is given, only paths inside it are acted on.

    Every action is journalled back (``tool="undo"``) so the undo is
    itself undoable and appears in the change listing.

    Returns ``{"reverted": [paths], "conflicts": [{"path", "reason"}],
    "skipped": [{"path", "reason"}], "undo_seqs": [int]}``.
    """
    result: dict = {"reverted": [], "conflicts": [], "skipped": [],
                    "undo_seqs": []}
    try:
        records = _read_journal(session_id)
    except Exception:
        return result

    def _agent_records() -> list[dict]:
        """Entries the agent wrote and that are still standing.

        Undo records are excluded: "undo the last change" must walk
        back through the agent's edits, not toggle the newest one on
        and off. They stay selectable through ``turn``.
        """
        return [r for r in records
                if r.get("tool") != "undo" and not r.get("undone")]

    if scope == "last":
        selected = _agent_records()[-1:]
    elif scope == "turn":
        wanted = {int(s) for s in (turn_seqs or [])}
        selected = [r for r in records if r["seq"] in wanted]
    elif scope == "session":
        selected = _agent_records()
    else:
        result["skipped"].append(
            {"path": "", "reason": f"unknown scope {scope!r}"})
        return result

    if not selected and records:
        # Three empty lists read exactly like "there was nothing to undo".
        # Say which of the two it is: everything already undone, or a
        # turn/seq selection that matched no record.
        undone = [r for r in records if r.get("undone")]
        if scope in ("last", "session") and undone:
            result["skipped"].append({
                "path": str(undone[-1].get("path", "")),
                "reason": ("every recorded change in this session has "
                           "already been undone"),
            })
        elif scope == "turn":
            result["skipped"].append({
                "path": "",
                "reason": ("none of the given change ids are in this "
                           "session's journal"),
            })
        return result

    ws_resolved: Optional[Path] = None
    if workspace is not None:
        try:
            ws_resolved = Path(workspace).expanduser().resolve()
        except OSError:
            ws_resolved = Path(workspace)

    sdir = _session_dir(session_id)
    undone_now: dict[int, int] = {}      # reverted seq -> undo record seq

    def _journal_undo(rec: dict, target: Path, *, pre_bytes: Optional[bytes],
                      deleted: bool) -> None:
        """Record what the undo itself did, so it can be undone."""
        try:
            src_seq = int(rec.get("seq", 0))
            if deleted:
                out = record_deletion(
                    session_id, tool="undo", path=target, pre_bytes=pre_bytes)
            else:
                out = _write_record(
                    session_id, sdir, target,
                    tool="undo",
                    pre_bytes=bytes(pre_bytes or b""),
                    post_hash=sha256_file_bytes(target) or "",
                    created=pre_bytes is None,
                    extra={"raw": True, "undo_of": src_seq},
                )
            if out:
                undone_now[src_seq] = int(out.get("seq", 0))
                result["undo_seqs"].append(int(out.get("seq", 0)))
        except Exception:
            pass

    # Newest first: for chained edits to one file the latest post_hash is
    # what is on disk NOW; each restore re-establishes the previous
    # entry's post state, so the chain unwinds cleanly.
    for rec in sorted(selected, key=lambda r: r["seq"], reverse=True):
        path_str = str(rec.get("path", ""))
        try:
            target = Path(path_str)

            if ws_resolved is not None:
                try:
                    inside = target.resolve().is_relative_to(ws_resolved)
                except OSError:
                    inside = False
                if not inside:
                    result["skipped"].append(
                        {"path": path_str, "reason": "outside workspace"})
                    continue

            if rec.get("undone"):
                _u = rec.get("undone") or {}
                result["skipped"].append({
                    "path": path_str,
                    "reason": (
                        "already undone at "
                        f"{_u.get('ts', '?')} (undo record "
                        f"#{_u.get('by_seq', '?')}) — nothing left to restore"),
                })
                continue

            if rec.get("dropped"):
                result["skipped"].append({
                    "path": path_str,
                    "reason": (
                        "the pre-image was dropped at the "
                        f"{MAX_ENTRIES_PER_SESSION}-entry session cap; this "
                        "change is recorded but can no longer be restored"),
                })
                continue

            if rec.get("truncated"):
                result["skipped"].append({
                    "path": path_str,
                    "reason": ("pre-image was stored truncated (>2 MB); "
                               "revert refused to avoid corrupting the file"),
                })
                continue

            if rec.get("lossy"):
                result["skipped"].append({
                    "path": path_str,
                    "reason": ("the pre-image could not be stored in the "
                               "file's own encoding, so restoring it would "
                               "change bytes the agent never touched; "
                               "revert refused"),
                })
                continue

            byte_exact = bool(rec.get("binary") or rec.get("raw"))
            exists = target.exists()
            if not exists:
                current = None
            elif byte_exact:
                current = sha256_file_bytes(target)
            else:
                current = sha256_file(target)

            # A deletion's post state is the file's ABSENCE.
            if rec.get("deleted"):
                if exists:
                    result["conflicts"].append({
                        "path": path_str,
                        "reason": ("the file exists again since the agent "
                                   "deleted it — not overwriting"),
                    })
                    continue
                pre_data = _read_pre_image(sdir, rec, result, path_str)
                if pre_data is None:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                _atomic_restore_bytes(target, pre_data)
                result["reverted"].append(path_str)
                # Inverse of "recreate" is "create": undoing it removes it.
                _journal_undo(rec, target, pre_bytes=None, deleted=False)
                continue

            if not exists:
                if rec.get("created"):
                    result["skipped"].append(
                        {"path": path_str,
                         "reason": "created file already absent"})
                else:
                    result["conflicts"].append(
                        {"path": path_str,
                         "reason": "file missing since the agent's edit"})
                continue

            if current != rec.get("post_hash"):
                result["conflicts"].append({
                    "path": path_str,
                    "reason": ("content changed since the agent's edit "
                               "(post-hash mismatch) — not overwriting"),
                })
                continue

            if rec.get("created"):
                agent_bytes = target.read_bytes()
                target.unlink()
                result["reverted"].append(path_str)
                # Inverse of "delete" is "recreate what was deleted".
                _journal_undo(rec, target, pre_bytes=agent_bytes, deleted=True)
                continue

            pre_data = _read_pre_image(sdir, rec, result, path_str)
            if pre_data is None:
                continue

            if byte_exact:
                agent_bytes = target.read_bytes()
                _atomic_restore_bytes(target, pre_data)
                result["reverted"].append(path_str)
                _journal_undo(rec, target, pre_bytes=agent_bytes, deleted=False)
                continue

            # Legacy text record (written before byte-exact pre-images).
            try:
                pre_text = pre_data.decode("utf-8")
            except UnicodeDecodeError:
                result["conflicts"].append(
                    {"path": path_str,
                     "reason": "pre-image is not decodable text"})
                continue
            if sha256_text(pre_text) != rec.get("pre_hash"):
                result["conflicts"].append(
                    {"path": path_str, "reason": "pre-image corrupt "
                     "(stored hash mismatch)"})
                continue
            agent_bytes = target.read_bytes()
            _atomic_restore(target, pre_text)
            result["reverted"].append(path_str)
            _journal_undo(rec, target, pre_bytes=agent_bytes, deleted=False)
        except Exception as exc:
            result["conflicts"].append(
                {"path": path_str, "reason": f"revert failed: {exc}"})

    if undone_now:
        _mark_undone(session_id, undone_now)
    return result


def _read_pre_image(
    sdir: Path, rec: dict, result: dict, path_str: str,
) -> Optional[bytes]:
    """Load and verify a record's stored pre-image.

    Returns ``None`` after filing the conflict, so the caller can just
    ``continue``. Verifying the stored hash first is what keeps a
    half-written or hand-edited store from being written over a user
    file."""
    pre_path = sdir / str(rec.get("pre_file", ""))
    try:
        data = pre_path.read_bytes()
    except OSError:
        result["conflicts"].append(
            {"path": path_str, "reason": "pre-image file missing"})
        return None
    if bool(rec.get("binary") or rec.get("raw")):
        if sha256_bytes(data) != rec.get("pre_hash"):
            result["conflicts"].append(
                {"path": path_str,
                 "reason": "pre-image corrupt (stored hash mismatch)"})
            return None
    return data


def _mark_undone(session_id: Any, undone: dict[int, int]) -> None:
    """Flag the reverted records so a second undo says "already undone".

    Without this the second run recomputes the hash of a file the UNDO
    wrote, finds it different from the agent's post-image, and reports
    "content changed since the agent's edit" — blaming the user for the
    agent's own undo. Never raises: the file has already been restored,
    and losing the bookkeeping must not turn that into an error.
    """
    try:
        records = _read_journal(session_id)
        stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
        changed = False
        for rec in records:
            by = undone.get(int(rec.get("seq", -1)))
            if by is not None and not rec.get("undone"):
                rec["undone"] = {"ts": stamp, "by_seq": by}
                changed = True
        if changed:
            _atomic_write(_journal_path(session_id), "".join(
                json.dumps(r, ensure_ascii=False) + "\n" for r in records))
    except Exception:
        pass
