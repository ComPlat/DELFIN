"""Undo journal for agent file edits — pre-image capture + safe revert.

The agent's write path (write_file / edit_file / apply_patch /
notebook_edit) reads each file's old content before overwriting it (for
the chat diff) and then throws it away. Calc workspaces frequently are
NOT git repos, so a single agent edit clobbers uncommitted user work
irrecoverably. This module keeps that pre-image:

Storage layout (per session, under ``~/.delfin/undo/<sid>/``)::

    journal.jsonl               one JSON record per captured change
    <seq:06d>-<prehash8>.txt    full pre-image (content-addressed)

Journal record schema::

    {seq, ts, tool, path, pre_hash, post_hash, pre_file,
     created: bool,      # True → the file did not exist before the write
     truncated: bool}    # True → pre-image exceeded MAX_PRE_IMAGE_BYTES

Guarantees:

* ``record_change`` never raises — a broken undo store must not block the
  tool path it observes (returns ``None`` on any failure).
* Pre-images and journal rewrites are atomic (temp file + ``os.replace``)
  and 0600; the per-session directory is 0700.
* Pre-images larger than ``MAX_PRE_IMAGE_BYTES`` (2 MB) are stored
  TRUNCATED and flagged; :func:`revert` refuses those entries (a partial
  restore would silently corrupt the file — losing undo for a huge file
  is better than "undoing" it into garbage).
* Per-session cap of ``MAX_ENTRIES_PER_SESSION`` (500) entries; oldest
  entries are pruned and their pre-image files unlinked.
* ``revert`` only ever touches a file whose CURRENT sha256 equals the
  recorded ``post_hash``. Anything else — the user or another process
  changed the file since the agent's edit — is reported as a conflict
  and the file is left untouched.
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

def record_change(
    session_id: Any,
    *,
    tool: str,
    path: Any,
    old_text: Optional[str],
    new_text: Optional[str],
) -> Optional[dict]:
    """Capture one file change: store the full pre-image + a journal line.

    ``old_text is None`` means the file did not exist before the write
    (``created=True``; revert deletes it). Returns
    ``{"seq": int, "pre_hash": str, "post_hash": str}`` on success and
    ``None`` on ANY failure — the journal must never break the write path
    it observes.
    """
    try:
        sdir = _session_dir(session_id)
        sdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(sdir, 0o700)
        except OSError:
            pass

        created = old_text is None
        pre_text = "" if created else str(old_text)
        post_text = "" if new_text is None else str(new_text)
        pre_hash = sha256_text(pre_text)
        post_hash = sha256_text(post_text)

        # Size guard: store a truncated pre-image, flagged so revert
        # refuses it (restoring a truncation would corrupt the file).
        raw = pre_text.encode("utf-8", errors="replace")
        truncated = len(raw) > MAX_PRE_IMAGE_BYTES
        stored_text = (
            raw[:MAX_PRE_IMAGE_BYTES].decode("utf-8", errors="ignore")
            if truncated else pre_text
        )

        records = _read_journal(session_id)
        seq = (records[-1]["seq"] + 1) if records else 1

        pre_file = f"{seq:06d}-{pre_hash[:8]}.txt"
        _atomic_write(sdir / pre_file, stored_text)

        rec = {
            "seq": seq,
            "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "tool": str(tool),
            "path": str(path),
            "pre_hash": pre_hash,
            "post_hash": post_hash,
            "pre_file": pre_file,
            "created": created,
            "truncated": truncated,
        }
        journal = _journal_path(session_id)
        line = json.dumps(rec, ensure_ascii=False) + "\n"
        with journal.open("a", encoding="utf-8") as fh:
            fh.write(line)
        _set_file_perms(journal)

        # Per-session cap: prune oldest entries (journal rewritten
        # atomically, orphaned pre-image files unlinked).
        records.append(rec)
        if len(records) > MAX_ENTRIES_PER_SESSION:
            drop = records[: len(records) - MAX_ENTRIES_PER_SESSION]
            keep = records[len(records) - MAX_ENTRIES_PER_SESSION:]
            _atomic_write(journal, "".join(
                json.dumps(r, ensure_ascii=False) + "\n" for r in keep))
            for old in drop:
                try:
                    (sdir / str(old.get("pre_file", ""))).unlink()
                except OSError:
                    pass

        return {"seq": seq, "pre_hash": pre_hash, "post_hash": post_hash}
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
        ``"last"``    — only the most recent journal entry.
        ``"turn"``    — the caller-supplied ``turn_seqs`` list (the caller
                        tracks turn boundaries; seqs outside the journal
                        are ignored).
        ``"session"`` — every entry in the session, newest first, so
                        chained edits to one file unwind step by step.

    Safety contract, per entry:

    * The file's CURRENT sha256 must equal the recorded ``post_hash``;
      any mismatch (user/other-process change since the agent's edit)
      is a conflict and the file is NEVER overwritten.
    * ``created`` entries are deleted on revert — again only when the
      current hash still matches.
    * ``truncated`` (oversize) entries are refused (skipped).
    * When ``workspace`` is given, only paths inside it are acted on.

    Returns ``{"reverted": [paths], "conflicts": [{"path", "reason"}],
    "skipped": [{"path", "reason"}]}``.
    """
    result: dict = {"reverted": [], "conflicts": [], "skipped": []}
    try:
        records = _read_journal(session_id)
    except Exception:
        return result

    if scope == "last":
        selected = records[-1:]
    elif scope == "turn":
        wanted = {int(s) for s in (turn_seqs or [])}
        selected = [r for r in records if r["seq"] in wanted]
    elif scope == "session":
        selected = list(records)
    else:
        result["skipped"].append(
            {"path": "", "reason": f"unknown scope {scope!r}"})
        return result

    ws_resolved: Optional[Path] = None
    if workspace is not None:
        try:
            ws_resolved = Path(workspace).expanduser().resolve()
        except OSError:
            ws_resolved = Path(workspace)

    sdir = _session_dir(session_id)

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

            if rec.get("truncated"):
                result["skipped"].append({
                    "path": path_str,
                    "reason": ("pre-image was stored truncated (>2 MB); "
                               "revert refused to avoid corrupting the file"),
                })
                continue

            exists = target.exists()
            current = sha256_file(target) if exists else None

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
                target.unlink()
                result["reverted"].append(path_str)
                continue

            pre_path = sdir / str(rec.get("pre_file", ""))
            try:
                # Byte-exact read (no universal-newline translation).
                pre_text = pre_path.read_bytes().decode("utf-8")
            except OSError:
                result["conflicts"].append(
                    {"path": path_str, "reason": "pre-image file missing"})
                continue
            if sha256_text(pre_text) != rec.get("pre_hash"):
                result["conflicts"].append(
                    {"path": path_str, "reason": "pre-image corrupt "
                     "(stored hash mismatch)"})
                continue

            _atomic_restore(target, pre_text)
            result["reverted"].append(path_str)
        except Exception as exc:
            result["conflicts"].append(
                {"path": path_str, "reason": f"revert failed: {exc}"})
    return result
