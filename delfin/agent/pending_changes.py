"""Pending-edit store for the diff-approval permission mode.

In ``diff_approval`` mode the agent's file mutations are NOT applied
immediately: each write is captured here as a pending unified diff that
the user reviews in the dashboard. Only an explicit approval applies
the staged content to the workspace; a rejection discards it. This
complements plan mode (which blocks writes entirely) — work proceeds,
but every applied change stays under user control.

Storage layout (per session, under ``~/.delfin/pending/<sid>/``)::

    pending.jsonl        one JSON record per staged change
    <seq:06d>-old.txt    full pre-image (absent for created files)
    <seq:06d>-new.txt    full post-image (the content approval applies)

Record schema::

    {seq, id, created_at, tool, path, note,
     status,              # "pending" | "applied" | "rejected"
     created: bool,       # True → the file did not exist when staged
     old_hash, new_hash,  # sha256 of pre-/post-image
     old_file, new_file,  # image file names (old_file "" for created)
     diff,                # unified diff, capped for display
     applied_at?, rejected_at?, undo_seq?, last_conflict?}

Guarantees (mirroring ``change_journal``):

* No public function raises — a broken store must not take down the
  tool path or the dashboard. ``stage`` reports failure as
  ``{"error": ...}``; the others return safe defaults.
* Journal rewrites and image files are atomic (temp + ``os.replace``)
  and 0600; the per-session directory is 0700. Session ids are
  sanitized before use as a directory name (no traversal).
* ``approve`` only writes when the target file's CURRENT content still
  hashes to the staged pre-image (created files: the target must not
  exist). Anything else is a conflict — reported with a reason, the
  file is left untouched, and the change stays pending.
* Approved changes are recorded in the undo journal
  (:func:`change_journal.record_change`), so session-level undo covers
  them exactly like directly applied edits.
* One pending change per file: staging a second change for a path that
  already has one pending is refused. The disk file does not move while
  a change is pending, so a second stage would be based on the
  un-applied state and conflict at approval time; refusing at stage
  time surfaces the problem when the model can still react.
* ``stage`` refuses when the pending queue is full — silently dropping
  changes that await a user decision would defeat the mode. Terminal
  records (applied/rejected) are pruned beyond a separate cap.
"""

from __future__ import annotations

import difflib
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Optional

from .change_journal import sha256_file, sha256_text

MAX_IMAGE_BYTES = 2 * 1024 * 1024     # per image; larger stages are refused
MAX_PENDING_PER_SESSION = 200         # un-reviewed queue cap (stage refuses)
MAX_TERMINAL_RECORDS = 300            # applied/rejected history cap (pruned)
MAX_DIFF_CHARS = 4_000                # diff text stored inside the record

_JOURNAL_NAME = "pending.jsonl"


# ---------------------------------------------------------------------------
# Helpers (conventions shared with change_journal)
# ---------------------------------------------------------------------------

def _safe_session_id(session_id: Any) -> str:
    """Sanitize a session id for use as a directory name (no traversal)."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", str(session_id or "") or "session")[:40]


def _pending_root() -> Path:
    return Path.home() / ".delfin" / "pending"


def _session_dir(session_id: Any) -> Path:
    return _pending_root() / _safe_session_id(session_id)


def _journal_path(session_id: Any) -> Path:
    return _session_dir(session_id) / _JOURNAL_NAME


def _set_file_perms(path: Path) -> None:
    """Best-effort 0600 on a journal/image file (per-user data)."""
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def _atomic_write(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` atomically (temp file + ``os.replace``).

    ``newline=""`` disables newline translation so images round-trip
    byte-exact (same pattern as ``change_journal._atomic_write``)."""
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


def _atomic_write_workspace(
    path: Path, text: str, *, default_mode: int = 0o644,
) -> None:
    """Atomically write a WORKSPACE file, preserving its current permission
    bits (an approved user file must not silently become 0600 like the
    internal store files). New files get ``default_mode`` best-effort."""
    import tempfile
    mode: Optional[int] = None
    try:
        mode = os.stat(path).st_mode & 0o7777
    except OSError:
        pass
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".approve.tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        try:
            os.chmod(tmp, mode if mode is not None else default_mode)
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


def _rewrite_journal(session_id: Any, records: list[dict]) -> None:
    _atomic_write(_journal_path(session_id), "".join(
        json.dumps(r, ensure_ascii=False) + "\n" for r in records))


def _unlink_images(session_id: Any, rec: dict) -> None:
    """Best-effort removal of a record's image files (terminal states)."""
    sdir = _session_dir(session_id)
    for key in ("old_file", "new_file"):
        name = str(rec.get(key, "") or "")
        if not name:
            continue
        try:
            (sdir / name).unlink()
        except OSError:
            pass


def _unified_diff(old: str, new: str, label: str, max_chars: int) -> str:
    """Unified diff in the same shape as the chat write-path diffs
    (``a/<label>`` → ``b/<label>``, 3 context lines), capped by chars."""
    diff = list(difflib.unified_diff(
        old.splitlines(keepends=False),
        new.splitlines(keepends=False),
        fromfile=f"a/{label}",
        tofile=f"b/{label}",
        lineterm="",
        n=3,
    ))
    if not diff:
        return "(no changes)"
    text = "\n".join(diff)
    if len(text) > max_chars:
        text = text[:max_chars] + "\n... (diff truncated)"
    return text


def _display_label(path_str: str) -> str:
    """Short label for diffs/render: home-relative when possible."""
    try:
        home = str(Path.home())
        if path_str.startswith(home + os.sep):
            return "~" + path_str[len(home):]
    except Exception:
        pass
    return path_str


def _match_id(rec: dict, change_id: Any) -> bool:
    """True when ``change_id`` names ``rec`` — accepts "3", 3, "000003"."""
    wanted = str(change_id).strip().lstrip("0") or "0"
    return str(rec.get("seq")) == wanted or str(rec.get("id")) == str(change_id).strip()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def stage(
    session_id: Any,
    *,
    tool: str,
    path: Any,
    old_text: Optional[str],
    new_text: str,
    note: str = "",
) -> dict:
    """Capture one file change as pending — BEFORE any application.

    ``old_text is None`` means the file does not exist yet (approval
    creates it). Returns the full record on success; on ANY failure or
    refusal returns ``{"error": reason}`` and never raises — the write
    path that calls this must not break on store trouble.
    """
    try:
        if not isinstance(new_text, str):
            return {"error": "new_text must be a string"}
        if old_text is not None and not isinstance(old_text, str):
            return {"error": "old_text must be a string or None"}
        if old_text == new_text:
            return {"error": "staged change is empty (old and new content are identical)"}
        path_str = str(path)
        if not path_str:
            return {"error": "path is required"}
        for what, text in (("old_text", old_text), ("new_text", new_text)):
            if text is not None and len(text.encode("utf-8", errors="replace")) > MAX_IMAGE_BYTES:
                return {"error": (
                    f"{what} exceeds {MAX_IMAGE_BYTES // (1024 * 1024)} MB — "
                    "too large to stage for approval; apply this change "
                    "outside diff-approval mode or split it up."
                )}

        records = _read_journal(session_id)
        pending = [r for r in records if r.get("status") == "pending"]
        if len(pending) >= MAX_PENDING_PER_SESSION:
            return {"error": (
                f"pending queue is full ({MAX_PENDING_PER_SESSION} staged "
                "changes) — the user must approve or reject some before "
                "more can be staged."
            )}
        for r in pending:
            if str(r.get("path", "")) == path_str:
                return {"error": (
                    f"a change for '{path_str}' is already pending "
                    f"(id {r.get('id')}) — the user must approve or reject "
                    "it before this file can be staged again."
                )}

        sdir = _session_dir(session_id)
        sdir.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(sdir, 0o700)
        except OSError:
            pass

        created = old_text is None
        pre_text = "" if created else str(old_text)
        seq = (records[-1]["seq"] + 1) if records else 1

        old_file = "" if created else f"{seq:06d}-old.txt"
        new_file = f"{seq:06d}-new.txt"
        if old_file:
            _atomic_write(sdir / old_file, pre_text)
        _atomic_write(sdir / new_file, new_text)

        rec = {
            "seq": seq,
            "id": str(seq),
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "tool": str(tool),
            "path": path_str,
            "note": str(note or "")[:500],
            "status": "pending",
            "created": created,
            "old_hash": sha256_text(pre_text),
            "new_hash": sha256_text(new_text),
            "old_file": old_file,
            "new_file": new_file,
            "diff": _unified_diff(
                pre_text, new_text, _display_label(path_str), MAX_DIFF_CHARS),
        }
        journal = _journal_path(session_id)
        with journal.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        _set_file_perms(journal)

        # History cap: prune the OLDEST terminal (applied/rejected)
        # records only — pending entries always survive.
        records.append(rec)
        terminal = [r for r in records if r.get("status") != "pending"]
        if len(terminal) > MAX_TERMINAL_RECORDS:
            drop_seqs = {r["seq"] for r in terminal[: len(terminal) - MAX_TERMINAL_RECORDS]}
            keep = [r for r in records if r["seq"] not in drop_seqs]
            _rewrite_journal(session_id, keep)
            for r in terminal:
                if r["seq"] in drop_seqs:
                    _unlink_images(session_id, r)

        return dict(rec)
    except Exception as exc:
        return {"error": f"stage failed: {exc}"}


def list_pending(session_id: Any) -> list[dict]:
    """Pending (un-reviewed) records for a session, oldest first."""
    try:
        return [dict(r) for r in _read_journal(session_id)
                if r.get("status") == "pending"]
    except Exception:
        return []


def get(session_id: Any, change_id: Any) -> Optional[dict]:
    """One record (any status) by id, or None. Never raises."""
    try:
        for rec in _read_journal(session_id):
            if _match_id(rec, change_id):
                return dict(rec)
    except Exception:
        pass
    return None


def approve(session_id: Any, change_id: Any, *, workspace: Any) -> dict:
    """Apply one pending change to the workspace, with conflict safety.

    Safety contract:

    * The target must resolve inside ``workspace``.
    * Modified files: the CURRENT sha256 must equal the staged
      ``old_hash`` — any mismatch (user/other-process edit since
      staging, or the file vanished) is a conflict; the file is never
      touched and the change stays pending with ``last_conflict`` set.
    * Created files: the target must NOT exist yet.
    * The write is atomic and preserves the file's permission bits.
    * On success the change is recorded in the undo journal
      (``change_journal.record_change``) so session undo covers it.

    Returns ``{"status": "applied"|"conflict"|"not_found"|
    "already_applied"|"already_rejected"|"error", ...}``. Never raises.
    """
    try:
        records = _read_journal(session_id)
        rec = next((r for r in records if _match_id(r, change_id)), None)
        if rec is None:
            return {"status": "not_found", "id": str(change_id)}
        if rec.get("status") == "applied":
            return {"status": "already_applied", "id": rec["id"],
                    "path": rec.get("path", "")}
        if rec.get("status") == "rejected":
            return {"status": "already_rejected", "id": rec["id"],
                    "path": rec.get("path", "")}

        def _conflict(reason: str) -> dict:
            rec["last_conflict"] = reason
            try:
                _rewrite_journal(session_id, records)
            except Exception:
                pass
            return {"status": "conflict", "id": rec["id"],
                    "path": rec.get("path", ""), "reason": reason}

        target = Path(str(rec.get("path", "")))
        try:
            ws = Path(workspace).expanduser().resolve()
        except OSError:
            ws = Path(workspace)
        try:
            inside = target.expanduser().resolve().is_relative_to(ws)
        except OSError:
            inside = False
        if not inside:
            return _conflict("target is outside the workspace")

        sdir = _session_dir(session_id)
        try:
            new_text = (sdir / str(rec.get("new_file", ""))).read_bytes().decode("utf-8")
        except OSError:
            return _conflict("staged post-image file is missing")
        if sha256_text(new_text) != rec.get("new_hash"):
            return _conflict("staged post-image is corrupt (hash mismatch)")

        if rec.get("created"):
            if target.exists():
                return _conflict(
                    "file appeared on disk since staging — not overwriting")
            target.parent.mkdir(parents=True, exist_ok=True)
            _atomic_write_workspace(target, new_text)
            old_for_undo: Optional[str] = None
        else:
            if not target.exists():
                return _conflict("file is missing since staging")
            if sha256_file(target) != rec.get("old_hash"):
                return _conflict(
                    "content changed since staging (hash mismatch) — "
                    "not overwriting; re-stage from the current content")
            try:
                old_for_undo = (
                    sdir / str(rec.get("old_file", ""))
                ).read_bytes().decode("utf-8")
            except OSError:
                return _conflict("staged pre-image file is missing")
            if sha256_text(old_for_undo) != rec.get("old_hash"):
                return _conflict("staged pre-image is corrupt (hash mismatch)")
            _atomic_write_workspace(target, new_text)

        # Undo-journal hand-off: an approved change must be revertible
        # exactly like a directly applied edit. Done here (not in the
        # wiring) so EVERY approval surface gets undo coverage.
        undo_seq: Optional[int] = None
        try:
            from . import change_journal as _cj
            undo = _cj.record_change(
                session_id, tool=str(rec.get("tool", "")),
                path=str(target), old_text=old_for_undo, new_text=new_text)
            if undo is not None:
                undo_seq = undo.get("seq")
        except Exception:
            pass

        rec["status"] = "applied"
        rec["applied_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        rec.pop("last_conflict", None)
        if undo_seq is not None:
            rec["undo_seq"] = undo_seq
        _rewrite_journal(session_id, records)
        _unlink_images(session_id, rec)

        return {"status": "applied", "id": rec["id"],
                "path": str(target), "undo_seq": undo_seq}
    except Exception as exc:
        return {"status": "error", "id": str(change_id),
                "reason": f"approve failed: {exc}"}


def reject(session_id: Any, change_id: Any) -> dict:
    """Discard one pending change (images removed, record kept for audit).

    Returns ``{"status": "rejected"|"not_found"|"already_applied"|
    "already_rejected"|"error", ...}``. Never raises."""
    try:
        records = _read_journal(session_id)
        rec = next((r for r in records if _match_id(r, change_id)), None)
        if rec is None:
            return {"status": "not_found", "id": str(change_id)}
        if rec.get("status") == "applied":
            return {"status": "already_applied", "id": rec["id"],
                    "path": rec.get("path", "")}
        if rec.get("status") == "rejected":
            return {"status": "already_rejected", "id": rec["id"],
                    "path": rec.get("path", "")}
        rec["status"] = "rejected"
        rec["rejected_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        rec.pop("last_conflict", None)
        _rewrite_journal(session_id, records)
        _unlink_images(session_id, rec)
        return {"status": "rejected", "id": rec["id"],
                "path": rec.get("path", "")}
    except Exception as exc:
        return {"status": "error", "id": str(change_id),
                "reason": f"reject failed: {exc}"}


def approve_all(session_id: Any, *, workspace: Any) -> dict:
    """Approve every pending change, oldest first.

    Returns ``{"applied": [...], "conflicts": [...], "errors": [...]}``
    with one entry per change (the per-change ``approve`` result).
    Never raises."""
    out: dict = {"applied": [], "conflicts": [], "errors": []}
    try:
        for rec in list_pending(session_id):
            res = approve(session_id, rec["id"], workspace=workspace)
            status = res.get("status", "error")
            if status == "applied":
                out["applied"].append(res)
            elif status == "conflict":
                out["conflicts"].append(res)
            else:
                out["errors"].append(res)
    except Exception as exc:
        out["errors"].append({"status": "error", "reason": str(exc)})
    return out


def reject_all(session_id: Any) -> dict:
    """Reject every pending change. Returns ``{"rejected": [ids],
    "errors": [...]}``. Never raises."""
    out: dict = {"rejected": [], "errors": []}
    try:
        for rec in list_pending(session_id):
            res = reject(session_id, rec["id"])
            if res.get("status") == "rejected":
                out["rejected"].append(res.get("id"))
            else:
                out["errors"].append(res)
    except Exception as exc:
        out["errors"].append({"status": "error", "reason": str(exc)})
    return out


def render_pending(
    session_id: Any,
    *,
    max_items: int = 10,
    max_diff_lines: int = 12,
) -> str:
    """Compact markdown summary of the pending queue for the dashboard.

    Never raises; returns a short notice when the queue is empty."""
    try:
        pending = list_pending(session_id)
    except Exception:
        pending = []
    if not pending:
        return "No pending changes — nothing awaits approval."
    lines: list[str] = [
        f"**Pending changes — {len(pending)} awaiting approval**", ""]
    for rec in pending[:max_items]:
        kind = "new file" if rec.get("created") else "edit"
        lines.append(
            f"**#{rec.get('id')}** · {rec.get('tool', '?')} ({kind}) · "
            f"`{_display_label(str(rec.get('path', '')))}` · "
            f"{rec.get('created_at', '')}"
        )
        if rec.get("note"):
            lines.append(f"  note: {rec['note']}")
        if rec.get("last_conflict"):
            lines.append(f"  last conflict: {rec['last_conflict']}")
        diff_lines = str(rec.get("diff", "")).splitlines()
        shown = diff_lines[:max_diff_lines]
        lines.append("```diff")
        lines.extend(shown)
        if len(diff_lines) > max_diff_lines:
            lines.append(f"... ({len(diff_lines) - max_diff_lines} more diff lines)")
        lines.append("```")
        lines.append("")
    if len(pending) > max_items:
        lines.append(f"... and {len(pending) - max_items} more (use /pending).")
    lines.append(
        "Approve: `/approve <id>` or `/approve all` · "
        "Reject: `/reject <id>` or `/reject all`")
    return "\n".join(lines)


__all__ = [
    "stage", "list_pending", "get", "approve", "reject",
    "approve_all", "reject_all", "render_pending",
    "MAX_IMAGE_BYTES", "MAX_PENDING_PER_SESSION",
    "MAX_TERMINAL_RECORDS", "MAX_DIFF_CHARS",
]
