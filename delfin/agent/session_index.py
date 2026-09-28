"""Full-text index over DELFIN's own session transcripts.

A finished session leaves behind distilled facts (memories) but the episode
itself -- "how did we solve the Fe complex problem?" -- was not findable.
This module indexes the two archives DELFIN already keeps, and nothing
else:

* ``~/.delfin/agent_sessions/<id>.json`` (``session_store``: title, dates,
  ``chat_messages`` / ``engine_messages``)
* ``~/.delfin/transcript_archive/<id>.jsonl`` (pre-compaction snapshots)

There are deliberately NO directory walks over HOME or anywhere else
(cluster operating rule). The index lives in
``~/.delfin/session_index.sqlite`` (0600), is incremental (mtime + size per
source file), and every text is scrubbed with ``memory_store._without_secrets``
BEFORE it is stored -- a credential that reached the index would be
searchable, which is worse than it sitting in a 0600 file.

When the Python sqlite3 of this machine has no FTS5 module, the same
interface answers via plain ``LIKE`` search: a slower index, not a missing
one. Both paths are proven by tests.
"""

from __future__ import annotations

import json
import os
import re
import sqlite3
import threading
import time
from dataclasses import dataclass
from pathlib import Path

from . import session_store

#: Hard caps so a search can never balloon: at most five hits by default,
#: at most 20 when the caller asks for more, snippets at most 240 chars.
DEFAULT_LIMIT = 5
MAX_LIMIT = 20
MAX_SNIPPET = 240

_write_lock = threading.Lock()


def _index_path() -> Path:
    """``~/.delfin/session_index.sqlite``. A function (not a constant) so
    the state-path redirect tables and tests can point it elsewhere."""
    return Path.home() / ".delfin" / "session_index.sqlite"


@dataclass
class Hit:
    """One search result: what session matched, and a short context."""
    session_id: str
    date: str
    title: str
    snippet: str


_fts5_cache: bool | None = None


def _fts5_available() -> bool:
    """Whether this interpreter's sqlite3 has the FTS5 extension."""
    global _fts5_cache
    if _fts5_cache is None:
        try:
            con = sqlite3.connect(":memory:")
            try:
                con.execute("CREATE VIRTUAL TABLE probe USING fts5(x)")
            finally:
                con.close()
            _fts5_cache = True
        except (sqlite3.Error, OSError):
            _fts5_cache = False
    return _fts5_cache


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


def _open_index() -> sqlite3.Connection:
    """Open (and create) the index with the schema this build needs.

    The file is created 0600 from the start: session transcripts are
    user-only, and the index makes them searchable, so it inherits the
    same rule (``state_paths.secure_file``).
    """
    p = _index_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(p), timeout=30)
    con.execute("PRAGMA journal_mode=WAL")
    if _fts5_available():
        con.execute(
            "CREATE VIRTUAL TABLE IF NOT EXISTS documents USING fts5("
            "session_id UNINDEXED, created_at UNINDEXED, title UNINDEXED, "
            "text)")
    else:
        # LIKE fallback: same columns, no virtual table. Search reads this
        # table with substring matching instead of the MATCH operator.
        con.execute(
            "CREATE TABLE IF NOT EXISTS documents ("
            "session_id TEXT, created_at REAL, title TEXT, text TEXT)")
        con.execute(
            "CREATE INDEX IF NOT EXISTS documents_session "
            "ON documents(session_id)")
    con.execute(
        "CREATE TABLE IF NOT EXISTS sources ("
        "session_id TEXT NOT NULL, kind TEXT NOT NULL, "
        "mtime REAL NOT NULL, size INTEGER NOT NULL, "
        "PRIMARY KEY (session_id, kind))")
    # Be strict about permissions on every open: an index created by an
    # older run (or a copy) must not stay wider than the transcripts are.
    try:
        os.chmod(p, 0o600)
    except OSError:
        pass
    return con


# ---------------------------------------------------------------------------
# Text extraction
# ---------------------------------------------------------------------------


def _message_text(msg) -> str:
    """One message's searchable text. Never raises."""
    if not isinstance(msg, dict):
        return ""
    content = msg.get("content", "")
    if not isinstance(content, str):
        try:
            content = json.dumps(content, ensure_ascii=False)
        except (TypeError, ValueError):
            content = ""
    return f"{msg.get('role', '')}: {content}"


def _session_json_documents(session_id: str):
    """(kind, created_at, title, text) from the session json, if present."""
    p = Path(session_store._SESSIONS_DIR) / f"{session_id}.json"
    if not p.exists():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, ValueError):
        return []
    if not isinstance(data, dict):
        return []
    parts: list[str] = []
    # engine_messages can carry content as structured lists (tool calls);
    # _message_text stringifies what is not already a string.
    for key in ("chat_messages", "engine_messages"):
        msgs = data.get(key) or []
        if isinstance(msgs, list):
            parts.extend(_message_text(m) for m in msgs)
    title = str(data.get("title") or "Untitled")
    created = float(data.get("created_at") or data.get("updated_at") or 0.0)
    stat = p.stat()
    return [(f"session:{stat.st_mtime}:{stat.st_size}",
             created, title, "\n".join(parts))]


def _archive_documents(session_id: str):
    """(kind, created_at, title, text) from the transcript archive."""
    p = session_store._transcript_archive_path() / f"{session_id}.jsonl"
    if not p.exists():
        return []
    try:
        stat = p.stat()
        lines = p.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    parts: list[str] = []
    first_created = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(rec, dict):
            continue
        if first_created is None:
            first_created = rec.get("compacted_at")
        msgs = rec.get("messages") or []
        if isinstance(msgs, list):
            parts.extend(_message_text(m) for m in msgs)
    return [(f"archive:{stat.st_mtime}:{stat.st_size}",
             float(first_created or 0.0), "Untitled",
             "\n".join(parts))]


# ---------------------------------------------------------------------------
# Indexing
# ---------------------------------------------------------------------------


def _scrub(text: str) -> str:
    """Searchable text without credentials (``memory_store._without_secrets``).

    Applied BEFORE storage, never after: a secret that reached the index
    would be searchable, which is worse than it sitting in a 0600 file.
    Never raises.
    """
    try:
        from .memory_store import _without_secrets
        return _without_secrets(text)
    except Exception:
        return text


def _source_state(session_id: str):
    """The recorded (mtime, size) of a session's source documents, as a
    ``{kind: (mtime, size)}`` dict -- or None when nothing is indexed.

    ``kind`` is ``"session"`` or ``"archive"`` (the stat pair in the row's
    id is a detail of change detection, not part of the interface).
    """
    try:
        con = _open_index()
        try:
            rows = con.execute(
                "SELECT kind, mtime, size FROM sources WHERE session_id = ?",
                (session_id,)).fetchall()
        finally:
            con.close()
    except sqlite3.Error:
        return None
    if not rows:
        return None
    return {kind: (float(mtime), int(size)) for kind, mtime, size in rows}


def index_session(session_id: str):
    """(Re)index one session's archived material. Never raises.

    Incremental: a source file whose mtime and size match the recorded
    state is not re-read. Changed sources replace their old documents --
    re-indexing never duplicates. Returns True when something was indexed,
    False when there was nothing to do or nothing to index.
    """
    if not session_id:
        return False
    try:
        with _write_lock:
            return _index_session_locked(session_id)
    except Exception:
        # The index is a convenience: its own failure must never break the
        # caller (least of all session end).
        return False


def _index_session_locked(session_id: str) -> bool:
    docs = []
    docs.extend(_session_json_documents(session_id))
    docs.extend(_archive_documents(session_id))
    if not docs:
        return False
    con = _open_index()
    try:
        con.execute("BEGIN")
        known = {
            row[0]: (float(row[1]), int(row[2])) for row in con.execute(
                "SELECT kind, mtime, size FROM sources WHERE session_id = ?",
                (session_id,)).fetchall()}
        changed = False
        for doc_id, created, title, text in docs:
            kind, _, _ = doc_id.partition(":")
            stat_mtime = float(doc_id.split(":")[1])
            stat_size = int(doc_id.split(":")[2])
            if known.get(kind) == (stat_mtime, stat_size):
                continue  # unchanged since the last pass
            # One changed source rewrites the session's whole entry: FTS5
            # has no delete-by-value on unindexed columns, and the
            # documents of one session belong together anyway.
            con.execute("DELETE FROM documents WHERE session_id = ?",
                        (session_id,))
            con.execute("DELETE FROM sources WHERE session_id = ?",
                        (session_id,))
            for doc_id2, created2, title2, text2 in docs:
                con.execute(
                    "INSERT INTO documents "
                    "(session_id, created_at, title, text) "
                    "VALUES (?, ?, ?, ?)",
                    (session_id, created2, title2, _scrub(text2)))
                kind2, _, _ = doc_id2.partition(":")
                m2 = float(doc_id2.split(":")[1])
                s2 = int(doc_id2.split(":")[2])
                con.execute(
                    "INSERT OR REPLACE INTO sources "
                    "(session_id, kind, mtime, size) VALUES (?, ?, ?, ?)",
                    (session_id, kind2, m2, s2))
            changed = True
            break
        con.commit()
        return changed
    except sqlite3.Error:
        con.rollback()
        raise
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------


def _snippet(text: str, term: str, width: int = MAX_SNIPPET) -> str:
    """A context window around the first occurrence of ``term`` in
    ``text``, capped at ``width`` characters (half of it context on each
    side). Falls back to the head of the text when nothing matches."""
    if not text:
        return ""
    pos = text.lower().find(term.lower())
    if pos < 0:
        return text[:width]
    half = max(width // 2 - len(term) // 2, 0)
    start = max(pos - half, 0)
    end = min(start + width, len(text))
    out = text[start:end]
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(text) else ""
    return f"{prefix}{out}{suffix}"


def _date_str(created: float) -> str:
    """UTC YYYY-MM-DD of a session's creation timestamp."""
    try:
        return time.strftime("%Y-%m-%d", time.gmtime(created or 0.0))
    except (ValueError, OverflowError, OSError):
        return ""


def _sanitize_query(query: str) -> str:
    """A query safe to hand to FTS5's MATCH (quoted phrase) or LIKE.

    The raw query never reaches SQL as syntax: the FTS5 path wraps every
    term in double quotes, the LIKE path percent-escapes.
    """
    terms = [t for t in re.split(r"\s+", (query or "").strip()) if t]
    if _fts5_available():
        return " ".join(f'"{t.replace(chr(34), "")}"' for t in terms)
    return " ".join(terms)


BACKFILL_BATCH = 10


def backfill_missing(batch: int = BACKFILL_BATCH) -> int:
    """Index sessions that exist on disk but are missing from the index.

    Self-healing for a measured gap (2026-09-29): every session-end index
    pass failed silently on the real machine (``index_session`` never
    raises, by design), leaving 25 finished sessions unsearchable and the
    index DB absent. This pass finds the gap the other way round: scan the
    sessions directory, index what the ``sources`` table does not know,
    bounded to ``batch`` sessions so a search never stalls on a machine
    with thousands of files. Corrupt files are skipped by ``index_session``
    itself (its never-raise contract). Returns how many sessions were
    newly indexed.
    """
    try:
        d = Path(session_store._SESSIONS_DIR)
        if not d.is_dir():
            return 0
        con = _open_index()
        try:
            known = {row[0] for row in con.execute(
                "SELECT DISTINCT session_id FROM sources").fetchall()}
        finally:
            con.close()
        missing = [p.stem for p in d.glob("*.json")
                   if p.stem not in known and p.stem != ""]
        if not missing:
            return 0
        n = 0
        for sid in missing[:max(1, int(batch))]:
            if index_session(sid):
                n += 1
        return n
    except Exception:
        # Same contract as everything else here: the index is a
        # convenience, its failure must never break the caller.
        return 0


def search(query: str, limit: int = DEFAULT_LIMIT) -> list[Hit]:
    """Search indexed sessions. Never raises; empty query -> no hits.

    Result count is capped (``limit`` at most ``MAX_LIMIT``), and every
    snippet is capped at ``MAX_SNIPPET`` characters -- a search must stay
    cheap no matter how large the transcripts were.

    Self-healing: when nothing is indexed yet, one bounded backfill pass
    runs first (see ``backfill_missing``) so a silently failed session-end
    index cannot make sessions unfindable forever.
    """
    if not (query or "").strip():
        return []
    limit = max(1, min(int(limit or DEFAULT_LIMIT), MAX_LIMIT))
    try:
        if backfill_missing() > 0:
            pass  # healed below; the search itself is the retry
        return _search_locked(query, limit)
    except Exception:
        return []


def _search_locked(query: str, limit: int) -> list[Hit]:
    con = _open_index()
    try:
        if _fts5_available():
            match = _sanitize_query(query)
            rows = con.execute(
                "SELECT session_id, created_at, title, text, "
                "snippet(documents, 3, '', '', '...', 12) "
                "FROM documents WHERE documents MATCH ? "
                "ORDER BY created_at DESC LIMIT ?",
                (match, limit)).fetchall()
            hits = [Hit(r[0], _date_str(r[1] or 0.0), str(r[2] or ""),
                        str(r[4] or "")[:MAX_SNIPPET]) for r in rows]
            if hits:
                return hits
            # AND-semantics of FTS5 can be stricter than the user meant;
            # fall through to the substring scan below, which ORs nothing
            # but also requires only one term to match.
        term = _sanitize_query(query).split(" ")[0] if query.strip() else ""
        if not term:
            return []
        like = f"%{term.replace('%', '')}%"
        rows = con.execute(
            "SELECT session_id, created_at, title, text FROM documents "
            "WHERE text LIKE ? OR title LIKE ? "
            "ORDER BY created_at DESC LIMIT ?",
            (like, like, limit)).fetchall()
        # One hit per session: a session with both json and archive
        # documents must not crowd the result list.
        seen: dict[str, Hit] = {}
        for session_id, created, title, text in rows:
            if session_id in seen:
                continue
            seen[session_id] = Hit(
                session_id, _date_str(created or 0.0), str(title or ""),
                _snippet(str(text or ""), term))
        return list(seen.values())[:limit]
    finally:
        con.close()
