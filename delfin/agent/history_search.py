"""Model-callable recall over a session's FULL conversation history.

DELFIN compacts long sessions aggressively (3-stage ladder). The full
pre-compaction transcripts are archived per session as append-only JSONL
(:func:`delfin.agent.session_store.archive_pre_compaction_transcript`),
but until now only the USER could browse them (``/session archive``,
``/session search`` in the dashboard). The model itself had no way to
page detail back in — after a few compactions, "what did we decide
earlier?" degrades into confabulation risk.

This module is the library behind the ``history_search`` /
``history_get`` agent tools:

- :func:`history_search` — BM25-ranked (substring fallback) search over
  BOTH the live message list (the engine's ``self.messages``, supplied
  by the caller — this module never imports the engine) AND every
  archived pre-compaction record for the session.
- :func:`history_get` — fetch the full text of one hit by its ``ref``.

It lives in its own module rather than in ``session_store`` because
session_store is the persistence layer; this is retrieval/ranking built
ON TOP of it and of ``memory_store``'s BM25. Keeping it separate avoids
coupling storage to ranking and keeps session_store import-cheap.

Guarantees:

- Never raises: corrupt JSONL lines / malformed records are skipped;
  unexpected failures return an empty result or an ``{"error": ...}``
  payload the model can read.
- Path-safe: session ids are validated (no separators, no ``..``)
  BEFORE any archive path is built.
- Cheap on big sessions: the archive is streamed line-by-line and the
  combined scan is bounded to the most recent :data:`SCAN_CAP` messages
  (a note is appended to the result when the cap truncated the scan).
- Deterministic: ties break newest-first; ordering is a total order.
"""

from __future__ import annotations

import json
import re
from collections import deque
from pathlib import Path
from typing import Any

from delfin.agent.memory_store import _tokenize, bm25_scores
from delfin.agent.session_store import _transcript_archive_dir

# Combined live+archive scan bound: only the most recent SCAN_CAP message
# records are tokenised/scored. Refs stay stable under the cap because
# they use absolute positions (live list index / archive line + message
# index), never positions within the capped window.
SCAN_CAP = 5000

# Per-message text bound for scanning/snippet extraction. history_get
# re-reads the message from its source, so the full text stays reachable
# even when scoring saw only this prefix.
_MAX_TEXT_CHARS = 20_000

_SNIPPET_CHARS = 200

# Same character family session_store uses when it sanitises ids for
# filenames (see save_handoff_brief) — but REJECTING instead of
# substituting: a substituted id would silently point at a different
# session's archive. No path separators, no leading dot, no "..".
_SESSION_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _valid_session_id(session_id: str) -> bool:
    sid = str(session_id or "")
    return bool(_SESSION_ID_RE.match(sid)) and ".." not in sid


def _archive_path(session_id: str) -> Path:
    """Archive file for ``session_id`` — same layout session_store writes."""
    return _transcript_archive_dir() / f"{session_id}.jsonl"


def _msg_text(m: dict[str, Any]) -> str:
    """Flatten a message's content to text, same rule the archiver uses."""
    c = m.get("content", "")
    if isinstance(c, str):
        return c
    try:
        return json.dumps(c, ensure_ascii=False)
    except (TypeError, ValueError):
        return str(c)


def _collect_entries(
    session_id: str,
    messages: list[dict[str, Any]] | None,
    cap: int,
) -> tuple[list[dict[str, Any]], bool]:
    """Gather scannable message entries, oldest -> newest, bounded to
    the most recent ``cap``. Returns ``(entries, capped)``.

    Archive first (chronologically older than the live list), then the
    live messages. The JSONL is streamed and each line parsed once;
    corrupt lines are skipped but still consume their line number, so
    ``record`` (the seq id inside archive refs) is a raw line index and
    stays stable as the append-only file grows.
    """
    seen = 0
    buf: deque[dict[str, Any]] = deque(maxlen=cap)

    p = _archive_path(session_id)
    if p.is_file():
        try:
            with p.open("r", encoding="utf-8", errors="replace") as fh:
                for seq, line in enumerate(fh):
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        continue   # torn/corrupt line — skip, keep seq
                    if not isinstance(rec, dict):
                        continue
                    msgs = rec.get("messages")
                    if not isinstance(msgs, list):
                        continue
                    ts = rec.get("compacted_at")
                    if not isinstance(ts, (int, float)):
                        ts = None
                    for i, m in enumerate(msgs):
                        if not isinstance(m, dict):
                            continue
                        seen += 1
                        buf.append({
                            "source": "archive",
                            "record": seq,
                            "index": i,
                            "role": str(m.get("role", "") or ""),
                            "ts": ts,
                            "text": _msg_text(m)[:_MAX_TEXT_CHARS],
                        })
        except OSError:
            pass

    for i, m in enumerate(messages or []):
        if not isinstance(m, dict):
            continue
        seen += 1
        ts = m.get("ts")
        if not isinstance(ts, (int, float)):
            ts = None
        buf.append({
            "source": "live",
            "index": i,
            "role": str(m.get("role", "") or ""),
            "ts": ts,
            "text": _msg_text(m)[:_MAX_TEXT_CHARS],
        })

    return list(buf), seen > len(buf)


def _snippet(text: str, terms: list[str]) -> str:
    """~200 chars around the earliest occurrence of any query term."""
    flat = " ".join(text.split()) or text
    low = flat.casefold()
    pos = -1
    for t in terms:
        t = str(t or "").casefold()
        if not t:
            continue
        i = low.find(t)
        if i >= 0 and (pos < 0 or i < pos):
            pos = i
    if pos < 0:
        pos = 0
    start = max(0, pos - _SNIPPET_CHARS // 2)
    end = min(len(flat), start + _SNIPPET_CHARS)
    start = max(0, end - _SNIPPET_CHARS)
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(flat) else ""
    return f"{prefix}{flat[start:end]}{suffix}"


def history_search(
    session_id: str,
    query: str,
    *,
    messages: list[dict[str, Any]] | None = None,
    max_results: int = 8,
) -> list[dict[str, Any]]:
    """Search this session's live messages + archived transcripts.

    Parameters
    ----------
    session_id : str
        The session whose transcript archive to search. Validated
        against path traversal; an invalid id yields a single
        ``{"error": ...}`` entry.
    query : str
        Ranked with BM25 (:func:`memory_store.bm25_scores`) when it
        yields query tokens; otherwise (short/stopword-only queries,
        e.g. "S1") a case-insensitive substring match.
    messages : list, optional
        The engine's LIVE message list — supplied by the caller so this
        module never imports the engine. ``None`` -> archive-only.
    max_results : int
        Maximum hits returned (default 8).

    Returns
    -------
    list of dict
        Hits sorted by score desc (ties newest-first), each::

            {"ref": "live:<i>" | "archive:<rec>:<i>",
             "source": "live" | "archive",
             "index": <message index>,           # + "record" for archive
             "role": ..., "ts": <float | None>,
             "snippet": ~200 chars around the best match,
             "score": <float>}

        When the scan was capped, a final ``{"note": ..., "capped":
        True}`` entry is appended. Never raises.
    """
    try:
        sid = str(session_id or "")
        if not _valid_session_id(sid):
            return [{"error": f"invalid session id: {sid!r}"}]
        try:
            max_results = int(max_results)
        except (TypeError, ValueError):
            max_results = 8
        if max_results <= 0:
            max_results = 8

        cap = max(1, int(SCAN_CAP))
        entries, capped = _collect_entries(sid, messages, cap)

        q = str(query or "")
        tokens = sorted(_tokenize(q))
        scored: list[tuple[float, int]] = []
        if entries and tokens:
            scores = bm25_scores(q, [e["text"] for e in entries])
            terms: list[str] = tokens
            scored = [(s, i) for i, s in enumerate(scores) if s > 0]
        elif entries:
            needle = q.strip()
            terms = [needle] if needle else []
            nlow = needle.casefold()
            for i, e in enumerate(entries):
                if not nlow or nlow in e["text"].casefold():
                    scored.append((1.0, i))
        else:
            terms = []

        # Score desc; ties newest-first (entries are chronological, so a
        # larger position index means a more recent message).
        scored.sort(key=lambda si: (-si[0], -si[1]))

        hits: list[dict[str, Any]] = []
        for s, i in scored[:max_results]:
            e = entries[i]
            if e["source"] == "live":
                hit: dict[str, Any] = {"ref": f"live:{e['index']}"}
            else:
                hit = {
                    "ref": f"archive:{e['record']}:{e['index']}",
                    "record": e["record"],
                }
            hit.update({
                "source": e["source"],
                "index": e["index"],
                "role": e["role"],
                "ts": e["ts"],
                "snippet": _snippet(e["text"], terms),
                "score": round(float(s), 4),
            })
            hits.append(hit)

        if capped:
            hits.append({
                "note": (
                    f"scan capped: only the most recent {cap} messages "
                    "(live + archived) were searched; older history was "
                    "not scanned"
                ),
                "capped": True,
            })
        return hits
    except Exception:
        return []


def _parse_ref(ref: str) -> tuple[str, int, int] | None:
    """Parse ``live:<i>`` / ``archive:<rec>:<i>`` -> (source, rec, idx)."""
    parts = str(ref or "").strip().split(":")
    try:
        if len(parts) == 2 and parts[0] == "live":
            idx = int(parts[1])
            if idx >= 0:
                return ("live", -1, idx)
        elif len(parts) == 3 and parts[0] == "archive":
            rec, idx = int(parts[1]), int(parts[2])
            if rec >= 0 and idx >= 0:
                return ("archive", rec, idx)
    except ValueError:
        pass
    return None


def _truncate(text: str, max_chars: int) -> tuple[str, bool]:
    """Head+tail truncation with an explicit omission marker."""
    if max_chars <= 0:
        max_chars = 4000
    if len(text) <= max_chars:
        return text, False
    head = max(1, int(max_chars * 0.6))
    tail = max(1, max_chars - head)
    omitted = len(text) - head - tail
    marker = f"\n... [truncated: {omitted} chars omitted] ...\n"
    return text[:head] + marker + text[-tail:], True


def history_get(
    session_id: str,
    ref: str,
    *,
    messages: list[dict[str, Any]] | None = None,
    max_chars: int = 4000,
) -> dict[str, Any]:
    """Fetch the full text of ONE :func:`history_search` hit by its ref.

    ``ref`` is the ``"live:<i>"`` / ``"archive:<rec>:<i>"`` string a
    search hit carries. Archive lookups stream the JSONL only up to the
    target line. Text longer than ``max_chars`` is returned head+tail
    with a ``[truncated: N chars omitted]`` marker.

    Returns ``{"ref", "source", "index", "role", "ts", "text",
    "truncated", "total_chars"}`` (+``"record"`` for archive refs) or
    ``{"error": ...}``. Never raises.
    """
    try:
        sid = str(session_id or "")
        if not _valid_session_id(sid):
            return {"error": f"invalid session id: {sid!r}"}
        parsed = _parse_ref(ref)
        if parsed is None:
            return {"error": (
                f"invalid ref: {str(ref)!r} — expected 'live:<index>' or "
                "'archive:<record>:<index>' as returned by history_search"
            )}
        source, rec_seq, idx = parsed
        try:
            max_chars = int(max_chars)
        except (TypeError, ValueError):
            max_chars = 4000

        if source == "live":
            if messages is None:
                return {"error": (
                    "live messages are not available in this context; "
                    "only archive:* refs can be resolved"
                )}
            if not (0 <= idx < len(messages)) or not isinstance(messages[idx], dict):
                return {"error": f"live message index {idx} out of range"}
            m = messages[idx]
            ts = m.get("ts")
            if not isinstance(ts, (int, float)):
                ts = None
            full = _msg_text(m)
            text, truncated = _truncate(full, max_chars)
            return {
                "ref": f"live:{idx}",
                "source": "live",
                "index": idx,
                "role": str(m.get("role", "") or ""),
                "ts": ts,
                "text": text,
                "truncated": truncated,
                "total_chars": len(full),
            }

        # archive: stream up to the target line only
        p = _archive_path(sid)
        if not p.is_file():
            return {"error": f"no transcript archive for session '{sid}'"}
        rec: Any = None
        try:
            with p.open("r", encoding="utf-8", errors="replace") as fh:
                for seq, line in enumerate(fh):
                    if seq < rec_seq:
                        continue
                    if seq > rec_seq:
                        break
                    line = line.strip()
                    if line:
                        try:
                            rec = json.loads(line)
                        except json.JSONDecodeError:
                            rec = None
                    break
        except OSError:
            return {"error": f"transcript archive for '{sid}' unreadable"}
        if not isinstance(rec, dict):
            return {"error": f"archive record {rec_seq} not found or unreadable"}
        msgs = rec.get("messages")
        if not isinstance(msgs, list) or not (0 <= idx < len(msgs)):
            return {"error": (
                f"message index {idx} out of range in archive record {rec_seq}"
            )}
        m = msgs[idx]
        if not isinstance(m, dict):
            return {"error": (
                f"archive record {rec_seq} message {idx} is malformed"
            )}
        ts = rec.get("compacted_at")
        if not isinstance(ts, (int, float)):
            ts = None
        full = _msg_text(m)
        text, truncated = _truncate(full, max_chars)
        return {
            "ref": f"archive:{rec_seq}:{idx}",
            "source": "archive",
            "record": rec_seq,
            "index": idx,
            "role": str(m.get("role", "") or ""),
            "ts": ts,
            "text": text,
            "truncated": truncated,
            "total_chars": len(full),
        }
    except Exception:
        return {"error": "history_get failed unexpectedly"}
