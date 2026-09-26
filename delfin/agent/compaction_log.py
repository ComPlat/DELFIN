"""One small append-only record per compaction, on disk.

Real operation is where compaction losses show, and nobody could see them:
``last_compaction_info`` holds only the LAST event, in RAM, and the archived
transcript holds everything but is unbounded prose. This module is the
missing middle — one JSONL line per compaction event:

* timestamp, session id
* tokens before / after, messages replaced
* which sections the working-state block carried (headings and counts
  only — never contents, which stay in the archived transcript)

Written best-effort: a broken log must never break compaction itself.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from .state_paths import ensure_dir, open_append

#: Resolved at import time so the suite can redirect it
#: (``USER_STATE_SINKS`` in state_paths.py) — see the sink table's header
#: comment for why every module-level state path must be listed there.
_LOG_DIR = Path.home() / ".delfin" / "compaction_log"

_MAX_RECORD_CHARS = 2000


def _log_path(session_id: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in session_id)
    return ensure_dir(_LOG_DIR) / f"{safe or 'session'}.jsonl"


def record_compaction(
    session_id: str,
    *,
    kind: str,
    tokens_before: int | None = None,
    tokens_after: int | None = None,
    messages_compacted: int = 0,
    messages_trimmed: int = 0,
    pinned_kept: int = 0,
    forced: bool = False,
    note: str = "",
    state_block: str = "",
    compacted_texts: list[str] | None = None,
) -> str | None:
    """Append one compaction event to the session's log; return the line.

    ``state_block`` is kept as its SECTION HEADINGS and per-section counts
    only (never contents — those live on in the archived transcript).
    Returns ``None`` (never raises) when nothing recordable happened.
    """
    sections = _state_block_sections(state_block)
    if kind == "summary_unavailable" and messages_compacted == 0:
        sections = None  # still record: the failure IS the observable fact
    record = {
        "ts": round(time.time(), 3),
        "session": (session_id or "")[:64],
        "kind": kind,
        "tokens_before": tokens_before,
        "tokens_after": tokens_after,
        "messages_compacted": messages_compacted,
        "messages_trimmed": messages_trimmed,
        "pinned_kept": pinned_kept,
        "forced": bool(forced),
        "note": (note or "")[:200],
        "state_block_sections": sections,
        "n_compacted_texts": len(compacted_texts or []),
    }
    line = json.dumps(record, ensure_ascii=False)
    if len(line) > _MAX_RECORD_CHARS:
        record["note"] = ""
        line = json.dumps(record, ensure_ascii=False)
    try:
        with open_append(_log_path(session_id or "")) as fh:
            fh.write(line + "\n")
        return line
    except Exception:
        return None


def read_compactions(session_id: str) -> list[dict]:
    """The parsed records for one session, oldest first; [] on any error."""
    try:
        path = _log_path(session_id or "")
        if not path.exists():
            return []
        out = []
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except Exception:
                continue
        return out
    except Exception:
        return []


def _state_block_sections(block: str) -> dict[str, int]:
    """Heading -> line count for each section of a working-state block.

    Only the fixed headings the block builder emits (working_state.py)
    are counted; anything else in the text is ignored, so no message
    contents can leak into the log.
    """
    if not block:
        return {}
    headings = (
        "Open tasks (task tool):",
        "Standing instructions:",
        "Recent denials (do not retry):",
        "Operator refusals (do not ask again):",
        "Last test outcomes:",
        "Recently worked on:",
        "Files changed this session (change journal):",
    )
    lines = block.splitlines()
    counts: dict[str, int] = {}
    current: str | None = None
    for ln in lines:
        if ln in headings:
            current = ln
            counts.setdefault(current, 0)
            continue
        if current is not None:
            if ln.startswith("  ") and ln.strip():
                counts[current] += 1
            else:
                current = None
    return counts
