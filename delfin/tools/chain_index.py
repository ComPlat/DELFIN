"""The chain index (JC2 target 2): a per-user pointer list of chains.

When a Pipeline creates a chain it records one entry here — chain name,
work dir, creation time and (when the run happens inside an agent
session) that session's id, passed in explicitly by the caller that
already holds it (the same dispatch boundary that fills the step
journal — ``permissions.task_session_id``; never discovered from the
environment here).

Design rules:

- The index is a POINTER, never the truth.  Every consumer re-reads
  each chain's own ``chain_state.json`` for its current state; the
  index only says where to look.
- Per-user, under the DELFIN state dir, resolved PER CALL from
  ``Path.home()`` — never an import-time constant (the state-paths
  watcher test refuses those).
- Appended atomically (one JSON line, fsync), bounded: on each write,
  entries whose chain state file is gone or whose chain is fully
  finished are pruned and the file is capped at the newest MAX_ENTRIES.
- Never a tree walk of HOME or of a workspace: only the recorded work
  dirs are stat'ed, when the bound check prunes.
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

#: The index file lives directly under ``~/.delfin`` (a per-user list,
#: not per workspace: an agent resumes across workspaces and a crashed
#: chain's directory is where it is, not where this session runs).
INDEX_FILENAME = "chain_index.jsonl"

#: The index stays small: older entries — finished chains, vanished
#: chains — are pruned on write and the file is capped at this many
#: newest entries.
MAX_ENTRIES = 100


class ChainIndexError(RuntimeError):
    """The chain index could not be written (never fatal for the run).

    Raised by :func:`record` so the caller can say, once and visibly,
    why a chain is not indexed — the run itself continues.
    """


def index_path() -> Path:
    """The per-user chain index file, resolved per call.

    Resolved per call (wave-15a rule: a store path is never an
    import-time constant) so a redirected home — a test, a benchmark
    run — moves it with the home.
    """
    return Path.home() / ".delfin" / INDEX_FILENAME


def record(
    name: str,
    work_dir: str | Path,
    *,
    agent_session_id: Optional[str] = None,
    created_at: Optional[float] = None,
) -> bool:
    """Append one chain entry to the per-user index.

    Best-effort by design: an index failure must never break the run
    that created the chain — a warning is logged, the run goes on.
    Returns True when the entry was written.

    Live-run A (operator, 2026-11-16): inside the agent sandbox HOME is
    read-only, so this fails on EVERY chain started through bash — and
    a silent failure made the restart offer empty.  The exception's
    reason is therefore returned in the message of the raised
    :class:`ChainIndexError` so the caller (pipeline.run) can say, once,
    why nothing was indexed.
    """
    entry = {
        "name": str(name),
        "work_dir": str(work_dir),
        "agent_session_id": str(agent_session_id) if agent_session_id else None,
        "created_at": float(created_at if created_at is not None else time.time()),
    }
    try:
        path = index_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        _rewrite(path, [entry])
        return True
    except OSError as exc:
        logger.warning("chain index: could not record %s (%s)", name, exc)
        raise ChainIndexError(
            f"chain '{name}' not indexed: {exc}") from exc


def load(*, agent_session_id: Optional[str] = None) -> List[Dict[str, Any]]:
    """The recorded entries, newest last; a torn line is skipped with a
    log line.

    When *agent_session_id* is given, only that session's entries are
    returned (the crash-recovery offer shows a session ITS chains).
    """
    path = index_path()
    if not path.is_file():
        return []
    entries: List[Dict[str, Any]] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        logger.warning("chain index: unreadable (%s)", exc)
        return []
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        try:
            entry = json.loads(line)
            if not isinstance(entry, dict):
                raise ValueError("not an object")
        except ValueError as exc:
            logger.warning(
                "chain index: skipping torn line %d in %s (%s)",
                lineno, path, exc,
            )
            continue
        entries.append(entry)
    if agent_session_id is not None:
        entries = [e for e in entries
                   if e.get("agent_session_id") == str(agent_session_id)]
    return entries


def _entry_state(entry: Dict[str, Any]) -> str:
    """The chain's CURRENT state, read from its own state file.

    The index is a pointer: the state file is the truth.  ``gone`` when
    the work dir or the state file vanished; ``finished`` when every
    step carries the finished state; otherwise ``unfinished``.
    """
    from delfin.tools import chain_state as _cs

    wd = entry.get("work_dir")
    if not wd:
        return "gone"
    state_file = Path(wd) / "chain_state.json"
    if not state_file.is_file():
        return "gone"
    try:
        chain = _cs.ChainState.load(state_file)
    except Exception as exc:  # noqa: BLE001 - a corrupt file is a fact
        logger.warning(
            "chain index: state file %s unreadable (%s)", state_file, exc,
        )
        return "unfinished"
    if chain.steps and all(s.state == "finished" for s in chain.steps):
        return "finished"
    return "unfinished"


def _rewrite(path: Path, new_entries: List[Dict[str, Any]]) -> None:
    """Prune + append atomically: rewrite the bounded entry list.

    Vanished chains (no state file any more) are pruned eagerly — they
    are pure noise.  Finished chains STAY until the bound evicts them:
    `delfin chain list` should still show a chain that just finished,
    and the offer logic filters finished ones by reading each chain's
    own state anyway (the index is a pointer).  The file holds at most
    MAX_ENTRIES newest entries.  Atomic via temp-file + os.replace, so
    a crash mid-write leaves the previous index readable.
    """
    kept: List[Dict[str, Any]] = []
    for e in load():
        if _entry_state(e) != "gone":
            kept.append(e)
    kept.extend(new_entries)
    kept = kept[-MAX_ENTRIES:]
    payload = "".join(
        json.dumps(e, sort_keys=True) + "\n" for e in kept
    )
    fd, tmp_name = tempfile.mkstemp(
        prefix=".chain_index.", suffix=".tmp", dir=str(path.parent)
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


__all__ = ["INDEX_FILENAME", "MAX_ENTRIES", "index_path", "record", "load"]
