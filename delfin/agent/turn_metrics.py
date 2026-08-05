"""Per-turn timing metrics for the DELFIN agent.

Records one entry per agent turn — provider/model/role, total wall time,
**time-to-first-token (ttft)**, output size and tool-call count — to
``~/.delfin/turn_metrics/<session>.jsonl``.

The point is to make transient backend stalls *visible after the fact*. When a
trivial turn takes a minute and a half (e.g. the field report where Azure
GPT-5.4 answered "Hallo" in 92.7s), the breakdown tells whether the time went
into **waiting for the first token** (a backend/queue stall — high ttft, tiny
output, no tools) versus heavy generation or many tool rounds. That distinction
is exactly what a single "turn took 92.7s" number cannot give.

Best-effort and dependency-free: never raises, caps file size, no-ops on IO
error. Mirrors :mod:`tool_trace`.
"""

from __future__ import annotations

import os

import json
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "turn_metrics"
_MAX_BYTES = 2 * 1024 * 1024    # trim the file when it grows past this
_KEEP_TAIL = 2000               # lines kept when trimming

# A turn that took this long to its FIRST token, while emitting little and
# calling no tools, is flagged as a likely backend stall in summaries.
_SLOW_TTFT_MS = 20_000


def _safe(session: str) -> str:
    s = "".join(c if (c.isalnum() or c in "-_.") else "_"
                for c in (session or "default"))
    return s[:80] or "default"


def metrics_path(session: str) -> Path:
    return _DIR / f"{_safe(session)}.jsonl"


def record(
    session: str,
    *,
    provider: str = "",
    model: str = "",
    role: str = "",
    total_ms: int = 0,
    ttft_ms: int | None = None,
    output_chars: int = 0,
    tool_calls: int = 0,
    stopped: bool = False,
    error: str = "",
    input_tokens: int = 0,
    output_tokens: int = 0,
    cached_tokens: int = 0,
) -> None:
    """Append one turn-timing entry to the session log. Never raises.

    The token counts are here because a claim about them could not be
    checked. The prompt is deliberately byte-stable and re-sent every turn,
    justified by prefix caching making that nearly free -- and a probe once
    found ``cached_tokens = 0`` on the KIT endpoint, which would make the
    justification wrong and real shortening the only lever that helps.

    Neither could be confirmed later, because the number reached no
    persistent store: the only consumer renders a cache line ONLY when the
    count is above zero, so on an endpoint that reports none it silently
    vanished. Recording it makes the question answerable from data instead
    of from one probe somebody remembers running.
    """
    try:
        _DIR.mkdir(parents=True, exist_ok=True)
        p = metrics_path(session)
        try:
            if p.exists() and p.stat().st_size > _MAX_BYTES:
                lines = p.read_text(encoding="utf-8").splitlines()[-_KEEP_TAIL:]
                p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except Exception:
            pass
        entry = {
            "ts": time.time(),
            "provider": str(provider or ""),
            "model": str(model or ""),
            "role": str(role or ""),
            "total_ms": int(total_ms),
            "ttft_ms": (int(ttft_ms) if ttft_ms is not None else None),
            "output_chars": int(output_chars),
            "tool_calls": int(tool_calls),
            "stopped": bool(stopped),
            "error": str(error or "")[:300],
            "input_tokens": int(input_tokens or 0),
            "output_tokens": int(output_tokens or 0),
            "cached_tokens": int(cached_tokens or 0),
        }
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")
        # 0600 AFTER the write: on the first append the file does
        # not exist yet, so a chmod before it silently did nothing.
        # These files carry raw tool output, commands and paths;
        # they were created at the process umask (observed 0664)
        # and a bug report bundles them, adding group-read.
        try:
            os.chmod(p, 0o600)
        except OSError:
            pass
    except Exception:
        pass


def read(session: str, *, last_n: int | None = None) -> list[dict]:
    """Return the turn entries for ``session`` (optionally the last N)."""
    try:
        p = metrics_path(session)
        if not p.exists():
            return []
        out = []
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except Exception:
                continue
        return out[-last_n:] if last_n else out
    except Exception:
        return []


def is_stall(entry: dict) -> bool:
    """A turn dominated by waiting for the first token — little output, no
    tools — i.e. the backend, not the agent, ate the time."""
    try:
        ttft = entry.get("ttft_ms")
        return (ttft is not None and int(ttft) >= _SLOW_TTFT_MS
                and int(entry.get("tool_calls") or 0) == 0
                and int(entry.get("output_chars") or 0) <= 400)
    except Exception:
        return False


def _percentile(sorted_vals: list[int], q: float) -> int:
    """Nearest-rank percentile of a pre-sorted list (stdlib only)."""
    if not sorted_vals:
        return 0
    k = int(round((q / 100.0) * (len(sorted_vals) - 1)))
    k = max(0, min(len(sorted_vals) - 1, k))
    return sorted_vals[k]


def aggregate_turn_stats(
    window_days: float = 7,
    *,
    dir_path: Path | None = None,
    now: float | None = None,
) -> dict:
    """Windowed roll-up across ALL session files — the telemetry consumer
    the per-turn logs have been missing.

    Returns ``{turns, avg_ttft_ms, p90_ttft_ms, stalls, stopped_count}``
    where ``stalls`` counts entries flagged by :func:`is_stall` and the
    ttft numbers cover only turns that recorded a first token.  Entries
    outside the last ``window_days`` are skipped (``window_days <= 0``
    disables the window); entries without a parseable ``ts`` stay
    visible.  Best-effort: never raises — corrupt lines/files skipped.
    """
    empty = {"turns": 0, "avg_ttft_ms": 0, "p90_ttft_ms": 0,
             "stalls": 0, "stopped_count": 0}
    try:
        base = Path(dir_path) if dir_path else _DIR
        cutoff: float | None = None
        if window_days and float(window_days) > 0:
            cutoff = ((now if now is not None else time.time())
                      - float(window_days) * 86400.0)
        turns = stalls = stopped = 0
        ttfts: list[int] = []
        for fp in sorted(base.glob("*.jsonl")):
            try:
                lines = fp.read_text(encoding="utf-8").splitlines()
            except Exception:
                continue
            for ln in lines:
                try:
                    e = json.loads(ln)
                except Exception:
                    continue
                if not isinstance(e, dict):
                    continue
                if cutoff is not None:
                    try:
                        if float(e["ts"]) < cutoff:
                            continue
                    except (KeyError, TypeError, ValueError):
                        pass                    # no/bad ts → keep visible
                turns += 1
                ttft = e.get("ttft_ms")
                if ttft is not None:
                    try:
                        ttfts.append(int(ttft))
                    except (TypeError, ValueError):
                        pass
                if is_stall(e):
                    stalls += 1
                if e.get("stopped"):
                    stopped += 1
        ttfts.sort()
        return {
            "turns": turns,
            "avg_ttft_ms": (int(round(sum(ttfts) / len(ttfts)))
                            if ttfts else 0),
            "p90_ttft_ms": _percentile(ttfts, 90),
            "stalls": stalls,
            "stopped_count": stopped,
        }
    except Exception:
        return empty


def format_summary(entries: list[dict], *, limit: int = 30) -> str:
    """One line per turn; flags likely backend stalls. Empty when no entries."""
    if not entries:
        return ""
    rows = []
    for e in entries[-limit:]:
        total = int(e.get("total_ms") or 0)
        ttft = e.get("ttft_ms")
        ttft_s = f"{int(ttft)/1000:.1f}s" if ttft is not None else "—"
        flag = "  ⚠ backend-stall" if is_stall(e) else ""
        rows.append(
            f"{e.get('model','?')}  total={total/1000:.1f}s  "
            f"ttft={ttft_s}  out={int(e.get('output_chars') or 0)}ch  "
            f"tools={int(e.get('tool_calls') or 0)}{flag}"
        )
    return "\n".join(rows)


__all__ = ["metrics_path", "record", "read", "is_stall", "format_summary",
           "aggregate_turn_stats"]
