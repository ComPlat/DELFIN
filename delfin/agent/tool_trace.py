"""Per-session tool-call trace for the DELFIN agent.

Records every tool the agent calls — name, input, truncated output, duration,
ok/error — to ``~/.delfin/tool_traces/<session>.jsonl`` so a turn can be
replayed / audited, and so bug reports ship the exact sequence of actions the
agent took (the missing piece when diagnosing a failed session).

Best-effort and dependency-free: never raises, caps per-entry output and
per-file size, and silently no-ops on any IO error.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

_DIR = Path.home() / ".delfin" / "tool_traces"
_OUT_CAP = 2000             # chars kept per input/output
_MAX_BYTES = 4 * 1024 * 1024   # trim the file when it grows past this
_KEEP_TAIL = 1500          # lines kept when trimming


def _safe(session: str) -> str:
    s = "".join(c if (c.isalnum() or c in "-_.") else "_"
                for c in (session or "default"))
    return s[:80] or "default"


def trace_path(session: str) -> Path:
    return _DIR / f"{_safe(session)}.jsonl"


def _to_text(v: Any) -> str:
    if isinstance(v, (dict, list)):
        try:
            return json.dumps(v, ensure_ascii=False)
        except Exception:
            return str(v)
    return str(v if v is not None else "")


def record(
    session: str,
    *,
    tool: str,
    tool_input: Any = "",
    output: Any = "",
    duration_ms: int = 0,
    ok: bool = True,
    error: str = "",
    round: int = 0,
) -> None:
    """Append one tool-call entry to the session trace. Never raises."""
    try:
        _DIR.mkdir(parents=True, exist_ok=True)
        p = trace_path(session)
        # Bound the file: when it gets big, keep only the recent tail.
        try:
            if p.exists() and p.stat().st_size > _MAX_BYTES:
                lines = p.read_text(encoding="utf-8").splitlines()[-_KEEP_TAIL:]
                p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        except Exception:
            pass
        entry = {
            "ts": time.time(),
            "tool": str(tool),
            "input": _to_text(tool_input)[:_OUT_CAP],
            "output": _to_text(output)[:_OUT_CAP],
            "duration_ms": int(duration_ms),
            "ok": bool(ok),
            "error": str(error or "")[:300],
        }
        with p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except Exception:
        pass


def read(session: str, *, last_n: int | None = None) -> list[dict]:
    """Return the trace entries for ``session`` (optionally the last N)."""
    try:
        lines = trace_path(session).read_text(encoding="utf-8").splitlines()
    except Exception:
        return []
    if last_n:
        lines = lines[-int(last_n):]
    out: list[dict] = []
    for ln in lines:
        try:
            out.append(json.loads(ln))
        except Exception:
            continue
    return out


def clear(session: str) -> None:
    try:
        trace_path(session).unlink()
    except Exception:
        pass


def format_summary(entries: list[dict], *, limit: int = 30) -> str:
    """One line per call: ✓/✗ tool(args) → output… (Nms). For /trace + reports."""
    if not entries:
        return "(no tool calls recorded)"
    rows = []
    for e in entries[-limit:]:
        mark = "✓" if e.get("ok", True) else "✗"
        inp = (e.get("input") or "")[:60]
        out = (e.get("error") or e.get("output") or "")[:70].replace("\n", " ")
        dur = e.get("duration_ms", 0)
        rows.append(f"  {mark} {e.get('tool','?')}({inp}) → {out} ({dur}ms)")
    head = f"Tool trace — {len(entries)} call(s)" + (
        f", showing last {limit}" if len(entries) > limit else "")
    return head + ":\n" + "\n".join(rows)


def format_panel_html(entries: list[dict], *, limit: int = 12) -> str:
    """Compact HTML feed of the most-recent tool calls for the dashboard's live
    tool-trace panel — newest first, ✓/✗ glyph, tool name, short input, dur.
    Empty string when there are no calls (panel hides itself)."""
    from html import escape as _esc
    if not entries:
        return ""
    rows = []
    for e in entries[-limit:][::-1]:          # newest first
        ok = bool(e.get("ok", True))
        glyph, colour = ("✓", "#2e7d32") if ok else ("✗", "#c62828")
        tool = _esc(str(e.get("tool", "?")))
        inp = _esc(_to_text(e.get("input", ""))[:48])
        dur = int(e.get("duration_ms") or 0)
        dur_s = f"{dur / 1000:.1f}s" if dur >= 1000 else f"{dur}ms"
        rows.append(
            "<div style='font-family:monospace;font-size:11px;white-space:nowrap;"
            "overflow:hidden;text-overflow:ellipsis;'>"
            f"<span style='color:{colour};'>{glyph}</span> "
            f"<b>{tool}</b> <span style='color:#888;'>{inp}</span> "
            f"<span style='color:#aaa;'>· {dur_s}</span></div>"
        )
    head = (f"<div style='font-size:11px;color:#546e7a;'>"
            f"<b>🔧 Tools</b> · {len(entries)} call(s)</div>")
    return (head + "<div style='max-height:118px;overflow-y:auto;"
            "border:1px solid #2a2a2a;border-radius:4px;padding:2px 6px;"
            "margin-top:2px;'>" + "".join(rows) + "</div>")


__all__ = ["record", "read", "clear", "trace_path", "format_summary",
           "format_panel_html"]
