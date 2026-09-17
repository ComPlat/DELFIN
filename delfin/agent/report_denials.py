"""Where does the agent lose turns to refusals? A denials report.

Aggregates two REAL sources — no third one exists:

* "delfin.agent.tool_trace"  (durable, on disk): every tool call whose
  ``ok`` is False and whose error text looks like a refusal/blocked action
  ("not on the auto-allow list", "escapes workspace sandbox", "refusing to
  overwrite", …) is one lost turn. Counted per (category, tool) across all
  session traces in ``~/.delfin/tool_traces``, listed and read through
  ``tool_trace.sessions()`` and ``tool_trace.read(..., root=...)``, with
  ``aggregate_tools()``
  supplying the per-tool call/error totals so a refusal count can be read
  against a denominator.
* ``delfin.agent.security_events`` (in-process ring buffer): the permission
  gate's own ``record()``ed kinds, summarized via ``counts()``, ``recent()``
  and ``known_kinds()``. NOTE: these live in the memory of the process whose
  gate produced them — a CLI run in a fresh process usually sees an empty
  buffer. The tool_trace half is the durable record.

Best-effort like its sources: a missing file, a corrupt line or a missing
import is skipped, never raised on. Secrets never appear in the output:
example details are truncated and scrubbed (key/token/password assignments
and secret paths are redacted before anything is kept).

CLI: ``python -m delfin.agent.report_denials`` prints the report.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from . import security_events, tool_trace

# Refusal categories: (label, regex over the recorded error text). Order
# matters — first match wins. Kept deliberately coarse: the goal is "where
# do turns die", not a taxonomy of the permission system.
_CATEGORIES: list[tuple[str, re.Pattern[str]]] = [
    ("allowlist", re.compile(
        r"not on the (auto-)?allow|requires approval|approval", re.I)),
    ("sandbox", re.compile(
        r"escapes workspace sandbox|outside the allowed workspace roots|"
        r"outside an allowed|not within the allowed", re.I)),
    ("read_only", re.compile(r"read-only|readonly", re.I)),
    ("write_guard", re.compile(
        r"refusing to overwrite|prior read_file|read_file first|"
        r"read it first", re.I)),
    ("secret", re.compile(r"\.ssh|\.env|credential|\.pem|\.key\b|secret", re.I)),
    ("denied", re.compile(r"denied|refus|blocked|permission", re.I)),
]

# Anything that is an approval-window timeout is NOT a refusal (the user was
# simply away) — excluded before categorisation.
_NOT_A_REFUSAL = re.compile(
    r"TIMED OUT|approval window|the user is away", re.I)

_MAX_EXAMPLES = 3
_DETAIL_CAP = 160          # chars kept per example, post-scrub

_SECRET_SUB = re.compile(
    r"(?i)(api[_-]?key|token|password|secret|bearer)\s*[\"']?\s*[:=]\s*\S+")
_SECRET_PATH = re.compile(r"(?i)\S*(\.ssh/|\.env|\.pem|\.key|credentials)\S*")


def _scrub(text: str) -> str:
    """Redact secret-like assignments and paths; never raises."""
    try:
        text = _SECRET_SUB.sub(r"\1=«redacted»", text or "")
        text = _SECRET_PATH.sub("«secret-path»", text)
        return text.replace("\n", " ")[:_DETAIL_CAP]
    except Exception:
        return ""


def _categorize(error: str) -> str | None:
    if _NOT_A_REFUSAL.search(error or ""):
        return None
    for label, pat in _CATEGORIES:
        if pat.search(error or ""):
            return label
    return None


def collect(
    *,
    dir_path: Path | None = None,
    include_security_events: bool = True,
    max_sessions: int = 0,
) -> dict:
    """Ranked refusals/blocked actions by kind and tool. Never raises.

    Returns::

        {
          "sessions": N, "entries": N, "refusals": N,
          "ranked": [ {category, tool, count, share,
                       examples: [str, …]}, … ],       # sorted by count
          "tool_totals": {tool: {calls, errors}},      # from aggregate_tools
          "security_events": {"total", "blocked"} | None,
        }
    """
    out: dict[str, Any] = {
        "sessions": 0, "entries": 0, "refusals": 0,
        "ranked": [], "tool_totals": {}, "security_events": None,
    }
    try:
        base = Path(dir_path) if dir_path else tool_trace._DIR
        # The trace's own listing and reading: this module used to walk
        # the directory and parse the JSONL itself, because read() could
        # not be pointed at another directory. It can now.
        sessions = sorted(tool_trace.sessions(root=base))
        if max_sessions:
            sessions = sessions[-int(max_sessions):]
    except Exception:
        return out

    groups: dict[tuple[str, str], dict] = {}
    for sess in sessions:
        entries = tool_trace.read(sess, root=base)
        out["sessions"] += 1
        for e in entries:
            out["entries"] += 1
            if e.get("ok", True):
                continue
            err = str(e.get("error") or "")
            cat = _categorize(err)
            if cat is None:
                continue
            tool = str(e.get("tool") or "?")
            out["refusals"] += 1
            g = groups.setdefault(
                (cat, tool), {"count": 0, "examples": []})
            g["count"] += 1
            if len(g["examples"]) < _MAX_EXAMPLES:
                detail = _scrub(err) or "(no error text)"
                inp = _scrub(str(e.get("input") or ""))
                g["examples"].append(
                    f"{inp} → {detail}" if inp else detail)

    ranked = [
        {"category": cat, "tool": tool, "count": g["count"],
         "share": (g["count"] / out["refusals"]) if out["refusals"] else 0.0,
         "examples": g["examples"]}
        for (cat, tool), g in groups.items()
    ]
    ranked.sort(key=lambda r: (-r["count"], r["category"], r["tool"]))
    out["ranked"] = ranked

    try:
        for row in tool_trace.aggregate_tools(dir_path=base):
            out["tool_totals"][row["tool"]] = {
                "calls": row["calls"], "errors": row["errors"]}
    except Exception:
        pass

    if include_security_events:
        try:
            c = security_events.counts()
            if c.get("total"):
                out["security_events"] = {
                    "total": c["total"], "blocked": c["blocked"],
                    "recent": [
                        {"kind": e.kind, "tool": e.tool, "blocked": e.blocked}
                        for e in security_events.recent(limit=10)],
                    "known_kinds": sorted(security_events.known_kinds()),
                }
        except Exception:
            pass
    return out


def format_text(data: dict) -> str:
    """Compact ranked table for ``collect()`` output."""
    try:
        head = (f"Denials report — {data.get('sessions', 0)} session(s), "
                f"{data.get('entries', 0)} tool call(s), "
                f"{data.get('refusals', 0)} refusal(s)/blocked")
        ranked = list(data.get("ranked") or [])
        if not ranked:
            return head + "\n(no refusals recorded)"
        lines = [head, ""]
        lines.append(f"{'count':>5}  {'share':>5}  {'category':<12} tool")
        totals = data.get("tool_totals") or {}
        for r in ranked:
            tool = str(r.get("tool") or "?")
            tot = totals.get(tool)
            denom = f"/{tot['calls']}" if tot else ""
            lines.append(
                f"{r.get('count', 0):>5}  "
                f"{100.0 * float(r.get('share') or 0.0):>4.0f}%  "
                f"{str(r.get('category') or '?'):<12} {tool}{denom}")
        lines.append("")
        lines.append("examples:")
        for r in ranked:
            for ex in (r.get("examples") or [])[:_MAX_EXAMPLES]:
                lines.append(f"  [{r.get('category')}/{r.get('tool')}] {ex}")
        se = data.get("security_events")
        if se:
            lines.append("")
            lines.append(
                f"security events (in-process buffer): {se.get('total')} total, "
                f"{se.get('blocked')} blocked")
            for e in se.get("recent") or []:
                mark = "⛔" if e.get("blocked") else "•"
                lines.append(f"  {mark} {e.get('kind')} ({e.get('tool')})")
        return "\n".join(lines)
    except Exception:
        return "(denials report: could not format)"


def main() -> None:
    print(format_text(collect()))


if __name__ == "__main__":
    main()
