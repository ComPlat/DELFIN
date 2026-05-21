"""Read-only inventory of DELFIN scheduled agent tasks.

Schedules live in ``~/.delfin/schedules.json``::

    {
      "schedules": [
        {
          "name": "weekly-recalc-triage",
          "cron": "0 9 * * MON",
          "command": "/skill recalc-failed",
          "enabled": true,
          "note": "Find failed jobs every Monday 09:00."
        }
      ]
    }

Execution is intentionally out of scope for the dashboard.  The user
runs a system cron (or systemd timer) that consumes this file and
launches DELFIN agent runs accordingly.  The dashboard only renders
what's there so users have a clear picture of their automation.
"""
from __future__ import annotations

import html as _html
import json
from dataclasses import dataclass
from pathlib import Path


_DEFAULT_PATH = Path.home() / ".delfin" / "schedules.json"


@dataclass(frozen=True)
class Schedule:
    """One scheduled agent task."""

    name: str
    cron: str
    command: str
    enabled: bool
    note: str
    source_path: str


def _read(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def discover_schedules(path: Path | None = None) -> list[Schedule]:
    """Return all parseable schedule entries from ``path``."""
    p = path or _DEFAULT_PATH
    data = _read(p)
    raw = data.get("schedules", []) or []
    out: list[Schedule] = []
    for entry in raw:
        if not isinstance(entry, dict):
            continue
        out.append(Schedule(
            name=str(entry.get("name", "") or ""),
            cron=str(entry.get("cron", "") or ""),
            command=str(entry.get("command", "") or ""),
            enabled=bool(entry.get("enabled", True)),
            note=str(entry.get("note", "") or ""),
            source_path=str(p),
        ))
    return out


def render_schedules_html(items: list[Schedule]) -> str:
    """Render schedule inventory as compact HTML."""
    if not items:
        return (
            '<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
            'font-weight:700;">Scheduled Agent Tasks</h4>'
            '<p style="margin:0;color:#9ca3af;font-size:11px;">'
            'No schedules configured. Create '
            '<code>~/.delfin/schedules.json</code> with a '
            '<code>{"schedules": [...]}</code> array.</p>'
        )

    rows: list[str] = []
    for entry in items:
        if entry.enabled:
            badge = ('<span style="background:#10b981;color:#fff;font-size:10px;'
                     'padding:1px 6px;border-radius:8px;font-weight:700;">ON</span>')
        else:
            badge = ('<span style="background:#9ca3af;color:#fff;font-size:10px;'
                     'padding:1px 6px;border-radius:8px;font-weight:700;">OFF</span>')
        name = _html.escape(entry.name or "(unnamed)")
        cron = _html.escape(entry.cron or "(no schedule)")
        cmd = _html.escape(entry.command[:200])
        note = entry.note
        note_html = (
            f'<div style="font-size:10px;color:#9ca3af;margin-top:2px;">'
            f'{_html.escape(note)}</div>' if note else ""
        )
        rows.append(
            '<li style="padding:6px 8px;border-bottom:1px solid #f3f4f6;'
            'font-size:11px;">'
            f'<div style="display:flex;gap:10px;align-items:center;">'
            f'{badge}'
            f'<span style="color:#111827;font-weight:600;">{name}</span>'
            f'<code style="background:#fef3c7;color:#92400e;padding:1px 6px;'
            f'border-radius:4px;">{cron}</code>'
            '</div>'
            f'<div style="color:#374151;font-size:11px;margin-top:2px;'
            f'padding-left:6px;">→ <code>{cmd}</code></div>'
            + note_html +
            '</li>'
        )

    return (
        f'<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
        f'font-weight:700;">Scheduled Agent Tasks ({len(items)})</h4>'
        f'<div style="border:1px solid #e5e7eb;border-radius:8px;'
        f'background:#fff;padding:6px 4px;">'
        f'<ul style="list-style:none;margin:0;padding:0;">'
        + "".join(rows) +
        '</ul></div>'
    )


__all__ = ["Schedule", "discover_schedules", "render_schedules_html"]
