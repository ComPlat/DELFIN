"""Helpers around Claude Code hooks (read + write ``settings.json``).

Hooks are configured in ``~/.claude/settings.json`` (user scope) or
``<repo>/.claude/settings.json`` (project scope).  This module provides
pure functions to discover, read, render, and (optionally) add/remove
hook entries — separated from any widget code so it stays testable.

The supported event types follow the Claude Code spec:
``PreToolUse``, ``PostToolUse``, ``Stop``, ``UserPromptSubmit``,
``Notification``.  Each event holds a list of ``{matcher, hooks: [...]}``.
"""
from __future__ import annotations

import html as _html
import json
from dataclasses import dataclass
from pathlib import Path


_DEFAULT_USER_PATH = Path.home() / ".claude" / "settings.json"
SUPPORTED_EVENTS: tuple[str, ...] = (
    "PreToolUse",
    "PostToolUse",
    "Stop",
    "UserPromptSubmit",
    "Notification",
)


@dataclass(frozen=True)
class HookEntry:
    """A single hook command attached to an event."""

    event: str          # one of SUPPORTED_EVENTS
    matcher: str        # regex/glob, e.g. "Edit|Write" or "*"
    hook_type: str      # "command" (only type currently spec'd)
    command: str        # shell command to run
    source_path: str    # which settings.json this came from


def project_settings_path(repo_dir: Path | str) -> Path:
    """Project-scoped Claude settings path."""
    return Path(repo_dir) / ".claude" / "settings.json"


def _read_settings(path: Path) -> dict:
    """Load a settings.json file safely.  Missing/corrupt → empty dict."""
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def discover_hooks(
    *,
    user_path: Path | None = None,
    project_path: Path | None = None,
) -> list[HookEntry]:
    """Return all hooks from user + project settings, in priority order.

    Project hooks come first (they take priority over user hooks).
    """
    out: list[HookEntry] = []
    paths: list[Path] = []
    if project_path is not None:
        paths.append(project_path)
    paths.append(user_path or _DEFAULT_USER_PATH)

    for path in paths:
        data = _read_settings(path)
        hooks_block = data.get("hooks", {}) or {}
        for event, entries in hooks_block.items():
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                matcher = str(entry.get("matcher", "*"))
                for hk in (entry.get("hooks") or []):
                    if not isinstance(hk, dict):
                        continue
                    out.append(HookEntry(
                        event=event,
                        matcher=matcher,
                        hook_type=str(hk.get("type", "command")),
                        command=str(hk.get("command", "")),
                        source_path=str(path),
                    ))
    return out


def add_hook(
    event: str,
    matcher: str,
    command: str,
    *,
    settings_path: Path,
    hook_type: str = "command",
) -> bool:
    """Append a hook to ``settings_path``, creating the file if needed.

    Returns True if the hook was newly added (False if already present).
    """
    if event not in SUPPORTED_EVENTS:
        raise ValueError(
            f"Unsupported hook event {event!r}. "
            f"Allowed: {', '.join(SUPPORTED_EVENTS)}"
        )
    data = _read_settings(settings_path)
    hooks_block = data.setdefault("hooks", {})
    event_list = hooks_block.setdefault(event, [])

    # Look for an existing matcher group; reuse it if found
    target_group = None
    for group in event_list:
        if isinstance(group, dict) and group.get("matcher") == matcher:
            target_group = group
            break
    if target_group is None:
        target_group = {"matcher": matcher, "hooks": []}
        event_list.append(target_group)

    new_hook = {"type": hook_type, "command": command}
    if new_hook in (target_group.get("hooks") or []):
        return False
    target_group.setdefault("hooks", []).append(new_hook)

    settings_path.parent.mkdir(parents=True, exist_ok=True)
    settings_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return True


def remove_hook(
    event: str,
    matcher: str,
    command: str,
    *,
    settings_path: Path,
) -> bool:
    """Remove a hook from ``settings_path``.

    Returns True if removed, False if not found.  Empty matcher groups
    and empty event arrays are pruned to keep the file clean.
    """
    data = _read_settings(settings_path)
    hooks_block = data.get("hooks") or {}
    event_list = hooks_block.get(event)
    if not isinstance(event_list, list):
        return False

    removed = False
    for group in list(event_list):
        if not isinstance(group, dict) or group.get("matcher") != matcher:
            continue
        kept = [
            h for h in (group.get("hooks") or [])
            if not (isinstance(h, dict) and h.get("command") == command)
        ]
        if len(kept) != len(group.get("hooks") or []):
            removed = True
        if kept:
            group["hooks"] = kept
        else:
            event_list.remove(group)

    if removed:
        if not event_list:
            hooks_block.pop(event, None)
        if not hooks_block:
            data.pop("hooks", None)
        settings_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return removed


def render_hooks_html(hooks: list[HookEntry]) -> str:
    """Render a compact read-only inventory of configured hooks."""
    if not hooks:
        return (
            '<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
            'font-weight:700;">Claude Hooks</h4>'
            '<p style="margin:0;color:#9ca3af;font-size:11px;">'
            'No hooks configured. Add them to '
            '<code>~/.claude/settings.json</code> or '
            '<code>&lt;repo&gt;/.claude/settings.json</code>.</p>'
        )

    by_event: dict[str, list[HookEntry]] = {}
    for hook in hooks:
        by_event.setdefault(hook.event, []).append(hook)

    sections: list[str] = []
    for event, entries in by_event.items():
        rows: list[str] = []
        for entry in entries:
            cmd = _html.escape(entry.command[:200])
            matcher = _html.escape(entry.matcher)
            src = _html.escape(entry.source_path)
            rows.append(
                '<li style="padding:6px 8px;border-bottom:1px solid #f3f4f6;'
                'font-size:11px;">'
                f'<code style="background:#fef3c7;color:#92400e;padding:1px 6px;'
                f'border-radius:4px;">{matcher}</code>'
                f' <span style="color:#374151;">→ {cmd}</span>'
                f'<div style="color:#9ca3af;font-size:10px;margin-top:2px;">'
                f'{src}</div>'
                '</li>'
            )
        sections.append(
            f'<div style="margin-bottom:10px;">'
            f'<div style="font-size:11px;font-weight:700;color:#3b82f6;'
            f'text-transform:uppercase;letter-spacing:0.4px;'
            f'padding:3px 8px;">{_html.escape(event)} '
            f'({len(entries)})</div>'
            f'<ul style="list-style:none;margin:0;padding:0;">'
            + "".join(rows) + '</ul></div>'
        )

    return (
        f'<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
        f'font-weight:700;">Claude Hooks ({len(hooks)})</h4>'
        f'<div style="border:1px solid #e5e7eb;border-radius:8px;'
        f'background:#fff;padding:10px 12px;">'
        + "".join(sections) +
        '</div>'
    )


__all__ = [
    "HookEntry",
    "SUPPORTED_EVENTS",
    "project_settings_path",
    "discover_hooks",
    "add_hook",
    "remove_hook",
    "render_hooks_html",
]
