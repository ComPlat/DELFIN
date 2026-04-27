"""Read-only inspector for Claude Code ``settings.json``.

Surfaces the most relevant fields (model, env, permissions.allow / deny,
output style) so users can see what's active without hunting through the
JSON.  We intentionally don't write here — settings are security-relevant
and the user should edit them in their editor.
"""
from __future__ import annotations

import html as _html
import json
from dataclasses import dataclass
from pathlib import Path


_DEFAULT_USER_PATH = Path.home() / ".claude" / "settings.json"


@dataclass(frozen=True)
class SettingsView:
    """Compact read-only view of a settings.json scope."""

    source_path: str
    model: str
    env: dict[str, str]
    allow: list[str]
    deny: list[str]
    output_style: str
    other_keys: list[str]    # top-level keys we don't render in detail


def _read_settings(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _extract(data: dict, source_path: str) -> SettingsView:
    permissions = data.get("permissions", {}) or {}
    allow = permissions.get("allow", []) or []
    deny = permissions.get("deny", []) or []
    env = data.get("env", {}) or {}
    handled = {"model", "env", "permissions", "outputStyle", "hooks"}
    other = sorted(k for k in data.keys() if k not in handled)
    return SettingsView(
        source_path=source_path,
        model=str(data.get("model", "")),
        env={str(k): str(v) for k, v in env.items()} if isinstance(env, dict) else {},
        allow=[str(x) for x in allow if isinstance(x, (str, int, float))],
        deny=[str(x) for x in deny if isinstance(x, (str, int, float))],
        output_style=str(data.get("outputStyle", "")),
        other_keys=other,
    )


def discover_settings(
    *,
    user_path: Path | None = None,
    project_path: Path | None = None,
) -> list[SettingsView]:
    """Return per-scope settings views.  Project first, then user."""
    paths: list[Path] = []
    if project_path is not None:
        paths.append(project_path)
    paths.append(user_path or _DEFAULT_USER_PATH)
    out: list[SettingsView] = []
    for path in paths:
        data = _read_settings(path)
        if not data:
            continue
        out.append(_extract(data, str(path)))
    return out


def render_settings_html(views: list[SettingsView]) -> str:
    """Render the inspector view; empty input → empty string."""
    if not views:
        return (
            '<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
            'font-weight:700;">Claude Settings</h4>'
            '<p style="margin:0;color:#9ca3af;font-size:11px;">'
            'No <code>settings.json</code> found at '
            '<code>~/.claude/</code> or <code>&lt;repo&gt;/.claude/</code>.</p>'
        )

    cards: list[str] = []
    for view in views:
        scope_label = "project" if "/.claude/" in view.source_path and \
            not view.source_path.startswith(str(Path.home())) else "user"
        rows: list[str] = []
        if view.model:
            rows.append(_kv_row("model", view.model))
        if view.output_style:
            rows.append(_kv_row("outputStyle", view.output_style))
        if view.env:
            env_text = ", ".join(f"{k}={v}" for k, v in sorted(view.env.items()))
            rows.append(_kv_row("env", env_text))
        if view.allow:
            rows.append(_list_row("permissions.allow", view.allow,
                                  chip_color="#10b981"))
        if view.deny:
            rows.append(_list_row("permissions.deny", view.deny,
                                  chip_color="#ef4444"))
        if view.other_keys:
            rows.append(_kv_row("other keys", ", ".join(view.other_keys)))

        if not rows:
            rows.append(
                '<div style="font-size:11px;color:#9ca3af;">'
                '(empty file)</div>'
            )

        cards.append(
            '<div style="border:1px solid #e5e7eb;border-radius:8px;'
            'background:#fff;padding:10px 12px;margin-bottom:8px;">'
            f'<div style="display:flex;align-items:center;gap:8px;'
            f'margin-bottom:6px;">'
            f'<span style="background:#dbeafe;color:#1e40af;font-size:10px;'
            f'font-weight:700;padding:1px 6px;border-radius:8px;'
            f'text-transform:uppercase;">{scope_label}</span>'
            f'<code style="font-size:11px;color:#6b7280;">'
            f'{_html.escape(view.source_path)}</code>'
            '</div>'
            + "".join(rows)
            + '</div>'
        )

    return (
        f'<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
        f'font-weight:700;">Claude Settings ({len(views)})</h4>'
        + "".join(cards)
    )


def _kv_row(label: str, value: str) -> str:
    return (
        '<div style="display:flex;gap:10px;font-size:11px;padding:3px 0;">'
        f'<span style="color:#6b7280;font-weight:600;min-width:140px;">'
        f'{_html.escape(label)}</span>'
        f'<span style="color:#111827;">{_html.escape(value)}</span>'
        '</div>'
    )


def _list_row(label: str, items: list[str], chip_color: str = "#3b82f6") -> str:
    chips = "".join(
        f'<span style="background:{chip_color}1a;color:{chip_color};'
        f'font-size:10px;padding:1px 7px;border-radius:8px;'
        f'margin:2px 4px 2px 0;display:inline-block;">'
        f'{_html.escape(item)}</span>'
        for item in items
    )
    return (
        '<div style="display:flex;gap:10px;font-size:11px;padding:3px 0;">'
        f'<span style="color:#6b7280;font-weight:600;min-width:140px;">'
        f'{_html.escape(label)}</span>'
        f'<span>{chips}</span>'
        '</div>'
    )


__all__ = [
    "SettingsView",
    "discover_settings",
    "render_settings_html",
]
