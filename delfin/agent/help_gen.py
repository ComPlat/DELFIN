"""Generated help for the dashboard slash-command surface.

The dashboard keeps ONE catalog of slash commands — a palette table of
``(category, command, summary, has_args)`` rows. The old ``/help``
handler was a second, hand-written copy of that knowledge and drifted
(missing commands, stale summaries). This module renders help text FROM
the palette rows verbatim, so the palette and ``/help`` can never
disagree, and exposes a coverage check that flags drift between the
palette and the dispatcher's registered slash prefixes.

Both helpers never raise: help is a read-only convenience surface and
must not be able to break the chat loop.
"""

from __future__ import annotations

from typing import Any


def _normalized_rows(palette_rows: Any) -> list[tuple[str, str, str, bool]]:
    """Coerce arbitrary input into well-formed palette rows.

    Malformed rows (too short, non-indexable, command without a leading
    slash) are skipped instead of failing the whole help page."""
    try:
        iterable = list(palette_rows or ())
    except TypeError:
        return []
    rows: list[tuple[str, str, str, bool]] = []
    for row in iterable:
        try:
            category = str(row[0])
            command = str(row[1])
            summary = str(row[2])
            has_args = bool(row[3]) if len(row) > 3 else False
        except Exception:
            continue
        if not command.startswith("/"):
            continue
        rows.append((category, command, summary, has_args))
    return rows


def generate_help(
    palette_rows: Any,
    *,
    category: str = "",
    search: str = "",
) -> str:
    """Grouped, aligned help text rendered from the palette rows.

    - Groups by category in first-appearance order (the palette's own
      curated order), commands aligned in one column across all groups.
    - ``has_args`` rows get a trailing ellipsis marker.
    - ``category`` filters by exact category name (case-insensitive);
      ``search`` filters by substring over command + summary + category
      (the palette's own filter semantics). Both may be combined.

    Never raises; renders the rows verbatim with no command knowledge of
    its own."""
    try:
        rows = _normalized_rows(palette_rows)
        cat = (category or "").strip().lower()
        query = (search or "").strip().lower()
        if cat:
            rows = [r for r in rows if r[0].lower() == cat]
        if query:
            rows = [r for r in rows
                    if query in f"{r[1]} {r[2]} {r[0]}".lower()]

        filters = []
        if cat:
            filters.append(f"category {category.strip()!r}")
        if query:
            filters.append(f"matching {search.strip()!r}")
        scope = " ".join(filters)
        if not rows:
            hint = "Try /help with no filter."
            if scope:
                return f"No commands {scope}. {hint}"
            return "No commands available."

        groups: dict[str, list[tuple[str, str, str, bool]]] = {}
        order: list[str] = []
        for row in rows:
            if row[0] not in groups:
                groups[row[0]] = []
                order.append(row[0])
            groups[row[0]].append(row)

        def _display(command: str, has_args: bool) -> str:
            return command + (" ..." if has_args else "")

        width = max(len(_display(r[1], r[3])) for r in rows)
        header = (f"Commands {scope} ({len(rows)}):" if scope
                  else f"Available commands ({len(rows)}):")
        lines = [header]
        for name in order:
            lines.append("")
            lines.append(f"{name}:")
            for _, command, summary, has_args in groups[name]:
                disp = _display(command, has_args)
                lines.append(f"  {disp:<{width}}  — {summary}")
        return "\n".join(lines)
    except Exception as exc:
        return f"Help unavailable: {exc}"


def coverage_gaps(palette_rows: Any, slash_prefixes: Any) -> dict:
    """Drift detector between the palette table and the dispatcher.

    A palette row is reachable only when its command's first token is a
    registered slash prefix; a registered prefix is discoverable only
    when at least one palette row documents it. Returns::

        {
          "ok": bool,                      # no drift in either direction
          "unregistered_commands": [...],  # rows the dispatcher ignores
          "unlisted_prefixes": [...],      # prefixes /help cannot show
        }

    Never raises (adds an ``"error"`` key on unexpected failure)."""
    try:
        rows = _normalized_rows(palette_rows)
        try:
            prefixes = {
                str(p) for p in (slash_prefixes or ())
                if str(p).startswith("/")
            }
        except TypeError:
            prefixes = set()
        covered: set[str] = set()
        unregistered: set[str] = set()
        for _, command, _, _ in rows:
            first = command.split()[0]
            covered.add(first)
            if first not in prefixes:
                unregistered.add(command)
        unlisted = sorted(prefixes - covered)
        result = {
            "ok": not unregistered and not unlisted,
            "unregistered_commands": sorted(unregistered),
            "unlisted_prefixes": unlisted,
        }
        return result
    except Exception as exc:
        return {
            "ok": False,
            "error": str(exc),
            "unregistered_commands": [],
            "unlisted_prefixes": [],
        }


__all__ = ["coverage_gaps", "generate_help"]
