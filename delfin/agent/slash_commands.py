"""User-defined slash commands as markdown templates.

Drop a ``<name>.md`` file in any of the search locations below and the
chat will accept ``/<name>`` as a slash command that expands to the
file's body. Lighter weight than a full Skill — no Skill frontmatter,
no Skill-runtime, just text substitution.

Search order (later wins on collision):
  1. ``~/.delfin/commands/`` — user-global, survives projects
  2. ``<workspace>/.delfin/commands/`` — project-scoped
  3. ``~/.claude/commands/`` — canonical CLI compatibility layer

Frontmatter (optional, all keys optional)::

    ---
    description: Short one-line summary shown in /help and the palette.
    argument-hint: <pr-number>      # what the user types after /name
    ---

Inside the body, the following placeholders are substituted before the
text is injected into the chat:

  ``$ARGUMENTS``        full text the user typed after the command
  ``$1``, ``$2``, …    individual whitespace-separated tokens
  ``$@``                alias for ``$ARGUMENTS``

Unknown placeholders are left as-is so a literal ``$foo`` in a template
doesn't accidentally vanish.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


_NAME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_\-]{0,40}$")
_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)


@dataclass(frozen=True)
class SlashTemplate:
    name: str
    description: str
    argument_hint: str
    body: str
    source: Path

    @property
    def slash(self) -> str:
        return f"/{self.name}"


def _search_dirs(workspace: Path | str | None) -> list[Path]:
    """Order: user-delfin → project-delfin → user-canonical-compat."""
    out: list[Path] = [
        Path.home() / ".delfin" / "commands",
    ]
    if workspace is not None:
        out.append(Path(workspace) / ".delfin" / "commands")
    out.append(Path.home() / ".claude" / "commands")
    return out


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    m = _FRONTMATTER_RE.match(text or "")
    if not m:
        return {}, text or ""
    block, body = m.group(1), text[m.end():]
    meta: dict[str, str] = {}
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, _, v = line.partition(":")
        meta[k.strip()] = v.strip().strip("'").strip('"')
    return meta, body


def _load_one(path: Path) -> SlashTemplate | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    if not text.strip():
        return None
    meta, body = _parse_frontmatter(text)
    name = path.stem
    if not _NAME_RE.match(name):
        return None
    description = (meta.get("description") or "").strip()
    argument_hint = (meta.get("argument-hint") or meta.get("argument_hint") or "").strip()
    return SlashTemplate(
        name=name,
        description=description or _first_nonempty_line(body) or name,
        argument_hint=argument_hint,
        body=body.strip(),
        source=path,
    )


def _first_nonempty_line(text: str, *, max_len: int = 100) -> str:
    for line in (text or "").splitlines():
        s = line.strip().lstrip("#").strip()
        if s:
            return s[:max_len]
    return ""


def discover_commands(workspace: Path | str | None = None) -> list[SlashTemplate]:
    """Return every command markdown found across the search layers.

    Project-scoped files override user-global files with the same name
    (last-wins because we iterate in priority order, project after user).
    """
    by_name: dict[str, SlashTemplate] = {}
    for d in _search_dirs(Path(workspace) if workspace else None):
        if not d.is_dir():
            continue
        try:
            paths = sorted(d.glob("*.md"))
        except OSError:
            continue
        for p in paths:
            tpl = _load_one(p)
            if tpl is None:
                continue
            by_name[tpl.name] = tpl
    return sorted(by_name.values(), key=lambda t: t.name)


def get_command(
    name: str, workspace: Path | str | None = None,
) -> SlashTemplate | None:
    target = (name or "").strip().lstrip("/").lower()
    if not target:
        return None
    for tpl in discover_commands(workspace):
        if tpl.name.lower() == target:
            return tpl
    return None


def expand_template(template: str, args: str) -> str:
    """Substitute ``$ARGUMENTS`` / ``$@`` / ``$1`` … in ``template``.

    Unknown ``$<word>`` placeholders are left untouched so users can
    write ``"$0.05"`` inside a template without it being eaten. Numeric
    positional placeholders beyond available tokens become empty
    strings (matching POSIX shell semantics).
    """
    template = template or ""
    args = (args or "").strip()
    tokens = args.split() if args else []

    def _sub(m: re.Match[str]) -> str:
        ref = m.group(1)
        if ref in ("ARGUMENTS", "@"):
            return args
        if ref.isdigit():
            idx = int(ref)
            if idx == 0:
                return args
            return tokens[idx - 1] if idx <= len(tokens) else ""
        return m.group(0)

    pattern = re.compile(r"\$(\{[A-Za-z_@][A-Za-z0-9_@]*\}|[A-Za-z0-9_@]+)")

    def _norm(m: re.Match[str]) -> str:
        raw = m.group(1)
        ref = raw[1:-1] if raw.startswith("{") else raw
        if ref in ("ARGUMENTS", "@"):
            return args
        if ref.isdigit():
            idx = int(ref)
            if idx == 0:
                return args
            return tokens[idx - 1] if idx <= len(tokens) else ""
        return m.group(0)

    return pattern.sub(_norm, template)


def render_command(
    name_or_template: "str | SlashTemplate",
    args: str = "",
    workspace: Path | str | None = None,
) -> str:
    """Resolve a command name (or pre-loaded ``SlashTemplate``), expand
    ``$ARGUMENTS`` and positional placeholders, and return the final
    text ready to send to the agent as a user message."""
    tpl = (
        name_or_template
        if isinstance(name_or_template, SlashTemplate)
        else get_command(name_or_template, workspace)
    )
    if tpl is None:
        raise KeyError(f"unknown slash command: {name_or_template!r}")
    return expand_template(tpl.body, args)


__all__ = [
    "SlashTemplate",
    "discover_commands",
    "get_command",
    "expand_template",
    "render_command",
]
