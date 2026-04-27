"""Skill registry for the DELFIN agent.

A *skill* is a small, named markdown prompt the user can invoke from the
dashboard via ``/skill <name>``.  The body of the skill is sent to the
agent as a user message — the same effect as if the user had typed the
prompt themselves, but kept reusable and discoverable.

Skill files live under ``delfin/agent/pack/skills/``.  Each ``<name>.md``
file follows a tiny convention so we don't need a YAML parser dependency
for the body:

    # <Title>
    > <Optional one-line description>

    <Skill body — full prompt sent to the agent>

The leading ``#`` is the title (used in the palette).  A line beginning
with ``>`` immediately after is treated as the short summary.  Anything
else is the skill body.

If a skill file is malformed (missing title), it is silently ignored.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# Skill record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Skill:
    """One named skill loaded from disk."""

    name: str          # filename stem, e.g. ``tune-control``
    title: str         # user-visible title (from the leading ``# ...``)
    description: str   # one-line summary (from a leading ``> ...`` line)
    body: str          # full prompt text sent to the agent
    source_path: str   # absolute path of the MD file (for debugging)

    @property
    def slash_command(self) -> str:
        """The full slash-command including the ``/skill`` prefix."""
        return f"/skill {self.name}"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

_DEFAULT_DIR = Path(__file__).resolve().parent / "pack" / "skills"


def parse_skill_text(name: str, text: str, source_path: str = "") -> Skill | None:
    """Parse skill markdown into a ``Skill`` record.

    Returns ``None`` if the file does not start with a ``# Title`` line.
    """
    if not text or not text.strip():
        return None

    lines = text.splitlines()
    title = ""
    description = ""
    body_start = 0

    # First non-blank line must be the title (# Heading)
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            title = stripped.lstrip("#").strip()
            body_start = i + 1
        break

    if not title:
        return None

    # Optional description on a leading ``>`` line
    while body_start < len(lines):
        candidate = lines[body_start].strip()
        if not candidate:
            body_start += 1
            continue
        if candidate.startswith(">"):
            description = candidate.lstrip(">").strip()
            body_start += 1
        break

    body = "\n".join(lines[body_start:]).strip()
    return Skill(
        name=name,
        title=title,
        description=description,
        body=body,
        source_path=source_path,
    )


def discover_skills(skills_dir: Path | None = None) -> list[Skill]:
    """Return all valid skills found in ``skills_dir``.

    Sorted alphabetically by name for deterministic UI ordering.
    """
    base = skills_dir if skills_dir is not None else _DEFAULT_DIR
    if not base.is_dir():
        return []

    out: list[Skill] = []
    for path in sorted(base.glob("*.md")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        skill = parse_skill_text(path.stem, text, source_path=str(path))
        if skill is not None:
            out.append(skill)
    return out


def find_skill(name: str, skills_dir: Path | None = None) -> Skill | None:
    """Find a single skill by name (case-insensitive)."""
    target = (name or "").strip().lower()
    if not target:
        return None
    for skill in discover_skills(skills_dir):
        if skill.name.lower() == target:
            return skill
    return None


def format_skill_message(skill: Skill) -> str:
    """Build the user-message expansion for a skill invocation.

    Adds a tiny header so the agent knows which skill is active and the
    user's intent, without obscuring the original skill body.
    """
    header = f"[Skill: {skill.title}]"
    return f"{header}\n\n{skill.body}".strip()


__all__ = [
    "Skill",
    "parse_skill_text",
    "discover_skills",
    "find_skill",
    "format_skill_message",
]
