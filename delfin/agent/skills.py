"""User-invocable skills (.delfin-native).

Skills are short, parameterised playbooks the user can invoke from
the chat input via ``/<skill-name>`` or that the agent can invoke via
the ``Skill`` tool. They live as Markdown files with a YAML
front-matter::

    ---
    name: security-review
    description: Run a security audit of pending changes
    ---

    # Security review

    1. Inspect every changed file for hardcoded secrets …

Search order (later wins on name collision):

  1. ``~/.delfin/skills/<name>/SKILL.md``        — user-global
  2. ``~/.delfin/skills/<name>.md``              — user-global flat
  3. ``<workspace>/.delfin/skills/<name>/SKILL.md`` — project-scoped
  4. ``<workspace>/.delfin/skills/<name>.md``    — project-scoped flat

The loader is intentionally permissive: missing front-matter is
tolerated (the file's first heading becomes the description). Only
file IO errors silently skip a candidate.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class Skill:
    name: str
    description: str
    body: str
    source: Path
    # Where this playbook applies, from the ``domains:`` front-matter key.
    # Empty means everywhere -- a user's own skill keeps working without
    # knowing that domains exist, and nothing can lose a skill it has today.
    domains: tuple[str, ...] = ()

    def applies_to(self, domain: str) -> bool:
        """Whether this skill belongs in a session working in ``domain``.

        A skill is offered where its TOOLS are. The curated set splits
        cleanly: the chemistry playbooks call ``search_docs`` and
        ``search_calcs``, the administrative ones call ``read_document``
        and ``compare_tables``, and neither role can reach the other's
        surface -- so offering both to both meant every session was told
        about playbooks whose first step it could not execute.
        """
        if not self.domains:
            return True
        d = (domain or "").strip().lower()
        if not d:
            return True
        return d in self.domains

    def to_summary(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "source": str(self.source),
            "domains": list(self.domains),
        }


_SKILL_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_\-]{0,63}$")


def _parse_frontmatter(text: str) -> tuple[dict, str]:
    """Extract YAML-style front-matter; return (meta, remaining_body)."""
    if not text.startswith("---"):
        return {}, text
    parts = text.split("\n", 1)
    if len(parts) < 2:
        return {}, text
    rest = parts[1]
    end = rest.find("\n---")
    if end < 0:
        return {}, text
    block = rest[:end]
    body = rest[end + 4 :].lstrip("\n")
    meta: dict[str, str] = {}
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            continue
        k, _, v = line.partition(":")
        meta[k.strip()] = v.strip().strip("\"").strip("'")
    return meta, body


def _parse_domains(raw: object) -> tuple[str, ...]:
    """Read a ``domains:`` value as a tuple of lowercase names.

    Accepts what people actually write: ``office``, ``office, code``, or
    the inline YAML list ``[office, code]``. Anything unparseable reads as
    "no domain given", which means universal -- a typo must not be able to
    hide a skill from every session at once.
    """
    if not isinstance(raw, str):
        return ()
    text = raw.strip().strip("[]")
    parts = [p.strip().strip("\"'").lower()
             for p in text.replace(";", ",").split(",")]
    return tuple(p for p in parts if p)


def _first_heading(body: str) -> str:
    for line in body.splitlines():
        line = line.strip()
        if line.startswith("#"):
            return line.lstrip("#").strip()
    return ""


# Built-in "curated set" skills shipped with the package.
_PACK_SKILLS_DIR = Path(__file__).resolve().parent / "pack" / "skills"


def _skill_dirs(workspace: Path | None) -> list[Path]:
    # Pack (built-in) skills first — the lowest-precedence curated base — then
    # user-global, then project-scoped. Later dirs override earlier ones on a
    # name collision (project > user > pack). Without the pack dir the 9
    # shipped chemistry skills were undiscoverable by the `skill` tool.
    out: list[Path] = [_PACK_SKILLS_DIR, Path.home() / ".delfin" / "skills"]
    if workspace is not None:
        out.append(Path(workspace) / ".delfin" / "skills")
    return out


def _candidates_in_dir(d: Path) -> Iterable[Path]:
    if not d.is_dir():
        return ()
    found: list[Path] = []
    try:
        for child in d.iterdir():
            if child.is_dir():
                p = child / "SKILL.md"
                if p.is_file():
                    found.append(p)
            elif child.suffix == ".md" and child.name != "SKILL.md":
                found.append(child)
    except OSError:
        pass
    return found


def _load_one(path: Path) -> Skill | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    if not text.strip():
        return None
    meta, body = _parse_frontmatter(text)
    name = (meta.get("name") or "").strip()
    if not name:
        # Folder-style: parent dir name; flat-style: stem
        if path.name == "SKILL.md":
            name = path.parent.name
        else:
            name = path.stem
    if not _SKILL_NAME_RE.match(name):
        return None
    description = (meta.get("description") or "").strip()
    if not description:
        description = _first_heading(body)
    return Skill(
        name=name,
        description=description,
        body=body.strip() or text.strip(),
        source=path,
        domains=_parse_domains(meta.get("domains")),
    )


def discover_skills(
    workspace: Path | str | None = None, domain: str = "",
) -> list[Skill]:
    """Return all skills found in user + project directories.

    On name collisions the later (project-scoped) skill wins, mirroring
    the canonical "project overrides user" semantics.

    With ``domain`` set, only the skills that apply there are returned.
    The default returns everything, so a caller that never learned about
    domains sees exactly what it saw before.
    """
    by_name: dict[str, Skill] = {}
    for d in _skill_dirs(Path(workspace) if workspace else None):
        for path in _candidates_in_dir(d):
            sk = _load_one(path)
            if sk is None:
                continue
            by_name[sk.name] = sk
    found = sorted(by_name.values(), key=lambda s: s.name)
    if not domain:
        return found
    return [s for s in found if s.applies_to(domain)]


def get_skill(
    name: str, workspace: Path | str | None = None, domain: str = "",
) -> Skill | None:
    """Look up a single skill by name.

    ``domain`` filters the same way as ``discover_skills``. Callers that
    need to tell "does not exist" apart from "does not apply here" should
    look it up WITHOUT a domain and ask ``applies_to`` themselves — the
    two cases deserve different messages.
    """
    for sk in discover_skills(workspace, domain=domain):
        if sk.name == name:
            return sk
    return None


def render_skill_invocation(skill: Skill, args: str = "") -> str:
    """Render the skill body into a ready-to-inject prompt block."""
    args = (args or "").strip()
    header = f"--- Skill: {skill.name} ---"
    arg_line = f"\nArguments: {args}\n" if args else "\n"
    return f"{header}{arg_line}\n{skill.body}\n--- end skill ---\n"


__all__ = [
    "Skill",
    "discover_skills",
    "get_skill",
    "render_skill_invocation",
]
