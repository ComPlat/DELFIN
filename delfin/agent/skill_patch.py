"""Propose a new version of a skill as a patch (alt -> neu).

Hermes-style skill patching, with DELFIN's one non-negotiable
difference: the ACTIVE skill is never modified. A patch only ever
produces a *proposal* of the next version, via
``skill_proposals.propose(..., base_version=<old version>)``
(Paket 1); nothing activates automatically.

The replacement itself reuses the edit engine's ``apply_edit`` — one
unique exact match, otherwise an error that names the lines (and the
near misses). Nothing here reimplements that.

Public surface (contract with Paket 1 / api_client):

* ``propose_skill_patch(name, old, new, reason, evidence, ...)`` —
  the tool executor's entry point.
* ``archive_previous_version(skill_source, base_version)`` — called by
  ``skill_proposals.accept()`` when a patched skill (a proposal with a
  ``base_version``) becomes active: it moves the previous skill text to
  ``~/.delfin/skills/<name>/_versions/<version>.md`` and never
  overwrites an existing archive entry.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Protocol

from . import skills
from .edit_engine import apply_edit

__all__ = [
    "SkillPatchError",
    "Proposals",
    "propose_skill_patch",
    "archive_previous_version",
]


class SkillPatchError(Exception):
    """A patch could not be turned into a proposal.

    The message is user-facing (English) and names the reason: unknown
    skill, ambiguous or missing old_string, missing evidence.
    """


class Proposals(Protocol):
    """What ``delfin.agent.skill_proposals`` (Paket 1) provides.

    A structural stand-in so tests (and this module, until Paket 1
    lands on this branch) can pass a fake with the same shape.
    """

    def propose(self, name: str, text: str, *, evidence, source: str,
                base_version: str = ""): ...


_VERSION_RE = re.compile(r"^version:\s*(\d+)\s*$", re.MULTILINE)


def _normalise_name(name: str) -> str:
    """Accept what a model passes: '/tune' or 'tune.md' mean 'tune'.

    Same normalisation ``_execute_skill`` applies to the ``skill``
    tool — a leading slash and a copied extension are not a name.
    """
    name = (name or "").strip().lstrip("/").strip()
    if name.lower().endswith(".md"):
        name = name[:-3]
    return name.strip()


def _skill_text(skill) -> str:
    """The skill file's full text, front-matter included."""
    try:
        return Path(skill.source).read_text(encoding="utf-8")
    except OSError as exc:
        raise SkillPatchError(
            f"cannot read the active skill file '{skill.source}': {exc}"
        ) from exc


def _bump_version(text: str) -> tuple[str, str, str]:
    """Return (base_version, new_version, new_text).

    ``base_version`` is the front-matter version being patched ("1"
    when none is declared — the first patch of an unversioned skill
    treats what is on disk as v1), ``new_version`` is that + 1, and
    ``new_text`` carries the bumped ``version:`` line.
    """
    m = _VERSION_RE.search(text)
    if m:
        base = m.group(1)
    else:
        base = "1"
        text = "---\n" + text.lstrip("-\n") if text.startswith("---") \
            else text
        # Insert a version into the front-matter block, or add one.
        if text.startswith("---"):
            end = text.find("\n---", 3)
            if end >= 0:
                text = (text[:end] + f"\nversion: {base}"
                        + text[end:])
            else:
                text = f"---\nversion: {base}\n{text}"
        else:
            text = f"---\nversion: {base}\n---\n\n{text}"
    new = str(int(base) + 1)
    new_text = _VERSION_RE.sub(lambda _: f"version: {new}", text, count=1)
    if not _VERSION_RE.search(new_text) or \
            _VERSION_RE.search(new_text).group(1) != new:  # pragma: no cover
        raise SkillPatchError("could not bump the skill version")
    return base, new, new_text


def _load_proposals():
    """The real proposal store (Paket 1), imported lazily.

    Importing at module level would fail while Paket 1 is not on this
    branch; the lazy import keeps this module importable and lets the
    caller inject a stand-in instead.
    """
    try:
        from . import skill_proposals
    except ImportError as exc:  # pragma: no cover - depends on Paket 1
        raise SkillPatchError(
            "skill_proposals (Paket 1) is not available: "
            f"{exc}. The patch tool needs the proposal store."
        ) from exc
    return skill_proposals


def propose_skill_patch(
    name: str,
    old: str,
    new: str,
    reason: str,
    evidence,
    *,
    workspace=None,
    proposals: Optional[Proposals] = None,
) -> dict:
    """Propose the next version of an active skill as a patch.

    Never modifies the active skill: the patched text goes to the
    proposal store as a proposal with ``base_version`` = the version
    on disk. Activation is a human's ``accept`` (Paket 1).

    Args:
        name: skill name ('/tune' and 'tune.md' normalise to 'tune').
        old, new: the replacement, exactly one unique ``old`` match in
            the skill text (``edit_engine.apply_edit`` semantics).
        reason: why the patch improves the skill (kept on record).
        evidence: non-empty list of evidence dicts (kind/ref/detail).
            Without evidence there is no proposal — an unproven
            success must not become a rule.
        workspace: workspace for skill discovery.
        proposals: proposal store to use; defaults to the real
            ``skill_proposals`` module (Paket 1), imported lazily.

    Returns:
        ``{"name", "base_version", "new_version", "status"}`` of the
        proposal that was created.

    Raises:
        SkillPatchError: unknown skill, no evidence, or the patch does
            not apply (ambiguous / missing ``old`` — the message names
            the matching and near-miss lines).
    """
    name = _normalise_name(name)
    if not name:
        raise SkillPatchError("skill name must be non-empty")
    evidence = list(evidence or [])
    if not evidence:
        raise SkillPatchError(
            "no evidence given: a patch without evidence cannot become "
            "a proposal. Cite a green test, a documented run or a "
            "verified recipe.")
    skill = skills.get_skill(name, workspace)
    if skill is None:
        raise SkillPatchError(
            f"no active skill named '{name}' — a patch can only revise "
            "an existing skill, it cannot create one. Have it written "
            "as a proposal first (skill_learning / Paket 2).")
    text = _skill_text(skill)
    result = apply_edit(text, old or "", new or "")
    if not result.applied:
        detail = result.error or "the old string did not apply"
        lines = ""
        if result.match_lines:
            lines = (" Matches on lines: "
                     + ", ".join(str(n) for n in result.match_lines) + ".")
        elif result.near_misses:
            lines = (" Near matches: "
                     + "; ".join(f"line {nm.line} ({nm.reason})"
                                 for nm in result.near_misses) + ".")
        raise SkillPatchError(
            f"patch to '{name}' did not apply: {detail}.{lines}")
    base, new_version, new_text = _bump_version(result.new_text)
    store = proposals if proposals is not None else _load_proposals()
    proposal = store.propose(
        name, new_text, evidence=evidence, source="skill_patch",
        base_version=base)
    return {
        "name": name,
        "base_version": base,
        "new_version": new_version,
        "status": getattr(proposal, "status", "pending"),
    }


def archive_previous_version(skill_source, base_version: str) -> Path:
    """Move the previous skill text under ``_versions/`` — never overwrite.

    Called by ``skill_proposals.accept()`` (Paket 1) when a proposal
    with a ``base_version`` becomes the active skill: the version that
    was active until then is preserved verbatim at

        ``<skill dir>/_versions/<base_version>.md``

    (a flat ``<name>.md`` skill archives next to itself, a folder
    skill inside its own folder). An existing archive file for that
    version is an error — accepting must never silently overwrite
    history; the caller decides whether a different version number is
    in order.

    Args:
        skill_source: path of the skill file that was active until now.
        base_version: its front-matter version; non-empty.

    Returns:
        The path of the archived copy.

    Raises:
        SkillPatchError: no version given, unreadable source, or an
            archive entry for that version already exists.
    """
    base_version = (base_version or "").strip()
    if not base_version:
        raise SkillPatchError(
            "cannot archive a skill version: no version given. A skill "
            "without a version cannot be replaced by a patched one.")
    src = Path(skill_source)
    try:
        text = src.read_text(encoding="utf-8")
    except OSError as exc:
        raise SkillPatchError(
            f"cannot read '{src}' to archive it: {exc}") from exc
    archive_dir = src.parent / "_versions"
    where = archive_dir / f"{base_version}.md"
    if where.exists():
        raise SkillPatchError(
            f"archive entry '{where}' already exists — the previous "
            "version is not overwritten. Rename the proposal's "
            "version instead.")
    archive_dir.mkdir(parents=True, exist_ok=True)
    where.write_text(text, encoding="utf-8")
    return where
