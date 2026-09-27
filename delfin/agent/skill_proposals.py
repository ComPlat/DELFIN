"""Skill proposals: a self-written skill is only a PROPOSAL.

The Hermes lesson (arXiv 2608.12851): an agent that activates its own
skills turns unsafe successes into permanent rules. So a proposal lives
in ``~/.delfin/skills/_proposals`` where ``discover_skills`` never looks
(skills.py only loads direct ``<name>/SKILL.md`` children of the skill
dirs), it is blocked when the safety check finds anything, and only a
human ``accept`` moves it into ``~/.delfin/skills/<name>/``.

Never overwritten, never deleted: a name conflict gets a new name, a
rejection moves the folder to ``rejected/``.
"""
from __future__ import annotations

import dataclasses
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from .memory_store import _atomic_write, _private_dir, _set_file_perms

def _proposals_dir() -> Path:
    # Resolved per call, not at import time: tests (and any redirected
    # HOME) must be able to move the store after this module is loaded.
    return Path.home() / ".delfin" / "skills" / "_proposals"


def __getattr__(name):
    # PEP 562: ``PROPOSALS_DIR`` reads the CURRENT home, like every
    # other path in memory_store. A module-level constant would pin the
    # first HOME it ever saw.
    if name == "PROPOSALS_DIR":
        return _proposals_dir()
    raise AttributeError(name)


_STATUS_PENDING = "pending"
_STATUS_BLOCKED = "blocked"
_STATUS_ACCEPTED = "accepted"
_STATUS_REJECTED = "rejected"


@dataclasses.dataclass
class Evidence:
    """The proof a proposal rests on. No evidence, no proposal."""

    kind: str  # "test" | "calc" | "job" | "recipe"
    ref: str   # test node id, calc folder, job id, ...
    detail: str = ""
    verified_at: str = ""


@dataclasses.dataclass
class Proposal:
    name: str
    text: str
    evidence: list[Evidence]
    source: str
    status: str  # pending | blocked | accepted | rejected
    findings: list[str]
    created: str
    base_version: str = ""


# --------------------------------------------------------------------------
# helpers


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _safety_findings(text: str) -> list[str]:
    """Ask skill_safety (package 3); a check that cannot run blocks.

    An absent or crashing module is never a silent pass: the unchecked
    proposal goes to status "blocked" with a finding the human can read.
    Until package 3 lands, that is every proposal -- the correct default
    is "not verified", not "verified by absence of a checker".
    """
    try:
        from . import skill_safety  # noqa: PLC0415
    except ImportError as exc:
        return [f"safety check unavailable: {exc}"]
    try:
        return [str(f) for f in skill_safety.check(text)]
    except Exception as exc:  # the check itself is the guard: it may not
        # crash its caller, and it may not pass what it could not read.
        return [f"safety check failed: {exc!r}"]


def _prop_dir(name: str, *, root: Path | None = None) -> Path:
    return (root or _proposals_dir()) / name


def _free_name(name: str, root: Path | None = None,
               extra: Path | None = None) -> str:
    """A name whose folder is free; conflicts get ``-2``, ``-3``, ..."""
    base = root or _proposals_dir()
    taken = lambda n: _prop_dir(n, root=base).exists() or (
        extra is not None and (extra / n).exists())
    if not taken(name):
        return name
    i = 2
    while taken(f"{name}-{i}"):
        i += 1
    return f"{name}-{i}"


def _write_proposal(proposal: Proposal, dirpath: Path) -> None:
    """Persist a proposal as ``SKILL.md`` + ``proposal.json`` (0600)."""
    _private_dir(dirpath)
    _atomic_write(dirpath / "SKILL.md", proposal.text)
    data = dataclasses.asdict(proposal)
    _atomic_write(dirpath / "proposal.json", json.dumps(data, ensure_ascii=False))
    _set_file_perms(dirpath / "proposal.json")


def _read_proposal(dirpath: Path) -> Proposal:
    data = json.loads((dirpath / "proposal.json").read_text(encoding="utf-8"))
    data["evidence"] = [Evidence(**e) for e in data.get("evidence", [])]
    return Proposal(**data)


def _iter_proposal_dirs(root: Path | None = None) -> list[Path]:
    base = root or _proposals_dir()
    if not base.is_dir():
        return []
    out = []
    for child in sorted(base.iterdir()):
        if child.is_dir() and (child / "proposal.json").is_file():
            out.append(child)
    return out


# --------------------------------------------------------------------------
# public API


def propose(name: str, text: str, *, evidence, source: str,
            base_version: str = "") -> Proposal:
    """Record a new skill proposal. Invisible until accepted by a human.

    Raises ValueError without evidence. A name conflict gets a new name
    (never an overwrite). Safety findings move the proposal to status
    "blocked" and are stored for the human to read.
    """
    evs = [e if isinstance(e, Evidence) else Evidence(**e) for e in (evidence or [])]
    if not evs:
        raise ValueError("a skill proposal needs evidence: no evidence, no rule")
    name = (name or "").strip()
    if not name or "/" in name or name.startswith("."):
        raise ValueError(f"invalid proposal name: {name!r}")
    findings = _safety_findings(text)
    proposal = Proposal(
        name=name,
        text=text,
        evidence=evs,
        source=source,
        status=_STATUS_BLOCKED if findings else _STATUS_PENDING,
        findings=findings,
        created=_now(),
        base_version=base_version,
    )
    skills_root = _proposals_dir().parent
    # a pending proposal must not shadow a skill that already exists, and
    # two proposals must not collide: both get a fresh name, nothing is
    # overwritten.
    proposal.name = _free_name(name, extra=skills_root)
    _write_proposal(proposal, _prop_dir(proposal.name))
    return proposal


def list_proposals(status: str | None = None) -> list[Proposal]:
    """All proposals on disk, optionally filtered by status.

    Includes ``rejected/``: a rejection is a decision to be findable,
    not a deletion. Accepted proposals are NOT listed -- they are live
    skills now, and the live skills tree is their home.
    """
    out = []
    roots = [_proposals_dir(), _proposals_dir() / "rejected"]
    for base in roots:
        for d in _iter_proposal_dirs(base):
            try:
                p = _read_proposal(d)
            except (OSError, ValueError, KeyError, TypeError):
                continue
            if status is None or p.status == status:
                out.append(p)
    return sorted(out, key=lambda p: (p.status, p.name))


def get_proposal(name: str) -> Proposal | None:
    d = _prop_dir(name)
    if not (d / "proposal.json").is_file():
        return None
    try:
        return _read_proposal(d)
    except (OSError, ValueError, KeyError, TypeError):
        return None


def accept(name: str, *, by: str) -> Path:
    """Move a pending proposal into the live skills directory.

    Only status "pending" is acceptable: "blocked" means the safety check
    found something and a human must read it, not wave it through here.
    A name already taken by a live skill gets a new name.
    """
    d = _prop_dir(name)
    if not (d / "proposal.json").is_file():
        raise ValueError(f"no proposal named {name!r}")
    proposal = _read_proposal(d)
    if proposal.status == _STATUS_BLOCKED:
        raise ValueError(
            f"proposal {name!r} is blocked by the safety check; "
            "review its findings before accepting")
    if proposal.status != _STATUS_PENDING:
        raise ValueError(f"proposal {name!r} is already {proposal.status}")

    skills_root = _proposals_dir().parent
    # Only the live skills tree counts here: the proposal's own folder in
    # _proposals is being moved away, so counting it would rename every
    # accepted proposal to "<name>-2" unconditionally.
    final_name = _free_name(name, root=skills_root)
    target = skills_root / final_name
    _private_dir(skills_root)
    os.replace(d, target)
    proposal.status = _STATUS_ACCEPTED
    proposal.name = final_name
    # keep the acceptance traceable inside the live skill folder
    _atomic_write(target / "proposal.json",
                  json.dumps(dataclasses.asdict(proposal), ensure_ascii=False))
    _set_file_perms(target / "proposal.json")
    return target / "SKILL.md" if (target / "SKILL.md").exists() else target


def reject(name: str, *, reason: str, by: str) -> Path:
    """Move a proposal to ``rejected/`` with the reason. Nothing is deleted."""
    if not reason or not reason.strip():
        raise ValueError("a rejection needs a reason")
    d = _prop_dir(name)
    if not (d / "proposal.json").is_file():
        raise ValueError(f"no proposal named {name!r}")
    proposal = _read_proposal(d)
    if proposal.status == _STATUS_ACCEPTED:
        raise ValueError(f"proposal {name!r} is already accepted; "
                         "rejecting a live skill is not this call")
    rejected_root = _proposals_dir() / "rejected"
    _private_dir(rejected_root)
    final = _free_name(name, root=rejected_root)
    target = rejected_root / final
    proposal.status = _STATUS_REJECTED
    proposal.findings = list(proposal.findings)
    os.replace(d, target)
    _write_proposal(proposal, target)
    (target / "REJECTED.txt").write_text(
        f"reason: {reason.strip()}\nby: {by}\n", encoding="utf-8")
    _set_file_perms(target / "REJECTED.txt")
    return target

