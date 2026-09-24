"""Tidy a typed-memory store by proposing, not by sweeping.

Input: a store directory. Output: a proposal — near-duplicate pairs to
fold together, model-written entries nobody has recalled to retire, and
the counts before and after. ``propose`` writes nothing. ``apply``
performs the proposal and never deletes: a retired file is moved to
``<store>/retired/``.

Why a proposal. A store is what the agent learned, and an unsupervised
pass over it can drop a fact without anyone noticing for weeks. The
opposite failure was measured on 2026-09-24: 358 of 371 notes had been
unreachable for months, because their store was keyed to paths that no
longer existed, and that went unnoticed as well. Both directions are
silent, so this one is made visible instead of automatic.

Reused rather than reimplemented: ``_jaccard`` and
``_merge_similarity_threshold`` decide what counts as a near-duplicate,
``list_typed_memories`` supplies the records with their usage and decay
metadata, and ``_AGENT_MEMORY_MAX_AGE_DAYS`` defines disuse. What is new
is when they are applied: merging happens today only as a memory is
WRITTEN, so look-alikes already in a store stay there, and
``prune_memories`` deletes on its own schedule without showing anyone
what it took.

Two rules the proposal never breaks. Memories of different types are not
merged -- the type is a classification, and folding across it loses that.
Memories the USER wrote are never retired: disuse is not a reason to drop
what somebody put there deliberately.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from delfin.agent.memory_store import (
    _AGENT_MEMORY_MAX_AGE_DAYS,
    _jaccard,
    _merge_similarity_threshold,
)

#: Written into a folded file between the two bodies, so a later reader
#: can see that two notes became one and where the seam is.
_SEAM = "\n\n---\n\n"

_DAY = 86400


@dataclass
class Note:
    """One memory file, as the proposal sees it."""
    path: Path
    name: str
    type: str
    source: str
    updated_at: int
    body: str


@dataclass
class Merge:
    keep: Note
    drop: Note
    similarity: float


@dataclass
class Proposal:
    store: Path
    merges: List[Merge] = field(default_factory=list)
    retire: List[Note] = field(default_factory=list)
    before: int = 0

    @property
    def after(self) -> int:
        return max(0, self.before - len(self.merges) - len(self.retire))

    def select(self, names) -> "Proposal":
        """A copy holding only the items the caller named.

        All-or-nothing takes every decision away at once, and a tidy pass
        is a list of separate judgements. A merge is named by EITHER of
        its two notes, because that is how a reader refers to it. A name
        that matches nothing selects nothing rather than everything: the
        failure of a typo has to be doing less, never more.
        """
        wanted = {str(n).strip() for n in (names or ()) if str(n).strip()}
        picked = Proposal(store=self.store, before=self.before)
        picked.merges = [m for m in self.merges
                         if m.keep.name in wanted or m.drop.name in wanted]
        picked.retire = [r for r in self.retire if r.name in wanted]
        return picked

    def render(self) -> str:
        """The report a person reads before deciding."""
        lines: List[str] = []
        if self.merges:
            lines.append(f"{len(self.merges)} to merge")
            for m in self.merges:
                lines.append(f"    {m.keep.name}  +  {m.drop.name}"
                             f"   (similarity {m.similarity:.2f})")
                lines.append(f"      keep: {_one_line(m.keep.body)}")
                lines.append(f"      fold: {_one_line(m.drop.body)}")
        if self.retire:
            lines.append(f"{len(self.retire)} to retire "
                         f"(model-written, not recalled in "
                         f"{_AGENT_MEMORY_MAX_AGE_DAYS} days)")
            for n in self.retire:
                lines.append(f"    {n.name}   {_one_line(n.body)}")
        if not lines:
            lines.append("nothing to tidy")
        lines.append("")
        lines.append(f"{self.before} notes before · {self.after} after")
        if self.merges or self.retire:
            lines.append("nothing changed — run again with --apply")
        return "\n".join(lines)


def _one_line(text: str, width: int = 68) -> str:
    flat = " ".join(str(text or "").split())
    return flat if len(flat) <= width else flat[:width - 1] + "…"


def _read(path: Path) -> Optional[Note]:
    """One note, or None when the file cannot be read as one.

    A damaged file is skipped rather than guessed at: a tidy pass that
    rewrites something it did not understand is the failure this whole
    module is shaped to avoid.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    if not text.startswith("---"):
        return None
    try:
        head, body = text.split("---", 2)[1], text.split("---", 2)[2]
    except IndexError:
        return None
    meta = {}
    for line in head.splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip()
    try:
        updated = int(meta.get("updated_at") or 0)
    except ValueError:
        updated = 0
    name = meta.get("name") or path.stem
    mtype = meta.get("type") or _type_from_name(path.stem)
    if not mtype:
        return None
    return Note(path=path, name=name, type=mtype,
                source=(meta.get("source") or "user").strip().lower(),
                updated_at=updated, body=body.strip())


def _type_from_name(stem: str) -> str:
    for t in ("feedback", "project", "reference", "user"):
        if stem.startswith(t + "_"):
            return t
    return ""


def hint(store: Path, *, now: Optional[int] = None) -> str:
    """One line for the session start, or "" when there is nothing to say.

    Measured 2026-09-24: proposing over 1000 notes costs 43 ms, and over
    the real stores (189 notes) it found one merge and no retirements. So
    it is cheap enough to ask on every start and quiet enough not to nag
    — which is the point, because a command nobody is told about is a
    command nobody runs.

    Never raises: a hint that can take the session start with it is not a
    hint.
    """
    try:
        p = propose(store, now=now)
        n = len(p.merges) + len(p.retire)
        if not n:
            return ""
        what = "note" if n == 1 else "notes"
        return (f"memory     {n} {what} could be tidied "
                f"(/tidy shows what, and changes nothing)")
    except Exception:
        return ""


def propose(store: Path, *, now: Optional[int] = None) -> Proposal:
    """What tidying this store would do. Writes nothing."""
    store = Path(store)
    moment = int(now if now is not None else time.time())
    notes: List[Note] = []
    try:
        files = sorted(store.glob("*.md"))
    except OSError:
        files = []
    for path in files:
        note = _read(path)
        if note is not None:
            notes.append(note)

    proposal = Proposal(store=store, before=len(notes))
    threshold = _merge_similarity_threshold()

    # Retire first, so a note on its way out is not also proposed for a
    # merge -- two actions on one file is how an apply loses a body.
    retiring = set()
    horizon = moment - _AGENT_MEMORY_MAX_AGE_DAYS * _DAY
    for note in notes:
        if note.source == "agent" and note.updated_at and (
                note.updated_at < horizon):
            proposal.retire.append(note)
            retiring.add(note.path)

    remaining = [n for n in notes if n.path not in retiring]
    paired = set()
    for i, first in enumerate(remaining):
        if first.path in paired:
            continue
        for second in remaining[i + 1:]:
            if second.path in paired or second.type != first.type:
                continue
            score = _jaccard(first.body, second.body)
            if score < threshold:
                continue
            older, newer = ((first, second)
                            if first.updated_at <= second.updated_at
                            else (second, first))
            proposal.merges.append(
                Merge(keep=older, drop=newer, similarity=score))
            paired.add(first.path)
            paired.add(second.path)
            break
    return proposal


def apply(proposal: Proposal) -> dict:
    """Carry out *proposal*. Returns what was done. Deletes nothing."""
    done = {"merged": 0, "retired": 0, "failed": 0}
    for merge in proposal.merges:
        try:
            keep_text = merge.keep.path.read_text(encoding="utf-8")
            if merge.drop.body and merge.drop.body not in keep_text:
                merge.keep.path.write_text(
                    keep_text.rstrip("\n") + _SEAM + merge.drop.body + "\n",
                    encoding="utf-8")
            _retire_file(proposal.store, merge.drop.path)
            done["merged"] += 1
        except OSError:
            done["failed"] += 1
    for note in proposal.retire:
        try:
            _retire_file(proposal.store, note.path)
            done["retired"] += 1
        except OSError:
            done["failed"] += 1
    return done


def _retire_file(store: Path, path: Path) -> None:
    """Move a file out of the store, dated. Never unlink.

    The store shrinks, the text stays readable, and a wrong call costs a
    move rather than a fact.
    """
    room = Path(store) / "retired"
    room.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%d", time.localtime())
    target = room / f"{stamp}_{path.name}"
    n = 2
    while target.exists():
        target = room / f"{stamp}_{path.stem}-{n}{path.suffix}"
        n += 1
    path.replace(target)


__all__ = ["Merge", "Note", "Proposal", "apply", "hint", "propose"]
