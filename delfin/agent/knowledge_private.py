"""Personal lessons: what an agent learned, kept private until its person shares it.

Input of :func:`propose`: lesson text (``Situation:/Lesson:/Why:`` lines,
or one sentence that becomes the lesson), at least one evidence pointer
(``kind:ref``, :func:`knowledge.valid_evidence`), a source ("agent",
"user", "swarm"), optional title, topic, session id and workspace.
Output: a :class:`Lesson` with a small integer id, status ``private``.

Storage: ``~/.delfin/knowledge/private/lessons.json``, 0600 in a 0700
folder (``memory_store._private_dir``/``_atomic_write``). Nothing here
writes into a repository or runs a command: a lesson reaches git only
through ``knowledge_share``, behind the person's own confirmation, and
that module reads this store, never the other way round (pinned by a
test over this file's imports and calls).

Dedup: a new lesson whose text is a near-duplicate of a stored one with
the same polarity words (``memory_store._jaccard`` >= 0.72,
``_inverts_polarity`` False) merges its evidence into the stored one and
keeps that id; its share hint is not repeated.

Hint: :func:`pending_hints` returns one line per agent-proposed private
lesson that has not been hinted yet and marks it hinted, so the end of a
turn names a lesson once: "This lesson could help others: /share 3
(public)". The label reads the cached visibility of the target
(``repo_visibility``, no subprocess); unknown counts as public.
"""
from __future__ import annotations

import json
import re
import threading
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

PRIVATE = "private"
QUEUED = "queued"          # waiting in the admin review queue
PREPARED = "prepared"      # a branch was prepared for the person's PR
REJECTED = "rejected"
WITHDRAWN = "withdrawn"
STATUSES = (PRIVATE, QUEUED, PREPARED, REJECTED, WITHDRAWN)

_LOCK = threading.Lock()
_MERGE_SIMILARITY = 0.72


@dataclass
class Lesson:
    id: int
    title: str
    situation: str
    lesson: str
    why: str
    evidence: list = field(default_factory=list)
    topic: str = ""
    source: str = "agent"
    session_id: str = ""
    workspace: str = ""
    created: float = 0.0
    status: str = PRIVATE
    hinted: bool = False
    note: str = ""            # the reviewer's reason, shown to the author
    prepared_layer: str = ""  # what the person asked the agent to prepare
    shared_ref: str = ""      # queue id or branch once shared

    @property
    def text(self) -> str:
        return (f"{self.title}\nSituation: {self.situation}\n"
                f"Lesson: {self.lesson}\nWhy: {self.why}")


def store_path() -> Path:
    return Path.home() / ".delfin" / "knowledge" / "private" / "lessons.json"


def _read() -> dict:
    try:
        data = json.loads(store_path().read_text(encoding="utf-8"))
        if isinstance(data, dict) and isinstance(data.get("lessons"), dict):
            return data
    except (OSError, ValueError):
        pass
    return {"next_id": 1, "lessons": {}}


def _write(data: dict) -> None:
    from .memory_store import _atomic_write, _private_dir, _set_file_perms
    p = store_path()
    _private_dir(p.parent)
    _atomic_write(p, json.dumps(data, ensure_ascii=False, indent=1))
    _set_file_perms(p)


def _from_row(row: dict) -> Lesson:
    known = Lesson.__dataclass_fields__
    return Lesson(**{k: v for k, v in row.items() if k in known})


def _split(text: str) -> tuple[str, str, str]:
    fields = {"situation": "", "lesson": "", "why": ""}
    found = False
    for m in re.finditer(r"(?is)\b(Situation|Lesson|Why):\s*(.*?)"
                         r"(?=\b(?:Situation|Lesson|Why):|\Z)", text or ""):
        fields[m.group(1).lower()] = " ".join(m.group(2).split())
        found = True
    if not found:
        fields["lesson"] = " ".join(str(text or "").split())
    return fields["situation"], fields["lesson"], fields["why"]


def _title_of(lesson: str) -> str:
    first = re.split(r"(?<=[.!?])\s", lesson.strip(), maxsplit=1)[0]
    return first[:100].rstrip(".") or "Lesson"


def propose(text: str, *, evidence, source: str = "agent",
            title: str = "", topic: str = "", session_id: str = "",
            workspace: str = "", prepared_layer: str = "") -> Lesson:
    """Store a private lesson; raise ValueError without a valid pointer."""
    from .knowledge import valid_evidence
    pointers = [str(e).strip() for e in (evidence or []) if str(e).strip()]
    if not pointers:
        raise ValueError("a lesson needs evidence: a pointer such as "
                         "test:<node>, session:<id>, job:<id>, commit:<sha>")
    bad = [p for p in pointers if not valid_evidence(p)]
    if bad:
        raise ValueError(f"evidence must be a pointer (kind:ref), not "
                         f"content: {bad[0][:60]!r}")
    situation, lesson, why = _split(text)
    if not lesson:
        raise ValueError("a lesson needs the lesson itself")
    from .memory_store import _inverts_polarity, _jaccard
    with _LOCK:
        data = _read()
        for row in data["lessons"].values():
            old = _from_row(row)
            if old.status == WITHDRAWN:
                continue
            new_text = f"{situation} {lesson}"
            old_text = f"{old.situation} {old.lesson}"
            if (_jaccard(new_text, old_text) >= _MERGE_SIMILARITY
                    and not _inverts_polarity(lesson, old.lesson)):
                for p in pointers:
                    if p not in old.evidence:
                        old.evidence.append(p)
                if prepared_layer:
                    old.prepared_layer = prepared_layer
                data["lessons"][str(old.id)] = asdict(old)
                _write(data)
                return old
        les = Lesson(
            id=int(data.get("next_id", 1)),
            title=(title or "").strip()[:100] or _title_of(lesson),
            situation=situation, lesson=lesson, why=why,
            evidence=pointers, topic=(topic or "").strip(),
            source=str(source or "agent"), session_id=str(session_id or ""),
            workspace=str(workspace or ""), created=time.time(),
            prepared_layer=str(prepared_layer or ""))
        data["lessons"][str(les.id)] = asdict(les)
        data["next_id"] = les.id + 1
        _write(data)
        return les


def get(lesson_id) -> Optional[Lesson]:
    row = _read()["lessons"].get(str(lesson_id))
    return _from_row(row) if isinstance(row, dict) else None


def list_lessons() -> list[Lesson]:
    rows = _read()["lessons"].values()
    return sorted((_from_row(r) for r in rows if isinstance(r, dict)),
                  key=lambda x: x.id)


def _update(lesson_id, **changes) -> Lesson:
    with _LOCK:
        data = _read()
        row = data["lessons"].get(str(lesson_id))
        if not isinstance(row, dict):
            raise KeyError(f"no lesson {lesson_id}")
        les = _from_row(row)
        for k, v in changes.items():
            setattr(les, k, v)
        data["lessons"][str(les.id)] = asdict(les)
        _write(data)
        return les


def edit(lesson_id, **fields) -> Lesson:
    """The person's own rewording; the text is taken as given. ``evidence``
    replaces the pointers (space or comma separated, each kind:ref)."""
    allowed = {"title", "situation", "lesson", "why", "topic"}
    changes = {k: " ".join(str(v).split()) for k, v in fields.items()
               if k in allowed and v is not None}
    if fields.get("evidence") is not None:
        from .knowledge import valid_evidence
        pointers = [x for x in re.split(r"[\s,]+", str(fields["evidence"]))
                    if x]
        bad = [x for x in pointers if not valid_evidence(x)]
        if not pointers or bad:
            raise ValueError("evidence must be kind:ref pointers, e.g. "
                             "commit:<sha> or test:<node>")
        changes["evidence"] = pointers
    les = get(lesson_id)
    if les is None:
        raise KeyError(f"no lesson {lesson_id}")
    if les.status in (REJECTED, WITHDRAWN):
        changes["status"] = PRIVATE
    return _update(lesson_id, **changes)


def set_status(lesson_id, status: str, *, note: str = "",
               shared_ref: str = "") -> Lesson:
    if status not in STATUSES:
        raise ValueError(f"unknown status {status!r}")
    changes = {"status": status}
    if note:
        changes["note"] = note
    if shared_ref:
        changes["shared_ref"] = shared_ref
    return _update(lesson_id, **changes)


def delete(lesson_id) -> bool:
    with _LOCK:
        data = _read()
        if data["lessons"].pop(str(lesson_id), None) is None:
            return False
        _write(data)
        return True


def target_label(les: Lesson) -> str:
    """"public" unless the lesson's target is known to be private."""
    from . import repo_visibility as rv
    from .knowledge import global_dir, is_delfin_repo, repo_root_of
    ws = les.workspace
    if not ws or is_delfin_repo(ws) or les.prepared_layer == "global":
        root = global_dir().parents[3]
    else:
        root = repo_root_of(ws)
    return rv.label(rv.visibility(root, probe=False))


def pending_hints() -> list[str]:
    """One line per agent-proposed private lesson not hinted yet; marks
    each hinted, so a lesson is named once."""
    lines = []
    for les in list_lessons():
        if les.hinted or les.status != PRIVATE or les.source == "user":
            continue
        _update(les.id, hinted=True)
        lines.append(f"This lesson could help others: /share {les.id} "
                     f"({target_label(les)}) -- {les.title}")
    return lines
