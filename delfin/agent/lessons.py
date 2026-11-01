"""Reviewed lessons carried across tasks.

A maintainer curates Markdown lesson files under a configured directory
(``agent.lessons_dir``, default unset = off). Each bullet is one lesson,
carrying a short scope tag; for a new task the loader picks the few whose
scope matches the task text, at most ``max_lessons`` within a character
budget, in deterministic order, with no model call for selection.

Lesson source format (a catalog a maintainer keeps in Git):

    ## orca-cards
    - [orca-cards,dispersion] Match the dispersion correction to the functional.
    - [orca-cards] The Hessian must be read from the same geometry.
    ## sandbox
    - [sandbox] The scratch TMPDIR must resolve inside the bwrap layer.

A bullet's scope is the ``[tag1,tag2]`` prefix (comma-separated); bullets
under a ``## <tag>`` section heading inherit that tag as well, so a
heading alone is a valid way to scope a group.

The explicit ``[tag]`` prefix and the behaviour of untrusted text are the
two reasons lessons cannot be treated as instructions:

- selection is DETERMINISTIC and stops at a hard budget — no model call,
  no ranking that could be steered by earlier lessons;
- each selected lesson is rendered as a bare data line and wrapped (by the
  caller) inside an ``untrusted`` notes fence; a lesson text that carries an
  instruction channel ("ignore all previous instructions", a role marker,
  a standing consent) is REFUSED via the same shape guard the memory store
  uses, so it never reaches the prompt;
- control characters are stripped;
- only files inside the configured directory are read; a ``source`` that
  walks out of it (``..``) is dropped.

No network access, no writes, no model call — this is a read-only selector.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from delfin.agent.repl_render import strip_control as _strip_control
from delfin.agent.untrusted import instruction_shapes

#: A lesson in a maintainer's catalog. ``text`` is the lesson body, ``scope``
#: its scope tags (a set; order is not significant), ``source`` the catalog
#: file name (basename only, never a path) and ``line`` the 1-based line of
#: the bullet in that file.
@dataclass(frozen=True)
class Lesson:
    text: str
    scope: frozenset[str] = field(default_factory=frozenset)
    source: str = ""
    line: int = 0


# ---------------------------------------------------------------------------
# Markdown parsing
# ---------------------------------------------------------------------------

_TAG_PREFIX_RE = re.compile(r"^\s*\[([^\]]+)\]\s*(.*)$")


def _section_tag(line: str) -> str | None:
    """The tag a ``## <tag>`` heading declares, or None."""
    m = re.match(r"^\s*##\s+([^\s#]+)", line)
    return m.group(1) if m else None


def parse_lessons_file(path: Path | str) -> list[Lesson]:
    """Parse one Markdown catalog file into lessons.

    A bullet is a line whose first non-space character is ``-``. Its scope
    is the union of an inline ``[tag1,tag2]`` prefix and the most recent
    ``## <tag>`` section heading. Control characters are stripped. The
    returned ``source`` is the file NAME only (the caller resolves the dir
    and names the file), so no lesson can carry a path that walks out of
    the configured lessons directory.
    """
    if isinstance(path, str):
        path = Path(path)
    p = path.resolve()
    source = p.name
    lessons: list[Lesson] = []
    current_heading: str | None = None
    try:
        lines = p.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError):
        return []
    for i, raw in enumerate(lines, start=1):
        heading = _section_tag(raw)
        if heading:
            current_heading = heading
            continue
        bullet = raw.lstrip()
        if not bullet.startswith("-"):
            continue
        body = bullet[1:].strip()
        inline_tags: set[str] = set()
        m = _TAG_PREFIX_RE.match(body)
        if m:
            tags = [t.strip() for t in m.group(1).split(",") if t.strip()]
            inline_tags = {t for t in tags}
            body = m.group(2).strip()
        if not body:
            continue
        scope = set(inline_tags)
        if current_heading:
            scope.add(current_heading)
        lessons.append(Lesson(
            text=_strip_control(body),
            scope=frozenset(scope),
            source=source,
            line=i,
        ))
    return lessons


def load_lessons_dir(lessons_dir: Path | str) -> list[Lesson]:
    """All lessons from the ``.md`` files directly inside ``lessons_dir``.

    Only the immediate children are read; subdirectories are ignored (a
    subclassed catalog convention is not assumed). Files are read in sorted
    name order so the result is deterministic across runs. A missing or
    unreadable directory yields an empty list — never raises.
    """
    if isinstance(lessons_dir, str):
        lessons_dir = Path(lessons_dir)
    try:
        if not lessons_dir.is_dir():
            return []
        files = sorted(lessons_dir.glob("*.md"))
    except OSError:
        return []
    try:
        root = lessons_dir.resolve()
    except OSError:
        return []
    out: list[Lesson] = []
    for f in files:
        try:
            # A symlink (or a name resolving elsewhere) must not pull in a
            # file outside the configured directory.
            if f.resolve().parent != root or not f.resolve().is_file():
                continue
            out.extend(parse_lessons_file(f))
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------------
# Selection (deterministic, no model call)
# ---------------------------------------------------------------------------

def select_lessons(
    task_text: str,
    items: list[Lesson],
    *,
    max_lessons: int = 8,
    max_chars: int = 1600,
    refused_cap: int = 16,
) -> dict:
    """Pick the lessons whose scope matches ``task_text``.

    The match is deterministic and read-only: a lesson is a candidate when
    any of its scope tags OR its source catalog name (extension-stripped)
    appears in the lowercased ``task_text``. Candidates are ranked by the
    number of matched tags (fewest-unrelated first), then by source file
    order within the catalog, then by line — so the result is repeated
    exactly every call.

    At most ``max_lessons`` lessons are returned, within ``max_chars`` of
    text (measured on the bodies alone, not the fence). A lesson whose body
    carries an injection shape (``instruction_shapes``) is refused and
    reported in ``refused`` — its scope match earns it nothing.

    Returns::

        {"lessons": [Lesson, ...],  # ≤ max_lessons, ≤ max_chars bodies
         "chars": int,              # total body characters selected
         "matched": int,            # ≤ max_lessons minus refused before cap
         "refused": [str, ...]}     # source names of refused lessons
    """
    needle = (task_text or "").lower()
    if not needle:
        return {"lessons": [], "chars": 0, "matched": 0, "refused": []}

    refused: list[str] = []
    candidates: list[tuple[int, int, int, Lesson]] = []
    for order, lesson in enumerate(items):
        tags = _matching_tags(lesson, needle)
        if not tags:
            continue
        if instruction_shapes(lesson.text):
            if lesson.source and lesson.source not in refused:
                refused.append(lesson.source)
            continue
        # Defense in depth: parse_lessons_file emits only basenames, but a
        # hand-constructed Lesson (or a future writer) must not carry a path
        # that points outside the configured lessons directory.
        if "/" in lesson.source or "\\" in lesson.source:
            continue
        # Rank: (matched-count DESC → -count ascending sort key, source
        # order, line). Source file name and line are capture-only ranking
        # ties; they never affect which lesson is refused.
        candidates.append((-len(tags), order, lesson.line, lesson))

    candidates.sort(key=lambda c: (c[0], c[1], c[2]))
    picked: list[Lesson] = []
    budget = 0
    for _neg, _order, _line, lesson in candidates:
        if len(picked) >= max_lessons:
            break
        cost = len(lesson.text) + 1  # +1 for the newline that separates them
        if budget + cost > max_chars and picked:
            break
        picked.append(lesson)
        budget += cost
    return {
        "lessons": picked,
        "chars": budget,
        "matched": len(picked),
        "refused": refused[:refused_cap],
    }


def _matching_tags(lesson: Lesson, needle: str) -> set[str]:
    """The lesson's tags present in the task text.

    A tag matches when it is a substring of the task text (the task says the
    tag verbatim) OR any ``-``/``_`` fragment of it appears in the task text
    (the task says a keyword the tag is built from — a tag ``orca-cards``
      matches a task about "orca"). So a maintainer's catalog about
      "the orca card mistake" surfaces on a task that only says "orca".
    """
    matched: set[str] = set()
    for tag in lesson.scope:
        t = tag.lower()
        if not t:
            continue
        if t in needle:
            matched.add(t)
            continue
        for frag in re.split(r"[-_\s]+", t):
            if len(frag) >= 3 and frag in needle:
                matched.add(t)
                break
    return matched
