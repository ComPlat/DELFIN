"""Extract fenced code blocks from DELFIN's markdown files.

The recipe tests read the manuals themselves, so a new recipe is checked the
day it is written.  This module owns the shared notion of a "block":

* a block starts at a line of three or more backticks, optionally followed
  by a language tag (bash, ini, yaml, python, text, bibtex, ...),
* ends at the next fence line at or above that length,
* and may be marked not-runnable by an HTML comment placed on the line
  directly above the fence::

    <!-- recipe: not-runnable — output example, nothing to run -->

Everything else — every block without such a marker — must be runnable by
the standard of its own language tag, and a block with NO tag must still
be recognizable as one of the shapes the tests know (SMILES, XYZ body,
slash-command help), otherwise it is a failure rather than silence.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# All markdown files that carry recipes.  Kept in one place so the tests and
# any future doc share the same list.
MANUAL_FILES = [
    "docs/USER_MANUAL.md",
    "README.md",
    "docs/SETTINGS_AND_SETUP.md",
    "docs/RETRY_LOGIC.md",
    "docs/methodology.md",
]

# <!-- recipe: not-runnable — reason -->
NOT_RUNNABLE_RE = re.compile(
    r"^\s*<!--\s*recipe:\s*not-runnable\s*[—-]\s*(?P<reason>.+?)\s*-->\s*$"
)

_FENCE_RE = re.compile(r"^(?P<fence>```+)\s*(?P<lang>\S*)\s*$")


@dataclass
class Block:
    """One fenced code block from a markdown file."""

    path: str
    line: int          # 1-based line of the opening fence
    lang: str          # tag after the opening fence, "" when untagged
    text: str          # the block body, fences stripped, no trailing newline
    not_runnable_reason: str | None = field(default=None)

    @property
    def where(self) -> str:
        return f"{self.path}:{self.line}"

    def lines(self) -> list[str]:
        return self.text.splitlines()


def blocks_of(md_text: str, path: str) -> list[Block]:
    """Return every fenced block in *md_text* in document order."""
    lines = md_text.splitlines()
    blocks: list[Block] = []
    i = 0
    pending_marker: str | None = None
    while i < len(lines):
        m = _FENCE_RE.match(lines[i])
        if m and m.group("fence").startswith("```"):
            fence_len = len(m.group("fence"))
            lang = m.group("lang")
            body: list[str] = []
            j = i + 1
            while j < len(lines):
                close = _FENCE_RE.match(lines[j])
                if close and len(close.group("fence")) >= fence_len and not close.group("lang"):
                    break
                body.append(lines[j])
                j += 1
            blocks.append(
                Block(
                    path=path,
                    line=i + 1,
                    lang=lang,
                    text="\n".join(body),
                    not_runnable_reason=pending_marker,
                )
            )
            pending_marker = None
            i = j + 1
            continue
        marker = NOT_RUNNABLE_RE.match(lines[i])
        if marker:
            pending_marker = marker.group("reason").strip()
        elif lines[i].strip():
            # A marker only counts directly above the fence; intervening
            # prose resets it so a stray comment never silences a recipe.
            pending_marker = None
        i += 1
    return blocks


def load_blocks(repo_root: str | Path) -> list[Block]:
    """Read every manual file under *repo_root* and return its blocks."""
    root = Path(repo_root)
    out: list[Block] = []
    for rel in MANUAL_FILES:
        p = root / rel
        out.extend(blocks_of(p.read_text(encoding="utf-8"), rel))
    return out


def repo_root() -> Path:
    """Repository root as seen from this test file (…/tests/..)."""
    return Path(__file__).resolve().parent.parent
