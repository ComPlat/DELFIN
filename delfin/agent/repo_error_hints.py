"""Which section of the ``git-errors`` skill a piece of text calls for.

Input: the text a shell command printed, or the text of a user message.
Output: the matching sections of the pack skill ``git-errors``, rendered
for the tool result (the section itself) or for the turn's steering tail
(a one-line pointer to the skill).

Semantics. Each section of ``pack/skills/git-errors.md`` has one trigger
pattern here, keyed by the section's heading; a test pins the two lists
to each other. A tool result carries a matched section at most once per
session (``KitToolPermissions.repo_error_hints_shown``), and at most two
sections per result. A user message that names one of the error texts,
or asks about checkouts, worktrees, the dashboard's branch or its ports,
gets a pointer to call the skill; nothing is added otherwise.

Cost (len/4, measured). The skill's entry in the ``skill`` tool's listing
is the only always-sent text: 19 tokens per request in a code session. A
matched section adds 85-235 tokens (the held-branch case: 205 with the
header) to the one result that matched; the pointer adds about 60 tokens
to one turn; the whole skill, loaded through the tool, is 1,851.

Rejected: loading the whole skill into every prompt (about 1,850 tokens
per request, for a situation most sessions never meet); a separate pattern
table in the skill's front matter (the loader reads one ``key: value`` per
line, and a regex written there would not survive its quote stripping).
"""

from __future__ import annotations

import re
from typing import Iterable, Optional

SKILL_NAME = "git-errors"

#: Most sections a single tool result carries.
MAX_SECTIONS_PER_RESULT = 2

#: How much of a long output is scanned: git prints its error last.
_SCAN_TAIL_CHARS = 20_000

#: (section heading, pattern). The heading is the section's ``## `` line
#: in the skill, character for character.
TRIGGERS: tuple[tuple[str, str], ...] = (
    ("Branch held by another worktree",
     r"is already (?:used by worktree|checked out) at '"),
    ("Detached HEAD",
     r"You are in 'detached HEAD' state|^HEAD detached (?:at|from) \S"),
    ("index.lock exists",
     r"Unable to create '[^'\n]*index\.lock': File exists"
     r"|Another git process seems to be running"),
    ("Push rejected or branches diverged",
     r"! \[rejected\][^\n]*\((?:fetch first|non-fast-forward)\)"
     r"|Updates were rejected because"
     r"|Your branch and '[^'\n]+' have diverged"
     r"|Need to specify how to reconcile divergent branches"),
    ("Untracked files would be overwritten",
     r"untracked working tree files would be (?:overwritten|removed) by"),
    ("Local changes would be overwritten",
     r"Your local changes to the following files would be overwritten by"
     r"|Please commit your changes or stash them"),
    ("Merge conflict",
     r"CONFLICT \([^)\n]+\): Merge conflict in"
     r"|Automatic merge failed; fix conflicts"
     r"|^<{7} \S"),
    ("Not a git repository",
     r"fatal: not a git repository"),
    ("Shallow clone",
     r"shallow update not allowed|\(grafted\b"
     r"|fatal: --unshallow on a complete repository"
     r"|unknown revision or path not in the working tree"),
    ("LFS pointer files",
     r"^version https://git-lfs\.github\.com/spec/v1"
     r"|git-lfs[^\n]{0,40}command not found"
     r"|should have been pointers"),
    ("Permission denied (publickey)",
     r"Permission denied \(publickey"),
    ("DELFIN's own checkouts and servers",
     r"\.delfin/worktrees/[^\s'\"]+[^\n]{0,80}(?:fatal|error|not a git)"
     r"|already open in another session"
     r"|held by a live session"
     r"|voila[^\n]{0,200}Address already in use"
     r"|Address already in use[^\n]{0,200}voila"),
    ("Contributing a change",
     r"protected branch hook declined|\bGH006\b"
     r"|\[remote rejected\][^\n]*protected"),
)

_COMPILED = tuple((h, re.compile(p, re.IGNORECASE | re.MULTILINE))
                  for h, p in TRIGGERS)

#: What a user writes when asking about DELFIN's own machinery, in the two
#: languages users write in. Matched only against user messages: these
#: words in a tool's output are ordinary content, not a question.
_USER_TOPICS = re.compile(
    r"\bworktrees?\b|\bcheckouts?\b|\barbeitskopie\w*"
    r"|\bdashboard\W{1,3}(?:\w+\W{1,3}){0,3}?(?:branch|zweig|port|server|kernel)\w*"
    r"|\b(?:branch|zweig|port|server|kernel)\w*\W{1,3}(?:\w+\W{1,3}){0,3}?(?:des |of the |vom )?dashboard"
    r"|\bvoila\b|\bwriter[- ]lock|\bsession lock",
    re.IGNORECASE)

#: A user asking for a pull request or a merge of one.
_USER_CONTRIBUTING = re.compile(
    r"\bpull[- ]requests?\b|\b(?:open|create|merge|erstell\w*|mach\w*)\s+"
    r"(?:\w+\s+){0,2}?(?:a |the |einen |den )?PR\b|\bgh pr\b",
    re.IGNORECASE)


def match(text: str) -> list[str]:
    """Headings whose trigger occurs in ``text``, in skill order."""
    body = str(text or "")[-_SCAN_TAIL_CHARS:]
    if not body.strip():
        return []
    return [h for h, rx in _COMPILED if rx.search(body)]


def user_topics(text: str) -> list[str]:
    """Headings a user message calls for: its error texts, plus the
    DELFIN section when it asks about checkouts, worktrees, the dashboard's
    branch or ports."""
    found = match(text)
    for heading, rx in (("DELFIN's own checkouts and servers", _USER_TOPICS),
                        ("Contributing a change", _USER_CONTRIBUTING)):
        if heading not in found and rx.search(str(text or "")):
            found.append(heading)
    return found


def sections(workspace: Optional[str] = None) -> dict[str, str]:
    """The skill's sections, heading -> text (heading line included).

    Read through the skill loader, so a project or user override of
    ``git-errors`` is what a hint quotes. Empty when no such skill loads.
    """
    from . import skills as _skills
    sk = _skills.get_skill(SKILL_NAME, workspace)
    if sk is None:
        return {}
    out: dict[str, str] = {}
    current: Optional[str] = None
    lines: list[str] = []
    for line in sk.body.splitlines():
        if line.startswith("## "):
            if current is not None:
                out[current] = "\n".join(lines).strip()
            current = line[3:].strip()
            lines = [line]
        elif current is not None:
            lines.append(line)
    if current is not None:
        out[current] = "\n".join(lines).strip()
    return out


def tool_hint(text: str, shown: Optional[set] = None,
              workspace: Optional[str] = None) -> str:
    """The matched sections for one tool result, or "".

    ``shown`` holds the headings already delivered in this session; the
    ones delivered now are added to it.
    """
    found = [h for h in match(text) if shown is None or h not in shown]
    if not found:
        return ""
    texts = sections(workspace)
    picked = [h for h in found if h in texts][:MAX_SECTIONS_PER_RESULT]
    if not picked:
        return ""
    if shown is not None:
        shown.update(picked)
    return (f"[{SKILL_NAME}] This output is a known git situation. From "
            f"the {SKILL_NAME} skill (all of it: skill(name="
            f"'{SKILL_NAME}')):\n\n" + "\n\n".join(texts[h] for h in picked))


def user_hint(text: str) -> str:
    """A one-line pointer for a user message that calls for the skill."""
    found = user_topics(text)
    if not found:
        return ""
    return (f"[{SKILL_NAME}] The message concerns: {'; '.join(found)}. "
            f"Before answering, call skill(name='{SKILL_NAME}') and follow "
            "the safe fix it gives; do not offer forced or destructive "
            "options in its place.")


def headings(body_lines: Iterable[str]) -> list[str]:
    """The ``## `` headings of a skill body, in order (for the parity test)."""
    return [ln[3:].strip() for ln in body_lines if ln.startswith("## ")]


__all__ = ["SKILL_NAME", "TRIGGERS", "match", "user_topics", "sections",
           "tool_hint", "user_hint", "headings"]
