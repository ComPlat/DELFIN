"""Safety review of skill proposals (Paket 3 of the learning wave).

A self-written skill is only ever a PROPOSAL. Before a human accepts
it, this module checks its text: an unsafe success must not become a
rule ("Skill Misevolution", arXiv 2608.12851 measured that agents left
unattended write unsafe skills in 21/21 configurations).

``check(text) -> list[str]`` returns one named finding per problem, an
empty list when the text is clean. Findings are REPORTED, never
silently repaired — a repaired unsafe proposal stays suspicious, and
the human sees what was found.

Design rule (wave contract): this module CALLS DELFIN's own live
classifiers — the functions in :mod:`delfin.agent.api_client` and
:mod:`delfin.agent.output_guard` that decide whether a bash command is
asked about, auto-allowed or refused at run time. Nothing here
re-implements, duplicates or widens them. Anything the live gate would
ask about or refuse is therefore a finding here.

What counts as a finding (phase 1, commands):

* a command on the deny list (``matches_bash_deny``) — refused live;
* a command the gate would not auto-allow (``matches_bash_auto_allow``
  is False) — asked about live, which a skill text cannot answer;
* a hidden payload (base64 into an interpreter) — not readable by any
  scanner, refused live.
"""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import Any

__all__ = ["check", "extract_commands"]


#: The review classifier is built once and reused; a fresh empty temp
#: dir would leak per call. This workspace never exists on disk, which
#: is exactly how a session's own fresh workspace behaves.
_SKILL_REVIEW_WORKSPACE = Path(
    os.environ.get("DELFIN_SKILL_REVIEW_WORKSPACE",
                   os.path.join(tempfile.gettempdir(),
                                "delfin-skill-review-workspace")))
_PERMS_CACHE: Any = None


#: A fenced code block; the tag names the language or is empty.
_FENCE_RE = re.compile(r"```([A-Za-z0-9_+-]*)[ \t]*\r?\n(.*?)```", re.DOTALL)

#: Inline code span: `...` (not ```` ``` ```` — fences are consumed first).
_INLINE_RE = re.compile(r"`([^`\n]+)`")

#: A shell-prompt line: ``$ cmd`` or ``% cmd`` at line start.
_PROMPT_LINE_RE = re.compile(r"^\s*[$%]\s*(.+?)\s*$", re.MULTILINE)

#: Line continuations and trailing comments: join and drop them so the
#: classifier sees the whole command, exactly as bash reads it.
_LINE_CONTINUATION_RE = re.compile(r"\\\s*\n\s*")
_TRAILING_COMMENT_RE = re.compile(r"\s+#\s.*$")

#: Languages whose fenced blocks hold shell commands.
_SHELL_FENCE_TAGS = {
    "", "bash", "sh", "shell", "zsh", "console", "shell-session",
}

#: A base64 payload piped into an interpreter — unreadable by any
#: scanner, and refused by the live gate for exactly that reason.
_BASE64_INTO_INTERPRETER_RE = re.compile(
    r"base64\b[^|;&]*\|\s*(?:sudo\s+)?(?:/\S*/)?"
    r"(?:bash|sh|zsh|python[^\s|;&]*)\b")

#: Program words that make a fenced block a command block even when the
#: tag says otherwise (```python with ``python -c`` inside is both).
_INTERPRETER_INVOCATION_RE = re.compile(
    r"^\s*(?:sudo\s+)?(?:/\S*/)?"
    r"(?:python[23]?(?:\.\d+)?|bash|sh|zsh|perl|ruby|node)\s+-\w*c",
    re.IGNORECASE)


def _normalise_block(body: str) -> str:
    """Join wrapped lines into one command string, drop comments."""
    joined = _LINE_CONTINUATION_RE.sub(" ", body)
    lines = []
    for line in joined.splitlines():
        line = _TRAILING_COMMENT_RE.sub("", line.rstrip()).strip()
        if line:
            lines.append(line)
    return lines


def extract_commands(text: str) -> list[tuple[str, str]]:
    """Every command candidate in a skill text, as ``(where, command)``.

    ``where`` names the place it was found ("code block", "prompt line",
    "inline code") so a finding can say what was looked at.
    """
    out: list[tuple[str, str]] = []
    consumed: list[tuple[int, int]] = []  # inline spans inside fences

    for m in _FENCE_RE.finditer(text):
        tag, body = m.group(1).lower(), m.group(2)
        start, end = m.span()
        consumed.append((start, end))
        lines = _normalise_block(body)
        if tag in _SHELL_FENCE_TAGS:
            out.extend(("code block", ln) for ln in lines)
        else:
            # A tagged non-shell block can still hold a shell command
            # (```text with ``bash -c '...'``); the interpreter-invocation
            # form is the one that runs code from its text.
            for ln in lines:
                if _INTERPRETER_INVOCATION_RE.match(ln):
                    out.append(("code block", ln))

    prose = text
    for start, end in consumed:
        prose = prose[:start] + "\n" * (prose[start:end].count("\n")) \
            + prose[end:]
    for m in _PROMPT_LINE_RE.finditer(prose):
        cmd = _TRAILING_COMMENT_RE.sub("", m.group(1).strip()).strip()
        if cmd:
            out.append(("prompt line", cmd))
    for m in _INLINE_RE.finditer(prose):
        cand = m.group(1).strip()
        # Only spans that parse as a command: a program word first.
        if cand and re.match(r"^[A-Za-z_][\w.-]*\s", cand + " "):
            out.append(("inline code", cand))
    return out


def _classifier():
    """A permissions object wired exactly like a live session's.

    ``KitToolPermissions`` with a workspace that does not exist is the
    same classification a session gets for its own (fresh) workspace:
    the deny list, the auto-allow list and the write-outside check all
    run unchanged. Built, cached, thrown away — nothing persists.
    """
    global _PERMS_CACHE
    if _PERMS_CACHE is None:
        from .api_client import KitToolPermissions
        _PERMS_CACHE = KitToolPermissions(
            workspace=_SKILL_REVIEW_WORKSPACE, mode="default")
    return _PERMS_CACHE


def check(text: str) -> list[str]:
    """Named findings in a draft skill text; empty list means clean.

    Phase 1 — commands: every command the text carries is judged by
    DELFIN's own live classifiers (imported, never copied):

    * ``matches_bash_deny``: refused at run time ⇒ finding;
    * ``matches_bash_auto_allow`` returning False: asked about at run
      time ⇒ finding — a skill text cannot answer a confirm dialog, so
      a command that needs one is not a rule the agent may apply
      unattended;
    * a base64 payload piped into an interpreter: unreadable by any
      scanner ⇒ finding.
    """
    findings: list[str] = []
    seen: set[str] = set()

    def _report(finding: str) -> None:
        if finding not in seen:
            seen.add(finding)
            findings.append(finding)

    for where, cmd in extract_commands(text):
        if _BASE64_INTO_INTERPRETER_RE.search(cmd):
            _report(
                f"{where}: hidden payload — a base64 blob piped into an "
                f"interpreter cannot be inspected: {cmd[:90]!r}")
            continue
        perms = _classifier()
        denied = perms.matches_bash_deny(cmd)
        if denied:
            _report(
                f"{where}: command on the deny list, refused at run time: "
                f"{cmd[:90]!r}")
            continue
        if not perms.matches_bash_auto_allow(cmd):
            _report(
                f"{where}: command the gate would ask about (not "
                f"auto-allowed), so a skill cannot apply it unattended: "
                f"{cmd[:90]!r}")
    return findings
