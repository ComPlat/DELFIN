"""What a turn hands the model besides the user's text, for every surface.

Input: the engine and the text of one turn. Output: the ``memory_context``
the engine puts into the system prompt (the project instruction files and
the greeting hint) and the user message (the text plus the verify hint).

One module because the three entry paths -- the terminal session
(``repl.run_turn``), the one-shot run (``cli.cmd_run``) and the dashboard
(``tab_agent``) -- must not each compose this themselves. They did: the
dashboard loaded ``DELFIN.MD`` / ``AGENTS.md`` and added both hints, the two
terminal paths sent none of it, so the same workspace and the same prompt
reached the model as two different requests, and the instruction file
``delfin-agent init`` scaffolds was never read by the program that wrote it.

A turn the harness starts on its own (``HarnessText``: a wake, a
continuation) gets the project rules but neither hint: both hints judge
what the user typed.

Rejected: loading the instruction files inside the engine. A sub-agent's
engine already receives them through ``subagents._inherited_context``
(capped for a child), so an engine-level load would send them twice to
every child.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

# Fact-question shapes: the user asks for exact ORCA / chemistry keyword
# names, which is where confident invention happens. German, because it
# recognises what the user writes.
_VERIFY_HINT_PATTERNS = re.compile(
    r"(?i)("
    r"welche[s]?\s+(?:exakten?\s+)?keyword|"        # "welche keyword(s)"
    r"keyword[- ]?name|"                              # "keyword-name"
    r"\bexakt(?:e|en)?\b.*\b(?:keyword|name|block)|"  # "exakte ... keyword"
    r"nutzt\s+orca|"                                  # "nutzt ORCA"
    r"in\s+orca|"                                     # "in ORCA"
    r"gibt\s+es\s+(?:ein|den|das)?\s*\b\w*\s*(?:keyword|block|flag|option)|"
    r"wie\s+hei[ßs]t\b|"                              # "wie heißt"
    r"%\w+\s*[- ]?block"                              # "%casscf-block"
    r")"
)

_VERIFY_HINT = (
    "(Factual question detected — look it up via search_docs/Read in "
    "the manual/code BEFORE answering; only mention keywords that "
    "are actually documented, none from memory.)"
)

_GREETING_HINT = (
    "The user's message is a greeting and nothing else: reply with ONE "
    "short sentence; do NOT call any tools."
)


def verify_hint(user_text: str, mode: str = "") -> str:
    """The ground-yourself line for a factual ORCA/keyword question, or "".

    ``mode`` is accepted for the dashboard's call shape and not read: the
    hint applies in every mode. Texts over 4000 characters are documents,
    not questions, and get nothing.
    """
    if not user_text or len(user_text) > 4000:
        return ""
    if _VERIFY_HINT_PATTERNS.search(user_text):
        return _VERIFY_HINT
    return ""


def greeting_hint(user_text: str) -> str:
    """The one-sentence-reply line for a bare greeting, or "".

    Judged by ``AgentEngine.is_bare_greeting`` -- the rule the engine uses
    to withhold the tool schemas -- so "hallo, prüfe die Tabelle" is work,
    not a greeting.
    """
    try:
        from .engine import AgentEngine
        if AgentEngine.is_bare_greeting(user_text or ""):
            return _GREETING_HINT
    except Exception:
        pass
    return ""


def project_instructions(engine: Any) -> str:
    """``DELFIN.MD`` / ``AGENTS.md`` from the workspace up, or "".

    Anchored at the permissions' workspace and its extra roots, as the
    dashboard always did; without a permissions object, at the process's
    working directory (``load_project_memory``'s default).
    """
    try:
        from .project_memory import load_project_memory
        perms = getattr(engine, "kit_permissions", None)
        workspace = getattr(perms, "workspace", None) if perms else None
        extra = list(getattr(perms, "extra_workspace_dirs", ()) or ()) \
            if perms else []
        return load_project_memory(
            cwd=Path(workspace) if workspace else None,
            extra_roots=extra)
    except Exception:
        return ""


def _is_harness(text: Any) -> bool:
    try:
        from .harness_rows import HarnessText
        return isinstance(text, HarnessText)
    except Exception:
        return False


def memory_context(engine: Any, user_text: str) -> str:
    """The project rules, then the greeting hint; "" when neither applies."""
    parts = [project_instructions(engine)]
    if not _is_harness(user_text):
        parts.append(greeting_hint(user_text))
    return "\n\n".join(p for p in parts if p)


def user_message(user_text: str) -> str:
    """``user_text`` with the verify hint appended when it applies.

    A ``HarnessText`` comes back unchanged (same object, same type), so the
    engine still marks the row as harness-written.
    """
    if _is_harness(user_text):
        return user_text
    hint = verify_hint(user_text)
    return f"{user_text}\n\n{hint}" if hint else user_text


__all__ = ["greeting_hint", "memory_context", "project_instructions",
           "user_message", "verify_hint"]
