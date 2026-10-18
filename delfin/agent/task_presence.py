"""Is there a task from the person in this session? (package B1)

Measured in wave 17 (2026-10-17): 25 sessions lost their first message
(a stale pause stopped the turn that carried it) and then IMPROVISED:
they searched ``~/.delfin`` memory, read old ``.gate/HANDOFF.md`` files
and old wave plans, and two hit severe refusals. The right behaviour is
to notice there is no task, say so, and ask the person.

``assess`` answers that from one history plus the workspace. It is PURE:
no I/O beyond reading the named workspace files for the stale check, no
engine, no model call, never raises (every error reads as "present" —
an assessment that throws must never block or rewrite a real turn).

How the person is recognised. Everything the harness tells the model
between the person's messages — operator notes, wake-ups, compaction
summaries, other sessions' messages — travels as role ``user`` but is
MARKED: ``HARNESS_KEY`` (``harness_rows.HARNESS_KEY``) on a history row,
or it arrived as ``HarnessText`` (a wake prompt, a session message).
Only an UNMARKED user row is the person's own words. Consequences:

* a history with only harness notes -> ``missing``;
* a tool result or file content that says "your task is ..." is not a
  user row at all, so hostile file text can never produce ``present``;
* a session message (harness-marked, rendered by
  ``session_messages.render``) is never the person.

``stale`` is the wave-17 trap itself: no person message in the history,
but task-shaped text (a HANDOFF, a plan) sits in the workspace — the
only thing the session could improvise from. ``missing`` is no task
anywhere; ``present`` means the person spoke a task in this session.

``prompt`` is the turn's own prompt argument as the caller received it.
A plain, non-empty ``str`` is a person typing NOW (a real first-turn
task must not be answered with "you have no task"). ``HarnessText`` is
the harness's own voice (a wake, a continuation) and never counts.
"""

from __future__ import annotations

import re

from .harness_rows import HARNESS_KEY, HarnessRow

#: The one plain instruction, said the same way on every surface.
NO_TASK_INSTRUCTION = (
    "You have no task from the person in this session. Say so and ask; "
    "do not search memory, other folders or old notes for one."
)

_TASK_FILE_NAMES = re.compile(
    r"(HANDOFF|PLAN|handoff|plan)", )

#: Text that reads as a task left for the session. Deliberately loose:
#: stale is a nudge ("an old note may look like your task"), not a
#: refusal, so a false positive costs one sentence, not a refusal.
_TASK_SHAPE = re.compile(
    r"(task|phase|implement|fix|build|write|test|commit|review|deliver)",
    re.IGNORECASE)

#: A history row that is the person's own words: role ``user`` and no
#: harness mark of either carrier.
def _is_person_row(row: object) -> bool:
    if not isinstance(row, dict):
        return False
    if str(row.get("role") or "") != "user":
        return False
    if row.get(HARNESS_KEY) or isinstance(row, HarnessRow):
        return False
    return bool(str(row.get("content") or "").strip())


def _person_spoke(history) -> bool:
    """True when at least one unmarked user row is in the history."""
    for row in (history or []):
        if _is_person_row(row):
            return True
    return False


def _taskish(text: str) -> bool:
    return bool(_TASK_SHAPE.search(str(text or "")))


def _workspace_task_notes(workspace) -> list[str]:
    """Task-shaped text sitting in the workspace, best-effort, bounded.

    Reads at most the few files a session actually improvises from
    (wave 17: HANDOFF.md, WELLE plans). A workspace that is None, not a
    path-like object, or unreadable contributes nothing — and any error
    reads as "no notes", never as a raise.
    """
    if workspace is None:
        return []
    try:
        base = workspace / ".gate"
        names = [n for n in ("HANDOFF.md",)
                 if (base / n).exists()]
        root_names = [n for n in list(workspace.glob("*.md"))
                      if _TASK_FILE_NAMES.search(n.name)]
        texts: list[str] = []
        for path in names and [base / n for n in names]:
            try:
                texts.append(str(path.read_text(encoding="utf-8",
                                                errors="replace"))[:4000])
            except Exception:
                pass
        for path in root_names[:5]:
            try:
                texts.append(str(path.read_text(encoding="utf-8",
                                                errors="replace"))[:4000])
            except Exception:
                pass
        return [t for t in texts if _taskish(t)]
    except Exception:
        return []


def assess(history, *, workspace, prompt=None):
    """Is there a task from the person? -> (status, reason).

    ``history`` is the session's message list (``engine.messages``);
    ``workspace`` the workspace root (or None); ``prompt`` the turn's
    own prompt as the caller received it. Returns one of

    * ``("present", reason)``  — the person's words are in the history
      (or the prompt is a plain non-empty string: they are typing now);
    * ``("stale", reason)``   — no person message, but task-shaped text
      sits in the workspace (the improvisation trap);
    * ``("missing", reason)`` — no task anywhere: ask the person.

    Never raises; an unexpected history shape reads as the safe answer
    for it (garbage history has no person's words in it).
    """
    try:
        # A plain typed prompt IS the person speaking now — but only on
        # an otherwise-empty history (the real first turn). Once turns
        # have run without a person row, a plain prompt is just as
        # likely a harness continuation that lost its type (wave 17:
        # the lost first message), so it does not rescue the session.
        if prompt is not None and type(prompt) is str \
                and str(prompt).strip() and not list(history or []):
            return ("present",
                    "the turn's prompt is a plain message (the person "
                    "is typing this turn)")
        if _person_spoke(history):
            return ("present",
                    "an unmarked user row (the person's own words) is "
                    "in this session's history")
        notes = _workspace_task_notes(workspace)
        if notes:
            return ("stale",
                    "no person message in this session's history, but "
                    "task-shaped text sits in the workspace (an old "
                    "HANDOFF or plan may look like a task)")
        return ("missing",
                "no person message in this session's history and no "
                "task text in the workspace")
    except Exception:
        # An assessment that throws must never block a real turn: read
        # as present so the hook stays out of the way.
        return ("present", "assessment failed; staying out of the way")


def instruction_for(status: str) -> str:
    """The one plain instruction for a missing/stale turn, or "".

    Nothing is refused by this module: it only informs.
    """
    if status == "present":
        return ""
    return NO_TASK_INSTRUCTION
