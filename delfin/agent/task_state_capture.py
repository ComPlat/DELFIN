"""The ONE writer of the per-session TaskState at session start.

``task_state.py`` renders ``Task:`` and ``Success criterion:`` and
``working_state._task_state_lines`` reads that store back into every
compaction block -- but nothing in the product wrote it: the
``state_commands`` writer API had zero call sites (its docstring claims
a repl wiring that does not exist; ``repl._OpenTasksBlock`` says the
same outright), so every reduction path -- compaction, fresh restart,
in-turn compaction -- lost what the person had asked and, when the
first message stated one, what "done" means.

This module is the writer. One writer, one place, called from the
shared turn entrances:

* ``repl.prepare_fresh_turn`` -- the headless CLI (``cli.cmd_run``) and
  the dashboard's turn path both route through it,
* ``repl.run_turn`` -- the terminal's attended path (and the worker's),
* so ``delfin-agent chat``/``run`` and the dashboard both capture.

Semantics:

* Only the FIRST real user message of the session writes (``first_task``
  is empty and the message is not harness transport). A later message is
  steering, not the session's task; it never overwrites the store here.
* The store is written ONLY when it does not already hold a task: the
  person's /state edit (state_commands) outranks the capture, and a
  restored session's existing store is kept -- the capture never
  demotes anything already stored.
* The success criterion is stored only when the message states one
  (``success criterion:`` / ``success criteria:``); nothing is invented.
* Best-effort and never raising: a capture that fails must not delay a
  turn.
"""
from __future__ import annotations

import re
from pathlib import Path

# Prefixes a task text must not carry. A wake or a system note is
# harness transport, not a person naming a task.
_HARNESS_PREFIXES = (
    "[System]", "[Message from", "[scheduled", "[watch]",
    "[Command results]", "[Working state", "[Conversation summary",
    "[Verify]", "task state (",
)

# The success-criterion phrase the capture recognises. Deliberately
# narrow: a criterion that is not stated is left empty, never invented.
_CRITERION_RE = re.compile(
    r"success criteri(?:a|on|ons)\s*[:\-–]\s*(.+)", re.IGNORECASE)

# Clip for the stored task; the store clips again at render time.
_MAX_TASK_CHARS = 300


def parse_success_criterion(prompt: str) -> str:
    """The success criterion the first message states, or "".

    Input: one user message. Output: the text after the first
    ``success criterion:`` marker, clipped; "" when the message states
    none. Never raises.
    """
    try:
        m = _CRITERION_RE.search(str(prompt or ""))
    except Exception:
        return ""
    if not m:
        return ""
    text = str(m.group(1) or "").strip()
    # A criterion runs to the end of its line: a multi-section message
    # must not swallow the sections after it.
    text = text.splitlines()[0].strip() if text else ""
    return text[:300].strip()


def _is_harness_transport(prompt: str) -> bool:
    text = str(prompt or "").lstrip()
    lowered = text.lower()
    return any(lowered.startswith(p.lower()) for p in _HARNESS_PREFIXES)


def capture_from_first_message(*, engine, workspace, prompt: str) -> bool:
    """Write the session's task (and stated success criterion) once.

    Input: the session's engine (``session_id`` is read), the workspace
    the store lives in, and the prompt about to be sent. Output: True
    when the store was written, False when it already held a task, the
    prompt is not a person's first message, or the write failed --
    never raises.

    Semantics: the store is the durable carrier every reduction path
    reads; this is its only product writer, so a session that reaches
    its first reduction with a stated task keeps it.
    """
    try:
        text = str(prompt or "")
        if not text.strip():
            return False
        if _is_harness_transport(text):
            return False
        from . import task_state as _ts
        workspace = Path(workspace) if workspace is not None else None
        if workspace is None:
            return False
        path = workspace / _ts.DEFAULT_FILENAME
        st = _ts.open(path)
        try:
            st.load()
        except Exception:
            st.begin()
        # The capture is the FIRST-message writer: a store that already
        # holds a task (a person's /state edit, or a restored session)
        # is never demoted.
        if str(getattr(st, "task", "") or "").strip():
            return False
        task = text.strip()[:_MAX_TASK_CHARS]
        st.set_person_instruction(task)
        criterion = parse_success_criterion(text)
        if criterion:
            st.set_success_criterion(criterion)
        st.save()
        return True
    except Exception:
        return False
