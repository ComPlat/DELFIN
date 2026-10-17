"""The person's window onto the persistent task state: /state.

Three commands over the TaskState store (delfin/agent/task_state.py),
wired into the terminal by repl.py (KG1) as one line:

    from .state_commands import dispatch as _state_dispatch

* ``/state``            -> show(st): the whole store, bounded.
* ``/state edit <f>…``  -> edit(st, field, text): the PERSON names a
  field and a value; goal edits are marked as the person's instruction.
* ``/state clear-wait`` -> clear_wait(st): empty the waiting-on field.

And one guard the terminal must call for every incoming message:

* ``apply_incoming(origin, text)`` — only origin == "person" writes the
  goal / person's-instruction fields. A session or operator message is
  at most recorded and never becomes the goal: a message from another
  session (or any text inside a file or tool output) must not be able to
  steer the store, because steering is the person's alone. The person's
  newest instruction replaces any stored plan.
"""
from __future__ import annotations

# Fields /state edit accepts. Deliberately narrow: a typos is refused
# with the list, not guessed. "goal" is the alias for the task field.
_EDITABLE = {
    "goal": "task",
    "success criterion": "success_criterion",
    "success_criterion": "success_criterion",
}

# Clip for a person-set value; render clips again.
_MAX_EDIT_CHARS = 300

# Whose message may steer. Nothing else — not "session", not "operator",
# not a scheduled or tool-emitted message.
_PERSON_ORIGINS = ("person", "user")


def show(st) -> str:
    """The whole current state, bounded and scrubbed (render() does the
    ceiling and the secret scrubbing). An empty store says so."""
    body = st.render()
    if not body:
        return "Task state: empty (nothing committed yet)."
    return body


def edit(st, field: str, text: str) -> str:
    """Set one field from the terminal. The goal alias routes through
    set_person_instruction so the edit is marked as the person's."""
    key = str(field or "").strip().lower()
    target = _EDITABLE.get(key)
    if target is None:
        known = ", ".join(sorted(_EDITABLE))
        return (f"Unknown field {field!r}. Editable: {known} "
                f"(clear-wait clears waiting on).")
    value = str(text or "").strip()[:_MAX_EDIT_CHARS]
    if not value:
        return f"Nothing to set for {field!r}: empty value."
    if target == "task":
        st.set_person_instruction(value)
    else:
        st.set_success_criterion(value)
    return f"{field!r} set (person's edit, marked as theirs)."


def clear_wait(st) -> str:
    """Empty the waiting-on field."""
    st.waiting_for = ""
    return "Waiting-on cleared."


def apply_incoming(origin: str, text: str, st=None) -> None:
    """Feed one incoming message to the state.

    Only a person's message writes the goal. A session or operator
    message may not set goal or person's-instruction — even when its
    text says "new goal:" — so the state can never be steered by anyone
    but the person.

    The caller binds the store: the terminal passes the session's
    TaskState (``st``); without one, nothing is written (the guard must
    never invent a store path on its own — a wrong path would write
    outside the session's workspace).
    """
    if st is None:
        # No bound store: refuse to write rather than guess a path.
        return
    if not str(text or "").strip():
        return
    if str(origin or "").strip().lower() not in _PERSON_ORIGINS:
        # A note at most; the store keeps the standing goal untouched.
        # Deliberately NOT recorded in the goal/instruction fields.
        return
    st.set_person_instruction(str(text).strip()[:_MAX_EDIT_CHARS])
    st.save()
