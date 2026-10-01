"""The append transformation for write_file(mode="append").

Three of six sessions in Welle 4 wrote large files in pieces with

    cat >> file <<EOF   ...   EOF

in bash -- every one of them a write the gate has to ask about, because
write_file could not append. ``write_file(path, content, mode="append")``
is the replacement, and the point of this module is that the handler
does NOT need a second write path: appending is transformed into the
content a plain write would carry (old text + separator + tail), so the
read-baseline check, the mtime check, the atomic swap, the undo journal
(pre-image + post hash) and the diff all run unchanged.

``old_text is None`` means the file does not exist yet; an append to a
not-yet-existing file is a create (the piecewise pattern starts with
file creation, and refusing it would force a mode switch mid-stream).
"""

from __future__ import annotations


class AppendRejected(ValueError):
    """The append was refused before anything was written.

    The message says which rule fired, so the caller can surface it
    instead of silently doing something else (like a no-op write).
    """


def build_new_text(old_text: str | None, tail: str) -> str:
    """The full text a plain write_file call must carry for an append.

    ``old_text``: the file's current text, or ``None`` when it does not
    exist. ``tail``: the text to append.

    Rules, each with its reason:
    - a tail that is empty or only whitespace is refused: it would be a
      journaled change that changes nothing, and every such record
      spends one of the session's capped pre-image slots for noise.
    - a missing final newline in ``old_text`` is supplied: without it,
      appending would glue the tail onto the last line, and the next
      read would show a line that was never written.
    - ``old_text is None`` (file does not exist) is allowed and yields
      the tail alone, so the piecewise pattern can start with append
      calls throughout.

    Raises :class:`AppendRejected` on a refused tail.
    """
    if not isinstance(tail, str) or not tail.strip():
        raise AppendRejected(
            "append content is empty or whitespace only — nothing to "
            "append; a no-op write would be journalled as a change.")
    if old_text is None:
        return tail
    if not isinstance(old_text, str):
        raise AppendRejected("old_text must be a string or None")
    if old_text == "":
        return tail
    if not old_text.endswith("\n"):
        old_text = old_text + "\n"
    return old_text + tail
