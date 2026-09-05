"""Read and write the user's text files without changing what they are.

Every write tool used ``Path.read_text`` and ``Path.write_text``. Both
are convenience wrappers with two behaviours that are wrong for a file
somebody else owns:

* ``read_text`` opens in universal-newline mode, so a CRLF file arrives
  as LF and is written back as LF. Editing one word in a Windows-authored
  ``.csv`` rewrote **every line ending in the file**, and the diff shown
  to the user was the one line that was meant to change.
* ``errors="replace"`` turns any byte that is not UTF-8 into U+FFFD, and
  the next write makes that permanent. ``café`` in a LaTeX file became
  ``caf<?>`` — on a line the diff displayed as unchanged context.

Both reproduced exactly as described.

And ``write_text`` truncates before it writes. A 2000-byte file whose
write raised partway through was left at **zero bytes**, with the
pre-image still only in memory and discarded by the error return — so
the undo journal had nothing either. That is the one case where the
journal is the last copy, and it is the case it does not cover.

The internal stores in this project already know all of this: the change
journal writes through a temporary file and ``os.replace``, with
``newline=""``. This module gives the user's own files the same
treatment, and adds the part the journal does not need — remembering how
the file was encoded and terminated so it can be put back that way.

The text handed out is always normalised to ``\\n`` with no BOM, because
that is what the model matches its edit strings against. The convention
travels alongside it and is re-applied on write.
"""

from __future__ import annotations

import dataclasses
import os
import tempfile
from pathlib import Path

# Tried in order. utf-8 first and strictly: a file that decodes as utf-8
# is utf-8, and guessing further would be worse than the guess it fixes.
# cp1252 before latin-1 because it is the Windows default the users of
# these files actually have, and it decodes the same byte range plus the
# printable 0x80-0x9F block that latin-1 turns into control characters.
_ENCODINGS: tuple[str, ...] = ("utf-8", "cp1252", "latin-1")

_BOM = "﻿"


@dataclasses.dataclass(frozen=True)
class TextFile:
    """A file's content plus the three things a naive write destroys."""

    text: str
    """Normalised: LF line endings, no BOM. What edit strings match."""

    encoding: str
    newline: str
    """The file's own line ending: ``\\n``, ``\\r\\n`` or ``\\r``."""

    bom: bool

    def restore(self, text: str) -> str:
        """Put ``text`` back into this file's own convention."""
        out = text
        if self.newline != "\n":
            out = out.replace("\n", self.newline)
        if self.bom:
            out = _BOM + out
        return out


def _detect_newline(raw: str) -> str:
    """The dominant terminator, deciding ties towards the first one seen.

    A mixed file is a real thing (a patch applied to a CRLF file, a
    generator that forgot). Rewriting it wholesale to the majority is
    still a change nobody asked for, but it is the smaller one, and the
    alternative — tracking per line — is not worth its complexity here.
    """
    crlf = raw.count("\r\n")
    cr = raw.count("\r") - crlf
    lf = raw.count("\n") - crlf
    if crlf and crlf >= lf and crlf >= cr:
        return "\r\n"
    if cr and cr > lf:
        return "\r"
    return "\n"


def read_text_file(path: Path | str) -> TextFile:
    """Decode ``path`` without losing a byte and without translating.

    Raises ``UnicodeDecodeError`` only if no candidate encoding fits,
    which for the last candidate cannot happen — latin-1 maps every
    byte. A caller that wants to refuse binary content should check the
    content, not rely on a decode failure.
    """
    raw_bytes = Path(path).read_bytes()
    text = ""
    encoding = _ENCODINGS[-1]
    for candidate in _ENCODINGS:
        try:
            text = raw_bytes.decode(candidate)
            encoding = candidate
            break
        except UnicodeDecodeError:
            continue
    bom = text.startswith(_BOM)
    if bom:
        text = text[len(_BOM):]
    newline = _detect_newline(text)
    if newline != "\n":
        text = text.replace(newline, "\n")
    else:
        # A file that is mostly LF can still hold a stray CR.
        text = text.replace("\r\n", "\n")
    return TextFile(text=text, encoding=encoding, newline=newline, bom=bom)


def write_text_file(
    path: Path | str, text: str, *, like: TextFile | None = None,
) -> list[str]:
    """Write ``text`` atomically, in ``like``'s convention if given.

    Returns notes for anything the caller should pass on — currently one
    case: the new text could not be encoded the way the file was, so it
    was written as UTF-8. Changing a file's encoding is a real change
    and saying nothing about it is how the old code lost data.

    Atomic because a failed write must leave the previous content
    intact. The temporary file is created in the target's own directory
    so ``os.replace`` stays within one filesystem, and the original
    file's permission bits are carried over — a mode-0600 file must not
    come back readable by the whole group.
    """
    target = Path(path)
    notes: list[str] = []
    payload = like.restore(text) if like is not None else text
    encoding = like.encoding if like is not None else "utf-8"
    try:
        data = payload.encode(encoding)
    except UnicodeEncodeError:
        if encoding != "utf-8":
            data = payload.encode("utf-8")
            notes.append(
                f"{target.name} was {encoding} and the new text cannot be "
                f"written in it; the file is now UTF-8.")
        else:
            raise

    target.parent.mkdir(parents=True, exist_ok=True)
    mode: int | None = None
    try:
        mode = target.stat().st_mode & 0o7777
    except OSError:
        mode = None

    fd, tmp = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp",
                               dir=str(target.parent))
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        if mode is not None:
            os.chmod(tmp, mode)
        os.replace(tmp, target)
    except BaseException:
        # Including KeyboardInterrupt: a half-written temporary file left
        # beside the user's file is litter, and the point of the whole
        # exercise is that the original is still there.
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return notes


__all__ = ["TextFile", "read_text_file", "write_text_file"]
