"""Reading and writing a JSON state file so a damaged one is not replaced.

Input: a path, and for a write the object to store. Output: the parsed
object, or for a write nothing; a damaged file is moved aside and its
new location returned.

The problem this exists for was measured on 2026-11-17 in two stores
independently. Each had one read that answered two different questions
with the same value:

    except (FileNotFoundError, json.JSONDecodeError):
        return {}          # "absent" and "unreadable" are not the same

A reader may conflate them -- a damaged store must not stop the agent
when its value can come from somewhere else. A WRITER may not: it
builds the new file from that empty answer and the old contents are
gone. Measured in ``memory_store``: two remembered facts, a file
truncated the way a full disk or a kill mid-write truncates it, one
``save_memory`` call, and the two facts were gone with no copy
anywhere. ``credentials`` lost two API keys the same way.

So the read here distinguishes the two cases and the caller chooses,
and the write never lands on top of something it could not read.

Rejected alternatives:

* A writer that simply refuses on a damaged file. The person then has no
  way forward that does not involve a terminal and a text editor, which
  a standing rule of this project rules out.
* Keeping one ``.bak`` slot. The second damage overwrites the evidence of
  the first, which is how a neighbouring store in this codebase loses the
  original. The quarantine name carries the date and is reserved with
  O_EXCL.
* Repairing the file's contents (salvaging the parsable prefix). A
  truncated JSON prefix is not a document, and guessing where a record
  ended invents data. The bytes are kept untouched for a person or a
  later tool to look at.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
from pathlib import Path
from typing import Any


class StateFileUnreadable(OSError):
    """The file exists but cannot be parsed, so a writer must not replace it."""


def read_json(path: str | Path, *, default: Any = None,
              strict: bool = False) -> Any:
    """The parsed file, or *default* where it is absent or unreadable.

    ``strict`` raises :class:`StateFileUnreadable` instead of returning
    *default* for a file that exists but cannot be read or parsed. An
    absent file is NOT an error in either mode: nothing was ever stored,
    so there is nothing to lose by writing.
    """
    p = Path(path)
    if not p.exists():
        return default
    try:
        text = p.read_text(encoding="utf-8")
    except OSError as exc:
        if strict:
            raise StateFileUnreadable(
                f"{p} cannot be read ({exc.__class__.__name__})") from None
        return default
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        if strict:
            raise StateFileUnreadable(f"{p} is not readable JSON: {exc}") \
                from None
        return default


def quarantine(path: str | Path) -> Path | None:
    """Move a file aside, keeping its bytes and its mode; return where.

    None when there was nothing to move. A rename, so the contents are
    never rewritten and the mode travels with the file -- whatever is in
    it stays as readable, or as private, as it was.

    The destination is reserved with O_EXCL before the rename because
    ``os.replace`` overwrites its destination without a word: two
    corruptions inside one second would otherwise leave only the second
    rescue, which is the single-slot loss the dated name exists to avoid.
    """
    p = Path(path)
    if not p.exists():
        return None
    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    for n in range(100):
        suffix = f".unreadable-{stamp}" + (f"-{n}" if n else "")
        aside = p.with_name(p.name + suffix)
        try:
            fd = os.open(aside, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            continue
        except OSError:
            return None
        os.close(fd)
        try:
            os.replace(p, aside)
        except OSError:
            try:
                aside.unlink()
            except OSError:
                pass
            return None
        return aside
    return None


def read_json_for_write(path: str | Path, *,
                        default: Any = None) -> tuple[Any, Path | None]:
    """What to modify, and where a damaged file was moved to.

    One place decides what a writer does about a damaged file, so several
    writers of one store cannot drift apart on it.
    """
    try:
        return read_json(path, default=default, strict=True), None
    except StateFileUnreadable:
        return default, quarantine(path)


def write_json_atomic(path: str | Path, data: Any, *,
                      mode: int = 0o600, indent: int = 2) -> None:
    """Write *data* as JSON so the file is never seen half-written.

    A temporary file in the same directory, then a rename: a reader
    either sees the previous file or the new one. The plain
    ``path.write_text`` this replaces truncates first, so a crash or a
    full disk between the truncation and the last byte leaves exactly
    the damaged file the rest of this module is about.

    *mode* is applied to the temporary file before any content is
    written, so a store holding secrets is never briefly world-readable.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
    try:
        payload = json.dumps(data, indent=indent, ensure_ascii=False)
        os.write(fd, payload.encode("utf-8"))
        os.fsync(fd)
    finally:
        os.close(fd)
    os.replace(tmp, p)
