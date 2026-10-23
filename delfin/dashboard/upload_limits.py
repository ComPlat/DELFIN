"""Size and count limits for what the dashboard accepts as an upload.

Input: an uploaded file's bytes, the chunks of a browser-side chunked
upload, or a ZIP archive that is unpacked into the calculations tree.
Output: the data is written, or :class:`UploadRefused` (a ``ValueError``)
is raised with a message that names the limit and the variable that sets
it; nothing partial is left behind.

Limits (bytes / entries; an unset, unparsable or non-positive value means
the default):

``DELFIN_UPLOAD_MAX_BYTES``        one uploaded file, default 2 GiB
``DELFIN_ZIP_MAX_UNPACKED_BYTES``  everything one ZIP unpacks to, default 8 GiB
``DELFIN_ZIP_MAX_ENTRIES``         files in one ZIP, default 20 000

The unpacked size is counted on the bytes actually written, not on the
sizes the archive declares (those are checked first as a fast refusal, but
a crafted archive can declare anything). Files are created with
``O_EXCL | O_NOFOLLOW``, so an existing name or a planted link at the
destination is never written through.

Every ``FileUpload`` of the dashboard writes through here: the Calculations
browser and the Literature tab create new files (:func:`write_new_file`),
the ORCA Builder's job directory and the agent's upload folder replace what
is there (:func:`write_replacing`); the name the browser sends is reduced
to a file name first (:func:`upload_name`). Uploaded files are 0600.

Rejected: a compression-ratio limit. An absolute bound on what is written
covers the same denial of service without refusing legitimately sparse
data (cube files and zero-padded binaries compress 100:1 and more).
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Iterable
from pathlib import Path
from typing import BinaryIO

UPLOAD_ENV = "DELFIN_UPLOAD_MAX_BYTES"
UNPACKED_ENV = "DELFIN_ZIP_MAX_UNPACKED_BYTES"
ENTRIES_ENV = "DELFIN_ZIP_MAX_ENTRIES"

DEFAULT_UPLOAD_BYTES = 2 * 1024 ** 3
DEFAULT_UNPACKED_BYTES = 8 * 1024 ** 3
DEFAULT_ENTRIES = 20_000

#: A chunked upload announces how many chunks follow; a placeholder list of
#: that length is allocated on the first chunk. Bounded so the announcement
#: alone cannot exhaust memory.
MAX_CHUNKS = 1_000_000

_COPY_BLOCK = 1024 * 1024


class UploadRefused(ValueError):
    """An upload or archive over one of the limits."""


def _limit(name: str, default: int) -> int:
    try:
        value = int(str(os.environ.get(name, "")).strip())
    except ValueError:
        return default
    return value if value > 0 else default


def max_upload_bytes() -> int:
    return _limit(UPLOAD_ENV, DEFAULT_UPLOAD_BYTES)


def max_unpacked_bytes() -> int:
    return _limit(UNPACKED_ENV, DEFAULT_UNPACKED_BYTES)


def max_entries() -> int:
    return _limit(ENTRIES_ENV, DEFAULT_ENTRIES)


def _size(n: int) -> str:
    for unit in ("B", "KiB", "MiB", "GiB"):
        if n < 1024 or unit == "GiB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n} B"


def check_upload_size(nbytes: int, name: str = "upload") -> None:
    limit = max_upload_bytes()
    if nbytes > limit:
        raise UploadRefused(
            f"{name} is {_size(nbytes)}, over the upload limit of "
            f"{_size(limit)} ({UPLOAD_ENV}).")


def check_chunk_plan(chunk_total: int, received_chars: int,
                     name: str = "upload") -> None:
    """Refuse a chunked upload that announces or delivers too much.

    *received_chars* counts the base64 text received so far for one file.
    """
    if chunk_total > MAX_CHUNKS:
        raise UploadRefused(
            f"{name} announces {chunk_total} chunks (at most {MAX_CHUNKS}).")
    encoded_limit = (max_upload_bytes() * 4) // 3 + 4
    if received_chars > encoded_limit:
        raise UploadRefused(
            f"{name} is over the upload limit of {_size(max_upload_bytes())} "
            f"({UPLOAD_ENV}).")


def check_archive(infos: Iterable, name: str = "archive") -> None:
    """Fast refusal from what a ZIP declares: entry count and total size."""
    files = [info for info in infos if not info.is_dir()]
    if len(files) > max_entries():
        raise UploadRefused(
            f"{name} holds {len(files)} files, over the limit of "
            f"{max_entries()} ({ENTRIES_ENV}).")
    declared = sum(max(0, int(getattr(info, "file_size", 0) or 0))
                   for info in files)
    if declared > max_unpacked_bytes():
        raise UploadRefused(
            f"{name} unpacks to {_size(declared)}, over the limit of "
            f"{_size(max_unpacked_bytes())} ({UNPACKED_ENV}).")


class Budget:
    """Bytes that may still be written while unpacking one archive."""

    def __init__(self, total: int | None = None, name: str = "archive"):
        self.left = max_unpacked_bytes() if total is None else int(total)
        self.name = name

    def spend(self, n: int) -> None:
        self.left -= n
        if self.left < 0:
            raise UploadRefused(
                f"{self.name} unpacks to more than "
                f"{_size(max_unpacked_bytes())} ({UNPACKED_ENV}).")


def write_new_file(destination: Path, source: BinaryIO,
                   budget: Budget | None = None) -> int:
    """Copy *source* into a file that must not exist yet; return its size.

    The file is created ``O_EXCL | O_NOFOLLOW``, owner-only (0600): on a
    shared host an upload is the uploader's until they share it.
    On any error -- the budget exceeded included -- the partial file is
    removed and the error re-raised.
    """
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(destination), flags, 0o600)
    written = 0
    try:
        with os.fdopen(fd, "wb") as sink:
            while True:
                block = source.read(_COPY_BLOCK)
                if not block:
                    break
                if budget is not None:
                    budget.spend(len(block))
                sink.write(block)
                written += len(block)
    except BaseException:
        try:
            os.unlink(str(destination))
        except OSError:
            pass
        raise
    return written


def write_new_bytes(destination: Path, payload: bytes) -> int:
    """:func:`write_new_file` for bytes already in memory."""
    import io

    return write_new_file(destination, io.BytesIO(bytes(payload)))


def upload_name(name: object) -> str:
    """The file name an upload is saved under: the last path component.

    The browser sends a bare name, but the value reaches the kernel from
    the page, and a page is not trusted to have sent one. ``../x``,
    ``/a/b/x`` and ``a\\x`` become ``x``; a name that is empty, ``.``,
    ``..`` or holds a NUL is refused.
    """
    raw = str(name or "").replace("\\", "/")
    base = raw.rsplit("/", 1)[-1].strip()
    if not base or base in (".", "..") or "\0" in base:
        raise UploadRefused(f"{raw!r} is not a file name.")
    return base


def write_replacing(destination: Path, payload: bytes) -> int:
    """Write *payload* to *destination*, replacing what is there; return its size.

    For the uploads whose destination may already exist (a job directory
    the ORCA Builder saves into again, the agent's upload folder). The
    bytes go to a new owner-only (0600) file beside the destination,
    created ``O_EXCL | O_NOFOLLOW``, which is then renamed over it: a link
    that sits at the destination is replaced, never written through, and
    a reader never sees half a file. A directory at the destination is an
    error (``IsADirectoryError``/``OSError``); so is a payload over
    ``DELFIN_UPLOAD_MAX_BYTES`` (:class:`UploadRefused`). Nothing is left
    behind on any error.
    """
    import secrets

    destination = Path(destination)
    payload = bytes(payload)
    check_upload_size(len(payload), destination.name)
    if destination.is_dir() and not destination.is_symlink():
        raise IsADirectoryError(str(destination))
    tmp = destination.with_name(
        f".{destination.name}.{secrets.token_hex(4)}.upload")
    write_new_bytes(tmp, payload)
    try:
        os.replace(str(tmp), str(destination))
    except BaseException:
        try:
            os.unlink(str(tmp))
        except OSError:
            pass
        raise
    return len(payload)


def remove_tree(path: Path) -> None:
    shutil.rmtree(str(path), ignore_errors=True)
