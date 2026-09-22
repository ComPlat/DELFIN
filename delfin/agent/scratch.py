"""Scratch directories that do not outlive the process that made them.

Measured on the login node on 2026-09-19: ``/tmp`` held 1071 directories
named ``delfin-gfnff-topo-*`` -- the GFN-FF bonding perception the
structure editor keeps for one molecule. It is dropped when the molecule
changes and by nothing at all when the session ends, so every closed tab
left one behind for good.

``atexit`` alone would not have caught them: a Voila kernel is killed,
and a killed process runs no handler. So a folder carries a stamp saying
which process owns it, and the next run sweeps the folders whose owner is
gone. The two halves cover each other -- the handler keeps the sweep from
having anything to do, the sweep covers the deaths the handler misses.

What the sweep must never do is take a directory it does not own. An
unstamped one was made by somebody else; one whose owner is running is
somebody's open session; one stamped on another host says nothing here,
because a pid means nothing across machines. All three are left alone,
and a folder is removed only on the positive answer that its owner is
gone.
"""

from __future__ import annotations

import atexit
import json
import os
import shutil
import socket
import tempfile
from pathlib import Path
from typing import List, Optional, Union

from delfin.agent.proc_identity import alive, process_start

#: The stamp inside a scratch directory. A dotfile, so a tool that lists
#: the calculation's own output does not show it.
OWNER_NAME = ".delfin-scratch-owner.json"

#: Paths this process made, for the clean way out.
_MINE: set = set()


def owned_dir(prefix: str, base: Union[str, Path, None] = None) -> Path:
    """A fresh scratch directory, stamped with this process's identity.

    A drop-in for ``tempfile.mkdtemp(prefix=...)`` that can be cleaned up
    afterwards by somebody who was not there when it was made.
    """
    folder = Path(tempfile.mkdtemp(prefix=prefix,
                                   dir=str(base) if base else None))
    pid = os.getpid()
    record = {"pid": pid, "proc_start": process_start(pid),
              "host": socket.gethostname()}
    try:
        (folder / OWNER_NAME).write_text(
            json.dumps(record), encoding="utf-8")
    except OSError:
        # An unstamped folder is only as bad as what we had before: the
        # sweep will leave it alone, and the exit handler still takes it.
        pass
    _MINE.add(str(folder))
    return folder


def release(folder: Union[str, Path, None]) -> None:
    """Remove a scratch directory now. Quiet when it is already gone."""
    if not folder:
        return
    shutil.rmtree(str(folder), ignore_errors=True)
    _MINE.discard(str(folder))


@atexit.register
def _release_mine() -> None:
    """The clean way out, for the exits that get one."""
    for path in list(_MINE):
        shutil.rmtree(path, ignore_errors=True)
    _MINE.clear()


def sweep(prefix: str, base: Union[str, Path, None] = None) -> List[Path]:
    """Remove stamped *prefix* directories whose owner is gone.

    Returns what was taken, so a caller can say so. Never raises: a
    tidy-up that brings a session down is worse than the mess.
    """
    room = Path(base) if base else Path(tempfile.gettempdir())
    taken: List[Path] = []
    try:
        entries = list(os.scandir(room))
    except OSError:
        return taken
    for entry in entries:
        try:
            if not entry.name.startswith(prefix) or not entry.is_dir():
                continue
        except OSError:
            continue
        folder = Path(entry.path)
        record = _stamp_of(folder)
        if record is None:
            continue          # not ours to judge
        if alive(record.get("pid"), str(record.get("proc_start") or ""),
                 str(record.get("host") or "")) is not False:
            continue          # running, or unanswerable -- leave it
        try:
            shutil.rmtree(folder, ignore_errors=True)
            taken.append(folder)
        except OSError:
            continue
    return taken


def _stamp_of(folder: Path) -> Optional[dict]:
    """The owner record, or None when there is none to read.

    None means "say nothing about this directory". An unreadable stamp is
    not evidence that the work is finished.
    """
    try:
        text = (folder / OWNER_NAME).read_text(encoding="utf-8")
    except OSError:
        return None
    try:
        record = json.loads(text)
    except (ValueError, TypeError):
        return None
    return record if isinstance(record, dict) else None
