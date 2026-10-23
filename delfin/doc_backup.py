"""Versioned copies of a document, kept where the user can find them.

The dashboard's file browser makes one of these before it saves an edit,
and the agent has to make the same one before it saves its own — into
the same folder, under the same names. Two mechanisms writing two
histories of the same file is worse than one: the user opens the folder
they know and sees half the story.

The rules are the browser's, not new ones:

* copies go into a ``Backups`` folder beside the document, so the
  document's own folder stays a list of documents;
* each save keeps its own copy (``name.bak.ext``, ``name.bak2.ext``, …)
  rather than replacing the previous one — otherwise the second save
  finds a backup already there and keeps nothing, and a safety net for
  one edit is not a history;
* nothing is ever overwritten, and a failure to copy is reported rather
  than silently skipped: the caller decides whether to proceed without
  one.

This module holds no UI dependency on purpose, so both sides can use it.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

BACKUP_DIR_NAME = "Backups"
MAX_BACKUP_VERSIONS = 999


def backup_dir_for(path) -> Path:
    """The folder the copies of *path* belong in."""
    return Path(path).parent / BACKUP_DIR_NAME


def versioned_backup_path(path, folder=None) -> Optional[Path]:
    """The next free ``name.bakN.ext``, or None when the numbering is used up.

    The first copy is ``name.bak.ext`` without a number, matching what the
    browser already writes, so both fill one sequence instead of starting
    two.
    """
    path = Path(path)
    folder = Path(folder) if folder is not None else backup_dir_for(path)
    for version in range(1, MAX_BACKUP_VERSIONS + 1):
        mark = "bak" if version == 1 else f"bak{version}"
        candidate = folder / f"{path.stem}.{mark}{path.suffix}"
        if not candidate.exists():
            return candidate
    return None


def make_backup(path, folder=None) -> Optional[Path]:
    """Copy *path* aside before it is changed. Returns the copy, or None.

    None means no copy was made — the file does not exist yet, the folder
    could not be created, or the numbering is exhausted. Callers report
    that rather than treating it as success, because "changed without a
    backup" and "changed with one" are different things to tell a user
    about their records.
    """
    path = Path(path)
    if not path.is_file():
        return None
    folder = Path(folder) if folder is not None else backup_dir_for(path)
    try:
        folder.mkdir(parents=True, exist_ok=True)
    except OSError:
        return None
    target = versioned_backup_path(path, folder)
    if target is None or target.exists():
        return None
    try:
        shutil.copy2(str(path), str(target))
    except OSError:
        return None
    return target


def list_backups(path, folder=None) -> list[Path]:
    """Existing copies of *path*, oldest first."""
    path = Path(path)
    folder = Path(folder) if folder is not None else backup_dir_for(path)
    if not folder.is_dir():
        return []
    found = []
    for candidate in folder.iterdir():
        if not candidate.is_file() or candidate.suffix != path.suffix:
            continue
        stem = candidate.stem
        if stem == f"{path.stem}.bak":
            found.append((1, candidate))
        elif stem.startswith(f"{path.stem}.bak"):
            tail = stem[len(f"{path.stem}.bak"):]
            if tail.isdigit():
                found.append((int(tail), candidate))
    return [c for _, c in sorted(found)]
