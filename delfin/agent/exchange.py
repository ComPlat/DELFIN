"""One directory where sessions working on the same repository exchange files.

Named `exchange` and not `handover`: that module already exists and means
something else -- handing a refused COMMAND to the user to run. This is
about files between sessions. The user-facing wording stays "handover
directory", which is what it is called in the refusal that offers it.

Three sessions worked together on one repository, each in its own
worktree, and exchanged FILENAMES. Every handover failed:

    File not found: reports/20261007-chemdarwin-design-loop-blueprint-….md

-- the report existed, in the author's worktree, which is not a path the
reader may open. The containment was right to refuse: a session that can
read another session's worktree can read its whole workspace, and the
boundary is the only thing that makes parallel sessions safe to run.

So the answer is not a wider boundary. It is ONE directory, named, inside
the repository's own ``.delfin``, shared by the repository and every
worktree of it, and nothing else.

Anchored on git's COMMON directory rather than the checkout root: the
root differs per worktree (that is what a worktree is), the common
directory does not. Three sessions in three worktrees of one repository
therefore agree on one path without being told it, and sessions in a
different repository never see it.

What this does NOT solve, stated because the field report shows the other
half: `chemdarwin_interface.py` was missing from a session's worktree,
and that is CODE. Merging code between sessions is a merge, not a drop
-- this is for the documents sessions write for each other (reports,
schemas, specifications, CSVs).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

#: The leaf. A name, not a dot-directory: it is meant to be found, and a
#: session that is told "put it in the handover directory" has to be able
#: to see it in a listing.
_LEAF = "handover"

#: Where it lives relative to the repository that owns it. Inside
#: ``.delfin`` because that directory is already DELFIN's and already
#: holds the worktrees these sessions run in; a new top-level directory
#: in somebody's repository is not ours to create.
_PARENT = ".delfin"


def directory_for(workspace: str | Path, *, create: bool = True) -> Optional[Path]:
    """The handover directory for the repository *workspace* belongs to.

    Input: any directory inside a checkout or one of its worktrees.
    Output: ``<repository>/.delfin/handover``, or None outside a
    repository and wherever the path cannot be established.

    Semantics: the SAME path for every worktree of one repository, and a
    different one for a different repository. ``create`` makes it, owner
    only; False answers where it would be without touching the disk.

    Never raises. A session must start on a host where this cannot be
    made -- it then simply has no shared directory, which is the state
    every session was in before.
    """
    try:
        from . import session_presence as _presence
        info = _presence.repository_of(str(workspace or ""))
        common = str(info.get("common_dir") or "")
        if not common:
            return None
        owner = Path(common).parent
        if not owner.is_dir():
            return None
        target = owner / _PARENT / _LEAF
        # Not the repository, not .delfin itself. Those hold the
        # worktrees, the session store and the settings; handing them
        # over as a shared root would hand over everything this is
        # supposed to keep separate. Checked rather than assumed,
        # because the whole value of this module is that the answer is
        # one leaf and not a parent of anything.
        if target == owner or target == owner / _PARENT:
            return None
        if create:
            target.mkdir(parents=True, exist_ok=True)
            try:
                os.chmod(target, 0o700)
            except OSError:
                pass
        return target.resolve() if target.exists() else target
    except Exception:
        return None


def contains(path: str | Path, workspace: str | Path) -> bool:
    """Whether *path* really lies in this workspace's handover directory.

    Resolved on both sides before comparing, so a symlink placed inside
    the directory is not a way out of it: the shared root is the one
    place two sessions can both write, which makes it the one place to
    plant one.

    Never raises; an unresolvable path is not inside.
    """
    try:
        room = directory_for(workspace, create=False)
        if room is None:
            return False
        here = Path(path).expanduser().resolve()
        room = room.resolve() if room.exists() else room
        return here == room or room in here.parents
    except Exception:
        return False


def describe(workspace: str | Path) -> str:
    """One line for the model, or "" where there is no such directory.

    The agent cannot use a directory it has not been told about, and the
    field report is three sessions failing to hand anything over while
    the boundary refusing them was working exactly as designed.
    """
    room = directory_for(workspace, create=False)
    if room is None:
        return ""
    return (f"Shared handover directory: {room} — readable and writable by "
            "every session working on this repository, including sessions "
            "in other worktrees. Put reports, schemas and data you want "
            "another session to read there, and tell it the path. It is "
            "the ONLY path outside your own workspace that another "
            "session can reach; their worktrees are not readable.")


#: The setting that turns it off. On by default: a shared directory that
#: has to be switched on is one nobody switches on, and the field report
#: is sessions failing to hand anything over at all.
_SETTING = "handover_dir"


def enabled() -> bool:
    """Whether the shared directory is granted. Never raises; on by
    default, including where the settings file cannot be read."""
    try:
        from delfin.user_settings import load_settings
        value = ((load_settings() or {}).get("agent") or {}).get(_SETTING, True)
    except Exception:
        return True
    if isinstance(value, str):
        return value.strip().lower() not in ("0", "false", "off", "no")
    return bool(value)


def with_exchange(workspace: str | Path, extra) -> tuple:
    """*extra* workspace roots, with this repository's handover directory.

    Input: the session's workspace and the extra roots it already has.
    Output: the same tuple plus the handover directory, once. Semantics:
    unchanged outside a repository, with the setting off, or where the
    directory is already granted.

    One function with two callers, because the terminal path and the KIT
    path each assemble their own extra-roots tuple and a second copy of
    this decision is how one surface comes to have a shared directory
    and the other not.
    """
    current = tuple(extra or ())
    if not enabled():
        return current
    # create=False, and that is not an optimisation. Granting a root runs
    # whenever a session's permissions are built -- which the test suite
    # does thousands of times, and which a `delfin-agent doctor` does
    # once. Creating here put a directory into the real checkout as a
    # side effect of constructing an object: it appeared in the
    # repository, the suite's guard against writing into the checkout saw
    # it, and 14 unrelated tests failed. The directory is made when
    # something writes to it.
    room = directory_for(workspace, create=False)
    if room is None:
        return current
    resolved = {Path(p).expanduser().resolve() for p in current
                if str(p).strip()}
    if room in resolved:
        return current
    # Said once per process, because a root outside the workspace is a
    # door and the security panel is where the user sees which doors are
    # open.
    _announce(room)
    return current + (room,)


_ANNOUNCED: set = set()


def _announce(room: Path) -> None:
    """Record that the shared root is active. Never raises."""
    key = str(room)
    if key in _ANNOUNCED:
        return
    _ANNOUNCED.add(key)
    try:
        from .api_client import _record_security_event
        _record_security_event(
            "handover_root", "workspace",
            f"{room} is readable and writable by every session working on "
            "this repository (agent.handover_dir)", blocked=False)
    except Exception:
        pass

