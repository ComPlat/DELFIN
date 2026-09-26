"""One temp directory shared between the bash cage and the file tools.

Today the two see different ``/tmp`` directories: bash runs inside
bwrap with ``--tmpfs /tmp`` (mcp_isolation.bwrap_argv), a private empty
tmpfs, while read_file/write_file resolve absolute paths against the
real filesystem. A file a command writes to /tmp is therefore invisible
to read_file -- and two Welle-4 sessions bridged that by asking for
read access to ALL of /tmp, which the read gate rightly refused.

The design this module computes for:

- One temp directory per session, under DELFIN's state root:
  ``<state>/tmp/<sanitized-session-id>/``. It is NOT under the
  workspace, so scratch files never land in the user's project nor in
  git, and NOT under ``/tmp``, so nothing outside the mapping ever
  sees host /tmp.
- The cage binds that directory AS ``/tmp`` (``--bind <dir> /tmp``
  replacing ``--tmpfs /tmp``), so a command's /tmp and the file tools'
  mapped path are the same files.
- The file tools map any path under ``/tmp`` to the session directory
  before their normal path checks run; every other path passes through
  untouched, so nothing else changes behaviour.

What is allowed to live there, who cleans up, and why nothing can
escape -- see the security notes in .gate/TMP-SICHERHEIT.md (the
operator's build decision); the mapping rules themselves are pinned
by tests/test_shared_tmp.py.
"""

from __future__ import annotations

import os
import re
from pathlib import Path, PurePosixPath

# The cage-side mount point this module maps.
CAGE_TMP = "/tmp"

# Session ids become directory names. Keep it conservative: letters,
# digits, dash, underscore, dot; everything else becomes a dash.
_SAFE_ID = re.compile(r"[^A-Za-z0-9._-]+")

_DEFAULT_STATE_ROOT: Path | None = None


def state_root() -> Path:
    """The DELFIN state root shared temp directories live under.

    ``~/.delfin`` by default, overridable via the ``DELFIN_STATE``
    environment variable (used by tests and by callers that keep state
    elsewhere). Resolution errors fall back to ``~/.delfin`` -- a
    missing HOME must not turn a temp-path question into a crash.
    """
    env = os.environ.get("DELFIN_STATE", "").strip()
    if env:
        return Path(env).expanduser()
    try:
        return Path.home() / ".delfin"
    except Exception:  # pragma: no cover - HOME unset
        return Path("/.delfin")


def _is_safe_dir(p: Path) -> bool:
    """A candidate location must be an existing directory reached
    through real directories only: every component is lstat'd and must
    be a plain directory, never a symlink or a file. A pre-placed
    symlink anywhere on the way disqualifies the candidate."""
    try:
        cur = Path(p.anchor) if p.is_absolute() else Path(".")
        for part in (p.parts[1:] if p.is_absolute() else p.parts):
            cur = cur / part
            if (cur.lstat().st_mode & 0o170000) != 0o040000:
                return False  # symlink, file, or anything but a plain dir
        return True
    except OSError:
        return False


def location_candidates(uid: int | None = None) -> list[Path]:
    """Where the shared temp tree may live, best first.

    Chain (reasoning in .gate/ORT.md):
    1. ``$DELFIN_TMP_ROOT`` -- explicit operator override, always wins.
    2. ``$XDG_RUNTIME_DIR`` -- per-user tmpfs, wiped at logout, no HOME
       quota (the operator's two objections to the HOME location).
    3. ``/run/user/<uid>`` -- the same directory when XDG is not
       exported (recomputed from the real uid when not given).
    4. The state root (``$DELFIN_STATE`` / ``~/.delfin``) -- always
       present, last resort.

    A candidate counts only if it is a real directory (no symlinks, no
    files) that already exists; the caller does not create these roots.
    The host's ``/tmp`` and ``$TMPDIR`` are never candidates -- on this
    cluster TMPDIR is a shared scratch filesystem, exactly what the
    operator rejected.
    """
    candidates: list[Path] = []
    env_root = os.environ.get("DELFIN_TMP_ROOT", "").strip()
    if env_root:
        # Explicit operator override wins unconditionally: it need not
        # exist yet -- Phase 2 creates it with the pre-placement checks
        # (owner, 0700, no symlinks). Everything below is only for the
        # implicit chain.
        return [Path(env_root)]
    xdg = os.environ.get("XDG_RUNTIME_DIR", "").strip()
    if xdg:
        p = Path(xdg)
        if p not in candidates:
            candidates.append(p)
    real_uid = os.getuid() if uid is None else uid
    run_user = Path(f"/run/user/{real_uid}")
    if run_user not in candidates:
        candidates.append(run_user)
    sr = state_root()
    if sr not in candidates:
        candidates.append(sr)
    return [c for c in candidates if _is_safe_dir(c)] or [sr]


def resolve_location(uid: int | None = None) -> Path:
    """The first usable location from :func:`location_candidates`."""
    return location_candidates(uid=uid)[0]


def session_dir(session_id: str, base: Path | None = None) -> Path:
    """This session's shared temp directory, under ``<base>/tmp/``.

    The id is sanitized to a single path component: any run of
    unsafe characters collapses to one dash, so
    ``session_dir("../../../etc")`` is a directory NEXT TO the other
    session directories, never an escape out of ``<base>/tmp``.
    Deterministic for the same id -- the cage and the file tools must
    agree on the directory without communicating.

    ``base`` defaults to :func:`resolve_location` -- the XDG runtime
    dir first, the state root only as a fallback (see .gate/ORT.md).
    """
    root = (Path(base) if base is not None else resolve_location()) / "tmp"
    safe = _SAFE_ID.sub("-", str(session_id)).strip("-.") or "session"
    return root / safe


def map_to_host(path: str, session_dir: Path) -> Path:
    """A cage-side path as the file tools must resolve it.

    Rules, in order:
    1. Normalize ``..`` and ``.`` FIRST (``os.path.normpath``), so a
       path that climbs out of /tmp is judged after it does -- the
       alternative would join ``<session_dir>/../..`` and hand the
       tools a path OUTSIDE the share.
    2. A path that (after normalization) is ``/tmp`` or under it maps
       componentwise into the session directory.
    3. Everything else -- other absolute paths, relative paths --
       passes through unchanged: the caller's normal resolution and
       permission checks stay exactly as they are. A normalized
       ``/tmp/../etc/passwd`` IS ``/etc/passwd``, so it passes through
       and the read gate judges it as the host path it is.
    """
    normalized = os.path.normpath(str(path))
    pure = PurePosixPath(normalized)
    cage_root = PurePosixPath(CAGE_TMP)
    if pure.is_absolute() and (pure == cage_root
                               or cage_root in pure.parents):
        remainder = pure.relative_to(CAGE_TMP)
        target = Path(session_dir)
        for part in remainder.parts:
            if part in ("..", "."):
                continue  # normpath already removed these; belt+braces
            target = target / part
        return target
    return Path(normalized)


def cage_mount(session_dir: Path) -> Path:
    """The directory the cage should bind as ``/tmp``.

    Trivial today, but it is the one place that names the contract:
    ``bwrap_argv`` replaces ``["--tmpfs", "/tmp"]`` with
    ``["--bind", str(cage_mount(...)), "/tmp"]``.
    """
    return Path(session_dir)


def bwrap_substitution(session_dir: Path) -> tuple[list[str], list[str]]:
    """(old argv fragment, new argv fragment) for ``bwrap_argv``.

    Returned as the exact lists so the cage build can assert the swap
    rather than string-splice it.
    """
    return (["--tmpfs", CAGE_TMP],
            ["--bind", str(Path(session_dir)), CAGE_TMP])
