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
import stat
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


class PrePlacedError(RuntimeError):
    """The session temp directory (or a path component of it) already
    exists but is not OURS: wrong owner, wrong permissions, a symlink,
    or a regular file. Raised by :func:`ensure_session_dir` -- the
    pre-placed object is never adopted, and the caller must not fall
    back to "just use it anyway"."""


def _lstat_mode(p: Path) -> int | None:
    """The lstat mode of ``p`` (does not follow symlinks), or None."""
    try:
        return os.lstat(p).st_mode
    except OSError:
        return None


def ensure_session_dir(session_dir_path: Path) -> Path:
    """Create -- or safely adopt -- one session's shared temp directory.

    The pre-placement attack this refuses: another user creates the
    directory (or a symlink at, or anywhere on the way to, that path)
    before we do, and thereby reads what the cage writes or points it
    outside. session_dir() always has the shape ``<base>/tmp/<safe>``,
    so the rules split at OUR namespace boundary:

    - EVERY component of the path is lstat'd and must be a plain
      directory -- a symlink or a file anywhere in the path is refused,
      never followed (the O_NOFOLLOW equivalent).
    - The last two components -- the ``tmp`` collection directory and
      the session directory itself, the parts DELFIN owns -- must, if
      they already exist, be owned by the current user with mode
      exactly 0700; anything looser is a pre-placed directory and is
      refused, untouched. If missing they are created with
      os.mkdir(..., 0o700) plus an explicit chmod (mkdir masks with the
      umask, so the chmod pins the exact mode).
    - Components above that boundary are provided by the environment
      (``/run/user/<uid>``, a private cage ``/tmp``) and may legitimately
      be 0755 or root-owned; only the no-symlink rule applies to them --
      we cannot refuse a system directory, but we can refuse to walk
      through anything that is not a real directory.

    Idempotent: a directory this function created (or one already
    meeting all the rules) passes unchanged on the second call.
    """
    p = Path(session_dir_path)
    parts = p.parts
    if not p.is_absolute():
        raise PrePlacedError(
            f"session temp dir must be absolute, got {p!s}")
    strict_count = 2  # the 'tmp' collection dir + the session dir itself
    cur = Path(parts[0])
    for idx, part in enumerate(parts[1:], start=1):
        cur = cur / part
        mode = _lstat_mode(cur)
        if mode is None:
            os.mkdir(cur, 0o700)
            os.chmod(cur, 0o700)  # umask-proof: pin the exact mode
            continue
        if not stat.S_ISDIR(mode):
            raise PrePlacedError(
                f"{cur!s} exists but is not a plain directory (symlink "
                f"or file): refusing to use or create beneath it")
        if idx > len(parts) - 1 - strict_count:
            # our namespace: 'tmp' collection dir and session dir
            st = os.lstat(cur)
            if st.st_uid != os.getuid():
                raise PrePlacedError(
                    f"{cur!s} is owned by uid {st.st_uid}, not "
                    f"{os.getuid()}: refusing to use a foreign directory")
            if stat.S_IMODE(mode) != 0o700:
                raise PrePlacedError(
                    f"{cur!s} has mode {oct(stat.S_IMODE(mode))}, not "
                    f"0700: refusing to adopt a pre-placed directory")
    return p


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
