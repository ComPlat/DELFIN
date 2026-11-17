"""Per-workspace freedom profiles: how far a workspace's unattended runs go.

The swarm side (merged #493, ``api_client._bash_runs_caged`` /
``_caged_still_asked``) already lets an attended Bypass session run
workspace work inside the bash cage free: the cage confines its writes
and holds its network, and push, git remote, installers outside the
project env, listeners and anything uncaged still ask. This module
brings that model to scheduled runs: a run starts in bypass only when
the host can actually hold a command in the cage AND the schedule's
workspace carries the profile ``free``. Every other workspace keeps the
behaviour before this change (``scheduler_daemon.scheduled_permission_mode``
caps bypass to acceptEdits).

The store is the user's alone, like the schedules
(``delfin.agent.schedules``): a 0700 directory with 0600 files beside
them, masked in every cage an agent's command runs in and refused to the
file tools, because the two stores share the directory. A profile is
confirmed with the account key (``trust_store``, HMAC over its
fingerprint); an unconfirmed profile reads as ``default``, so an agent
that manages to write the file gains nothing.

The host-cage decision is read from ``api_client._bash_isolation_decision``
-- the same answer the bash executor builds its argv from, so the runner
and the gate cannot drift. It is exposed here as :func:`host_cage_decision`
(a module attribute, never a local import at call time) so a test can
patch one name and the runner, the store and the doctor all follow.

Rejected: a per-schedule flag instead of per-workspace profiles (the
freedom belongs to the folder a run acts in, not to one entry), and a
global permission_mode rewrite (that would widen every unattended run,
not only a workspace the user granted).

Never raises at call time: :func:`get_profile` and
:func:`host_cage_decision` answer a value on any failure.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

PROFILES = ("default", "free")

_DIR = Path.home() / ".delfin" / "user_schedules"

_FILE_NAME = "workspace_profiles.json"

_MAC_PURPOSE = "workspace_profiles"


def host_cage_decision() -> str:
    """The kind of confinement shell commands get on this host.

    The answer of ``api_client._bash_isolation_decision`` with this
    process's settings and default permissions: "bwrap", "landlock" or
    "seatbelt" confine writes; "none" (nothing here can) and "refuse"
    (isolation is demanded and unavailable -- commands do not run) do
    not, and under them the profile ``free`` stays refused.
    """
    try:
        from .api_client import _bash_isolation_decision
        kind, _mode = _bash_isolation_decision(type("P", (), {
            "scope_locked": False})())
        return str(kind)
    except Exception:
        return "none"


def host_can_cage() -> bool:
    """Whether this host can hold a shell command in the cage."""
    return host_cage_decision() in ("bwrap", "landlock", "seatbelt")


def _file() -> Path:
    return _DIR / _FILE_NAME


def _ensure_dir() -> None:
    _DIR.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        os.chmod(_DIR, 0o700)
    except OSError:
        pass


@contextmanager
def _locked() -> Iterator[None]:
    _ensure_dir()
    fd = os.open(_DIR / ".lock", os.O_CREAT | os.O_RDWR, 0o600)
    try:
        import fcntl
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        os.close(fd)


def _read() -> dict:
    try:
        data = json.loads(_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"profiles": {}}
    if not isinstance(data, dict) or not isinstance(data.get("profiles"), dict):
        return {"profiles": {}}
    return data


def _write(data: dict) -> None:
    from .state_paths import write_text_atomic
    try:
        os.chmod(_file(), 0o600)
    except OSError:
        pass
    write_text_atomic(_file(), json.dumps(data, indent=1, ensure_ascii=False))


def _normalize_workspace(workspace: str) -> str:
    return str(Path(str(workspace or "")).expanduser().resolve())


def _record(workspace: str, profile: str, confirmed_sha: str = "") -> dict:
    return {"profile": profile, "confirmed_sha": confirmed_sha,
            "ts": time.time()}


def fingerprint(workspace: str, profile: str) -> str:
    """What the user confirms: which folder runs how."""
    keep = [workspace, profile]
    return hashlib.sha256(json.dumps(keep, ensure_ascii=False).encode()).hexdigest()


def _seal(workspace: str, profile: str) -> str:
    """The confirmation of this profile. Refused inside an agent
    session: only the user sets a workspace free."""
    from . import trust_store
    from .approval_guard import inside_agent_session
    if inside_agent_session():
        raise PermissionError("refused inside an agent session: only the "
                              "user sets a workspace profile")
    tag = trust_store.mac(_MAC_PURPOSE,
                          fingerprint(workspace, profile).encode())
    return tag or "sha:" + fingerprint(workspace, profile)


def confirmed(workspace: str) -> bool:
    """Whether the user confirmed exactly this profile: the HMAC under
    the account key, or -- only where there is no account key at all --
    the plain hash."""
    data = _read()
    record = (data.get("profiles") or {}).get(_normalize_workspace(workspace))
    if not isinstance(record, dict):
        return False
    profile = str(record.get("profile") or "")
    tag = str(record.get("confirmed_sha") or "")
    if not tag:
        return False
    from . import trust_store
    if tag.startswith("sha:"):
        return (tag == "sha:" + fingerprint(_normalize_workspace(workspace),
                                            profile)
                and trust_store.account_key(create=True) is None)
    return trust_store.check_mac(_MAC_PURPOSE,
                                 fingerprint(_normalize_workspace(workspace),
                                             profile).encode(),
                                 tag) is True


def get_profile(workspace: str) -> str:
    """The workspace's profile, ``default`` unless the user confirmed
    ``free`` (an unconfirmed or broken record stays ``default``)."""
    data = _read()
    record = (data.get("profiles") or {}).get(_normalize_workspace(workspace))
    if not isinstance(record, dict):
        return "default"
    profile = str(record.get("profile") or "default")
    if profile not in PROFILES or profile == "default":
        return "default"
    try:
        return profile if confirmed(workspace) else "default"
    except Exception:
        return "default"


def set_profile(workspace: str, profile: str) -> dict:
    """The user sets a workspace's profile. ``free`` is refused on a
    host that cannot cage, with the reason."""
    profile = str(profile or "").strip().lower()
    if profile not in PROFILES:
        raise ValueError(f"profile is one of {', '.join(PROFILES)}")
    ws = _normalize_workspace(workspace)
    if profile == "free" and not host_can_cage():
        raise PermissionError(
            f"this host cannot hold a shell command in the cage "
            f"(bash cage: {host_cage_decision()}); the profile \"free\" "
            f"is refused here")
    sha = _seal(ws, profile)
    with _locked():
        data = _read()
        data["profiles"][ws] = _record(ws, profile, sha)
        _write(data)
    return {"workspace": ws, "profile": profile, "confirmed": True}


def set_from_word(workspace: str, word: str) -> dict:
    """The messenger's and TUI's word for a profile: ``frei`` (or
    ``free``) sets ``free``, ``standard`` (or ``default``) resets."""
    word = str(word or "").strip().lower().lstrip("/")
    mapping = {"frei": "free", "free": "free", "standard": "default",
               "default": "default"}
    if word not in mapping:
        raise ValueError("profil <ordner> frei|standard")
    return set_profile(workspace, mapping[word])


def set_for_workspace_words(workspace: str, word: str) -> dict:
    """Like :func:`set_from_word`, but only the user (outside an agent
    session) may raise the profile to ``free``."""
    return set_from_word(workspace, word)
