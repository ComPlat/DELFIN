"""The agent's shell command inside a container (``agent.bash_isolation = "container"``).

Settings (user file only; the agent's write gates refuse it and no tool
changes it):

``agent.bash_isolation = "container"``
    selects this backend for ``bash``, ``bash_background`` and
    ``run_tests`` (all three build their argv through
    ``api_client._bash_isolation_argv``).
``agent.container.runtime``
    ``auto`` (default) | ``apptainer`` | ``singularity`` | ``podman`` |
    ``docker``. ``auto`` takes the first of ``RUNTIMES`` that is installed
    and has the image. A named runtime is never swapped for another.
``agent.container.image``
    OCI runtimes: an image reference present locally (never pulled at run
    time, ``--pull=never``); default ``DEFAULT_OCI_IMAGE``. Apptainer: a
    ``.sif`` file or sandbox directory, or a URI (``docker://...``) the
    runtime resolves itself; default ``default_sif()``.

Output: ``session_argv(cmd, run_cwd, perms, extra_write) -> list | None``.
A list is the full argv: the session's writable roots (workspace roots,
``extra_write``, the git metadata directories) bound read-write; read-only
and confirm-write roots, git hooks and approval answer rooms bound
read-only; credential directories of ``sandbox._HOME_SECRET_DIRS`` that lie
inside a bound root covered by an empty read-only directory, credential
files by ``/dev/null``; nothing else of the host's filesystem. Network
none, unless ``agent.sandbox_network.mode`` is ``open``. The user's uid
(Apptainer by default, Podman ``--userns=keep-id``, Docker ``--user``),
no capabilities, no-new-privileges, a read-only root filesystem (OCI
``--read-only``; a SIF image is read-only by construction).

``None`` means no container can run (no runtime, no image, unknown
runtime name): the caller falls back to bubblewrap/Landlock/Seatbelt and
a notice is recorded once. It is never an open shell and never a broken
run.

Rejected alternatives:
- the egress proxy inside the container. The proxy is a Unix socket and
  the existing forwarder runs in the supervising process's network
  namespace; a container with its own namespace would need a forwarder
  inside the image. ``proxy`` therefore means network none here, and the
  notice says so.
- wrapping the runtime in the bubblewrap process cage: a runtime needs
  user namespaces of its own and does not start inside the cage.
- pulling a missing image on first use: a first command would wait for a
  download of unknown size, from a registry the network policy may not
  allow. The image is a setup step the user does once.

Cost: one ``shutil.which`` per runtime and, for an OCI runtime, one
``<runtime> image inspect`` per process (cached); nothing for a session
whose setting is not ``container``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

RUNTIMES = ("apptainer", "singularity", "podman", "docker")
_APPTAINER = ("apptainer", "singularity")
DEFAULT_OCI_IMAGE = "docker.io/library/python:3.11-slim"
_SIF_SOURCE = "docker://python:3.11-slim"


def default_sif() -> Path:
    """Where the Apptainer image is looked for when none is configured."""
    return Path.home() / ".local" / "share" / "delfin" / "containers" / "python_3.11-slim.sif"


@dataclass(frozen=True)
class Choice:
    runtime: str        # "" when no container can run
    binary: str
    image: str
    reason: str         # why not, when runtime is ""


_CHOICE: Choice | None = None
_ANNOUNCED: set[str] = set()
_EMPTY_DIR: str = ""


def reset() -> None:
    """Forget the cached choice and the once-per-process notices."""
    global _CHOICE
    _CHOICE = None
    _ANNOUNCED.clear()


def _settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        cfg = ((load_settings() or {}).get("agent") or {}).get("container") or {}
        return dict(cfg) if isinstance(cfg, dict) else {}
    except Exception:
        return {}


def _configured_mode() -> str:
    """The isolation mode in force: the per-process override, else the setting."""
    try:
        from . import api_client as ac
        if ac._BASH_ISOLATION_OVERRIDE:
            return ac._BASH_ISOLATION_OVERRIDE.strip().lower()
    except Exception:
        pass
    try:
        from delfin.user_settings import load_settings
        return str(((load_settings() or {}).get("agent") or {})
                   .get("bash_isolation", "auto") or "auto").strip().lower()
    except Exception:
        return "auto"


def _image_present(runtime: str, binary: str, image: str) -> bool:
    if runtime in _APPTAINER:
        if "://" in image:
            return True
        return Path(image).expanduser().exists()
    try:
        done = subprocess.run([binary, "image", "inspect", image],
                              capture_output=True, timeout=20)
        return done.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def _missing_image_reason(runtime: str, image: str) -> str:
    if runtime in _APPTAINER:
        return (f"the image {image} is not a .sif file on this host (build it "
                f"once: {runtime} pull {image} {_SIF_SOURCE})")
    return (f"the image {image} is not present for {runtime}, or its daemon "
            f"does not answer (pull it once: {runtime} pull {image})")


def choose() -> Choice:
    """The runtime and image a command would use, cached per process."""
    global _CHOICE
    if _CHOICE is not None:
        return _CHOICE
    cfg = _settings()
    pref = str(cfg.get("runtime") or "auto").strip().lower()
    if pref != "auto" and pref not in RUNTIMES:
        _CHOICE = Choice("", "", "", f"agent.container.runtime {pref!r} is not "
                         f"one of auto, {', '.join(RUNTIMES)}")
        return _CHOICE
    names = RUNTIMES if pref == "auto" else (pref,)
    found = [(n, shutil.which(n)) for n in names]
    found = [(n, b) for n, b in found if b]
    if not found:
        what = ", ".join(names)
        _CHOICE = Choice("", "", "", f"no container runtime ({what}) is installed")
        return _CHOICE
    configured = str(cfg.get("image") or "").strip()
    first_reason = ""
    for name, binary in found:
        image = configured or (str(default_sif()) if name in _APPTAINER
                               else DEFAULT_OCI_IMAGE)
        if _image_present(name, binary, image):
            _CHOICE = Choice(name, binary, image, "")
            return _CHOICE
        first_reason = first_reason or _missing_image_reason(name, image)
    _CHOICE = Choice("", "", "", first_reason)
    return _CHOICE


def _empty_dir() -> str:
    """An empty, read-only-mounted source that covers a secret directory."""
    global _EMPTY_DIR
    if not _EMPTY_DIR or not os.path.isdir(_EMPTY_DIR):
        _EMPTY_DIR = tempfile.mkdtemp(prefix="delfin-empty-")
        os.chmod(_EMPTY_DIR, 0o500)
    return _EMPTY_DIR


def build_argv(runtime: str, binary: str, image: str, command: list[str], *,
               cwd: str, write_roots, read_roots, mask_dirs=(), mask_files=(),
               network: str = "none", uid: int | None = None,
               gid: int | None = None) -> list[str]:
    """The runtime's argv for ``command``. Pure apart from the existence
    checks: a root that does not exist is not mounted (Docker would create
    it as root, Podman would refuse to start)."""
    uid = os.getuid() if uid is None else uid
    gid = os.getgid() if gid is None else gid
    seen: set[str] = set()
    mounts: list[tuple[str, str, str]] = []
    for mode, roots in (("rw", write_roots), ("ro", read_roots)):
        for r in roots:
            r = str(r)
            if r in seen or not os.path.exists(r):
                continue
            seen.add(r)
            mounts.append((r, r, mode))
    if mask_dirs:
        empty = _empty_dir()
        mounts += [(empty, str(d), "ro") for d in mask_dirs]
    mounts += [("/dev/null", str(f), "ro") for f in mask_files]
    closed = network != "open"

    if runtime in _APPTAINER:
        argv = [binary, "exec", "--containall", "--cleanenv", "--no-privs",
                "--no-mount", "hostfs,bind-paths", "--pwd", cwd]
        if closed:
            argv += ["--net", "--network", "none"]
        for src, dst, mode in mounts:
            argv += ["--bind", f"{src}:{dst}:{mode}"]
        return argv + [image, *command]

    argv = [binary, "run", "--rm", "--pull=never", "--read-only",
            "--cap-drop=ALL", "--security-opt=no-new-privileges",
            "--env", "HOME=/tmp", "--workdir", cwd]
    if runtime == "podman":
        # keep-id maps the user's uid into the container unchanged, so a
        # file written in a root belongs to the user on the host; the
        # read-only root comes with tmpfs on /tmp, /run and /var/tmp.
        argv += ["--userns=keep-id"]
    else:
        argv += ["--init", "--user", f"{uid}:{gid}",
                 "--tmpfs", "/tmp:rw,nosuid,nodev,mode=1777"]
    if closed:
        argv += ["--network=none"]
    for src, dst, mode in mounts:
        argv += ["-v", f"{src}:{dst}:{mode}"]
    return argv + [image, *command]


def _announce(key: str, detail: str) -> None:
    if key in _ANNOUNCED:
        return
    _ANNOUNCED.add(key)
    try:
        from . import api_client as ac
        ac._record_security_event("isolation", "bash", detail, blocked=False)
    except Exception:
        pass


def _inside(path: Path, roots) -> bool:
    return any(path == r or r in path.parents for r in roots)


def session_argv(cmd: str, run_cwd, perms, extra_write=()) -> list[str] | None:
    """The container argv for one command of a session, or None (see module)."""
    choice = choose()
    if not choice.runtime:
        _announce("fallback", (
            "container isolation was requested (agent.bash_isolation = "
            f"\"container\") but {choice.reason}; commands fall back to "
            "bubblewrap, Landlock or Seatbelt, whichever holds here"))
        return None
    from . import api_client as ac

    def _resolved(items) -> list[Path]:
        out = []
        for p in items:
            try:
                out.append(Path(p).resolve())
            except (OSError, RuntimeError):
                continue
        return out

    try:
        writable = _resolved(perms.all_workspace_roots())
    except Exception:
        writable = _resolved([getattr(perms, "workspace", ".")])
    writable += _resolved(extra_write)
    git_roots = _resolved(ac._git_metadata_roots(getattr(perms, "workspace", ".")))
    writable += git_roots
    try:
        readable = _resolved(perms.all_readable_roots())
    except Exception:
        readable = []
    read_only = [r for r in readable if r not in writable]
    # Nested read-only mounts inside a writable root, after it: git hooks
    # (a planted hook would run on the user's next git command) and the
    # rooms where approvals are answered.
    read_only += [g / "hooks" for g in git_roots]
    try:
        from . import approval_guard as _ag
        read_only += [room for room in _resolved(_ag.answer_rooms())
                      if _inside(room, writable)]
    except Exception:
        pass
    bound = writable + read_only
    mask_dirs, mask_files = [], []
    try:
        from .sandbox import _HOME_SECRET_DIRS
    except Exception:
        _HOME_SECRET_DIRS = ()
    for rel in _HOME_SECRET_DIRS:
        try:
            p = (Path.home() / rel).resolve()
        except (OSError, RuntimeError):
            continue
        if not _inside(p, bound):
            continue
        if p.is_dir():
            mask_dirs.append(str(p))
        elif p.is_file():
            mask_files.append(str(p))
    try:
        from . import egress_proxy as _ep
        net_setting = _ep.network_mode()
    except Exception:
        net_setting = "none"
    network = "open" if net_setting == "open" else "none"
    note = ""
    if net_setting == "proxy":
        note = ("; the egress proxy is not reachable from inside a container, "
                "so the allowed domains are not either (set "
                "agent.sandbox_network.mode = \"open\" to give it the network)")
    _announce("active", (
        f"commands run in a {choice.runtime} container (image {choice.image}): "
        "only the session's roots are mounted, as the user's uid, without "
        f"capabilities; network {network}{note}"))
    return build_argv(choice.runtime, choice.binary, choice.image,
                      ["/bin/bash", "-c", cmd],
                      cwd=str(Path(run_cwd).resolve()),
                      write_roots=[str(p) for p in writable],
                      read_roots=[str(p) for p in read_only],
                      mask_dirs=mask_dirs, mask_files=mask_files,
                      network=network)
