"""Who a phone can belong to: the account owner or one person of people mode.

People mode (the Windows launcher, ``~/.config/delfin/launcher-people.conf``
in the account's home) gives each person ``~/people/<name>`` as HOME: own
sessions, agent service (``<home>/.run``), messenger activation records and
pairing requests. Everyone is the same Unix user, so the separation drawn
here is an application-level boundary: the gateway and the bridge never
route one person's phone into another person's sessions, but a process of
the account can read every person's files.

Input: a person name as the bridge reports it (the directory a pairing
request was found in). Output: that person's home, or None when the name is
not a person of this account. The owner is the empty name ``""``: the
account's own home and the sessions outside ``~/people``.
"""
from __future__ import annotations

import os
import re
import stat
from pathlib import Path

#: The account owner: the person an admin pairs with ``pi_messenger.sh --pair``.
OWNER = ""

#: The launcher's person names (tools/windows-launcher/remote_bootstrap.sh:
#: lower-case letters, digits, ``_`` and ``-``, not starting with ``-``).
_NAME = re.compile(r"[a-z0-9_][a-z0-9_-]{0,31}")

PEOPLE_CONF = Path(".config") / "delfin" / "launcher-people.conf"


def valid(name: object) -> bool:
    return isinstance(name, str) and (name == OWNER or bool(_NAME.fullmatch(name)))


def account_home() -> Path:
    """The account's home: the HOME of the bridge (the service user)."""
    return Path.home()


def mode_on(account: Path | None = None) -> bool:
    return ((account or account_home()) / PEOPLE_CONF).is_file()


def home(name: str, account: Path | None = None, *, uid: int | None = None) -> Path | None:
    """The person's home, or None. A person's home must be a real directory
    (never a link) of ``uid`` (default: this process's) under
    ``<account>/people``, and people mode must be on."""
    account = account or account_home()
    uid = os.getuid() if uid is None else uid
    if name == OWNER:
        return account
    if not valid(name) or not mode_on(account):
        return None
    path = account / "people" / name
    try:
        info = os.lstat(path)
    except OSError:
        return None
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != uid:
        return None
    return path


def names(account: Path | None = None) -> list[str]:
    """The owner, then every person of people mode, sorted."""
    account = account or account_home()
    found = [OWNER]
    if mode_on(account):
        try:
            entries = sorted(os.listdir(account / "people"))
        except OSError:
            entries = []
        found += [n for n in entries if n != OWNER and home(n, account) is not None]
    return found


def runtime_socket(name: str, account: Path | None = None) -> Path | None:
    """A person's agent-service socket (``<home>/.run/delfin-<uid>/agent.sock``,
    the launcher's DELFIN_RUNTIME_DIR); None for the owner (the bridge's
    default) and for an unknown name."""
    if name == OWNER:
        return None
    base = home(name, account)
    return None if base is None else base / ".run" / f"delfin-{os.getuid()}" / "agent.sock"


def runtime_dir(name: str, account: Path | None = None) -> Path | None:
    sock = runtime_socket(name, account)
    return None if sock is None else sock.parent


def door_root(name: str, account: Path | None = None) -> Path | None:
    """The private runtime directory a person's pairing door must lie in:
    the owner's is this process's own (the bridge runs with the owner's
    runtime directory), a person's ``<home>/.run/delfin-<uid>``."""
    if name == OWNER:
        from delfin.agent.service import rundir
        return rundir.runtime_dir()
    return runtime_dir(name, account)


def activation_dir(name: str, account: Path | None = None) -> Path | None:
    base = home(name, account)
    return None if base is None else base / ".delfin" / "messenger_sessions"


def pairing_dir(name: str, account: Path | None = None) -> Path | None:
    base = home(name, account)
    return None if base is None else base / ".delfin" / "messenger_pairing"


#: A host with the starter unit (pi_messenger.sh: delfin-person-start.path)
#: names its request directory here; the bridge, which runs without network
#: in its own sandbox, then only leaves a request file, and the starter (the
#: service user, no sandbox of the bridge's, NoNewPrivileges) runs the start.
START_REQUESTS_ENV = "DELFIN_MESSENGER_PERSON_STARTS"

#: How long the bridge waits for a started service to listen.
START_WAIT_S = 30.0

#: The start command, as the agent-service unit runs it: the account's
#: agent.env (the sandbox settings every terminal sources), then the given
#: command. ``$1`` is agent.env, the rest the command.
_WITH_AGENT_ENV = 'f="$1"; shift; if [ -f "$f" ]; then . "$f"; fi; exec "$@"'


def service_env(name: str, account: Path | None = None) -> dict | None:
    """The environment a person's own windows run with (the Windows
    launcher's people mode, remote_bootstrap.sh): their HOME, DELFIN_PERSON,
    their runtime directory, the account's settings as defaults. None for
    anything but a person of people mode."""
    account = account or account_home()
    base = home(name, account) if name != OWNER else None
    if base is None:
        return None
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("DELFIN_MESSENGER_") and k not in {"DELFIN_PERSON", "DELFIN_RUNTIME_DIR"}}
    env.update(HOME=str(base), DELFIN_PERSON=name, DELFIN_RUNTIME_DIR=str(base / ".run"))
    # As the launcher exports it: the account's file, whether or not it exists yet.
    env["DELFIN_SETTINGS_DEFAULTS"] = str(account / ".delfin_settings.json")
    return env


def start_service(name: str, *, run=None, account: Path | None = None) -> bool:
    """Start person ``name``'s agent service the way their own windows do:
    the launcher's start command (``python -m delfin.agent.service.launcher``,
    which starts the supervisor detached and waits for its hello), in the
    person's environment, same Unix user, after the account's agent.env.
    True when the launcher said the service runs. Refused (False, nothing
    run) for anything but a real, owned ``~/people/<name>``."""
    import subprocess
    from delfin.agent.service import launcher
    account = account or account_home()
    env = service_env(name, account)
    if env is None:
        return False
    runtime = Path(env["DELFIN_RUNTIME_DIR"])
    runtime.mkdir(mode=0o700, exist_ok=True)
    if os.path.islink(runtime) or not runtime.is_dir():
        return False
    cmd = ["/bin/bash", "-c", _WITH_AGENT_ENV, "delfin-person-service",
           str(account / "delfin" / "agent.env"), *launcher.service_command()]
    run = run or subprocess.run
    try:
        done = run(cmd, env=env, cwd=str(runtime), stdin=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   start_new_session=True, timeout=90, check=False)
    except (OSError, subprocess.SubprocessError):
        return False
    return getattr(done, "returncode", 1) == 0


def transient_argv(name: str, *, account: Path | None = None, python: str | None = None) -> list | None:
    """``systemd-run`` for person ``name``'s agent service as a transient
    unit of its own (``delfin-person-<name>``: its own cgroup, stopped like
    the account's unit, KillMode=mixed): the supervisor, as
    delfin-agent-service runs it, as the account's user, in the person's
    environment after the account's agent.env. None for anything but a
    real ~/people/<name> owned by the account's user."""
    import shutil
    import sys
    account = account or account_home()
    try:
        owner = os.stat(account)
    except OSError:
        return None
    if name == OWNER or not valid(name) or home(name, account, uid=owner.st_uid) is None:
        return None
    base = account / "people" / name
    env = {"HOME": str(base), "DELFIN_PERSON": name, "DELFIN_RUNTIME_DIR": str(base / ".run"),
           "DELFIN_SETTINGS_DEFAULTS": str(account / ".delfin_settings.json"),
           "PATH": os.environ.get("PATH", "/usr/local/bin:/usr/bin:/bin")}
    argv = [shutil.which("systemd-run") or "/usr/bin/systemd-run", "--quiet", "--collect",
            f"--unit=delfin-person-{name}", f"--uid={owner.st_uid}", f"--gid={owner.st_gid}",
            f"--working-directory={base}", "--property=KillMode=mixed",
            "--property=TimeoutStopSec=60"]
    argv += [f"--setenv={k}={v}" for k, v in env.items()]
    argv += ["/bin/bash", "-c", _WITH_AGENT_ENV, "delfin-person-service",
             str(account / "delfin" / "agent.env"), python or sys.executable, "-P",
             "-m", "delfin.agent.service.supervisor"]
    return argv


def start_unit(name: str, *, account: Path | None = None, python: str | None = None,
               run=None) -> bool:
    """The starter unit's start: :func:`transient_argv`, run once."""
    import subprocess
    argv = transient_argv(name, account=account, python=python)
    if argv is None:
        return False
    try:
        done = (run or subprocess.run)(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL, timeout=60, check=False)
    except (OSError, subprocess.SubprocessError):
        return False
    return getattr(done, "returncode", 1) == 0


def request_start(name: str, *, run=None, account: Path | None = None) -> bool:
    """Start ``name``'s service: directly, or, where the starter unit is
    installed (START_REQUESTS_ENV), by a request file it acts on."""
    if not valid(name) or name == OWNER or home(name, account) is None:
        return False
    where = os.environ.get(START_REQUESTS_ENV, "").strip()
    if not where:
        return start_service(name, run=run, account=account)
    path = Path(where) / f"start-{name}"
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0), 0o600)
        os.close(fd)
    except OSError:
        return False
    return True


def start_requested(directory: Path, *, start=None, uid: int | None = None) -> list[str]:
    """The starter unit's work: each ``start-<name>`` file of ``uid`` (the
    owner of ``directory``, the bridge's user) in ``directory`` is removed,
    and ``name``'s service started when it is a person. Returns the names
    started."""
    start = start or start_service
    if uid is None:
        try:
            uid = os.stat(directory).st_uid
        except OSError:
            return []
    done = []
    try:
        entries = sorted(os.listdir(directory))
    except OSError:
        return done
    for entry in entries:
        path = Path(directory) / entry
        try:
            info = os.lstat(path)
            os.unlink(path)
        except OSError:
            continue
        name = entry[len("start-"):] if entry.startswith("start-") else ""
        if info.st_uid != uid or not valid(name) or name == OWNER:
            continue
        if start(name):
            done.append(name)
    return done


def main(argv: list[str] | None = None) -> int:
    """``python -m delfin.messenger.people start-requested DIR PYTHON`` (the
    starter unit, as root from the root-owned gateway venv): each request
    becomes a transient unit of the request directory owner's user, running
    PYTHON (that user's DELFIN interpreter)."""
    import pwd
    import sys
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3 or args[0] != "start-requested":
        print("usage: python -m delfin.messenger.people start-requested DIR PYTHON", file=sys.stderr)
        return 2
    directory, python = Path(args[1]), args[2]
    try:
        uid = os.stat(directory).st_uid
        account = Path(pwd.getpwuid(uid).pw_dir)
    except (OSError, KeyError):
        return 1
    if uid == 0:
        return 1                                   # never a person service as root
    started = start_requested(directory, uid=uid,
                              start=lambda name: start_unit(name, account=account, python=python))
    for name in started:
        print(f"started the agent service of {name}")
    return 0


def workspaces_root(name: str, owner_root: str, account: Path | None = None) -> Path | None:
    """Where ``neu NAME`` puts a person's session folders: the owner's
    configured root (DELFIN_MESSENGER_WORKSPACES), or ``<home>/DELFIN-Handy``
    of the person when the owner's root is configured at all."""
    if not owner_root:
        return None
    if name == OWNER:
        return Path(owner_root)
    base = home(name, account)
    return None if base is None else base / "DELFIN-Handy"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
