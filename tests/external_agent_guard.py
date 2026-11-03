"""The test suite reaches no external agent CLI and no remote model endpoint.

Input: the host as the suite finds it -- an external agent CLI
(``claude``, ``codex``) possibly on PATH, provider keys and the calling
agent's own session variables possibly in the environment.

Output, installed once per session by ``tests/conftest.py``:

* Discovery. ``shutil.which`` answers ``None`` for an agent CLI name
  whenever the answer would be a real host binary, so every product
  site that asks "is the CLI installed?" (``CLIClient``,
  ``CodexCLIClient``, the dashboard agent tab's provider list, the
  settings tab, ``delfin.dashboard``) sees none. Product code resolves
  the binaries at seven separate ``shutil.which`` calls, not through one
  helper, so the discovery function itself is the single place.
* Spawning. ``subprocess.Popen._execute_child`` (every ``subprocess``
  and ``asyncio`` child), ``os.execv``/``os.execve`` and
  ``os.posix_spawn``/``os.posix_spawnp`` refuse a program, or an argument,
  whose real path is a host agent CLI or lies in its install directory
  (``~/.local/share/<cli>/versions``, an npm package directory). A shell
  that resolves the bare name finds a tripwire script first on PATH
  instead; the script logs its argv and exits 127.
* Network. ``socket.socket.connect``/``connect_ex`` refuse an AF_INET or
  AF_INET6 destination that is not loopback. AF_UNIX is not touched.
* Credentials. Provider keys, tokens and the calling agent's own
  session variables (``CLAUDE*``, ``CODEX*``, ``ANTHROPIC_*``,
  ``OPENAI_*``, ``KIT_*``, ``OLLAMA_*``) are removed from ``os.environ``
  for the session and restored at its end; the per-test credential store
  redirect stays in ``conftest._ci_parity_no_secrets``.

A refusal raises :class:`ExternalAgentTripwire` (a ``PermissionError``,
so the product's own "not installed / unreachable" handling runs) AND is
recorded; the per-test fixture fails the test at teardown with the
command or address, so a product ``except Exception`` cannot hide it.

Opt-out, still bounded: ``@pytest.mark.fake_agent_cli`` (or the
``fake_agent_cli`` fixture) puts stub ``claude``/``codex`` executables
first on PATH, so discovery by name finds them and a test may start them;
each stub appends its argv to ``<stub dir>/calls.log`` and exits 0. Real
binaries and non-loopback connections stay refused under the marker.

Rejected: removing the CLI's directories from PATH (``~/.local/bin``,
``/usr/local/bin`` and the node bin directory also hold tools tests use);
a per-test ``monkeypatch`` of each ``shutil.which`` call site (seven
sites today, and the eighth would be missed); a silent skip (the
2026-11-03 incident was a test that passed while 13 prompts went out).
Measured cost: one ``realpath`` per spawned argument that contains a
``/`` and one ``ipaddress`` parse per connect.
"""

from __future__ import annotations

import ipaddress
import os
import shutil
import socket
import stat
import subprocess
import tempfile
import threading

#: Program names of external agent CLIs.
AGENT_CLI_NAMES = ("claude", "codex")

#: Environment prefixes removed for the session: provider credentials,
#: provider endpoints, and the variables an attended agent session hands
#: its children (one of them names the CLI's own executable).
STRIPPED_PREFIXES = ("CLAUDE", "CODEX", "ANTHROPIC_", "OPENAI_", "KIT_",
                     "OLLAMA_")

MARKER = "fake_agent_cli"


class ExternalAgentTripwire(PermissionError):
    """A test tried to reach an external agent CLI or a remote endpoint."""


_STATE: dict = {
    "installed": False,
    "real": frozenset(),          # realpaths of host agent CLI binaries
    "prefixes": (),               # their install directories
    "shim_dir": "",
    "shim_log": "",
    "saved_env": {},
    "originals": {},
}
_TRIPS: list[str] = []
_LOCK = threading.Lock()


def _where() -> str:
    """The innermost frames of the product (``delfin/``) and the test that
    led to the refusal: what a reader needs to replace it with a fake."""
    import traceback
    frames = [f for f in traceback.extract_stack()[:-3]
              if (os.sep + "delfin" + os.sep in f.filename
                  or os.path.basename(f.filename).startswith("test_"))
              and not f.filename.endswith("external_agent_guard.py")]
    return " <- ".join(f"{os.path.basename(f.filename)}:{f.lineno} {f.name}"
                       for f in reversed(frames[-4:]))


def _record(what: str) -> None:
    where = _where()
    with _LOCK:
        _TRIPS.append(f"{what}  [at {where}]" if where else what)


def take_trips() -> list[str]:
    """Every refusal since the last call, in-process and from the shim."""
    with _LOCK:
        out = list(_TRIPS)
        _TRIPS.clear()
    log = _STATE.get("shim_log")
    if log and os.path.exists(log):
        try:
            with open(log, encoding="utf-8", errors="replace") as fh:
                lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
            open(log, "w").close()
        except OSError:
            lines = []
        out.extend(f"spawned through the shell: {ln}" for ln in lines)
    return out


def _host_binaries(path_value: str) -> tuple[frozenset, tuple]:
    """Real paths of every agent CLI on ``path_value``, and the install
    directories they resolve into (not a directory on PATH itself)."""
    dirs = [d for d in path_value.split(os.pathsep) if d]
    on_path = {os.path.realpath(d) for d in dirs}
    real: set[str] = set()
    prefixes: set[str] = set()
    candidates = [os.path.join(d, n) for d in dirs for n in AGENT_CLI_NAMES]
    for var in ("CLAUDE_CODE_EXECPATH",):
        if os.environ.get(var):
            candidates.append(os.environ[var])
    for cand in candidates:
        if not (os.path.isfile(cand) and os.access(cand, os.X_OK)):
            continue
        rp = os.path.realpath(cand)
        real.add(rp)
        real.add(os.path.abspath(cand))
        parent = os.path.dirname(rp)
        if parent not in on_path:
            prefixes.add(parent + os.sep)
    return frozenset(real), tuple(sorted(prefixes))


def is_host_agent_cli(path: str) -> bool:
    """Whether ``path`` is, or lies inside the install of, a host agent CLI."""
    try:
        rp = os.path.realpath(path)
    except (OSError, ValueError):
        return False
    if rp in _STATE["real"]:
        return True
    return any(rp.startswith(p) for p in _STATE["prefixes"])


def _write_stub(directory: str, name: str, body: str) -> None:
    target = os.path.join(directory, name)
    with open(target, "w") as fh:
        fh.write(body)
    os.chmod(target, stat.S_IRWXU)


def make_stub_dir(directory: str) -> str:
    """Stub agent CLIs for ``@pytest.mark.fake_agent_cli``: log argv, exit 0."""
    log = os.path.join(directory, "calls.log")
    for name in AGENT_CLI_NAMES:
        _write_stub(directory, name,
                    "#!/bin/sh\n"
                    f"printf '%s\\n' \"$0 $*\" >> '{log}'\n"
                    "exit 0\n")
    return directory


# --------------------------------------------------------------------------
# the wrappers
# --------------------------------------------------------------------------

def _which(cmd, mode=os.F_OK | os.X_OK, path=None):
    original = _STATE["originals"]["which"]
    name = os.path.basename(str(cmd)) if cmd else ""
    if name not in AGENT_CLI_NAMES:
        return original(cmd, mode, path)
    if os.sep in str(cmd):
        if is_host_agent_cli(str(cmd)):
            _record(f"resolved a host agent CLI by path: {cmd}")
            return None
        return original(cmd, mode, path)
    found = original(cmd, mode, path)
    if found is None:
        return None
    if (os.path.dirname(os.path.abspath(found)) == _STATE["shim_dir"]
            or is_host_agent_cli(found)):
        return None
    return found


def _refuse_program(argv, executable=None, env=None) -> None:
    if isinstance(argv, (str, bytes, os.PathLike)):
        argv = [argv]
    items = [os.fsdecode(a) for a in (argv or []) if isinstance(a, (str, bytes, os.PathLike))]
    if executable is not None:
        items.insert(0, os.fsdecode(executable))
    if not items:
        return
    prog = items[0]
    if os.sep not in prog:
        search = None
        if env is not None:
            search = env.get("PATH") if hasattr(env, "get") else None
        resolved = _STATE["originals"]["which"](prog, path=search)
        if resolved and os.path.dirname(os.path.abspath(resolved)) == _STATE["shim_dir"]:
            _record(f"spawned an agent CLI by name: {' '.join(items)[:300]}")
            raise ExternalAgentTripwire(
                f"test tripwire: spawning the external agent CLI {prog!r}")
        prog = resolved or prog
    for candidate in [prog] + [a for a in items[1:8] if os.sep in a]:
        if is_host_agent_cli(candidate):
            _record(f"spawned a host agent CLI: {' '.join(items)[:300]}")
            raise ExternalAgentTripwire(
                f"test tripwire: spawning the external agent CLI {candidate!r}")


def _execute_child(self, args, executable, *rest, **kw):
    # positional order (3.11-3.13): ..., cwd, env, ...; env is rest[4]
    env = rest[4] if len(rest) > 4 else kw.get("env")
    shell = rest[7] if len(rest) > 7 else kw.get("shell", False)
    if not shell:
        _refuse_program(args, executable, env)
    return _STATE["originals"]["_execute_child"](self, args, executable, *rest, **kw)


def _wrap_exec(key):
    original = _STATE["originals"][key]

    def wrapper(path, args, *rest, **kw):
        _refuse_program([path] + list(args or [])[1:])
        return original(path, args, *rest, **kw)
    wrapper.__name__ = getattr(original, "__name__", key)
    return wrapper


def _is_loopback(host) -> bool:
    if isinstance(host, bytes):
        host = host.decode("ascii", "replace")
    host = str(host).strip("[]").split("%", 1)[0]
    if host in ("", "localhost") or host.endswith(".localhost"):
        return True
    try:
        ip = ipaddress.ip_address(host)
    except ValueError:
        try:
            infos = _STATE["originals"]["getaddrinfo"](host, None)
        except OSError:
            return True       # unresolvable: the connect cannot leave either
        return bool(infos) and all(
            _is_loopback(info[4][0]) for info in infos)
    if ip.is_loopback or ip.is_unspecified:
        return True
    mapped = getattr(ip, "ipv4_mapped", None)
    return bool(mapped and mapped.is_loopback)


def _check_address(sock, address) -> None:
    if sock.family not in (socket.AF_INET, socket.AF_INET6):
        return
    if not isinstance(address, tuple) or not address:
        return
    if _is_loopback(address[0]):
        return
    _record(f"network connection to {address[0]}:{address[1] if len(address) > 1 else '?'}")
    raise ExternalAgentTripwire(
        f"test tripwire: network connection to non-loopback {address!r}")


def _connect(self, address):
    _check_address(self, address)
    return _STATE["originals"]["connect"](self, address)


def _connect_ex(self, address):
    _check_address(self, address)
    return _STATE["originals"]["connect_ex"](self, address)


# --------------------------------------------------------------------------
# install / uninstall
# --------------------------------------------------------------------------

def install(is_credential) -> None:
    """``is_credential``: the suite's rule for a secret variable name."""
    if _STATE["installed"]:
        return
    path_value = os.environ.get("PATH", "")
    _STATE["real"], _STATE["prefixes"] = _host_binaries(path_value)

    shim = tempfile.mkdtemp(prefix="delfin-agent-cli-tripwire-")
    log = os.path.join(shim, "spawned.log")
    open(log, "w").close()
    for name in AGENT_CLI_NAMES:
        _write_stub(shim, name,
                    "#!/bin/sh\n"
                    "echo \"test tripwire: the external agent CLI $0 was "
                    "started by a test; refused\" >&2\n"
                    f"printf '%s\\n' \"$0 $*\" >> '{log}'\n"
                    "exit 127\n")
    _STATE["shim_dir"] = os.path.realpath(shim)
    _STATE["shim_log"] = log

    saved = {}
    for name in list(os.environ):
        upper = name.upper()
        if upper.startswith(STRIPPED_PREFIXES) or is_credential(upper):
            saved[name] = os.environ.pop(name)
    saved["PATH"] = path_value
    _STATE["saved_env"] = saved
    os.environ["PATH"] = _STATE["shim_dir"] + os.pathsep + path_value

    o = _STATE["originals"]
    o["which"] = shutil.which
    o["_execute_child"] = subprocess.Popen._execute_child
    o["connect"] = socket.socket.connect
    o["connect_ex"] = socket.socket.connect_ex
    o["getaddrinfo"] = socket.getaddrinfo
    shutil.which = _which
    subprocess.Popen._execute_child = _execute_child
    socket.socket.connect = _connect
    socket.socket.connect_ex = _connect_ex
    for key in ("execv", "execve", "posix_spawn", "posix_spawnp"):
        if hasattr(os, key):
            o[key] = getattr(os, key)
            setattr(os, key, _wrap_exec(key))
    _STATE["installed"] = True


def uninstall() -> None:
    if not _STATE["installed"]:
        return
    o = _STATE["originals"]
    shutil.which = o.pop("which")
    subprocess.Popen._execute_child = o.pop("_execute_child")
    socket.socket.connect = o.pop("connect")
    socket.socket.connect_ex = o.pop("connect_ex")
    o.pop("getaddrinfo", None)
    for key in ("execv", "execve", "posix_spawn", "posix_spawnp"):
        if key in o:
            setattr(os, key, o.pop(key))
    saved = _STATE["saved_env"]
    os.environ["PATH"] = saved.pop("PATH", os.environ.get("PATH", ""))
    os.environ.update(saved)
    _STATE["saved_env"] = {}
    shutil.rmtree(_STATE["shim_dir"], ignore_errors=True)
    _STATE["installed"] = False

