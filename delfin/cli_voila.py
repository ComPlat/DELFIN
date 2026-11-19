"""``delfin-voila`` – launch the DELFIN Dashboard via Voila.

Usage::

    delfin-voila               # start on default port 8866
    delfin-voila --port 9000   # custom port
    delfin-voila --open-browser
"""

import argparse
import atexit
import importlib.util
import importlib.resources
import ipaddress
import os
import re
import socket
import shutil
import stat
import shlex
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path


def _get_voila_static_root() -> str | None:
    """Return the installed Voilà static root when available."""
    spec = importlib.util.find_spec("voila")
    if spec is None or not spec.origin:
        return None
    static_root = Path(spec.origin).resolve().parent / "static"
    if static_root.is_dir():
        return str(static_root)
    return None


def _is_safe_notebook_source(path: Path) -> bool:
    """Reject a discovered notebook an attacker could have planted.

    ``_find_notebook`` walks cwd, its parents and $HOME. On a shared host a
    world-writable parent (e.g. ``/tmp``) could hold a hostile
    ``delfin_dashboard.ipynb`` that Voilà would then TRUST and EXECUTE as the
    launching user — arbitrary code execution via a search-path hijack. Only
    accept a search-path notebook that is owned by the current user and is not
    world-writable; anything else falls through to the packaged copy that ships
    with the install (installed by pip, implicitly trusted).
    """
    try:
        st = path.stat()
    except OSError:
        return False
    if hasattr(os, "geteuid"):
        if st.st_uid != os.geteuid():
            return False
        if st.st_mode & 0o002:  # world-writable file
            return False
    return True


def _find_notebook() -> str:
    """Locate the dashboard notebook, preferring a trusted checkout copy."""
    search_roots = [Path.cwd(), *Path.cwd().parents, Path.home()]
    rel_candidates = [
        Path("delfin_dashboard.ipynb"),
        Path("software/delfin/delfin_dashboard.ipynb"),
    ]
    for root in search_roots:
        for rel in rel_candidates:
            candidate = (root / rel).resolve()
            if candidate.is_file() and _is_safe_notebook_source(candidate):
                return str(candidate)

    ref = importlib.resources.files("delfin.notebooks").joinpath(
        "delfin_dashboard.ipynb"
    )
    return str(ref)


def _latest_vscode_ipc_socket() -> str | None:
    """Return the newest live VS Code CLI IPC socket if available."""
    runtime_dir = Path(f"/run/user/{os.getuid()}")
    try:
        candidates = sorted(
            [
                path for path in runtime_dir.glob("vscode-ipc-*.sock")
                if path.is_socket()
            ],
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
    except Exception:
        return None

    for path in candidates:
        try:
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
                client.settimeout(0.15)
                client.connect(str(path))
            return str(path)
        except Exception:
            continue
    return None


def _terminal_gone(signum, frame):
    """SIGTERM or SIGHUP reach the launcher: stop exactly as Ctrl+C does."""
    raise KeyboardInterrupt


def _install_stop_signals() -> dict:
    """A closed terminal (SIGHUP) or a plain kill (SIGTERM) used to end the
    launcher without its cleanup, leaving the server and everything under
    it running: one server stayed up five days. Both now take the Ctrl+C
    path. Returns the handlers they replaced."""
    import signal
    previous: dict = {}
    for name in ("SIGTERM", "SIGHUP"):
        sig = getattr(signal, name, None)
        if sig is None:
            continue
        try:
            previous[sig] = signal.signal(sig, _terminal_gone)
        except (ValueError, OSError):
            pass
    return previous


def _restore_signals(previous: dict) -> None:
    import signal
    for sig, handler in (previous or {}).items():
        try:
            signal.signal(sig, handler)
        except (ValueError, OSError, TypeError):
            pass


def _end_agent_processes() -> None:
    """End every detached process started under this launcher: background
    shells, MCP servers, daemons (the lifeline ledger)."""
    try:
        from delfin.agent import lifeline
        ended = lifeline.end_children()
    except Exception:
        return
    if ended:
        print(f"Stopped {len(ended)} background process(es) started from "
              "this dashboard.")


def _is_vscode_session(env: dict[str, str] | None = None) -> bool:
    """Return True when the launcher runs inside a VS Code terminal session."""
    env = env or os.environ
    term_program = str(env.get("TERM_PROGRAM") or "")
    browser = str(env.get("BROWSER") or "")
    ipc_hook = str(env.get("VSCODE_IPC_HOOK_CLI") or "").strip()
    return term_program == "vscode" or "browser.sh" in browser or bool(ipc_hook)


def _launch_is_interactive() -> bool:
    """True when a person typed the launch: stdin and stdout are terminals."""
    try:
        return bool(sys.stdin and sys.stdin.isatty() and sys.stdout and sys.stdout.isatty())
    except Exception:
        return False


def _decide_open_browser(args, env: dict[str, str] | None = None,
                         interactive: bool | None = None) -> bool:
    """Whether this launch opens the dashboard in a browser.

    The explicit switches win. Otherwise a tab opens only for a launch a
    person typed into a VS Code terminal. A launch from a script, a test,
    an agent or a service inherits the same VS Code variables, and on
    2026-09-16 every such launch -- browser drives, dashboard tests --
    opened a tab on the developer's desktop. Those are never interactive
    (no terminal on stdin/stdout), and ``DELFIN_NO_BROWSER`` switches the
    automatic open off for a whole environment.
    """
    env = env if env is not None else os.environ
    if getattr(args, "open_browser", None) is True:
        return True
    if getattr(args, "no_browser", None) is True:
        return False
    if str(env.get("DELFIN_NO_BROWSER") or "").strip().lower() in ("1", "true", "yes", "on"):
        return False
    if interactive is None:
        interactive = _launch_is_interactive()
    if not interactive:
        return False
    return _is_vscode_session(env)


def _prepare_voila_env(open_browser: bool) -> dict[str, str]:
    env = os.environ.copy()
    if not open_browser:
        return env

    if not _is_vscode_session(env):
        return env

    current_hook = str(env.get("VSCODE_IPC_HOOK_CLI") or "").strip()
    if current_hook and Path(current_hook).is_socket():
        return env

    socket_path = _latest_vscode_ipc_socket()
    if socket_path:
        env["VSCODE_IPC_HOOK_CLI"] = socket_path
    return env


def _get_local_ips() -> list[str]:
    """Return non-loopback IPv4 addresses of this machine."""
    ips: list[str] = []
    try:
        for info in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
            addr = info[4][0]
            if not addr.startswith("127."):
                ips.append(addr)
    except Exception:
        pass
    if not ips:
        # Fallback: connect to a public IP to discover the default route address.
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
                s.connect(("8.8.8.8", 80))
                ips.append(s.getsockname()[0])
        except Exception:
            pass
    return list(dict.fromkeys(ips))  # dedupe, preserve order


def _detect_login_node() -> str | None:
    """Try to identify the login-node DNS name from the SSH connection."""
    # SSH_CONNECTION = "<client_ip> <client_port> <server_ip> <server_port>"
    ssh_conn = os.environ.get("SSH_CONNECTION", "")
    parts = ssh_conn.split()
    if len(parts) >= 3:
        server_ip = parts[2]
        try:
            fqdn = socket.getfqdn(server_ip)
            if fqdn and fqdn != server_ip:
                # Return the short hostname (e.g. "uc3" from "uc3.scc.kit.edu")
                # but keep the domain for external access.
                return fqdn
        except Exception:
            pass
        # Fallback: try reverse DNS
        try:
            hostname = socket.gethostbyaddr(server_ip)[0]
            if hostname:
                return hostname
        except Exception:
            pass
    return None


def _voila_is_available() -> bool:
    """Return True when the current Python can import voila."""
    return importlib.util.find_spec("voila") is not None


def _wait_for_port(host: str, port: int, timeout: float = 10.0) -> bool:
    """Poll until a TCP port becomes reachable or timeout expires."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(0.2)
            if sock.connect_ex((host, port)) == 0:
                return True
        time.sleep(0.1)
    return False


def _open_browser_url(url: str, env: dict[str, str]) -> None:
    """Best-effort browser open for VS Code remote sessions."""
    browser = str(env.get("BROWSER") or "").strip()
    if browser:
        try:
            subprocess.Popen(
                [*shlex.split(browser), url],
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return
        except Exception:
            pass

    try:
        subprocess.Popen(
            [sys.executable, "-m", "webbrowser", "-t", url],
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        pass


def _select_port(requested_port: int) -> int:
    """Return the first free TCP port starting at requested_port."""
    port = requested_port
    max_port = port + 100
    while port <= max_port:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            if sock.connect_ex(("127.0.0.1", port)) != 0:
                return port
        if port == requested_port:
            print(f"Port {port} is in use, searching for a free port...")
        port += 1

    raise RuntimeError(
        f"Error: no free port found in range {requested_port}–{max_port}."
    )


def _stage_notebook_under_root(notebook: str, root_dir: str) -> str:
    """Copy packaged notebooks into root_dir so Voila can serve them safely."""
    notebook_path = Path(notebook).resolve()
    root_path = Path(root_dir).resolve()

    try:
        notebook_path.relative_to(root_path)
        return str(notebook_path)
    except ValueError:
        pass

    staged_dir = root_path / "delfin_voila_runtime"
    staged_path = staged_dir / notebook_path.name
    staged_dir.mkdir(parents=True, exist_ok=True)

    if (
        not staged_path.exists()
        or notebook_path.stat().st_mtime > staged_path.stat().st_mtime
    ):
        shutil.copy2(notebook_path, staged_path)

    return str(staged_path)


def _trust_notebook(notebook: str) -> None:
    """Mark the notebook as trusted so Voilà skips the security warning."""
    try:
        subprocess.run(
            [sys.executable, "-m", "jupyter", "trust", notebook],
            capture_output=True,
            check=False,
        )
    except Exception:
        pass


def _refuse_insecure_no_token(ip: str) -> str:
    """Return a human-readable refusal for disabling auth in Voilà mode."""
    return (
        "Error: --no-token is no longer allowed for delfin-voila.\n"
        "The dashboard can start kernels, run jobs, and drive the DELFIN "
        "agent, so token authentication is mandatory even on 127.0.0.1 "
        "on shared HPC/login nodes.\n"
        "Use the generated token file or put the service behind an "
        "authenticated SSH tunnel/reverse proxy."
    )


#: Frontend extensions the browser is not asked to load.  See the comment at
#: the ``--VoilaConfiguration.extension_denylist`` argument for the measurement.
DEFAULT_EXTENSION_DENYLIST = ('nglview-js-widgets', 'jupyterlab-plotly')


def _extension_denylist() -> list:
    """Which frontend extensions to leave out, honouring an override.

    ``DELFIN_VOILA_EXTENSION_DENYLIST`` is a comma-separated list; empty means
    leave everything in.
    """
    raw = os.environ.get('DELFIN_VOILA_EXTENSION_DENYLIST')
    if raw is None:
        return list(DEFAULT_EXTENSION_DENYLIST)
    return [name.strip() for name in raw.split(',') if name.strip()]


def _validate_access_token(token: str) -> str:
    """Return a validated dashboard token or raise ValueError."""
    value = str(token or "").strip()
    if not value:
        raise ValueError("dashboard token must not be empty")
    if any(ch.isspace() or ord(ch) < 32 for ch in value):
        raise ValueError("dashboard token must not contain whitespace/control characters")
    # token_urlsafe(32) yields 43 chars. Reject weak hand-written tokens so a
    # shared login node cannot be protected by something guessable like "test".
    if len(value) < 32:
        raise ValueError("dashboard token must be at least 32 characters")
    return value


def _is_loopback_bind(ip: str) -> bool:
    """Return True for loopback-only bind addresses."""
    value = str(ip or "").strip().lower()
    if value in {"localhost", "127.0.0.1", "::1"}:
        return True
    try:
        return ipaddress.ip_address(value).is_loopback
    except ValueError:
        return False


def _secure_cache_dir(name: str) -> Path:
    """Create a private per-user DELFIN cache directory."""
    import tempfile

    cache_home = os.environ.get("XDG_CACHE_HOME", "").strip()
    bases = [
        Path(cache_home).expanduser() / "delfin" if cache_home else Path.home() / ".cache" / "delfin",
        Path(tempfile.gettempdir()) / f"delfin-{os.getuid()}",
    ]
    last_error: OSError | None = None
    for base in bases:
        path = base / name
        try:
            path.mkdir(parents=True, exist_ok=True)
            try:
                os.chmod(base, 0o700)
                os.chmod(path, 0o700)
            except OSError:
                pass
            return path
        except OSError as exc:
            last_error = exc
            continue
    raise RuntimeError(f"could not create private cache directory: {last_error}")


def _create_token_file(token: str) -> Path:
    """Write *token* to a 0600 file, falling back from bad runtime dirs."""
    import tempfile

    candidates: list[Path] = []
    runtime_raw = os.environ.get("XDG_RUNTIME_DIR")
    if runtime_raw:
        candidates.append(Path(runtime_raw))
    candidates.append(_secure_cache_dir("tokens"))

    last_error: OSError | None = None
    for directory in candidates:
        try:
            directory.mkdir(parents=True, exist_ok=True)
            try:
                os.chmod(directory, 0o700)
            except OSError:
                pass
            fd, path_str = tempfile.mkstemp(
                prefix=f"voila-token-{os.getpid()}-", dir=str(directory)
            )
            token_file = Path(path_str)
            try:
                os.write(fd, token.encode("utf-8"))
            finally:
                os.close(fd)
            try:
                os.chmod(token_file, 0o600)
            except OSError:
                pass
            return token_file
        except OSError as exc:
            last_error = exc
            continue
    raise RuntimeError(f"could not create secure token file: {last_error}")


#: Set to 1 to accept a DELFIN_VOILA_ROOT_DIR of /, /home, /tmp or $HOME.
BROAD_ROOT_ENV = "DELFIN_VOILA_ALLOW_BROAD_ROOT"


def _broad_roots() -> set:
    """Directories too wide to serve: every user's files, or all of one's own."""
    return {
        Path("/").resolve(),
        Path("/home").resolve(),
        Path("/tmp").resolve(),
        Path.home().resolve(),
    }


def _checked_root_override(raw: str) -> Path | None:
    """The DELFIN_VOILA_ROOT_DIR override, or None when it is not taken.

    Taken only when it names an existing directory that is not one of
    :func:`_broad_roots` (unless BROAD_ROOT_ENV=1). A refused override is
    reported on stderr and the default root is used instead, so the
    dashboard still starts.
    """
    try:
        path = Path(raw).expanduser().resolve()
    except (OSError, RuntimeError):
        path = None
    reason = ""
    if path is None or not path.is_dir():
        reason = "is not an existing directory"
    elif (path in _broad_roots()
          and os.environ.get(BROAD_ROOT_ENV, "").strip() != "1"):
        reason = (f"would serve far more than the dashboard needs "
                  f"(set {BROAD_ROOT_ENV}=1 to allow it)")
    if reason:
        print(f"  \u26a0  DELFIN_VOILA_ROOT_DIR={raw!r} {reason}; "
              "using the default root instead.", file=sys.stderr)
        return None
    return path


def _default_voila_root() -> Path:
    """Return a narrow default root_dir for files served by Voilà."""
    override = os.environ.get("DELFIN_VOILA_ROOT_DIR", "").strip()
    if override:
        checked = _checked_root_override(override)
        if checked is not None:
            return checked

    cwd = Path.cwd().resolve()
    home = Path.home().resolve()
    if cwd not in _broad_roots() and home in cwd.parents:
        return cwd
    return _secure_cache_dir("voila-root").resolve()


def _kernel_transport_args() -> list:
    """The server's kernels talk over Unix sockets, not loopback TCP.

    Over TCP every kernel channel is a port on 127.0.0.1, which any user
    of the machine can connect to. The key in the connection file signs
    messages, it does not encrypt them, and the IOPub channel publishes to
    any subscriber: measured 2026-10-20, a client without the key read a
    kernel's printed output over TCP (ipykernel 7.3, jupyter_client 8.10).
    With ``ipc`` the sockets are files beside the connection file, in the
    Jupyter runtime directory, which is the user's own (mode 0700), so no
    other user can reach them. ipykernel 7 also stops warning about an
    unencrypted TCP kernel. Windows has no ipc transport and keeps TCP.
    ``DELFIN_KERNEL_TRANSPORT=tcp`` restores TCP (a kernel on another
    machine, a remote debugger).

    A socket path is at most 107 bytes on Linux and 103 on macOS
    (sockaddr_un), and each one is ``<runtime>/kernel-<uuid>-ipc-<port>``;
    under a deep runtime directory no kernel could start. Then the server
    keeps TCP and says so, instead of a dashboard that never loads."""
    if os.name == "nt":
        return []
    transport = os.environ.get("DELFIN_KERNEL_TRANSPORT", "ipc").strip().lower()
    if transport != "ipc":
        return []
    try:
        from jupyter_core.paths import jupyter_runtime_dir
        runtime = os.path.realpath(jupyter_runtime_dir())
    except Exception:
        return []
    try:
        import zmq
        limit = int(getattr(zmq, "IPC_PATH_MAX_LEN", 0) or 0)
    except Exception:
        limit = 0
    limit = limit or (103 if sys.platform == "darwin" else 107)

    def fits(path: str) -> bool:
        longest = (len(os.fsencode(path)) + len("/kernel-") + 36
                   + len("-ipc-") + 5)
        return longest <= limit

    if fits(runtime) and _is_private_dir(runtime):
        return ["--KernelManager.transport=ipc"]
    moved = _private_runtime_dir(fits)
    if moved is None:
        print(f"Note: kernels talk over loopback TCP, which other users of "
              f"this machine can reach: {runtime} is too long for a Unix "
              "socket or readable by others, and neither $XDG_RUNTIME_DIR "
              "nor /tmp gave a private directory. Set JUPYTER_RUNTIME_DIR "
              "to a short private directory (mode 0700) to keep them "
              "private.", file=sys.stderr)
        return []
    # Read by the server, which is started with this process's
    # environment (main() copies it into the server's env).
    os.environ["JUPYTER_RUNTIME_DIR"] = moved
    print(f"Kernels use private sockets under {moved} (the default "
          f"runtime directory is too long for a socket or not private).",
          file=sys.stderr)
    return ["--KernelManager.transport=ipc"]


def _is_private_dir(path: str) -> bool:
    """A directory, not a symlink, owned by this user, no group or other
    permissions: where connection files (the kernels' keys) and sockets
    may live."""
    try:
        st = os.lstat(path)
    except OSError:
        return False
    return (stat.S_ISDIR(st.st_mode) and st.st_uid == os.getuid()
            and not st.st_mode & 0o077)


#: Where the fixed and the one-off runtime directories are made: short on
#: every Unix, unlike $TMPDIR on macOS (/var/folders/...).
#: ``DELFIN_PRIVATE_TMP`` replaces it (the test suite points it at a
#: directory of its own, so no test server writes the user's real
#: /tmp/delfin-k-<uid>); the same checks apply to it.
_TMP_ROOT = "/tmp"
PRIVATE_TMP_ENV = "DELFIN_PRIVATE_TMP"


def _private_runtime_dir(fits) -> str | None:
    """A short private directory for the kernels' connection files and
    sockets, or None.

    In order:
    - ``$XDG_RUNTIME_DIR/delfin-jupyter`` when ``$XDG_RUNTIME_DIR`` is this
      user's private directory (on Linux a tmpfs per login, mode 0700);
    - ``/tmp/delfin-k-<uid>``, made 0700 by this user. It is the same
      directory on every start, which the login store in it needs: a
      restart keeps the logins (a fresh directory per start ended them);
    - a new directory from ``mkdtemp`` under /tmp, made 0700 with a name
      nobody can take first, removed when this process exits; only when
      the fixed name exists and is not this user's private directory.
    Every candidate is checked after it exists: a symlink, another owner
    or open permissions is refused, not repaired.

    The rules are the established ones for a per-user socket directory:
    the XDG Base Directory Specification's runtime directory (owned by the
    user, mode 0700), and in a shared /tmp the pattern tmux and ssh-agent
    use (a fixed per-user name, accepted only when it is the user's own
    0700 directory and not a symlink; CWE-377/CWE-379). /tmp itself must
    carry the sticky bit or belong to the user: without it, any user
    could rename or replace a directory in it, so neither the fixed name
    nor a mkdtemp directory is used there.
    """
    xdg = os.environ.get("XDG_RUNTIME_DIR", "").strip()
    if xdg and _is_private_dir(xdg):
        candidate = os.path.join(xdg, "delfin-jupyter")
        try:
            os.mkdir(candidate, 0o700)
        except FileExistsError:
            pass
        except OSError:
            candidate = ""
        if candidate and _is_private_dir(candidate) and fits(candidate):
            return candidate
    root = os.environ.get(PRIVATE_TMP_ENV, "").strip() or _TMP_ROOT
    base = root if os.path.isdir(root) else tempfile.gettempdir()
    try:
        bst = os.lstat(base)
    except OSError:
        return None
    if not (stat.S_ISDIR(bst.st_mode)
            and (bst.st_mode & stat.S_ISVTX or bst.st_uid == os.getuid())):
        return None
    fixed = os.path.join(base, f"delfin-k-{os.getuid()}")
    try:
        os.mkdir(fixed, 0o700)
    except FileExistsError:
        pass
    except OSError:
        fixed = ""
    if fixed and _is_private_dir(fixed) and fits(fixed):
        return fixed
    try:
        made = tempfile.mkdtemp(prefix="delfin-k-", dir=base)
    except OSError:
        return None
    if not (_is_private_dir(made) and fits(made)):
        shutil.rmtree(made, ignore_errors=True)
        return None
    atexit.register(shutil.rmtree, made, True)
    return made


def _file_allowlist() -> str:
    from delfin.dashboard.web_security import FILE_ALLOWLIST

    return FILE_ALLOWLIST


def _dashboard_default_url(notebook: str, root_dir: str) -> str:
    """URL that renders the dashboard, relative to *root_dir*.

    Used as ``ServerApp.default_url`` so a plain ``localhost:PORT/`` request is
    redirected straight to the rendered dashboard — the user experience is
    unchanged even though we now run under an authenticating jupyter-server.
    The notebook is always staged under ``root_dir`` first, so the relative
    path resolves cleanly.
    """
    from urllib.parse import quote

    rel = Path(notebook).resolve().relative_to(Path(root_dir).resolve()).as_posix()
    return "/voila/render/" + quote(rel)


def _login_logo_data_uri() -> str:
    """Return the packaged DELFIN logo as a base64 ``data:`` URI, or ''.

    Mirrors the dashboard logo loader (``dashboard.__init__._load_logo_data_uri``)
    so the branded login page shows the same mark as the running dashboard.
    """
    import base64

    try:
        logo = importlib.resources.files("delfin").joinpath("logo", "DELFIN_logo.png")
        if logo.is_file():
            data = logo.read_bytes()
            if data:
                return "data:image/png;base64," + base64.b64encode(data).decode("ascii")
    except Exception:
        pass
    return ""


def _render_login_template(template: str, logo_uri: str) -> str:
    """``template`` with the DELFIN logo inlined as a ``data:`` URI, or
    without its logo tile when there is no logo (no broken image)."""
    if logo_uri:
        return template.replace("__DELFIN_LOGO_URI__", logo_uri)
    import re

    return re.sub(
        r"<!--LOGO_TILE_START-->.*?<!--LOGO_TILE_END-->",
        "",
        template,
        flags=re.DOTALL,
    ).replace("__DELFIN_LOGO_URI__", "")


def _prepare_login_template_dir() -> str | None:
    """Materialise the DELFIN-branded sign-in and sign-out pages for
    jupyter_server.

    Reads every shipped template in ``voila_assets/login_templates``
    (``login.html``, ``logout.html`` and the card style they include),
    inlines the DELFIN logo as a ``data:`` URI, and writes the results into
    a private 0700 directory that is handed to jupyter_server via
    ``extra_template_paths`` (searched before the packaged templates, so
    ours win). Returns the directory path, or ``None`` if ``login.html`` is
    unavailable (jupyter_server then falls back to its own pages — auth is
    unaffected either way). ``logout.html`` used to be jupyter_server's own:
    an unstyled "Successfully logged out. Proceed to the login page."
    """
    try:
        folder = importlib.resources.files("delfin.voila_assets").joinpath(
            "login_templates")
        templates = {f.name: f.read_text(encoding="utf-8")
                     for f in folder.iterdir()
                     if f.name.endswith(".html") and f.is_file()}
    except Exception:
        return None
    if "login.html" not in templates:
        return None

    logo_uri = _login_logo_data_uri()
    try:
        target_dir = _secure_cache_dir("voila-login")
        for name, template in templates.items():
            (target_dir / name).write_text(
                _render_login_template(template, logo_uri), encoding="utf-8")
        return str(target_dir)
    except Exception:
        return None


#: The dashboard's own Voila page template (``voila_assets/jupyter_data/
#: nbconvert/templates/delfin``): Voila's ``lab`` template without the two
#: Font Awesome stylesheets it links from a CDN, which the page CSP refuses.
#: Font Awesome itself is bundled in Voila's ``voila-style.js`` and served
#: from the dashboard's own origin.
PAGE_TEMPLATE = "delfin"


def _page_template_args(env: dict) -> list[str]:
    """Select the dashboard's page template; return the server arguments.

    Voila finds templates only in the Jupyter data directories, so the
    shipped data directory is put first on ``JUPYTER_PATH`` in *env* (an
    existing ``JUPYTER_PATH`` is kept behind it). The template lives under
    ``nbconvert/templates``, which Voila searches too, and not under
    ``voila/``: Voila takes its data directory (the themes) from the first
    ``voila`` directory on the path, and that must stay its own.

    Returns ``[]`` (Voila's ``lab`` template, as before) when the shipped
    directory cannot be found as a real directory.
    """
    try:
        root = Path(str(importlib.resources.files("delfin.voila_assets")
                        .joinpath("jupyter_data")))
    except Exception:
        return []
    if not (root / "nbconvert" / "templates" / PAGE_TEMPLATE
            / "index.html.j2").is_file():
        return []
    existing = env.get("JUPYTER_PATH", "")
    parts = [str(root)] + [p for p in existing.split(os.pathsep)
                           if p and p != str(root)]
    env["JUPYTER_PATH"] = os.pathsep.join(parts)
    return [f"--VoilaConfiguration.template={PAGE_TEMPLATE}"]


def _prepare_published_dir() -> str | None:
    """Create the one directory the server publishes, and return it.

    The drawing editor is a browser application: it has to be fetched over
    HTTP, so it has to sit somewhere the server will hand out. Until now the
    only such place was Voilà's root directory — which is whatever directory
    the dashboard was started in — so thirty megabytes of Ketcher were copied
    into it, one copy per launch directory, and a user found them in their
    archive folder next to the calculations.

    A static route of its own fixes that: the editor is loaded from where it
    is kept, in ``~/.delfin/served``, and nothing is written into the launch
    directory at all. Returns the directory, or ``None`` if it cannot be made
    (the editor then falls back to the copy under the root, as before).

    ``/static/`` is served WITHOUT a token — jupyter_server needs it for the
    login page — which is why what is published gets a directory of its own
    rather than a route into ``~/.delfin``, where the credentials and the
    agent's memory live. It holds the editor and nothing else, and the editor
    is a public open-source bundle. ``/voila/files/``, which reaches the
    launch directory, requires the token (``web_security``).
    """
    try:
        from delfin.dashboard import ketcher

        published = ketcher.published_root()
        published.mkdir(parents=True, exist_ok=True)
        return str(published)
    except Exception:
        return None


# Benign lines jupyter_server emits that only confuse users without signalling
# any real problem: an internal "ServerApp.token config is deprecated" warning
# (jupyter_server reads its OWN deprecated trait during startup) and a routine
# login-cookie clear. jupyter_server force-shows some of these past Python's
# warning filters (it resets them at startup), so the only reliable way to keep
# them out of the terminal is to filter the child's stderr.
_BENIGN_STDERR_SUBSTRINGS = (
    "ServerApp.token config is deprecated",
    "Clearing invalid/expired login cookie",
    "cannot be longer than the websocket_ping_interval",
    "Setting websocket_ping_timeout",
    # Routine shutdown chatter on Ctrl+C — we print our own "Stopping…/Stopped."
    "Shutting down",
    "received signal",
)


#: A server warning with no text: jupyter_server logs ``str(e)`` of the
#: ``WebSocketClosedError`` raised for each kernel message that reaches a
#: browser connection which has closed, and that exception has no text.
_EMPTY_WARNING_RE = re.compile(rb"^(\[W [^\]]*\])\s*$")

_CLOSED_CONNECTION_NOTE = (
    b" kernel output reached a browser connection that had closed (page "
    b"reloaded, tab closed or connection lost); repeats directly after "
    b"it are not shown\n")


def _forward_filtered_stderr(stream) -> None:
    """Forward the child's stderr to ours, dropping only benign noise lines.

    Matches an explicit allowlist of known-harmless messages and passes
    EVERYTHING else through unchanged — genuine warnings and errors are never
    hidden (fail-loud preserved). A ``warnings.warn`` header is printed with an
    indented source snippet on the next line; that snippet is dropped too when
    its header was dropped.

    A run of warnings without text (``_EMPTY_WARNING_RE``) is printed once,
    with what it means, instead of one empty header per kernel message.
    """
    drop_indented_continuation = False
    in_empty_run = False
    out = getattr(sys.stderr, "buffer", None)
    try:
        for raw in iter(stream.readline, b""):
            empty = _EMPTY_WARNING_RE.match(raw)
            if empty:
                if in_empty_run:
                    continue
                in_empty_run = True
                raw = empty.group(1) + _CLOSED_CONNECTION_NOTE
            else:
                in_empty_run = False
            text = raw.decode("utf-8", "replace")
            if drop_indented_continuation and text[:1].isspace():
                drop_indented_continuation = False
                continue
            drop_indented_continuation = False
            if any(sub in text for sub in _BENIGN_STDERR_SUBSTRINGS):
                if "Warning:" in text:
                    drop_indented_continuation = True
                continue
            if out is not None:
                out.write(raw)
                out.flush()
            else:
                sys.stderr.write(text)
                sys.stderr.flush()
    except Exception:
        # A daemon log-pump must never crash the launcher; go silent on this
        # stream rather than raise in a background thread.
        pass


# -- a server that outlives its terminal --------------------------------------
#
# A launch typed at a terminal starts this same launcher again in a session
# of its own (``_start_detached``): no terminal, so no SIGHUP when the
# terminal closes or the SSH connection drops; stdin from /dev/null; its
# output appended to an owner-only log. The detached launcher is the one
# that runs the server and stops it -- the lifeline of everything the
# dashboard starts, exactly as the foreground launcher was -- so
# ``delfin-voila stop`` (SIGTERM) takes the same Ctrl+C path. The launch
# itself waits until the server answers, prints where it is and returns.
#
# Rejected: tmux (what --keep did): one more program to have, and keys to
# learn, for what a session of its own does. A systemd user unit: not on
# every host (macOS, a container, a cluster without lingering), and the
# Pi's system unit already runs the launcher in the foreground. A daemon
# library (python-daemon): a dependency for one detached start the
# standard library does in a few lines, the way the agent service starts
# (agent/service/launcher.py).

#: How long a detached start waits for its server to answer.
DETACHED_START_TIMEOUT_S = 120.0

#: How long ``delfin-voila stop`` waits for the launcher to end. The
#: launcher gives its server 15 s, then the background processes 5 s.
STOP_TIMEOUT_S = 60.0

#: How long ``stop`` waits for kept agent sessions to reach the service.
_HANDOVER_WAIT_S = 30.0

#: ``DELFIN_VOILA_FOREGROUND=1`` keeps every launch in the foreground.
FOREGROUND_ENV = "DELFIN_VOILA_FOREGROUND"

#: What ``--keep`` says, once, now that it starts no tmux.
KEEP_NOTE = ("--keep: no tmux any more -- the dashboard keeps running after "
             "this terminal closes by default, and every window it opens "
             "starts with Keep dashboard on.")

#: The detached child imports DELFIN from where this process did, and never
#: from the launch folder it runs in: ``python -c`` puts the working folder
#: first on sys.path, and a json.py there would be imported in place of the
#: standard one.
_CHILD_BOOT = (
    "import os, sys\n"
    "sys.path[:] = [p for p in sys.path if p not in ('', '.')]\n"
    "root = os.environ.pop('DELFIN_VOILA_IMPORT_ROOT', '')\n"
    "if root:\n"
    "    import importlib.util as u\n"
    "    spec = u.find_spec('delfin')\n"
    "    where = os.path.dirname(os.path.dirname(spec.origin)) if spec and spec.origin else ''\n"
    "    if os.path.realpath(where or '/nonexistent') != os.path.realpath(root):\n"
    "        sys.path.insert(0, root)\n"
    "from delfin.cli_voila import main\n"
    "raise SystemExit(main(sys.argv[1:]))\n"
)

_TRUTHY = ("1", "true", "yes", "on")


def _stdin_is_a_terminal() -> bool:
    try:
        return bool(sys.stdin is not None and sys.stdin.isatty())
    except Exception:
        return False


def _single_server_applies(env) -> bool:
    """One server per user and host, except a branch preview (its own
    server beside the real one) and a SLURM job (``--job`` keeps its own
    list, and a job on a node is no login of the user's)."""
    return _announces_where(env) and not env.get("SLURM_JOB_ID")


def _should_detach(args, env=None) -> tuple:
    """(detach, why): whether this launch starts the server detached.

    Only a launch a person typed at a terminal does. A caller that owns
    the dashboard's lifetime -- a script, a service manager (systemd,
    launchd: neither gives a terminal), a SLURM job, tmux or screen, a
    branch preview -- keeps it in the foreground, as before.
    """
    from delfin.dashboard import server_record as _sr

    e = os.environ if env is None else env
    if getattr(args, "foreground", False):
        return False, "--foreground"
    if _sr.is_detached(e):
        return False, "this is the detached launcher already"
    if str(e.get(FOREGROUND_ENV, "")).strip().lower() in _TRUTHY:
        return False, f"{FOREGROUND_ENV} is set"
    if os.name == "nt":
        return False, "Windows: the server runs in this window"
    if not _announces_where(e):
        return False, "a branch preview closes with its branch"
    if e.get("SLURM_JOB_ID"):
        return False, "inside a SLURM job, which owns its lifetime"
    if e.get("TMUX"):
        return False, "inside tmux, which outlives the terminal already"
    if e.get("STY"):
        return False, "inside screen, which outlives the terminal already"
    if not _stdin_is_a_terminal():
        return False, "no terminal: a script or a service manager owns it"
    return True, "typed at a terminal"


#: Options with a value whose value must not reach the child's command line.
_SECRET_OPTIONS = ("--token",)
#: Options the launch itself handles (the browser), not the child.
_PARENT_OPTIONS = ("--open-browser", "--no-browser", "--foreground")


def _detached_child_argv(argv) -> list:
    """The detached launcher's arguments: these, without the token (it goes
    in the environment: /proc/<pid>/cmdline is readable by every user)
    and the browser switches (the launch opens the browser), plus
    ``--foreground --no-browser``."""
    out: list = []
    skip = False
    for arg in list(argv or []):
        if skip:
            skip = False
            continue
        if arg in _SECRET_OPTIONS:
            skip = True
            continue
        if any(arg.startswith(opt + "=") for opt in _SECRET_OPTIONS):
            continue
        if arg in _PARENT_OPTIONS:
            continue
        out.append(arg)
    return out + ["--foreground", "--no-browser"]


def _detached_child_env(token: str) -> dict:
    """The detached launcher's environment: this one, without the lifeline
    of the terminal it was typed in (it is the lifeline of what it starts)
    and without a server token it might have inherited (typed in a
    dashboard's own terminal), with the token it is to use."""
    import delfin
    from delfin.dashboard import server_record as _sr

    env = {k: v for k, v in os.environ.items()
           if k not in ("DELFIN_LIFELINE_PID", "DELFIN_LIFELINE_TICKS",
                        "JUPYTER_TOKEN", "JPY_API_TOKEN", "DELFIN_VOILA_PORT")}
    env["DELFIN_VOILA_TOKEN"] = token
    env[_sr.DETACHED_ENV] = "1"
    env["DELFIN_LAUNCH_CWD"] = os.getcwd()
    env["PYTHONUNBUFFERED"] = "1"
    env["DELFIN_VOILA_IMPORT_ROOT"] = str(Path(delfin.__file__).resolve().parents[1])
    return env


def _detached_child_command(child_argv) -> list:
    return [sys.executable, "-c", _CHILD_BOOT, *list(child_argv)]


def _log_tail(path, lines: int = 15) -> str:
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    return "\n".join(text.splitlines()[-lines:])


def _launch_token(args) -> str:
    """The token for a new detached server: --token, DELFIN_VOILA_TOKEN, or
    a new one. Raises ValueError for a weak one."""
    import secrets

    given = args.token or os.environ.get("DELFIN_VOILA_TOKEN")
    if given:
        return _validate_access_token(given)
    return secrets.token_urlsafe(32)


def _start_detached(args, argv) -> int:
    """Start the launcher again in a session of its own; return once its
    server answers (0) or it gave up (its exit code, and its log)."""
    from delfin.agent import state_paths
    from delfin.dashboard import server_record as _sr

    if args.no_token:
        print(_refuse_insecure_no_token(args.ip or "127.0.0.1"), file=sys.stderr)
        return 2
    try:
        token = _launch_token(args)
    except ValueError as exc:
        print(f"Error: insecure --token: {exc}.", file=sys.stderr)
        return 2
    # A login node keeps its limits: the load watcher samples every DELFIN
    # process of the user there by its command line, a detached one too.
    try:
        from delfin import hpc_guard
        if hpc_guard.on_login_node():
            hpc_guard.ensure_watcher()
    except Exception:
        pass
    log = _sr.log_file()
    state_paths.ensure_dir(log.parent)
    fd = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    try:
        state_paths.secure_file(log)
        os.write(fd, (f"\n--- delfin-voila, detached, "
                      f"{time.strftime('%Y-%m-%d %H:%M:%S')} ---\n").encode())
        print("Starting the DELFIN dashboard in the background ...", flush=True)
        child = subprocess.Popen(
            _detached_child_command(_detached_child_argv(argv)),
            stdin=subprocess.DEVNULL, stdout=fd, stderr=fd,
            start_new_session=True, close_fds=True, cwd=os.getcwd(),
            env=_detached_child_env(token))
    finally:
        os.close(fd)
    record = _sr.wait_for_record(pid=int(getattr(child, "pid", 0) or 0),
                                 timeout=DETACHED_START_TIMEOUT_S,
                                 poll=child.poll)
    if record is None:
        code = child.poll()
        if code == 0:
            # It found a server another launch had started a moment
            # earlier, and left: that one is the answer.
            running = _sr.running()
            if running is not None:
                return _use_running(running, args, argv)
        if code is None:
            print(f"The dashboard server has not answered within "
                  f"{DETACHED_START_TIMEOUT_S:.0f} s and is still starting; "
                  f"`delfin-voila status` shows it once it does. Its log: "
                  f"{log}\n{_log_tail(log)}", file=sys.stderr)
            return 1
        print(f"The dashboard server did not start (exit {code}). Its log: "
              f"{log}\n{_log_tail(log)}", file=sys.stderr)
        return int(code) or 1
    _print_ready(record, fresh=True)
    _open_for(args, record)
    return 0


def _open_for(args, record) -> None:
    """Open the dashboard in a browser where the launch would have."""
    if not _decide_open_browser(args):
        return
    env = _prepare_voila_env(open_browser=True)
    _open_browser_url(_direct_url(record, with_token=True), env)


def _direct_url(record, *, with_token: bool) -> str:
    scheme = str(record.get("scheme") or "http")
    ip = str(record.get("ip") or "127.0.0.1")
    host = "localhost" if ip in ("127.0.0.1", "0.0.0.0", "localhost") else ip
    url = f"{scheme}://{host}:{int(record.get('port') or 0)}{record.get('url') or '/'}"
    if with_token and record.get("token"):
        url += f"?token={record['token']}"
    return url


def _windows(n) -> str:
    if n is None:
        return "windows unknown"
    return "no window" if n == 0 else f"{n} window{'s' * (n != 1)}"


def _print_kept(record, kernel_list) -> None:
    from delfin.dashboard import server_record as _sr

    kept = _sr.kept_dashboards(record, kernel_list)
    if not kept:
        return
    print(f"  Kept dashboards ({len(kept)}):")
    for row in kept:
        state = ("gone (its kernel ended)" if row["state"] == "gone" else
                 f"last activity {_sr.ago(row['last_activity'])}, "
                 f"{_windows(row['windows'])}")
        print(f"    {row['name']}  {state}")
        if row["url"] and row["state"] != "gone":
            print(f"      {row['url']}")


def _print_ready(record, *, fresh: bool = False) -> None:
    """Where the dashboard is, what it holds, how to end it."""
    from delfin.dashboard import server_record as _sr

    if sys.stdout.isatty():
        print(f"  🔗 Open (you are logged in directly): "
              f"{_direct_url(record, with_token=True)}")
    else:
        # Stdout may be a log other users can read: no token in it.
        print(f"  Open: {_direct_url(record, with_token=False)}")
        print(f"  (the login token is in {_sr.record_file()}, readable only "
              "by you)")
    if record.get("detached"):
        if fresh:
            print("  It runs in the background: closing this terminal ends "
                  "nothing.")
        if record.get("log"):
            print(f"  Log:    {record['log']}")
    ssh = os.environ.get("SSH_CONNECTION", "")
    if ssh and _is_loopback_bind(str(record.get("ip") or "127.0.0.1")):
        node = _detect_login_node() or socket.gethostname()
        user = os.environ.get("USER", os.environ.get("LOGNAME", "<user>"))
        port = int(record.get("port") or 0)
        print(f"  From your computer:  ssh -N -L {port}:localhost:{port} "
              f"{user}@{node}")
        print(f"  then open http://localhost:{port} there.")
    _print_kept(record, _sr.kernels(record, timeout=3.0))
    print("  Status: delfin-voila status    End it: delfin-voila stop\n")


#: The options that only shape a NEW server.
_NEW_SERVER_OPTIONS = ("--port", "--ip", "--dark", "--token", "--resume",
                       "--agent-terminal", "--certfile", "--keyfile",
                       "--allow-remote-bind", "--keep")


def _new_server_options(argv, record) -> list:
    """The options given that a running server cannot take any more."""
    items = list(argv or [])
    named = []
    for arg in items:
        name = arg.split("=", 1)[0]
        if name in _NEW_SERVER_OPTIONS and name not in named:
            named.append(name)
    if "--port" in named:
        wanted = None
        for i, arg in enumerate(items):
            try:
                if arg == "--port" and i + 1 < len(items):
                    wanted = int(items[i + 1])
                elif arg.startswith("--port="):
                    wanted = int(arg.split("=", 1)[1])
            except ValueError:
                pass
        if wanted == int(record.get("port") or 0):
            named.remove("--port")
    return named


def _use_running(record, args, argv) -> int:
    """A server already runs on this host: open it, start nothing."""
    port = int(record.get("port") or 0)
    if getattr(args, "strict_port", False) and port != int(args.port):
        print(f"Error: a dashboard server already runs on this host on port "
              f"{port}; --strict-port asked for {args.port}. "
              "`delfin-voila stop` ends it.", file=sys.stderr)
        return 1
    print(f"The DELFIN dashboard already runs on this host (port {port}).")
    ignored = _new_server_options(argv, record)
    if ignored:
        print(f"  {', '.join(ignored)} shape a new server; this one keeps "
              "its own (delfin-voila stop, then start again).")
    _print_ready(record)
    _open_for(args, record)
    return 0


def _await_hand_over(ids, timeout: float) -> tuple:
    """(handed, pending): which of the kept agent sessions *ids* left their
    record (resumed in the agent service) within *timeout*."""
    from delfin.dashboard import kept_agent_sessions as _kas

    wanted = list(ids)
    deadline = time.monotonic() + timeout
    while True:
        left = {str(r.get("session_id")) for r in _kas.records()}
        pending = [i for i in wanted if i in left]
        if not pending or time.monotonic() >= deadline:
            return [i for i in wanted if i not in pending], pending
        time.sleep(0.25)


def _status() -> int:
    from delfin.dashboard import server_record as _sr

    record = _sr.running()
    if record is None:
        print("No dashboard server runs on this host. Start one with: "
              "delfin-voila")
        return 1
    try:
        up = _sr.span(time.time() - float(record.get("started_at") or 0))
    except (TypeError, ValueError):
        up = "?"
    how = "detached" if record.get("detached") else "in a terminal"
    print(f"Dashboard server: running on {record.get('host')}, port "
          f"{record.get('port')}, up {up} (pid {record.get('pid')}, {how})")
    kernel_list = _sr.kernels(record)
    if kernel_list is None:
        print("  The server did not answer; its windows are not listed.")
    else:
        windows = sum(int(k.get("connections") or 0) for k in kernel_list
                      if isinstance(k, dict))
        n = len(kernel_list)
        print(f"  {n} dashboard{'s' * (n != 1)} open, {_windows(windows)}")
    if sys.stdout.isatty():
        print(f"  🔗 Open: {_direct_url(record, with_token=True)}")
    else:
        print(f"  Open: {_direct_url(record, with_token=False)}  (the login "
              f"token is in {_sr.record_file()}, readable only by you)")
    if not _sr.kept_dashboards(record, kernel_list):
        print("  No kept dashboard.")
    _print_kept(record, kernel_list)
    agents = _sr.kept_agent_sessions()
    print(f"  {len(agents)} kept agent session{'s' * (len(agents) != 1)} "
          "(they go on in the agent service when the server stops)")
    if record.get("log"):
        print(f"  Log: {record['log']}")
    print("  End it: delfin-voila stop")
    return 0


def _stop() -> int:
    """End this host's server the way Ctrl+C did, and say what ended."""
    import signal

    from delfin.dashboard import kept_agent_sessions as _kas
    from delfin.dashboard import server_record as _sr

    record = _sr.running()
    if record is None:
        print("No dashboard server runs on this host.")
        return 0
    pid = int(record["pid"])
    kernel_list = _sr.kernels(record, timeout=5.0)
    kept = [r["name"] for r in _sr.kept_dashboards(record, kernel_list)
            if r["state"] != "gone"]
    agents = [str(r["session_id"]) for r in _sr.kept_agent_sessions()
              if not _kas.holder_gone(r)]
    try:
        up = _sr.span(time.time() - float(record.get("started_at") or 0))
    except (TypeError, ValueError):
        up = "?"
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    except PermissionError as exc:
        print(f"Could not stop the dashboard server (pid {pid}): {exc}",
              file=sys.stderr)
        return 1
    deadline = time.monotonic() + STOP_TIMEOUT_S
    while _sr.alive(record) and time.monotonic() < deadline:
        time.sleep(0.2)
    if _sr.alive(record):
        print(f"The dashboard server (pid {pid}) has not ended within "
              f"{STOP_TIMEOUT_S:.0f} s; it may still be ending its kernels. "
              f"Its log: {record.get('log') or '-'}", file=sys.stderr)
        return 1
    n = len(kernel_list) if kernel_list is not None else len(kept)
    line = (f"Stopped the dashboard server (pid {pid}, port "
            f"{record.get('port')}, up {up}): {n} dashboard{'s' * (n != 1)} "
            "ended")
    if kept:
        line += f" (kept: {', '.join(kept)})"
    print(line + ".")
    if agents:
        handed, pending = _await_hand_over(agents, _HANDOVER_WAIT_S)
        if handed:
            print(f"  {len(handed)} kept agent session"
                  f"{'s' * (len(handed) != 1)} went on in the agent service "
                  "(delfin-agent tui, or a new dashboard, opens them).")
        if pending:
            print(f"  {len(pending)} kept agent session"
                  f"{'s' * (len(pending) != 1)} not handed to the agent "
                  "service yet; saved, and tried again by the next dashboard "
                  f"(log: {_kas.record_dir() / 'hand_over.log'}).")
    return 0


def _server_command(name: str, rest) -> int:
    parser = argparse.ArgumentParser(
        prog=f"delfin-voila {name}",
        description=("Show this host's dashboard server, what it keeps and "
                     "how to reach it." if name == "status" else
                     "End this host's dashboard server the way Ctrl+C did: "
                     "its dashboards end, kept agent sessions go on in the "
                     "agent service."))
    parser.parse_args(list(rest))
    return _status() if name == "status" else _stop()


def _record_when_ready(proc, args, *, token: str, url: str, scheme: str,
                       resume_path: str):
    """Write this host's server record once the server answers, from a
    thread (the launcher waits for its server meanwhile); register its
    removal for when the launcher ends. Returns the thread."""
    from delfin.dashboard import server_record as _sr

    path = _sr.record_file()
    atexit.register(_sr.withdraw, path)
    ip = str(args.ip or "127.0.0.1")
    probe = "127.0.0.1" if ip in ("0.0.0.0", "127.0.0.1", "localhost") else ip

    def _run():
        if not _wait_for_port(probe, args.port, timeout=DETACHED_START_TIMEOUT_S):
            return
        try:
            from delfin import __version__ as version
        except Exception:
            version = ""
        record = _sr.own_record()
        record.update({
            "server_pid": int(getattr(proc, "pid", 0) or 0),
            "port": int(args.port), "ip": ip, "scheme": scheme,
            "url": url, "resume_path": resume_path, "token": token,
            "log": str(_sr.log_file()) if _sr.is_detached() else "",
            "detached": _sr.is_detached(), "started_at": time.time(),
            "version": str(version),
        })
        try:
            _sr.write(record)
        except OSError as exc:
            print(f"Note: the server record could not be written ({exc}); "
                  "a second delfin-voila will not find this server.",
                  file=sys.stderr)

    thread = threading.Thread(target=_run, name="delfin-server-record",
                              daemon=True)
    thread.start()
    return thread


def _announces_where(env) -> bool:
    """False for a branch preview (branch_preview.preview_env sets the
    markers): the note says where *the* dashboard runs, and a preview
    closes with its branch."""
    return not (env.get("DELFIN_PREVIEW_BRANCH") and env.get("DELFIN_PREVIEW_WORKTREE"))


def _parse_args(argv):
    """The launcher's options (``status`` and ``stop`` are read before)."""
    parser = argparse.ArgumentParser(
        prog="delfin-voila",
        description=(
            "Launch the DELFIN Dashboard with Voila. Typed at a terminal, "
            "the server runs in the background and outlives the terminal; "
            "a second delfin-voila opens the one that runs. "
            "`delfin-voila status` shows it, `delfin-voila stop` ends it."),
    )
    from delfin import __version__ as _delfin_version
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {_delfin_version}",
        help="Print the DELFIN version and exit.",
    )
    parser.add_argument(
        "--strict-port",
        action="store_true",
        help="Fail if the requested port is occupied instead of choosing another (SSH launchers).",
    )
    parser.add_argument(
        "--agent-terminal",
        metavar="COMMAND",
        default="",
        help=(
            "Offer COMMAND (a coding-agent CLI, by name on PATH or as a "
            "path) in a terminal panel of the Agent tab. Loopback only; the "
            "terminal runs that program and nothing else, and never a shell "
            "or an interpreter."
        ),
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8866,
        help="Port to serve on (default: 8866)",
    )
    browser_group = parser.add_mutually_exclusive_group()
    browser_group.add_argument(
        "--open-browser",
        action="store_true",
        default=None,
        help="Ask Voila to open a browser window",
    )
    browser_group.add_argument(
        "--no-browser",
        action="store_true",
        default=None,
        help="Disable automatic browser launch",
    )
    parser.add_argument(
        "--dark",
        action="store_true",
        help="Use dark theme",
    )
    parser.add_argument(
        "--ip",
        default=None,
        help=(
            "IP to bind to (default: 127.0.0.1 for local-only access; "
            "set 0.0.0.0 only when direct network access is required)"
        ),
    )
    parser.add_argument(
        "--allow-remote-bind",
        action="store_true",
        help=(
            "Allow binding to a non-loopback IP. Not recommended; prefer SSH "
            "tunnels or an authenticated HTTPS reverse proxy on shared servers."
        ),
    )
    parser.add_argument(
        "--certfile",
        default=None,
        help=(
            "TLS certificate (PEM) for HTTPS. Recommended together with "
            "--keyfile when serving to a network so the token is not sent in "
            "clear."
        ),
    )
    parser.add_argument(
        "--keyfile",
        default=None,
        help="TLS private key (PEM) matching --certfile.",
    )
    parser.add_argument(
        "--token",
        default=None,
        help=(
            "Authentication token for dashboard access. "
            "Auto-generated if not set."
        ),
    )
    parser.add_argument(
        "--no-token",
        action="store_true",
        help="Refused for security; token authentication is mandatory.",
    )
    parser.add_argument(
        "--resume",
        default="",
        metavar="SID",
        help=(
            "Auto-load a saved agent session on dashboard boot. Pass a "
            "specific session ID or 'latest' to resume the most recent. "
            "Equivalent to setting the DELFIN_RESUME_SESSION env var."
        ),
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help=(
            "Arm 'Keep dashboard' for every window this dashboard opens, "
            "so a closed window or a dropped connection costs nothing. "
            "Surviving the terminal needs no flag any more (no tmux)."
        ),
    )
    parser.add_argument(
        "--foreground",
        action="store_true",
        help=(
            "Run the server in this terminal, as before: Ctrl+C or closing "
            "the terminal ends it. For debugging; a launch without a "
            "terminal (a script, systemd, launchd, a SLURM job) and one "
            "inside tmux or screen stays in the foreground anyway. "
            "DELFIN_VOILA_FOREGROUND=1 makes it the default."
        ),
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List your dashboard jobs (node, time left, kept sessions).")
    parser.add_argument(
        "--attach", nargs="?", const="latest", default=None, metavar="JOB",
        help="Reconnect to a running dashboard job (the newest without JOB).")
    parser.add_argument(
        "--stop", nargs="?", const="latest", default=None, metavar="JOB",
        help="End a dashboard job (JOB, or 'all').")
    parser.add_argument(
        "--job", action="store_true",
        help=("Run the dashboard as a SLURM job on a compute node and tunnel "
              "to it (reconnects to a running one). For more than the one "
              "agent session a login node allows."))
    parser.add_argument(
        "--new", action="store_true",
        help="With --job: start another dashboard job although one runs.")
    return parser.parse_args(argv)


def main(argv=None):
    raw = list(argv) if argv is not None else sys.argv[1:]
    if raw and raw[0] in ("status", "stop"):
        return _server_command(raw[0], raw[1:])
    args = _parse_args(raw)
    from delfin.dashboard import server_record as _server_record
    if (args.job or args.list or args.attach is not None
            or args.stop is not None):
        # A dashboard in a SLURM job (issue #640): any number of agent
        # sessions, reached through a tunnel, managed by these flags. On a
        # login node the dashboard itself runs here as before, with ONE
        # agent session (agent_sessions.login_node_refusal) whose tests and
        # builds go to SLURM.
        from delfin import dashboard_job
        return dashboard_job.main(args, raw)
    if getattr(args, "keep", False):
        # It started tmux once; a dashboard outlives its terminal without
        # it now. What stays is the arming: the sessions this dashboard
        # opens start with Keep dashboard on. Said by the launch a person
        # typed, not again by the detached launcher it starts.
        if not _server_record.is_detached():
            print(KEEP_NOTE, file=sys.stderr)
        os.environ["DELFIN_KEEP_SESSIONS"] = "1"
    # This process carries the user's environment, which may hold an
    # exported provider key; the commands the dashboard runs are processes
    # of the same user. Unreadable to them, and said once if a key is there.
    try:
        from delfin.agent import process_guard as _process_guard
        _process_guard.protect("dashboard launcher")
        from delfin.user_settings import load_settings as _load_settings, take_security_notices
        _load_settings()                    # applies a pending security update
        for _notice in take_security_notices():
            print(_notice, file=sys.stderr)
        _exported = _process_guard.exported_provider_keys()
        # Unconditional, because the environment is not where the worst
        # case lives: an alias sets a key for one command and leaves the
        # line in a shell file for good, so nothing ever exports it and
        # the old `if _exported:` never looked.
        #
        # Do it rather than ask for it: an exported key goes into the
        # 0600 store and the export is commented out, with a copy of the
        # file kept. Only what could not be decided is left to say.
        _done = _process_guard.put_exported_keys_away(_exported)
        if _done:
            print(_done, file=sys.stderr)
        elif _exported:
            print("Warning: "
                  + _process_guard.exported_key_advice(_exported),
                  file=sys.stderr)
    except Exception:
        pass
    # Record the REAL shell cwd the user launched from, BEFORE Voila resets the
    # kernel's cwd to the notebook's directory (inside the delfin checkout).
    # ONLY the agent tab reads this — to decide where to build (launch dir =
    # workspace). Everything else keeps using ctx.repo_dir unchanged.
    os.environ.setdefault("DELFIN_LAUNCH_CWD", os.getcwd())
    if args.resume:
        # Propagate to the spawned voila subprocess via env; the agent-tab
        # reads DELFIN_RESUME_SESSION at create_tab() time.
        os.environ["DELFIN_RESUME_SESSION"] = args.resume

    # Check that voila is installed in the current Python environment.
    if not _voila_is_available():
        print(
            "Error: voila is not installed.\n"
            f"Install it with:  {sys.executable} -m pip install voila",
            file=sys.stderr,
        )
        sys.exit(1)

    # One server per user and host: a running one is opened, not joined
    # by a second (server_record). The detached launcher itself does not
    # ask: the launch that started it just did, and the lock below
    # settles a race between two.
    _single = _single_server_applies(os.environ)
    if _single and not _server_record.is_detached():
        _running = _server_record.running()
        if _running is not None and _server_record.port_answers(_running):
            return _use_running(_running, args, raw)
    _detach, _why_not = _should_detach(args)
    if _detach:
        return _start_detached(args, raw)
    _server_lock = None
    if _single:
        _server_lock = _server_record.Lock()
        if not _server_lock.acquire():
            # Another launcher on this host holds it: starting, or running.
            _other = _server_record.wait_for_record(timeout=60.0)
            if _other is not None:
                return _use_running(_other, args, raw)
            print("Error: another dashboard launcher holds this host's lock "
                  f"({_server_record.lock_file()}) and its server does not "
                  "answer. `delfin-voila status` shows it.", file=sys.stderr)
            return 1

    # Find a free port, starting from the requested one.
    try:
        selected_port = _select_port(args.port)
        if args.strict_port and selected_port != args.port:
            raise RuntimeError(f"Error: requested port {args.port} is occupied; refusing to change the SSH tunnel target.")
        args.port = selected_port
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    if args.ip is None:
        args.ip = "127.0.0.1"
    if not _is_loopback_bind(args.ip) and not args.allow_remote_bind:
        print(
            "Error: refusing to bind delfin-voila to a non-loopback address "
            f"({args.ip}). On HPC/login nodes, use the default 127.0.0.1 with "
            "an SSH tunnel. If you have a secured HTTPS reverse proxy and "
            "really need this, pass --allow-remote-bind.",
            file=sys.stderr,
        )
        sys.exit(2)

    # -- Agent terminal (opt-in) ----------------------------------------------
    # A terminal is a process under the user's account, outside the agent's
    # sandbox. Only on a loopback bind -- never with --allow-remote-bind --
    # only when the program exists, and never a shell or an interpreter
    # (agent_terminal.refusal): no command, no terminal.
    _term_argv: list[str] = []
    if args.agent_terminal:
        from delfin.dashboard import agent_terminal as _at
        if not _is_loopback_bind(args.ip):
            print("Error: --agent-terminal is only offered on a loopback "
                  f"bind, not on {args.ip}.", file=sys.stderr)
            sys.exit(2)
        _wanted = os.path.expanduser(args.agent_terminal)
        _found = (_wanted if os.path.isfile(_wanted)
                  and os.access(_wanted, os.X_OK) else shutil.which(_wanted))
        if not _found:
            print(f"Error: --agent-terminal: {args.agent_terminal!r} is not "
                  "an executable on PATH or a path to one.", file=sys.stderr)
            sys.exit(2)
        _found = os.path.realpath(os.path.abspath(_found))
        _why = _at.refusal(_found) or _at.refusal(args.agent_terminal)
        if _why:
            print(f"Error: --agent-terminal: {_why}.", file=sys.stderr)
            sys.exit(2)
        _term_argv = [_found]

    # -- Token-based authentication ----------------------------------------
    # Auto-generate a token by default. On shared multi-user hosts (HPC login
    # nodes), even 127.0.0.1 is reachable by other local users, so a token is
    # required regardless of bind address. --no-token is refused above.
    _token = ""
    _token_file: Path | None = None
    if args.no_token:
        print(_refuse_insecure_no_token(args.ip), file=sys.stderr)
        sys.exit(2)
    elif args.token or os.environ.get("DELFIN_VOILA_TOKEN"):
        # DELFIN_VOILA_TOKEN: set by a dashboard job's launcher, so the token
        # is in no command line; not passed on to the kernels.
        try:
            _token = _validate_access_token(
                args.token or os.environ.pop("DELFIN_VOILA_TOKEN"))
        except ValueError as exc:
            print(f"Error: insecure --token: {exc}.", file=sys.stderr)
            sys.exit(2)
    else:
        import atexit
        import secrets
        _token = secrets.token_urlsafe(32)
        # Write token to 0600 file instead of stdout — stdout may be captured
        # by journald/log files where it would be readable by other users on
        # shared hosts. XDG_RUNTIME_DIR is per-user (mode 0700) on Linux; we
        # fall back to a private DELFIN cache dir if it is absent or unwritable.
        _token_file = _create_token_file(_token)

        def _cleanup_token_file(p: Path = _token_file) -> None:
            try:
                p.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass

        atexit.register(_cleanup_token_file)
        print(
            f"\n  🔑 Auto-generated access token (mode 0600):\n"
            f"     {_token_file}\n"
            f"     Read it with:  cat {_token_file}\n"
        )

    open_browser = _decide_open_browser(args)
    notebook = _find_notebook()
    root_dir = str(_default_voila_root())
    notebook = _stage_notebook_under_root(notebook, root_dir)
    env = _prepare_voila_env(open_browser=open_browser)
    env.setdefault("DELFIN_VOILA_ROOT_DIR", root_dir)
    env["DELFIN_VOILA_PORT"] = str(args.port)

    # Quieter startup (cosmetic only — no endpoint/auth change). Suppress the
    # third-party DeprecationWarning spam (traitlets) and the Node url.parse
    # deprecation so the terminal shows the DELFIN banner, not a wall of noise.
    # Scoped to the Voilà subprocess env; DELFIN's own runtime is untouched.
    env.setdefault("PYTHONWARNINGS", "ignore::DeprecationWarning")
    _node_opts = env.get("NODE_OPTIONS", "").strip()
    env["NODE_OPTIONS"] = (
        f"{_node_opts} --no-deprecation".strip() if _node_opts else "--no-deprecation"
    )

    # Trust the notebook so Voilà doesn't warn about untrusted content.
    _trust_notebook(notebook)

    # Launch Voilà as a JUPYTER-SERVER EXTENSION, not standalone `voila`.
    # WHY (verified 2026-07-08, voila 0.5.12 / jupyter_server 2.20): standalone
    # `voila <nb>` does NOT enforce the token — the dashboard page, the kernel
    # REST API AND the kernel websocket are all reachable WITHOUT auth, so on a
    # shared host any local user reaching the port could drive the agent tab
    # (arbitrary CLI = RCE as the launching user). Under `jupyter server` with
    # the voila extension, jupyter_server enforces the token instead: an
    # unauthenticated GET -> 302 login and the kernel websocket -> 403.
    # `default_url` sends the browser straight to the rendered dashboard, so
    # `http://localhost:PORT/?token=…` behaves exactly as before for the user.
    dashboard_url = _dashboard_default_url(notebook, root_dir)

    # Coming back to a kept session.
    #
    # A session the user chose to keep is reached through a SECOND,
    # one-cell notebook rather than the dashboard's own address: that one
    # renders and executes the whole thing, which against a live kernel
    # would run all nineteen tabs a second time. The one-cell notebook
    # re-displays what the kernel already holds.
    #
    # The path is handed to the kernels through the environment because
    # only the launcher knows where the notebook was staged relative to
    # the server root; and the kernel manager is swapped for one that
    # recognises `?session=<name>` and hands back the running kernel
    # instead of starting a fresh one. An ordinary request goes through
    # it untouched.
    try:
        from delfin.dashboard import resume_server as _resume

        _resume_nb = _resume.stage_resume_notebook(root_dir)
        _trust_notebook(_resume_nb)
        resume_url_path = _dashboard_default_url(_resume_nb, root_dir)
        kernel_manager_class = (
            "delfin.dashboard.resume_server.ResumeAwareMappingKernelManager")
    except Exception:
        # A dashboard that starts is worth more than one that can be
        # resumed: without these two the default path is exactly what it
        # was, and keeping a session simply is not offered.
        resume_url_path = ""
        kernel_manager_class = ""

    # Kernels inherit the server's environment, and Voila copies it into
    # every kernel it starts, so this is how a dashboard learns its own
    # come-back address without the launcher and the kernel having to
    # agree on paths twice.
    if resume_url_path:
        env[_resume.RESUME_PATH_ENV] = resume_url_path
    # ... and the dashboard's own address, for the header's "New
    # dashboard" (dashboards_menu).
    from delfin.dashboard.session import DASHBOARD_PATH_ENV as _dash_path_env
    env[_dash_path_env] = dashboard_url
    # Voilà only needs the `voila` extension. Explicitly DISABLE the others that
    # auto-load from entry points (jupyterlab/notebook/lsp/terminals): each adds
    # HTTP endpoints = attack surface, boot time and log noise the dashboard
    # never uses. notebook_shim is left alone (harmless redirect shim).
    _extensions = (
        "{'voila': True, 'jupyterlab': False, 'notebook': False, "
        "'jupyter_lsp': False, 'jupyter_server_terminals': "
        + ("True" if _term_argv else "False") + ", "
        "'delfin.dashboard.server_guard': True}"
    )
    if _term_argv:
        import json as _json
        from delfin.dashboard import agent_terminal as _ct
        env[_ct.ARGV_ENV] = _json.dumps(_term_argv)
        env[_ct.CWD_ENV] = os.environ.get("DELFIN_LAUNCH_CWD") or os.getcwd()
        env[_ct.ENABLED_ENV] = "1"
    else:
        env.pop("DELFIN_AGENT_TERMINAL", None)
    import json

    from delfin.dashboard.web_security import NO_TYPESETTING_ROUTE

    # Server-side login sessions (logout, idle timeout) and refused logins
    # as security events: delfin.dashboard.identity. If it does not import
    # with this jupyter-server, the server starts with jupyter-server's own
    # provider and the terminal says what is missing.
    _identity_args = []
    try:
        from delfin.dashboard import identity as _identity
        _identity_args = [
            "--ServerApp.identity_provider_class=" + _identity.CLASS_PATH]
    except Exception as _identity_exc:  # noqa: BLE001
        print("Warning: dashboard logins have no logout or idle timeout on "
              f"this server ({_identity_exc}); they end only after "
              "DELFIN_DASHBOARD_LOGIN_DAYS.", file=sys.stderr)

    cmd = [
        sys.executable,
        "-m",
        "jupyter",
        "server",
        f"--ServerApp.jpserver_extensions={_extensions}",
        *_identity_args,
        *(["--ServerApp.kernel_manager_class="
           + kernel_manager_class] if kernel_manager_class else []),
        # A kernel whose last window closed is ended by the server after
        # a grace, unless the session was kept. The flags switch the
        # culler on; the rule itself lives in the kernel manager above.
        *(_resume.cull_config_args() if kernel_manager_class else []),
        *_kernel_transport_args(),
        f"--port={args.port}",
        *(["--ServerApp.port_retries=0"] if args.strict_port else []),
        f"--ServerApp.ip={args.ip}",
        f"--ServerApp.root_dir={root_dir}",
        f"--ServerApp.default_url={dashboard_url}",
        # Quieter startup: drop the [I …] link/load/serving/302 flood; WARN
        # still surfaces genuine warnings and errors (fail-loud preserved).
        "--ServerApp.log_level=WARN",
        # The kernel-to-browser channel is rate limited for a notebook, where a
        # runaway print loop can hang the page. The defaults are 1000 messages
        # and 1 MB per three seconds, and over them the server *stops sending*
        # -- so a widget silently keeps its old value and the tab looks stuck
        # with nothing in the log but "IOPub message rate exceeded". A dashboard
        # is the opposite case: a file browser showing a photograph, a PDF page
        # or a five-thousand-row grid sends several MB in one message on
        # purpose, and dropping it is worse than a busy browser. The window
        # stays, so a genuine runaway is still bounded, an order of magnitude
        # further out.
        "--ZMQChannelsWebsocketConnection.iopub_msg_rate_limit=10000",
        "--ZMQChannelsWebsocketConnection.iopub_data_rate_limit=200000000",
        # Info-disclosure: never surface Python tracebacks to the browser
        # (they leak paths, usernames, code). Off by default; set explicitly.
        "--VoilaConfiguration.show_tracebacks=False",
        # What the UI loads from /voila/files/: images, PDFs and the drawing
        # editor's launch copy (web_security.FILE_ALLOWLIST). The route itself
        # is put behind authentication by the server_guard extension, which
        # also serves anything active from user data as a download.
        "--VoilaConfiguration.file_allowlist=" + _file_allowlist(),
        "--VoilaConfiguration.preheat_kernel=False",
        "--VoilaConfiguration.default_pool_size=0",
        # Shrink the RCE surface: no server terminals, unless the agent
        # terminal was asked for -- and then only through the manager that
        # runs the one configured program and nothing else
        # (agent_terminal.py).
        *(["--ServerApp.terminals_enabled=True",
            "--TerminalsExtensionApp.terminal_manager_class="
            "delfin.dashboard.agent_terminal.OneCommandTerminalManager"]
          if _term_argv else ["--ServerApp.terminals_enabled=False"]),
        # Ctrl+C stops immediately instead of prompting "Shut down (y/[n])?".
        "--ServerApp.answer_yes=True",
        "--ServerApp.websocket_ping_interval=30000",
        "--ServerApp.websocket_ping_timeout=30000",
        # MathJax's default entry point is a loader shim, and it never stops
        # loading: it removes and re-appends its own script tag for as long as
        # the page is open. Measured on the running dashboard, five idle
        # seconds: 4.90 s of the main thread busy, against 0.06 s once the
        # request is answered. That is 98% of one core held by a page doing
        # nothing, which is what every later click had to queue behind.
        #
        # This asks for nothing at all rather than for a working MathJax,
        # because nothing is what the dashboard has today: window.MathJax is
        # undefined here and no output is typeset. Pointing at a real bundle
        # would fix the spin and *start* typesetting where there is a network
        # -- a $ in a job log becoming mathematics is a change in what the
        # dashboard shows, and this is a speed fix. The rendered document is
        # identical either way: 7626 elements, 15 script tags.
        #
        # The URL is an empty script on the dashboard's own origin
        # (web_security.NO_TYPESETTING_ROUTE, installed by server_guard), so
        # the page CSP's script-src needs no data:. Voila appends
        # ?config=..., which the route ignores. It goes through
        # tornado_settings because the dashboard runs Voila as a
        # jupyter-server extension, where the Voila app's own mathjax_url
        # trait is never read.
        "--ServerApp.tornado_settings="
        + json.dumps({"mathjax_url": NO_TYPESETTING_ROUTE}),
        # Frontend extensions the browser fetches, parses and evaluates at
        # every load, for widgets this dashboard never makes.
        #
        # nglview: 4.55 MB. The structure viewers are py3Dmol and there is no
        # `import nglview` in the source; the four mentions of the name are
        # Python-side package probes for the Settings tab, which this does not
        # touch. Measured over four interleaved loads each way: blocking 1581
        # -> 1227 ms, transfer 14.7 -> 10.2 MB.
        #
        # jupyterlab-plotly: 192.9 ms of self time in a CPU profile of the
        # load, 186.4 of it in one function -- the extension evaluating
        # itself. It is only reachable by a plotly figure rendered as a widget
        # in this page, and nothing can make one: the agent runs its Python in
        # subprocesses (sys.executable -c ...), whose output is a file or
        # text, and Voila has no editable cells. `delfin/tools/` contains no
        # plotly, FigureWidget or to_html. The Settings probe that reports
        # whether the plotly *package* is installed is a different question
        # and still answers it.
        #
        # DELFIN_VOILA_EXTENSION_DENYLIST overrides the list, so a machine
        # that does want either of them back can have it without editing this
        # file.
        f"--VoilaConfiguration.extension_denylist={_extension_denylist()!r}",
    ]

    import json

    # DELFIN-branded login page: drop our own login.html (logo inlined) onto
    # extra_template_paths so jupyter_server renders it instead of its generic
    # "Password or token" page. Appearance only — the LoginHandler still
    # validates the token/password server-side, so auth is unchanged.
    _login_dir = _prepare_login_template_dir()
    if _login_dir:
        cmd.append(
            "--ServerApp.extra_template_paths=" + json.dumps([_login_dir])
        )

    # The page template: Voila's lab template without its CDN links to
    # Font Awesome (refused by the page CSP; the font is bundled locally).
    cmd.extend(_page_template_args(env))

    # The drawing editor, published from where it is kept instead of being
    # copied into the directory the dashboard was started in. The dashboard
    # only takes this route when the variable says the route is there, so a
    # dashboard started any other way still finds its editor.
    _published_dir = _prepare_published_dir()
    if _published_dir:
        cmd.append(
            "--ServerApp.extra_static_paths=" + json.dumps([_published_dir])
        )
        env["DELFIN_KETCHER_URL"] = "/static/ketcher"

    # Token authentication — now actually ENFORCED by jupyter_server. Pass via
    # env var, not CLI: /proc/PID/cmdline is world-readable on shared hosts, so
    # a CLI token would leak to any local user via `ps`. Because the token is
    # env-supplied, jupyter_server does NOT echo it to stdout.
    if _token:
        env["JUPYTER_TOKEN"] = _token
    else:
        # Unreachable (main refuses to run without a token) — keep the
        # fail-closed default explicit anyway.
        cmd.append("--ServerApp.token=")

    # TLS: honour explicit certs so a networked bind does not send the token in
    # clear. Without certs we keep working but warn (SSH tunnel stays the
    # recommended path); loopback (the default) needs no TLS.
    if args.certfile and args.keyfile:
        cmd.append(f"--ServerApp.certfile={args.certfile}")
        cmd.append(f"--ServerApp.keyfile={args.keyfile}")
    elif not _is_loopback_bind(args.ip):
        print(
            "  ⚠  Remote bind without TLS: dashboard traffic (including the "
            "access token) is UNENCRYPTED. Prefer an SSH tunnel, or pass "
            "--certfile/--keyfile for HTTPS.",
            file=sys.stderr,
        )

    # We ALWAYS open the browser ourselves (below), with the token embedded so
    # the runner lands directly in the dashboard. Letting jupyter open it too
    # would race us and could open a token-less URL (-> login page), so disable
    # jupyter's own auto-open unconditionally.
    cmd.append("--no-browser")

    if args.dark:
        cmd.append("--VoilaConfiguration.theme=dark")

    scheme = "https" if (args.certfile and args.keyfile) else "http"
    bind_display = "localhost" if args.ip == "127.0.0.1" else args.ip
    print(f"Starting DELFIN Dashboard on {scheme}://{bind_display}:{args.port}")

    # The token is what keeps everyone else out: opening the port WITHOUT it
    # lands on the login page, so no other local user on this shared host can
    # get in. On an interactive terminal (YOUR own private session) show the
    # ready-to-open URL with the token embedded, so you are logged in with a
    # single click. When stdout is redirected/captured (service logs, job
    # launchers) isatty() is False and we fall back to the token-file path — the
    # raw token is never written to a shared log.
    # Point straight at the dashboard render URL (not "/"). Going through "/"
    # relies on a redirect to default_url that can DROP the token when the
    # browser still holds a stale login cookie for a reused port (localhost:PORT
    # accumulates cookies across sessions, each with a fresh cookie-secret) —
    # the token then never reaches the target and you land on /login. Hitting
    # the render URL directly consumes the token at the target handler, which
    # sets a fresh cookie regardless of any stale one.
    _show_url = bool(_token) and sys.stdout.isatty()
    if _show_url and args.ip != "0.0.0.0":
        print(
            f"  🔗 Open (you are logged in directly): "
            f"{scheme}://{bind_display}:{args.port}{dashboard_url}?token={_token}"
        )
    elif _token and _token_file is not None:
        print(f"  Append ?token=$(cat {_token_file}) to the URL.")

    # Show all reachable URLs so remote users know exactly what to type.
    if args.ip == "0.0.0.0":
        for ip in _get_local_ips():
            suffix = f"{dashboard_url}?token={_token}" if _show_url else ""
            print(f"  -> {scheme}://{ip}:{args.port}{suffix}")

    print()
    if args.ip == "0.0.0.0":
        hostname = socket.gethostname().split(".")[0]
        ips = _get_local_ips()
        target = ips[0] if ips else hostname
        # Detect the login node from SSH_CONNECTION (the IP the user connected to).
        login_node = _detect_login_node()
        user = os.environ.get("USER", os.environ.get("LOGNAME", "<user>"))
        tunnel_target = f"{user}@{login_node}" if login_node else "<user>@<login-node>"
        print(
            "Tip: If you cannot connect directly, set up an SSH tunnel.\n"
            "     Run this on your LOCAL machine (PowerShell / MobaXterm / Terminal):\n"
            "\n"
            f"       ssh -N -L {args.port}:{target}:{args.port} {tunnel_target}\n"
            "\n"
            f"     Then open http://localhost:{args.port} in your browser."
        )
    # Leave the note that says where this is running, so a later login
    # on another node can find its way back instead of guessing which
    # node it was. Best-effort: a dashboard that starts is worth more
    # than one that can be found again.
    try:
        from delfin.agent import where as _where

        if not _announces_where(os.environ):
            raise RuntimeError("a branch preview leaves the note alone")
        _where.announce_dashboard(
            port=args.port, token=_token, url=dashboard_url,
            resume_path=resume_url_path.split("?")[0] if resume_url_path else "")
        atexit.register(_where.withdraw_dashboard)
        _back = _where.reconnect_command(_where.dashboard())
        if _back:
            print(f"  ↩ Come back with:  {_back}")
            print("     From anywhere:     delfin-agent where\n")
    except Exception:
        pass
    if _server_record.is_detached():
        print("Detached from the terminal; `delfin-voila stop` ends it.\n")
    else:
        print("Press Ctrl+C to stop.\n")

    # Pipe stderr so we can drop the few benign jupyter noise lines that survive
    # Python's warning filters; stdout stays inherited (banner is untouched).
    # This process is the lifeline of everything the dashboard starts: the
    # server dies with it however it dies (parent death signal), kernels and
    # daemons watch it, and what runs detached is on its ledger.
    from delfin.agent import lifeline as _lifeline
    env.update(_lifeline.claim_root())
    # The private runtime directory _kernel_transport_args chose, when the
    # default one was too long for a socket or not private: env was copied
    # from os.environ before the choice was made.
    if os.environ.get("JUPYTER_RUNTIME_DIR"):
        env["JUPYTER_RUNTIME_DIR"] = os.environ["JUPYTER_RUNTIME_DIR"]
    proc = subprocess.Popen(
        cmd, env=env, stderr=subprocess.PIPE,
        preexec_fn=_lifeline.parent_death_signal if os.name == "posix" else None,
    )
    _previous_signals = _install_stop_signals()
    if _single:
        _record_when_ready(
            proc, args, token=_token, url=dashboard_url, scheme=scheme,
            resume_path=resume_url_path.split("?")[0] if resume_url_path else "")
    stderr_stream = getattr(proc, "stderr", None)
    stderr_thread = None
    if stderr_stream is not None:
        stderr_thread = threading.Thread(
            target=_forward_filtered_stderr, args=(stderr_stream,), daemon=True
        )
        stderr_thread.start()
    # Auto-open straight INTO the dashboard: open the render URL with the token
    # embedded (not bare "/", which would land on the login page). Only the
    # runner gets this URL; anyone else opening the port still hits the login.
    if open_browser and _token:
        browser_host = "localhost" if args.ip in {"127.0.0.1", "0.0.0.0"} else args.ip
        direct_url = (
            f"{scheme}://{browser_host}:{args.port}{dashboard_url}?token={_token}"
        )
        if _wait_for_port("127.0.0.1", args.port):
            _open_browser_url(direct_url, env)
    try:
        rc = proc.wait()
        if stderr_thread is not None:
            stderr_thread.join(timeout=2)
        _end_agent_processes()
        sys.exit(rc)
    except KeyboardInterrupt:
        # Ctrl+C, SIGTERM or the terminal closing: the end of everything,
        # kept sessions included. A second Ctrl+C must not cut it short.
        try:
            import signal as _signal
            _signal.signal(_signal.SIGINT, _signal.SIG_IGN)
        except Exception:
            pass
        print("\nStopping the dashboard and everything it started...")
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        if stderr_thread is not None:
            stderr_thread.join(timeout=2)
        _end_agent_processes()
        print("Stopped.")
    finally:
        _restore_signals(_previous_signals)
        if _single:
            from delfin.dashboard import server_record as _sr_end
            _sr_end.withdraw(_sr_end.record_file())
            if _server_lock is not None:
                _server_lock.release()


if __name__ == "__main__":
    main()
