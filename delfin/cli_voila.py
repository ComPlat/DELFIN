"""``delfin-voila`` – launch the DELFIN Dashboard via Voila.

Usage::

    delfin-voila               # start on default port 8866
    delfin-voila --port 9000   # custom port
    delfin-voila --open-browser
"""

import argparse
import importlib.util
import importlib.resources
import ipaddress
import os
import socket
import shutil
import shlex
import subprocess
import sys
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


def _is_vscode_session(env: dict[str, str] | None = None) -> bool:
    """Return True when the launcher runs inside a VS Code terminal session."""
    env = env or os.environ
    term_program = str(env.get("TERM_PROGRAM") or "")
    browser = str(env.get("BROWSER") or "")
    ipc_hook = str(env.get("VSCODE_IPC_HOOK_CLI") or "").strip()
    return term_program == "vscode" or "browser.sh" in browser or bool(ipc_hook)


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


def _default_voila_root() -> Path:
    """Return a narrow default root_dir for files served by Voilà."""
    override = os.environ.get("DELFIN_VOILA_ROOT_DIR", "").strip()
    if override:
        return Path(override).expanduser().resolve()

    cwd = Path.cwd().resolve()
    home = Path.home().resolve()
    forbidden = {
        Path("/").resolve(),
        Path("/home").resolve(),
        Path("/tmp").resolve(),
        home,
    }
    if cwd not in forbidden and home in cwd.parents:
        return cwd
    return _secure_cache_dir("voila-root").resolve()


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


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="delfin-voila",
        description="Launch the DELFIN Dashboard with Voila.",
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
    args = parser.parse_args(argv)
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

    # Find a free port, starting from the requested one.
    try:
        args.port = _select_port(args.port)
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

    # -- Token-based authentication ----------------------------------------
    # Auto-generate a token by default. On shared multi-user hosts (HPC login
    # nodes), even 127.0.0.1 is reachable by other local users, so a token is
    # required regardless of bind address. --no-token is refused above.
    _token = ""
    _token_file: Path | None = None
    if args.no_token:
        print(_refuse_insecure_no_token(args.ip), file=sys.stderr)
        sys.exit(2)
    elif args.token:
        try:
            _token = _validate_access_token(args.token)
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

    if args.open_browser is True:
        open_browser = True
    elif args.no_browser is True:
        open_browser = False
    else:
        open_browser = _is_vscode_session()
    notebook = _find_notebook()
    root_dir = str(_default_voila_root())
    notebook = _stage_notebook_under_root(notebook, root_dir)
    env = _prepare_voila_env(open_browser=open_browser)
    env.setdefault("DELFIN_VOILA_ROOT_DIR", root_dir)
    env["DELFIN_VOILA_PORT"] = str(args.port)

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
    cmd = [
        sys.executable,
        "-m",
        "jupyter",
        "server",
        "--ServerApp.jpserver_extensions={'voila': True}",
        f"--port={args.port}",
        f"--ServerApp.ip={args.ip}",
        f"--ServerApp.root_dir={root_dir}",
        f"--ServerApp.default_url={dashboard_url}",
        # Info-disclosure: never surface Python tracebacks to the browser
        # (they leak paths, usernames, code). Off by default; set explicitly.
        "--VoilaConfiguration.show_tracebacks=False",
        "--VoilaConfiguration.file_allowlist=.*\\.(png|jpg|gif|svg|js|css|html|ico|pdf)",
        "--VoilaConfiguration.preheat_kernel=False",
        "--VoilaConfiguration.default_pool_size=0",
        # Shrink the RCE surface: the dashboard never needs server terminals.
        "--ServerApp.terminals_enabled=False",
        "--ServerApp.websocket_ping_interval=30000",
        "--ServerApp.websocket_ping_timeout=30000",
    ]

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

    if open_browser:
        cmd.append("--ServerApp.open_browser=True")
    else:
        cmd.append("--no-browser")

    if args.dark:
        cmd.append("--VoilaConfiguration.theme=dark")

    scheme = "https" if (args.certfile and args.keyfile) else "http"
    bind_display = "localhost" if args.ip == "127.0.0.1" else args.ip
    # Never print the raw token. Terminals, notebooks, service logs, and job
    # launchers may retain stdout on shared systems.
    print(f"Starting DELFIN Dashboard on {scheme}://{bind_display}:{args.port}")
    if _token and _token_file is not None:
        print(f"  Append ?token=$(cat {_token_file}) to the URL.")

    # Show all reachable URLs so remote users know exactly what to type.
    if args.ip == "0.0.0.0":
        for ip in _get_local_ips():
            print(f"  -> {scheme}://{ip}:{args.port}")

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
    print("Press Ctrl+C to stop.\n")

    proc = subprocess.Popen(cmd, env=env)
    if open_browser and _is_vscode_session(env):
        browser_host = "localhost" if args.ip in {"127.0.0.1", "0.0.0.0"} else args.ip
        if _wait_for_port("127.0.0.1", args.port):
            _open_browser_url(f"http://{browser_host}:{args.port}", env)
    try:
        sys.exit(proc.wait())
    except KeyboardInterrupt:
        print("\nStopping...")
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        print("Stopped.")


if __name__ == "__main__":
    main()
