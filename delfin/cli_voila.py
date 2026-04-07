"""``delfin-voila`` – launch the DELFIN Dashboard via Voila.

Usage::

    delfin-voila               # start on default port 8866
    delfin-voila --port 9000   # custom port
    delfin-voila --open-browser
"""

import argparse
import importlib.util
import importlib.resources
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


def _find_notebook() -> str:
    """Locate the dashboard notebook, preferring the checkout copy."""
    search_roots = [Path.cwd(), *Path.cwd().parents, Path.home()]
    rel_candidates = [
        Path("delfin_dashboard.ipynb"),
        Path("software/delfin/delfin_dashboard.ipynb"),
    ]
    for root in search_roots:
        for rel in rel_candidates:
            candidate = (root / rel).resolve()
            if candidate.is_file():
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
        "--token",
        default=None,
        help=(
            "Authentication token for dashboard access. "
            "Auto-generated if not set and --ip is 0.0.0.0. "
            "Use --no-token to explicitly disable."
        ),
    )
    parser.add_argument(
        "--no-token",
        action="store_true",
        help="Disable token authentication (NOT recommended for network access)",
    )
    args = parser.parse_args(argv)

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

    # -- Token-based authentication ----------------------------------------
    # Auto-generate a token when binding to non-localhost (security).
    _token = ""
    if args.no_token:
        _token = ""
        if args.ip != "127.0.0.1":
            print(
                "\n  ⚠  WARNING: Token auth disabled on non-localhost binding!\n"
                "     The DELFIN Agent has FULL CLI ACCESS.\n"
                "     Anyone who can reach this port can execute commands.\n",
                file=sys.stderr,
            )
    elif args.token:
        _token = args.token
    elif args.ip != "127.0.0.1":
        # Auto-generate token for network-facing deployments
        import secrets
        _token = secrets.token_urlsafe(32)
        print(f"\n  🔑 Auto-generated access token (required in URL):\n")
        print(f"     {_token}\n")

    if args.open_browser is True:
        open_browser = True
    elif args.no_browser is True:
        open_browser = False
    else:
        open_browser = _is_vscode_session()
    notebook = _find_notebook()
    root_dir = str(
        Path(
            os.environ.get("DELFIN_VOILA_ROOT_DIR", str(Path.home().resolve()))
        ).resolve()
    )
    notebook = _stage_notebook_under_root(notebook, root_dir)
    env = _prepare_voila_env(open_browser=open_browser)
    env.setdefault("DELFIN_VOILA_ROOT_DIR", root_dir)
    env["DELFIN_VOILA_PORT"] = str(args.port)

    # Trust the notebook so Voilà doesn't warn about untrusted content.
    _trust_notebook(notebook)

    cmd = [
        sys.executable,
        "-m",
        "voila",
        notebook,
        f"--port={args.port}",
        f"--Voila.ip={args.ip}",
        "--show_tracebacks=True",
        f"--Voila.root_dir={root_dir}",
        "--VoilaConfiguration.file_allowlist=.*\\.(png|jpg|gif|svg|js|css|html|ico)",
        "--VoilaConfiguration.preheat_kernel=False",
        "--VoilaConfiguration.default_pool_size=0",
        # XSRF handled per-mode: disabled for localhost, enabled for network
        "--ServerApp.websocket_ping_interval=30000",
        "--ServerApp.websocket_ping_timeout=30000",
    ]

    # Token authentication for Voilà/Jupyter
    if _token:
        cmd.append(f"--ServerApp.token={_token}")
    elif args.ip == "127.0.0.1":
        # Localhost: no token needed, disable for convenience
        cmd.append("--ServerApp.token=")
        cmd.append("--ServerApp.disable_check_xsrf=True")

    voila_static_root = _get_voila_static_root()
    if voila_static_root:
        cmd.append(f"--Voila.static_root={voila_static_root}")

    if open_browser:
        cmd.append("--Voila.open_browser=True")
    else:
        cmd.append("--no-browser")

    if args.dark:
        cmd.append("--theme=dark")

    bind_display = "localhost" if args.ip == "127.0.0.1" else args.ip
    token_suffix = f"?token={_token}" if _token else ""
    print(f"Starting DELFIN Dashboard on http://{bind_display}:{args.port}{token_suffix}")

    # Show all reachable URLs so remote users know exactly what to type.
    if args.ip == "0.0.0.0":
        for ip in _get_local_ips():
            print(f"  -> http://{ip}:{args.port}{token_suffix}")

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
