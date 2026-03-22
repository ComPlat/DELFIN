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
import subprocess
import sys
from pathlib import Path


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


def _prepare_voila_env(open_browser: bool) -> dict[str, str]:
    env = os.environ.copy()
    if not open_browser:
        return env

    browser = str(env.get("BROWSER") or "")
    term_program = str(env.get("TERM_PROGRAM") or "")
    if term_program != "vscode" and "browser.sh" not in browser:
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


def _voila_is_available() -> bool:
    """Return True when the current Python can import voila."""
    return importlib.util.find_spec("voila") is not None


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
        help="Explicitly ask Voila to open a browser window",
    )
    browser_group.add_argument(
        "--no-browser",
        action="store_true",
        help="Default behaviour (kept for backwards compatibility)",
    )
    parser.add_argument(
        "--dark",
        action="store_true",
        help="Use dark theme",
    )
    parser.add_argument(
        "--ip",
        default="0.0.0.0",
        help="IP to bind to (default: 0.0.0.0 = all interfaces, use 127.0.0.1 for local only)",
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

    open_browser = bool(args.open_browser)
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
        "--VoilaConfiguration.preheat_kernels=True",
        "--Voila.tornado_settings=disable_check_xsrf=True",
    ]

    if open_browser:
        cmd.append("--Voila.open_browser=True")
    else:
        cmd.append("--no-browser")

    if args.dark:
        cmd.append("--theme=dark")

    bind_display = "localhost" if args.ip == "127.0.0.1" else args.ip
    print(f"Starting DELFIN Dashboard on http://{bind_display}:{args.port}")

    # Show all reachable URLs so remote users know exactly what to type.
    if args.ip == "0.0.0.0":
        for ip in _get_local_ips():
            print(f"  -> http://{ip}:{args.port}")

    print()
    if args.ip == "0.0.0.0":
        print(
            "Tip: From another machine, use one of the URLs above.\n"
            "     If it doesn't connect, check your firewall:\n"
            f"       sudo ufw allow {args.port}/tcp"
        )
    print("Press Ctrl+C to stop.\n")

    proc = subprocess.Popen(cmd, env=env)
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
