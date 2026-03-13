"""``delfin-voila`` – launch the DELFIN Dashboard via Voila.

Usage::

    delfin-voila              # start on default port 8866
    delfin-voila --port 9000  # custom port
    delfin-voila --no-browser # don't auto-open a browser
"""

import argparse
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
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not automatically open a browser window",
    )
    parser.add_argument(
        "--dark",
        action="store_true",
        help="Use dark theme",
    )
    args = parser.parse_args(argv)

    # Check that voila is installed
    if shutil.which("voila") is None:
        print(
            "Error: voila is not installed.\n"
            "Install it with:  pip install 'delfin-complat[dashboard]'",
            file=sys.stderr,
        )
        sys.exit(1)

    notebook = _find_notebook()
    root_dir = str(Path.home().resolve())
    env = _prepare_voila_env(open_browser=not args.no_browser)
    env.setdefault("DELFIN_VOILA_ROOT_DIR", root_dir)

    cmd = [
        sys.executable,
        "-m",
        "voila",
        notebook,
        f"--port={args.port}",
        "--show_tracebacks=True",
        f"--Voila.root_dir={root_dir}",
        "--ServerApp.websocket_ping_interval=30000",
        "--ServerApp.websocket_ping_timeout=30000",
        "--VoilaConfiguration.file_allowlist=.*\\.(png|jpg|gif|svg)",
        "--VoilaConfiguration.file_allowlist=.*\\.(js|css|html)",
    ]

    if args.no_browser:
        cmd.append("--no-browser")
    else:
        cmd.append("--Voila.open_browser=True")

    if args.dark:
        cmd.append("--theme=dark")

    print(f"Starting DELFIN Dashboard on http://localhost:{args.port}")
    print("Press Ctrl+C to stop.\n")

    try:
        proc = subprocess.run(
            cmd,
            env=env,
            check=False,
        )
        sys.exit(proc.returncode)
    except KeyboardInterrupt:
        print("\nStopped.")


if __name__ == "__main__":
    main()
