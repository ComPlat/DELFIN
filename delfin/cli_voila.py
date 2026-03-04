"""``delfin-voila`` – launch the DELFIN Dashboard via Voila.

Usage::

    delfin-voila              # start on default port 8866
    delfin-voila --port 9000  # custom port
    delfin-voila --no-browser # don't auto-open a browser
"""

import argparse
import importlib.resources
import shutil
import subprocess
import sys


def _find_notebook() -> str:
    """Locate the bundled ``delfin_dashboard.ipynb``."""
    ref = importlib.resources.files("delfin.notebooks").joinpath(
        "delfin_dashboard.ipynb"
    )
    # For editable installs the path is already a real file;
    # for wheel installs importlib.resources may return a traversable.
    path = str(ref)
    return path


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

    cmd = [
        sys.executable,
        "-m",
        "voila",
        notebook,
        f"--port={args.port}",
        "--show_tracebacks=True",
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
        proc = subprocess.run(cmd, check=False)
        sys.exit(proc.returncode)
    except KeyboardInterrupt:
        print("\nStopped.")


if __name__ == "__main__":
    main()
