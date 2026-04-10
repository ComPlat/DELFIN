"""Entry point: ``python -m delfin.doc_server`` runs the MCP server."""

from __future__ import annotations

import sys


def main() -> None:
    from .server import run_server

    run_server(sys.argv[1:])


if __name__ == "__main__":
    main()
