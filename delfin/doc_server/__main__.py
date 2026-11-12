"""Entry point: ``python -m delfin.doc_server`` runs the MCP server."""

from __future__ import annotations

import sys


def main() -> None:
    from delfin.hpc_guard import exit_on_login_node
    exit_on_login_node("the DELFIN docs server", "it starts with an agent session inside a SLURM job")
    from .server import run_server

    run_server(sys.argv[1:])


if __name__ == "__main__":
    main()
