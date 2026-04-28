"""Convenience entry: ``python -m tests.agent_dryrun ...`` runs the report."""

from .report import main


if __name__ == "__main__":
    raise SystemExit(main())
