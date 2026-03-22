"""Workflow wrapper for hyperpolarizability calculations."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class HyperpolWorkflow:
    """Hyperpolarizability calculations (beta tensor)."""

    name = "hyperpol"
    description = "Hyperpolarizability calculations (beta tensor)"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.hyperpol import run_single_hyperpol_workflow

        return run_single_hyperpol_workflow(**kwargs)

    def run_cli(self, argv: List[str]) -> int:
        from delfin.hyperpol import run_cli

        return run_cli(argv)


register(HyperpolWorkflow())
