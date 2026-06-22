"""Workflow wrapper for TADF-xTB screening."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class TadfXtbWorkflow:
    """TADF screening via xTB excited-state calculations."""

    name = "tadf_xtb"
    description = "TADF screening via xTB excited-state calculations"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.tadf_xtb import run_single_tadf_xtb

        return run_single_tadf_xtb(**kwargs)

    def run_cli(self, argv: List[str]) -> int:
        from delfin.tadf_xtb import run_cli

        return run_cli(argv)


register(TadfXtbWorkflow())
