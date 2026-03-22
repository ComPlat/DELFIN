"""Workflow wrapper for OCCUPIER (adaptive spin-state / redox potential)."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class OccupierWorkflow:
    """OCCUPIER: automated spin-state identification and redox potential calculation."""

    name = "occupier"
    description = "Adaptive spin-state identification and redox potential calculation"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.occupier import run_OCCUPIER

        return run_OCCUPIER()

    def run_cli(self, argv: List[str]) -> int:
        from delfin.occupier import run_OCCUPIER

        try:
            run_OCCUPIER()
            return 0
        except Exception:
            return 1


register(OccupierWorkflow())
