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

        # F6 (OC1 audit): forward the intended run folder. run_OCCUPIER takes
        # work_dir (delfin/occupier.py:680) precisely so a stage runs in its
        # own folder instead of the process CWD — dropping it silently re-runs
        # OCCUPIER where the caller did not ask for it.
        work_dir = config.get("work_dir") if config else None
        if work_dir is None:
            work_dir = kwargs.get("work_dir")
        if work_dir is None:
            run_dir = config.get("run_dir") if config else None
            if run_dir:
                work_dir = run_dir
        return run_OCCUPIER(work_dir=work_dir)

    def run_cli(self, argv: List[str]) -> int:
        from delfin.occupier import run_OCCUPIER

        try:
            run_OCCUPIER()
            return 0
        except Exception:
            return 1


register(OccupierWorkflow())
