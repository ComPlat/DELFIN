"""Workflow wrapper for CO2 Coordinator."""

from __future__ import annotations

import sys
from typing import Any, Dict, List

from delfin.workflows.registry import register


class Co2CoordinatorWorkflow:
    """CO2 Coordinator: automated CO2 placement, distance/rotation scans."""

    name = "co2_coordinator"
    description = "CO2 coordination: placement, distance scans, rotation scans"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.co2.CO2_Coordinator6 import main as co2_main

        return co2_main()

    def run_cli(self, argv: List[str]) -> int:
        from delfin.co2.CO2_Coordinator6 import main as co2_main

        # CO2 main() reads from its own CONTROL.txt via read_control_file()
        try:
            co2_main()
            return 0
        except SystemExit as e:
            return e.code if isinstance(e.code, int) else 1
        except Exception:
            return 1


register(Co2CoordinatorWorkflow())
