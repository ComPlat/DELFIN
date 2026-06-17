"""Multiwfn wavefunction analysis adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec


class MultiwfnBondOrdersAdapter(StepAdapter):
    name = "multiwfn_bond_orders"
    description = "Bond order analysis via Multiwfn"
    produces_geometry = False
    category = "analysis"
    params = (
        ParamSpec("input_file", "path", required=True,
                  description="Wavefunction file (.molden/.wfn/.wfx/.fch/.gbw)"),
        ParamSpec("method", "str", default="mayer", enum=("mayer", "wiberg", "fuzzy"),
                  description="Bond-order method"),
    )
    consumes = ("molden", "wfn")
    data_keys = (
        DataKeySpec("bond_orders", "list", "", "Per-bond order list"),
        DataKeySpec("output", "str", "", "Raw Multiwfn output"),
        DataKeySpec("method", "str"),
    )
    requires_binaries = ("Multiwfn",)

    def validate_params(self, **kwargs: Any) -> None:
        if "input_file" not in kwargs:
            raise ValueError("'input_file' parameter is required (path to .molden/.wfn/.gbw)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis

        input_file = kwargs["input_file"]
        method = kwargs.get("method", "mayer")

        result_dict = bond_order_analysis(input_file, method=method)

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=result_dict,
        )


register(MultiwfnBondOrdersAdapter())
