"""Workflow wrapper for IMAG (imaginary frequency elimination)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from delfin.workflows.registry import register


class ImagWorkflow:
    """IMAG: iterative elimination of imaginary frequencies."""

    name = "imag"
    description = "Iterative imaginary frequency elimination via ORCA re-optimization"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        """IMAG on one finished calculation, asked for explicitly (no IMAG switch is read).

        ``input_file`` is the calculation's .inp or .out; the other one and the
        .hess are found beside it under the same name.
        """
        from delfin.imag import _pipeline_run_orca, eliminate_imaginary_modes

        given = Path(kwargs["input_file"])
        if given.suffix not in (".inp", ".out"):
            raise ValueError(f"IMAG needs the calculation's .inp or .out, not {given.name}")
        return eliminate_imaginary_modes(
            label=kwargs.get("step_name", "imag"),
            input_path=given.with_suffix(".inp"),
            output_path=given.with_suffix(".out"),
            config=config,
            run_orca=_pipeline_run_orca,
            pal=kwargs.get("cores"),
        )

    def run_cli(self, argv: List[str]) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("input_file", help="The calculation's .inp or .out (the .hess beside it)")
        parser.add_argument("--cores", type=int, default=None)
        parser.add_argument("--max-rounds", type=int, default=None)
        args = parser.parse_args(argv)
        config: Dict[str, Any] = {}
        if args.max_rounds:
            config["IMAG_max_rounds"] = args.max_rounds
        try:
            result = self.run(config=config, input_file=args.input_file, cores=args.cores)
        except Exception:  # noqa: BLE001
            return 1
        return 0 if result.resolved else 1


register(ImagWorkflow())
