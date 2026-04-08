"""Workflow wrapper for IMAG (imaginary frequency elimination)."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class ImagWorkflow:
    """IMAG: iterative elimination of imaginary frequencies."""

    name = "imag"
    description = "Iterative imaginary frequency elimination via ORCA re-optimization"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.imag import run_IMAG

        return run_IMAG(
            input_file=kwargs["input_file"],
            hess_file=kwargs["hess_file"],
            charge=kwargs["charge"],
            multiplicity=kwargs["mult"],
            solvent=kwargs.get("solvent", config.get("solvent", "")),
            metals=kwargs.get("metals", config.get("metals", [])),
            config=config,
            main_basisset=kwargs.get("main_basisset", config.get("main_basisset", "def2-SVP")),
            metal_basisset=kwargs.get("metal_basisset", config.get("metal_basisset", "")),
            broken_sym=kwargs.get("broken_sym", config.get("broken_sym", False)),
            step_name=kwargs.get("step_name", "imag"),
            pal_override=kwargs.get("cores"),
        )

    def run_cli(self, argv: List[str]) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("input_file", help="Path to input XYZ/inp file")
        parser.add_argument("hess_file", help="Path to .hess file")
        parser.add_argument("--charge", type=int, required=True)
        parser.add_argument("--mult", type=int, required=True)
        parser.add_argument("--solvent", type=str, default="")
        parser.add_argument("--metals", type=str, nargs="*", default=[])
        parser.add_argument("--basis", type=str, default="def2-SVP")
        parser.add_argument("--metal-basis", type=str, default="")
        args = parser.parse_args(argv)

        from delfin.imag import run_IMAG
        try:
            run_IMAG(
                input_file=args.input_file,
                hess_file=args.hess_file,
                charge=args.charge,
                multiplicity=args.mult,
                solvent=args.solvent,
                metals=args.metals,
                config={"PAL": 1},
                main_basisset=args.basis,
                metal_basisset=args.metal_basis,
                broken_sym=False,
            )
            return 0
        except Exception:
            return 1


register(ImagWorkflow())
