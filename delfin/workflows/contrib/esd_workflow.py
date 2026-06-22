"""Workflow wrapper for ESD (Excited-State Dynamics)."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class EsdWorkflow:
    """ESD: excited-state dynamics (ISC, IC, fluorescence, phosphorescence)."""

    name = "esd"
    description = "Excited-state dynamics: ISC/RISC rates, fluorescence, phosphorescence"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.esd_module import execute_esd_jobs

        return execute_esd_jobs(
            config=config,
            charge=kwargs.get("charge", int(config.get("charge", 0))),
            solvent=kwargs.get("solvent", config.get("solvent", "")),
            metals=kwargs.get("metals", config.get("metals", [])),
            main_basisset=kwargs.get("main_basisset", config.get("main_basisset", "def2-SVP")),
            metal_basisset=kwargs.get("metal_basisset", config.get("metal_basisset", "")),
        )

    def run_cli(self, argv: List[str]) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("--charge", type=int, default=0)
        parser.add_argument("--solvent", type=str, default="")
        parser.add_argument("--metals", type=str, nargs="*", default=[])
        parser.add_argument("--basis", type=str, default="def2-SVP")
        parser.add_argument("--metal-basis", type=str, default="")
        args = parser.parse_args(argv)

        from delfin.config import read_control_file
        config = read_control_file()

        from delfin.esd_module import execute_esd_jobs
        try:
            execute_esd_jobs(
                config=config,
                charge=args.charge,
                solvent=args.solvent,
                metals=args.metals,
                main_basisset=args.basis,
                metal_basisset=args.metal_basis,
            )
            return 0
        except Exception:
            return 1


register(EsdWorkflow())
