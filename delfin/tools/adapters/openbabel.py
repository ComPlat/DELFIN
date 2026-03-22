"""Open Babel adapter — format conversion, conformer generation, property calculation.

Requires ``openbabel`` Python bindings (``pip install openbabel-wheel`` or
system Open Babel with Python bindings).

Usage::

    run_step("obabel_convert", geometry="mol.xyz", output_format="mol2")
    run_step("obabel_convert", input_file="mol.sdf", output_format="xyz")
    run_step("obabel_conformers", geometry="mol.xyz", n_conformers=50)
    run_step("obabel_from_smiles", smiles="CCO", output_format="xyz")
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class ObabelConvertAdapter(StepAdapter):
    """Convert between molecular file formats using Open Babel."""

    name = "obabel_convert"
    description = "Convert molecular file formats via Open Babel (xyz, sdf, mol2, pdb, …)"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "output_format" not in kwargs:
            raise ValueError("'output_format' parameter is required (e.g. 'xyz', 'mol2', 'sdf')")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        input_file = kwargs.get("input_file")
        if input_file:
            input_path = Path(input_file).resolve()
        elif geometry:
            input_path = geometry
        else:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry or input_file is required")

        if not input_path.is_file():
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"Input file not found: {input_path}")

        output_format = kwargs["output_format"]
        out_path = work_dir / f"converted.{output_format}"

        try:
            from openbabel import openbabel as ob

            conv = ob.OBConversion()
            in_fmt = input_path.suffix.lstrip(".")
            conv.SetInFormat(in_fmt)
            conv.SetOutFormat(output_format)

            mol = ob.OBMol()
            conv.ReadFile(mol, str(input_path))

            if kwargs.get("add_hydrogens", False):
                mol.AddHydrogens()
            if kwargs.get("gen3d", False):
                builder = ob.OBBuilder()
                builder.Build(mol)

            conv.WriteFile(mol, str(out_path))

        except ImportError:
            # Fallback to obabel CLI
            import subprocess
            cmd = ["obabel", str(input_path), "-O", str(out_path)]
            if kwargs.get("add_hydrogens"):
                cmd.append("-h")
            if kwargs.get("gen3d"):
                cmd.append("--gen3d")
            proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
            if proc.returncode != 0:
                return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                         error=f"obabel failed: {proc.stderr[:500]}")

        geom_out = out_path if out_path.suffix == ".xyz" and out_path.is_file() else None

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=geom_out,
            data={"input_format": input_path.suffix.lstrip("."),
                  "output_format": output_format,
                  "output_file": str(out_path)},
            artifacts={"converted": out_path} if out_path.is_file() else {},
        )


class ObabelConformersAdapter(StepAdapter):
    """Generate conformers using Open Babel's conformer search."""

    name = "obabel_conformers"
    description = "Conformer generation via Open Babel (genetic algorithm / systematic)"
    produces_geometry = True

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        n_conformers = kwargs.get("n_conformers", 50)
        out_path = work_dir / "conformers.xyz"

        try:
            import subprocess
            cmd = [
                "obabel", str(geometry), "-O", str(out_path),
                "--conformer",
                "--nconf", str(n_conformers),
                "--score", "energy",
                "--writeconformers",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
            if proc.returncode != 0:
                return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                         error=f"obabel conformer search failed: {proc.stderr[:500]}")
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(e))

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_path if out_path.is_file() else None,
            data={"n_conformers_requested": n_conformers},
        )


class ObabelFromSmilesAdapter(StepAdapter):
    """Generate 3D structure from SMILES via Open Babel."""

    name = "obabel_from_smiles"
    description = "SMILES to 3D geometry via Open Babel (alternative to RDKit)"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "smiles" not in kwargs:
            raise ValueError("'smiles' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        smiles = kwargs["smiles"]
        output_format = kwargs.get("output_format", "xyz")
        out_path = work_dir / f"structure.{output_format}"

        try:
            import subprocess
            cmd = [
                "obabel", f"-:{smiles}", "-O", str(out_path),
                "--gen3d",
            ]
            if kwargs.get("add_hydrogens", True):
                cmd.append("-h")
            proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
            if proc.returncode != 0:
                return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                         error=f"obabel failed: {proc.stderr[:500]}")
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(e))

        geom_out = out_path if out_path.suffix == ".xyz" and out_path.is_file() else None

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=geom_out,
            data={"smiles": smiles},
        )


register(ObabelConvertAdapter())
register(ObabelConformersAdapter())
register(ObabelFromSmilesAdapter())
