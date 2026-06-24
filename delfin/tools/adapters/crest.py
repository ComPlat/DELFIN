"""CREST conformer search adapter (native crest binary, full feature surface)."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, List, Optional, Sequence

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec

# CREST GFN level → CLI flag (gfn2 is the default, so no flag).
_GFN_FLAGS = {"gfn2": [], "gfn1": ["--gfn1"], "gfnff": ["--gfnff"], "gfn2//gfnff": ["--gfn2//gfnff"]}


def build_crest_cmd(
    geom_name: str,
    *,
    charge: int = 0,
    uhf: int = 0,
    cores: int = 1,
    ewin: float = 6.0,
    method: str = "gfn2",
    solvent: Optional[str] = None,
    extra_args: Sequence[str] = (),
) -> List[str]:
    """Assemble the ``crest`` command line (pure, so it is unit-testable)."""
    cmd: List[str] = ["crest", geom_name, "--chrg", str(int(charge)),
                      "--uhf", str(int(uhf)), "-T", str(int(cores)),
                      "--ewin", str(ewin)]
    cmd += _GFN_FLAGS.get(str(method).lower(), [])
    if solvent:
        cmd += ["--alpb", solvent]
    cmd += [str(a) for a in extra_args]
    return cmd


def _count_xyz_frames(path: Path) -> Optional[int]:
    """Number of structures in a multi-frame XYZ (ensemble) file."""
    try:
        lines = path.read_text().splitlines()
        natom = int(lines[0].strip())
        block = natom + 2
        if block > 0 and len(lines) % block == 0:
            return len(lines) // block
    except (ValueError, IndexError, OSError):
        return None
    return None


def _parse_crest_outputs(work_dir: Path) -> dict:
    """Conformer count + lowest energy from a finished CREST run."""
    data: dict = {}
    n = _count_xyz_frames(work_dir / "crest_conformers.xyz")
    if n is not None:
        data["n_conformers"] = n
    # The absolute lowest-conformer energy is in crest_best.xyz's comment line.
    try:
        comment = (work_dir / "crest_best.xyz").read_text().splitlines()[1]
        for tok in comment.split():
            try:
                data["lowest_energy_Eh"] = float(tok)
                break
            except ValueError:
                continue
    except (IndexError, OSError):
        pass
    return data


class CrestConformersAdapter(StepAdapter):
    name = "crest_conformers"
    description = "CREST conformer/isomer search (GFN0/1/2/FF, ALPB solvation)"
    produces_geometry = True
    category = "semiempirical"
    params = (
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1)"),
        ParamSpec("method", "str", default="gfn2", enum=tuple(_GFN_FLAGS),
                  description="GFN level for the search"),
        ParamSpec("solvent", "str", description="ALPB implicit solvent name"),
        ParamSpec("ewin", "float", default=6.0, unit="kcal/mol",
                  description="Energy window for retained conformers"),
        ParamSpec("extra_args", "list", description="Extra CLI args passed to CREST "
                  "(e.g. --nci, --tautomerize, --protonate, --qcg)"),
    )
    consumes = ("geometry",)
    produces = ("ensemble",)   # geometry is added automatically (produces_geometry)
    data_keys = (
        DataKeySpec("n_conformers", "int", "", "Number of conformers found"),
        DataKeySpec("lowest_energy_Eh", "float", "Eh", "Lowest conformer energy"),
    )
    requires_binaries = ("crest",)

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        local_geom = self._copy_geometry_to_workdir(geometry, work_dir, "input.xyz")
        cmd = build_crest_cmd(
            local_geom.name, charge=int(kwargs["charge"]),
            uhf=int(kwargs.get("mult", 1)) - 1, cores=cores,
            ewin=kwargs.get("ewin", 6.0), method=kwargs.get("method", "gfn2"),
            solvent=kwargs.get("solvent", ""), extra_args=kwargs.get("extra_args", []),
        )

        from delfin.qm_runtime import run_tool
        env = {"OMP_NUM_THREADS": str(cores), "OMP_STACKSIZE": "1G"}
        result = run_tool(cmd[0], cmd[1:], cwd=str(work_dir), env=env, capture_output=True)

        out_path = work_dir / "crest.out"
        out_path.write_text(result.stdout or "")

        best_xyz = work_dir / "crest_best.xyz"
        ensemble = work_dir / "crest_conformers.xyz"
        artifacts: dict[str, Path] = {}
        if ensemble.is_file():
            artifacts["ensemble"] = ensemble

        if result.returncode == 0 and best_xyz.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=best_xyz, output_file=out_path, artifacts=artifacts,
                data=_parse_crest_outputs(work_dir),
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            output_file=out_path, error=f"CREST failed (rc={result.returncode})",
        )


register(CrestConformersAdapter())
