"""ORCA calculation adapters (single-point, optimization, frequencies).

Supports the full ORCA feature set used by DELFIN workflows:

Simple usage::

    run_step("orca_sp", geometry="mol.xyz", charge=0, method="B3LYP", basis="def2-SVP")

Full-featured usage (matching existing DELFIN workflow capabilities)::

    run_step("orca_opt", geometry="mol.xyz",
        charge=2, mult=3,
        method="B3LYP", basis="def2-SVP",
        ri="RIJCOSX", aux_basis="def2/J",
        dispersion="D4",
        solvent="water", solvent_model="CPCM",
        metal_basis={"Fe": "def2-TZVP", "Cu": "def2-TZVP"},
        opt_level="TIGHTOPT",
        moread="previous.gbw",
        broken_sym="%scf BrokenSym 5,3 {4,5,6} {7,8,9} end",
        scf_maxiter=300,
        scf_extra="SOSCF true",
        extra_blocks="%geom TolE 1e-5 end",
        base_name="S0",
    )
"""

from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools.adapters.xtb import _read_xyz_coords


# ──────────────────────────────────────────────────────────────────────
#  Input builder
# ──────────────────────────────────────────────────────────────────────

def _assemble_bang(
    job_type: str,
    method: str,
    basis: str,
    *,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    solvent: str = "",
    solvent_model: str = "CPCM",
    relativity: str = "",
    opt_level: str = "",
    extra_keywords: str = "",
) -> str:
    """Assemble the ORCA `!` method line from named parameters.

    If the caller already supplies a complete ``bang`` line, this function
    is bypassed entirely.
    """
    parts = ["!"]

    # Method + basis
    parts.append(method)
    parts.append(basis)

    # Job type (SP, OPT, FREQ, etc.)
    if job_type:
        if opt_level and job_type == "OPT":
            parts.append(opt_level)  # e.g. TIGHTOPT replaces OPT
        else:
            parts.append(job_type)

    # RI acceleration
    if ri:
        parts.append(ri)
    if aux_basis:
        parts.append(aux_basis)

    # Dispersion
    if dispersion:
        parts.append(dispersion)

    # Solvation
    if solvent:
        parts.append(f"{solvent_model}({solvent})")

    # Relativity
    if relativity:
        parts.append(relativity)

    # Extra keywords (e.g. "MOREAD", "TDDFT NRoots 10")
    if extra_keywords:
        parts.append(extra_keywords)

    return " ".join(parts)


def _build_full_orca_input(
    bang_line: str,
    charge: int,
    mult: int,
    cores: int,
    maxcore: int,
    coords: str,
    *,
    base_name: str = "",
    moread: str = "",
    metal_basis: Optional[Dict[str, str]] = None,
    broken_sym: str = "",
    scf_maxiter: Optional[int] = None,
    scf_extra: str = "",
    freq_temp: Optional[float] = None,
    output_blocks: str = "",
    extra_blocks: str = "",
    xyzfile_mode: bool = False,
    xyzfile_path: str = "",
) -> str:
    """Build a full ORCA input with all supported features."""
    parts = []

    # --- Method line ---
    if not bang_line.startswith("!"):
        bang_line = "! " + bang_line
    # Add MOREAD to method line if gbw file specified
    if moread and "MOREAD" not in bang_line.upper():
        bang_line = bang_line.rstrip() + " MOREAD"
    parts.append(bang_line.rstrip() + "\n")

    # --- %base ---
    if base_name:
        parts.append(f'%base "{base_name}"\n')

    # --- %pal + %maxcore ---
    parts.append(f"%maxcore {maxcore}\n")
    parts.append(f"%pal nprocs {cores} end\n")

    # --- %moinp (MOREAD) ---
    if moread:
        moread_path = moread if moread.endswith(".gbw") else moread + ".gbw"
        parts.append(f'%moinp "{moread_path}"\n')

    # --- %scf ---
    scf_lines = []
    if scf_maxiter is not None:
        scf_lines.append(f"  maxiter {scf_maxiter}")
    if scf_extra:
        for line in scf_extra.strip().splitlines():
            scf_lines.append(f"  {line.strip()}")
    if scf_lines:
        parts.append("%scf\n" + "\n".join(scf_lines) + "\nend\n")

    # --- %basis (per-atom) ---
    if metal_basis:
        basis_lines = []
        for element, bs in metal_basis.items():
            basis_lines.append(f'  NewGTO {element} "{bs}" end')
        parts.append("%basis\n" + "\n".join(basis_lines) + "\nend\n")

    # --- %freq ---
    if freq_temp is not None:
        parts.append(f"%freq\n  Temp {freq_temp}\nend\n")

    # --- Broken symmetry ---
    if broken_sym:
        stripped = broken_sym.strip()
        if stripped:
            parts.append(stripped + "\n")

    # --- %output ---
    if output_blocks:
        parts.append(output_blocks.rstrip() + "\n")

    # --- Extra blocks (anything else) ---
    if extra_blocks:
        parts.append(extra_blocks.rstrip() + "\n")

    # --- Geometry ---
    if xyzfile_mode and xyzfile_path:
        parts.append(f"* xyzfile {charge} {mult} {xyzfile_path}\n")
    else:
        parts.append(f"*xyz {charge} {mult}\n{coords}\n*\n")

    return "".join(parts)


# ──────────────────────────────────────────────────────────────────────
#  Energy extraction
# ──────────────────────────────────────────────────────────────────────

def _extract_energy(out_path: Path) -> Optional[float]:
    """Try to extract the final electronic energy from an ORCA output."""
    try:
        from delfin.energies import find_gibbs_energy
        e = find_gibbs_energy(str(out_path))
        if e is not None:
            return e
    except Exception:
        pass
    try:
        text = out_path.read_text()
        for line in reversed(text.splitlines()):
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])
    except Exception:
        pass
    return None


def _extract_gibbs(out_path: Path) -> Optional[float]:
    """Extract Gibbs free energy if available (freq calc)."""
    try:
        text = out_path.read_text()
        for line in reversed(text.splitlines()):
            if "Final Gibbs free energy" in line or "Final Gibbs free enthalpy" in line:
                return float(line.split()[-2])
    except Exception:
        pass
    return None


# ──────────────────────────────────────────────────────────────────────
#  Base adapter
# ──────────────────────────────────────────────────────────────────────

class _OrcaBase(StepAdapter):
    """Shared logic for all ORCA adapters."""

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def _run(
        self,
        work_dir: Path,
        job_type: str,
        geometry: Optional[Path],
        cores: int,
        start: float,
        **kwargs: Any,
    ) -> StepResult:
        # Pop internal pipeline keys that shouldn't leak into ORCA logic
        kwargs.pop("_prev_artifacts", None)

        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        maxcore = kwargs.get("maxcore", 1000)
        timeout = kwargs.get("timeout", None)
        calc_name = kwargs.get("base_name", "calc")

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        # --- Assemble bang line ---
        if "bang" in kwargs:
            bang_line = kwargs["bang"]
        else:
            bang_line = _assemble_bang(
                job_type,
                method=kwargs.get("method", "B3LYP"),
                basis=kwargs.get("basis", "def2-SVP"),
                ri=kwargs.get("ri", ""),
                aux_basis=kwargs.get("aux_basis", ""),
                dispersion=kwargs.get("dispersion", ""),
                solvent=kwargs.get("solvent", ""),
                solvent_model=kwargs.get("solvent_model", "CPCM"),
                relativity=kwargs.get("relativity", ""),
                opt_level=kwargs.get("opt_level", ""),
                extra_keywords=kwargs.get("extra_keywords", ""),
            )

        # --- Handle MOREAD: copy GBW file into work_dir ---
        moread = kwargs.get("moread", "")
        if moread:
            moread_src = Path(moread)
            if not moread_src.suffix:
                moread_src = moread_src.with_suffix(".gbw")
            if moread_src.is_file():
                dest = work_dir / moread_src.name
                if not dest.exists():
                    shutil.copy2(moread_src, dest)
                moread = moread_src.stem  # just the name without .gbw

        coords = _read_xyz_coords(geometry)

        inp_content = _build_full_orca_input(
            bang_line, charge, mult, cores, maxcore, coords,
            base_name=kwargs.get("base_name", ""),
            moread=moread,
            metal_basis=kwargs.get("metal_basis"),
            broken_sym=kwargs.get("broken_sym", ""),
            scf_maxiter=kwargs.get("scf_maxiter"),
            scf_extra=kwargs.get("scf_extra", ""),
            freq_temp=kwargs.get("freq_temp"),
            output_blocks=kwargs.get("output_blocks", ""),
            extra_blocks=kwargs.get("extra_blocks", ""),
        )

        inp_path = work_dir / f"{calc_name}.inp"
        out_path = work_dir / f"{calc_name}.out"
        inp_path.write_text(inp_content)

        # --- Smart recalc: skip if input unchanged and output complete ---
        from delfin.smart_recalc import should_skip
        if should_skip(inp_path, out_path):
            logger.info("Smart recalc: skipping %s (unchanged input, complete output)", calc_name)
            success = True
        else:
            # --- Run ORCA with intelligent error recovery if config provided ---
            config = kwargs.get("config") or kwargs.get("_prev_artifacts", {}).get("_config")
            if config and config.get("enable_auto_recovery", "").lower() in ("yes", "true", "1"):
                from delfin.orca import run_orca_with_intelligent_recovery
                success = run_orca_with_intelligent_recovery(
                    str(inp_path), str(out_path), timeout=timeout,
                    working_dir=work_dir, config=config,
                )
            else:
                from delfin.orca import run_orca
                success = run_orca(str(inp_path), str(out_path), timeout=timeout, working_dir=work_dir)

        # --- Collect results ---
        data: dict[str, Any] = {}
        if out_path.is_file():
            energy = _extract_energy(out_path)
            if energy is not None:
                data["energy_Eh"] = energy
            gibbs = _extract_gibbs(out_path)
            if gibbs is not None:
                data["gibbs_Eh"] = gibbs

        # Look for output geometry (base_name.xyz or calc.xyz)
        xyz_path = work_dir / f"{calc_name}.xyz"
        geom_out = xyz_path if xyz_path.is_file() else None

        # Collect artifacts
        artifacts: dict[str, Path] = {}
        for ext, key in [(".gbw", "gbw"), (".hess", "hess"), (".densitiesinfo", "densities")]:
            f = work_dir / f"{calc_name}{ext}"
            if f.is_file():
                artifacts[key] = f

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir,
            start,
            geometry=geom_out,
            output_file=out_path if out_path.is_file() else None,
            data=data,
            artifacts=artifacts,
            error=None if success else "ORCA calculation failed",
        )


# ──────────────────────────────────────────────────────────────────────
#  Concrete adapters
# ──────────────────────────────────────────────────────────────────────

class OrcaSpAdapter(_OrcaBase):
    name = "orca_sp"
    description = "ORCA single-point energy calculation"
    produces_geometry = False

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        return self._run(work_dir, "SP", geometry, cores, time.monotonic(), **kwargs)


class OrcaOptAdapter(_OrcaBase):
    name = "orca_opt"
    description = "ORCA geometry optimization"
    produces_geometry = True

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        return self._run(work_dir, "OPT", geometry, cores, time.monotonic(), **kwargs)


class OrcaFreqAdapter(_OrcaBase):
    name = "orca_freq"
    description = "ORCA frequency calculation"
    produces_geometry = False

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        return self._run(work_dir, "FREQ", geometry, cores, time.monotonic(), **kwargs)


class OrcaOptFreqAdapter(_OrcaBase):
    name = "orca_opt_freq"
    description = "ORCA combined geometry optimization + frequency calculation"
    produces_geometry = True

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        return self._run(work_dir, "OPT FREQ", geometry, cores, time.monotonic(), **kwargs)


class OrcaTddftAdapter(_OrcaBase):
    name = "orca_tddft"
    description = "ORCA TD-DFT excited state calculation"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        super().validate_params(**kwargs)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        # Inject TDDFT + NRoots into extra_keywords
        nroots = kwargs.pop("nroots", 10)
        extra_kw = kwargs.get("extra_keywords", "")
        if "TDDFT" not in extra_kw.upper():
            extra_kw = f"TDDFT NRoots {nroots} " + extra_kw
        kwargs["extra_keywords"] = extra_kw.strip()
        return self._run(work_dir, "SP", geometry, cores, time.monotonic(), **kwargs)


register(OrcaSpAdapter())
register(OrcaOptAdapter())
register(OrcaFreqAdapter())
register(OrcaOptFreqAdapter())
register(OrcaTddftAdapter())
